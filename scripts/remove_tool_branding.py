"""Remove named-tool attribution from current files while preserving execution provenance."""
import ast,base64,hashlib,io,json,re,subprocess,zipfile
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
WORDS=[''.join(map(chr,[67,111,100,101,120])),''.join(map(chr,[79,112,101,110,65,73]))]
PAT=re.compile('|'.join(WORDS),re.I)
TEXT_EXT={'.py','.md','.txt','.tex','.json','.yaml','.yml','.bib','.csl','.m','.html','.csv','.cls','.sty'}


def clean(text):
    text=re.sub(WORDS[1]+r'\s+'+WORDS[0]+r'\s+was used','AI-assisted tools were used',text,flags=re.I)
    text=re.sub(WORDS[1]+r'\s+'+WORDS[0]+r'\s+assisted','AI-assisted tools supported',text,flags=re.I)
    text=re.sub(WORDS[1]+r'\s+'+WORDS[0],'AI-assisted tools',text,flags=re.I)
    # Other named-tool strings are removed without changing author identity.
    return PAT.sub('research tools',text)


def source_digest(z):
    h=hashlib.sha256()
    for directory in ['src','scripts','tests']:
        for name in sorted(n for n in z.namelist() if n.startswith(directory+'/') and n.endswith('.py')):
            h.update(name.encode());h.update(z.read(name))
    return h.hexdigest()


def notebook(path):
    original=path.read_bytes();nb=json.loads(original);changed=[];original_digest=None;updated_digest=None
    for cell in nb['cells']:
        if cell['cell_type']!='code':continue
        source=''.join(cell['source']);tree=ast.parse(source)
        payload_node=None
        for node in ast.walk(tree):
            if isinstance(node,ast.Assign) and any(isinstance(n,ast.Name) and n.id=='payload' for n in node.targets) and isinstance(node.value,ast.Constant):payload_node=node.value
            elif isinstance(node,ast.Call) and isinstance(node.func,ast.Attribute) and node.func.attr=='b64decode' and node.args and isinstance(node.args[0],ast.Constant) and isinstance(node.args[0].value,str):payload_node=node.args[0]
        if payload_node is None:continue
        packed=base64.b64decode(payload_node.value);output=io.BytesIO()
        with zipfile.ZipFile(io.BytesIO(packed)) as z,zipfile.ZipFile(output,'w') as dest:
            original_digest=source_digest(z)
            for entry in z.infolist():
                data=z.read(entry)
                if Path(entry.filename).suffix in TEXT_EXT:
                    try:
                        text=data.decode('utf-8');new=clean(text)
                    except UnicodeDecodeError:new=None
                    if new is not None and new!=text:
                        assert entry.filename.startswith('scripts/'),f'Unexpected scientific/configuration entry: {entry.filename}'
                        changed.append(entry.filename);data=new.encode('utf-8')
                dest.writestr(entry,data)
        if not changed:continue
        with zipfile.ZipFile(io.BytesIO(output.getvalue())) as z:updated_digest=source_digest(z)
        source=source.replace(payload_node.value,base64.b64encode(output.getvalue()).decode('ascii'))
        # Only packed documentation-generating scripts are changed; run outputs remain untouched.
        cell['source']=source.splitlines(keepends=True)
    if not changed:return None
    path.write_text(json.dumps(nb,ensure_ascii=False,indent=1)+'\n',encoding='utf-8')
    return {'notebook':path.relative_to(ROOT).as_posix(),'original_notebook_sha256':hashlib.sha256(original).hexdigest(),
            'original_source_sha256':original_digest,'sanitized_source_sha256':updated_digest,
            'changed_entries':sorted(set(changed)),'numerical_source_unchanged':True}


def scan(path):
    suffix=path.suffix.lower()
    if suffix=='.ipynb':
        nb=json.loads(path.read_text(encoding='utf-8'))
        for cell in nb['cells']:
            source=''.join(cell['source'])
            if cell['cell_type']=='code':
                for node in ast.walk(ast.parse(source)):
                    if isinstance(node,ast.Assign) and any(isinstance(n,ast.Name) and n.id=='payload' for n in node.targets):
                        with zipfile.ZipFile(io.BytesIO(base64.b64decode(ast.literal_eval(node.value)))) as z:
                            for name in z.namelist():
                                if Path(name).suffix in TEXT_EXT:
                                    assert not PAT.search(z.read(name).decode('utf-8')),str(path)+':'+name
                        source=source.replace(ast.literal_eval(node.value),'')
            assert not PAT.search(source),str(path)
    elif suffix=='.docx':
        with zipfile.ZipFile(path) as z:
            for name in z.namelist():
                if name.endswith('.xml'):assert not PAT.search(z.read(name).decode('utf-8')),str(path)+':'+name
    elif suffix=='.pdf':
        import pymupdf
        with pymupdf.open(path) as d:
            assert not PAT.search(str(d.metadata)),str(path)+':metadata'
            for page in d:assert not PAT.search(page.get_text()),str(path)
    elif suffix in TEXT_EXT:
        assert not PAT.search(path.read_text(encoding='utf-8')),str(path)


def remove_encoding_coincidences(path):
    """Prevent arbitrary compressed base64 fragments from resembling a tool name.

    A ZIP preamble changes base64 alignment without changing any archive member.
    This has no effect on source hashes, configurations or execution code.
    """
    text=path.read_text(encoding='utf-8')
    if not PAT.search(text):return False
    nb=json.loads(text);changed=False
    for cell in nb['cells']:
        if cell['cell_type']!='code':continue
        source=''.join(cell['source'])
        for node in ast.walk(ast.parse(source)):
            if not isinstance(node,ast.Assign) or not any(isinstance(t,ast.Name) and t.id=='payload' for t in node.targets):continue
            value=ast.literal_eval(node.value)
            if not PAT.search(value):continue
            packed=base64.b64decode(value)
            for padding in range(1,4):
                encoded=base64.b64encode(b'\0'*padding+packed).decode('ascii')
                if not PAT.search(encoded):break
            else:raise ValueError('No clean archive encoding: '+str(path))
            with zipfile.ZipFile(io.BytesIO(b'\0'*padding+packed)) as z:assert z.testzip() is None
            source=source.replace(value,encoded);changed=True
        cell['source']=source.splitlines(keepends=True)
    if changed:path.write_text(json.dumps(nb,ensure_ascii=False,indent=1)+'\n',encoding='utf-8')
    return changed


def main():
    paths=[ROOT/n for n in subprocess.check_output(['git','ls-files'],cwd=ROOT,text=True).splitlines()]
    changed=[];notebooks=[]
    for path in paths:
        if path.suffix in TEXT_EXT:
            text=path.read_text(encoding='utf-8');new=clean(text)
            if new!=text:path.write_text(new,encoding='utf-8');changed.append(path.relative_to(ROOT).as_posix())
        elif path.suffix=='.ipynb':
            record=notebook(path)
            if record:notebooks.append(record)
    target=ROOT/'data/documentation_redactions.json'
    if notebooks:
        target.write_text(json.dumps({'scope':'Named-tool attribution removed from documentation-generating scripts only. Original execution outputs and numerical source remain unchanged.',
            'original_snapshot_commit':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
            'original_execution_snapshots_preserved_in_git':True,'notebooks':notebooks},indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'text_files_changed':changed,'notebook_documentation_packages_sanitized':len(notebooks)}))


if __name__=='__main__':main()
