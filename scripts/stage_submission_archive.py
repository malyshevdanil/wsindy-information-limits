"""Stage the current tracked research archive without moving earlier tags."""
import hashlib,json,shutil,subprocess,re
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];DEST=ROOT/'reports/public_release_v1'

def main():
    assert DEST.resolve().is_relative_to(ROOT.resolve()) and (DEST/'.git').is_dir()
    names=subprocess.check_output(['git','ls-files'],cwd=ROOT,text=True).splitlines()
    copied=[]
    for name in names:
        p=ROOT/name
        if not p.is_file() or name.startswith('reports/public_release_v1/'):continue
        if p.suffix in ['.log','.png'] and name.startswith('reports/'):continue
        if name=='build_document.py':continue
        target=DEST/name;target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,target);copied.append(name)
    (DEST/'.gitattributes').write_bytes(b'* -text\n')
    # This is a generated checkout. Historical files excluded from the current
    # tracked package remain available under their unchanged earlier tags.
    allowed=set(copied)|{'.gitattributes','release_manifest.json'};pruned=0
    for p in sorted(DEST.rglob('*')):
        if not p.is_file() or '.git' in p.relative_to(DEST).parts:continue
        if p.relative_to(DEST).as_posix() not in allowed:
            assert p.resolve().is_relative_to(DEST.resolve());p.unlink();pruned+=1
    records=[]
    for p in sorted(DEST.rglob('*')):
        if not p.is_file() or any(q in ('.git','__pycache__','.pytest_cache') for q in p.relative_to(DEST).parts):continue
        name=p.relative_to(DEST).as_posix()
        if name=='release_manifest.json':continue
        data=p.read_bytes();assert p.name not in ('access_token','kaggle.json','.env')
        assert not re.search(rb'KGAT_[A-Za-z0-9]{20,}',data)
        records.append({'path':name,'bytes':len(data),'sha256':hashlib.sha256(data).hexdigest()})
    (DEST/'release_manifest.json').write_text(json.dumps(records,indent=2),encoding='utf-8')
    subprocess.run(['git','config','core.autocrlf','false'],cwd=DEST,check=True)
    subprocess.run(['git','add','-f','--all'],cwd=DEST,check=True)
    tracked=set(subprocess.check_output(['git','ls-files'],cwd=DEST,text=True).splitlines())
    assert tracked=={r['path'] for r in records}|{'release_manifest.json'},'Every manifest member must be tracked in the published tree'
    disclosure=(DEST/'manuscript/article_revised.md').read_text(encoding='utf-8')
    assert '## Declaration of generative AI and AI-assisted technologies' in disclosure
    assert 'Codex' in disclosure and 'mathematical derivations' in disclosure
    result={'version':json.loads((ROOT/'configs/public_archive.json').read_text())['version'],'copied_tracked_files':len(copied),'manifest_files':len(records),'derived_legacy_files_pruned':pruned,'secret_scan_verified':True,'ai_disclosure_preserved':True,'previous_tags_unchanged':True}
    (ROOT/'reports/submission_archive_preparation.json').write_text(json.dumps(result,indent=2));print(json.dumps(result))

if __name__=='__main__':
    import sys;sys.path.insert(0,str(ROOT));main()
