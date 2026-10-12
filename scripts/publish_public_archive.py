"""Stage, publish and anonymously verify an authorized versioned GitHub archive."""
import argparse,hashlib,json,os,re,subprocess,time,zipfile
from pathlib import Path
import requests
ROOT=Path(__file__).resolve().parents[1]
CHECKOUT=ROOT/'reports/public_release_v1'
OWNER='malyshevdanil';REPO='wsindy-information-limits'


def config():return json.loads((ROOT/'configs/public_archive.json').read_text())


def git(*args,check=True):
    process=subprocess.run(['git',*args],cwd=CHECKOUT,capture_output=True,text=True,timeout=1200,
        env={**os.environ,'GIT_TERMINAL_PROMPT':'0','GCM_INTERACTIVE':'never'})
    if check and process.returncode:raise RuntimeError('Git operation failed: '+args[0])
    return process.stdout.strip()


def prepare():
    assert CHECKOUT.resolve().is_relative_to(ROOT.resolve()) and (CHECKOUT/'.git').is_dir()
    for name in ('reproducibility_code','reproducibility_data'):
        with zipfile.ZipFile(ROOT/'dist'/(name+'.zip')) as archive:
            assert archive.testzip() is None
            for entry in archive.infolist():
                path=(CHECKOUT/entry.filename).resolve()
                if not path.is_relative_to(CHECKOUT.resolve()):raise ValueError('Unsafe release member')
                assert not any(part in ('.git','.kaggle','author_repos','author_sources') for part in Path(entry.filename).parts)
                path.parent.mkdir(parents=True,exist_ok=True);path.write_bytes(archive.read(entry))
    for source in (ROOT/'manuscript/figures/publication').glob('*'):
        target=CHECKOUT/source.relative_to(ROOT);target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(source.read_bytes())
    for name in ('article_revised','supplement','weak_efficiency_theory','independent_records_theory','heterogeneous_records_theory'):
        for suffix in ('.docx','.pdf','.tex'):
            source=ROOT/'manuscript'/(name+suffix);(CHECKOUT/'manuscript'/source.name).write_bytes(source.read_bytes())
    for name in ('article_ejs.tex','article_ejs.pdf'):
        source=ROOT/'manuscript'/name;(CHECKOUT/'manuscript'/name).write_bytes(source.read_bytes())
    (CHECKOUT/'.gitattributes').write_bytes(b'* -text\n')
    records=[]
    for path in sorted(CHECKOUT.rglob('*')):
        if not path.is_file() or '.git' in path.relative_to(CHECKOUT).parts:continue
        name=path.relative_to(CHECKOUT).as_posix()
        if name=='release_manifest.json':continue
        data=path.read_bytes()
        assert path.name not in ('access_token','kaggle.json','.env')
        assert not re.search(rb'KGAT_[A-Za-z0-9]{20,}',data),f'Unexpected token-shaped data: {name}'
        records.append(dict(path=name,bytes=len(data),sha256=hashlib.sha256(data).hexdigest()))
    (CHECKOUT/'release_manifest.json').write_text(json.dumps(records,indent=2),encoding='utf-8')
    git('config','core.autocrlf','false')
    git('add','--all')
    tracked=git('ls-files').splitlines()
    assert set(tracked)=={r['path'] for r in records}|{'release_manifest.json'}
    (ROOT/'reports/public_archive_preparation.json').write_text(json.dumps(dict(version=config()['version'],
        staged_files=len(tracked),secret_scan=True,sha256_manifest=True),indent=2))
    print(json.dumps(dict(staged_files=len(tracked),version=config()['version'])))


def public_get(path):
    for attempt in range(6):
        try:
            response=requests.get('https://api.github.com/repos/'+OWNER+'/'+REPO+path,
                params={'verify':time.time_ns()},timeout=(20,60),headers={'Accept':'application/vnd.github+json'})
            if response.status_code==200:return response.json()
        except requests.RequestException:pass
        time.sleep(2)
    raise RuntimeError('Anonymous GitHub verification failed') from None


def publish():
    cfg=config();version=cfg['version']
    assert re.fullmatch(r'v\d+\.\d+\.\d+',version)
    git('config','user.name',OWNER);git('config','user.email',OWNER+'@users.noreply.github.com')
    if git('status','--porcelain'):
        git('-c','commit.gpgsign=false','commit','-m','Validate paired information limits and independent-record weak GMM')
    head=git('rev-parse','HEAD')
    tags=git('tag','--list',version).splitlines()
    if tags:assert git('rev-parse',version)==head,'Never move an existing release tag'
    else:git('tag',version)
    git('-c','http.version=HTTP/1.1','-c','http.postBuffer=104857600','push','--atomic','origin','main',version)
    print(json.dumps(dict(pushed=True,tag=version,commit=head)))


def verify():
    cfg=config();version=cfg['version']
    repository=public_get('');assert repository['private'] is False
    ref=public_get('/git/ref/tags/'+version);commit=ref['object']['sha'];assert ref['object']['type']=='commit'
    tree=public_get('/git/trees/'+commit+'?recursive=1')
    assert not tree['truncated']
    blobs={r['path']:r['sha'] for r in tree['tree'] if r['type']=='blob'}
    tracked=git('ls-files').splitlines();assert set(blobs)==set(tracked)
    for name in tracked:
        data=(CHECKOUT/name).read_bytes();actual=hashlib.sha1(b'blob '+str(len(data)).encode()+b'\0'+data).hexdigest()
        assert actual==blobs[name],f'Public byte mismatch: {name}'
    for r in json.loads((CHECKOUT/'release_manifest.json').read_text()):
        assert hashlib.sha256((CHECKOUT/r['path']).read_bytes()).hexdigest()==r['sha256']
    result=dict(published=True,repository=cfg['repository'],version_url=cfg['version_url'],tag=version,commit=commit,
        tracked_files=len(tracked),anonymous_access_verified=True,all_public_blob_hashes_verified=True,
        sha256_manifest_verified=True,source_and_data_verified=True,doi=None)
    (ROOT/'reports/public_archive_validation.json').write_text(json.dumps(result,indent=2))
    print(json.dumps(result))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('action',choices=['prepare','publish','verify'])
    globals()[parser.parse_args().action]()
