"""Upload an already validated archive through GitHub's Git data API.

Used when the large git smart-HTTP upload stalls. Existing public tags are
never rewritten. Credentials remain exclusively in process memory.
"""
import argparse,base64,concurrent.futures,hashlib,json,os,subprocess,time
from pathlib import Path
import requests

ROOT=Path(__file__).resolve().parents[1];CHECKOUT=ROOT/'reports/public_release_v1'
BASE='https://api.github.com/repos/malyshevdanil/wsindy-information-limits'


def git(*args):
    p=subprocess.run(['git',*args],cwd=CHECKOUT,capture_output=True,text=True,timeout=120,
        env={**os.environ,'GIT_TERMINAL_PROMPT':'0','GCM_INTERACTIVE':'never'})
    if p.returncode:raise RuntimeError('Local archive git operation failed: '+args[0])
    return p.stdout.strip()


def main(metadata_only=False):
    cfg=json.loads((ROOT/'configs/public_archive.json').read_text());version=cfg['version']
    assert (CHECKOUT/'.git').is_dir() and CHECKOUT.resolve().is_relative_to(ROOT.resolve())
    p=subprocess.run(['git','credential','fill'],input='protocol=https\nhost=github.com\n\n',
        capture_output=True,text=True,timeout=60,env={**os.environ,'GIT_TERMINAL_PROMPT':'0','GCM_INTERACTIVE':'never'})
    if p.returncode:raise RuntimeError('Credential lookup failed')
    c=dict(line.split('=',1) for line in p.stdout.splitlines() if '=' in line)
    headers={'Authorization':'Bearer '+c['password'],'Accept':'application/vnd.github+json'}
    def call(method,path,body=None,allowed=(200,201)):
        for attempt in range(5):
            try:
                r=requests.request(method,BASE+path,headers=headers,json=body,timeout=(20,90))
                if r.status_code in allowed:return r
                transient=False
                if r.status_code in (400,422):
                    try:
                        transient=any(s in r.json().get('message','').lower() for s in ('malformed request','request timed out'))
                    except ValueError:pass
                if r.status_code not in (429,500,502,503,504) and not transient:
                    try:
                        detail=r.json()
                        print(json.dumps(dict(api_error_message=detail.get('message'),errors=detail.get('errors'))),flush=True)
                    except ValueError:pass
                    raise RuntimeError('GitHub API operation failed: '+method+' '+path+' HTTP '+str(r.status_code))
                print(json.dumps(dict(network_retry=attempt+1,operation=path,http_status=r.status_code)),flush=True)
            except requests.RequestException as error:
                print(json.dumps(dict(network_retry=attempt+1,operation=path,error_type=type(error).__name__)),flush=True)
            time.sleep(2**attempt)
        raise RuntimeError('GitHub API request timed out: '+path)
    if metadata_only:
        tag=call('GET','/git/ref/tags/'+version).json()['object']['sha']
        assert tag==json.loads((ROOT/'reports/public_archive_validation.json').read_text())['commit']
    else:
        assert call('GET','/git/ref/tags/'+version,allowed=(404,)).status_code==404
        if not git('branch','--list','unpublished-smart-http-v1.2'):
            git('branch','unpublished-smart-http-v1.2','HEAD')
        if git('tag','--list',version):git('tag','-d',version)
    ref=call('GET','/git/ref/heads/main').json();parent=ref['object']['sha']
    old_tree=call('GET','/git/trees/'+parent+'?recursive=1').json();assert not old_tree['truncated']
    old={r['path']:r for r in old_tree['tree'] if r['type']=='blob'}
    local={}
    for name in git('ls-files').splitlines():
        data=(CHECKOUT/name).read_bytes();sha=hashlib.sha1(b'blob '+str(len(data)).encode()+b'\0'+data).hexdigest()
        local[name]=sha
    pending=[(name,sha) for name,sha in local.items() if name not in old or old[name]['sha']!=sha]
    done_file=ROOT/'reports/api_uploaded_blobs.json';done=json.loads(done_file.read_text()) if done_file.exists() else {}
    def upload(item):
        name,sha=item
        if done.get(name)==sha:return name,sha
        # A successful server upload may outlive a lost client response.
        with requests.get(BASE+'/git/blobs/'+sha,headers=headers,stream=True,timeout=(20,45)) as check:
            if check.status_code==200:return name,sha
            if check.status_code!=404:raise RuntimeError('Blob lookup HTTP '+str(check.status_code))
        response=call('POST','/git/blobs',dict(content=base64.b64encode((CHECKOUT/name).read_bytes()).decode(),encoding='base64'))
        assert response.json()['sha']==sha
        return name,sha
    print(json.dumps(dict(archive_blobs_to_upload=len(pending),previously_uploaded=len(done))),flush=True)
    with concurrent.futures.ThreadPoolExecutor(max_workers=8) as pool:
        for i,future in enumerate(concurrent.futures.as_completed([pool.submit(upload,item) for item in pending]),1):
            name,sha=future.result();done[name]=sha
            done_file.write_text(json.dumps(done,indent=2))
            if i%25==0 or i==len(pending):print(json.dumps(dict(blobs_ready=i,total=len(pending))),flush=True)
    entries=[dict(path=name,mode='100644',type='blob',sha=sha) for name,sha in pending]
    entries += [dict(path=name,mode=old[name]['mode'],type='blob',sha=None) for name in old if name not in local]
    # Large recursive tree construction can time out even after blob uploads.
    # Keep each transaction small, using the preceding tree as its base.
    tree=old_tree['sha']
    for start in range(0,len(entries),60):
        tree=call('POST','/git/trees',dict(base_tree=tree,tree=entries[start:start+60])).json()['sha']
        print(json.dumps(dict(tree_entries_ready=min(start+60,len(entries)),total=len(entries))),flush=True)
    assert tree==git('rev-parse','HEAD^{tree}'),'API tree differs from checked archive'
    message='Attach verified Zenodo DOI and final submission metadata' if metadata_only else 'Validate nuisance initial-state theory, heterogeneous records and EJS manuscript'
    commit=call('POST','/git/commits',dict(message=message,
        tree=tree,parents=[parent])).json()['sha']
    # The update is a fast-forward from the freshly read public main.
    call('PATCH','/git/refs/heads/main',dict(sha=commit,force=False))
    if metadata_only:
        assert call('GET','/git/ref/tags/'+version).json()['object']['sha']==tag
    else:
        assert call('GET','/git/ref/tags/'+version,allowed=(404,)).status_code==404
        call('POST','/git/refs',dict(ref='refs/tags/'+version,sha=commit))
    # Preserve the unpublished local history before aligning the derived checkout.
    git('fetch','origin','main','--tags')
    git('reset','--soft',commit)
    result=dict(published=True,version=version,commit=commit,tree=tree,blobs_uploaded=len(pending),transport='GitHub Git data API')
    result['immutable_archive_tag_unchanged']=metadata_only
    name='public_metadata_update.json' if metadata_only else 'api_archive_publication.json'
    (ROOT/'reports'/name).write_text(json.dumps(result,indent=2));print(json.dumps(result))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--metadata',action='store_true');main(parser.parse_args().metadata)
