"""Publish the authorized GitHub archive release using existing GCM credentials.

Credentials remain in process memory and are never included in reports or logs.
"""
import json,os,subprocess
from pathlib import Path
import requests

ROOT=Path(__file__).resolve().parents[1]


def main():
    cfg=json.loads((ROOT/'configs/public_archive.json').read_text())
    verified=json.loads((ROOT/'reports/public_archive_validation.json').read_text())
    assert verified['tag']==cfg['version'] and verified['all_public_blob_hashes_verified']
    owner='malyshevdanil';repo='wsindy-information-limits';base=f'https://api.github.com/repos/{owner}/{repo}'
    process=subprocess.run(['git','credential','fill'],input='protocol=https\nhost=github.com\n\n',
        capture_output=True,text=True,timeout=60,env={**os.environ,'GIT_TERMINAL_PROMPT':'0','GCM_INTERACTIVE':'never'})
    if process.returncode:raise RuntimeError('Existing GitHub credential lookup failed')
    credential=dict(line.split('=',1) for line in process.stdout.splitlines() if '=' in line)
    session=requests.Session();session.headers.update({'Authorization':'Bearer '+credential['password'],
        'Accept':'application/vnd.github+json','X-GitHub-Api-Version':'2022-11-28'})
    user=session.get('https://api.github.com/user',timeout=30)
    assert user.status_code==200 and user.json()['login']==owner
    response=session.get(base+'/releases/tags/'+cfg['version'],timeout=30)
    if response.status_code==404:
        notes='''Theory and reproducible benchmarks for boundary-free weak ODE estimation with nuisance initial states.

Includes exact weak-class and local-information results; full polynomial covariance; independent-instrument heterogeneous-record simulations with 120/480 tests, 400 repetitions and record-bootstrap inference; matched 128-path likelihood at 5% noise; author-method comparisons and all frozen Kaggle execution snapshots.

The heterogeneous-record benchmark uses known OU covariance and conditions on a separate frozen training sample. At 5% noise its feasible 480-test scaled RMSE remains above the local reference; at 0.5% and 800 records joint coverage is 92.25%. The archive preserves these limitations and every planned result. It makes no claim of journal acceptance or submission.

The manuscript uses the official IMS EJS class. Creator ORCID: https://orcid.org/0009-0007-0947-3866. Code and data are accompanied by SHA256 manifests and verified replay scripts.'''
        response=session.post(base+'/releases',json=dict(tag_name=cfg['version'],target_commitish=verified['commit'],
            name=cfg['version']+' - nuisance initial states and record-level inference',body=notes,
            draft=False,prerelease=False),timeout=60)
    if response.status_code not in (200,201):
        raise RuntimeError('GitHub release request failed with HTTP '+str(response.status_code))
    data=response.json();assert data['tag_name']==cfg['version'] and not data['draft']
    result=dict(published=True,tag=data['tag_name'],url=data['html_url'],release_id=data['id'],
        created_at=data['created_at'],published_at=data['published_at'],zenodo_doi_pending=True)
    (ROOT/'reports/github_release_validation.json').write_text(json.dumps(result,indent=2))
    print(json.dumps(result))


if __name__=='__main__':main()
