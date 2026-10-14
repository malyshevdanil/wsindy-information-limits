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
        notes='''The information cost of boundary-free weak estimation under correlated noise.

Includes exact adjoint-range characterizations and initial-state information cost; a constructive full-library small-noise efficiency theorem under a checkable jet condition; its exact Van der Pol verification; a deterministic finite-span covariance reference; and paired boundary and testing benchmarks.

Independent small-noise Wald validation retains 2000 null observations and 1000 observations for each of two alternatives. The size tolerance passes; the full power tolerance check fails at half a reference unit. All five guarded-estimator noise levels and their matched covariance reference are retained. Moderate-noise bootstrap validation is incomplete following a recorded compute-budget decision; unsubmitted frozen notebooks are protocol artifacts, not completed results.

The main manuscript uses the official IMS EJS class. Core proofs are in a shortened submission supplement; independent-record, fission, cross-fitting and incidental-state extensions are preserved separately. Source and original-byte outputs carry SHA-256 manifests. This archive makes no claim of journal submission or acceptance. Creator ORCID: https://orcid.org/0009-0007-0947-3866.'''
        response=session.post(base+'/releases',json=dict(tag_name=cfg['version'],target_commitish=verified['commit'],
            name=cfg['version']+' - constructive efficiency and information cost',body=notes,
            draft=False,prerelease=False),timeout=60)
    if response.status_code not in (200,201):
        raise RuntimeError('GitHub release request failed with HTTP '+str(response.status_code))
    data=response.json();assert data['tag_name']==cfg['version'] and not data['draft']
    result=dict(published=True,tag=data['tag_name'],url=data['html_url'],release_id=data['id'],
        created_at=data['created_at'],published_at=data['published_at'],zenodo_doi_pending=not bool(cfg.get('doi')),doi=cfg.get('doi'))
    (ROOT/'reports/github_release_validation.json').write_text(json.dumps(result,indent=2))
    print(json.dumps(result))


if __name__=='__main__':main()
