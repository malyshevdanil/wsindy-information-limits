"""Verify an automatically registered Zenodo DOI through public endpoints."""
import json,re,subprocess
from pathlib import Path
import requests
ROOT=Path(__file__).resolve().parents[1]


def main():
    cfg=json.loads((ROOT/'configs/public_archive.json').read_text())
    # Windows requests connections have repeatedly stalled on this host;
    # IPv4 curl retrieved and verified the preceding public record.
    curl=['curl.exe','-4','--http1.1','--connect-timeout','15','--max-time','35','-sS']
    r=subprocess.run(curl+['-I','https://zenodo.org/badge/latestdoi/1414450592'],capture_output=True,text=True,timeout=40)
    if r.returncode:
        print(json.dumps(dict(doi_verified=False,transport_error='curl',exit_code=r.returncode)));return
    locations=re.findall(r'^location:\s*(.*?)\s*$',r.stdout,re.I|re.M)
    location=locations[-1] if locations else '';match=re.search(r'(?:doi.org/)?(10\.5281/zenodo\.([0-9]+))',location)
    if not match:
        print(json.dumps(dict(doi_verified=False,badge_redirect_missing=True)));return
    doi,record_id=match.group(1),match.group(2)
    response=subprocess.run(curl+['--fail','https://zenodo.org/api/records/'+record_id],capture_output=True,text=True,encoding='utf-8',timeout=40)
    if response.returncode:
        print(json.dumps(dict(doi_verified=False,transport_error='curl',exit_code=response.returncode)));return
    record=json.loads(response.stdout);metadata=record['metadata']
    assert record['doi']==doi
    if metadata.get('version')!=cfg['version']:
        print(json.dumps(dict(doi_verified=False,current_version_pending=True,
            observed_version=metadata.get('version'),expected_version=cfg['version'])));return
    assert any(c.get('orcid')=='0009-0007-0947-3866' for c in metadata['creators'])
    assert metadata['title']==json.loads((ROOT/'.zenodo.json').read_text())['title']
    result=dict(doi_verified=True,doi=doi,url='https://doi.org/'+doi,record_url='https://zenodo.org/records/'+record_id,
        version=metadata['version'],creator_orcid_verified=True,public_metadata_verified=True)
    (ROOT/'reports/zenodo_current_public_record.json').write_text(response.stdout,encoding='utf-8')
    (ROOT/'reports/zenodo_archive_validation.json').write_text(json.dumps(result,indent=2))
    print(json.dumps(result))


if __name__=='__main__':main()
