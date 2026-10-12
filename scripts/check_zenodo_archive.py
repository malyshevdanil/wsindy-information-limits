"""Verify an automatically registered Zenodo DOI through public endpoints."""
import json,re
from pathlib import Path
import requests
ROOT=Path(__file__).resolve().parents[1]


def main():
    cfg=json.loads((ROOT/'configs/public_archive.json').read_text())
    r=requests.get('https://zenodo.org/badge/latestdoi/1414450592',allow_redirects=False,timeout=(20,45))
    location=r.headers.get('Location','');match=re.search(r'(?:doi.org/)?(10\.5281/zenodo\.([0-9]+))',location)
    if not match:
        print(json.dumps(dict(doi_verified=False,badge_http_status=r.status_code)));return
    doi,record_id=match.group(1),match.group(2)
    response=requests.get('https://zenodo.org/api/records/'+record_id,timeout=(20,45))
    if response.status_code!=200:
        print(json.dumps(dict(doi_verified=False,record_http_status=response.status_code)));return
    record=response.json();metadata=record['metadata']
    assert record['doi']==doi and metadata['version']==cfg['version']
    assert any(c.get('orcid')=='0009-0007-0947-3866' for c in metadata['creators'])
    assert 'nuisance initial states' in metadata['title']
    result=dict(doi_verified=True,doi=doi,url='https://doi.org/'+doi,record_url='https://zenodo.org/records/'+record_id,
        version=metadata['version'],creator_orcid_verified=True,public_metadata_verified=True)
    (ROOT/'reports/zenodo_archive_validation.json').write_text(json.dumps(result,indent=2))
    print(json.dumps(result))


if __name__=='__main__':main()
