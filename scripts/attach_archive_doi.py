"""Attach a publicly verified DOI without changing the immutable archive tag."""
import hashlib,json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]


def main():
    verified=json.loads((ROOT/'reports/zenodo_archive_validation.json').read_text())
    assert verified['doi_verified'] and verified['creator_orcid_verified']
    p=ROOT/'configs/public_archive.json';cfg=json.loads(p.read_text())
    assert cfg['version']==verified['version'];cfg['doi']=verified['doi'];p.write_text(json.dumps(cfg,indent=2))
    p=ROOT/'manuscript/article_revised.md';before=p.read_text(encoding='utf-8')
    after=before.replace('a DOI has not been assigned.',
        'the code and data are preserved in the [Zenodo archive]('+verified['url']+').')
    assert before.count('a DOI has not been assigned.')==1
    p.write_text(after,encoding='utf-8')
    p=ROOT/'CITATION.cff';p.write_text(p.read_text(encoding='utf-8')+'doi: "'+verified['doi']+'"\n',encoding='utf-8')
    p=ROOT/'README.md';p.write_text(p.read_text(encoding='utf-8')+'\nThe verified v1.2.0 code and data archive is preserved at [Zenodo DOI '+verified['doi']+']('+verified['url']+').\n',encoding='utf-8')
    p=ROOT/'manuscript/submission/declarations_to_complete.md'
    p.write_text(p.read_text(encoding='utf-8').replace('A DOI has not been assigned.',
        'The v1.2.0 archive DOI is '+verified['url']+'.'),encoding='utf-8')
    report=dict(doi=verified['doi'],immutable_tag_unchanged=True,manuscript_change='Data availability DOI sentence only',
        manuscript_before_sha256=hashlib.sha256(before.encode()).hexdigest(),manuscript_after_sha256=hashlib.sha256(after.encode()).hexdigest())
    (ROOT/'reports/doi_metadata_attachment.json').write_text(json.dumps(report,indent=2));print(json.dumps(report))


if __name__=='__main__':main()
