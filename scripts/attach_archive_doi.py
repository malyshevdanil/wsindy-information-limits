"""Attach a publicly verified DOI without changing the immutable archive tag."""
import hashlib,json,re
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]


def main():
    verified=json.loads((ROOT/'reports/zenodo_archive_validation.json').read_text())
    assert verified['doi_verified'] and verified['creator_orcid_verified']
    p=ROOT/'configs/public_archive.json';cfg=json.loads(p.read_text())
    assert cfg['version']==verified['version'];cfg['doi']=verified['doi'];p.write_text(json.dumps(cfg,indent=2))
    p=ROOT/'manuscript/article_revised.md';before=p.read_text(encoding='utf-8')
    start=before.index('## Data availability\n\n')+len('## Data availability\n\n')
    stop=before.index('\n\n',start)
    paragraph=(f'The current code and synthetic-data archive, {cfg["version"]}, is publicly '
        f'available in the [versioned reproducibility repository]({cfg["version_url"]}) '
        f'and its [Zenodo snapshot]({verified["url"]}). Configurations, seed ranges, '
        'source fingerprints, coefficient estimates, covariance summaries and SHA256 '
        'manifests permit exact input reconstruction and paired audit. Full synthetic '
        'time series are regenerated from the frozen seed protocol. Pinned download '
        'specifications retrieve original author sources. The separate research archive '
        'preserves earlier extensions. The incomplete moderate-noise bootstrap study is '
        'identified as such; unsubmitted notebooks document protocols only. '
        f'The [preceding {cfg["previous_version"]} snapshot](https://doi.org/{cfg["previous_doi"]}) '
        'remains unchanged. No external experimental dataset is used.')
    after=before[:start]+paragraph+before[stop:]
    p.write_text(after,encoding='utf-8')
    p=ROOT/'CITATION.cff';text=re.sub(r'^doi:.*\n?', '', p.read_text(encoding='utf-8'),flags=re.M)
    p.write_text(text.rstrip()+'\ndoi: "'+verified['doi']+'"\n',encoding='utf-8')
    p=ROOT/'README.md';text=p.read_text(encoding='utf-8')
    marker='\nCurrent reproducibility archive: '
    text=text.split(marker)[0].rstrip()
    p.write_text(text+marker+f'[{cfg["version"]}]({cfg["version_url"]}), '
        f'[Zenodo DOI {verified["doi"]}]({verified["url"]}).\n',encoding='utf-8')
    p=ROOT/'manuscript/submission/declarations_to_complete.md'
    p.write_text(p.read_text(encoding='utf-8').replace('A DOI has not been assigned.',
        'The current archive DOI is '+verified['url']+'.'),encoding='utf-8')
    p=ROOT/'manuscript/submission/cover_letter_draft.md'
    text=p.read_text(encoding='utf-8').replace(
        'The current versioned archive identifier and factual submission declarations will be finalized before submission.',
        f'The code and data archive is {cfg["version"]}, DOI {verified["doi"]}. '
        'Factual submission declarations will be finalized before submission.')
    p.write_text(text,encoding='utf-8')
    report=dict(doi=verified['doi'],immutable_tag_unchanged=True,manuscript_change='Data availability DOI sentence only',
        manuscript_before_sha256=hashlib.sha256(before.encode()).hexdigest(),manuscript_after_sha256=hashlib.sha256(after.encode()).hexdigest())
    (ROOT/'reports/doi_metadata_attachment.json').write_text(json.dumps(report,indent=2));print(json.dumps(report))


if __name__=='__main__':main()
