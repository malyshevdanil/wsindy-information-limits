"""Rebuild numerical summaries and manuscript text from the public ZIP pair."""
import hashlib,json,subprocess,sys,zipfile
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
def main():
    target=ROOT/'reports/release_check_nuisance_v1'
    target.mkdir(exist_ok=False)
    for name in ['reproducibility_code','reproducibility_data']:
        with zipfile.ZipFile(ROOT/f'dist/{name}.zip') as archive:
            assert archive.testzip() is None
            for entry in archive.infolist():assert (target/entry.filename).resolve().is_relative_to(target.resolve())
            archive.extractall(target)
    commands=[*[["scripts/unpack_kaggle_outputs.py",f"efficiency_match_{i}"] for i in range(4)],
        ['scripts/analyze_revision.py','--stage','revision_fast'],['scripts/analyze_full_adaptive.py'],
        ['scripts/analyze_reviewer_validation.py'],['scripts/prepare_efficiency_manuscript.py'],
        ['scripts/analyze_efficiency.py'],['scripts/analyze_efficiency_followup.py'],
        ['-m','scripts.analyze_information_followup'],['scripts/prepare_efficiency_followup_manuscript.py'],
        ['scripts/analyze_independent_records.py'],['scripts/prepare_independent_records_manuscript.py'],
        ['-m','scripts.analyze_record_bootstrap'],['scripts/prepare_nuisance_manuscript.py'],['-m','scripts.analyze_heterogeneous_records'],['scripts/prepare_heterogeneous_results.py'],['-m','scripts.analyze_moderate_likelihood'],['scripts/prepare_moderate_likelihood_results.py'],['scripts/validate_nuisance_runs.py'],['scripts/validate_independent_records.py'],
        ['scripts/validate_efficiency.py'],['scripts/validate_efficiency_followup.py']]
    for args in commands:subprocess.run([sys.executable,*args],cwd=target,check=True,stdout=subprocess.DEVNULL)
    files=['reports/efficiency_results.json','reports/efficiency_followup_results.json','reports/information_followup_results.json',
        'reports/independent_records_results.json','reports/record_bootstrap_results.json','reports/heterogeneous_records_results.json','reports/moderate_likelihood_results.json','manuscript/heterogeneous_records_theory.md','manuscript/independent_records_theory.md',
        'manuscript/article_revised.md','manuscript/supplement.md','manuscript/weak_efficiency_theory.md','reports/efficiency_followup_RU.md']
    records=[]
    for name in files:
        a=(ROOT/name).read_bytes();b=(target/name).read_bytes();assert a==b,f'Reproduction discrepancy: {name}'
        records.append(dict(path=name,sha256=hashlib.sha256(a).hexdigest(),byte_identical=True))
    result=dict(completed=True,fresh_code_and_data_zip_extraction=True,all_statistics_intervals_and_manuscript_text_identical=True,files=records)
    (ROOT/'reports/efficiency_followup_reproduction.json').write_text(json.dumps(result,indent=2),encoding='utf-8');print(json.dumps(result))
if __name__=='__main__':main()
