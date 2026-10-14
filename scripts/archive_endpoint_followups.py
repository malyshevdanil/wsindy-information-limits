"""Preserve complete audit/paired outputs byte-for-byte, without raw fit scratch files."""
import hashlib,json,shutil
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
DEST=ROOT/'data/endpoint_followups'
STAGES=[f'boundary_ablate_v2_{i}' for i in range(4)]+['endpoint_vdp5','endpoint_lorenz']+[f'support_wald_v2_{i}' for i in range(4)]


def main():
    DEST.mkdir(parents=True,exist_ok=True);files=[]
    for stage in STAGES:
        source=ROOT/f'results/kaggle_{stage}'
        assert json.loads((source/'completion.json').read_text())['status']=='complete'
        target=DEST/stage;target.mkdir(exist_ok=True)
        for name in ['manifest.json','completion.json','rows.json','quadrature_gates.json','basis.json','oracle_geometry.json','information.json','author_provenance.json']:
            path=source/name
            if not path.exists():continue
            shutil.copyfile(path,target/name)
            assert path.read_bytes()==(target/name).read_bytes()
            files.append({'path':str((target/name).relative_to(DEST)).replace('\\','/'),'bytes':path.stat().st_size,'sha256':hashlib.sha256(path.read_bytes()).hexdigest()})
    quarantine=DEST/'numerical_amendment_v1';quarantine.mkdir(exist_ok=True)
    diagnostics={
        'high_noise_preflight_gates.json':ROOT/'results/failed_support_wald_2/results_support_wald_2/quadrature_gates.json',
        'support_wald_2.log':ROOT/'reports/progress_support_wald_2.log',
        'support_wald_3.log':ROOT/'reports/progress_support_wald_3.log'}
    for name,path in diagnostics.items():
        assert path.exists();shutil.copyfile(path,quarantine/name)
        data=path.read_bytes();files.append({'path':str((quarantine/name).relative_to(DEST)).replace('\\','/'),'bytes':len(data),'sha256':hashlib.sha256(data).hexdigest()})
    (DEST/'manifest.json').write_text(json.dumps({'scope':'complete endpoint ablation/generalization and primary amended Wald study; all paired records retained; failed v1 preflight diagnostics quarantined',
        'stages':STAGES,'files':files,'frozen_notebooks':'../../kaggle/<stage>/; embedded source is authoritative for each run',
        'observations':'regenerate using frozen configurations; SHA256 observation hashes preserved in rows',
        'previous_archive':'../endpoint_confirmation/; preceding benchmark is separate'},indent=2)+'\n',encoding='utf-8')
    (DEST/'README.md').write_text('# Endpoint and targeted-selection follow-ups\n\n'
        'Complete paired outputs for 500 ablation records, 128 moderate-noise Van der Pol records, 256 Lorenz records and 1600 paired Wald/MSTLS outcomes. '
        'The original 1000 MSTLS outcomes are reused; 600 refined-grid outcomes are new. All scientific configurations and embedded run source are frozen in the corresponding Kaggle notebooks.\n\n'
        'The primary Wald version corrects a deterministic pilot-model quadrature defect before inspection of the original noisy Wald results. '
        'Original versions are retained separately. Counts, null-calibration failure at 5% and failed moderate-noise endpoint recovery criteria are preserved. '
        'Author-fit scratch files and synthetic time-series arrays are omitted; seeds, source, estimates and paired SHA256 hashes suffice for reconstruction. '
        'The SHA256 manifest hashes original file bytes, with no newline conversion.\n\n'
        'For local postprocessing, copy each complete stage directory to `results/kaggle_<stage>/`, '
        'run `scripts/analyze_boundary_ablation.py`, `scripts/analyze_boundary_followups.py vdp5`, '
        '`scripts/analyze_boundary_followups.py lorenz`, `scripts/analyze_boundary_followups.py wald`, '
        'then `scripts/report_endpoint_followups.py`. Those scripts require the complete paired samples and verify the frozen notebook source. '
        '`scripts/validate_endpoint_followups.py` checks compact output hashes, notebook packages and publication files without simulations.\n',encoding='utf-8')
    print(json.dumps({'files':len(files),'bytes':sum(f['bytes'] for f in files),'stages':len(STAGES)}))


if __name__=='__main__':main()
