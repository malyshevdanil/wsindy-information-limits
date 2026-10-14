"""Read-only collection and fingerprint/paired-input audit of the direct screen."""
import argparse, json, sys, zipfile
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT))


def collect(stage):
    from kaggle.api.kaggle_api_extended import KaggleApi
    from scripts.collect_stopped_runs import redact, snapshot
    api = KaggleApi(); api.authenticate()
    ref = 'malyshevdanil/wsindy-sampling-saturation-'+stage.replace('_', '-')
    response = api.kernels_status(ref)
    status = getattr(response.status, 'name', str(response.status))
    if status not in ('COMPLETE', 'ERROR'):
        print(stage+': '+status+'; no submission performed'); return
    dest = ROOT/f'results/kaggle_{stage}'
    dest.mkdir(exist_ok=True)
    if status == 'COMPLETE':
        from scripts.download_kaggle_archive import download
        download(stage)
        with zipfile.ZipFile(ROOT/f'results/kaggle_{stage}_download/wsindy_{stage}_outputs.zip') as z:
            assert z.testzip() is None
            assert all((dest/n).resolve().is_relative_to(dest.resolve()) for n in z.namelist())
            z.extractall(dest)
        folder = dest
    else:
        api.kernels_output(ref, str(dest), file_pattern='(^|/)results_'+stage+'/', quiet=True)
        folder = dest/('results_'+stage)
    for log in dest.glob('*.log'): log.write_text(redact(log.read_text(encoding='utf-8')), encoding='utf-8')
    manifest = json.loads((folder/'manifest.json').read_text())
    verified = snapshot(stage, manifest)
    completion = json.loads((folder/'completion.json').read_text())
    gate = json.loads((folder/'quadrature_gate.json').read_text())
    report = {'stage': stage, 'api_status': status, 'snapshot': verified,
              'completion': completion, 'quadrature_gate': gate,
              'full_2000_path_experiment_completed': False}
    if completion['status'] == 'complete':
        p = manifest['protocol']; rows = json.loads((folder/'rows.json').read_text())
        old = {}
        for shard in range(4):
            for r in json.loads((ROOT/f'results/kaggle_efficiency_power_{shard}/rows.json').read_text()):
                if r['weight_model'] == 'ou': old[r['seed']] = r
        matched = all(r['input_sha256'] == old[r['seed']]['input_sha256'] for r in rows)
        report['paired_inputs_verified'] = matched
        if not matched:
            report['statistical_promotion'] = False
            report['reason'] = 'Archived paired input fingerprint mismatch; no coefficient-error summaries promoted.'
            # Preserve the numerical gate and execution audit, but do not
            # inspect/promote the outcome-dependent completion summary.
            report['completion'] = {'status': completion['status'], 'paths': completion['paths']}
            (ROOT/f'reports/{stage}_validation.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
            print(stage+': input replay mismatch; outcome summaries withheld'); return
        groups = []
        for method in ('feasible_boundary_weak', 'oracle_boundary_weak'):
            subset = [r for r in rows if r['method'] == method]
            assert [r['seed'] for r in subset] == list(range(*[p['screen_seed_indices'][0], p['screen_seed_indices'][1]+1]))
            assert all(r['input_sha256'] == old[r['seed']]['input_sha256'] for r in subset)
            err = np.array([r['coefficients'] for r in subset])-np.array([1., -1., 1., -1.])
            assert np.isfinite(err).all()
            groups.append({'method': method, 'paths': len(subset),
                           'relative_rmse': float(np.sqrt(np.mean(np.sum(err*err, axis=1)))/2),
                           'relative_bias': float(np.linalg.norm(err.mean(axis=0))/2),
                           'relative_dispersion': float(np.sqrt(np.mean(np.sum((err-err.mean(axis=0))**2, axis=1)))/2)})
        report.update(groups=groups, paired_inputs_verified=True,
                      caveat='32-path screen; threshold result is provisional and not a 2000-path validation.')
    path = ROOT/f'reports/{stage}_validation.json'
    path.write_text(json.dumps(report, indent=2), encoding='utf-8')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('stage', choices=['boundary_screen', 'boundary_screen_v2', 'boundary_screen_v3', 'boundary_screen_v4'])
    collect(parser.parse_args().stage)
