"""Preserve complete numerical summaries with hashes; omit redundant raw fit files."""
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def main():
    status = json.loads((ROOT/'reports/endpoint_confirmation_status.json').read_text())
    assert len(status) == 5
    assert all(r['status'] == 'COMPLETE' and r['snapshot']['source_sha256_verified']
               and r['completion']['status'] == 'complete' for r in status)
    destination = ROOT/'data/endpoint_confirmation'
    files = {}
    for record in status:
        stage = record['stage']
        source = ROOT/f'results/kaggle_{stage}'
        target = destination/stage
        target.mkdir(parents=True, exist_ok=True)
        names = ['manifest.json', 'completion.json', 'rows.json']
        names += (['basis.json', 'quadrature_gate.json'] if stage.startswith('boundary')
                  else ['information.json', 'delta_grids.json', 'author_provenance.json'])
        for name in names:
            data = (source/name).read_bytes()
            (target/name).write_bytes(data)
            files[f'{stage}/{name}'] = hashlib.sha256(data).hexdigest()
    for name in ['endpoint_confirmation_results.json', 'support_extra_term_results.json',
                 'endpoint_confirmation_status.json']:
        data = (ROOT/'reports'/name).read_bytes()
        (destination/name).write_bytes(data)
        files[name] = hashlib.sha256(data).hexdigest()
    (destination/'sha256.json').write_text(json.dumps(files, indent=2)+'\n', encoding='utf-8')
    (destination/'README.md').write_text(
        '# Complete endpoint and added-term study outputs\n\n'
        'All four endpoint shards (1968 paired records, four branches) and all 1000 '
        'complete author WSINDy fits passed the frozen source/configuration audit. '
        'The rows here preserve every coefficient result, including zero inclusions.\n\n'
        'Frozen protocols are `configs/boundary_confirmation.json` and '
        '`configs/support_extra_term.json`; execution sources are embedded in the '
        'five frozen notebooks under `kaggle/boundary_confirm_*` and `kaggle/support_extra`. '
        'Observed time series can be regenerated from those protocols. '
        'Full fit artifacts remain in the downloaded Kaggle output ZIPs; this compact '
        'archive contains the numerical results and provenance needed for analysis.\n\n'
        'To rerun local analysis, copy each stage directory to '
        '`results/kaggle_<stage>/`, then run `scripts/analyze_endpoint_confirmation.py`, '
        '`scripts/report_endpoint_confirmation.py` and `scripts/analyze_support_extra_term.py`. '
        'Those scripts verify the archived notebook source and require complete samples. '
        'Research simulations run on Kaggle.\n', encoding='utf-8')
    print(f'Archived {len(files)} complete summary files with SHA256 hashes')


if __name__ == '__main__':
    main()
