"""Freeze the weak component only after the independent amplitude pilot."""
import json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]


def prepare():
    main=json.loads((ROOT/'configs/incidental_main.json').read_text())
    p=json.loads((ROOT/'configs/crossfit_records.json').read_text())
    p.update(protocol='matched_incidental_crossfit_fission_v1',levels=main['levels'],
        record_counts=main['record_counts'],smoke_record_counts=[50,800],
        monte_carlo_repetitions=main['monte_carlo_repetitions'],seed_base=main['seed_base'],
        extra_record_seed_base=4400000,bootstrap_seed_base=main['weak_bootstrap_seed_base'],
        cache_kernel_source='malyshevdanil/wsindy-sampling-saturation-crossfita-geometry',
        matched_likelihood_protocol='incidental_main.json')
    target=ROOT/'configs/incidental_weak.json'
    if target.exists():assert json.loads(target.read_text())==p
    else:target.write_text(json.dumps(p,indent=2))


if __name__=='__main__':prepare()
