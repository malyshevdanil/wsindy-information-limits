"""Freeze a higher-order cache revision; retain seeds and the 5% criterion."""
import json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]


def freeze(path,value):
    if path.exists():assert json.loads(path.read_text())==value
    else:path.write_text(json.dumps(value,indent=2))


def main():
    c=json.loads((ROOT/'configs/crossfit_cache.json').read_text())
    c.update(degree=7,node_power=10,status='repair probe required before geometry',
        revision_reason='Degree-five validation failed at axis probes; same support and 5% loss threshold, more polynomial terms and design nodes')
    freeze(ROOT/'configs/crossfit_cache_v3.json',c)
    for family,original in (('repair120','crossfit120_records.json'),('repair480','crossfit_records.json')):
        p=json.loads((ROOT/'configs'/original).read_text())
        p.update(protocol=p['protocol']+'_cache_v3',cache_kernel_source='malyshevdanil/wsindy-sampling-saturation-'+family+'-geometry',
            previous_protocol=original,cache_protocol_file='crossfit_cache_v3.json')
        freeze(ROOT/'configs'/(family+'_records.json'),p)
    freeze(ROOT/'configs/cache_repair_probe.json',dict(protocol='degree_seven_cache_repair_probe_v1',
        families=[120,480],levels=[.005,.05],initial=[2.25,.25],threshold=.05,
        validation='zero, four prescribed offsets, and all eight axis endpoints; unchanged parameter support'))
    freeze(ROOT/'configs/cache_repair_probe_A.json',dict(protocol='degree_seven_cache_repair_large_noise_probe_v1',
        families=[480],levels=[.1,.2],initial=[2.25,.25],threshold=.05,
        validation='zero, four prescribed offsets, and all eight axis endpoints; 5% amplitude covered by shared repair probe'))
    original=json.loads((ROOT/'configs/incidental_weak.json').read_text())
    for index,level in enumerate(original['levels']):
        family=f'repairA{index}';p={**original}
        p.update(protocol=original['protocol']+'_cache_v3_level_'+str(index),levels=[level],
            cache_kernel_source='malyshevdanil/wsindy-sampling-saturation-'+family.lower()+'-geometry',
            previous_protocol='incidental_weak.json',cache_protocol_file='crossfit_cache_v3.json',
            bootstrap_amplitude_index=index,bootstrap_seed_encoding='base + repetition*1000000 + amplitude_index*100000 + n')
        freeze(ROOT/'configs'/(family+'_records.json'),p)


if __name__=='__main__':main()
