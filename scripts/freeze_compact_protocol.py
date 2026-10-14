"""Freeze the user's final 120-window budget on fresh physical seed ranges."""
import hashlib,json
from datetime import datetime,timezone
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]


def freeze(name,value):
    path=ROOT/'configs'/name
    assert not path.exists(), 'Frozen configurations must never be overwritten'
    path.write_text(json.dumps(value,indent=2))
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    p=json.loads((ROOT/'configs/incidental_main.json').read_text())
    p.update(protocol='compact_incidental_comparison_v3',previous_protocol='incidental_trimmed.json',
        levels=[.05,.1],record_counts=[100,400,800],monte_carlo_repetitions=200,
        seed_base=7100000,parametric_bootstrap_draws=0,bootstrap_correction_cells=[],
        methods=['profile_likelihood','oracle_likelihood','crossfit_fission_gmm'],weak_windows=120,
        weak_families=['repairAa','repairAb'],status='frozen fresh-seed compact comparison; no MLE bootstrap correction',
        amendment='User-directed final budget: 120 windows, 5/10% noise, n=100/400/800, 200 repetitions')
    hashes={'incidental_compact.json':freeze('incidental_compact.json',p)}
    base=json.loads((ROOT/'configs/repair120_records.json').read_text())
    for prefix,levels,counts,seed,extra,oracle in [('repairBa',[.005,.05],[100,800],7200000,7400000,True),('repairAa',[.05],[100,400,800],7100000,7300000,False),('repairAb',[.1],[100,400,800],7100000,7300000,False)]:
        config={**base,'protocol':'compact_120_'+prefix,'previous_protocol':'repair120_records.json',
            'levels':levels,'record_counts':counts,'monte_carlo_repetitions':200,
            'seed_base':seed,'extra_record_seed_base':extra,'include_oracle_fission_comparison':oracle,
            'include_split_diagnostic':False,'compare_frozen_bootstrap':False,
            'bootstrap_seed_base':820000000 if prefix=='repairBa' else 810000000,
            'bootstrap_amplitude_index':0 if prefix!='repairAb' else 1,
            'cache_kernel_source':'malyshevdanil/wsindy-sampling-saturation-'+prefix.lower()+'-geometry',
            'instrument_cache':'degree-seven Legendre map, 1024 fixed uniform Sobol nodes; 5% trace-variance gate',
            'smoke_record_counts':[100,800],'matched_incidental_main':prefix!='repairBa',
            'budget_comparison':'200 repetitions, 120 tests, fivefold cross-fitting, no external training sample or split diagnostic'}
        if prefix=='repairAa':
            config.update(geometry_source_stage='repairBa_geometry',geometry_protocol_file='repairBa_records.json',
                cache_kernel_source='malyshevdanil/wsindy-sampling-saturation-repairba-geometry')
        hashes[prefix+'_records.json']=freeze(prefix+'_records.json',config)
    probe=dict(families=[120],levels=[.1],initial=[2.25,.25],threshold=.05)
    hashes['cache120_main_probe.json']=freeze('cache120_main_probe.json',probe)
    decision=dict(decision_time_utc=datetime.now(timezone.utc).isoformat(),user_directed=True,
        previous_protocols_preserved=True,main_numerical_outcomes_used_to_select_scope=False,
        pilot_and_deterministic_checks_used=True,protocol_sha256=hashes,
        physical_seeds_A=[7100000,7100199],physical_seeds_B=[7200000,7200199],
        planned_mle_bootstrap_refits=0,planned_profile_fits=1200,planned_oracle_fits=1200,
        replicas=200,weak_windows=120,known_covariance=True,
        reuse='A 5%-noise instruments reuse the fully validated B cache; physical A and B data have disjoint seeds',
        excluded_scope=['20%-noise main comparison','480-window new MC','MLE bootstrap correction','external 1000-record training diagnostic'],
        publication_rule='No empirical claim of absence of a crossing until matched comparison is complete; analytic crossing is labelled a local order-of-magnitude forecast')
    (ROOT/'reports/compact_protocol_decision.json').write_text(json.dumps(decision,indent=2))
    print('Compact 120-window protocols frozen on fresh physical seeds')


if __name__=='__main__':main()
