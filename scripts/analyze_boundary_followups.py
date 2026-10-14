"""Complete-sample generalization and calibrated paired Wald/MSTLS summaries."""
import argparse,json,sys
from pathlib import Path
import numpy as np
from scipy.stats import norm
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from scripts.analyze_support_extra_term import wilson


def verified(stage):
    from scripts.collect_stopped_runs import snapshot
    folder=ROOT/f'results/kaggle_{stage}';m=json.loads((folder/'manifest.json').read_text())
    audit=snapshot(stage,m);c=json.loads((folder/'completion.json').read_text());assert c['status']=='complete'
    return folder,m,audit


def general(system):
    stage='endpoint_'+system;folder,m,audit=verified(stage);p=m['protocol'];settings=p[system]
    gates=json.loads((folder/'quadrature_gates.json').read_text());assert all(g['relative_bias']<=p['clean_mean_quadrature_bias_max'] for g in gates)
    rows=json.loads((folder/'rows.json').read_text());methods=sorted({r['method'] for r in rows});n=settings['paths']
    indexed={method:{r['seed']:r for r in rows if r['method']==method} for method in methods}
    assert len(rows)==n*len(methods) and all(sorted(indexed[method])==list(range(n)) for method in methods)
    assert all(len({indexed[method][seed]['input_sha256'] for method in methods})==1 for seed in range(n))
    truth=np.array(settings.get('truth',[1.,-1.,1.,-1.]));scale=np.linalg.norm(truth)
    boot_indices=np.random.default_rng(p['bootstrap_seed']).integers(0,n,(p['bootstrap_draws'],n));groups=[];boot={};rmses={}
    for method in methods:
        errors=np.array([indexed[method][i]['coefficients'] for i in range(n)])-truth;assert np.isfinite(errors).all()
        loss=np.sum(errors**2,axis=1)/scale**2;rmses[method]=float(np.sqrt(loss.mean()));boot[method]=np.sqrt(loss[boot_indices].mean(axis=1))
        groups.append({'method':method,'paths':n,'relative_rmse':rmses[method],'relative_bias':float(np.linalg.norm(errors.mean(axis=0))/scale),
                       'relative_dispersion':float(np.sqrt(np.mean(np.sum((errors-errors.mean(axis=0))**2,axis=1)))/scale),
                       'rmse_ci95':np.quantile(boot[method],[.025,.975]).tolist()})
    gap=rmses['matched_boundary_free']-rmses['known_initial_mle'];bgap=boot['matched_boundary_free']-boot['known_initial_mle'];valid=bgap>0
    fraction=(rmses['matched_boundary_free']-rmses['universal_boundary'])/gap if gap>0 else None
    bf=(boot['matched_boundary_free'][valid]-boot['universal_boundary'][valid])/bgap[valid]
    report={'completed':True,'system':system,'protocol':p,'audit':audit,'quadrature_gates':gates,'groups':groups,
            'fraction_gap_closed':fraction,'fraction_ci95':np.quantile(bf,[.025,.975]).tolist() if len(bf) else None,
            'pass_frozen_point_rule':fraction>=.5 if fraction is not None else None,
            'oracle_geometry':json.loads((folder/'oracle_geometry.json').read_text()),
            'limitations':['Common known-initial experiment; unknown-initial MLE is a contextual reference.',
                           'First-order/full-Hermite geometry is a fixed-clean-design surrogate, not actual estimator variance decomposition.']}
    (ROOT/f'reports/boundary_generalization_{system}_results.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    print(json.dumps({'system':system,'groups':groups,'fraction_gap_closed':fraction,'fraction_ci95':report['fraction_ci95']},indent=2))


def wald():
    p=json.loads((ROOT/'configs/support_wald_v2.json').read_text());rows=[];audits=[];gates=[]
    for i in range(4):
        folder,m,audit=verified(f'support_wald_v2_{i}');assert m['protocol']==p
        audits.append(audit);g=json.loads((folder/'quadrature_gates.json').read_text());assert len(g)==8
        assert all(v['candidate_absolute_bias']<=p['clean_mean_quadrature_coefficient_absolute_error_max'] for v in g)
        gates.append(g);rows.extend(json.loads((folder/'rows.json').read_text()))
    assert len(rows)==3200
    refs=json.loads((ROOT/'configs/support_wald_references.json').read_text());groups=[];calibration=[]
    indices=np.random.default_rng(p['bootstrap_seed']).integers(0,100,(p['bootstrap_draws'],100))
    for level in [.005,.05]:
        local=next(r for r in refs['information_references'] if r['level']==level and r['label']=='full_library_boundary_free')
        for di,multiplier in enumerate(p['multipliers']):
            vectors={}
            for method in ['weak_gls_wald','author_mstls']:
                sample=sorted([r for r in rows if r['level']==level and r['delta_index']==di and r['method']==method],key=lambda r:r['seed'])
                assert [r['seed'] for r in sample]==list(range(100))
                assert all(r['multiplier']==multiplier for r in sample)
                vectors[method]=np.array([r['included'] for r in sample],dtype=float)
                count=int(vectors[method].sum());ci=wilson(count,100)
                g={'level':level,'delta_index':di,'multiplier':multiplier,'delta':sample[0]['delta'],'method':method,
                   'included':count,'paths':100,'rate':count/100,'wilson_ci95':ci}
                if method=='weak_gls_wald':
                    error=np.array([r['candidate_coefficient']-r['delta'] for r in sample]);g.update(candidate_bias=float(error.mean()),candidate_rmse=float(np.sqrt(np.mean(error**2))),median_se=float(np.median([r['standard_error'] for r in sample])))
                groups.append(g)
                if di==0:calibration.append({'level':level,'method':method,'nominal_alpha':p['alpha'] if method=='weak_gls_wald' else None,
                    'null_rate':count/100,'wilson_ci95':ci,'calibration_gate_pass':ci[0]<=p['alpha']<=ci[1] if method=='weak_gls_wald' else None})
            for seed in range(100):
                paired=[r for r in rows if r['level']==level and r['delta_index']==di and r['seed']==seed]
                assert len(paired)==2 and len({r['input_sha256'] for r in paired})==1
            difference=vectors['weak_gls_wald']-vectors['author_mstls'];boot=difference[indices].mean(axis=1)
            delta=abs(multiplier)*refs['delta_references'][str(level)];ncp=delta*np.sqrt(local['information']);critical=norm.ppf(1-p['alpha']/2)
            groups.append({'level':level,'delta_index':di,'method':'paired_difference','wald_minus_mstls':float(difference.mean()),
                           'difference_ci95':np.quantile(boot,[.025,.975]).tolist(),
                           'null_local_gaussian_wald_power_reference':float(norm.cdf(-critical-ncp)+norm.sf(critical-ncp))})
    report={'completed':True,'protocol':p,'audits':audits,'quadrature_gates':gates,'calibration':calibration,'groups':groups,
            'limitations':['A Wald test at alpha=.05 and the parent equal-prior two-point error reference are distinct criteria.',
                           'A failed null-calibration gate precludes claiming calibrated power efficiency.',
                           'Exact residual Hermite covariance at a pilot is not exact finite-noise variance of a random-design estimator.',
                           'Finite-amplitude and full-library sparse-selection outcomes do not identify variance-loss ratios.']}
    (ROOT/'reports/support_wald_results.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    print(json.dumps({'calibration':calibration,'groups':groups},indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('mode',choices=['vdp5','lorenz','wald']);args=parser.parse_args()
    wald() if args.mode=='wald' else general(args.mode)
