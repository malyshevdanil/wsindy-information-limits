"""Queue and audit deterministic forecasting before further main MC jobs."""
import argparse,json,sys,time,zipfile
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))


def main(bias_check=False,collect_only=False):
    hold=ROOT/'reports/kaggle_submission_hold.json'
    if not collect_only and hold.exists() and json.loads(hold.read_text()).get('enabled'):
        print('Automatic Kaggle submissions held; use --collect-only for finished outputs',flush=True);return
    stage='inc_bias_check' if bias_check else 'inc_forecast'
    config_file='incidental_bias_check.json' if bias_check else 'incidental_forecast.json'
    report=ROOT/('reports/incidental_finite_n_bias.json' if bias_check else 'reports/incidental_forecast_exact.json')
    while not report.exists():
        try:
            from kaggle.api.kaggle_api_extended import KaggleApi
            api=KaggleApi();api.authenticate()
            try:response=api.kernels_status('malyshevdanil/wsindy-sampling-saturation-'+stage.replace('_','-'))
            except Exception:
                if collect_only:raise
                api.kernels_status('malyshevdanil/wsindy-sampling-saturation-efficiency-match')
                response=api.kernels_push(str(ROOT/'kaggle'/stage))
                print(stage+' waiting for CPU slot' if getattr(response,'error',None) else stage+' submitted',flush=True)
                time.sleep(45);continue
            status=getattr(response.status,'name',str(response.status)).split('.')[-1].upper()
            if collect_only and status!='COMPLETE':
                print(stage+' is '+status+'; no submission performed',flush=True);return
            if status=='CANCEL_ACKNOWLEDGED':
                response=api.kernels_push(str(ROOT/'kaggle'/stage))
                print(stage+' waiting for CPU slot' if getattr(response,'error',None) else stage+' resubmitted frozen notebook',flush=True)
            elif status in ('ERROR','FAILED','CANCELED','CANCELLED'):
                (ROOT/f'reports/{stage}_failure.json').write_text(json.dumps(dict(stage=stage,status=status)))
                return
            elif status=='COMPLETE':
                from scripts.download_kaggle_archive import download
                download(stage);folder=ROOT/f'results/kaggle_{stage}';folder.mkdir(exist_ok=True)
                with zipfile.ZipFile(ROOT/f'results/kaggle_{stage}_download/wsindy_{stage}_outputs.zip') as z:
                    assert z.testzip() is None
                    assert all((folder/n).resolve().is_relative_to(folder.resolve()) for n in z.namelist())
                    z.extractall(folder)
                from scripts.validate_nuisance_runs import snapshot
                verified=snapshot(stage,config_file)
                config=json.loads((ROOT/'configs'/config_file).read_text())
                assert json.loads((folder/'completion.json').read_text())['status']=='complete'
                if bias_check:
                    checks=json.loads((folder/'finite_difference_checks.json').read_text())
                    assert max(c['relative_error'] for c in checks)<=config['relative_finite_difference_tolerance']
                    rows=json.loads((folder/'bias_check.json').read_text())['rows']
                    previous=json.loads((ROOT/'reports/incidental_curvature_validation.json').read_text())['results']['results'][-1]['bias_coefficient']
                    assert np.allclose(rows[-1]['incidental_coefficient'],previous,rtol=1e-5,atol=1e-8)
                    report.write_text(json.dumps(dict(completed=True,execution_snapshot=verified,protocol=config,rows=rows,finite_difference_checks=checks,
                        caveat='Second-order small-noise and first-order 1/n expansion; higher-noise remainders and local optimizer branch effects require Monte Carlo checks.'),indent=2))
                    print('Finite-n joint curvature verified',flush=True);return
                checks=json.loads((folder/'quadrature_checks.json').read_text())
                assert max(c['relative_difference'] for c in checks)<=config['quadrature_relative_tolerance']
                rows=json.loads((folder/'forecast_variance.json').read_text())['rows']
                assert len(rows)==len(config['initial_quadrature_orders'])*len(config['families'])*len(config['levels'])*2
                curvature=json.loads((ROOT/'reports/incidental_curvature_validation.json').read_text())['results']['results'][-1]
                coefficient=np.array(curvature['bias_coefficient']);I=np.array(curvature['mean_information'])
                from scripts.forecast_incidental_crossing import crossing
                groups=[]
                for row in rows:
                    V=np.array(row['covariance_of_sqrt_n_estimate']);assert np.linalg.eigvalsh(V).min()>0
                    assert np.isclose(np.trace(V)/4,row['relative_variance'])
                    if row['order']!=config['initial_quadrature_orders'][-1]:continue
                    level=row['level'];mle=level**2*np.trace(np.linalg.inv(I))/4;bias=level**4*(coefficient@coefficient)/4
                    groups.append({**row,'mle_leading_relative_variance':float(mle),'mle_bias_squared':float(bias),
                        'n_star':crossing(row['relative_variance'],mle,bias)})
                report.write_text(json.dumps(dict(completed=True,execution_snapshot=verified,protocol=config,groups=groups,quadrature_checks=checks,
                    caveat='Oracle deterministic instruments require true per-record initial states. These are a variance benchmark, not the new feasible cross-fitted variance. Crossing assumes small-noise MLE bias/variance and negligible weak finite-n bias.'),indent=2))
                import matplotlib
                matplotlib.use('Agg')
                import matplotlib.pyplot as plt
                fig,ax=plt.subplots(figsize=(7,4.5))
                for windows in config['families']:
                    for method in ('oracle_original','oracle_fission'):
                        g=[r for r in groups if r['windows']==windows and r['method']==method]
                        ax.plot([100*r['level'] for r in g],[max(1,r['n_star']) for r in g],marker='o',label=f'{windows} tests, {method}')
                ax.set_yscale('log');ax.axhline(800,color='black',linestyle=':',label='Largest planned n')
                ax.set(xlabel='Noise amplitude (%)',ylabel='Leading-order crossing n*',title='Oracle weak variance benchmark; feasible weights pending')
                ax.grid(alpha=.25);ax.legend(fontsize=8);fig.tight_layout()
                for extension in ('png','pdf'):fig.savefig(ROOT/f'reports/figures/incidental_crossing_exact.{extension}',dpi=200)
                plt.close(fig);print('Deterministic forecast verified; submission hold is managed separately',flush=True);return
        except (AssertionError,FileNotFoundError) as error:
            (ROOT/f'reports/{stage}_failure.json').write_text(json.dumps(dict(stage=stage,error_type=type(error).__name__)))
            return
        except Exception as error:print(json.dumps(dict(retry=True,error_type=type(error).__name__)),flush=True)
        time.sleep(45)


if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--bias-check',action='store_true');a.add_argument('--collect-only',action='store_true')
    args=a.parse_args();main(args.bias_check,args.collect_only)
