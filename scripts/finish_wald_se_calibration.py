"""Read-only collection and full-sample analysis of frozen SE diagnostics."""
import argparse,json,sys,time,zipfile,hashlib
from pathlib import Path
import numpy as np
from scipy.stats import norm
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from scripts.collect_stopped_runs import snapshot
from scripts.download_kaggle_archive import download


def write(path,value):path.write_text(json.dumps(value,indent=2,allow_nan=False),encoding='utf-8')


def wilson(k,n):
    z=norm.ppf(.975);p=k/n;c=(p+z*z/(2*n))/(1+z*z/n)
    h=z*np.sqrt(p*(1-p)/n+z*z/(4*n*n))/(1+z*z/n)
    return [float(c-h),float(c+h)]


def collect(stage,api):
    folder=ROOT/'results'/f'secal_{stage}'
    if not (folder/'completion.json').exists():
        download(stage,api=api)
        source=ROOT/'results'/f'kaggle_{stage}_download'/f'wsindy_{stage}_outputs.zip'
        folder.mkdir(parents=True,exist_ok=True)
        with zipfile.ZipFile(source) as archive:
            assert archive.testzip() is None
            for name in archive.namelist():
                target=(folder/name).resolve()
                assert target.is_relative_to(folder.resolve())
            archive.extractall(folder)
    manifest=json.loads((folder/'manifest.json').read_text())
    proof=snapshot(stage,manifest)
    protocol=json.loads((ROOT/'configs/wald_se_calibration.json').read_text())
    assert manifest['protocol']==protocol
    completion=json.loads((folder/'completion.json').read_text());assert completion['status']=='complete'
    expected=1 if stage=='secal_smoke' else 20
    assert completion['paths']==expected and completion['rows']==3*expected
    rows=json.loads((folder/'rows.json').read_text());assert len(rows)==3*expected
    assert len({(r['seed'],r['method']) for r in rows})==len(rows)
    refs={r['seed']:r for r in json.loads((ROOT/'configs/wald_se_references.json').read_text())['rows']}
    for r in rows:
        ref=refs[r['seed']];assert r['input_sha256']==ref['input_sha256']
        assert abs(r['candidate_coefficient']-ref['candidate_coefficient'])<=protocol['coefficient_replay_absolute_tolerance']
        if r['method']=='nominal':assert abs(r['standard_error']/ref['standard_error']-1)<=protocol['se_replay_relative_tolerance']
    proof.update(stage=stage,completion=completion)
    return rows,proof,folder


def analyze():
    allrows=[];proofs=[];diagnostics=[]
    bundle=ROOT/'data/wald_se_calibration'
    if (bundle/'manifest.json').exists():
        inventory=json.loads((bundle/'manifest.json').read_text())
        for record in inventory['files']:
            file=(bundle/record['path']).resolve()
            assert file.is_relative_to(bundle.resolve())
            assert file.stat().st_size==record['bytes']
            assert hashlib.sha256(file.read_bytes()).hexdigest()==record['sha256']
    protocol=json.loads((ROOT/'configs/wald_se_calibration.json').read_text())
    decision=json.loads((ROOT/'reports/wald_se_protocol_decision.json').read_text())
    assert hashlib.sha256((ROOT/'configs/wald_se_calibration.json').read_bytes()).hexdigest()==decision['protocol_sha256']
    assert hashlib.sha256((ROOT/'configs/wald_se_references.json').read_bytes()).hexdigest()==decision['reference_sha256']
    for i in range(5):
        stage=f'secal_{i}';folder=bundle/stage if (bundle/stage).exists() else ROOT/'results'/f'secal_{stage}'
        manifest=json.loads((folder/'manifest.json').read_text());proofs.append(snapshot(stage,manifest))
        assert manifest['protocol']==protocol
        completion=json.loads((folder/'completion.json').read_text());assert completion['status']=='complete' and completion['paths']==20 and completion['rows']==60
        rows=json.loads((folder/'rows.json').read_text());assert len(rows)==60
        for method in ['nominal','design_noise_sandwich','parametric_bootstrap']:
            assert sorted(r['seed'] for r in rows if r['method']==method)==[s for s in range(100) if s%5==i]
        allrows.extend(rows)
        diagnostics.extend(json.loads(p.read_text()) for p in sorted(folder.glob('diagnostics_s*.json')))
    assert len(allrows)==300 and len(diagnostics)==100
    methods=['nominal','design_noise_sandwich','parametric_bootstrap'];summary={}
    ref={r['seed']:r for r in json.loads((ROOT/'configs/wald_se_references.json').read_text())['rows']}
    for method in methods:
        rows=sorted((r for r in allrows if r['method']==method),key=lambda r:r['seed'])
        assert [r['seed'] for r in rows]==list(range(100))
        assert all(r['input_sha256']==ref[r['seed']]['input_sha256'] for r in rows)
        assert all(abs(r['candidate_coefficient']-ref[r['seed']]['candidate_coefficient'])<=1e-8 for r in rows)
        assert all(np.isfinite(r['standard_error']) and r['standard_error']>0 for r in rows)
        assert all(r['rejected']==(abs(r['candidate_coefficient']/r['standard_error'])>norm.ppf(.975)) for r in rows)
        if method=='nominal':assert all(abs(r['standard_error']/ref[r['seed']]['standard_error']-1)<=1e-8 for r in rows)
        k=sum(r['rejected'] for r in rows);interval=wilson(k,100)
        summary[method]={'paths':100,'rejections':k,'rejection_rate':k/100,'wilson95':interval,
            'compatible_with_5_percent':interval[0]<=.05<=interval[1],
            'median_se':float(np.median([r['standard_error'] for r in rows])),
            'median_se_ratio_to_nominal':float(np.median([r['standard_error']/ref[r['seed']]['standard_error'] for r in rows]))}
    assert summary['nominal']['rejections']==23
    assert all(len(d['bootstrap_targets'])==399 and d['score_draws']==999 for d in diagnostics)
    boot_rows={r['seed']:r for r in allrows if r['method']=='parametric_bootstrap'}
    assert sorted(d['seed'] for d in diagnostics)==list(range(100))
    assert all(abs(np.std(d['bootstrap_targets'],ddof=1)-boot_rows[d['seed']]['standard_error'])<1e-12 for d in diagnostics)
    result={'protocol':'vdp5_null_paired_wald_se_v1','complete_paths':100,'methods':summary,
        'all_source_snapshots_verified':True,'diagnostic_previously_seen_null_paths':True,
        'independent_confirmation':False,'score_seconds_total':sum(d['score_seconds'] for d in diagnostics),
        'bootstrap_seconds_total':sum(d['bootstrap_seconds'] for d in diagnostics),
        'mean_bootstrap_bias':float(np.mean([d['bootstrap_mean'] for d in diagnostics])),
        'mean_score_influence':float(np.mean([d['score_influence_mean'] for d in diagnostics]))}
    write(ROOT/'reports/wald_se_calibration_results.json',result)
    lines=['# Диагностика SE Wald при 5% шума','',
        'Все 100 исходных нулевых путей сохранены. Коэффициент и порог z0.975 одинаковы; меняется только SE. Критерий заморожен до новых расчётов: 95% Wilson-интервал должен содержать 5%.','',
        '| SE | Ложные срабатывания | Wilson 95% | Совместимость с 5% | Медиана отношения SE к номинальному |',
        '|---|---:|---|---|---:|']
    for method,s in summary.items():
        lo,hi=s['wilson95'];lines.append(f"| {method} | {s['rejections']}/100 | [{100*lo:.1f}%; {100*hi:.1f}%] | {'да' if s['compatible_with_5_percent'] else 'нет'} | {s['median_se_ratio_to_nominal']:.3f} |")
    lines+=['','Это диагностика на ранее просмотренных путях. Совместимость с 5% на 100 путях требует независимого подтверждения и сама по себе не устанавливает калибровку мощности.',
        'Сэндвич включает случайный полиномиальный дизайн, но фиксирует наблюдённые веса; его meat оценён по 999 модельным OU-репликам. Бутстреп использует 399 реплик с полным повторным построением пилота и весов. Оба используют нулевую ограниченную core-five пилотную модель и известную ковариацию OU. Поправки смещения и подстройки критического уровня отсутствуют.']
    (ROOT/'reports/wald_se_calibration_RU.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    print(json.dumps(summary),flush=True)


def main(smoke=False):
    from kaggle.api.kaggle_api_extended import KaggleApi
    stages=['secal_smoke'] if smoke else [f'secal_{i}' for i in range(5)]
    done={}
    while len(done)<len(stages):
        statuses=[]
        try:
            api=KaggleApi();api.authenticate()
            for stage in stages:
                if stage in done:continue
                meta=json.loads((ROOT/f'kaggle/{stage}/kernel-metadata.json').read_text())
                response=api.kernels_status(meta['id']);state=getattr(response.status,'name',str(response.status))
                statuses.append({'stage':stage,'status':state})
                if state=='COMPLETE':
                    _,proof,folder=collect(stage,api);done[stage]=proof
                    if smoke:
                        diagnostic=json.loads(next(folder.glob('diagnostics_s*.json')).read_text())
                        timing={k:diagnostic[k] for k in ('score_seconds','bootstrap_seconds','total_seconds')}
                        write(ROOT/'reports/wald_se_smoke_validation.json',{'execution_validated':True,'proof':proof,'timing':timing,'scientific_calibration_assessed':False})
                        print('Smoke execution/replay validated: '+json.dumps(timing),flush=True)
                elif state in ('ERROR','CANCELED','CANCEL_ACKNOWLEDGED'):
                    write(ROOT/'reports/wald_se_job_status.json',{'statuses':statuses,'failed_stage':stage,'automatic_resubmission':False})
                    print('Failed stage '+stage+' preserved; no automatic resubmission',flush=True);return
            write(ROOT/'reports/wald_se_job_status.json',{'statuses':statuses,'completed':list(done),'automatic_resubmission':False})
        except AssertionError:
            write(ROOT/'reports/wald_se_job_status.json',{'audit_failed':True,'automatic_resubmission':False,'completed':list(done)})
            print('Output audit failed; stop collection without statistical promotion',flush=True);return
        except Exception as error:
            print('Read-only status/collection retry: '+type(error).__name__,flush=True)
        if len(done)<len(stages):time.sleep(45)
    if not smoke:analyze()


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--smoke',action='store_true');parser.add_argument('--analyze-only',action='store_true');a=parser.parse_args()
    if a.analyze_only:analyze()
    else:main(a.smoke)
