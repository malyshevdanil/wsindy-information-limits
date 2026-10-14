"""Freeze private CPU execution notebooks for the efficiency protocol."""
import argparse,base64,io,json,zipfile
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]


def build(stage):
    folder=ROOT/'kaggle'/stage;folder.mkdir(parents=True,exist_ok=True)
    if (folder/f'wsindy_{stage}.ipynb').exists():raise FileExistsError('Preserve execution snapshots')
    buffer=io.BytesIO()
    with zipfile.ZipFile(buffer,'w',zipfile.ZIP_DEFLATED) as archive:
        for directory in ('src','scripts','tests','configs'):
            for p in sorted((ROOT/directory).rglob('*')):
                if stage=='guarded_geometry_reference' and directory=='configs' and p.name!='guarded_geometry_reference.json':continue
                if stage.startswith('iq_') and directory=='configs' and p.name not in ('independent_validation.json','independent_validation_v2.json','support_wald_v2.json'):
                    continue
                if stage.startswith('secal_') and directory=='configs' and p.name not in ('wald_se_calibration.json','wald_se_references.json','support_wald_v2.json','support_extra_term.json'):
                    continue
                if stage.startswith('boundary_ablate_') and directory=='configs' and p.name not in ('boundary_ablation.json','boundary_ablation_references.json'):
                    continue
                if stage.startswith('support_wald_') and directory=='configs' and p.name not in ('support_wald.json','support_wald_v2.json','support_extra_term.json','support_wald_references.json','author_pins.json'):
                    continue
                if (stage.startswith('boundary_general_') or stage in ('endpoint_vdp5','endpoint_lorenz')) and directory=='configs' and p.name not in ('boundary_generalization.json','boundary_vdp5_references.json'):
                    continue
                if p.suffix in ('.py','.json','.m'):archive.write(p,p.relative_to(ROOT).as_posix())
        if stage.startswith('efficiency_refine'):
            p=ROOT/'configs/refinement_pilots.npz';archive.write(p,p.relative_to(ROOT).as_posix())
        cache_source=None
        if stage.startswith('crossfit_records_'):
            cache_source=json.loads((ROOT/'configs/crossfit_records.json').read_text()).get('cache_kernel_source')
        if stage.startswith('crossfit120_records_'):
            cache_source=json.loads((ROOT/'configs/crossfit120_records.json').read_text()).get('cache_kernel_source')
        if stage.startswith('crossfitA_records_'):
            cache_source=json.loads((ROOT/'configs/incidental_weak.json').read_text()).get('cache_kernel_source')
        if stage.startswith('repair') and '_records_' in stage:
            cache_source=json.loads((ROOT/'configs'/(stage.split('_')[0]+'_records.json')).read_text()).get('cache_kernel_source')
        if stage.startswith('crossfit_records_') and not cache_source:
            p=ROOT/'configs/crossfit_instruments.npz';archive.write(p,p.relative_to(ROOT).as_posix())
        if stage.startswith(('efficiency_fisher','efficiency_noise')) or stage=='moderate_likelihood':
            p=ROOT/'configs/information_pilots.npz';archive.write(p,p.relative_to(ROOT).as_posix())
    payload=base64.b64encode(buffer.getvalue()).decode()
    setup="import os, sys, subprocess, io, base64, zipfile\nfrom pathlib import Path\nfor key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'\nroot=Path('/kaggle/working/wsindy_research');root.mkdir(exist_ok=True)\n"+f"payload='{payload}'\nwith zipfile.ZipFile(io.BytesIO(base64.b64decode(payload))) as archive:archive.extractall(root)\n"
    if cache_source:
        cache_name=stage.split('_')[0] if stage.startswith('repair') else 'crossfitA' if stage.startswith('crossfitA_records_') else 'crossfit120' if stage.startswith('crossfit120_records_') else 'crossfit'
        setup+=f"import hashlib, json, shutil\nexpected=json.loads((root/'configs/{cache_name}_cache_geometry.json').read_text())['cache_sha256']\ncandidates=[p for p in Path('/kaggle/input').rglob('crossfit_instruments.npz') if hashlib.sha256(p.read_bytes()).hexdigest()==expected]\nif len(candidates)!=1:raise RuntimeError('Unique verified model-instrument cache missing')\nshutil.copyfile(candidates[0],root/'configs/{cache_name}_instruments.npz')\n"
    followup=stage.startswith('efficiency_power_') or (stage.startswith('efficiency_match') or stage=='efficiency_matched_wendy')
    if (stage.startswith('efficiency_match') or stage=='efficiency_matched_wendy' or stage=='support_extra' or stage.startswith('support_wald_')):
        setup+="subprocess.run(['apt-get','update','-qq'],check=True,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)\nsubprocess.run(['apt-get','install','-y','-qq','--no-install-recommends','octave','octave-statistics','octave-signal'],check=True,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)\n"
    runner='run_efficiency_followup.py' if followup else 'run_efficiency_validation.py'
    if stage.startswith(('efficiency_fisher','efficiency_noise')):runner='run_information_followup.py'
    if stage=='efficiency_lorenz':runner='run_lorenz_efficiency.py'
    if stage.startswith('efficiency_refine') or stage=='efficiency_fixture':runner='run_efficiency_checks.py'
    if stage.startswith('independent_records_'):runner='run_independent_records.py'
    if stage.startswith('record_bootstrap_'):runner='run_record_bootstrap.py'
    if stage.startswith('heterogeneous_records_'):runner='run_heterogeneous_records.py'
    if stage=='moderate_likelihood':runner='run_moderate_likelihood.py'
    if stage=='crossfit_geometry':runner='run_crossfit_geometry.py'
    if stage=='crossfit_cache_probe':runner='run_crossfit_cache_probe.py'
    if stage=='crossfit_polynomial_probe':runner='run_crossfit_polynomial_probe.py'
    if stage=='crossfit_poly_probe_v2':runner='run_crossfit_polynomial_probe.py'
    if stage=='crossfit_geometry_v2':runner='run_crossfit_geometry_v2.py'
    if stage=='crossfit120_geometry':runner='run_crossfit_geometry_v2.py'
    if stage=='crossfitA_geometry':runner='run_crossfit_geometry_v2.py'
    if stage=='incidental_smoke' or stage.startswith('incidental_pilot_'):runner='run_incidental_pilot.py'
    if stage=='incidental_curvature':runner='run_incidental_curvature.py'
    if stage=='inc_main_smoke' or stage.startswith(('inc_main_l','inc_trim_l','inc_short_l')):runner='run_incidental_main.py'
    if stage=='inc_forecast':runner='run_incidental_forecast.py'
    if stage=='inc_bias_check':runner='run_incidental_bias_check.py'
    if stage.startswith('boundary_screen'):runner='run_boundary_weak_screen.py'
    if stage.startswith('boundary_confirm_'):runner='run_boundary_confirmation.py'
    if stage.startswith('boundary_ablate_'):runner='run_boundary_ablation.py'
    if stage=='support_distance':runner='run_support_distance_screen.py'
    if stage=='support_extra':runner='run_support_extra_term.py'
    if stage.startswith('support_wald_'):runner='run_support_wald.py'
    if stage.startswith('secal_'):runner='run_wald_se_calibration.py'
    if stage.startswith('iq_'):runner='run_independent_validation.py'
    if stage=='guarded_geometry_reference':runner='run_guarded_geometry_reference.py'
    if stage.startswith('boundary_general_'):runner='run_boundary_generalization.py'
    if stage in ('endpoint_vdp5','endpoint_lorenz'):runner='run_boundary_generalization.py'
    if stage.startswith('crossfit_records_'):runner='run_crossfit_records.py'
    if stage.startswith('crossfit120_records_'):runner='run_crossfit_records.py'
    if stage.startswith('crossfitA_records_'):runner='run_crossfit_records.py'
    if stage.startswith('repair') and '_records_' in stage:runner='run_crossfit_records.py'
    if stage.startswith('repair') and stage.endswith('_geometry'):runner='run_crossfit_geometry_v2.py'
    if stage in ('cache_fix_probe','cache_fix_a','cache_fix_v4','cache120a'):runner='run_cache_repair_probe.py'
    run=f"subprocess.run([sys.executable,'scripts/{runner}','--stage','{stage}','--output','/kaggle/working/results_{stage}'],cwd=root,check=True)\n"
    finish=f"import shutil\nshutil.make_archive('/kaggle/working/wsindy_{stage}_outputs','zip','/kaggle/working/results_{stage}')\n"
    cells=[]
    for i,source in enumerate((setup,run,finish)):
        compile(source,f'{stage}-{i}','exec');cells.append({'cell_type':'code','id':f'{stage}-{i}','metadata':{},'execution_count':None,'outputs':[],'source':source.splitlines(keepends=True)})
    notebook={'nbformat':4,'nbformat_minor':5,'metadata':{'kernelspec':{'display_name':'Python 3','language':'python','name':'python3'}},'cells':cells}
    name=f'wsindy_{stage}.ipynb';(folder/name).write_text(json.dumps(notebook),encoding='utf-8')
    meta={'id':'malyshevdanil/wsindy-sampling-saturation-'+stage.replace('_','-'),'title':('WSINDy Sampling Saturation '+stage.replace('_',' ').title())[:50],'code_file':name,'language':'python','kernel_type':'notebook','is_private':True,'enable_gpu':False,'enable_internet':(stage.startswith('efficiency_match') or stage=='efficiency_matched_wendy' or stage=='support_extra' or stage.startswith('support_wald_')),'dataset_sources':[],'competition_sources':[],'kernel_sources':[]}
    if cache_source:meta['kernel_sources']=[cache_source]
    meta['id']=meta['id'].lower()
    (folder/'kernel-metadata.json').write_text(json.dumps(meta,indent=2))
    print(stage+' frozen')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage');build(p.parse_args().stage)
