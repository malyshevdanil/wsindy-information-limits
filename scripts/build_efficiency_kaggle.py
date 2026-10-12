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
                if p.suffix in ('.py','.json','.m'):archive.write(p,p.relative_to(ROOT).as_posix())
        if stage.startswith('efficiency_refine'):
            p=ROOT/'configs/refinement_pilots.npz';archive.write(p,p.relative_to(ROOT).as_posix())
        if stage.startswith(('efficiency_fisher','efficiency_noise')) or stage=='moderate_likelihood':
            p=ROOT/'configs/information_pilots.npz';archive.write(p,p.relative_to(ROOT).as_posix())
    payload=base64.b64encode(buffer.getvalue()).decode()
    setup="import os, sys, subprocess, io, base64, zipfile\nfrom pathlib import Path\nfor key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'\nroot=Path('/kaggle/working/wsindy_research');root.mkdir(exist_ok=True)\n"+f"payload='{payload}'\nwith zipfile.ZipFile(io.BytesIO(base64.b64decode(payload))) as archive:archive.extractall(root)\n"
    followup=stage.startswith('efficiency_power_') or (stage.startswith('efficiency_match') or stage=='efficiency_matched_wendy')
    if (stage.startswith('efficiency_match') or stage=='efficiency_matched_wendy'):
        setup+="subprocess.run(['apt-get','update','-qq'],check=True,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)\nsubprocess.run(['apt-get','install','-y','-qq','--no-install-recommends','octave','octave-statistics','octave-signal'],check=True,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)\n"
    runner='run_efficiency_followup.py' if followup else 'run_efficiency_validation.py'
    if stage.startswith(('efficiency_fisher','efficiency_noise')):runner='run_information_followup.py'
    if stage=='efficiency_lorenz':runner='run_lorenz_efficiency.py'
    if stage.startswith('efficiency_refine') or stage=='efficiency_fixture':runner='run_efficiency_checks.py'
    if stage.startswith('independent_records_'):runner='run_independent_records.py'
    if stage.startswith('record_bootstrap_'):runner='run_record_bootstrap.py'
    if stage.startswith('heterogeneous_records_'):runner='run_heterogeneous_records.py'
    if stage=='moderate_likelihood':runner='run_moderate_likelihood.py'
    run=f"subprocess.run([sys.executable,'scripts/{runner}','--stage','{stage}','--output','/kaggle/working/results_{stage}'],cwd=root,check=True)\n"
    finish=f"import shutil\nshutil.make_archive('/kaggle/working/wsindy_{stage}_outputs','zip','/kaggle/working/results_{stage}')\n"
    cells=[]
    for i,source in enumerate((setup,run,finish)):
        compile(source,f'{stage}-{i}','exec');cells.append({'cell_type':'code','id':f'{stage}-{i}','metadata':{},'execution_count':None,'outputs':[],'source':source.splitlines(keepends=True)})
    notebook={'nbformat':4,'nbformat_minor':5,'metadata':{'kernelspec':{'display_name':'Python 3','language':'python','name':'python3'}},'cells':cells}
    name=f'wsindy_{stage}.ipynb';(folder/name).write_text(json.dumps(notebook),encoding='utf-8')
    meta={'id':'malyshevdanil/wsindy-sampling-saturation-'+stage.replace('_','-'),'title':('WSINDy Sampling Saturation '+stage.replace('_',' ').title())[:50],'code_file':name,'language':'python','kernel_type':'notebook','is_private':True,'enable_gpu':False,'enable_internet':(stage.startswith('efficiency_match') or stage=='efficiency_matched_wendy'),'dataset_sources':[],'competition_sources':[],'kernel_sources':[]}
    (folder/'kernel-metadata.json').write_text(json.dumps(meta,indent=2))
    print(stage+' frozen')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage');build(p.parse_args().stage)
