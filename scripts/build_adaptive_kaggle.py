"""Private Kaggle notebooks for pinned full author pipelines, without credentials."""
import argparse
import base64
import io
import json
from pathlib import Path
import zipfile

ROOT = Path(__file__).resolve().parents[1]


def build(stage):
    folder = ROOT / 'kaggle' / stage
    folder.mkdir(parents=True, exist_ok=True)
    files = [p for name in ('src', 'scripts', 'configs') for p in (ROOT/name).rglob('*')
             if p.suffix in ('.py', '.m', '.json') and '__pycache__' not in p.parts]
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, 'w', zipfile.ZIP_DEFLATED) as archive:
        for p in files:
            info = zipfile.ZipInfo(p.relative_to(ROOT).as_posix(), (2026,10,11,0,0,0))
            info.compress_type = zipfile.ZIP_DEFLATED
            archive.writestr(info, p.read_bytes())
    setup = """import os, sys, subprocess, io, base64, zipfile, json
from pathlib import Path
for name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[name] = '1'
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2'
os.environ['TF_USE_LEGACY_KERAS'] = '1'
root = Path('/kaggle/working/wsindy_research')
root.mkdir(exist_ok=True)
""" + f"with zipfile.ZipFile(io.BytesIO(base64.b64decode('{base64.b64encode(buffer.getvalue()).decode()}'))) as archive:\n    archive.extractall(root)\n"
    probe = """import importlib.metadata, shutil
info = {'packages': {p: importlib.metadata.version(p) for p in ('numpy','scipy','tensorflow')},
        'uid': os.getuid(), 'octave_before': shutil.which('octave')}
try: info['tf_keras'] = importlib.metadata.version('tf-keras')
except importlib.metadata.PackageNotFoundError: info['tf_keras'] = None
import tensorflow as tf
info['gpu'] = [d.name for d in tf.config.list_physical_devices('GPU')]
print(json.dumps(info), flush=True)
subprocess.run(['apt-get','update','-qq'], check=True)
subprocess.run(['apt-get','install','-y','-qq','--no-install-recommends','octave','octave-statistics','octave-signal'], check=True, stdout=subprocess.DEVNULL)
info['octave'] = subprocess.check_output(['octave','--version'],text=True).splitlines()[0]
Path('/kaggle/working/environment.json').write_text(json.dumps(info,indent=2))
print(json.dumps(info), flush=True)
"""
    run = probe
    if stage != 'adaptive_probe':
        run += f"subprocess.run([sys.executable, 'scripts/run_full_adaptive.py', '--stage', '{stage}', '--output', '/kaggle/working/results_{stage}'], cwd=root, check=True)\n"
    finish = f"import shutil\nshutil.make_archive('/kaggle/working/wsindy_{stage}_outputs', 'zip', '/kaggle/working/{'results_'+stage if stage!='adaptive_probe' else 'wsindy_research'}')\n"
    cells=[]
    for i,source in enumerate((setup,run,finish)):
        compile(source, f'{stage}-cell-{i}', 'exec')
        cells.append({'cell_type':'code','id':f'{stage}-{i}','metadata':{},'execution_count':None,'outputs':[],
                      'source':source.splitlines(keepends=True)})
    title='WSINDy Sampling Saturation '+stage.replace('_',' ').title()
    notebook={'nbformat':4,'nbformat_minor':5,'metadata':{'kernelspec':{'display_name':'Python 3','language':'python','name':'python3'}},'cells':cells}
    name=f'wsindy_{stage}.ipynb'
    (folder/name).write_text(json.dumps(notebook),encoding='utf-8')
    meta={'id':'malyshevdanil/'+title.lower().replace(' ','-'),'title':title,'code_file':name,
          'language':'python','kernel_type':'notebook','is_private':True,'enable_gpu':False,
          'enable_internet':True,'dataset_sources':[],'competition_sources':[],'kernel_sources':[]}
    (folder/'kernel-metadata.json').write_text(json.dumps(meta,indent=2),encoding='utf-8')
    print(json.dumps({'stage':stage,'files':len(files),'archive_bytes':len(buffer.getvalue())}))


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('stage', choices=['adaptive_probe','adaptive_smoke','adaptive_octave_smoke','adaptive_classical','adaptive_wm',*[f'adaptive_wm_{i}' for i in range(4)]])
    build(parser.parse_args().stage)
