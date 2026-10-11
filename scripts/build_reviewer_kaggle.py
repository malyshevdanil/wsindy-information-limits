"""Freeze the reviewer calculations in a private CPU notebook."""
import base64
import io
import json
import zipfile
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
folder=ROOT/'kaggle/reviewer';folder.mkdir(parents=True,exist_ok=True)
buffer=io.BytesIO()
with zipfile.ZipFile(buffer,'w',zipfile.ZIP_DEFLATED) as archive:
    for directory in ('src','scripts','tests','configs'):
        for p in sorted((ROOT/directory).glob('*')):
            if p.suffix in ('.py','.json'):
                archive.write(p,p.relative_to(ROOT).as_posix())
payload=base64.b64encode(buffer.getvalue()).decode()
setup=f'''import os, sys, subprocess, io, base64, zipfile
from pathlib import Path
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):
    os.environ[key]='1'
root=Path('/kaggle/working/wsindy_research')
root.mkdir(exist_ok=True)
payload='{payload}'
with zipfile.ZipFile(io.BytesIO(base64.b64decode(payload))) as archive:
    archive.extractall(root)
'''
run="subprocess.run([sys.executable,'scripts/run_reviewer_validation.py','--output','/kaggle/working/results_reviewer'],cwd=root,check=True)\n"
finish="import shutil\nshutil.make_archive('/kaggle/working/wsindy_reviewer_outputs','zip','/kaggle/working/results_reviewer')\n"
cells=[]
for i,source in enumerate((setup,run,finish)):
    compile(source,f'reviewer-{i}','exec')
    cells.append({'cell_type':'code','id':f'reviewer-{i}','metadata':{},'execution_count':None,'outputs':[],'source':source.splitlines(keepends=True)})
notebook={'nbformat':4,'nbformat_minor':5,'metadata':{'kernelspec':{'display_name':'Python 3','language':'python','name':'python3'}},'cells':cells}
(folder/'wsindy_reviewer.ipynb').write_text(json.dumps(notebook),encoding='utf-8')
(folder/'kernel-metadata.json').write_text(json.dumps({'id':'malyshevdanil/wsindy-sampling-saturation-reviewer','title':'WSINDy Sampling Saturation Reviewer',
    'code_file':'wsindy_reviewer.ipynb','language':'python','kernel_type':'notebook','is_private':True,
    'enable_gpu':False,'enable_internet':False,'dataset_sources':[],'competition_sources':[],'kernel_sources':[]},indent=2))
print('Reviewer notebook frozen')
