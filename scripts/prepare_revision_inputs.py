"""Prepare private Kaggle input data; no credentials or external mutations."""
import hashlib
import json
from pathlib import Path
import shutil
ROOT=Path(__file__).resolve().parents[1]
out=ROOT/'kaggle/revision_inputs'
out.mkdir(parents=True,exist_ok=True)
files=sorted((ROOT/'results/kaggle_extended').glob('residual_*.npz'))
if len(files)!=6:raise ValueError('Six input files required')
hashes={}
for file in files:
    shutil.copy2(file,out/file.name)
    hashes[file.name]=hashlib.sha256(file.read_bytes()).hexdigest()
(out/'input_manifest.json').write_text(json.dumps(hashes,indent=2))
(out/'dataset-metadata.json').write_text(json.dumps({'title':'WSINDy Revision Covariance Inputs',
    'id':'malyshevdanil/wsindy-revision-covariance-inputs','licenses':[{'name':'CC0-1.0'}]},indent=2))
print(json.dumps({'files':len(files),'directory':str(out)}))
