"""Download a portable official Tectonic release into this workspace."""
import hashlib
import io
import json
from pathlib import Path
import requests
import zipfile
ROOT=Path(__file__).resolve().parents[1]
folder=ROOT/'tools/tectonic'
folder.mkdir(parents=True,exist_ok=True)
url='https://github.com/tectonic-typesetting/tectonic/releases/download/tectonic%400.17.0/tectonic-0.17.0-x86_64-pc-windows-msvc.zip'
r=requests.get(url,timeout=(30,180));r.raise_for_status()
with zipfile.ZipFile(io.BytesIO(r.content)) as archive:
    for member in archive.infolist():
        if not (folder/member.filename).resolve().is_relative_to(folder.resolve()):raise ValueError('Unsafe member')
    archive.extractall(folder)
(folder/'download_manifest.json').write_text(json.dumps({'url':url,'sha256':hashlib.sha256(r.content).hexdigest()},indent=2))
print('Portable Tectonic installed')
