"""Build a private self-contained notebook. This command does not contact Kaggle."""
import argparse
import base64
import io
import json
from pathlib import Path
import re
import zipfile

ROOT = Path(__file__).resolve().parents[1]


def cell(kind, source):
    result = {'cell_type': kind, 'metadata': {}, 'source': source.splitlines(keepends=True)}
    if kind == 'code':
        result.update({'execution_count': None, 'outputs': []})
    return result


def build(username, stage):
    if username != 'YOUR_KAGGLE_USERNAME' and not re.fullmatch('[A-Za-z0-9_-]+', username):
        raise ValueError('Invalid Kaggle username')
    output = ROOT / 'kaggle' / stage
    output.mkdir(parents=True, exist_ok=True)
    files = [ROOT / 'requirements.txt']
    for folder, pattern in [('src', '*.py'), ('scripts', '*.py'), ('tests', '*.py'), ('configs', '*.json')]:
        files += sorted((ROOT / folder).glob(pattern))
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, 'w', compression=zipfile.ZIP_DEFLATED) as archive:
        for path in files:
            name = path.relative_to(ROOT).as_posix()
            if stage == 'revision' and path.suffix == '.npz':
                name = 'inputs/covariance/' + path.name
            info = zipfile.ZipInfo(name, date_time=(2026, 10, 11, 0, 0, 0))
            info.compress_type = zipfile.ZIP_DEFLATED
            archive.writestr(info, path.read_bytes())
    encoded = base64.b64encode(buffer.getvalue()).decode('ascii')
    unpack = (
        'import os, sys, subprocess, base64, io, zipfile\n'
        'from pathlib import Path\n'
        "for name in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):\n"
        "    os.environ[name] = '1'\n"
        "root = Path('/kaggle/working/wsindy_research')\n"
        "root.mkdir(parents=True, exist_ok=True)\n"
        f"payload = '{encoded}'\n"
        'with zipfile.ZipFile(io.BytesIO(base64.b64decode(payload))) as archive:\n'
        '    for entry in archive.infolist():\n'
        '        destination = (root / entry.filename).resolve()\n'
        '        if not destination.is_relative_to(root.resolve()):\n'
        "            raise ValueError('Archive path outside workspace')\n"
        '    archive.extractall(root)\n'
        "print('Source unpacked:', root)\n"
    )
    verify = (
        'import importlib.metadata\n'
        "print({p: importlib.metadata.version(p) for p in ('numpy', 'scipy', 'matplotlib', 'pytest')})\n"
        "subprocess.run([sys.executable, '-m', 'pytest', '-q', 'tests'], cwd=root, check=True)\n"
    )
    run = (
        f"subprocess.run([sys.executable, 'scripts/run_experiments.py', '--config', 'configs/{stage}.json', "
        f"'--output', '/kaggle/working/results_{stage}'], cwd=root, check=True)\n"
        f"subprocess.run([sys.executable, 'scripts/summarize_results.py', '/kaggle/working/results_{stage}'], "
        'cwd=root, check=True)\n'
    )
    if stage == 'extended':
        run = ("subprocess.run([sys.executable, 'scripts/run_extended_validation.py', '--output', "
               "'/kaggle/working/results_extended'], cwd=root, check=True)\n")
    if stage == 'replication':
        run = ("subprocess.run([sys.executable, 'scripts/run_replication_comparison.py', '--output', "
               "'/kaggle/working/results_replication'], cwd=root, check=True)\n")
    if stage in ('revision', 'revision_fast'):
        run = ("import shutil\n"
               "input_dir = root / 'inputs/covariance'\n"
               "input_dir.mkdir(parents=True, exist_ok=True)\n"
               "for file in Path('/kaggle/input').rglob('residual_*.npz'):\n"
               "    shutil.copy2(file, input_dir / file.name)\n"
               f"subprocess.run([sys.executable, 'scripts/{'run_revision_fast.py' if stage=='revision_fast' else 'run_revision_validation.py'}', '--output', "
               f"'/kaggle/working/results_{stage}'], cwd=root, check=True)\n")
    if stage == 'reweighting':
        run = ("subprocess.run([sys.executable, 'scripts/run_reweighting_validation.py', '--output', "
               "'/kaggle/working/results_reweighting'], cwd=root, check=True)\n")
    if stage == 'geometry':
        run = ("subprocess.run([sys.executable, 'scripts/run_geometry_validation.py', '--output', "
               "'/kaggle/working/results_geometry'], cwd=root, check=True)\n")
    if stage == 'covariance':
        run += ("subprocess.run([sys.executable, 'scripts/run_scalar.py', '--output', "
                "'/kaggle/working/results_covariance/scalar', '--repeats', '1000'], cwd=root, check=True)\n")
    archive_output = (
        'import shutil\n'
        f"shutil.make_archive('/kaggle/working/wsindy_{stage}_outputs', 'zip', '/kaggle/working/results_{stage}')\n"
        "print('Saved manifest, records and archive in /kaggle/working')\n"
    )
    title = f"WSINDy Sampling Saturation {stage.replace('_',' ').title()}"
    notebook = {'nbformat': 4, 'nbformat_minor': 5,
                'metadata': {'kernelspec': {'display_name': 'Python 3', 'language': 'python', 'name': 'python3'},
                             'language_info': {'name': 'python'}},
                'cells': [cell('markdown', f'# {title}\n\nPrivate CPU research run. Synthetic data only. '
                              'No credentials are embedded. Correctness tests run before experiments. '
                              'The smoke stage is an installation check; the pilot is exploratory.\n'),
                          cell('code', unpack), cell('code', verify), cell('code', run),
                          cell('code', archive_output)]}
    for i, entry in enumerate(notebook['cells']):
        entry['id'] = f'wsindy-{stage}-{i}'
    name = f'wsindy_{stage}.ipynb'
    (output / name).write_text(json.dumps(notebook, ensure_ascii=False, indent=1), encoding='utf-8')
    slug = title.lower().replace(' ', '-')
    metadata = {'id': f'{username}/{slug}', 'title': title, 'code_file': name,
                'language': 'python', 'kernel_type': 'notebook', 'is_private': True,
                'enable_gpu': False, 'enable_internet': False,
                'dataset_sources': [], 'competition_sources': [], 'kernel_sources': []}
    if stage in ('revision', 'revision_fast'):
        metadata['dataset_sources'] = [f'{username}/wsindy-revision-covariance-inputs']
    (output / 'kernel-metadata.json').write_text(json.dumps(metadata, indent=2), encoding='utf-8')
    print(json.dumps({'directory': str(output), 'stage': stage, 'private': True,
                      'packed_files': len(files), 'archive_bytes': len(buffer.getvalue())}))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--username', default='YOUR_KAGGLE_USERNAME')
    parser.add_argument('--stage', choices=['smoke', 'pilot', 'covariance', 'extended', 'replication', 'revision', 'revision_fast', 'reweighting', 'geometry'], default='smoke')
    arguments = parser.parse_args()
    build(arguments.username, arguments.stage)
