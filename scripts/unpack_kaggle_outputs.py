"""Extract validated downloaded archives into their explicit stage directories."""
import argparse
from pathlib import Path
import zipfile
ROOT=Path(__file__).resolve().parents[1]
parser=argparse.ArgumentParser();parser.add_argument('stage',choices=['moderate_likelihood',*[f'heterogeneous_records_{i}' for i in range(4)],*[f'record_bootstrap_{i}' for i in range(4)],*[f'independent_records_{i}' for i in range(4)],'revision','revision_fast','reweighting','geometry','reviewer','efficiency_refine_0','efficiency_refine_1','efficiency_fixture','efficiency_fisher_v2','efficiency_noise_v2','efficiency_fisher','efficiency_noise','efficiency_lorenz',*[f'efficiency_match_{i}' for i in range(4)],'efficiency_match','efficiency_matched_wendy',*[f'efficiency_power_{i}' for i in range(4)],'efficiency_scalar','efficiency_scalar_v2','efficiency_scalar_v3','efficiency_geometry','efficiency_geometry_v2','efficiency_vdp_0','efficiency_vdp_1','efficiency_vdp_polish_0','efficiency_vdp_polish_1','efficiency_vdp_final_0','efficiency_vdp_final_1','adaptive_smoke','adaptive_octave_smoke','adaptive_classical','adaptive_wm',*[f'adaptive_wm_{i}' for i in range(4)]])
stage=parser.parse_args().stage
out=(ROOT/f'results/kaggle_{stage}').resolve();out.mkdir(parents=True,exist_ok=True)
with zipfile.ZipFile(ROOT/f'results/kaggle_{stage}_download/wsindy_{stage}_outputs.zip') as archive:
    if archive.testzip() is not None:raise ValueError('CRC failure')
    for entry in archive.infolist():
        if not (out/entry.filename).resolve().is_relative_to(out):raise ValueError('Unsafe output member')
    archive.extractall(out)
print(f'{stage}: extracted')
