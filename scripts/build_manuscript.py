"""Compatibility entry point; native equations are produced by Pandoc."""
from pathlib import Path
import subprocess
import sys
import pypandoc
ROOT=Path(__file__).resolve().parents[1]
if __name__=='__main__':
    subprocess.run([sys.executable,str(ROOT/'scripts/build_publication.py')],check=True,cwd=ROOT)
    report=ROOT/'reports/research_progress_RU.md'
    pypandoc.convert_file(str(report),'docx',format='markdown+tex_math_dollars+tex_math_single_backslash',outputfile=str(report.with_suffix('.docx')))
