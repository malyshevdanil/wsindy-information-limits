"""Add primary-source metadata for the incidental-parameter comparison."""
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
ENTRIES=[
    '@article{neyman1948,author={Neyman, Jerzy and Scott, Elizabeth L.},title={Consistent Estimates Based on Partially Consistent Observations},journal={Econometrica},year={1948},volume={16},number={1},pages={1--32},doi={10.2307/1914288}}',
    '@article{lancaster2000,author={Lancaster, Tony},title={The incidental parameter problem since 1948},journal={Journal of Econometrics},year={2000},volume={95},number={2},pages={391--413},doi={10.1016/S0304-4076(99)00044-5}}',
    '@article{hahn2004,author={Hahn, Jinyong and Newey, Whitney K.},title={Jackknife and Analytical Bias Reduction for Nonlinear Panel Models},journal={Econometrica},year={2004},volume={72},number={4},pages={1295--1319},doi={10.1111/j.1468-0262.2004.00533.x}}',
    '@incollection{arellano2007,author={Arellano, Manuel and Hahn, Jinyong},title={Understanding Bias in Nonlinear Panel Models: Some Recent Developments},booktitle={Advances in Economics and Econometrics: Theory and Applications, Ninth World Congress},editor={Blundell, Richard and Newey, Whitney and Persson, Torsten},publisher={Cambridge University Press},year={2007},volume={3},pages={381--409},doi={10.1017/CBO9780511607547.013}}'
]


def main():
    p=ROOT/'manuscript/publication_references.bib';text=p.read_text(encoding='utf-8')
    for entry in ENTRIES:
        key=entry.split('{',1)[1].split(',',1)[0]
        if '{'+key+',' not in text:text+='\n'+entry+'\n'
    p.write_text(text,encoding='utf-8')


if __name__=='__main__':main()
