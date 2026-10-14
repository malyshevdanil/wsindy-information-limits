"""Recompute paired uncertainty from byte-preserved original binary outcomes."""
import json,hashlib,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from src.paired_binary import paired_difference_interval


def main():
    bundle=ROOT/'data/endpoint_followups';manifest=json.loads((bundle/'manifest.json').read_text())
    for f in manifest['files']:
        assert hashlib.sha256((bundle/f['path']).read_bytes()).hexdigest()==f['sha256']
    rows=[]
    for i in range(4):rows+=json.loads((bundle/f'support_wald_v2_{i}/rows.json').read_text())
    old=json.loads((ROOT/'reports/support_wald_results.json').read_text());groups=[]
    for level in [.005,.05]:
        for di in range(8):
            vectors={}
            for method in ['weak_gls_wald','author_mstls']:
                sample=sorted((r for r in rows if r['level']==level and r['delta_index']==di and r['method']==method),key=lambda r:r['seed'])
                assert [r['seed'] for r in sample]==list(range(100))
                vectors[method]=np.array([r['included'] for r in sample])
            oldgroup=next(g for g in old['groups'] if g['level']==level and g.get('delta_index')==di and g['method']=='paired_difference')
            interval=paired_difference_interval(vectors['weak_gls_wald'],vectors['author_mstls'])
            assert abs(interval['difference']-oldgroup['wald_minus_mstls'])<1e-12
            groups.append({'level':level,'delta_index':di,'magnitude':abs(old['protocol']['multipliers'][di]),**interval,
                'archived_bootstrap_ci95':oldgroup['difference_ci95']})
    result={'method':groups[0]['method'],'all_1600_paired_outcomes_retained':True,'counts_unchanged':True,
        'confidence_intervals_pointwise':True,'post_study_uncertainty_reanalysis':True,'groups':groups}
    (ROOT/'reports/support_wald_paired_intervals.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    p=ROOT/'manuscript/supplement.md';text=p.read_text(encoding='utf-8')
    start=text.index('Table S19.');end=text.index('## S24.',start)
    original=text[start:end];lines=original.splitlines()
    new=[];index=0
    for line in lines:
        if line.startswith('| 0.5% |') or line.startswith('| 5% |'):
            g=groups[index];index+=1;parts=line.split('|');lo,hi=g['ci95']
            parts[-2]=f" {100*g['difference']:.0f} ({100*lo:.1f}\u2013{100*hi:.1f}) "
            line='|'.join(parts)
        new.append(line)
    assert index==16
    section='\n'.join(new)+'\n\n'
    section=section.replace('difference intervals use 5000 paired bootstrap resamples.',
        'difference intervals are conservative pointwise 95% intervals obtained from Bonferroni-projected Clopper\u2013Pearson limits for the two discordant cells. S23.2 defines them; they remain nondegenerate at extreme samples.')
    text=text[:start]+section+text[end:];p.write_text(text,encoding='utf-8')
    print(json.dumps(groups),flush=True)


if __name__=='__main__':main()
