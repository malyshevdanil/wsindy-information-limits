"""Insert proved cross-fitting extension, without inserting pending MC results."""
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'manuscript'


def main():
    path=BASE/'article_revised.md';text=path.read_text(encoding='utf-8')
    section=r'''**Corollary 2 (record-wise cross-fitting; assumptions and proof in S16.4).** Split independent records into a fixed number of balanced folds before observing the data. For each held-out fold, fit a coefficient pilot on its complement and construct the independent-view instruments from that pilot. If the pilots converge to a deterministic limit, their score maps satisfy the uniform continuity and envelope conditions (A6), and the quadrature defects satisfy (A5), the pooled estimator has the same first-order Gaussian limit and consistent record-score sandwich as the deterministic limiting instrument map. Every record enters the final estimating equation, with no external training sample.

The key mean identity holds for every instrument pilot: conditional on a fold complement, each held-out score is centered. Consistent complement pilots therefore change the normalized score sum by $o_P(1)$. Independence is used separately within each fold. Cross-fitting follows established sample-reuse methodology [@chernozhukov2018]; its role here is to preserve the corrected weak moment identity with heterogeneous initial states. The known covariance assumption remains in force. A record bootstrap can refit all complement pilots and retain every physical record's copies within its original fold; its finite-sample calibration requires numerical assessment.

'''
    if '**Corollary 2 (record-wise cross-fitting;' not in text:
        text=text.replace('### 7.2. Model weights',section+'### 7.2. Model weights')
        text=text.replace('Its local information reference is the sum of the individual nuisance-adjusted informations.',
            'Its local information reference is the sum of the individual nuisance-adjusted informations. Corollary 2 reuses all records through fold-wise independent pilots under explicit continuity conditions.')
    text=text.replace('Independent instruments preserve corrected weak moment identities at fixed noise and permit record-level sandwich inference without fitting every initial state in the estimating equation.',
        'Independent instruments preserve corrected weak moment identities at fixed noise. Fold-wise cross-fitting reuses all records and permits record-level sandwich inference without fitting every initial state in the estimating equation.')
    path.write_text(text,encoding='utf-8')
    proof=(BASE/'heterogeneous_records_theory.md').read_text(encoding='utf-8')
    extension=proof[proof.index('## S16.'):proof.index('## References')]
    path=BASE/'supplement.md';text=path.read_text(encoding='utf-8')
    start=text.index('## S16.');end=text.index('## S17.',start)
    text=text[:start]+extension+'\n'+text[end:];path.write_text(text,encoding='utf-8')


if __name__=='__main__':main()
