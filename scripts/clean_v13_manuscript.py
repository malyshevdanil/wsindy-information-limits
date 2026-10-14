"""Apply review corrections without inserting results of unfinished experiments."""
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]


def main():
    for name in ('article_revised.md', 'heterogeneous_records_theory.md',
                 'independent_records_theory.md', 'supplement.md'):
        path = ROOT/'manuscript'/name
        text = path.read_text(encoding='utf-8')
        text = text.replace('General dual influences act through the $L^2(P)$ closure of integrated Gaussian observations, through Gaussian observations.',
                            'General dual influences are defined by $L^2(P)$ limits of integrated Gaussian observations.')
        text = text.replace('Lemma 2', 'Lemma S14.1').replace('; Authorship', '; authorship')
        text = text.replace('Under (A1)–(A5), conditional on the independent training data,',
                            'With known sensor covariance $C_{\\nu_0}$ and under (A1)–(A5), conditional on the independent training data,')
        text = text.replace('Independent instruments preserve corrected weak moment identities at fixed noise.',
                            'With known sensor covariance, independent instruments preserve corrected weak moment identities at fixed noise.')
        text = text.replace('Its bootstrap at 50 records overcovers and greatly expands the confidence region, motivating the richer model-weighted heterogeneous-record experiment.',
                            'At 50 records, joint coverage is 52.4% for Wald regions and 99.9% for the studentized record bootstrap; the latter expands mean region volume by a factor of 248 (S17). This calibration diagnostic motivates assessment of richer model-weighted heterogeneous-record equations.')
        text = text.replace('For a fixed moment family with correct marginal correction, the independent-record variance under weights $W$ is $P_W\\Omega_{\\mathrm{true}}P_W^\\top/n$. The same sandwich describes the leading small-noise single-record variance using the order-one residual covariance.',
                            'For heterogeneous independent records, correctly centered instrumented scores have asymptotic coefficient covariance $B^{-1}V_\\psi B^{-\\top}/n$ (Theorem 3). For a single record, the leading small-noise covariance under weights $W$ is $P_W\\Omega_{\\mathrm{true},1}P_W^\\top$, using the order-one residual covariance.')
        text = text.replace('Replication allows consistent estimation of both OU range and variances under the shared-trajectory construction.',
                            'The pair-difference calibration in Section 7.3 requires equal trajectories. With heterogeneous initial states, covariance is supplied in Theorem 3; full covariance calibration requires a separate argument.')
        text = text.replace('Model trajectories at 25 initial-state grid points',
                            'Model trajectories at a 5-by-5 instrument grid (25 initial-state points)')
        text = text.replace('The 9-by-9 sensitivity-information grid has',
                            'A separate 9-by-9 sensitivity-information grid (81 benchmark points) has')
        text = text.replace('The tagged snapshot identifies the version supporting this article;',
                            'The v1.2.0 snapshot preserves the preceding numerical benchmark. The current cross-fitting and incidental-state extension is under validation and will have a separate versioned archive;')
        if name == 'article_revised.md':
            marker = '### 8.5.'
            paragraph = ('On these same 128 records at 5% noise, unknown-initial local likelihood has RMSE 0.0172, author WENDy 0.0241 and weak GLS 0.0256. The weak estimates have approximately 40% and 49% higher RMSE than likelihood, respectively. This finite-noise comparison lies within the same nuisance-initial-state information class; local small-noise efficiency alone does not predict the moderate-noise ordering.\n\n')
            if paragraph not in text:
                text = text.replace(marker, paragraph+marker, 1)
            marker = '### 7.4.'
            qv = ('For pure OU measurement error and a continuously differentiable deterministic trajectory, quadratic variation removes the finite-variation signal: $[Y]_T=2\\sigma^2T/\\tau$ coordinatewise. Thus infill realized quadratic variation divided by $2T$ consistently estimates $\\sigma^2/\\tau$ from one record, without fitting an initial state. This identifies the microergodic ratio; full covariance calibration for heterogeneous records remains open.\n\n')
            # Insert within covariance subsection, before the numerical section.
            if qv not in text:
                at = text.index('## 8.')
                text = text[:at]+qv+text[at:]
        path.write_text(text, encoding='utf-8')


if __name__ == '__main__':
    main()
