"""Remove project-history prose and place disclosures in the main article."""
from pathlib import Path
import json
import re
ROOT=Path(__file__).resolve().parents[1]


def prepare():
    p=ROOT/'manuscript/article_revised.md';text=p.read_text(encoding='utf-8')
    if 'Efficient weak-form estimation' in text:
        print('Efficiency revision uses prepare_efficiency_manuscript.py and analyze_efficiency.py')
        return
    if '### 4.3. Optimal combination within a finite window family' not in text:
        section=r'''### 4.3. Optimal combination within a finite window family

For scalar exponential dynamics, let $w_k=-\phi'_k-a\phi_k$. Define the normalizing moments and residual covariance by

$$
D_k=\int\phi_k x_a\,dt,\qquad
C_{k\ell}=\iint w_k(s)C(s,t)w_\ell(t)\,ds\,dt.
$$

A combined Gaussian influence $\ell=\sum_k\alpha_kw_k$ is normalized by $\alpha^TD=1$. If the residual covariance matrix $C$ is positive definite, minimizing $\alpha^TC\alpha$ under this constraint gives

$$
\alpha_* = \frac{C^{-1}D}{D^TC^{-1}D},\qquad
V_*=(D^TC^{-1}D)^{-1},\qquad
\operatorname{Eff}_{\mathcal F}=\frac{D^TC^{-1}D}{\mathcal I_T(a)}.
$$

This follows either from a Lagrange multiplier or covariance Cauchy–Schwarz. With a singular matrix, the inverse is restricted to its retained positive eigenspace and the constraint must remain feasible. The weights use the generating coefficient and known covariance: this is an oracle local optimum within the specified span, rather than a demonstrated adaptive estimator or a global optimum over all compactly supported functions. Every member still annihilates initial-state sensitivity, so $V_*\geq\mathcal I_{a\mid u}^{-1}$. Section 7.8 reports declared families, spectral cutoffs and quadrature checks.

'''
        text=text.replace('### 4.4. Initial-state information lost by weak influences',section+'### 4.4. Initial-state information lost by weak influences')
    start=text.index('The completed exploratory studies use') if 'The completed exploratory studies use' in text else -1
    if start>=0:
        end=text.index('<!-- REVISION_RESULTS_START -->',start)
        text=text[:start]+'''The numerical experiments use exponential decay, damped, Van der Pol and Duffing oscillators. Fixed windows and nested grids distinguish numerical resolution from independent-noise averaging. They evaluate Gaussian likelihood, sparse weak regression, polynomial residual and coefficient covariance, covariance reweighting, heterogeneous-record aggregation and complete adaptive pipelines. Known marginal variance is used where explicitly stated. Settings and provenance are detailed in the supplement.

'''+text[end:]
    text=text.replace('An additional fresh-path comparison','An additional comparison')
    text=text.replace('and 16 new independent paths','and 16 independent paths')
    text=text.replace('The earlier grid diagnostic did not establish general accuracy preservation and is documented as a negative supplementary result. Historical replication is exploratory. ','')
    disclosures='''## Declaration of generative AI and AI-assisted technologies

AI-assisted tools supported literature searches, mathematical drafting, source development, experiment orchestration, numerical interpretation and manuscript preparation. Numerical results are derived from archived execution outputs, and figures are produced by plotting code. The tool is not an author. The human authors must review and edit the manuscript and accept responsibility for its content before submission; this statement does not assert that their review has already occurred.

## Data availability

Code, synthetic data, frozen execution notebooks and numerical outputs are prepared as versioned reproducibility archives. A public version-specific repository link must be verified and inserted before submission.

'''
    if '## Declaration of generative AI and AI-assisted technologies' not in text:
        text=text.replace('## References',disclosures+'## References')
    archive_path=ROOT/'configs/public_archive.json'
    if archive_path.exists():
        archive=json.loads(archive_path.read_text())
        text=text.replace('Code, synthetic data, frozen execution notebooks and numerical outputs are prepared as versioned reproducibility archives. A public version-specific repository link must be verified and inserted before submission.',
            'Code, synthetic observations, numerical outputs and frozen execution notebooks are available in the public [reproducibility repository, version '+archive['version']+']('+archive['version_url']+'). The repository includes protocol configurations, source fingerprints, pinned author-code download specifications, compatibility adapters, postprocessing scripts and SHA256 manifests. Original author sources are downloaded from their pinned repositories rather than redistributed. No external experimental dataset is used. The tagged snapshot identifies the version supporting this article; a DOI has not been assigned.')
    p.write_text(text,encoding='utf-8')
    p=ROOT/'manuscript/supplement.md';text=p.read_text(encoding='utf-8')
    if '`configs/revision.json` was frozen' in text:
        begin=text.index('`configs/revision.json` was frozen')
        end=text.index('The new seed ranges',begin)
        text=text[:begin]+'''`configs/revision.json` was fixed before the corresponding numerical outputs; it is not a public preregistration. Numerical experiments execute in CPU Kaggle notebooks. Batched deterministic Van der Pol integration uses the recorded solver and tolerances; a correctness check compares batched and individual trajectories. Duplicate executions are not counted as independent evidence. Six synthetic residual arrays are used for covariance-error calibration. Notebook source fingerprints identify the executed implementation; local manuscript edits do not regenerate those execution snapshots.

'''+text[end:]
    text=text.replace('The new seed ranges','Seed ranges').replace('Earlier exploratory ranges are disjoint. ','Other experiment ranges are disjoint. ')
    text=text.replace('Seven original correctness tests','Seven correctness tests').replace('The five original completed notebooks remain immutable snapshots. ','Completed notebooks remain immutable snapshots. ')
    ai_start=text.index('AI-assisted tools supported') if 'AI-assisted tools supported' in text else -1
    if ai_start>=0:
        ai_end=text.index('## S5.',ai_start)
        text=text[:ai_start]+text[ai_end:]
    if '## S5. Exploratory negative diagnostic and replication' in text:
        start=text.index('## S5. Exploratory negative diagnostic and replication');end=text.index('## S6.',start)
        text=text[:start]+r'''## S5. Aggregation and approximation diagnostics

Coefficient-level aggregation diagnostics use the same 400 blocks, random seeds, initial conditions and sensor errors as Table 3. A separate execution preserves the coefficient matrices without changing the estimators. It replays the existing experiment and is not an additional set of independent records. Error norms are checked against the original output before combining these results.

For coefficient error vectors $e_m$ normalized by $\|\xi^*\|_F$, define $B=\|\bar e\|$, $S^2=n^{-1}\sum_m\|e_m-\bar e\|^2$ and $R^2=n^{-1}\sum_m\|e_m\|^2$. The empirical identity $R^2=B^2+S^2$ separates squared bias of the sample mean from dispersion; finite-sample $B^2$ includes sampling error of that mean. An adjusted estimate subtracts $\operatorname{tr}(\widehat{\operatorname{Cov}}e)/n$ using the unbiased sample covariance. Medians alone do not identify this decomposition.

The relative RMS linearization diagnostic in Table 2 divides the RMS difference between actual and linearized coefficient displacements by the RMS linear displacement. It can exceed one. All 2000 paths pass the singular-value guard at each amplitude, so guarded and complete empirical samples coincide; this is not proof of finite population moments.

'''+text[end:]
    p.write_text(text,encoding='utf-8')


if __name__=='__main__':prepare()
