"""Condense the current submission, preserving removed research in an archive."""
import re
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'manuscript'

def main():
    p=BASE/'article_revised.md';s=p.read_text(encoding='utf-8')
    related='''## 2. Related estimation methods

ODE estimation has long used trajectory smoothing and integral matching: spline least squares [@varah1982], nonparametric gradient matching [@brunel2008], measurement-error regression [@liang2008], smooth-and-match estimation [@gugushvili2012], and direct integral estimation [@dattner2015]. Generalized profiling can be efficient under its approximation and independent-sampling assumptions [@ramsay2007; @qi2010]; one-step likelihood refinement of an ODE pilot is also established [@dattner2018]. These provide full-observation references. Our question is the information attainable by restricted weak equations under physical temporal correlation.

WENDy treats weak parameter estimation as an errors-in-variables problem, retaining coupled feature and response perturbations and iteratively updating their linearized covariance [@bortz2023]. WENDy-MLE uses an approximate Gaussian weak-residual likelihood with parameter-dependent covariance and analytic derivatives [@rummel2025]. Its log determinant and covariance-derivative terms distinguish its finite-noise objective from frozen-weight GLS. Polynomial Hermite correction and exact residual covariance here extend that coupled calculation to all Gaussian orders for supplied OU noise. They do not alone correct the nonlinear coefficient estimator. Gaussian correction of noisy polynomial features has classical precedents [@fuller1987; @stefanski1989], and WENDy finite-noise bias and uncertainty remain active topics [@chawla2026].

Test design also has substantial precedents. Tran and Bortz choose local radii through quadrature-error estimates [@tran2026], while Fourier Weak SINDy selects sinusoidal frequencies spectrally [@chen2026]. WmSINDy combines weak regression with dynamics and noise reconstruction [@lopez2026]. Ensemble-SINDy quantifies feature inclusion by bootstrap aggregation [@fasel2022]. Our author-method checks preserve their distinct structural inputs and are documented in the research archive. The new class characterization and constructive full-library result complement these methods; the supplied-candidate Wald comparison establishes no general superiority for structure selection.

### 2.1. Complementarity with orthogonality-condition estimation

Brunel et al. [@brunel2014] combine a nonparametric proxy, orthogonality conditions and GMM in an independent-observation, increasing-point experiment, including delay equations. They recognize boundary-information loss and adapt moments to known boundary values. Our deterministic polynomial setting is narrower. It instead treats fixed-horizon physical OU error, identifies the complete adjoint influence classes, gives their sharp initial-state Schur-complement relation and an explicit OU optimal test, and constructs an attainable full-library estimator under a checkable jet condition. The small-noise and sequential dense-span limits differ from increasing-point asymptotics. The results complement each other.

'''
    s=s[:s.index('## 2. Related estimation methods')]+related+s[s.index('## 3. Observation model'):]
    s=s.replace('Let $\\mathcal Q_h$ be composite trapezoidal quadrature with weights $w_n$.',
        'Let $\\mathcal Q_h$ denote composite trapezoidal quadrature with weights $w_n$ in the formulas; matched endpoint confirmations use composite Simpson quadrature as specified in their frozen protocols.')
    s=s.replace('S20 and S23 give the pinned pipeline and intermediate-grid protocol.',
        'Research archive S20 and supplement S23 give the pinned pipeline and intermediate-grid protocol.')
    s=re.sub(r'\n## Declaration of generative AI and AI-assisted technologies\n.*?(?=\n## Data availability)', '\n',s,flags=re.S)
    p.write_text(s,encoding='utf-8')
    p=BASE/'supplement.md';full=p.read_text(encoding='utf-8')
    theory=(BASE/'self_instrumented_gls_theory.md').read_text(encoding='utf-8').split('## S29.',1)[0].strip()
    start=full.index('## S28.');end=full.index('## S29.')
    full=full[:start]+theory+'\n\n'+full[end:]
    parts=re.split(r'(?m)(?=^## S\d+\.)',full);prefix=parts[0];keep=[];removed=[]
    for part in parts[1:]:
        number=int(re.match(r'## S(\d+)\.',part)[1])
        (keep if number==1 or number>=22 else removed).append(part)
    archive=BASE/'research_extensions_archive.md'
    if removed and not archive.exists():
        head=re.sub(r'title:.*', 'title: "Research archive: supplementary extensions and historical benchmarks"',prefix)
        archive.write_text(head+'This archive preserves earlier research extensions and supporting benchmarks. Its S-labels retain the original cross-references. These extensions are outside the submitted mathematical claims.\n\n'+''.join(removed)+'\n\n## References\n',encoding='utf-8')
    intro='The submission supplement contains the proofs and checks supporting the current paper. Historical S-labels are retained for reproducibility. Sections S2–S21, including independent-record GMM, fission, cross-fitting and incidental-state studies, are preserved separately in [Research archive](research_extensions_archive.pdf). References to those omitted labels below point to that archive.\n\n'
    p.write_text(prefix+intro+''.join(keep),encoding='utf-8')

if __name__=='__main__':main()
