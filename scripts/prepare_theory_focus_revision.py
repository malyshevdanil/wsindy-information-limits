"""Focus the manuscript on class geometry and constructive attainability."""
import re
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'manuscript'

def main():
    p=BASE/'article_revised.md';s=p.read_text(encoding='utf-8')
    abstract='''Boundary-free weak equations discard information about known initial states under correlated measurement error. For deterministic ODEs observed with stationary Ornstein–Uhlenbeck noise, we identify their closed adjoint influence class and sharp nuisance-adjusted information bound. The restriction raises the local RMSE reference by 29% for Van der Pol and almost doubles it for Lorenz. For polynomial libraries satisfying a nonzero trajectory-jet condition, we construct a feasible full-library estimator with smoothed-data instruments, vanishing regularization and bounded output. Pilot rank, covariance-weight consistency and uniform integrability are derived; the estimator attains finite-span covariance and reaches the weak-class bound in sequential small-noise and dense-span limits. Exact symbolic calculation verifies the condition for the eighteen-coefficient Van der Pol library. Known-initial traces restore the full closed influence class. This established adaptation, due to Brunel and coauthors, recovers about 97% and 95% of the empirical small-noise gaps in Van der Pol and Lorenz. Matched covariance benchmarks support constructive attainability while showing excess risk and incomplete boundary recovery at moderate noise.'''
    a=s.index('## Abstract\n\n')+len('## Abstract\n\n');b=s.index('**Keywords:**',a);s=s[:a]+abstract+'\n\n'+s[b:]
    a=s.index('**Proposition 3 (');b=s.index('### 5.5.',a)
    prop='''**Proposition 3 (guarded fixed-family weak GLS).** Fix a finite deterministic boundary-free test family, a correctly specified polynomial vector field linear in its coefficients, known initial state and supplied independent-coordinate OU covariance shape. Suppose the true coefficient is interior to a compact final box; clean and nearby pilot trajectories exist in a common bounded region; the clean pilot and final designs have full column rank; and the first-order residual covariance $\\Sigma_1$ is positive definite. Quadrature weights have bounded total mass, converge to their integrals, and clean quadrature and pilot-model numerical defects are $o(\\epsilon)$ with the local derivative control in (G4) of S27.2.

Use corrected least squares for the pilot, clip it to its supplied box, evaluate the full Hermite covariance on its integrated trajectory, and apply the vanishing-ridge, guarded and box-projected GLS routine described below. Then its pilot is root-$\\epsilon$ consistent, its normalized weights converge, and

$$
\\epsilon^{-1}(\\widehat\\theta_{K,\\epsilon}-\\theta)\\Rightarrow N(0,V_K),\\qquad
\\epsilon^{-2}\\operatorname{Cov}(\\widehat\\theta_{K,\\epsilon})\\to V_K.
$$

The squared scaled errors are uniformly integrable. These conclusions are derived from the stated structural and numerical conditions; an asymptotic expansion, weight consistency and integrability are not assumed. S27.2 proves the proposition. For nested dense families satisfying these conditions separately, the subsequent family limit gives the bound of Theorem 1. Theorem 3 supplies a full-library smoothed-data construction whose clean ranks follow from a jet condition instead of being assumed.

The archived solvers in Tables 1–6 use fixed relative ridge and an unbounded final estimate. Their local matrix geometry and empirical recovery are separate from the risk guarantee above. For fixed design $D$ and residual shape $\\Sigma_1$, their regularized linear map has

$$
W_\\rho=[\\Sigma_1+\\rho\\operatorname{tr}(\\Sigma_1)I/q]^{-1},\\quad
P_\\rho=(D^\\top W_\\rho D)^{-1}D^\\top W_\\rho,\\quad
V_\\rho=P_\\rho\\Sigma_1P_\\rho^\\top.
$$

Here $q$ is the number of moments, and $V_0$ is ideal GLS covariance. The exact identity

$$
V_\\rho-V_0=(P_\\rho-P_0)\\Sigma_1(P_\\rho-P_0)^\\top\\succeq0
$$

quantifies ridge cost without any estimator-limit assumption. S27.3 evaluates it on the actual archived spans. This matrix check cannot establish moments of an unbounded nonlinear final estimate. Table S27 separates those archived benchmarks, the constructive estimator and supplied-candidate testing.

'''
    s=s[:a]+prop+s[b:]
    a=s.index('**Theorem 3 (');algorithm='''The following procedure estimates every supplied library coefficient. Let $B$ be the compact coefficient box and use uniform grids with an even number of intervals. Smooth the linear interpolant of $Y$ with

$$
K_b(t)=b^{-1}(315/256)[1-(t/b)^2]_+^4,\\quad b=\\epsilon^{1/10}/4,
\\qquad w(t)=[1-((t-T/2)/(T/2-1/2))^2]_+^4.
$$

Reflect the interpolant at the endpoints. The proxy $z$ supplies the $m$ feature tests $w g_j(z)$ and their analytic derivatives. Their Hermite-corrected weak equations give a full-library least-squares pilot, clipped to $B$. A fixed deterministic-window pilot on a declared correct structural family gives the alternative construction in Proposition 3.

Append the cubic splines $B_3(t/a-j)$, $a=T/2^k$, $j=2,\\ldots,2^k-2$, where $B_3$ is the normalized cardinal cubic spline supported on $[-2,2]$. Whiten the combined $H^1$ Gram matrix, retaining eigenvalues above $10^{-10}\\epsilon^{1/10}$ times its largest eigenvalue. The exact Hermite covariance (Proposition 4), evaluated at the proxy and pilot, gives $\\widehat\\Sigma_1=\\widehat\\Omega/\\epsilon^2$. Add ridge $\\epsilon^2\\operatorname{tr}(\\widehat\\Sigma_1)I/q$, then solve corrected GLS and project the result to $B$. Relative minimum-eigenvalue guards use $10^{-8}\\epsilon^2$; singular, nonfinite or failed calculations return the box midpoint. S28 specifies the column scaling and covariance guard. No generating support or pilot ODE solution is supplied to this full-library procedure.

The technical restriction $T>1$ makes the fixed half-unit interior margins of $w$ nonempty. Proportional margins and a rescaled kernel give an equivalent construction on any positive horizon. Write $\\mathcal I_{e,C}$ for the initial-state-adjusted information of Theorem 1 computed with the unit-amplitude covariance $C$ of $\\eta$.

'''
    s=s[:a]+algorithm+s[a:]
    s=s.replace('Use the full-library smoothed-data pilot, Hermite-corrected moments, cubic-spline tests, vanishing relative ridge and bounded final estimate specified in S28.0, with $h^2=o(\\epsilon)$.','For the full-library procedure just defined, let $h^2=o(\\epsilon)$.')
    s=s.replace('The formulas and proof are in exact arithmetic. ','').replace('\\mathcal I_{\\theta\\mid u}','\\mathcal I_{e,C}')
    start=s.index('### 7.5. Prespecified');end=s.index('### 7.7.',start)
    diagnostic=s[start:end];figure=re.search(r'!\[Figure 5\..*?\]\([^\n]+\)',diagnostic,re.S)
    brief='''### 7.5. Direction-specific information cost

For the prespecified negative cubic perturbation $\\delta x^3$ in the second Van der Pol equation, profiling all seventeen other library coefficients gives local full-information and boundary-free two-point reference magnitudes 0.01928 and 0.01941 at 0.5% noise. The initial-state information cost for this direction is only 0.68%, despite its larger aggregate coefficient-estimation cost. Supplied-candidate Wald testing and library-wide sparse selection have different inputs and targets; their detailed comparison and calibration studies are supplementary (S23, S24 and S29).

'''
    s=s[:start]+brief+s[end:];s=s.replace('### 7.7.','### 7.6.').replace('Table 8.','Table 7.')
    s=s.replace('This supports the fixed-family small-noise covariance prediction and identifies finite-noise excess risk.','This supports the fixed-family small-noise covariance prediction and identifies finite-noise excess risk. At 5% noise the total eighteen-coefficient RMSE is 0.5813, or root-mean-square coordinate error 0.1370. Full-library precision is therefore limited at this amplitude even near its information reference. Each level has only 100 paths; the intervals leave the apparent lower-noise trend unresolved.')
    a=s.index('Information loss is direction specific.');b=s.index('\n\n',a)
    s=s[:a]+'''Information loss is direction specific: initial-state projection raises the local cubic-candidate detection reference by only 0.68%, compared with the larger aggregate coefficient-estimation cost. The supplementary testing study distinguishes prespecified candidates from unrestricted structure discovery. The main empirical attainability check instead matches the full-library estimator to its own finite-span covariance, with uncertainty reported for the limited 100-path samples. At moderate noise this large library remains imprecise, and the local theorem does not establish finite-noise inferential calibration.'''+s[b:]
    p.write_text(s,encoding='utf-8')
    p=BASE/'supplement.md';x=p.read_text(encoding='utf-8');x=x.replace('### S26.2. Proof of Proposition 3','### S26.2. A conditional GLS consequence')
    x=x.replace('Theorem 3 does not transfer its','Theorem 3 does not transfer its').replace('Five-level Table 8','Five-level Table 7').replace('For Table 8,','For Table 7,')
    x=x.replace('Original Table 7; independent size/power Table S26; incomplete 5% follow-up','Original Table S19; independent size/power Table S26; incomplete 5% follow-up')
    if figure:
        caption=figure[0].replace('Figure 5.','Figure S14.');anchor='### S23.2.'
        pos=x.index(anchor);x=x[:pos]+caption+'\n\n'+x[pos:]
    x=x.replace('with Table 7','with the original candidate comparison').replace('as Table 7','as the original candidate comparison')
    x=x.replace('None\nof the reported Monte Carlo outcomes uses it.','The fixed-ridge estimates in main Tables 1–6 use a different solver. The original_guarded rows of the five-level geometry study exercise this variant; the main full-library covariance check uses the distinct smoothed-data procedure of Theorem 3.')
    x=x.replace('**Proposition S27.1 (guarded fixed-family weak GLS).**','**Proposition 3 (guarded fixed-family weak GLS; detailed assumptions).**')
    x=x.replace('\\mathcal I_{\\theta\\mid u}','\\mathcal I_{e,C}').replace('The mathematical statement treats its formulas in exact arithmetic; the\nindependent numerical study checks their floating-point implementation.','The independent numerical study checks the specified implementation.')
    p.write_text(x,encoding='utf-8')
    p=BASE/'self_instrumented_gls_theory.md';x=p.read_text(encoding='utf-8').replace('\\mathcal I_{\\theta\\mid u}','\\mathcal I_{e,C}')
    x=x.replace('The mathematical statement treats its formulas in exact arithmetic; the\nindependent numerical study checks their floating-point implementation.','The independent numerical study checks the specified implementation.')
    x=x.replace('Take $h^2=o(\\epsilon)$.','Use even uniform grids and take $h^2=o(\\epsilon)$. Write $C$ for the unit-amplitude covariance of $\\eta$ and $\\mathcal I_{e,C}$ for the information of Theorem 1 computed with that shape.')
    p.write_text(x,encoding='utf-8')

if __name__=='__main__':main()
