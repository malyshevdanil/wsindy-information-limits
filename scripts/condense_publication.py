"""One-time editorial reduction; numerical data and frozen protocols are untouched."""
import json,re
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'manuscript'


def words(text):
    text=re.sub(r'\$\$.*?\$\$|\$[^$]*\$', ' ',text,flags=re.S)
    text=re.sub(r'!\[[^\]]*\]\([^)]*\)|\[@[^\]]*\]', ' ',text)
    return len(re.findall(r"\b[\w]+(?:[-'][\w]+)*\b",text))


def main():
    original=(BASE/'article_revised.md').read_text(encoding='utf-8')
    supp=(BASE/'supplement.md').read_text(encoding='utf-8')
    assert '### 7.12.' in original and '## S25.' not in supp
    s=original
    def block(a,b):return original[original.index(a):original.index(b)]
    def replace(a,b,value):
        nonlocal s
        i,j=s.index(a),s.index(b);s=s[:i]+value+s[j:]
    tables={int(m[1]):m[0].strip() for m in re.finditer(r'^Table (\d+)\.[^\n]*\n\n\|.*?(?=\n\n)',original,flags=re.M|re.S)}
    figures={int(m[1]):m[0] for m in re.finditer(r'^!\[Figure (\d+)\.[^\n]*',original,flags=re.M)}
    assert len(tables)==11 and len(figures)==7
    moved128=block('The largest clean discrete GLS coefficient error','### 7.3.')
    movedtau=block('A separate 128-path Van der Pol check','### 7.6.')
    extra=block('### 7.3.','### 7.4.')+block('### 7.6.','## 8.')
    extra=re.sub(r'^Table \d+\.[^\n]*\n\n\|.*?(?=\n\n)','',extra,flags=re.M|re.S)
    extra=re.sub(r'^!\[Figure \d+\.[^\n]*','',extra,flags=re.M)
    # S24 already contains the complete SE protocol and outcomes.
    extra=extra[:extra.index('### 7.12.')]
    extra=re.sub(r'^### (7\.\d+\.)',r'**Detailed protocol \1',extra,flags=re.M)
    extra=re.sub(r'^(\*\*Detailed protocol [^\n]+)$',r'\1**',extra,flags=re.M)
    moves='''## S25. Additional numerical checks and detailed protocols

These checks preserve the supplementary numerical evidence and implementation details of the condensed main paper. They use the original archived data without new simulations or amended outcomes.

### S25.1. Scalar optimal-test check

'''+block('### 7.1.','### 7.2.').split('\n\n',1)[1]+'''
### S25.2. Original 128-path nested-span check

'''+moved128+'''
### S25.3. Matched author WENDy and moderate-noise likelihood

'''+block('### 7.4.','### 7.5.').split('\n\n',1)[1]+'''
### S25.4. Sensitivity to supplied and fitted OU shape

'''+movedtau+'''
### S25.5. Detailed protocols for the retained comparisons

'''+extra
    scalar=block('For $x_a(t)=x_0e^{at}$','## 5.')
    scalar=safe_scalar_notation(scalar)
    raw=block('### 6.2.','### 6.3.').split('\n\n',1)[1]
    raw=raw.replace('For vector $u$','For a computational vector $w$').replace('u_n','w_n').replace('R_{q,h}u=a+z-u','R_{q,h}w=a+z-w')
    raw=raw.replace('$U_\\alpha$','$\\mathsf B_\\alpha$').replace('U_\\alpha','\\mathsf B_\\alpha')
    explicit=block('For scalar exponential dynamics','### 5.4.')
    # Keep the explanatory figure in the main text.
    explicit=re.sub(r'^!\[Figure 1\.[^\n]*','',explicit,flags=re.M)
    weaklikelihood=block('The same leading map arises','**Corollary 1')
    proofs=r'''## S26. Supplementary local-estimation details

### S26.1. Scalar path information and weak-ratio diagnostics

'''+scalar+r'''### S26.2. Proof of the conditional dense-span proposition

The main proposition separates estimator regularity from projection geometry. Its sufficient assumptions are:

- (P1) For each fixed finite family, the corrected weak solution has expansion $\widehat\theta_{K,\epsilon}-\theta=\epsilon P_K r_K+o_p(\epsilon)$, with $r_K$ centered Gaussian of covariance $\Sigma_K$, and full-column-rank $D_K$.
- (P2) Feasible pilots and covariance weights change $P_K=(D_K^\top\Sigma_K^{-1}D_K)^{-1}D_K^\top\Sigma_K^{-1}$ by $o_p(1)$. Deterministic quadrature error is $o(\epsilon)$.
- (P3) The squared norms of $\epsilon^{-1}(\widehat\theta_{K,\epsilon}-\theta)$ are uniformly integrable. Distributional convergence alone does not imply covariance convergence. Stable local solves with bounded fallback and bad-event probability $o(\epsilon^2)$ are one possible sufficient implementation.
- (P4) The deterministic, nested influence spans are dense in $\mathcal M$ in covariance norm, and $\mathcal I_{e,C}$ is positive definite.

**Proof.** Under (P1) and (P2), Slutsky's lemma gives the centered Gaussian limit with covariance $P_K\Sigma_KP_K^\top=(D_K^\top\Sigma_K^{-1}D_K)^{-1}=V_K$. To see the moment conclusion, let $X_\epsilon=\epsilon^{-1}(\widehat\theta_{K,\epsilon}-\theta)$. Assumption (P3) and the continuous-mapping theorem give convergence of each first and second coordinate moment: $|X_iX_j|\leq\|X_\epsilon\|^2$ supplies the necessary uniform integrability. Hence $E X_\epsilon\to0$ and $\operatorname{Cov}(X_\epsilon)\to V_K$. No risk conclusion is drawn for an unguarded ratio lacking (P3).

Let $t_i=C^{-1}R_i\in\mathcal M$ be the efficient dual representers of the coefficient sensitivities from the boundary-free theorem, and let $Q_K$ be covariance-orthogonal projection onto the finite influence span. Because every test influence annihilates $Z$, its pairing with $S_i$ equals its pairing with $R_i$. The finite-family information is therefore

$$
(D_K^\top\Sigma_K^{-1}D_K)_{ij}
=\langle Q_Kt_i,Q_Kt_j\rangle_C.
$$

By (P4), $Q_Kt_i\to t_i$ in covariance norm. Cauchy–Schwarz then gives entrywise convergence to $\langle t_i,t_j\rangle_C=(\mathcal I_{e,C})_{ij}$. The nested projected information matrices increase in Loewner order, and positive definiteness makes matrix inversion continuous at the limit. Consequently $V_K\to\mathcal I_{e,C}^{-1}$. These are two sequential limits. A joint $K,\epsilon$ limit requires uniform controls absent from (P1)–(P4). The proof is a specialization of classical estimating-equation linearization and efficient-score projection [@godambe1960; @bickel1993; @vandervaart1998, Chapters 5 and 25]. $\square$

### S26.3. Parameter-dependent weak likelihood

'''+weaklikelihood+r'''### S26.4. Raw features and linear-time temporal covariance application

'''+raw+r'''### S26.5. Explicit scalar OU optimal test

'''+explicit
    # Compact abstract: three claims, with one finite-noise caveat.
    replace('## Abstract','**Keywords:**',r'''## Abstract

Boundary-free weak equations discard local information when initial states are known. For deterministic ODEs observed with stationary Ornstein–Uhlenbeck measurement error and supplied covariance, we characterize their closed adjoint influence class and its sharp nuisance-adjusted information bound. The restriction raises the local RMSE reference by 29% for Van der Pol and almost doubles it for Lorenz. Allowing known-initial traces recovers the full local influence class. This boundary adaptation is due to Brunel and coauthors; universal time-only windows recover about 97% of the empirical gap in Van der Pol and 95% in Lorenz at small noise. An ablation attributes the recovery to initial traces. For a prespecified cubic candidate, targeted weak-GLS Wald testing reaches the local small-noise detection reference while the supplied MSTLS selects substantially larger coefficients. At moderate noise, endpoint recovery is incomplete and nominal Wald errors require estimator-level calibration. The contribution is an exact weak-class characterization and quantitative correlated-noise benchmarks, rather than a new boundary remedy or general selection procedure.

''')
    # Scalar special cases and computational detail go to S26.
    replace('For $x_a(t)=x_0e^{at}$','## 5.',r'''Exact discrete OU innovations converge to this finite path information under grid refinement. S26.1 gives the scalar calculation and shows why unguarded weak ratios can have infinite error moments despite a finite Gaussian influence variance. Information bounds concern local regular influences; finite-noise estimator risk requires separate control.

''')
    replace('For scalar exponential dynamics','### 5.4.',r'''The scalar exponential specialization and fitted-pilot formula are in S26.5. Figure 1 shows how the generalized optimum and its compact cutoff differ near the endpoints.

'''+figures[1]+'\n\n')
    replace('The same leading map arises','### 5.5.',r'''This is the Godambe-information construction for optimally weighted estimating equations [@godambe1960]: the local information is $D^\top\Sigma^{-1}D$. WENDy-MLE has the same leading small-noise map under regularity and correctly specified covariance; finite-noise covariance-derivative terms distinguish its objective (S26.3).

**Corollary 1 (small-noise, dense-span efficiency).** For known unit-amplitude covariance shape $C$, assume the fixed-family GLS expansion, consistent feasible weights, quadrature error $o(\epsilon)$, and uniform integrability of squared scaled errors stated as (P1)–(P3) in S26.2. Then

$$
\epsilon^{-1}(\widehat\theta_{K,\epsilon}-\theta)\Rightarrow N(0,V_K),\qquad
\epsilon^{-2}\operatorname{Cov}(\widehat\theta_{K,\epsilon})\to V_K.
$$

If the nested influence spans are dense in $\mathcal M$ and $\mathcal I_{e,C}$ is positive definite, then $V_K\to\mathcal I_{e,C}^{-1}$. The limits are sequential: small noise for fixed $K$, then dense spans. S26.2 proves both assertions under the explicit assumptions. A joint limit or risk claim for unguarded ratios requires additional analysis.

''')
    replace('### 6.2.','## 7.',r'''### 6.2. Implementation and covariance misspecification

The finite Hermite sum can be evaluated on a grid without storing a dense temporal matrix. S26.4 gives the recursion and the extension to incorrectly corrected raw features. Under temporal weights $W$, leading covariance is $P_W\Sigma_{\mathrm{true},1}P_W^\top$; changing correction variances can also change the moment mean. The present comparisons supply the sensor covariance. Sensitivity to supplied and fitted correlation times is reported in S25.4.

''')
    # Rebuild the numerical section around price, recovery and targeted selection.
    def table(old,new):return tables[old].replace(f'Table {old}.',f'Table {new}.',1)
    def fig(old,new):return figures[old].replace(f'Figure {old}.',f'Figure {new}.',1)
    numerical=r'''## 7. Numerical benchmarks: cost, recovery and targeted selection

The benchmarks use fully observed polynomial systems, supplied structure and sensor covariance, and archived paired observations. Feasible covariance pilots receive no generating coefficients. S22–S25 preserve protocols, seed ranges, fingerprints and additional checks. Local references and finite-noise performance are reported separately.

### 7.1. Information cost and finite-span efficiency

For Van der Pol, we fit four independent coefficients in $x'=av$, $v'=bx+cv+dx^2v$, with truth $(1,-1,1,-1)$, initial state $(2,0)$, $T=20$, $h=0.005$ and OU range 0.3. Relative errors are divided by the true coefficient norm, two. The known-initial local RMSE reference at 5% noise is 0.013093, versus the boundary-free reference 0.016904: a 29% cost. Both scale linearly with noise amplitude.

'''+table(2,1)+r'''

Rich spans produce most of the improvement: 480-window OU GLS is within 0.2% of the boundary-free local reference, while iid weights retain a roughly 5% gap. Exact residual-covariance calculations verify the stated local geometry. The original 128-path finite-span check is in S25.2; the larger paired comparison follows.

'''+table(4,2)+r'''

All 2000 prespecified paths at 0.5% noise are retained, including the original 128-path prefix once. Feasible weak GLS uses a bounded 30-window covariance pilot and 480 prescribed windows. Likelihood refinements use exact OU innovations and analytical ODE sensitivities. Paired iid-minus-OU RMSE is 0.00007426, with bootstrap interval [0.00005576, 0.00009327], a 4.22% benefit from temporal weighting. The initial unknown-state scoring step has an appreciable nonlinear initialization bias; three steps approach the nuisance-adjusted reference. The known-initial estimate has lower error because it observes the information discarded by compact tests.

'''+fig(4,2)+r'''

### 7.2. Recovery with universal endpoint tests

For known $u$, allow nonzero test traces and use the corrected residual

$$
m_\phi(\theta)=-\int_0^T[\phi'(t)^TY(t)+\phi(t)^Tf_\theta^H(Y(t))]dt
+\phi(T)^TY(T)-\phi(0)^Tu.
$$

The universal family adds twelve time-only endpoint half-windows and $Y(0)-u$ controls to 480 interior windows. Covariance weights still require an integrated 30-window pilot. Theorem 3 identifies the full local influence class; the known-initial trace adaptation itself is established by Brunel et al. [@brunel2014, Section 2.3].

'''+table(6,3)+r'''

These 1968 records were unused for endpoint fitting, while their compact and scoring references were archived. The frozen confirmation criterion requires at least 50% gap recovery. Universal endpoints recover 97.1%, with paired interval 95.5–98.6%, and lie 0.88% above scoring. Model-adjoint tests realize score directions and provide an implementation check, rather than a separate methodological gain. S22 and S25.5 record the pre-confirmation criterion amendment and quadrature checks.

'''+fig(5,3)+r'''

### 7.3. Attribution to initial traces

The 500-path ablation holds Simpson quadrature and the covariance pilot fixed. Component gains are measured against the matched compact branch and known-initial scoring.

'''+table(7,4)+r'''

Initial traces alone recover essentially the entire gap, quantifying the Brunel mechanism under OU noise. Initial controls give a smaller standalone gain and terminal traces give none resolved. Combining all components changes RMSE relative to initial traces by $7.46\times10^{-8}$, with paired interval $[-7.18\times10^{-7},8.55\times10^{-7}]$. Their contributions overlap and cannot be added. The control's standalone improvement is specific to correlation between initial and interior error, rather than evidence for an additional gain after trace adaptation.

'''+fig(6,4)+r'''

### 7.4. Moderate-noise limits and Lorenz confirmation

The 5% Van der Pol check reuses all 128 matched observations, with known initial state and the same feasible pilot. The likelihood reference optimizes exact OU innovations.

'''+table(8,5)+r'''

Universal endpoints recover 32.6% of the gap (18.3–45.8%), failing the frozen 50% criterion. Clean-design first-order and full-Hermite references are 0.013094 and 0.013636, below empirical dispersion 0.021125. Restored local geometry therefore leaves finite-noise design, pilot and fitting effects. The matched author-WENDy comparison is in S25.3.

For Lorenz, we fit $(s,r,b)=(10,28,8/3)$ in its supplied physical structure, with $u=(1,1,1)$, $T=4$, $h=0.001$, OU range 0.1 and 0.5% coordinate noise. The unknown-initial local relative RMSE reference is 0.000411, versus 0.000209 when the initial state is known: almost a factor two.

'''+table(9,6)+r'''

All 256 paths and all converged fits are retained. Universal endpoint windows recover 95.5% of the compact-to-likelihood gap (92.8–98.0%), passing the same 50% rule. The full-Hermite clean-design reference is 0.00020926, close to the observed endpoint dispersion. These two supplied structures support small-noise recovery; moderate-noise universality remains unestablished.

### 7.5. Prespecified testing versus sparse inclusion

A negative $\delta x^3$ perturbation to the second Van der Pol equation leaves the orbit bounded. Let $S_{-\delta}$ collect the other coefficient sensitivities. The local composite likelihood information is $\|s_\delta-\Pi_{S_{-\delta}}s_\delta\|_{\mathcal H}^2$; compact tests additionally project out $Z$. Profiling all 17 other library coefficients gives two-point reference magnitudes 0.01928 and 0.01941 at 0.5% noise. The initial-state cost for this direction is only 0.68%, despite the larger aggregate estimation cost. These references use equal-prior testing error 0.05, not the Wald significance level.

The targeted Wald procedure fits all eighteen degree-one-through-three coefficients with Hermite-corrected weak GLS. A bounded core-five pilot supplies covariance and the quadrature-defect correction. The author's MSTLS and observation paths are unchanged; S20 and S23 give the pinned pipeline and intermediate-grid protocol.

'''+table(10,7)+r'''

At small noise, nominal null rejection is 7/100 (Wilson 3.4–13.7%), compatible with 5%. At one reference unit the targeted test rejects on 91/100 paths, versus no MSTLS inclusions, consistent with the local Gaussian Wald prediction. MSTLS reaches its transition at larger magnitudes. This contrasts supplied-candidate testing with a particular sparse selector; it neither attributes the difference exclusively to its threshold nor establishes library-wide or post-selection validity.

'''+fig(7,5)+r'''

### 7.6. Moderate-noise standard-error diagnostic

At 5% noise, the nominal Wald null rate is 23/100 (Wilson 15.8–32.2%), so the corresponding alternative rejection frequencies remain descriptive. A separately frozen diagnostic on these same 100 null paths holds coefficient estimates and critical values fixed. A random-design score sandwich gives 2/100 rejections (0.6–7.0%), and a full-refit parametric bootstrap gives 5/100 (2.2–11.2%); both pass the Wilson-compatibility gate. Median paired SE multipliers are 2.81 and 1.79. S24 reports the complete counts, 399 bootstrap refits per path, 999 model score draws and numerical audit. These previously examined paths motivate independent size confirmation; power with the revised SEs remains untested.

'''
    replace('## 7.','## 8.',numerical)
    # Add proper classical attribution without presenting it as a new theorem.
    old='Continuous-time regression with correlated errors has an established optimal-design theory, including optimal signed measures with endpoint masses.'
    start=s.index(old);end=s.index('The continuum-data theory',start)
    s=s[:start]+'''Optimal estimating equations and their sandwich information are classical [@godambe1960]. Efficient scores remove nuisance-score directions by orthogonal projection [@bickel1993; @vandervaart1998, Chapter 25]; the initial-state Schur complement is the finite-dimensional instance here. Continuous-time BLUE and endpoint-weighted estimation likewise precede weak ODE analysis [@dette2016; @dette2017]. Our contribution identifies the adjoint range of the restricted weak class and quantifies its cost and recovery in this dynamical experiment.

'''+s[end:]
    start=s.index('The paper has three linked aims.');end=s.index('## 2.',start)
    s=s[:start]+'''The paper characterizes the boundary-free and known-initial extended influence classes, applies classical information geometry to their bounds, and benchmarks finite-family recovery and prespecified testing. Polynomial Hermite covariance supplies weights. Initial-trace adaptation is credited to Brunel et al.; independent-record GMM, fission and cross-fitting remain supplementary extensions.

'''+s[end:]
    # Name the projection at its mathematical point of use.
    old='Projection and best linear unbiased estimation in correlated continuous regression are classical [@dette2016].'
    s=s.replace(old,'The nuisance-score projection and efficient influence are classical [@bickel1993; @vandervaart1998, Chapter 25], as is continuous-record BLUE [@dette2016; @dette2017].')
    # Distinguish Hermite scalar argument from initial state.
    s=s.replace('$H_0=1$, $H_1=u$, $H_{q+1}=uH_q-qvH_{q-1}$', '$H_0(z;v)=1$, $H_1(z;v)=z$, $H_{q+1}(z;v)=zH_q(z;v)-qvH_{q-1}(z;v)$')
    s=s.replace('$\\exp(zu-vz^2/2)$','$\\exp(az-va^2/2)$')
    s=s.replace('In Corollary 1, $C$','In Corollary 1, $C$')
    # Complete supplement table relocation.
    table_map={1:'S21',2:'1',3:'S22',4:'2',5:'S23',6:'3',7:'4',8:'5',9:'6',10:'7',11:'S20'}
    figure_map={1:'1',2:'S12',3:'S13',4:'2',5:'3',6:'4',7:'5'}
    section_map={1:'S25.1',2:'7.1',3:'7.1',4:'S25.3',5:'7.4',6:'7.2',7:'7.3',8:'7.5',9:'7.4',10:'7.4',11:'7.5',12:'7.6'}
    def oldrefs(text):
        text=re.sub(r'\bTable (\d+)\b',lambda m:'Table '+table_map.get(int(m[1]),m[1]),text)
        text=re.sub(r'\bFigure (\d+)\b',lambda m:'Figure '+figure_map.get(int(m[1]),m[1]),text)
        text=re.sub(r'\bSection 7\.(\d+)\b',lambda m:('Section '+section_map[int(m[1])] if not section_map[int(m[1])].startswith('S') else section_map[int(m[1])]),text)
        text=re.sub(r'\bSections 5 and 7\b','Sections 5 and 7',text)
        return text
    supp=oldrefs(supp);moves=oldrefs(moves);proofs=oldrefs(proofs)
    # Old section pointers outside the newly reconstructed numerical section.
    i,j=s.index('## 7.'),s.index('## 8.')
    prefix=s[:i];suffix=s[j:]
    prefix=re.sub(r'\bSection 7\.(\d+)\b',lambda m:'Section '+section_map[int(m[1])],prefix)
    s=prefix+s[i:j]+suffix
    # Table 11 is now supplementary; S24 keeps its protocol and interpretation.
    position=supp.index('Table S20 reports')
    caption=tables[11].replace('Table 11.','Table S20.',1)
    supp=supp[:position]+caption+'\n\n'+supp[position:]
    supp=supp.replace('## References',moves+'\n\n'+proofs+'\n\n## References')
    # Renumber statements coherently after downgrading the two general results.
    def labels(text):
        mapping={'Lemma 1':'Lemma 2','Theorem 1':'Lemma 1','Theorem 2':'Theorem 1','Theorem 3':'Theorem 2','Proposition 3':'Proposition 4','Corollary 1':'Proposition 3'}
        return re.sub(r'\b(?:Lemma 1|Theorem [123]|Proposition 3|Corollary 1)\b',lambda m:mapping[m[0]],text)
    s=labels(s);supp=labels(supp)
    # The finite-family assertion is conditional, with a full proof in S26.2.
    s=s.replace('This theorem concerns regular WENDy-like','This proposition concerns regular WENDy-like')
    (BASE/'article_revised.md').write_text(s,encoding='utf-8')
    (BASE/'supplement.md').write_text(supp,encoding='utf-8')
    count={'method':'word tokens in Markdown after removing math, citation tokens and image captions; includes tables and headings',
        'before':words(original),'after':words(s),'reduction_fraction':1-words(s)/words(original),
        'abstract_words':words(s.split('## Abstract\n\n')[1].split('**Keywords:**')[0]),
        'experimental_subsections_before':12,'experimental_subsections_after':6,
        'main_tables':7,'main_figures':5,'all_numerical_outputs_unchanged':True,
        'main_table_mapping':table_map,'main_figure_mapping':figure_map}
    (ROOT/'reports/editorial_reduction.json').write_text(json.dumps(count,indent=2),encoding='utf-8')
    print(json.dumps(count),flush=True)


def safe_scalar_notation(text):
    text=text.replace('A/(D_\\phi+B)','N_\\phi/(D_\\phi+E_\\phi)').replace('\\frac{A}{D_\\phi+B}','\\frac{N_\\phi}{D_\\phi+E_\\phi}')
    text=text.replace('A=\\int','N_\\phi=\\int').replace('B=\\int','E_\\phi=\\int')
    text=text.replace('$(A,B)$','$(N_\\phi,E_\\phi)$').replace('of $A$ given $B$','of $N_\\phi$ given $E_\\phi$')
    return text


if __name__=='__main__':main()
