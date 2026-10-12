"""Add the completed matched moderate-noise likelihood reference without a new table."""
import json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'manuscript'


def main():
    data=json.loads((ROOT/'reports/moderate_likelihood_results.json').read_text());known=next(r for r in data['groups'] if r['method']=='known_initial');unknown=next(r for r in data['groups'] if r['method']=='unknown_initial')
    p=BASE/'article_revised.md';text=p.read_text(encoding='utf-8');at=text.index('### 8.5.')
    paragraph=f"A separate bounded local Gaussian path-likelihood optimization at 5% uses these identical 128 observations and exact stationary OU innovations. With the initial state known, relative RMSE is {known['relative_rmse']:.5f} (95% path-bootstrap interval [{known['rmse_bootstrap_95'][0]:.5f}, {known['rmse_bootstrap_95'][1]:.5f}]); jointly fitting its two coordinates gives {unknown['relative_rmse']:.5f} ([{unknown['rmse_bootstrap_95'][0]:.5f}, {unknown['rmse_bootstrap_95'][1]:.5f}]). Analytical ODE sensitivities supply the Jacobian. All fits meet the optimizer's stopping conditions, and {known['bound_active_paths']+unknown['bound_active_paths']} activate a declared parameter bound. These local likelihood references quantify full-path accuracy at moderate noise separately from the weak component (S19).\n\n"
    text=text[:at]+paragraph+text[at:];p.write_text(text,encoding='utf-8')
    p=BASE/'supplement.md';supp=p.read_text(encoding='utf-8');s='''## S19. Matched moderate-noise path likelihood

The frozen `moderate_likelihood.json` protocol fits 128 known-initial and 128 unknown-initial paths at 5%, matching the observation SHA256 values in the full author WENDy comparison. Each fit minimizes the exact stationary OU innovation sum of squares at known covariance. DOP853 integrates state and coefficient/initial sensitivities; analytic Jacobians drive bounded nonlinear least squares, with all three stopping tolerances 1e-9 and a maximum 100 function evaluations. Initial coefficient estimates are the feasible OU weak coefficients already archived; unknown initial coordinates start from the noisy first observation. Bounds are [0.5,1.5], [-1.5,-0.5], [0.5,1.5], [-1.5,-0.5] for the four coefficients and [-3,3] for each initial coordinate. The fit is a local likelihood optimization within declared bounds. The optimizer success flag and all active-bound masks, evaluation counts and optimality statistics are retained for every path.

'''
    for r in data['groups']:
        s+=f"{r['method'].replace('_',' ').capitalize()}: RMSE {r['relative_rmse']:.6f}, bias {r['relative_bias']:.6f}, dispersion {r['relative_dispersion']:.6f}, maximum function evaluations {r['maximum_nfev']}, active-bound paths {r['bound_active_paths']}.\n\n"
    for r in data['paired_author_comparisons']:
        s+=f"Author WENDy minus {r['method'].replace('_',' ')} likelihood RMSE: {r['author_minus_likelihood_rmse']:.6f}, paired 95% interval [{r['bootstrap_95'][0]:.6f}, {r['bootstrap_95'][1]:.6f}].\n\n"
    s+='Intervals use 5000 paired path-bootstrap resamples with seed 3500000. The information scale is a local unbiased-covariance reference; moderate-noise RMSE also includes nonlinear bias. Every planned fit is retained.\n\n'
    supp=supp.replace('## References',s+'## References');p.write_text(supp,encoding='utf-8');print('Moderate-noise likelihood references added')


if __name__=='__main__':main()
