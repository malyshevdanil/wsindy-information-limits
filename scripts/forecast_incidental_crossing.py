"""Algebraic pre-main crossing forecast from validated curvature and past data.

Uses no main MC outcomes. Historical dispersion is a labelled sensitivity
scenario, not the variance of the new cross-fitted estimator.
"""
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1]


def crossing(weak_variance,mle_variance,bias_squared):
    assert bias_squared>0
    return max(0.,float((weak_variance-mle_variance)/bias_squared))


def main():
    validation=json.loads((ROOT/'reports/incidental_curvature_validation.json').read_text())
    assert validation['completed'] and validation['quadrature_convergence_verified']
    c=validation['results']['results'][-1];a=np.array(c['bias_coefficient']);I=np.array(c['mean_information'])
    bias_constant=float(a@a/4);mle_constant=float(np.trace(np.linalg.inv(I))/4)
    levels=np.linspace(.005,.2,100);curves=[]
    for inflation in (.01,.05,.1,.3):
        # Fractional increase in sqrt(n) RMSE over the efficient MLE SD.
        multiplier=(1+inflation)**2
        curves.append(dict(label=f'Weak SD +{100*inflation:g}%',kind='assumed relative SD inflation',
            values=[dict(level=float(s),n_star=crossing(multiplier*s*s*mle_constant,s*s*mle_constant,s**4*bias_constant)) for s in levels],
            at_planned_levels=[dict(level=s,n_star=crossing(multiplier*s*s*mle_constant,s*s*mle_constant,s**4*bias_constant)) for s in [.05,.1,.2]]))
    previous=json.loads((ROOT/'reports/heterogeneous_records_results.json').read_text())
    for method in ('oracle120','fission480'):
        g=next(g for g in previous['groups'] if g['n']==800 and g['level']==.05 and g['method']==method)
        weak_constant=800*g['relative_dispersion']**2/.05**2
        curves.append(dict(label='Historical '+method,kind='past dispersion; assumes linear noise scaling; new cross-fitting variance pending',
            values=[dict(level=float(s),n_star=crossing(s*s*weak_constant,s*s*mle_constant,s**4*bias_constant)) for s in levels],
            at_planned_levels=[dict(level=s,n_star=crossing(s*s*weak_constant,s*s*mle_constant,s**4*bias_constant)) for s in [.05,.1,.2]]))
    report=dict(completed=True,scope='Leading-order algebra plus explicitly labelled variance sensitivity scenarios',
        main_outcomes_used=False,formula='n_star=(A_weak(sigma)-A_MLE(sigma))_+ / ||b(sigma)||^2; relative coefficient metric divides by ||theta||^2=4',
        bias_coefficient=a.tolist(),relative_bias_squared_constant=bias_constant,mle_variance_constant=mle_constant,
        curves=curves,assumptions=['Negligible finite-n bias of weak coefficient estimate','MLE bias sigma^2*a and variance sigma^2*I^-1/n',
        'Fixed-noise higher-order corrections to MLE dispersion are pending','Historical scenarios are not forecasts of the new cross-fitted estimator'],
        extrapolation='The predicted crossing beyond n=800 is an extrapolation, not an observed win')
    (ROOT/'reports/incidental_crossing_forecast.json').write_text(json.dumps(report,indent=2))
    fig,ax=plt.subplots(figsize=(7,4.5))
    for curve in curves:ax.semilogy([100*r['level'] for r in curve['values']],[r['n_star'] for r in curve['values']],label=curve['label'],linestyle='--' if curve['kind'].startswith('past') else '-')
    ax.axhline(800,color='black',linestyle=':',label='Largest planned n = 800')
    ax.set(xlabel='Noise amplitude (%)',ylabel='Predicted crossing n*',title='Pre-main sensitivity forecast: leading-order bias')
    ax.grid(alpha=.25);ax.legend(fontsize=8);fig.tight_layout()
    folder=ROOT/'reports/figures';folder.mkdir(exist_ok=True)
    for suffix in ('png','pdf'):fig.savefig(folder/f'incidental_crossing_forecast.{suffix}',dpi=200)
    plt.close(fig)
    print(json.dumps({c['label']:c['at_planned_levels'] for c in curves},indent=2))


if __name__=='__main__':main()
