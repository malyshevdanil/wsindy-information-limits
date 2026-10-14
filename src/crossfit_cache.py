"""A bounded polynomial instrument map; coefficients fixed before simulations."""
import itertools
import numpy as np
from scipy.stats import qmc


def design_nodes(config,anchor):
    """Deterministic design, independent of physical Monte Carlo data."""
    unit=qmc.Sobol(4,scramble=True,seed=config['node_seed']).random_base2(config['node_power'])
    distribution=config.get('node_distribution','uniform')
    if distribution=='uniform':coordinates=2*unit-1
    elif distribution=='chebyshev':coordinates=np.cos(np.pi*unit)
    else:raise ValueError('Unknown cache design distribution')
    return np.asarray(anchor)+config['radius']*coordinates


def powers(degree):
    return np.array([p for p in itertools.product(range(degree+1),repeat=4) if sum(p)<=degree],dtype=int)


def features(q,anchor,radius,exponents):
    z=np.clip((np.asarray(q)-anchor)/radius,-1.,1.)
    shape=z.shape[:-1];z=z.reshape(-1,4)
    # Legendre basis is better conditioned than raw powers on the design cube.
    values=[np.polynomial.legendre.legvander(z[...,j],int(np.max(exponents[:,j]))) for j in range(4)]
    result=np.ones(z.shape[:-1]+(len(exponents),))
    for j in range(4):result*=values[j][...,exponents[:,j]]
    return result.reshape(shape+(len(exponents),))
