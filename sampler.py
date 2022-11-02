# Routines for sampling different types of distributions
# ---------------------------------------------------------

import numpy as np
import scipy.integrate as integrate


def sample_mono(val, n_samples):
    """ Sample delta function. """

    s = val * np.ones(n_samples)

    return s

def sample_uniform(val_range, n_samples):
    """ Sample values from uniform distribution. """

    s = np.random.uniform(val_range[0], val_range[1], size=n_samples)

    return s

def sample_discrete(vals, n_samples, weights=None):
    """ Sample values from a 1D discrete distribution. """

    vals = np.atleast_1d(vals)

    if weights is not None:
        weights = np.atleast_1d(weights)
        weights = weights / weights.sum()

    s = np.random.choice(vals, p=weights, size=n_samples)

    return s

def sample_inv_trans(f, x, n_samples):
    """ Randomly sample from a tabulated distribution representing f(x), by
    means of inverse transform sampling. """

    if any(f<0):
        raise ValueError('Distribution must be non-negative!')

    # Inverse transform sampling
    cdf = np.zeros_like(f)
    cdf[1:] = integrate.cumtrapz(f, x=x)

    r = np.random.rand(n_samples) * cdf[-1]
    s = np.interp(r, cdf, x)

    return s

def sample_acc_rej(f, lims, fmax, n_samples, quiet=True):
    """ Randomly sample from a distribution f, using the acceptance-rejection method.
    'f' should be callable, taking the dependent variables input. The dimensionality of f
    is inferred from the tuple 'lims', which should contain the ranges to sample from, i.e.

         lims = ([x0_min,x1_min, ...], [x0_max,x1_max, ...])

    'fmax' is the maximum number that f can take. """

    # Sample points xp and fp in the hypercube bounded by 'lims' and 'fmax', and keep only
    # the points that fall below f(*xp). Loop until we have collected enought samples.
    ndims = len(lims[0])

    s = np.zeros((n_samples, ndims))
    n_to_sample = n_samples
    n_sampled = 0

    while n_to_sample > 0:

        # Sample points
        xp = np.random.uniform(low=lims[0], high=lims[1], size=(n_to_sample,ndims))
        fp = np.random.uniform(low=0, high=fmax, size=n_to_sample)

        # Accept/reject
        x_acc = xp[fp < f(*xp.T)]

        # Save accepted points
        i0 = n_sampled

        n_sampled += len(x_acc)
        i1 = n_sampled

        s[i0:i1] = x_acc

        # How many points to sample next time
        n_to_sample -= len(x_acc)

        if not quiet:
            print('{:.2f}%'.format(n_sampled/n_samples*100))

    return s

def sample_sphere(n_samples):
    """ Pick Cartesian coordinates of points that are uniformly distributed over the surface
    of a unit sphere. Output has the shape (3, n_samples). """

    # Uniformly sample cosine of polar angle
    u = np.random.uniform(-1, 1, size=n_samples)

    # Uniformly sample the azimuthal angle
    theta = np.random.uniform(0, 2*np.pi, size=n_samples)

    # Cartesian coordinates
    x = np.sqrt(1-u**2) * np.cos(theta)
    y = np.sqrt(1-u**2) * np.sin(theta)
    z = u

    # Return array
    v = np.vstack((x,y,z))

    return v

def sample_mcmc(F, x0, dx, n_samples):
    """
    Sample F(x) by means of Markov Chain Monte Carlo sampling. 
    'F' should be callable, accepting a vector 'x' as input. 
    'x0' is the startin point for the Markov chain.
    'dx' is a suitable step size.
    """

    # Parameters for controlling how to step to new points
    mean_step = np.zeros_like(x0)
    cov_step = dx * np.eye(len(x0))
    
    # Array to store the resulting Markov chain
    chain = np.zeros((n_samples,len(x0)))

    for i in range(n_samples):
        
        # Evaluate function at the present point
        f0 = F(x0)
        
        # Step to next point
        x1 = x0 + np.random.multivariate_normal(mean_step, cov_step)

        # Evaluate function at the new point
        f1 = F(x1)

        # Calculate the Markov ratio
        r = f1/f0

        # Accept new point unconditionally if r >= 1,
        # but if r < 1, accept it with probability r
        if r < 1:
            u = np.random.rand()

            if u > r:
                x1 = x0

        # Store new point and roll values for next step
        chain[i] = x1
        x0 = x1

    return chain

        
