# Module for relativistic kinematics calculations.
# -------------------------------------------------

import numpy as np
import constants as const


def four_vector(E, p):
    E = np.atleast_1d(E)

    if p.ndim == 1:
        p = p.reshape(3,1)

    return np.vstack((E, p))


def get_four_momentum(vel, m):
    """
    Calculate the four-momentum (keV/c) for a particle with velocity 'vel' (m/s)
    and mass 'm' (keV/c**2).

    'vel' can also be a (3,N) array with velocities and 'm' can be a
    length-N array with masses.
    """

    vel = vel / const.c

    v2 = vel[0]**2 + vel[1]**2 + vel[2]**2

    gamma = 1/np.sqrt(1 - v2)

    E = gamma*m                           # total energy (keV)
    p = gamma*m*vel                       # three-momentum (keV/c)

    P = four_vector(E, p)

    return P


def mult_four_vectors(P1, P2):
    """
    Multiply two 4-vectors.
    Vectorized input works if 'P1' and/or 'P2' are (4,N) arrays.
    """

    prod = P1[0]*P2[0] - np.sum(P1[1:]*P2[1:], axis=0)

    return prod


def boost(P, beta):
    """
    Lorentz transformation of four-momentum P.

    All vector input should be 2D, i.e. shape (4,N) for four-vectors
    and (3,N) for three-vectors, where N is either equal to 1 or to
    the number of particles.
    """

    # Evaluate three-momentum parallel and perpendicular to beta
    E = P[0]
    p3 = P[1:]

    beta_mag = np.linalg.norm(beta, axis=0)
    i_nonzero = beta_mag != 0

    beta_dir = np.zeros_like(beta)
    beta_dir[:,i_nonzero] = beta[:,i_nonzero] / beta_mag[i_nonzero]  # avoid division by zero

    p3_par_mag = np.sum(p3*beta_dir, axis=0)
    p3_par = p3_par_mag * beta_dir
    p3_perp = p3 - p3_par

    # Gamma factor
    gamma = 1 / np.sqrt(1-beta_mag**2)

    # Lorentz boost
    E_prime = gamma*E - gamma*beta_mag*p3_par_mag
    p3_par_mag_prime = -gamma*beta_mag*E + gamma*p3_par_mag     # Only the parallel component of p3 is affected

    p3_prime = p3_par_mag_prime * beta_dir + p3_perp

    # Final result
    P_prime = four_vector(E_prime, p3_prime)

    return P_prime
