# Module for relativistic scattering calculations.
# --------------------------------------------------

import numpy as np
import relkin
from constants import c       # speed of light
import fusreact

sigma_diff = fusreact.sigma_diff   # could be changed by the user if required


def two_body_event(Pa, Pb, m1, m2, u1):
    """
    Calculate a possible 4-vector (along direction 'u') of particle 1 in
    the two-body event

           a + b --> 1 + 2

    'Pa' and 'Pb' are the incoming 4-momenta (keV/c).
    'm1' and 'm2' are the masses of the products (keV/c**2).

    All vector input should be 2D, i.e. shape (4,N) for four-vectors
    and (3,N) for three-vectors, where N is either equal to 1 or to
    the number of events.
    """

    # Normalize emission directions
    u1 = u1 / np.linalg.norm(u1, axis=0)

    # Total initial 4-momentum
    Ptot = Pa + Pb
    Mtot = np.sqrt(relkin.mult_four_vectors(Ptot, Ptot))      # total invariant mass

    # Calculate magnitude of 3-momentum
    p1_mag = invariant_to_momentum(m2**2, Mtot, Ptot, m1, u1)

    # Energy
    E1 = np.sqrt(p1_mag**2 + m1**2)

    # 4-vector
    P1 = relkin.four_vector(E1, p1_mag*u1)

    return P1


def three_body_event(Pa, Pb, m1, m2, m3, u1):
    """
    Calculate a possible 4-vector (direction 'u') of particle 1 in
    the three-body event

         a + b --> 1 + 2 + 3

    'pa' and 'pb' are the incoming 4-momenta (keV/c). 'm1', 'm2' and 'm3' are
    the masses of the products (keV/c**2).

    All vector input should be 2D, i.e. shape (4,N) for four-vectors
    and (3,N) for three-vectors, where N is either equal to 1 or to
    the number of events.
    """


def invariant_to_momentum(inv, Mtot, Ptot, m1, u1):
    """ Solution to the four-momentum conservation equation. """

    Cm  = Mtot**2 + m1**2 - inv
    p0u = np.sum(Ptot[1:] * u1, axis=0)    # projection of CM momentum on emission direction
    E0  = Ptot[0]

    first_term = Cm * p0u

    sq_term1 = Cm**2 * E0**2
    sq_term2 = 4*m1**2 * E0**2 * (E0**2 - p0u**2)
    second_term = np.sqrt(sq_term1 - sq_term2)

    p1 = (first_term + second_term) / (2 * (E0**2 - p0u**2))   # magnitude of three-momentum

    return p1


def dalitz_sample(Mtot, m1, m2, m3):
    """
    Sample final states from the dalitz plot resulting from the input
    masses (keV/c**2). The sampling is done uniformly with the shooting method.
    """


def get_reactivity(ma, mb, Pa, Pb, P1, reaction='d-d', product_state='GS', extra_out=False):
    """
    Calculate the reactivity (m**3/sr/s) for a reaction with
    given four-momenta (keV/c) for the reactants ('Pa' and 'Pb') and the
    product of interest ('P1').

    All vector input should be 2D, i.e. shape (4,N) for four-vectors
    and (3,N) for three-vectors, where N is either equal to 1 or to
    the number of events.
    """

    # CM quanteties
    Ptot = Pa + Pb
    Mtot = np.sqrt(relkin.mult_four_vectors(Ptot, Ptot))
    Tcm = Mtot - (ma + mb)

    # Relative velocity
    beta_a = Pa[1:] / Pa[0]
    Pb_a = relkin.boost(Pb, beta_a)   # boost particle b into the rest frame of particle b
    vrel = Pb_a[1:] / Pb_a[0]
    vrel_mag = np.linalg.norm(vrel, axis=0) * c      # m/s

    # Jacobian for transforming between CMS and LAB
    gamma_cm = Ptot[0] / Mtot
    beta_cm =  Ptot[1:] / Ptot[0]
    P1_cm = relkin.boost(P1, beta_cm)
    p1cm_mag = np.linalg.norm(P1_cm[1:], axis=0)
    p1_mag = np.linalg.norm(P1[1:], axis=0)

    u_dot_beta = np.sum(P1[1:]*beta_cm, axis=0) / p1_mag
    jacobian = p1_mag**2 / (gamma_cm * p1cm_mag * (p1_mag - P1[0]*u_dot_beta))

    # Angular differential cross section in the CMS
    p_rel = Pa[1:] - Pb[1:]
    p_rel_mag = np.linalg.norm(p_rel, axis=0)

    nonzero = p_rel_mag > 0
    u_rel = np.zeros_like(p_rel)
    u_rel[:,nonzero] = p_rel[:,nonzero] / p_rel_mag[nonzero]   # unit vector along relative motion of reactants

    costheta = np.sum(P1_cm[1:]*u_rel, axis=0) / p1cm_mag      # cosine of emission angle of particle 1 in the CMS

    sigma = sigma_diff(Tcm, costheta, reaction=reaction, product_state=product_state)

    # Reactivity (m**3/sr/s)
    sigmav = sigma * vrel_mag * np.abs(jacobian)

    if extra_out:
        return sigmav, jacobian, costheta
    else:
        return sigmav
