# Collection of useful functions for working with particles
# moving around in a tokamak.
# --------------------------------------------------------------

import numpy as np
import constants as const


def get_speed(E, m):
    """
    Calculate the speed (in m/s) of a particle with kinetic energy E (keV)
    and mass 'm' (keV/c**2). Relativistic kinematics is used.
    """

    E = E + m                              # total energy (keV)
    v = np.sqrt(1 - (m/E)**2)*const.c      # speed (m/s)

    return v


def get_energy(v, m):
    """
    Calculate the kinetic energy 'E' (in keV) of a particle with velocity 'v' (m/s)
    and mass 'm' (keV/c**2).

    'v' can be a length-N array with the values of the speed, or a (3,N) array
    with the velocity components.

    Relativistic kinematics is used.
    """

    v = np.atleast_1d(v)

    c = const.c         # speed of light (m/s)

    v = np.linalg.norm(v, ord=2, axis=0)
    gamma = 1/np.sqrt(1-(v/c)**2)

    E = (gamma-1)*m     # keV

    return E


def get_vpar_vperp(E, p, m):
    """
    Calculate velocity components (in m/s) parallel and perpendicular
    to the plasma magnetic field, for a particle of mass 'm' (keV/c**2)
    with given values of energy 'E' (keV) and pitch p = v_par/v.
    """

    v = get_speed(E, m)
    #p = max(-1.0,min(1.0,p)) # Clamping, to avoid weird stuff
    v_par  = v*p
    v_perp = v*np.sqrt(1-p**2)

    return v_par, v_perp


def get_EP(v_par, v_perp, m):
    """
    Calculate energy E (keV) and pitch p = v_par/v
    for given values of the velocity components (in m/s)
    parallel and perpendicular to the plasma magnetic field,
    for a particle with mass 'm' (keV/c**2).
    """

    v = np.sqrt(v_par**2 + v_perp**2)

    E = get_energy(v, m)
    p = v_par/v

    return E, p


def get_basis_vectors(B,ortho_e = np.array([[1],[0],[0]])):
    """
    Construct a set of orthogonal basis vectors that aligns with
    a given magnetic field vector 'B' (or a (3,N) array of N vectors).

    The components of the basis vectors are evaluated in cylindrical coordinates.

    (The unit of 'B' is not important, only its direction is used by the function.)
    """

    if np.ndim(B) == 1:
        B = B.reshape(3,1)    # reshape needed for broadcasting to work

    # Normalize the B-field vector
    B_norm = np.linalg.norm(B, ord=2, axis=0)
    b = B/B_norm

    #!# cosine of csi, the angle in between b and preferred orthogonal direction ortho_e.
    #!# Both are assumed to be expressed in cylindrical coordinates, 
    #!# shaped as (3,N) arrays, and normalized to 1.
    cos_csi = np.sum(b*ortho_e, axis=0)

    # Construct basis vectors that align with B
    # (only works if the phi component of B is non-zero)
    e1 = ortho_e - cos_csi*b
    e1 = e1 / np.linalg.norm(e1, ord=2, axis=0)
    e2 = b
    e3 = np.cross(e1,e2, axis=0)

    return e1, e2, e3


def sample_gyro_angle(N):
    """
    Sample 'N' random gyro angles in the interval [0,2*pi).
    """

    return 2*np.pi*np.random.rand(N)


def get_velocity(v_par, v_perp, e1, e2, e3, theta_g):
    """
    Calculate the velocity components (m/s) for given values of
    the speeds parallel and perpendicular to the plasma magnetic field.

    'e1', 'e2' and 'e3' are orthogonal basis vectors, such that 'e2' is
    parallel to the magnetic field vector.

    'theta_g' is the gyro radius.
    """

    v_par = v_par * e2
    v_perp = v_perp*(-np.sin(theta_g)*e1 + np.cos(theta_g)*e3)

    return v_par + v_perp


def get_Larmor_radius(v_perp, B, Z, m):
    """
    Calculate the Larmor radius for a particle with speed v_perp (m/s)
    relative to the plasma magnetic field 'B' (T).

    'v_perp' and 'B' can be length-N arrays with the relevant magnitudes.
    Alternatively, B could also be a (3,N) array with vector components.

    'Z' and 'm' are the atomic charge number and mass (keV/c**2) of the particle(s).
    """

    B = np.atleast_2d(B)

    q = Z*const.e

    B = np.linalg.norm(B, ord=2, axis=0)

    m = m * 1000 * const.e / const.c**2     # keV/c**2 -> kg
    r_L = m*v_perp / (q*B)

    return abs(r_L)


def add_gyration(E, p, m, B, theta_g=None, flr=False, Z=None, Rg=None, zg=None):
    """
    Calculate the cylindrical velocity vector of a particle in a tokamak.
    'theta_g' is the gyro angle; if 'None' it is sampled randomly from [0,2*pi).

    'E' is the particle kinetic energy (keV). 'p' = v_par/v is the pitch.
    'm' is the mass (keV/c**2).
    'B' is the magnetic field vector at the gyrocenter of the particle (T).

    If 'flr' = True, the function also returns the position of the particle
    (cylindrical coordinates). In this case, the gyrocenter position ('Rg','zg')
    needs to be provided (in meters), as well as the atomic charge number 'Z'.

    Vectorized input works if 'B' is given as a (3,N) array and all other
    input are length-N arrays (or scalars).
    """

    if np.ndim(B) == 1:
        B = B.reshape(3,1)    # reshape needed for broadcasting to work

    E = np.atleast_1d(E)

    # Local coordinate system that aligns with B
    e1, e2, e3 = get_basis_vectors(B)

    # Compute the velocity vector in the local coordinate system.
    v_par, v_perp = get_vpar_vperp(E, p, m)    # speeds

    if theta_g == None:
        theta_g = sample_gyro_angle(len(E))

    v = get_velocity(v_par, v_perp, e1, e2, e3, theta_g)    # velocity vectors

    # Compute particle position, if required
    if flr:
        # Larmor radius
        r_L = get_Larmor_radius(v_perp, B, Z, m)

        # The gyro-center is at (Rg,phi_g,zg). Actual position should be
        # displaced by one Larmor radius in the (e1,e3) plane.
        R_hat = np.array([1,0,0]).reshape(3,1)    # reshape needed for broadcasting to work
        z_hat = np.array([0,0,1]).reshape(3,1)

        x = Rg*R_hat + zg*z_hat + r_L*(np.cos(theta_g)*e1 + np.sin(theta_g)*e3)

        return v, x

    else:
        return v
