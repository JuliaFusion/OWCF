# Wrapper functions for calculating neutron spectra using the
# 'sampler' and 'relscatt' modules.
# ---------------------------------------------------------------------------

import numpy as np
from scipy.interpolate import interp1d, RectBivariateSpline

import constants
import relkin
import relscatt
import sampler
import tokapart
import vcone


phi_hat = np.array([0,1.0,0]).reshape(3,1)    # toroidal basis vector


# Class for holding particle info
# --------------------------------


class Particle:

    def __init__(self, name):

        name = name.lower()
        self.name = name

        if name == 'n':
            self.long_name = 'neutron'
            self.u = constants.mn
            self.q = 0.0
        elif name == 'p':
            self.long_name = 'proton'
            self.u = constants.mp
            self.q = constants.e
        elif name == 'g':
            self.long_name = 'gamma'
            self.u = 0.0
            self.q = 0.0
        elif name == 'd':
            self.long_name = 'deuteron'
            self.u = constants.md
            self.q = constants.e
        elif name == 't':
            self.long_name = 'triton'
            self.u = constants.mt
            self.q = constants.e
        elif name == '3he':
            self.long_name = 'helium-3'
            self.u = constants.m3He
            self.q = 2*(constants.e) 
        elif name == '4he':
            self.long_name = 'helium-4'
            self.u = constants.m4He
            self.q = 2*(constants.e)
        elif name == '9be':
            self.long_name = 'beryllium-9'
            self.u = constants.m9Be   
            self.q = 4*(constants.e)
        elif name == '12c_gs': #
            self.long_name = 'carbon-12-gs'
            self.u = constants.m12C
            self.q = 6*(constants.e)
        elif name == '12c_1l': #
            self.long_name = 'carbon-12-1l'
            self.u = constants.m12C + 4439.82/constants.u_keV
            self.q = 6*(constants.e)
            self.E_g = 4438.94 # emitted gamma-ray energy (see Valentini25NF)
        elif name == '12c_2l': #
            self.long_name = 'carbon-12-2l'
            self.u = constants.m12C + 7654.07/constants.u_keV
            self.q = 6*(constants.e)
            self.E_g = 3213.79 # emitted gamma-ray energy (see Valentini25NF)
        elif name == 'proj': # In this case, the Particle class is simply used as a data placeholder to enable the workflow of forward.py 
            self.long_name = 'projection'
            self.u = 0.0
            self.q = 0.0
        else:
            raise ValueError('Invalid particle name')


    def __repr__(self):
        return 'Particle: {}'.format(self.long_name)


    @property
    def m(self):
        return self.u * constants.u_keV


# Class for holding reactant info and samples
# ----------------------------------------------

class Reactant:

    def __init__(self, name, n_samples=0, B_dir=[0,1,0]):

        self.particle = Particle(name)
        self._B_dir = B_dir

        # Initialize without MC sample
        self.P = None
        self.n_samples = int(n_samples)

    def __repr__(self):
        return 'Reactant: {}'.format(self.particle.long_name)

    @property
    def m(self):
        return self.particle.m
    
    @property
    def q(self):
        return self.particle.q

    @property
    def B_dir(self):
        return self._B_dir

    @B_dir.setter
    def B_dir(self, B):
        self._B_dir = np.array(B, 'd')

    # Handling the sample vectors
    @property
    def v(self):
        return self._v

    @v.setter
    def v(self, v_):
        if v_.shape != (3,self.n_samples):
            raise ValueError('Wrong shape of v (wrong number of samples?)')
        else:
            self._v = v_
            self.P = relkin.get_four_momentum(v_, self.m)

    def set_energy_pitch(self, E, p):
        self.E = E
        self.pitch = p

        self.v = tokapart.add_gyration(E, p, self.m, self.B_dir)


    # Convenience functions for sampling various distributions

    def sample_maxwellian_dist(self, T, pitch_range=[-1,1], v_rot=np.reshape(np.array([0.0,0.0,0.0]),(3,1))):
        """ Sample reactant four-momenta from a Maxwellian distribution with a given temperature T
        (with uniform pitch distribution in a given range pitch_range). Units in keV.

        v_rot: A plasma rotation velocity vector in (R, phi, z) can also be supplied (m/s). Usually,
        this plasma rotation vector points purely in the toroidal direction [0.0, 1.0, 0.0] or 
        [0.0,-1.0,0.0] but not always. Rotation along the magnetic field lines might also exist, which
        is not purely in the toroidal direction when a poloidal magnetic field is taken into account."""

        m = self.m

        # The (kinetic) energy is distributed as a chi2 variable with 3 degrees of freedom
        E = np.random.chisquare(3, size=self.n_samples) * 0.5 * T

        # Sample pitch values
        pitch = sampler.sample_uniform(pitch_range, self.n_samples)

        # Store results
        self.set_energy_pitch(E, pitch)

        # Add rotation
        if not (np.sum(v_rot) == 0.0):
            self.v = self.v + v_rot

    def sample_mono_dist(self, E0, pitch_range=[-1,1]):
        """ Sample four-momenta of reactant particles ('a' or 'b') with a mono-energetic
        distribution (with uniform pitch distribution in a given range).
        Units in keV. """

        # Sample energies
        E = sampler.sample_mono(E0, self.n_samples)

        # Sample pitch values
        pitch = sampler.sample_uniform(pitch_range, self.n_samples)

        # Store results
        self.set_energy_pitch(E, pitch)

    def sample_box_dist(self, Ehigh, Elow=0.0, pitch_range=[-1,1]):
        """ Sample four-momenta of reactant particles ('a' or 'b') with a box
        distribution (with uniform pitch distribution in a given range).
        Units in keV. """

        # Sample energies
        E = sampler.sample_uniform([Ehigh, Elow], self.n_samples)

        # Sample pitch values
        pitch = sampler.sample_uniform(pitch_range, self.n_samples)

        # Store results
        self.set_energy_pitch(E, pitch)

    def sample_E_dist(self, Ep, fp, pitch_range=[-1,1], quiet=True):
        """ Sample four-momenta of reactant particles ('a' or 'b') with a tabulated
        energy distribution fp given at the points Ep (with uniform pitch distribution
        in a given range). Units in keV. """

        # Function to sample from
        f = interp1d(Ep, fp, bounds_error=False, fill_value=0)

        # Sampling space
        lims = ([Ep[0]], [Ep[-1]])
        fmax = fp.max()

        # Sample energies
        s = sampler.sample_acc_rej(f, lims, fmax, self.n_samples, quiet=quiet)
        E = s[:,0]

        # Sample pitch values
        pitch = sampler.sample_uniform(pitch_range, self.n_samples)

        # Store results
        self.set_energy_pitch(E, pitch)

    def sample_EP_dist(self, Ep, pitch_p, fp, quiet=True):
        """ Sample four-momenta of reactant particles ('a' or 'b') with a tabulated
        (E,pitch) distribution fp given at the points (Ep, pitch_p). Units in keV. """

        # Function to sample from
        rbs = RectBivariateSpline(Ep, pitch_p, fp)

        def f(x,y): return rbs(x, y, grid=False)

        # Sampling space
        lims = ([Ep[0], pitch_p[0]], [Ep[-1], pitch_p[-1]])
        fmax = fp.max()

        # Sample energy and pitch
        s = sampler.sample_acc_rej(f, lims, fmax, self.n_samples, quiet=quiet)
        E = s[:,0]
        pitch = s[:,1]

        # Store results
        self.set_energy_pitch(E, pitch)


# Calculation of spectrum from sampled four-momenta.
# -----------------------------------------------------

class SpectrumCalculator:
    """
    Calculate product spectrum from reactions between particles with
    four-momenta 'Pa' and 'Pb'. 'u1' is the product emission direction
    (if None, each event gets a random direction).
    """

    def __init__(self, reaction='d-d', product_state='gs', B_dir=[0,1,0], n_samples=1e6): #

        self._B_dir = B_dir
        self._n_samples = int(n_samples)
        self.product_state = product_state.lower() #        
        self.reaction = reaction
        self.weights = None

        # 4*pi emission by default
        self.u1 = None

    def __repr__(self):
        return 'Spectrum calculator: {}'.format(self.reaction)

    @property
    def reaction(self):
        return self._reaction

    @reaction.setter
    def reaction(self, r):

        r = r.lower()

        a,b = r.split('-')

        self.reactant_a = a
        self.reactant_b = b

        if r == 'd-d':
            # D(d,n)3He reaction
            self.product_1 = 'n'
            self.product_2 = '3he'

        elif all(p in [a,b] for p in ['p','t']):
            # T(p,g)4He
            self.product_1 = 'g'
            self.product_2 = '4he'

        elif all(p in [a,b] for p in ['d','t']):
            # T(d,n)4He reaction
            self.product_1 = 'n'
            self.product_2 = '4he'

        elif all(p in [a,b] for p in ['d','3he']):
            # 3He(d,p)4He reaction
            self.product_1 = 'p' # Since the proton is a charged particle, as of now, only an energy spectrum for 4*pi emission can be returned
            self.product_2 = '4he'

        elif all(p in [a,b] for p in ['4he','9be']): #
            # 9Be(4He,12C)n reaction
            self.product_1 = '12c_'+self.product_state
            self.product_2 = 'n'

        elif b == 'proj': # In this case, the SpectrumCalculator class is just being used to fascilitate the computation of projected velocities in forward.py
            self.product_1 = 'n' # Then the products don't matter
            self.product_2 = 'n' # Then the products don't matter
        else:
            raise ValueError('Invalid reaction')

        self._reaction = r

    @property
    def reactant_a(self):
        return self._reactant_a

    @reactant_a.setter
    def reactant_a(self, name):
        self._reactant_a = Reactant(name, n_samples=self.n_samples, B_dir=self.B_dir)

    @property
    def reactant_b(self):
        return self._reactant_b

    @reactant_b.setter
    def reactant_b(self, name):
        self._reactant_b = Reactant(name, n_samples=self.n_samples, B_dir=self.B_dir)

    @property
    def product_1(self):
        return self._product_1

    @product_1.setter
    def product_1(self, p):
        self._product_1 = Particle(p)

    @property
    def product_2(self):
        return self._product_2

    @product_2.setter
    def product_2(self, p):
        self._product_2 = Particle(p)

    @property
    def ma(self):
        return self.reactant_a.m

    @property
    def mb(self):
        return self.reactant_b.m
    
    @property
    def qa(self):
        return self.reactant_a.q
    
    @property
    def qb(self):
        return self.reactant_b.q

    @property
    def m1(self):
        return self._product_1.m

    @property
    def m2(self):
        return self._product_2.m

    @property
    def B_dir(self):
        return self._B_dir

    @B_dir.setter
    def B_dir(self, B):

        self._B_dir = np.array(B, 'd')

        self.reactant_a.B_dir = B
        self.reactant_b.B_dir = B

    @property
    def n_samples(self):
        return self._n_samples

    @n_samples.setter
    def n_samples(self, n):
        self._n_samples = n

        self.reactant_a.n_samples = n
        self.reactant_b.n_samples = n

    @property
    def weights(self):
        if self._weights is None:
            return np.ones(self.n_samples)/self.n_samples
        else:
            return self._weights

    @weights.setter
    def weights(self, w):
        self._weights = w

    def __call__(self, bins=None, bin_width=25.0, normalize=False):
        """
        Compute reactant spectrum. The units of the spectrum depends on the
        units of the weights attribute. The 'normalize' keyword
        can be used to normalize the spectrum to unit sum.

        If 'bins' is given it is passed as a keyword argument to np.histogram.
        If bins=None (default), the bins will be constructed so as to cover the full
        range of calculated product energies.

        The bin width of the spectrum histogram can be specified with 'bin_width' (keV).
        (only effective if bins=None).

        EXPERIMENTAL: Additionally, both the 'bins' and the 'bin_width' attributes can be length-2 lists,
        in which case the spectrum is resolved both in energy and cos(theta), where
        theta is the angle between the B-field and the emission direction.
        """

        # Reactant momenta
        Pa = self.reactant_a.P
        Pb = self.reactant_b.P

        # Emission directions
        if self.u1 is None:
            u1 = sampler.sample_sphere(self.n_samples)   # random emission directions
        elif self.product_state == 'gs': #
            u1 = np.asarray(self.u1) / np.linalg.norm(self.u1, ord=2, axis=0) # norm=1 is forced, since it has to be a unit vector
            if u1.ndim == 1:
                u1 = np.array(u1).reshape(3,1)   # same emission direction for all particles
        else:   # If it's the case of a 2-step gamma reaction, a dedicated function is called to sample the excited nuclei direction. 
                # In this case, N.B. that 'u1' is swapped with the new direction array 'u_decay_g' so that the code works in both 1- and
                # 2-step cases with the same workflow. 'u1' becomes the random direction of the nuclei, while 'u_decay_g' the direction of
                # the gamma to the detector.
            u1, fourpi_fraction = self._sample_u_exNucl(Pa, Pb)    # The fraction of domain for the extracted scattering angle is needed as
                                                                   # an additional weight factor, which is multiplied to the rest later.
            if self.u1.shape[1] == 1:
                u_decay_g = self.u1 / np.linalg.norm(self.u1, ord=2, axis=0)
                if u_decay_g.ndim == 1:
                    u_decay_g = np.array(u_decay_g).reshape(3,1)   # same emission direction for all particles
            else:
                u_decay_g = self.u1 / np.linalg.norm(self.u1, ord=2, axis=0)

            self.weights *= fourpi_fraction # Additional sampling happens within a restricted domain, which is different depending on the
                                            # energy of the CM. Therefore, an array with values in [0,1] is factored in with the rest of the 
                                            # weights.

        # Compute product energies
        P1= relscatt.two_body_event(Pa, Pb, self.m1, self.m2, u1)
        
        # Compute reactivity
        sigmav = relscatt.get_reactivity(self.ma, self.mb, Pa, Pb, P1, 
                                         reaction=self.reaction, product_state=self.product_state)   #

        # Bin all events
        weights_tot = self.weights * sigmav        

        # If a relativistic Doppler shift has to be added on top
        if self.product_state != 'gs':
            beta_1 = P1[1:] / P1[0]
            gamma_1 = (1 - np.sum(beta_1**2,axis=0))**-0.5
            E_s = self.product_1.E_g  # energy at source
            E_r = E_s / gamma_1 / (1 - np.sum(beta_1*u_decay_g,axis=0)) # energy at receiver

            P1 = relkin.four_vector(E_r, E_r*u_decay_g)    # P1 variable now refers to a different quantity completely! We leave the excited nucleus alone
                                                            # and only care about the gamma from here on.

        result = self._make_spec_hist(P1, bins, bin_width, weights_tot, normalize)

        return result

    def _sample_u_exNucl(self, Pa, Pb):   # Here the magic happens. The sampling in the acceptance cone for the excited nuclei direction is made.
        p_tot = Pa[1:] + Pb[1:]
        ortho_e = np.ones(p_tot.shape) * 2**-0.5    # A direction which is as little aligned as possible to p_tot has to be picked. This is not 
                                                    # much elegant, but works every time.
        np.put_along_axis(ortho_e,np.argmax(p_tot,axis=0,keepdims=True),0,axis=0)   # See above comment.
        e1, e2, e3 = tokapart.get_basis_vectors(p_tot, ortho_e)    # The tokapart function comes in handy. A minor change has been made to 
                                                                   # make it more versatile and use whatever orthogonal direction to define the 
                                                                   # 3D basis.

        m1_sq,m2_sq = self.m1**2, self.m2**2
        E_tot_sq = (Pa[0] + Pb[0])**2
        p_tot_mag_sq = np.linalg.norm(p_tot, axis=0)**2
        Mtot_sq = E_tot_sq - p_tot_mag_sq
        Cm_sq = (Mtot_sq + m1_sq - m2_sq)**2

        cos_csi_sq_min = (4*m1_sq*E_tot_sq-Cm_sq) / (4*m1_sq*p_tot_mag_sq)    # This comes from kinematics, so energy and momentum conservation.
        if np.count_nonzero(cos_csi_sq_min>1) > 0: # just a sanity check, it should not be possible
            raise Exception("Minimum cosine of csi is ill-defined!")
        cos_csi_min = np.sqrt(cos_csi_sq_min,dtype=complex)    # N.B. Cast as complex.
        cos_csi_min[np.abs(cos_csi_min.imag) > 0] = -1     # If the minimum of the squared cosine is negative, it means that every angle is acceptable! 
        
        cos_csi = np.random.uniform(cos_csi_min.real,1,cos_csi_min.size)   # The cosine of csi is finally sampled in the proper region of its domain. 
        phi = np.random.rand(cos_csi.size)*2*np.pi    # The scattering is assumed isotropic in the azimuthal angle, and accounts for downward direction
        f_o_b = np.random.choice([+1,-1], cos_csi.size)    # Forward or backward? This boolean-like quantity assign the same proability to either of the directions

        sin_csi = (1-cos_csi**2)**0.5
        u_exc = sin_csi*np.sin(phi)*e1 + f_o_b*cos_csi*e2 + sin_csi*np.cos(phi)*e3 # The direction versor is built, and later normalized to 1 just to be sure

        return u_exc / np.linalg.norm(u_exc, ord=2, axis=0), 2*np.pi*(1 - cos_csi_min.real) / 4/np.pi

    def _make_spec_hist(self, P1, bins, bin_width, weights, normalize):

        E1 = P1[0] - np.sqrt(relkin.mult_four_vectors(P1,P1),dtype=complex).real #! kinetic energy (keV).  This expression is valid for the 2-step photons as well

        # How to bin events
        if bins is None:

            return_bins = True

            if type(bin_width) == list:
                # Set up for computing 2D spectrum with default bins
                bins_E = np.arange(E1.min()-5*bin_width[0], E1.max()+6*bin_width[0], bin_width[0])
                bins_A = np.arange(-1, 1+bin_width[1], bin_width[1])

                bins = [bins_E, bins_A]

            else:
                # Set up for computing energy spectrum with default bins
                bins = np.arange(E1.min()-5*bin_width, E1.max()+6*bin_width, bin_width)

        else:
            return_bins = False

        # Bin in 1D or 2D
        if type(bins) == list:
            # Bin in both energy and emission direction
            A1 = np.dot(P1[1:].T, self.B_dir)
            p1 = np.linalg.norm(P1[1:], ord=2, axis=0)
            b = np.linalg.norm(self.B_dir, ord=2, axis=0)
            A1 = A1 / (p1*b)      # cosine of emission angles w.r.t. B-field

            spec = np.histogram2d(E1, A1, bins=bins, weights=weights)[0]
        else:
            spec = np.histogram(E1, bins=bins, weights=weights)[0]

        # Possibly normalize to unit sum
        if normalize:
            spec = spec / spec.sum()

        # Return results
        if return_bins:

            if spec.ndim == 2:
                # Return spectrum and bin edges of the 2D histogram
                return spec, bins[0], bins[1]

            else:
                # Return spectrum and bin centers of the 1D histogram
                E_axis = 0.5*(bins[1:] + bins[:-1])
                return spec, E_axis

        else:
            # Return only the spectrum, since the bins were supplied as input
            return spec
