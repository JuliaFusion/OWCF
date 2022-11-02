# Tools for working with TRANSP ion distributions.
# ------------------------------------------------------


import os

from netCDF4 import Dataset
import numpy as np
from scipy.interpolate import griddata, RectBivariateSpline, interp1d

import sampler


class Thermal:

    def __init__(self, transp_output, ion='D'):

        self.transp_output = transp_output
        self.step = transp_output.step
        self.ion = ion.upper()
        self.pop = 'BULK'      # type of population

        # Time window
        self.t0 = transp_output.t0
        self.t1 = transp_output.t1

        # Volume elements
        dV = transp_output.get_variable('DVOL')
        dV[1] = dV[1]/1e6     # m^-3

        self.dV = dV
        self.V = dV[1].sum()

        # Temperature and density profiles
        T = transp_output.get_variable('TI')

        if self.ion == 'D':
            n = transp_output.get_variable('ND')
        if self.ion == 'T':
            n = transp_output.get_variable('NT')
        if self.ion == '3HE':
            n = transp_output.get_variable('NHE3')

        T[1] = T[1] / 1e3    # keV
        n[1] = n[1] * 1e6    # m^-3

        self.n_max = n[1].max()

        self.T = T
        self.n = n

        self.n_interp = interp1d(n[0], n[1], fill_value='extrapolate')

        # Rotation
        omega = transp_output.get_variable('OMEGA')
        if omega != -1:
            self.omega = omega
        else:
            self.omega = None

        # Numer of particles in each volume element
        self.N_vol = self.n[1]*self.dV[1]
        self.N_vol_interp = interp1d(n[0], self.N_vol, fill_value='extrapolate')
        self.N_vol_max = self.N_vol.max()

        # Total number of particles in the distribution
        self.N = np.sum(self.N_vol)


    def __repr__(self):
        return 'TRANSP thermal {} dist from {}'.format(self.ion, self.transp_output.out_file)


    def get_density(self, R, Z):

        X = self.transp_output.get_rho(R, Z)

        #n = np.interp(X, self.n[0], self.n[1], right=0)
        n = self.n_interp(X)
        n[X>1.0] = 0

        return n


    def get_N_vol(self, R, Z):

        X = self.transp_output.get_rho(R, Z)

        N_vol = self.N_vol_interp(X)
        N_vol[X>1.0] = 0.0

        return N_vol


    def get_temperature(self, R, Z):

        X = self.transp_output.get_rho(R, Z)

        T = np.interp(X, self.T[0], self.T[1], right=0)

        return T


    def get_rotation(self, R, Z):

        if self.omega is None:
            vrot =  np.zeros_like(R)
        else:
            X = self.transp_output.get_rho(R, Z)

            omega = np.interp(X, self.omega[0], self.omega[1], right=0)     # rad/s
            vrot = omega * R       # m/s

        return vrot


    def sample(self, n_samples=None, R=None, Z=None):
        """
        Sample from the distribution. Must input either 'n_samples',
        (in which case the positions are sampled), or both 'R' and 'Z',
        (in which case energy and pitch are sampled at the correspoinding positions).
        """

        if (R is None) and (Z is None):
            sample_pos = True

            # Draw random positions
            Rmin = self.transp_output.R.min()
            Rmax = self.transp_output.R.max()
            Zmin = self.transp_output.Z.min()
            Zmax = self.transp_output.Z.max()
            lims = ([Rmin, Zmin], [Rmax, Zmax])

            s = sampler.sample_acc_rej(self.get_N_vol, lims, self.N_vol_max, n_samples)

            R = s[:,0]
            Z = s[:,1]
        else:
            sample_pos = False
            n_samples = len(R)

        # The (kinetic) energy is distributed as a chi2 variable with 3 degrees of freedom
        T = self.get_temperature(R, Z)
        E = np.random.chisquare(3, size=n_samples) * 0.5 * T

        # Pitch distribution is isotropic
        p = sampler.sample_uniform([-1,1], n_samples)

        if sample_pos:
            return R, Z, E, p
        else:
            return E, p


    def __call__(self, R, Z, E, p):
        """ Evaluate the distribution at the given point(s). """

        T = self.get_temperature(R,Z)
        n = self.get_density(R,Z)

        f = n*np.sqrt(E/np.pi)*T**(-1.5)*np.exp(-E/T)

        return f



class FBM:

    def __init__(self, transp_output, ion='D', pop='NBI', jdotb=None):

        self.transp_output = transp_output

        step = transp_output.step

        if step is None:
            raise ValueError('Invalid time step!')

        self.step = step

        self.ion = ion.upper()
        self.pop = pop.upper()

        suffix = self.ion + '_' + self.pop

        # Set default sampling method for EP distributions
        self.EP_sampling_method = 'acc_rej'

        # Open the FBM netCDF file
        fbm_cdf = transp_output._openfbm()

        # Time window
        self.t0 = transp_output.t0
        self.t1 = transp_output.t1

        # Extract the distribution data
        self.X  = fbm_cdf.variables['X2D'][:]
        self.TH = fbm_cdf.variables['TH2D'][:]
        self.R  = fbm_cdf.variables['R2D'][:] / 100
        self.Z  = fbm_cdf.variables['Z2D'][:] / 100

        self.index = np.arange(len(self.X))

        self.dV = fbm_cdf.variables['BMVOL'][:] / 1e6    # m^-3
        self.V = np.sum(self.dV)

        # Read fast ion distribution
        self.F = fbm_cdf.variables['F_{}'.format(suffix)][:] / 2.0 * 1000.0 * 1e6

        self.E  = fbm_cdf.variables['E_{}'.format(suffix)][:] / 1000.0
        self.A  = fbm_cdf.variables['A_{}'.format(suffix)][:]
        self.EB = fbm_cdf.variables['EB_{}'.format(suffix)][:] / 1000.0
        self.AB = fbm_cdf.variables['AB_{}'.format(suffix)][:]
        self.dE = np.diff(self.EB)
        self.dA = np.diff(self.AB)

        if jdotb is None:
            self.jdotb = int(fbm_cdf.variables['nsjdotb'][:])
        else:
            self.jdotb = int(jdotb)

        if self.jdotb == -1:
            # Negate and flip pitch axis in order to always have
            # pitch values given relative to the B-field.
            self.F = self.F[:, ::-1, :]
            self.A = -self.A[::-1]

        # Make interpolants for the (E,p) distributions in each spatial zone
        self.F_Ep = []
        self.F_Ep_max = np.zeros_like(self.X)

        for i in range(len(self.X)):
            F = RectBivariateSpline(self.E, self.A, self.F[i].T, kx=1, ky=1)
            self.F_Ep.append(F)

            # Also store maximum distribution value at each point,
            # in order to know suitable sampling regions later.
            self.F_Ep_max[i] = self.F[i].max()

        self.Emin = self.E[0]
        self.Emax = self.E[-1]
        self.pmin = -1
        self.pmax = 1

        # Calculate fast ion density at each grid point
        dE = self.dE
        dA = self.dA
        self.n = np.sum( np.sum(self.F*dE, axis=-1) *dA, axis=-1 )   # m^-3
        self.n_max = self.n.max()

        self.Rmin = self.R.min()
        self.Rmax = self.R.max()
        self.Zmin = self.Z.min()
        self.Zmax = self.Z.max()

        # Number of particles in each volume element
        self.N_vol = self.n * self.dV
        self.N_vol_max = self.N_vol.max()

        # Total number of particles in the distribution
        self.N = np.sum(self.n * self.dV)

        # Integrate distribution over RZ plane
        # (this factor can be handy for normalizing MC calculations)
        self.C_RZ = np.sum(self.n * self.dV/(2*np.pi*self.R))

        # Prepare for inverse CDF sampling
        dAdE = dA[None,:,None] * dE[None,None,:]          # phase space volume elements
        Ptab = self.F * dAdE
        self.Ptab = Ptab.reshape(len(self.X), len(self.E)*len(self.A))   # 1D tables of probabilities
        self.Ctab = np.cumsum(self.Ptab, axis=-1)

        Egrid, Agrid = np.meshgrid(self.E, self.A)
        self.Etab = Egrid.flatten()
        self.Atab = Agrid.flatten()

        self.i_tab = np.arange(len(self.Etab))


    def __repr__(self):
        return 'TRANSP FBM {} {} dist from {}'.format(self.ion, self.pop, self.transp_output.fbm_files[self.step-1])


    def get_density(self, R, Z):
        X = self.transp_output.get_rho(R, Z)
        n = griddata((self.R, self.Z), self.n, (R,Z), method='nearest', fill_value=0)
        return np.where(X<1.0, n, 0.0)


    def get_N_vol(self, R, Z):
        X = self.transp_output.get_rho(R, Z)
        N_vol = griddata((self.R, self.Z), self.N_vol, (R,Z), method='nearest', fill_value=0)
        return np.where(X<1.0, N_vol, 0.0)


    def get_F_Ep(self, i_vol, E, p):
        F = np.zeros_like(E, dtype='d')

        valid_E = (E < self.Emax) & (E > self.Emin)
        valid_p = (p < self.pmax) & (p > self.pmin)
        valid = valid_E & valid_p

        F[valid] = self.F_Ep[i_vol](E[valid], p[valid], grid=False)

        return F


    def get_spatial_index(self, R, Z):

        # Map each R,Z value to the closest spatial grid point.
        i_spatial = griddata((self.R, self.Z), self.index, (R,Z), method='nearest')

        # Points outside the plasma should not get a valid index
        X = self.transp_output.get_rho(R, Z)
        i_spatial[X>1] = self.index[-1] + 1

        # Map each spatial grid point to the corresponding R,Z values
        in_vol = np.empty(len(self.R), dtype='object')
        for i in range(len(in_vol)):
            in_vol[i] = []

        for i, i_s in enumerate(i_spatial):
            if i_s > self.index[-1]:
                continue   # don't include points outside the plasma
            else:
                in_vol[i_s].append(i)

        return i_spatial, in_vol


    def sample(self, n_samples=None, R=None, Z=None):
        """
        Sample from the distribution. Must input either 'n_samples',
        (in which case the positions are sampled), or both 'R' and 'Z',
        (in which case energy and pitch are sampled at the corresponding positions).
        """

        if (R is None) and (Z is None):
            print('Sampling positions...')
            sample_pos = True

            # Draw random positions
            def fun(R,Z): return self.get_density(R,Z)*R
            lims = ([self.Rmin,self.Zmin], [self.Rmax,self.Zmax])
            max_val = self.n_max * 3.2    # a bit of a hack...

            s = sampler.sample_acc_rej(fun, lims, max_val, n_samples)

            R = s[:,0]
            Z = s[:,1]

        else:
            R,Z = np.atleast_1d(R,Z)
            sample_pos = False
            n_samples = len(R)

        # Check in which volumes the sampled/requested positions reside,
        # and sample E and p from the appropriate distributions.
        print('Sampling energy and pitch...')

        i_spatial, in_vol = self.get_spatial_index(R, Z)

        E = np.zeros(n_samples)
        p = np.zeros(n_samples)

        for i,ip in enumerate(in_vol):

            if len(ip) == 0:
                continue    # no points here, moving on...

            # Draw samples
            if self.F_Ep_max[i] == 0.0:
                # If the distribution is zero in this volume element
                # we return samples with zero energy.
                E[ip] = 0.0
                p[ip] = 0.0
            else:
                if self.EP_sampling_method == 'acc_rej':
                    def fun(E,p): return self.get_F_Ep(i,E,p)   # function to sample from
                    lims = ([self.Emin, self.pmin], [self.Emax,self.pmax])  # sampling region
                    
                    s = sampler.sample_acc_rej(fun, lims, self.F_Ep_max[i], len(ip))
                    
                    E[ip] = s[:,0]
                    p[ip] = s[:,1]

                elif self.EP_sampling_method == 'inv_cdf':

                    # Make 1D table of probabilities
                    r = np.random.rand(len(ip)) * self.Ctab[i][-1]
                    s = np.interp(r, self.Ctab[i], self.i_tab)
                    
                    E[ip] = np.interp(s, self.i_tab, self.Etab)
                    p[ip] = np.interp(s, self.i_tab, self.Atab)
                    
        if self.EP_sampling_method == 'inv_cdf':
            # Add jitter to pitch samples
            dp = np.interp(p, self.A, self.dA) / 2.0
            p = np.random.uniform(low=p-dp, high=p+dp) 

        print('Done!')

        if sample_pos:
            return R, Z, E, p
        else:
            return E, p


    def __call__(self, R, Z, E, p):
        """ Evaluate the distribution at the given point(s). """

        R, Z, E, p = np.atleast_1d(R, Z, E, p)

        F = np.zeros_like(R, dtype='d')

        # Only consider points inside the last closed flux surface
        X = self.transp_output.get_rho(R, Z)
        valid = X < 1

        R, Z, E, p = R[valid], Z[valid], E[valid], p[valid]
        F_valid = F[valid]

        # Find the closest spatial grid point
        i_spatial, in_vol = self.get_spatial_index(R, Z)

        # Loop over all volumes that contain requested points
        for i in np.unique(i_spatial):
            ip = in_vol[i]
            F_valid[ip] = self.get_F_Ep(i, E[ip], p[ip])

        F[valid] = F_valid

        return F
