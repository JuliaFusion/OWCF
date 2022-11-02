# A module for loading and retrieving results from TRANSP simulations.
# ----------------------------------------------------------------------

import os

from netCDF4 import Dataset
import numpy as np
from scipy.interpolate import griddata

import sampler



class TranspOutput:
    """
    A class for storing the filenames of all relvant output files from 
    a TRANSP simulation, and methods for reading the output variables.

    If an AC file output step is given, it is also possible to load the 
    flux surfaces for the corresponding time window.
    """
    
    def __init__(self, runid, step=None, **kwargs):

        runid = runid.upper()
        self.runid = runid

        if 'out_file' in kwargs.keys():
            # Output netCDF file provided explicitly
            self.out_file = kwargs['out_file']
            self.fbm_files = kwargs.get('fbm_files', [])
            self.fbm_files.sort()

        else:
            # No output file provided, try to find it in the database
            shot = runid[:-3]
            seq = runid[-3:]
            self.out_file = ( '/common/transp_shared/Data/result/JET/{}/{}/{}.CDF'
                             .format(shot, seq, seq) )
            
            top_dir = kwargs.get('top_dir', None)
            if top_dir is None:
                transp_dir = os.environ['TRANSP']
                shot_dir = os.path.join(transp_dir, shot)
            else:
                shot_dir = top_dir
            
            run_dir = os.path.join(shot_dir, runid)
            
            # Find associated FBM files
            self.fbm_files = []
            
            for D in os.listdir(run_dir):
                if D[:4] == 'DATA':
                    fbm_path = os.path.join(run_dir,D,'{}_fi_{}.cdf'.format(runid,D[-1]))
                    if os.path.isfile(fbm_path):
                        self.fbm_files.append(fbm_path)

            self.fbm_files.sort()
            
        if len(self.fbm_files) > 0:
            # Get AC output times
            
            fbm_times = np.zeros((len(self.fbm_files),2))
            
            for i,f in enumerate(self.fbm_files):
                fbm_cdf = Dataset(f, mode='r')
            
                t_mid = fbm_cdf.variables['TIME'][0]
                dt    = fbm_cdf.variables['DT_AVG'][0]
                
                fbm_times[i,0] = np.round(t_mid - dt, 3)
                fbm_times[i,1] = np.round(t_mid + dt, 3)

                fbm_cdf.close()

            self.fbm_times = fbm_times

        self.step = step


    def __repr__(self):
        return 'TRANSP output from: {}'.format(self.out_file)

        
    @property
    def step(self):
        return self._step

    @step.setter
    def step(self, s):
        if s is None:
            self.t0 = None
            self.t1 = None
            self.tmid = None
        else:
            i = s - 1
            n_steps = len(self.fbm_files)
            if (s < 1) or (s > n_steps):
                raise ValueError('Only {} steps available!'.format(n_steps))
            else:
                self.t0 = self.fbm_times[i,0]
                self.t1 = self.fbm_times[i,1]
                self.tmid = 0.5*(self.t0 + self.t1)

        self._step = s
    
    
    def get_variable(self, name):

        cdf = self._opencdf()
        var = self._readvar(cdf, name)
            
        cdf.close()

        # Average over FBM time widow, if set
        if (var is not None) and (self.step is not None):
            t = var[0]
            t_window = (t>=self.t0) & (t<=self.t1)

            if len(var) == 2:
                var = np.mean(var[1][t_window])

            elif len(var) == 3:
                var = [var[1], np.mean(var[2][t_window],axis=0)]

        return var


    def get_rho(self, R, Z):
        """ Sqrt of normalized toroidal flux. """
        
        if not hasattr(self, 'X'):
            self.load_flux()
            
        return griddata((self.R, self.Z), self.X, (R,Z), fill_value=2)


    def load_flux(self, n_first_surf=5):

        cdf = self._opencdf()

        XB = self._readvar(cdf, 'XB')
        t = XB[0]
        XB = XB[2]
        nt = XB.shape[0]
        nx = XB.shape[1]
        n_surf = n_first_surf * (np.arange(nx)+1)          # number of points on each flux surface
        n_points = np.sum(n_surf)                       # total number of points

        R = np.zeros((nt, n_points))
        Z = np.zeros((nt, n_points))
        X = np.zeros((nt, n_points))
        theta = np.zeros((nt, n_points))

        RMC = cdf.variables['RMC00'][:]
        YMC = cdf.variables['YMC00'][:]
        RMS = np.zeros_like(RMC)
        YMS = np.zeros_like(YMC)
        i = 0
        more_moments = True

        while more_moments:

            # Loop over flux surfaces and compute R and Z for different theta
            i0 = 0
            for ix in range(nx):

                # Poloidal angle points
                i1 = i0 + n_surf[ix]
                th = np.linspace(0.0, 2*np.pi, n_surf[ix]+1)[:-1]     # exclude 2*pi...    

                # Evaluate R and Z coordinates
                R[:,i0:i1] = R[:,i0:i1] + np.outer(RMC[:,ix],np.cos(i*th)) + np.outer(RMS[:,ix],np.sin(i*th))
                Z[:,i0:i1] = Z[:,i0:i1] + np.outer(YMC[:,ix],np.cos(i*th)) + np.outer(YMS[:,ix],np.sin(i*th))

                X[:,i0:i1] = np.outer(XB[:,ix],np.ones(len(th)))
                theta[:,i0:i1] = np.outer(np.ones(nt),th)

                # Change start index for next surface
                i0 = i1 

            # Read next set of sin and cos moments, if available
            i += 1
            cos_name = 'MC{0:02d}'.format(i)
            sin_name = 'MS{0:02d}'.format(i)
            try:
                RMC = cdf.variables['R' + cos_name][:]
                YMC = cdf.variables['Y' + cos_name][:]
                RMS = cdf.variables['R' + sin_name][:]
                YMS = cdf.variables['Y' + sin_name][:]
            except:
                more_moments = False

        cdf.close()

        # Average over FBM time widow, if set
        if self.step is not None:
            t_window = (t>=self.t0) & (t<self.t1)

            R = np.mean(R[t_window], axis=0)
            Z = np.mean(Z[t_window], axis=0)
            X = np.mean(X[t_window], axis=0)
            theta = np.mean(theta[t_window], axis=0)
        
        self.R =  R/100.0 
        self.Z = Z/100.0
        self.X = X


    def _opencdf(self):
        return Dataset(self.out_file, mode='r')


    def _openfbm(self):

        if self.step is None:
            raise ValueError('No FBM step set!')

        i = self.step - 1
        return Dataset(self.fbm_files[i], mode='r')


    def _readvar(self, cdf, name):
        name = name.upper()

        try:
            cdfvar = cdf.variables[name]
        except:
            return None

        # Coordinate axes
        dims = cdfvar.dimensions

        # Case 1D
        if len(dims) == 1:
            t = cdf.variables[dims[0]][:]
            f = cdfvar[:]                

            var = [t,f]

        # Case 2D
        elif len(dims) == 2:
            t = cdf.variables[dims[0]][:]
            x = cdf.variables[dims[1]][:]
            f = cdfvar[:]

            x = x[0]    # assuming x bins never change over time...

            var = [t,x,f]

        return var
