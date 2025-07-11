# ---------------------------------------------------------------------------------------------
# Forward modelling of fusion product energy spectrum for a given set
# of Monte-Carlo samples (fuel ions), synthetic diagnostic viewing cone and bulk
# plasma distribution. Toroidal symmetry is assumed.
#
# Written by Jacob Eriksson, Andrea Valentini and Henrik Järleblad. Last maintained 2025-07-11.
# ---------------------------------------------------------------------------------------------


import numpy as np

import spec
import tokapart
import vcone
import constants

class Forward(object):
    """
    A class representing a forward model for fusion product energy spectrum diagnostics in toroidally symmetric
    magnetic confinement fusion reactors. Please note, compared to the old version of the forward.py script,
    this new approach is used simply to be able to pre-load and prepare as much as possible, prior to using the 
    forward model in a distributed (parallel computing) for-loop. The new version is also cleaner.

    To compute a fusion product energy spectrum with the Forward class, the user can either:
    - Call the __init__() method, followed by the Forward.calc() method (old approach)
    - Call the __init__() method, followed by the Forward.add_gyration() method, manually modify the 
      outputs as needed and then run the Forward.compute_spectrum() method. This new approach gives the 
      user more control over what happens in between each step in the Forward.calc() method.
    """

    def getReactionForm(self, fusion_reaction):
        """
        A very short, but necessary, Python implementation of a function from the
        OWCF/misc/availReacts.jl script. This is to avoid having to import a Julia
        function into Python.
        """
        if "-" in fusion_reaction:
            return 2
        if "(" in fusion_reaction:
            return 1
        return 3

    def __init__(self, viewing_cone_name, fusion_reaction, bulk_dist):
        """
        Initialize the Forward object. The inputs are as follows:
        - viewing_cone_name - A String. Should be either any of the valid synthetic diagnostic names in OWCF/vcone.py 
                              or the path to an output file of the OWCF/extra/createCustomLOS.jl script. The path to
                              an output file of the LINE21 code is also accepted - String
        - fusion_reaction - A String. The fusion reaction String should be specified using one of the following forms:
                                (1) 'a(b,c)d' 
                                (2) 'a(b,c)d-l' 
                                (3) 'b' 
                            where a is thermal ion, b is fast ion, c is fusion product particle of interest, 
                            d is fusion product particle of disinterest and l is the nuclear energy state of c. 
                            l can be GS, 1L or 2L, corresponding to Ground State (GS), 1st excited energy level (1L) 
                            and 2nd excited energy level (2L). For lists of available fusion reactions and particle species, 
                            please see OWCF/misc/availReacts.jl and OWCF/misc/species_func.jl. The reaction forms imply:
                                (1) The standard fusion reaction form. The nuclear energy level 'l' of the fusion product 
                                    particle of interest 'c' is automatically assumed to be GS (if relevant).
                                (2) The advanced fusion reaction form. The nuclear energy level 'l' of the fusion product 
                                    particle of interest 'c' is taken to be 'l'.
                                (3) The projected velocity reaction form. No fusion reaction is computed. Instead, the 
                                    (energy) spectrum is returned as a velocity spectrum from the velocity vectors of 
                                    the ion 'b', projected onto the diagnostic line-of-sight (LOS), using the (E,p,R,z) 
                                    points that are inside the LOS.
                            PLEASE NOTE! Specify alpha particles as '4he' or '4He' (NOT 'he4' or 'He4'). The same goes 
                            for helium-3 (specify as '3he', NOT 'he3'). Etc.
        - bulk_dist - A String equal to "" or a Thermal object as defined in OWCF/transp_dists.py
        """
        if not viewing_cone_name=="":
            self.viewing_cone = vcone.ViewingCone(viewing_cone_name)
        else:
            self.viewing_cone = ""

        self.reaction = fusion_reaction

        self.bulk_dist = bulk_dist

    @property
    def reaction(self):
        return self._reaction

    @reaction.setter
    def reaction(self, r):
        self.reaction_form = self.getReactionForm(r)

        if self.reaction_form==1:
            self._reaction = r
            self.product_of_interest_energy_state = 'GS' # Ground state
            reactants, products = r.split(',')
            self.thermal_ion, self.fast_ion = reactants.split('(')
            self.product_of_interest, self.product_of_disinterest = products.split(')')
            a = self.fast_ion
            b = self.thermal_ion
            self.projVel = False
        elif self.reaction_form==2:
            self._reaction, self.product_of_interest_energy_state = r.split('-')
            reactants, products = self._reaction.split(',')
            self.thermal_ion, self.fast_ion = reactants.split('(')
            self.product_of_interest, self.product_of_disinterest = products.split(')')
            a = self.fast_ion
            b = self.thermal_ion
            self.projVel = False
        else:
            self._reaction = r
            self.product_of_interest_energy_state = "" # Ground state
            self.thermal_ion = ""
            self.fast_ion = r
            self.product_of_interest = ""
            self.product_of_disinterest = ""
            a = r
            b = 'proj'
            self.projVel = True
        self.reaction_short = a+"-"+b # Only the reactants, on the form [fast ion]-[thermal ion]. If projected velocities: [fast ion]-proj

        if self.projVel:
            self.spectrum_calculator = ""
        else:
            self.spectrum_calculator = spec.SpectrumCalculator(self.reaction_short, self.product_of_interest_energy_state)

    @property
    def bulk_dist(self):
        return self._bulk_dist
    
    @bulk_dist.setter
    def bulk_dist(self, bd):
        self._bulk_dist = bd
        if not bd=="":
            self.bulk_dist_name = "TRANSP"
        else:
            self.bulk_dist_name = "Custom"

    def __repr__(self):
        if not self.viewing_cone =="":
            return 'Forward model with \n ---> {} \n ---> Fusion reaction: {} \n ---> Bulk plasma distribution: {}'.format(self.viewing_cone, self.reaction, self.bulk_dist_name)
        else:
            return 'Forward model with \n ---> Spherical 4*pi emission \n ---> Fusion reaction: {} \n ---> Bulk plasma distribution: {}'.format(self.reaction, self.bulk_dist_name)


    def calc(self, E, p, R, z, weights, bulk_dist, Ed_bins, B_vec, n_repeat=1, reaction='D(D,n)3He', 
             bulk_temp='', bulk_dens='', flr=False, v_rot=np.zeros((3,1))):
        """
        Compute the (total) expected energy spectrum of fusion product particles of interest, given
        the fusion reaction specified by the 'reaction' input variable and the diagnostic 
        viewing cone set via the set_view() function above. The fusion reaction is
        specified via the 'reaction' input variable. The 'reaction' input variable should be 
        specified using one of the following forms:
            (1) 'a(b,c)d' 
            (2) 'a(b,c)d-l' 
            (3) 'b' 
        where a is thermal ion, b is fast ion, c is fusion product particle of interest, 
        d is fusion product particle of disinterest and l is the nuclear energy state of c. 
        l can be GS, 1L or 2L, corresponding to Ground State (GS), 1st excited energy level (1L) 
        and 2nd excited energy level (2L). For lists of available fusion reactions and particle species, 
        please see OWCF/misc/availReacts.jl and OWCF/misc/species_func.jl. The reaction forms imply:
            (1) The standard fusion reaction form. The nuclear energy level 'l' of the fusion product 
                particle of interest 'c' is automatically assumed to be GS (if relevant).
            (2) The advanced fusion reaction form. The nuclear energy level 'l' of the fusion product 
                particle of interest 'c' is taken to be 'l'.
            (3) The projected velocity reaction form. No fusion reaction is computed. Instead, the 
                (energy) spectrum is returned as a velocity spectrum from the velocity vectors of 
                the ion 'b', projected onto the diagnostic line-of-sight (LOS), using the (E,p,R,z) 
                points that are inside the LOS.
        PLEASE NOTE! Specify alpha particles as '4he' or '4He' (NOT 'he4' or 'He4'). The same goes 
        for helium-3 (specify as '3he', NOT 'he3'). Etc.

        The E, p, R and z inputs are the energy (E)(in keV), pitch (p), major radius (R)(in meters)
        and vertical position (z)(in meters) coordinates of the fast ions. len(E)==len(p)==
        len(R)==len(z)==len(weights) must hold. PLEASE NOTE! The E, p, R and z inputs are NOT
        grid points, but points! I.e. the total energy spectrum of fusion product particles of
        interest is computed from the combined spectra produced by each (E[i],p[i],R[i],z[i])
        point for i in range(len(E)).

        The array weights ´weights´ contains the statistical weights of each MC sample.
        The sum of the weights should be equal to the total number of particles in the
        sampled distribution.

        The thermal ion temperature and density distributions are specified via the bulk_dist
        input object. This is an object of the Thermal class, defined in OWCF/transp_dists.py.
        The bulk_dist input object can also be specified as ''. The 'bulk_temp' and 'bulk_dens'
        inputs will then be used to model the thermal ion temperature and density (see below).

        Ed_bins specifies the diagnostic measurement bin edges (in keV).

        B_vec is a (3,N) array containing the (R,phi,z) components of the magnetic field 
        vector at all of the N (E,p,R,z) points (in Teslas).

        n_repeat is the number of gyro angles sampled for each (E,p) point.

        bulk_temp will be used as the (mean) Maxwellian bulk plasma temperature, if bulk_dist 
        is not specified (''). It should be an array equal in length to E,p,R,z and weights.

        bulk_dens will be used as the bulk plasma density, if bulk_dist is not specified ('').
        It should be an array equal in length to E,p,R,z and weights.

        flr should be set to true if finite Larmor radius effects should be taken into account.
        The computation will be more costly, but signals will be more realistic.

        v_rot is a (3,N) array containing the plasma rotation velocity vector (in m/s) in (R,phi,z) coordinates at all of the 
        N (E,p,R,z) points. If left unspecified (or np.sum(v_rot)==0), np.zeros(np.shape(E)) will be used.

        All input arrays are assumed to be numpy arrays. 
        """

        if (np.sum(v_rot) == 0.0) and np.shape(v_rot)[1]==1: # If not specified
            v_rot = v_rot.repeat(len(E), axis=1)

        reaction = reaction.lower() # Turn reaction into lowercase
        E = np.array([max(0.0, x) for x in E]) # Enforce a minimum of 0 for the energy. In case some sampling algorithm has sampled energy values slightly below 0. 2025-06-05
        p = np.array([max(-1.0,x) for x in [min(1.0,x) for x in p]]) # Clamping the pitch between -1.0 and 1.0. In case some sampling algorithm has sampled pitch values slightly outside of the [-1.0,1.0] range. 2025-06-05

        # Check fusion reaction input. Deduce thermal ion, fast ion and the nucleus energy level of the fusion product particle of interest.
        reaction_form = self.getReactionForm(reaction)

        if reaction_form==1:
            energy_state = 'GS'
            reactants, products = reaction.split(',')
            thermal_ion, fast_ion = reactants.split('(')
            a = fast_ion
            b = thermal_ion
        elif reaction_form==2:
            fusion_reaction, energy_state = reaction.split('-')
            reactants, products = fusion_reaction.split(',')
            thermal_ion, fast_ion = reactants.split('(')
            a = fast_ion
            b = thermal_ion
        else:
            energy_state = 'GS'
            a = reaction
            b = 'proj'
        reaction_short = a+"-"+b # Only the reactants, on the form [fast ion]-[thermal ion]

        projVel = False # By default, assume that the user does not wish to compute projected velocities
        if b == 'proj': # But if b is 'proj'...
            projVel = True # ... Then the calc() function should simply compute, bin and return a spectrum of binned projected velocities

        if bulk_dist=='' and (not projVel):
                if (not len(bulk_temp)==len(E)) or (not len(bulk_dens)==len(E)):
                    raise Exception("bulk_temp and/or bulk_dens was not equal in length to the E,p,R,z and weights input. Please correct and re-try.")
        # Repeat inputs, to prepare for gyro-angle sampling
        E = E.repeat(n_repeat) # Repeat n_repeat number of times, for gyro-angle sampling
        p = p.repeat(n_repeat) # ...
        R = R.repeat(n_repeat)
        z = z.repeat(n_repeat)
        weights = (1/n_repeat) * weights.repeat(n_repeat) # Re-scale weights
        B_vec = B_vec.repeat(n_repeat, axis=1)
        if bulk_dist=='' and (not projVel):
            bulk_temp = bulk_temp.repeat(n_repeat)
            bulk_dens = bulk_dens.repeat(n_repeat)
        v_rot = v_rot.repeat(n_repeat, axis=1) # Repeat n_repeat number of times, for gyro-angle sampling        
        spec_calc = spec.SpectrumCalculator(reaction=reaction_short,product_state=energy_state)
        m = spec_calc.ma # The mass of the fast ion (in keV/c**2)

        # Expand to gyro-radius points, if FLR effects are to be included
        if flr:
            q_ecu = (spec_calc.qa)/(constants.e) # The charge of the fast ion, in elemental charge units
            v, x = tokapart.add_gyration(E, p, m, B_vec, flr=flr, Z=q_ecu, Rg=R, zg=z) # Add gyro-motion to the guiding centre, via randomness. Need mass in keV/c**2 Return actual gyro-orbit positions (x)
            R = np.squeeze(x[0,:]) # New R-coordinates, corresponding to gyro-radius points
            z = np.squeeze(x[2,:]) # New z-coordinates, corresponding to gyro-radius points
        else:
            v = tokapart.add_gyration(E, p, m, B_vec) # Add gyro-motion to the guiding centre, via randomness. Don't return actual gyro-orbit positions (x)

        v = v + v_rot # Add plasma rotation to fast-ion velocities

        # Map (E,p,R,z) points to viewing cone voxels, if viewing cone has been specified
        if self.viewing_cone == "":
            R_vc = R # No viewing cone specified. All R points are accepted (/inside the universe)
            z_vc = z # No viewing cone specified. All z points are accepted (/inside the universe)
            v_vc = v # No viewing cone specified. All v vectors are accepted (/inside the universe)
            v_rot_vc = v_rot # No viewing cone specified. All v_rot vectors are accepted (/inside the universe)
            weights_vc = weights # No viewing cone specified. All weights are accepted
            bulk_temp_vc = bulk_temp
            bulk_dens_vc = bulk_dens
            B_vec_vc = B_vec # No viewing cone specified. All magnetic field vectors are accepted
            omega = 4*np.pi # Spherical emission, in all directions
            dphi = 2*np.pi # All toroidal angles are included
        else:
            i_voxels, i_points = self.viewing_cone.map_points_voxels(R, z) # Map the R,z sample points to the viewing cone voxels
            if len(i_points) == 0:
                #print('No points inside viewing cone!')
                return np.zeros(len(Ed_bins)-1)
            R_vc = R[i_points] # Extract R points within diagnostic viewing cone
            z_vc = z[i_points] # Extract z points within diagnostic viewing cone
            v_vc = v[:,i_points] # Extract v vectors with origins within diagnostic viewing cone
            v_rot_vc = v_rot # Just for script functionality
            if not (np.sum(v_rot) == 0.0):
                v_rot_vc = v_rot[:,i_points] # Extract v_rot vectors with origins within diagnostic viewing cone
            weights_vc = weights[i_points] # Extract weights corresponding to R,z points within diagnostic viewing cone
            if bulk_dist=='' and (not projVel):
                bulk_temp_vc = bulk_temp[i_points]
                bulk_dens_vc = bulk_dens[i_points]
            B_vec_vc = B_vec[:,i_points] # Extract magnetic field vectors corresponding to R,z points within diagnostic viewing cone
            omega = self.viewing_cone.OMEGA[i_voxels] # Only specific solid angles, in steradian
            dphi = self.viewing_cone.dPHI # What is the incremental toroidal angle that one diagnostic viewing cone voxel occupies?
            u1 = self.viewing_cone.U[:,i_voxels] # What is the emission direction of the viewing cone?
            spec_calc.u1 = u1 # Specify it for the spectrum calculator. Otherwise, 4*pi emission will be assumed
        
        if projVel: # If we simply want the projected fast-ion velocities...
            spec_weights = weights_vc * omega * (dphi/(2*np.pi)) # Use weighting without bulk density
            if spec_calc.u1 is None:
                v_vc = np.linalg.norm(v_vc, axis=0) # Just compute the speed of all velocity vectors
                return np.histogram(v_vc, bins=Ed_bins, weights=spec_weights)[0] # Bin all velocities
            else:
                u1 = np.array(spec_calc.u1) # The vector parallel to the diagnostic sightline (pointing towards the detector)
                if u1.ndim == 1:
                    u1 = np.array(u1).reshape(3,1)   # same direction towards the detector for all samples (could be the case for a really small stagnation orbits for example)
                u1_norm = np.linalg.norm(u1,axis=0)
                proj_speeds = np.einsum("ij,ij->j",u1/u1_norm,v_vc) # Compute the projected velocities
                return np.histogram(proj_speeds, bins=Ed_bins, weights=spec_weights)[0] # Bin all projected velocities

        # Set necessary spectrum computation variable values
        spec_calc.n_samples = len(weights_vc) # The actual number of samples to be counted towards the MC simulation
        spec_calc.reactant_a.v = v_vc # The fast ion velocities (inside the viewing cone)
        spec_calc.B_dir = B_vec_vc # The exact magnetic field at all points of interest

        # Set bulk temperature and density
        if not bulk_dist=='': # If there is a bulk distribution specified
            T_bulk = bulk_dist.get_temperature(R_vc, z_vc) # Load its bulk temperature for the R,z points
            n_bulk = bulk_dist.get_density(R_vc, z_vc) # Load its density for the R,z points
        else: # If there is not a bulk distribution specified
            T_bulk = bulk_temp_vc # It must have been manually specified
            n_bulk = bulk_dens_vc # It must have been manually specified

        spec_calc.reactant_b.sample_maxwellian_dist(T_bulk, v_rot=v_rot_vc)

        spec_calc.weights = weights_vc * omega * n_bulk * (dphi/(2*np.pi))

        s = spec_calc(bins=Ed_bins)

        return s
    
    def add_gyration(self, E, p, R, z, B_vec, particle_species=None, n_gyro=20):
        """
        Add gyration to the guiding-centre points defined by the (E,p,R,z) input vectors.
        For a given index 'i', the corresponding (E,p,R,z) point is given by (E[i], p[i], R[i], z[i]).
        For every (E,p,R,z) point, sample n_gyro number of gyro-orbit points. Return both the 
        velocity vectors v at the gyro-orbit points x, and the gyro-orbit points themselves (x).

        The input variables are as follows:
        - E - The energy of the guiding-centre points. In keV - Vector (N,)
        - p - The pitch (v_||/v) of the guiding-centre points - Vector (N,)
        - R - The major radius coordinate of the guiding-centre points. In meters - Vector (N,)
        - z - The vertical coordinates of the guiding-centre points. In meters - Vector (N,)
        - B_vec - The magnetic field at the guiding-centre points. Array with shape (3,N) where
                  N is the number of (E,p,R,z) points. In Tesla - Array
        - particle_species - The particle species of interest. See OWCF/spec.py/Particle for more info.
                             If left unspecified, use self.fast_ion - String
        - n_gyro - The guiding-centre (E,p,R,z) point will be sampled into this many gyro-orbit points - Int

        The output variable are as follows:
        - v - The velocity vectors (v_R, v_phi, v_z) of the particle at x (see below). In meters/second - Array (3, N*n_gyro)
        - x - The gyro-orbit points (R, phi, z) of the particle. In (meters, radians, meters) - Array (3, N*n_gyro)
        """

        if particle_species==None:
            particle_species = self.fast_ion

        E = np.array([max(0.0, x) for x in E]) # Enforce a minimum of 0 for the energy. In case some sampling algorithm has sampled energy values slightly below 0. 2025-06-05
        p = np.array([max(-1.0,x) for x in [min(1.0,x) for x in p]]) # Clamping the pitch between -1.0 and 1.0. In case some sampling algorithm has sampled pitch values slightly outside of the [-1.0,1.0] range. 2025-06-05

        # Repeat inputs, to prepare for gyro-angle sampling
        E = E.repeat(n_gyro) # Repeat n_gyro number of times, for gyro-angle sampling
        p = p.repeat(n_gyro) # ...
        R = R.repeat(n_gyro)
        z = z.repeat(n_gyro)
        B_vec = B_vec.repeat(n_gyro, axis=1)
        prtcl = spec.Particle(particle_species)

        q_ecu = prtcl.q/(constants.e) # The charge of the particle, in elemental charge units
        v, x = tokapart.add_gyration(E, p, prtcl.m, B_vec, flr=True, Z=q_ecu, Rg=R, zg=z) # Add gyro-motion to the guiding centre, via randomness. Need mass in keV/c**2 Return actual gyro-orbit positions (x)

        return v, x
    
    def get_indices_of_points_in_viewing_cone(self, R, z, unique=False):
        """
        Given the (R,z) points as input, get the indices of the points that are inside the poloidal projection of the self.viewing_cone object.
        Please note! len(R)==len(z) must hold.

        The input variables are as follows:

        - R - The major radius coordinate points of interest. In meters - Vector (N,)
        - z - The vertical coordinate points of interest. In meters - Vector (N,)

        The input keyword arguments are as follows:

        - unique - If set to true, only unique indices will be returned as output. Please see the 'i_points' documentation below for more info.

        The output variable is as follows:
        - i_points - The indices of the input (R,z) points that were found to be inside the poloidal projection of the self.viewing_cone object. 
                     Please note! len(i_points) > len(R) might be true. Since each (R,z) points might be mapped to several viewing cone voxels.
                     If unique indices are wanted, set the input keyword argument 'unique' to true - Vector (M,)
        """

        if not len(R)==len(z):
            raise ValueError("Input R and z were not equal in shape (length). Please correct and re-try!")

        if self.viewing_cone == "":
            return range(len(R))
        else:
            i_voxels, i_points = self.viewing_cone.map_points_voxels(R, z) # Map the R,z points to the voxels of the viewing cone model
            # PLEASE NOTE! len(i_voxels)==len(i_points) is always true. len(i_points) > len(R) is likely to be true. This is because a single (R,z) point is likely to be mapped to several voxels
            if len(i_points) == 0:
                return range(0)
            if unique:
                return np.unique(i_points)
            return i_points
        
    def compute_spectrum(self, Ed_bins, R, z, v, weights, B_vecs, bulk_temp=np.array([-1.0]), bulk_dens=np.array([-1.0]), v_rot=np.zeros((3,1)), verbose=False):
        """
        Given the diagnostic viewing cone, fusion reaction and bulk (thermal) plasma distribution already set via the Forward.__init__() method,
        compute the synthetic diagnostic measurement spectrum using the following inputs:

        - Ed_bins - The edges of the bins into which to bin the computed measurements. The number of edges is one greater than the number of bins - Vector (nEd+1,)
        - R - The major radius coordinates to be used as the origin of the synthetic measurements - Vector (N,)
        - z - The vertical coordinates to be used as the origin -||- - Vector (N,)
        - v - The velocity vectors of the fast ion (one of the fusion reactants) - Array (3,N)
        - weights - The weights to be used for binning the computed measurements - Vector (N,)
        - B_vecs - The magnetic field at the input (R,z) points - Array (3,N)

        The keyword argument inputs are as follows:
        
        - bulk_temp - If no bulk (thermal) plasma distribution has been set via the Forward.__init__() method (i.e. it is ""), then the bulk plasma temperature 
                      at the input (R,z) points needs to be specified manually, via the 'bulk_temp' keyword argument input. In keV - Vector (N,)
        - bulk_dens - If no bulk (thermal) plasma distribution has been -||-, then the bulk plasma density at the input (R,z) points -||-, via the 'bulk_dens'
                      keyword argument input. In m^-3 - Vector (N,)
        - v_rot - The plasma rotation. In meters/second - Array (3,N)
        - verbose - If set to true, the function will be talkative! - Boolean
        """

        if not len(R)==len(z):
            raise ValueError("Input R and z were not equal in shape (length). Please correct and re-try!")

        if not np.shape(v)[1]==len(R):
            raise ValueError("Number of samples in input v is {} but number of samples in input R is {}. To compute a synthetic diagnostic spectrum, shape(v)[1] must match len(R). Please correct and re-try!".format(np.shape(v)[1], len(R)))
        
        if not np.shape(v)[1]==len(weights):
            raise ValueError("Number of samples in input v is {} but number of input weights is {}")

        if not np.shape(v)==np.shape(B_vecs):
            raise ValueError("Shape of input v is {} but shape of input B_vecs is {}. To compute synthetic diagnostic spectrum, shapes must match. Please correct and re-try!".format(np.shape(v), np.shape(B_vecs)))
        
        if (not self.bulk_dist_name == "TRANSP") and (not self.projVel):
            if (not len(bulk_temp)==len(R)) or np.sum(bulk_temp)==-1.0:
                raise ValueError("Forward object set to use custom bulk plasma tempemperature values. But input bulk_temp was not specified (correctly). Please correct and re-try!")
            
            if (not len(bulk_dens)==len(R)) or np.sum(bulk_dens)==-1.0:
                raise ValueError("Forward object set to use custom bulk plasma density values. But input bulk_dens was not specified (correctly). Please correct and re-try!")

        # Add the plasma rotation to the velocity vectors (if valid)
        if not np.sum(v_rot)==0.0: # If the plasma rotation is specified
            if not np.shape(v)==np.shape(v_rot):
                raise ValueError("The shape of plasma rotation input keyword argument 'v_rot' was not equal to the shape of input 'v'. Please correct and re-try.")
        else: # If the plasma rotation is left unspecified
            v_rot = np.zeros(np.shape(v)) # Just set zero plasma rotation everywhere
        
        if self.viewing_cone == "":
            # No viewing cone specified. All input points are accepted
            omega = 4*np.pi
            dphi = 2*np.pi
            self.spectrum_calculator.u1 = None # No viewing cone, so no vector pointing towards any detector
            self.spectrum_calculator.n_samples = len(R) # The number of samples to be counted towards the MC simulation
        else:
            i_voxels, i_points = self.viewing_cone.map_points_voxels(R, z) # Map the R,z points to the voxels of the viewing cone model
            # PLEASE NOTE! len(i_voxels)==len(i_points) is always true. But len(i_points) > len(R) is likely to be true. This is because a single (R,z) point is likely to be mapped to several voxels
            if len(i_points) == 0:
                return np.zeros(len(Ed_bins)-1) # No points in viewing cone
            R = R[i_points] # The R points inside the diagnostic viewing cone
            z = z[i_points] # The z points -||-
            v = v[:,i_points] # The v vectors -||-
            weights = weights[i_points] # The weights -||-
            B_vecs = B_vecs[:,i_points] # The B_vecs -||-
            if (not self.bulk_dist_name == "TRANSP") and (not self.projVel):
                bulk_temp = bulk_temp[i_points]
                bulk_dens = bulk_dens[i_points]
            v_rot = v_rot[:,i_points] # The v_rot vectors -||-
            omega = self.viewing_cone.OMEGA[i_voxels] # The solid angle values, in steradian, for the (R,z) points inside the viewing cone
            dphi = self.viewing_cone.dPHI # The incremental phi values of the (R,phi,z) grid used to create the model of the viewing cone
            self.spectrum_calculator.u1 = self.viewing_cone.U[:,i_voxels] # The direction towards the detector at the (R,z) points inside the viewing cone
            self.spectrum_calculator.n_samples = len(i_points) # The actual number of samples to be counted towards the MC simulation

        # Add plasma rotation (can be zero)
        v = v + v_rot

        if self.projVel:
            verbose and print("compute_spectrum(): self.projVel==True. Computing projected velocities... ")
            spec_weights = weights * self.omega * (self.dphi/(2*np.pi)) # Use weighting without bulk density
            if self.spectrum_calculator.u1 is None: # If there is no synthetic diagnostic viewing cone
                verbose and print("compute_spectrum(): No diagnostic viewing cone has been set. Binning speed of all input velocities... ")
                v = np.linalg.norm(v, axis=0) # Just compute the speed of all velocity vectors
                return np.histogram(v, bins=Ed_bins, weights=spec_weights)[0] # Bin all velocities
            else:
                u1 = np.array(self.spectrum_calculator.u1) # The vector parallel to the diagnostic sightline (pointing towards the detector)
                if u1.ndim == 1:
                    u1 = np.array(u1).reshape(3,1)   # same direction towards the detector for all samples (could be the case for a really small stagnation orbits for example)
                u1_norm = np.linalg.norm(u1,axis=0)
                proj_speeds = np.einsum("ij,ij->j",u1/u1_norm,v) # Compute the projected velocities
                return np.histogram(proj_speeds, bins=Ed_bins, weights=spec_weights)[0] # Bin all projected velocities
            
        # If not self.projVel, i.e. compute energy spectrum for fusion product particle of interest
        self.spectrum_calculator.reactant_a.v = v # The fast ion velocities (inside the viewing cone)
        self.spectrum_calculator.B_dir = B_vecs # The exact magnetic field at all points of interest

        if self.bulk_dist_name == "TRANSP":
            T_bulk = self.bulk_dist.get_temperature(R, z) # Load the bulk temperature for the R,z points
            n_bulk = self.bulk_dist.get_density(R, z) # Load the bulk density for the R,z points
        else:
            T_bulk = bulk_temp # The bulk plasma temperature for the R,z points, in keV
            n_bulk = bulk_dens # The bulk plasma temperature for the R,z points, in m^-3
        self.spectrum_calculator.reactant_b.sample_maxwellian_dist(T_bulk, v_rot=v_rot) # Sample velocities from a Maxwellian with mean T_bulk for thermal plasma distribution
        self.spectrum_calculator.weights = weights * omega * n_bulk * (dphi/(2*np.pi)) # The total Monte-Carlo weights, including the solid angle, bulk plasma density and fractional toroidal angle

        return self.spectrum_calculator(bins=Ed_bins) # Compute and return the computed spectrum
            