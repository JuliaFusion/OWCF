# Forward modelling of neutron energy spectrom for a given set
# of Monte-Carlo samples (fuel ions), and a given viewing cone.
# -------------------------------------------------------------


#from turtle import xcor
import numpy as np

import spec
import tokapart
import vcone

class Forward(object):
    """
    A class representing a forward model. Please note, compared to the old version of the forward.py script
    (where there was no Forward class), this approach is used simply to be able to pre-load and prepare
    as much as possible, prior to using the forward model in a distributed (parallel computing) for-loop
    in Julia, and sending it as input into a Julia function.
    """

    def __init__(self, viewing_cone_name):
        if not viewing_cone_name=="":
            self.viewing_cone = vcone.ViewingCone(viewing_cone_name)
        else:
            self.viewing_cone = ""

    def __repr__(self):
        if not self.viewing_cone =="":
            return 'Forward model with {}'.format(self.viewing_cone)
        else:
            return 'Forward model with spherical 4*pi emission'


    def calc(self, E, p, R, Z, weights, bulk_dist, Ed_bins, B_vec,
            n_repeat=1, reaction='d-d', bulk_temp='', bulk_dens=''):
        """
        Calculate spectrum for fast ions (specified by E, p, R, Z)
        reacting with the given bulk distribution, for a given viewing cone
        (which is set with the set_view function above).

        The array weights ´weights´ contains the statistical weights of each MC sample.
        The sum of the weights should be equal to the total number of particles in the
        sampled distribution.

        Ed_bins specifies the neutron bin edges (keV).

        B_vec is an (3,N) array containing the components of the (R,phi,Z) components
        of the magnetic field vector at all the N particle locations.

        n_repeat is the number of gyro angles sampled for each energy-pitch value.

        reaction is the fusion reaction to simulate ('D-D', 'D-T' etc). The species the the left of the '-' is always the 
        thermal species and the species to the right of the '-' is always the fast ion. However, as can be seen in the code,
        the order is flipped once it's been provided as input to the calc() function. 
        
        reaction = reaction[::-1]
        
        This is because the DRESS Python code framework assumes the order to be the other way around. 

        bulk_temp will be used as the bulk plasma temperature, if bulk_dist is not specified.
        It should be an array equal in length to E,p,R,Z and weights.

        bulk_dens will be used as the bulk plasma density, if bulk_dist is not specified.
        It should be an array equal in length to E,p,R,Z and weights.
        """
        
        # Here, the ordering in the input variable reaction is a-b, where 'a' is the thermal species and 'b' is the fast-ion species.
        a,b = reaction.split('-') # This annoying scheme needs to be performed because the OWCF uses the ordering 'thermal-fast' while the DRESS code uses 'fast-thermal'
        reaction = b+"-"+a # So, D-T becomes T-D. And proj-D becomes D-proj. So here, the ordering has been reversed.
        a,b = reaction.split('-') # So here, D-T becomes T-D. And proj-D becomes D-proj. 'a' will be the fast-ion species, and 'b' will be the thermal species.

        calcProjVel = False # By default, assume that the user does not wish to compute projected velocities
        if b == 'proj': # But if b is 'proj'...
            calcProjVel = True # ... Then the calc() function should simply compute, bin and return a spectrum of binned projected velocities

        if bulk_dist=='' and (not calcProjVel):
                if (not len(bulk_temp)==len(E)) or (not len(bulk_dens)==len(E)):
                    raise Exception("bulk_temp and/or bulk_dens was not equal in length to the E,p,R,Z and weights input. Please correct and re-try.")
        E = E.repeat(n_repeat)
        p = p.repeat(n_repeat)
        R = R.repeat(n_repeat)
        Z = Z.repeat(n_repeat)
        weights = weights.repeat(n_repeat)
        B_vec = B_vec.repeat(n_repeat, axis=1)
        if bulk_dist=='' and (not calcProjVel):
            bulk_temp = bulk_temp.repeat(n_repeat)
            bulk_dens = bulk_dens.repeat(n_repeat)
        
        spec_calc = spec.SpectrumCalculator(reaction=reaction)

        # Map orbit points to viewing cone voxels, if viewing cone has been specified
        if self.viewing_cone == "":
            E_vc = E # No viewing cone specified. All energy points are accepted (/inside the universe)
            p_vc = p # No viewing cone specified. All pitch points are accepted (/inside the universe)
            R_vc = R # No viewing cone specified. All R points are accepted (/inside the universe)
            Z_vc = Z # No viewing cone specified. All Z points are accepted (/inside the universe)
            weights_vc = weights # No viewing cone specified. All weights are accepted
            bulk_temp_vc = bulk_temp
            bulk_dens_vc = bulk_dens
            B_vec_vc = B_vec # No viewing cone specified. All magnetic field vectors are accepted
            omega = 4*np.pi # Spherical emission, in all directions
            dphi = 2*np.pi # All toroidal angles are included
        else:
            i_voxels, i_points = self.viewing_cone.map_points_voxels(R, Z) # Map the R,z sample points to the viewing cone voxels
            if len(i_points) == 0:
                #print('No points inside viewing cone!')
                return np.zeros(len(Ed_bins)-1)

            E_vc = E[i_points] # Extract energy points within diagnostic viewing cone
            p_vc = p[i_points] # Extract pitch points within diagnostic viewing cone
            R_vc = R[i_points] # Extract R points within diagnostic viewing cone
            Z_vc = Z[i_points] # Extract Z points within diagnostic viewing cone
            weights_vc = weights[i_points] # Extract weights corresponding to R,Z points within diagnostic viewing cone
            if bulk_dist=='' and (not calcProjVel):
                bulk_temp_vc = bulk_temp[i_points]
                bulk_dens_vc = bulk_dens[i_points]
            B_vec_vc = B_vec[:,i_points] # Extract magnetic field vectors corresponding to R,Z points within diagnostic viewing cone
            omega = self.viewing_cone.OMEGA[i_voxels] # Only specific angles, in steradian
            dphi = self.viewing_cone.dPHI # What is the incremental toroidal angle that one diagnostic viewing cone voxel occupies?
            u1 = self.viewing_cone.U[:,i_voxels] # What is the emission direction of the viewing cone?
            spec_calc.u1 = u1 # Specify it for the spectrum calculator. Otherwise, 4*pi emission will be assumed

        m = spec_calc.ma # The mass of the fast ion
        v = tokapart.add_gyration(E_vc, p_vc, m, B_vec_vc) # Add gyro-motion to the guiding centre, via randomness

        if calcProjVel: # If we simply want the projected fast-ion velocities...
            spec_weights = weights_vc * omega * (dphi/(2*np.pi)) / n_repeat # Use weighting without bulk density
            if spec_calc.u1 is None:
                v = np.linalg.norm(v, axis=0) # Just compute the speed of all velocity vectors
                return np.histogram(v, bins=Ed_bins, weights=spec_weights)[0] # Bin all velocities
            else:
                u1 = np.array(spec_calc.u1) # The vector parallel to the diagnostic sightline (pointing towards the detector)
                if u1.ndim == 1:
                    u1 = np.array(u1).reshape(3,1)   # same direction towards the detector for all samples (could be the case for a really small stagnation orbits for example)
                u1_norm = np.linalg.norm(u1,axis=0)
                #proj_speeds = np.linalg.norm((u1*v)/u1_norm,axis=0) # The magnitude of the projection of the ion velocities onto the vector parallel to the diagnostic sightline
                proj_speeds = np.einsum("ij,ij->j",u1/u1_norm,v) # Compute the projected velocities
                return np.histogram(proj_speeds, bins=Ed_bins, weights=spec_weights)[0] # Bin all projected velocities


        # Compute spectrum
        spec_calc.n_samples = len(E_vc)
        spec_calc.reactant_a.v = v

        if not bulk_dist=='': # If there is a bulk distribution specified
            T_bulk = bulk_dist.get_temperature(R_vc, Z_vc) # Load its bulk temperature for the R,Z points
            n_bulk = bulk_dist.get_density(R_vc, Z_vc) # Load its density for the R,Z points
        else: # If there is not a bulk distribution specified
            T_bulk = bulk_temp_vc # It must have been manually specified
            n_bulk = bulk_dens_vc # It must have been manually specified

        spec_calc.reactant_b.sample_maxwellian_dist(T_bulk)

        spec_calc.weights = weights_vc * omega * n_bulk * (dphi/(2*np.pi)) / n_repeat
        #spec_calc.weights = weights_vc * R_vc * omega * n_bulk * dphi / n_repeat

        s = spec_calc(bins=Ed_bins)

        return s
