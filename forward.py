# Forward modelling of neutron energy spectrom for a given set
# of Monte-Carlo samples (fuel ions), and a given viewing cone.
# -------------------------------------------------------------


#from turtle import xcor
import numpy as np

import spec
import tokapart
import vcone
import constants

"""
A very short, but necessary, Python implementation of a function from the
OWCF/misc/availReacts.jl script. This is to avoid having to import a Julia
function into Python.
"""
def getReactionForm(fusion_reaction):
    if "-" in fusion_reaction:
        return 2
    if "(" in fusion_reaction:
        return 1
    return 3

class Forward(object):
    """
    A class representing a forward model. Please note, compared to the old version of the forward.py script
    (where there was no Forward class), this approach is used simply to be able to pre-load and prepare
    as much as possible, prior to using the forward model in a distributed (parallel computing) for-loop.
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


    def calc(self, E, p, R, z, weights, bulk_dist, Ed_bins, B_vec, n_repeat=1, reaction='D(T,n)4He', 
             bulk_temp='', bulk_dens='', flr=False, v_rot=np.reshape(np.array([0.0,0.0,0.0]),(3,1))):
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
        N (E,p,R,z) points.

        All input arrays are assumed to be numpy arrays. 
        """

        # Check fusion reaction input. Deduce thermal ion, fast ion and the nucleus energy level of the fusion product particle of interest.
        reaction_form = getReactionForm(reaction)

        if reaction_form==1:
            energy_state = 'GS'
            reactants, products = reaction.split(',')
            thermal_ion, fast_ion = reactants.split('(')
            a = fast_ion
            b = thermal_ion
        elif reaction_form==2:
            fusion_reaction, energy_state = reaction.split('-')
            reactants, products = fusion_reaction.split(',')
            thermal_ion, fast_ion = fusion_reaction.split('(')
            a = fast_ion
            b = thermal_ion
        else:
            energy_state = 'GS'
            a = fast_ion
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
        weights = weights.repeat(n_repeat)
        B_vec = B_vec.repeat(n_repeat, axis=1)
        if bulk_dist=='' and (not projVel):
            bulk_temp = bulk_temp.repeat(n_repeat)
            bulk_dens = bulk_dens.repeat(n_repeat)
        if not (np.sum(v_rot) == 0.0):
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

        if not (np.sum(v_rot) == 0.0):
            v = v + v_rot # Add plasma rotation to fast-ion velocities

        # Map orbit points to viewing cone voxels, if viewing cone has been specified
        if self.viewing_cone == "":
            #E_vc = E # No viewing cone specified. All energy points are accepted (/inside the universe)
            #p_vc = p # No viewing cone specified. All pitch points are accepted (/inside the universe)
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
            
            #E_vc = E[i_points] # Extract energy points within diagnostic viewing cone
            #p_vc = p[i_points] # Extract pitch points within diagnostic viewing cone
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
            spec_weights = weights_vc * omega * (dphi/(2*np.pi)) / n_repeat # Use weighting without bulk density
            if spec_calc.u1 is None:
                v_vc = np.linalg.norm(v_vc, axis=0) # Just compute the speed of all velocity vectors
                return np.histogram(v_vc, bins=Ed_bins, weights=spec_weights)[0] # Bin all velocities
            else:
                u1 = np.array(spec_calc.u1) # The vector parallel to the diagnostic sightline (pointing towards the detector)
                if u1.ndim == 1:
                    u1 = np.array(u1).reshape(3,1)   # same direction towards the detector for all samples (could be the case for a really small stagnation orbits for example)
                u1_norm = np.linalg.norm(u1,axis=0)
                #proj_speeds = np.linalg.norm((u1*v)/u1_norm,axis=0) # The magnitude of the projection of the ion velocities onto the vector parallel to the diagnostic sightline
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

        spec_calc.weights = weights_vc * omega * n_bulk * (dphi/(2*np.pi)) / n_repeat
        #spec_calc.weights = weights_vc * n_bulk * (dphi/(2*np.pi)) / n_repeat
        #spec_calc.weights = weights_vc * R_vc * omega * n_bulk * dphi / n_repeat

        s = spec_calc(bins=Ed_bins)

        return s
