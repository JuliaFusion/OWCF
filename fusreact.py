# Cross sections and reactivities.
# ----------------------------------

import os

import numpy as np
from scipy.special import legendre
import pandas as pd
from scipy.interpolate import interp1d
from scipy.interpolate import RegularGridInterpolator

import constants


# Load Legendre coefficients for DD differential cross section
# -----------------------------------------------------------------
dd_coeff_endf = np.loadtxt('fit_data/ddn3he_legendre_endf.txt').T
# Same, but for DT neutrons
dt_coeff_endf = pd.read_csv("fit_data/dtn4He_legendre_Drosg.txt",sep="\t",header=None)

# Convert from LAB to COM (equal masses) and energies to keV
dd_coeff_endf[0,:] = dd_coeff_endf[0,:]/2000

# Number of coefficients in the expansion
n_dd_coeff = dd_coeff_endf.shape[0] - 1

# Also load necessary data for alpha-9Be reactions
alpha9Be_1L_dataframe = pd.read_csv('fit_data/AlphaBe9_1L.txt', header=None)
alpha9Be_2L_dataframe = pd.read_csv('fit_data/AlphaBe9_2L.txt', header=None)

def sigma_tot(E, reaction='d-d'):
    """
    Returns the fusion cross section (in m**2) for a vector of CM energies (in keV).
    """

    E = np.atleast_1d(E).astype('d')
    sigma = np.zeros_like(E)

    reaction = reaction.lower()
    a,b = reaction.split('-')

    if reaction == 'd-d':
        # Constants for the D(d,n)3He reaction
        B_G = 31.3970 * np.ones_like(E)

        A1  = 5.3701e4
        A2  = 3.3027e2
        A3  = -1.2706e-1
        A4  = 2.9327e-5
        A5  = -2.5151e-9

        # Calculate the S factor
        S = A1+E*(A2+E*(A3+E*(A4+E*A5)))

    elif all(p in [a,b] for p in ['d','t']):
        # Constants for the T(d,n)4He reaction
        # (low and high energy parametrization)
        B_G  = 34.3827 * np.ones_like(E)
        A1_l = 6.927e4
        A2_l = 7.454e8
        A3_l = 2.05e6
        A4_l = 5.2002e4
        A5_l = 0.0
        B1_l = 6.38e1
        B2_l = -9.95e-1
        B3_l = 6.981e-5
        B4_l = 1.728e-4
        A1_h = -1.4714e6
        A2_h = 0.0
        A3_h = 0.0
        A4_h = 0.0
        A5_h = 0.0
        B1_h = -8.4127e-3
        B2_h = 4.7983e-6
        B3_h = -1.0748e-9
        B4_h = 8.5184e-14

        I_l = E<=550.0
        I_h = E>550.0
        El = E[I_l]
        Eh = E[I_h]

        # Calculate the S factor
        S = np.zeros(np.shape(E))
        S[I_l] = (A1_l + El*(A2_l + El*(A3_l + El*(A4_l + El*A5_l))))/\
                (1.0 + El*(B1_l + El*(B2_l + El*(B3_l + El*B4_l))))
        S[I_h] = (A1_h + Eh*(A2_h + Eh*(A3_h + Eh*(A4_h + Eh*A5_h))))/\
                (1.0 + Eh*(B1_h + Eh*(B2_h + Eh*(B3_h + Eh*B4_h))))

    elif all(p in [a,b] for p in ['d','3he']):
        # Constants for the 3He(d,p)4He reaction
        # (low and high energy parametrization)
        B_G  = 68.7508 * np.ones_like(E)
        A1_l = 5.7501e6
        A2_l = 2.5226e3
        A3_l = 4.5566e1
        A4_l = 0.0
        A5_l = 0.0
        B1_l = -3.1995e-3
        B2_l = -8.5530e-6
        B3_l = 5.9014e-8
        B4_l = 0.0
        A1_h = -8.3993e5
        A2_h = 0.0
        A3_h = 0.0
        A4_h = 0.0
        A5_h = 0.0
        B1_h = -2.6830e-3
        B2_h = 1.1633e-6
        B3_h = -2.1332e-10
        B4_h = 1.4250e-14

        I_l = E<=900.0
        I_h = E>900.0
        El = E[I_l]
        Eh = E[I_h]

        # Calculate the S factor
        S = np.zeros(np.shape(E))
        S[I_l] = (A1_l + El*(A2_l + El*(A3_l + El*(A4_l + El*A5_l))))/\
                (1.0 + El*(B1_l + El*(B2_l + El*(B3_l + El*B4_l))))
        S[I_h] = (A1_h + Eh*(A2_h + Eh*(A3_h + Eh*(A4_h + Eh*A5_h))))/\
                (1.0 + Eh*(B1_h + Eh*(B2_h + Eh*(B3_h + Eh*B4_h))))

    elif all(p in [a,b] for p in ['p','t']) or all(p in [a,b] for p in ['h','t']):
        # T(p,g)4He reaction. Parameterization provided by Massimo Nocente.

        # Convert to LAB energies in MeV
        Ecm = E / 1000.0                      # CM kinetic energy (MeV)
        mp = constants.mp * constants.u_keV / 1000.0
        mt = constants.mt * constants.u_keV / 1000.0
        E = Ecm + mp + mt                     # total CM energy (MeV)
        E = (E**2 - mp**2 - mt**2) / (2*mt)   # total LAB energy of the projectile (MeV)
        E = E - mp                            # LAB kinetic energy of the projectile (MeV)
        #E = (mp + mt) * Ecm / mt              # classical limit

        # Parameterization coefficients
        A0 = np.zeros_like(E)
        A1 = np.zeros_like(E)
        A2 = np.zeros_like(E)
        B1 = np.zeros_like(E)
        B2 = np.zeros_like(E)
        B_G = np.zeros_like(E)

        low = E<=0.5
        high = E>0.5

        A0[low] = 2.087
        A1[low] = 4.495e1
        A2[low] = 0.0
        B1[low] = -3.977e-1
        B2[low] = 0.0
        B_G[low] = 1.094

        A0[high] = 1.928e2
        A1[high] = -1.643e1
        A2[high] = 1.066e2
        B1[high] = -2.011e-1
        B2[high] = 9.213e-2
        B_G[high] = 2.510

        # Calculate the S factor
        S1 = A0 + A1*E + A2*E**2
        S2 = 1 + B1*E + B2*E**2
        S = S1 / S2 * 1e-3            # to make units mb later

        # Put cross section to zero if the energy is above the
        # upper limit of the parameterization
        out_of_range = E > 6.0
        S[out_of_range] = 0.0

    else:
        print('Unknown reaction')
        return

    # Calculate the cross section
    nonzero = (E > 0)
    sigma[nonzero] = S[nonzero]/E[nonzero]*np.exp(-B_G[nonzero]/np.sqrt(E[nonzero]))     # mb

    return sigma * 1e-31     # m**2


def sigma_diff(E, costheta, reaction='d-d', product_state='gs', anisotropic=True):
    """
    Returns the differential fusion cross section (in m**2/sr) for a vector of CM energies (in keV).
    """

    E = np.atleast_1d(E)
    costheta = np.atleast_1d(costheta)

    reaction = reaction.lower()
    a,b = reaction.split('-')

    if reaction == 'd-d' and anisotropic:
        # Compute the angular dependence of the DD reaction
        coeff = dd_coeff_endf
        A     = interp1d(coeff[0],coeff[1:], fill_value=0, bounds_error=False)(E)
        prob  = 0.5*np.ones_like(E)

        for i in range(n_dd_coeff):
            l = i+1
            P = legendre(l)(costheta)
            prob += (2*l+1)/2 * A[i] * P

        sigma = prob*sigma_tot(E, reaction=reaction)/(2*np.pi)

    elif all(p in [a,b] for p in ['p','t']) and anisotropic:
        # T(p,g)4He reaction. Parameterization provided by Massimo Nocente.

        # First convert to LAB energies in MeV
        Ecm = E

        mp = constants.mp * constants.u_keV / 1000.0
        mt = constants.mt * constants.u_keV / 1000.0
        E = E / 1000.0 + mp + mt              # total CM energy (MeV)
        E = (E**2 - mp**2 - mt**2) / (2*mt)   # total LAB energy of the projectile (MeV)
        E = E - mp                            # LAB kinetic energy of the projectile (MeV)

        # Compute angular dependence term
        a = 0.01 + 0.02*E
        sintheta = np.sqrt(1 - costheta**2)
        W = (sintheta + a*sintheta*costheta)**2
        C = 8.0*np.pi*(a**2 + 5.0) / 15.0

        # Compute cross section
        sigma = sigma_tot(Ecm, reaction=reaction) * W / C
        #sigma = sigma_tot(Ecm, reaction=reaction) / (4*np.pi)
        #sigma = 1.0
    elif all(p in [a,b] for p in ['4he','9be']) and product_state == '1l': 
        # 9Be(4He,ng)12C reaction. Cross-sectional patchwork provided by Massimo Nocente.
        # 1->0 transition of 12C
        df = alpha9Be_1L_dataframe

        # First extract quantities of interest from the dataframe
        dims = tuple(map(int,df[0][0].split(' '))) # tuple of relevant dimensions
        energies = np.asarray( df[0][1].split('   '), dtype=np.float64 ).ravel() # energies of CM in MeV
        cosines = df[0][2:].str.split('\t', expand=True).to_numpy(dtype=np.float64)[:,0].ravel() # cosine values of scattering angle in CM frame 
        diff_cross = df[0][2:].str.split('\t', expand=True).to_numpy(dtype=np.float64)[:,1:] # corresponding differential cross section values in mb/sr

        if not anisotropic: # sigma_tot = INT( 2*pi*diff_cross * dcos ), across all angles
            diff_cross = np.tile( np.einsum('ij,i->j',diff_cross,np.convolve(np.diff(cosines),[0.5,0.5])) / 2, (cosines.size,1) ) 

        # Define linear interpolator, given the patchwork
        # N.B: if the later inputted energy exceeds the upper limit, the interpolator returns zero values, across all angles, yet no warnings 
        interp = RegularGridInterpolator((energies,cosines), diff_cross.transpose(), bounds_error=False, fill_value=0)

        # Interpolate over requested values; zero values will give no contribution to the spectra
        E /= 1000 # MeV
        sigma = interp((E,-1*costheta)) * 1e-31 # m**2/sr; the -1 factor is due to the carbon-12 being followed, and not the neutron (for which we have the emission probability, technically)

    elif all(p in [a,b] for p in ['4he','9be']) and product_state == '2l': 
        # 9Be(4He,ng)12C reaction. Cross-sectional patchwork provided by Massimo Nocente.
        # 2->1 transition of 12C
        df = alpha9Be_2L_dataframe

        # First extract quantities of interest from the dataframe
        dims = tuple(map(int,df[0][0].split(' '))) # tuple of relevant dimensions
        energies = np.asarray( df[0][1].split('   '), dtype=np.float64 ).ravel() # energies of CM in MeV
        cosines = df[0][2:].str.split('\t', expand=True).to_numpy(dtype=np.float64)[:,0].ravel() # cosine values of scattering angle in CM frame 
        diff_cross = df[0][2:].str.split('\t', expand=True).to_numpy(dtype=np.float64)[:,1:] # corresponding differential cross section values in mb/sr

        if not anisotropic: # sigma_tot = INT( 2*pi*diff_cross * dcos ), across all angles
            diff_cross = np.tile( np.einsum('ij,i->j',diff_cross,np.convolve(np.diff(cosines),[0.5,0.5])) / 2, (cosines.size,1) ) 

        # Define linear interpolator, given the patchwork
        # N.B: if the later inputted energy exceeds the upper limit, the interpolator returns zero values, across all angles, yet no warnings 
        interp = RegularGridInterpolator((energies,cosines), diff_cross.transpose(), bounds_error=False, fill_value=0)

        # Interpolate over requested values; zero values will give no contribution to the spectra
        E /= 1000 # MeV
        sigma = interp((E,-1*costheta)) * 1e-31 # m**2/sr; the -1 factor is due to the carbon-12 being followed, and not the neutron (for which we have the emission probability, technically)
        
    elif all(p in [a,b] for p in ['d','t']): 
        # t(d,n)4he reaction. Tabulated data from Drosg2015, IAEA Nuclear Data Section, Vienna.
        df = dt_coeff_endf

        # Define energy,cosCM interpolating grid
        low_bound = 0.9 * np.min(E)*(constants.md+constants.mt)/constants.mt
        max_bound = 1.1 * np.max(E)*(constants.md+constants.mt)/constants.mt
        x, y = np.linspace(low_bound, max_bound,1000), np.linspace(-1,1,300)

        sigma_totale = sigma_tot(constants.mt/(constants.md+constants.mt)*x,'d-t')
        A = interp1d(df[0].to_numpy(dtype=np.float64),
                     df.iloc[:, 2:].to_numpy(dtype=np.float64), 
                     axis=0,fill_value=0, bounds_error=False)(x*10**-3)
        obj_f = (sigma_totale/2/A[:,0])[:,np.newaxis] * np.sum([A[:,l][:,np.newaxis]*
                                                                legendre(l)(y)[np.newaxis,:]
                                                                for l in range(df.shape[1]-2)],axis=0)
        interpolatore = RegularGridInterpolator((x,y), 
                                                np.nan_to_num(obj_f),
                                                fill_value=0, bounds_error=False)

        # Interpolate over requested values; zero values will give no contribution to the spectra
        E_deuterium_beam = E * (constants.md+constants.mt)/constants.mt
        sigma = interpolatore((E_deuterium_beam,costheta)) /2/np.pi # m**2/sr;

    else:
        # Assume that the other cross sections are isotropic
        sigma = sigma_tot(E, reaction=reaction)/(4*np.pi)

    return sigma


def reactivity(T, reaction='d-d'):
    """ Evaluate the parameterized thermonuclear reactivity (m**3/s) for temperature T (keV). """

    reaction = reaction.lower()
    a,b = reaction.split('-')

    if reaction == 'd-d':
        # Constants for the DD reaction
       	B_G  = 31.3970
        mrc2 = 937814

       	C1 =  5.43360e-12
       	C2 =  5.85778e-3
       	C3 =  7.68222e-3
       	C4 =  0.0
       	C5 = -2.96400e-6
        C6 =  0.0
        C7 =  0.0

    elif all(p in [a,b] for p in ['d','t']):
        # Constants for the DT reaction
       	B_G  = 34.3827
        mrc2 = 1124656

       	C1 =  1.17302e-9
       	C2 =  1.51361e-2
       	C3 =  7.51886e-2
       	C4 =  4.60643e-3
        C5 =  1.35000e-2
        C6 = -1.06750e-4
        C7 =  1.36600e-5


    # Evaluate the expression for the reactivity, Eq. (12) in the paper by Bosch and Hale.
    theta  = T / (1 - T*(C2+T*(C4+T*C6))/(1+T*(C3+T*(C5+T*C7))) )
    xi     = (B_G**2/4/theta)**(1.0/3.0)
    sigmav = C1*theta*np.sqrt(xi/mrc2/T**3) * np.exp(-3*xi)

    sigmav = sigmav*1e-6      # convert from cm**3/s to m**3/s

    return sigmav
