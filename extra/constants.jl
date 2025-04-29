############################################### constants.jl ######################################################

#### Description:
# This file contains (physical) constants that the OWCF uses in various scripts. To avoid unnecessary warning prints,
# the constants have not been prefixed with the 'const' declaration. NEVERTHELESS, ALL CONSTANTS IN THIS SCRIPTS 
# SHOULD BE LEFT UNCHANGED. TO BE CLEAR, DO NOT CHANGE ANYTHING IN THIS FILE. Thank you!

#### Further info:
# As of the current OWCF version, the constants included in this file are as follows:
#   - OWCF_ICRF_STREAMLINES_SUPPORTED_COORD_SYSTEMS
#       These are the coordinate systems that are currently supported by the OWCF for computing streamlines that show
#       the direction of phase-space flow during heating by electromagnetic waves in the ion cyclotron range of 
#       frequencies (ICRF).
#   - OWCF_AMU_TO_KEV
#       This constant is used to convert values in atomic mass units (amu) to values in kiloelectronvolt (keV)
#   - OWCF_AMU_TO_KG
#       This constant is used to convert values in amu to values in kilogram (kg)
#   - OWCF_ELECTRON_MASS_AMU
#       This constant is the value of the electron (rest) mass, in amu
#   - OWCF_NEUTRON_MASS_AMU
#       This constant is the value of the neutron (rest) mass, in amu
#   - OWCF_PROTON_MASS_AMU
#       This constant is the value of the proton (rest) mass, in amu
#   - OWCF_e0
#       This constant is the value of the elementary charge, in Coulomb
#   - OWCF_ϵ0
#       This constant is the value of the vacuum permittivity, in Coulomb per Volt per meter
#   - OWCF_kB
#       This constant is the value of the Boltzmann constant, in Joules per Kelvin

#### Other:
# All the OWCF constants have been named with the prefix 'OWCF_' to clearly avoid conflicts with constants in other
# packages and/or user-defined constants.

# Script written by Henrik Järleblad. Last maintained 2025-04-25.
###################################################################################################################

OWCF_ICRF_STREAMLINES_SUPPORTED_COORD_SYSTEMS = "(E,p), (vpara,vperp), (E,mu,Pphi), (E,mu,Pphi,sigma), (E,Lambda,Pphi_n) and (E,Lambda,Pphi_n,sigma)"

OWCF_AMU_TO_KEV = 931494.0954 # keV. To be multiplied with atomic mass units (amu) to convert to keV
OWCF_AMU_TO_KG  = 1.660539040e-27 # kg. To be multiplied with atomic mass units (amu) to convert to kg
OWCF_ELECTRON_MASS_AMU = 5.48579909070e-4 # atomic mass units (amu)
OWCF_NEUTRON_MASS_AMU = 1.00866491588 # amu
OWCF_PROTON_MASS_AMU = 1.007276466879 # amu

OWCF_e0 = 1.602176634e-19 # Elementary charge, Coulomb
OWCF_ϵ0 = 8.8541878128e-12 # Permittivity of free space, C/(V*m)
OWCF_kB = 1.380649e-23 # Boltzmann constant, J/K