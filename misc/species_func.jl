################################### species_func.jl ##########################################

#### Description:
# This file stores functions that can be used to return the mass/charge of an ion (or electron)
# given a letter input. The most common letter input codes are the following:
#
# "H" or "proton" -> Proton
# "D" or "2H" -> Deuteron (deuterium)
# "T" or "3H" -> Triton (tritium)
# "3he" -> Helium-3on (helium-3 nucleus)
# "e" -> Electron
# "4he" or "alpha" -> Helium-4on (alpha particle)
# 
# The rest of the letter input codes can be found (without explanation) in the constant 
# 'OWCF_ChemElemSymbol_to_Z' below this file description. Also, PLEASE NOTE, do NOT use "p" to 
# refer to a proton. This could be confused with the lowercase symbol of the chemical element 
# Phosphorus. Instead, please use "H", "h" or "proton".

# If the functions detect an input equal to a lowercase "N", i.e. "n", they will assume you mean 
# a neutron, and return a guiding-centre particle (GCParticle) object with mass equal to the 
# mass of a neutron and zero charge.

# PLEASE NOTE! The ion mass (m_i) that will be returned for an ion with nuclear charge number (Z) 
# is the standard atomic mass (m_A) MINUS the mass of all electrons needed to make the ion neutral 
# (Z x m_e, where m_e is the electron rest mass) PLUS the total binding energy of all removed 
# electrons (B_e(Z)). Thus, the formula for m_i is:
#
#   m_i = m_A - Z*m_e + B_e(Z)
#
# where B_e(Z) can be approximated (in eV) by:
#
#   B_e(Z) = 14.4381*(Z^2.39) + 1.55468*(10^-6)*(Z^5.35)
#
# The above expressions and the values in OWCF_ATOMIC_MASS_TABLE are from the following paper:
#
#   Wang Meng et al 2017 Chinese Phys. C 41 030003

#### Further info:
# getSpeciesAmu() can be used to return the particle mass in atomic mass units (2.0141017778 for deuterium for example)
# getSpeciesMass() can be used to return the particle mass in kg (3.344496968928368e-27 kg for deuterium for example)
# getSpeciesEcu() can be used to return the particle charge in elementary charge units (1.0 for deuterium for example)
# getSpeciesCharge() can be used to return the particle charge in Coulomb (1.60217733e-19 C for deuterium for example)
# getGCP() can be used to return the guiding-centre particle object for that ion.
#
# The guiding-centre particle (GCParticle) structs are defined in the GuidingCenterOrbits.jl Julia package.

#### Other:
# Prior to loading this collection of functions (for example via 'include(misc/species_func.jl)' when standing
# in the OWCF-folder), you should have performed the usual activation of the OWCF Julia environment.
# This is done as follows
#
# folderpath_OWCF = "/path/to/the/OWCF/folder/"
# using Distributed
# @everywhere begin
#     using Pkg
#     cd(folderpath_OWCF)
#     Pkg.activate(".")
# end
# 
# This is performed in e.g. every OWCF start template file. It is also possible to skip the Distributed.jl
# package. This is then done as 
#
# folderpath_OWCF = "/path/to/the/OWCF/folder/"
# using Pkg
# cd(folderpath_OWCF)#
# Pkg.activate(".")

# Script written by Henrik Järleblad. Last maintained 2025-01-29.
###################################################################################################

using GuidingCenterOrbits
const OWCF_e0 = 1.6021766208e-19 # Coulomb
const OWCF_u_keV = 931494.0954 # keV
const OWCF_u_kg  = 1.660539040e-27 # kg
const OWCF_electron_mass_amu = 5.48579909070e-4 # atomic mass units (amu)
const OWCF_neutron_mass_amu = 1.00866491588 # amu
const OWCF_proton_mass_amu = 1.007276466879 # amu
# A dictionary with chemical element (and electron) symbols (lowercase) as keys and nuclear charge number (Z) as values
const OWCF_ChemElemSymbol_to_Z = Dict("e"=>-1,"h"=>1,"proton"=>1,"d"=>1,"t"=>1,"he"=>2,"alpha"=>2,"li"=>3,"be"=>4,"b"=>5,"c"=>6,"n"=>7,"o"=>8,"f"=>9,"ne"=>10,
                                      "na"=>11,"mg"=>12,"al"=>13,"si"=>14,"p"=>15,"s"=>16,"cl"=>17,"ar"=>18,"k"=>19,"ca"=>20,"sc"=>21,"ti"=>22,"v"=>23,"cr"=>24,
                                      "mn"=>25,"fe"=>26) # Useful fusion stops at iron (Fe)
# A dictionary with chemical element (and electron) symbols (lowercase) as keys and atomic mass number (A) as values
const OWCF_ChemElemSymbol_to_A = Dict("e"=>OWCF_electron_mass_amu,"h"=>1,"proton"=>1,"d"=>2,"t"=>3,"he"=>4,"alpha"=>4,"li"=>7,"be"=>9,"b"=>11,"c"=>12,"n"=>14,
                                      "o"=>16,"f"=>19,"ne"=>20,"na"=>23,"mg"=>24,"al"=>27,"si"=>28,"p"=>31,"s"=>32,"cl"=>35,"ar"=>40,"k"=>39,"ca"=>40,"sc"=>45,
                                      "ti"=>48,"v"=>51,"cr"=>52,"mn"=>55,"fe"=>56) # Useful fusion stops at iron (Fe)
const OWCF_SPECIES = [key for key in keys(OWCF_ChemElemSymbol_to_Z)] # Only the keys of the OWCF_ChemElemSymbol_to_Z dictionary 
##############################################
# The rows of the matrix below correspond to atomic mass number (A) and the columns correspond to nuclear charge number (Z), 
# in accordance with Wang Meng et al 2017 Chinese Phys. C 41 030003. The values of the matrix elements are written in micro atomic 
# mass units (μu), and the whole matrix is multiplied by 10^-6 (as you can see below) to have the elements be accessible in atomic 
# mass units (u).
const OWCF_ATOMIC_MASS_TABLE = (1.0e-6) .*[1007825.03224 NaN NaN NaN NaN NaN NaN NaN; 
                                           2014101.77811 NaN NaN NaN NaN NaN NaN NaN; 
                                           3016049.28199 3016029.32265 3030780 NaN NaN NaN NaN NaN; 
                                           4026430 4002603.25413 4027190 NaN NaN NaN NaN NaN; 
                                           5035310 5012057 5012540 5039870 NaN NaN NaN NaN; 
                                           6044960 6018885.89 6015122.8874 6019726 6050800 NaN NaN NaN;
                                           7052750 7027991 7016003.437 7016928.72 7029712 NaN NaN NaN;
                                           NaN 8033934.39 8022486.25 8005305.10 8024607.3 8037643 NaN NaN;
                                           NaN 9043950 9026790.19 9012183.07 9013329.6 9031037.2 NaN NaN;
                                           NaN 10052820 10035483 10013534.70 10012936.862 10016853.22 10041650 NaN;
                                           NaN NaN 11043723.6 11021661.08 11009305.167 11011432.60 11026090 NaN;
                                           NaN NaN 12052610 12026922.1 12014352.6 12000000.0 12018613.2 12034262] 
##############################################
const number_pattern = r"[0-9]+" # To be able to isolate atomic mass numbers in Strings
const letter_pattern = r"[a-z]+" # To be able to isolate chemical element symbols in Strings

# ---------------------------------------------------------------------------------------------------------------------------------
"""
    B_e(Z)

Function description to be written here.
"""
function OWCF_B_e(Z)
    return 14.4381*(Z^2.39) + 1.55468*(10^-6)*(Z^5.35)
end

"""
    getGCP(species_identifier)
    getGCP(-||-; E=0.0, p=0.0, R=0.0, z=0.0)

Return a guiding-centre particle struct, as defined in the Julia package 
GuidingCenterOrbits.jl in the file particles.jl. The particle species is 
determined by the species_identifier. Please see the keys of the OWCF_ChemElemSymbol_to_Z
dictionary above for a list of all available species. Example of valid inputs:
    - getGCP("H") # Hydrogen nucleus
    - getGCP("proton") # Proton (same as H)
    - getGCP("D") # Deuterium
    - getGCP("T") # Tritium
    - getGCP("3he") # Helium-3
    - getGCP("4he") # Helium-4
    - getGCP("alpha") # Alpha particle (same as 4he)
    - getGCP("e") # Electron
Example of non-valid (!!!) inputs:
    - getGCP("p") # Proton
DO NOT use "p" for a proton, as this may be confused with the lowercase chemical element 
symbol for Phosphorus. Use "H" or "proton" instead. In addition, if a lowercase "N" is 
provided as the 'species_identifier' input, i.e. "n", the function will assume that you mean 
a neutron and return a guiding-centre particle struct with mass equal to the neutron mass 
and zero charge. The keyword arguments are:
- E: The particle energy in keV - Float64
- p: The particle pitch (v_||/v) - Float64
- R: The particle major radius coordinate in meters - Float64
- z: The particle vertical coordinate in meters - Float64
- verbose: If true, the function will talk a lot!
"""
function getGCP(species_identifier::AbstractString; E=0.0, p=0.0, R=0.0, z=0.0, verbose=false)
    if species_identifier=="n"
        @warn "getGCP() input was $(species_identifier). Assuming a neutron."
        return (GuidingCenterOrbits.GCParticle(E,p,R,z,OWCF_neutron_mass_amu*OWCF_u_kg,0)) # Neutron mass, zero charge
    end
    if species_identifier=="g"
        @warn "getGCP() input was $(species_identifier). Assuming a gamma photon."
        return (GuidingCenterOrbits.GCParticle(E,p,R,z,0.0,0)) # Zero mass, zero charge
    end
    species_identifier_lowercase = lowercase(species_identifier)
    species_A = match(number_pattern,species_identifier_lowercase) # The atomic mass number, e.g. RegexMatch("10") of "3he". nothing if only e.g. "he"
    species_ChemElemSymbol = match(letter_pattern,species_identifier_lowercase) # The chemical element symbol
    if !(species_ChemElemSymbol in OWCF_SPECIES)
        error("getGCP() input $(species_ChemElemSymbol) did not match any of the particle species available in the OWCF. Please choose from $(OWCF_SPECIES). Please correct and re-try.")
    end
    if isnothing(species_A) # No atomic mass number included with the species_identifier
        species_A = OWCF_ChemElemSymbol_to_A[species_identifier_lowercase] # Determine it via look-up Dictionary
    else
        species_A = parse(Int64,species_A.match) # Otherwise, parse it
    end
    species_Z = OWCF_ChemElemSymbol_to_Z[species_ChemElemSymbol] # Determine the nuclear charge number
    if species_Z==-1 # If electron
        verbose && println("getGCP() input was $(species_identifier). Deduced electron!")
        return (GuidingCenterOrbits.GCParticle(E,p,R,z,OWCF_electron_mass_amu*OWCF_u_kg,species_Z))
    end
    if species_A>size(OWCF_ATOMIC_MASS_TABLE,1)
        m_A = species_A # Use the atomic mass number as an approximation. This will be refined in future versions of the OWCF
    else
        m_A = OWCF_ATOMIC_MASS_TABLE[species_A,species_Z] # The atomic mass (in amu)
    end
    species_B_e = (1.0e-3)*OWCF_B_e(species_Z) # eV to keV
    m_N = m_A - species_Z*OWCF_electron_mass_amu + species_B_e/OWCF_u_keV # Please see Wang Meng et al 2017 Chinese Phys. C 41 030003

    return GuidingCenterOrbits.GCParticle(E,p,R,z,m_N*OWCF_u_kg,species_Z)
end

"""
    getSpeciesMass(species_identifier)

Function description to be written here.
"""
function getSpeciesMass(species_identifier::AbstractString)
    return getGCP(species_identifier).m
end

"""
    getSpeciesAmu(species_identifier)

Function description to be written here.
"""
function getSpeciesAmu(species_identifier::AbstractString)
    return getSpeciesMass(species_identifier) / GuidingCenterOrbits.mass_u
end

"""
    getSpeciesEcu(species_identifier)

Function description to be written here.
"""
function getSpeciesEcu(species_identifier::AbstractString)
    return getGCP(species_identifier).q
end

"""
    getSpeciesCharge(species_identifier)

Function description to be written here.
"""
function getSpeciesCharge(species_identifier::AbstractString)
    return getSpeciesEcu(species_identifier) * GuidingCenterOrbits.e0
end