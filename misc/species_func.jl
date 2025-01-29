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

const OWCF_ChemElemSymbol_to_Z = Dict("e"=>-1,"h"=>1,"he"=>2,"li"=>3,"be"=>4,"b"=>5,"c"=>6,"n"=>7,"o"=>8,"f"=>9,"ne"=>10,"na"=>11,"mg"=>12,) # A dictionary with chemical element (and electron) symbols (lowercase) as keys and nuclear charge number (Z) as values
const OWCF_SPECIES = [key for key in keys(OWCF_ChemElemSymbol_to_Z)] # Only the keys of the OWCF_ChemElemSymbol_to_Z dictionary 
const OWCF_e0 = 1.6021766208e-19 # Coulomb
const OWCF_u_keV = 931494.0954 # keV
const OWCF_u_kg  = 1.660539040e-27 # kg
const OWCF_neutron_mass_amu = 1.00866491588 # atomic mass units (amu)
const OWCF_proton_mass_amu = 1.007276466879 # amu
##############################################
# The rows of the matrix below correspond to atomic mass number (A) and the columns correspond to nuclear charge number (Z), 
# in accordance with Wang Meng et al 2017 Chinese Phys. C 41 030003. The values of the matrix elements are written in micro atomic 
# mass units (μu), and the whole matrix is multiplied by 10^-6 (as you can see below) to have the elements be accessible in atomic 
# mass units (u).
const OWCF_ATOMIC_MASS_TABLE = (1.0e-6) .*[] 
##############################################

function getSpeciesAmu(species_identifier::AbstractString)
    return getSpeciesMass(species_identifier) / GuidingCenterOrbits.mass_u
end

function getSpeciesMass(species_identifier::AbstractString)
    return (getGCP(species_identifier)(0.0,0.0,0.0,0.0)).m # 0.0, 0.0, 0.0, 0.0 is just to activate the function
end

function getSpeciesEcu(species_identifier::AbstractString)
    return (getGCP(species_identifier)(0.0,0.0,0.0,0.0)).q # 0.0, 0.0, 0.0, 0.0 is just to activate the function
end

function getSpeciesCharge(species_identifier::AbstractString)
    return getSpeciesEcu(species_identifier) * GuidingCenterOrbits.e0
end

function getGCP(species_identifier::AbstractString)
    if species_identifier=="n"
        return (GuidingCenterOrbits.GCParticle())
    end
    if lowercase(species_identifier)=="d"
        q = 1 # Ecu
        m = ((OWCF_neutron_mass_amu+OWCF_proton_mass_amu)*OWCF_u_keV)/OWCF_u_keV

        return GuidingCenterOrbits.GCDeuteron
    elseif lowercase(species_identifier)=="t"
        return GuidingCenterOrbits.GCTriton
    elseif lowercase(species_identifier)=="3he"
        return GuidingCenterOrbits.GCHelium3
    elseif lowercase(species_identifier)=="p"
        return GuidingCenterOrbits.GCProton
    elseif lowercase(species_identifier)=="e"
        return GuidingCenterOrbits.GCElectron
    elseif (lowercase(species_identifier)=="alpha" || lowercase(species_identifier)=="4he")
        return GuidingCenterOrbits.GCAlpha
    elseif lowercase(species_identifier)==""
    else
        error("getGCP got unknown species as input. Please re-try.")
    end
end
