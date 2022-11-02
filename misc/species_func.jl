################################### species_func.jl ##########################################

#### Description:
# This script stores functions that can be used to return the mass/charge of an ion (or electron)
# given a letter input. The letter input codes are the following:
# "D" -> Deuteron (deuterium)
# "T" -> Triton (tritium)
# "3he" -> Helium-3on (helium-3 nucleus)
# "p" -> Proton
# "e" -> Electron
# "4he" or "alpha" -> Helium-4on (alpha particle)
#
# getSpeciesAmu() can be used to return the particle mass in atomic mass units (2.0141017778 for deuterium for example)
# getSpeciesMass() can be used to return the particle mass in kg (3.344496968928368e-27 kg for deuterium for example)
# getSpeciesEcu() can be used to return the particle charge in elementary charge units (1.0 for deuterium for example)
# getSpeciesCharge() can be used to return the particle charge in Coulomb (1.60217733e-19 C for deuterium for example)
# getGCP() can be used to return the guiding-centre particle object for that ion.
#
# The guiding-centre particle objects are defined in GuidingCenterOrbits.jl.
#
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
#
# In future versions of the OWCF, as the capabilities grows, it is envisioned how the number of
# letter inputs for particles grows and possible also the number of pertaining utility functions.
# They will then be written and stored in this script.

# Script written by Henrik Järleblad. Last maintained 2022-08-27.
###################################################################################################

using GuidingCenterOrbits

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
    if lowercase(species_identifier)=="d"
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
    else
        error("getSpeciesMass got unknown species as input. Please re-try.")
    end
end
