################################### temp_n_dens.jl ##########################################

#### Description:
# This script is comprised of a few functions, getAnalyticalTemp(T_axis::Float64, rho_p::Float64), 
# getAnalyticalDens(n_axis::Float64, rho_p::Float64), getTempProfileFromTRANSP(timepoint, ...
# filepath_thermal_distr::String, species::String) and getDensProfileFromTRANSP(timepoint, 
# filepath_thermal_distr::String, species::String;). They are all related to the temperature 
# and density profiles of the specific plasma species.
#
# Given a specified bulk temperature on axis and a normalized flux coordinate, getAnalyticalTemp() 
# returns the bulk plasma temperature at the normalized flux coordinate. The analytical bulk 
# temperature profile is adapted from a 4th degree polynomial fit of the bulk temperature profile of 
# JET shot No 96100.
#
# Given a specified bulk density on axis and a normalized flux coordinate, getAnalyticalDens()
# returns a bulk plasma density at the normalized flux coordinate. The analytical bulk density
# profile is an inverse exponential function based on the bulk density profile of JET shot #96100.
# It is given by
#
# /                   2               \
# | 2 -  ---------------------------  | * n_axis
# \      2 - ρ_pol* exp(8(ρ_pol-1))   /
#
# which tends towards n_axis as ρ_pol -> 0.0 and towards 0.0 as ρ_pol -> 1.0.
#
# The getTempProfileFromTRANSP() and getDensProfileFromTRANSP() functions return an interpolation 
# object representing the temperature and density profiles, respectively, at the specified timepoint,
# for the specified TRANSP .cdf pulse file, for the specified species (electrons are also valid).

#### Outputs
# -

#### Saved files
# -

### Other
# -

# Script written by Henrik Järleblad. Last maintained 2024-06-25.
###################################################################################################

using NetCDF
using Interpolations
include("rewriteReacts.jl")
include("load_TRANSP_interp_object.jl")

"""
    getAnalyticalTemp(T_axis, ρ_pol)

Input T_axis should be in keV. If input ρ_pol >= 1.0, then 0.0 will be returned (vacuum SOL assumed).
SOL = scrape-off layer.
"""
function getAnalyticalTemp(T_axis::T where {T<:Real}, ρ_pol::R where {R<:Real})

    if ρ_pol >= 1.0
        return 0.0 # Consider temperature to be zero outside LCFS
    end

    if ρ_pol <= 0.0
        return T_axis # Consider it to be the magnetic axis
    end

    return T_axis * (1.0+0.27947736204726326*ρ_pol-4.489653472449832*(ρ_pol^2)+6.117271127479773*(ρ_pol^3)-2.9026931867548087*(ρ_pol^4))
end

"""
    getAnalyticalDens(n_axis, ρ_pol)

Input n_axis should be in m^-3. If input ρ_pol >= 1.0, then 0.0 will be retured (vacuum SOL assumed).
SOL = scrape-off layer.
"""
function getAnalyticalDens(n_axis::T where {T<:Real}, ρ_pol::R where {R<:Real})

    if ρ_pol >= 1.0
        return 0.0
    end

    if ρ_pol <= 0.0
        return n_axis # Consider it to be the magnetic axis
    end

    return n_axis*(2-2/(2-ρ_pol*exp(8*(ρ_pol-1))))
end

"""
    getTempProfileFromTRANSP(timepoint::Float64, filepath_thermal_distr::String, species::String)
    getTempProfileFromTRANSP(-||-; kwargs... )

Given a timepoint (seconds, including decimals), load the TRANSP .cdf pulse file at filepath_thermal_distr and 
extract the temperature data at that timepoint, for the specified species. Return an interpolation object that 
represents the temperature profile at that point.

--- Input:
timepoint - The timepoint of interest (seconds) - Float64
filepath_thermal_distr - The file path to the TRANSP .cdf pulse file - String
species - The particle species of interest ("D","T","e" etc) - String
verbose - If true, the function will talk a lot! - Bool
"""
function getTempProfileFromTRANSP(timepoint::Float64,filepath_thermal_distr::String, species::String; verbose::Bool=false)
    T_identifier = ""
    if lowercase(species) in OWCF_SPECIES
        if species=="e"
            T_identifier = "TE"
        else
            T_identifier = "TI"
        end
    else
        error("No valid OWCF particle species recognized")
    end
    verbose && println("getTempProfileFromTRANSP(): T_identifier="*T_identifier)
    verbose && println("getTempProfileFromTRANSP(): Creating interpolation object... ")
    return getTRANSPInterpObject(T_identifier, timepoint, filepath_thermal_distr; verbose=verbose, temperature_eV=true)
end

"""
    getTempProfileFromTRANSP(timepoint::String, filepath_thermal_distr::String,species::String)
    getTempProfileFromTRANSP(-||-; kwargs... )

Given a timepoint (seconds, including decimals), load the TRANSP .cdf pulse file at filepath_thermal_distr and 
extract the temperature data at that timepoint. Return an interpolation object that represents the thermal 
temperature profile at that point.

--- Input:
timepoint - The timepoint of interest (seconds). Format "XX,YYYY" - String
filepath_thermal_distr - The file path to the TRANSP .cdf pulse file - String
verbose - If true, the function will talk a lot! - Bool
"""
function getTempProfileFromTRANSP(timepoint::String, filepath_thermal_distr::String,species::String; kwargs...)
    return getTempProfileFromTRANSP(parse(Float64,replace(timepoint,","=>".")), filepath_thermal_distr,species::String; kwargs...)
end

"""
    getTempProfileFromTRANSP(filepath_FI_distr::String, filepath_thermal_distr::String, species::String)
    getTempProfileFromTRANSP(-||-; kwargs... )

Given a TRANSP-NUBEAM output file 'filepath_FI_distr' (.cdf file format), load the TRANSP output file 'filepath_thermal_distr'
(.cdf file format) and extract the temperature data at the timepoint given by the 'TIME' key value loaded from the 'filepath_FI_distr'
file, for a specific particle species 'species'. Return an interpolation object that represents the temperature profile as a function of rho_pol 0 to 1.

--- Input:
filepath_FI_distr - The file path to the TRANSP-NUBEAM .cdf output file - String
filepath_thermal_distr - The file path to the TRANSP .cdf pulse file - String
species - The species of interest, e.g. "D", "T", "e" - String
"""
function getTempProfileFromTRANSP(filepath_FI_distr::String, filepath_thermal_distr::String, species::String; verbose=false, kwargs... )
    timepoint = Float64(ncread(filepath_FI_distr,"TIME")[1])
    verbose && println("getTempProfileFromTRANSP(): Timepoint loaded from $(filepath_FI_distr) input argument: $(timepoint) seconds... ")
    return getTempProfileFromTRANSP(timepoint,filepath_thermal_distr,species; verbose=verbose)
end

"""
    getDensProfileFromTRANSP(timepoint::Float64, filepath_thermal_distr::String, species::String)
    getDensProfileFromTRANSP(-||-; kwargs... )

Given a timepoint (seconds, including decimals), load the TRANSP .cdf pulse file at filepath_thermal_distr and 
extract the density data at that timepoint, for a specific species 'species'. Return an interpolation object 
that represents the density profile at that point.

--- Input:
timepoint - The timepoint of interest (seconds) - Float64
filepath_thermal_distr - The file path to the TRANSP .cdf pulse file - String
species - The species of interest, e.g. "D", "T", "e" - String
verbose - If true, the function will talk a lot! - Bool
"""
function getDensProfileFromTRANSP(timepoint::Float64,filepath_thermal_distr::String, species::String; verbose::Bool=false)
    verbose && println("getDensProfileFromTRANSP(): OWCF species label: "*species)
    species_TRANSP = thermalSpecies_OWCFtoTRANSP(species) # Convert from OWCF convention (e.g. 3he for Helium-3) to TRANSP convention (e.g. HE3 for Helium-3)
    verbose && println("getDensProfileFromTRANSP(): TRANSP species label: "*species_TRANSP)
    
    verbose && println("getDensProfileFromTRANSP(): Creating interpolation object... ")
    return getTRANSPInterpObject("N"*species_TRANSP, timepoint, filepath_thermal_distr; verbose=verbose, density_cm=true)
end

"""
    getDensProfileFromTRANSP(timepoint::String, filepath_thermal_distr::String, species::String)
    getDensProfileFromTRANSP(-||-; kwargs... )

Given a timepoint (seconds, including decimals), load the TRANSP .cdf output file 'filepath_thermal_distr' and 
extract the density data at that timepoint, for a specific species 'species'. Return an interpolation object 
that represents the density profile at that point.

--- Input:
timepoint - The timepoint of interest (seconds). Format "XX,YYYY" - String
filepath_thermal_distr - The file path to the TRANSP .cdf pulse file - String
species - The species of interest, e.g. "D", "T", "e" - String
"""
function getDensProfileFromTRANSP(timepoint::String, filepath_thermal_distr::String, species::String; kwargs...)
    return getDensProfileFromTRANSP(parse(Float64,replace(timepoint,","=>".")), filepath_thermal_distr, species; kwargs...)
end

"""
    getDensProfileFromTRANSP(filepath_FI_distr::String, filepath_thermal_distr::String, species::String)
    getDensProfileFromTRANSP(-||-; kwargs... )

Given a TRANSP-NUBEAM output file 'filepath_FI_distr' (.cdf file format), load the TRANSP output file 'filepath_thermal_distr'
(.cdf file format) and extract the density data at the timepoint given by the 'TIME' key value loaded from the 'filepath_FI_distr'
file, for a specific particle species 'species'. Return an interpolation object that represents the density profile as a function of rho_pol 0 to 1.

--- Input:
filepath_FI_distr - The file path to the TRANSP-NUBEAM .cdf output file - String
filepath_thermal_distr - The file path to the TRANSP .cdf pulse file - String
species - The species of interest, e.g. "D", "T", "e" - String
"""
function getDensProfileFromTRANSP(filepath_FI_distr::String, filepath_thermal_distr::String, species::String; verbose=false, kwargs... )
    timepoint = Float64(ncread(filepath_FI_distr,"TIME")[1])
    verbose && println("getDensProfileFromTRANSP(): Timepoint loaded from $(filepath_FI_distr) input argument: $(timepoint) seconds... ")
    return getDensProfileFromTRANSP(timepoint,filepath_thermal_distr,species; verbose=verbose)
end