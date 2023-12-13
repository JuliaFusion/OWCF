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
#
# Prior to loading this collection of functions (for example via 'include(misc/temp_n_dens.jl)' when standing
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

#### Inputs:
# PLEASE NOTE! In your script, you have to define the folderpath_OWCF variable, which is a string 
# with the file path to the OWCF folder on your computer, before loading this file with 
# include(folderpath_OWCF*"/misc/temp_n_dens.jl")

#### Outputs
# -

#### Saved files
# -

### Other
# -

# Script written by Henrik Järleblad. Last maintained 2023-11-23.
###################################################################################################

using NetCDF
using Interpolations
include(folderpath_OWCF*"/misc/rewriteReacts.jl")
include(folderpath_OWCF*"/misc/species_func.jl")

"""
    getAnalyticalTemp(T_axis, ρ_pol)

Input T_axis should be in keV. If input ρ_pol >= 1.0, then 0.0 will be returned (vacuum SOL assumed).
SOL = scrape-off layer.
"""
function getAnalyticalTemp(T_axis::Float64, ρ_pol::Float64)

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
function getAnalyticalDens(n_axis::Float64, ρ_pol::Float64)

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
    verbose && println("getTempProfileFromTRANSP(): Loading time... ")
    time_array_TRANSP = ncread(filepath_thermal_distr,"TIME")
    idx_min = argmin(abs.(time_array_TRANSP .- timepoint)) # Closest TRANSP timepoint to timepoint specified in start file
    verbose && println("getTempProfileFromTRANSP(): Time index of interest: $(idx_min)")
    verbose && println("getTempProfileFromTRANSP(): Loading poloidal flux... ")
    PLFLX_array = ncread(filepath_thermal_distr,"PLFLX")[:,idx_min] # Load poloidal flux data
    
    T_identifier = ""
    if lowercase(species) in VALID_OWCF_SPECIES
        if species=="e"
            global T_identifier = "TE"
        else
            global T_identifier = "TI"
        end
    else
        error("No valid OWCF species recognized")
    end
    verbose && println("getTempProfileFromTRANSP(): T_identifier="*T_identifier)
    temp_TRANSP = ncread(filepath_thermal_distr,T_identifier)[:,idx_min] # LOAD CORRESPONDING TRANSP TEMPERATURE DATA

    verbose && println("getTempProfileFromTRANSP(): Poloidal flux on-axis: $(PLFLX_array[1])")
    verbose && println("getTempProfileFromTRANSP(): Poloidal flux at separatrix: $(PLFLX_array[end])")
    PLFLX_axis = PLFLX_array[1]
    PLFLX_bdry = PLFLX_array[end]
    rho_array = sqrt.((PLFLX_array .- PLFLX_axis) ./ (PLFLX_bdry - PLFLX_axis)) # Compute the normalized poloidal flux points, rho
    
    verbose && println("getTempProfileFromTRANSP(): Creating interpolation object... ")
    temp_TRANSP_itp = Interpolations.interpolate((rho_array,),(1.0e-3) .*temp_TRANSP, Gridded(Linear())) # eV to keV
    return Interpolations.extrapolate(temp_TRANSP_itp,0) # Add what happens in extrapolation scenario. 0 means we assume vacuum scrape-off layer (SOL)
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
    verbose && println("getDensProfileFromTRANSP(): Loading time... ")
    time_array_TRANSP = ncread(filepath_thermal_distr,"TIME")
    idx_min = argmin(abs.(time_array_TRANSP .- timepoint)) # Closest TRANSP timepoint to timepoint specified in start file
    verbose && println("getDensProfileFromTRANSP(): Time index of interest: $(idx_min)")
    verbose && println("getDensProfileFromTRANSP(): Loading poloidal flux... ")
    PLFLX_array = ncread(filepath_thermal_distr,"PLFLX")[:,idx_min] # Load poloidal flux data
    
    verbose && println("getDensProfileFromTRANSP(): OWCF species label: "*species)
    species_TRANSP = thermalSpecies_OWCFtoTRANSP(species) # Convert from OWCF convention (e.g. 3he for Helium-3) to TRANSP convention (e.g. HE3 for Helium-3)
    verbose && println("getDensProfileFromTRANSP(): TRANSP species label: "*species_TRANSP)
    dens_TRANSP = ncread(filepath_thermal_distr,"N"*species_TRANSP)[:,idx_min] # LOAD CORRESPONDING TRANSP DENSITY DATA

    verbose && println("getDensProfileFromTRANSP(): Poloidal flux on-axis: $(PLFLX_array[1])")
    verbose && println("getDensProfileFromTRANSP(): Poloidal flux at separatrix: $(PLFLX_array[end])")
    PLFLX_axis = PLFLX_array[1]
    PLFLX_bdry = PLFLX_array[end]
    rho_array = sqrt.((PLFLX_array .- PLFLX_axis) ./ (PLFLX_bdry - PLFLX_axis)) # Compute the normalized poloidal flux points, rho
    
    verbose && println("getDensProfileFromTRANSP(): Creating interpolation object... ")
    dens_TRANSP_itp = Interpolations.interpolate((rho_array,),(1.0e6) .*dens_TRANSP, Gridded(Linear())) # cm^-3 to m^-3
    return Interpolations.extrapolate(dens_TRANSP_itp,0) # Add what happens in extrapolation scenario. 0 means we assume vacuum scrape-off layer (SOL)
end
"""
    getDensProfileFromTRANSP(timepoint::String, filepath_thermal_distr::String, species::String)
    getDensProfileFromTRANSP(-||-; kwargs... )

    Given a timepoint (seconds, including decimals), load the TRANSP .cdf pulse file at filepath_thermal_distr and 
    extract the density data at that timepoint, for a specific species 'species'. Return an interpolation object 
    that represents the density profile at that point.

--- Input:
timepoint - The timepoint of interest (seconds). Format "XX,YYYY" - String
filepath_thermal_distr - The file path to the TRANSP .cdf pulse file - String
species - The species of interest, e.g. "D", "T", "e" - String
verbose - If true, the function will talk a lot! - Bool
"""
function getDensProfileFromTRANSP(timepoint::String, filepath_thermal_distr::String, species::String; kwargs...)
    return getDensProfileFromTRANSP(parse(Float64,replace(timepoint,","=>".")), filepath_thermal_distr, species; kwargs...)
end