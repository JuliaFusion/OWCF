################################### load_TRANSP_interp_object.jl ##########################################

#### Description:
# This script is comprised of functions that will load from TRANSP and return interpolation objects.
# These interpolation objects will be functions of the normalized flux, rho, i.e. 
#
# rho = sqrt((psi - psi_axis)/(psi_bdry-psi_axis))
#
# where psi is the poloidal flux, psi_axis is the poloidal flux on-axis and psi_bdry is the poloidal flux 
# at the separatrix.
#
# Prior to loading this script (for example via 'include(misc/load_TRANSP_interp_object.jl)' when standing
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
# -

#### Outputs
# -

#### Saved files
# -

### Other
# -

# Script written by Henrik Järleblad. Last maintained 2024-08-13.
###################################################################################################

using NetCDF
using Interpolations

"""
    getTRANSPInterpObject(TRANSP_identifier, timepoint, filepath_TRANSP)
    getTRANSPInterpObject(-||-; verbose=false, density_cm=false)

- TRANSP_identifier: One of the variable identifiers used by TRANSP. E.g. "ND, "TI" or "OMEGA". Has to be a function of TRANSP'S "X" and "TIME". Provided as a String type.
- timepoint: The timepoint of interest in seconds, provided as a Float64 type
- filepath_TRANSP: The file path to a TRANSP .cdf output file, provided as a String type

Return an interpolation (actually extrapolation) object for the quantity identified by TRANSP_identifier.
This could be e.g. density, temperature etc. The interpolation object will be a function of normalized poloidal flux (going from 0 to 1).
- If verbose, the function will talk a lot. 
- If density_cm, the function will assume that the user tries to access bulk plasma density data from TRANSP. This is usually given in cm^-3. Since 
the OWCF works with m^-3, make a conversion, prior to creating and returning the interpolation object.
- If temperature_eV, the function will assume that the user tries to access bulk plasma temperature data from TRANSP. This is usually given in eV.
Since the OWCF works with keV, make a conversion, prior to creating and returning the interpolation object.
"""
function getTRANSPInterpObject(TRANSP_identifier::String, timepoint::Float64, filepath_TRANSP::String; verbose::Bool=false, density_cm::Bool=false, temperature_eV::Bool=false)
    verbose && println("getTRANSPInterpObject(): Loading time... ")
    time_array_TRANSP = ncread(filepath_TRANSP,"TIME")
    idx_min = argmin(abs.(time_array_TRANSP .- timepoint)) # Closest TRANSP timepoint to timepoint specified in start file
    verbose && println("getTRANSPInterpObject(): Time index of interest: $(idx_min).     Timepoint: $(time_array_TRANSP[idx_min]) s")
    verbose && println("getTRANSPInterpObject(): Loading poloidal flux... ")
    PLFLX_array = ncread(filepath_TRANSP,"PLFLX")[:,idx_min] # Load poloidal flux data

    q_TRANSP = ncread(filepath_TRANSP,TRANSP_identifier)[:,idx_min] # LOAD CORRESPONDING TRANSP DATA
    
    verbose && println("getTRANSPInterpObject(): Poloidal flux on-axis: $(PLFLX_array[1])")
    verbose && println("getTRANSPInterpObject(): Poloidal flux at separatrix: $(PLFLX_array[end])")
    PLFLX_axis = PLFLX_array[1]
    PLFLX_bdry = PLFLX_array[end]
    rho_array = sqrt.((PLFLX_array .- PLFLX_axis) ./ (PLFLX_bdry - PLFLX_axis)) # Compute the normalized poloidal flux points, rho

    verbose && println("getTRANSPInterpObject(): Creating interpolation (actually extrapolation) object... ")
    if density_cm
        verbose && println("getTRANSPInterpObject(): Density assumed. Converting from cm to m... ")
        q_TRANSP = (1.0e6) .*q_TRANSP
    end
    if temperature_eV
        verbose && println("getTRANSPInterpObject(): Temperature assumed. Converting from eV to keV... ")
        q_TRANSP = (1.0e-3) .*q_TRANSP
    end
    q_TRANSP_itp = Interpolations.interpolate((rho_array,), q_TRANSP, Gridded(Linear())) # Interpolation object
    return Interpolations.extrapolate(q_TRANSP_itp,0) # Add what happens in extrapolation scenario. 0 means we assume vacuum scrape-off layer (SOL)
end