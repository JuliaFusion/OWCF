################################### temp_n_dens.jl ##########################################

#### Description:
# This script is comprised of two functions, getAnalyticalTemp(T_axis::Float64, rho_p::Float64) 
# and getAnalyticalDens(n_axis::Float64, rho_p::Float64). 
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
# T_axis - Float64
# ρ_pol - Float64
#
# n_axis - Float64
# ρ_pol - Float64

#### Outputs
# T_at_ρ_pol - Float64
#
# n_at_ρ_pol - Float64

#### Saved files
# -

### Other
# -

# Script written by Henrik Järleblad. Last maintained 2022-09-10.
###################################################################################################

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