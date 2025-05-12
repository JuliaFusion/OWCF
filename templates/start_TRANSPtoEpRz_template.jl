################################ start_TRANSPtoEpRz_template.jl #########################################
# This file contains all the inputs that the script OWCF/extra/TRANSP_fastion_data_to_EpRz.jl needs to 
# load fast-ion distribution data from a TRANSP-NUBEAM output file in .cdf file format, and interpolate 
# it onto a regular (E,p,R,z) grid, where E is the fast-ion energy in keV, p is the pitch (v_||/v), R 
# is the major radius in meters and z is the vertical coordinate in meters. This file also executes the 
# script OWCF/extra/TRANSP_fastion_data_to_EpRz.jl after the inputs are defined.
#
# In TRANSP-NUBEAM output files, the fast-ion distribution data is given on an irregular spiral grid with 
# the magnetic flux surfaces as the abscissas (please see the OWCF/extra/TRANSP_spiral_grid_example.png).
# To be able to use the fast-ion distribution in the OWCF, we need to interpolate it onto a regular 
# (E,p,R,z) grid. This is done via Delaunay tesselation interpolation.
#
# The total number of fast ions is conserved. The TRANSP fast-ion distribution data is given in units of 
# keV^-1 m^-3. Therefore, the resulting f(E,p,R,z) fast-ion distribution is saved in the same units.
# To get the total number of fast ions in the output data from the TRANSP_fastion_data_to_EpRz.jl script,
# integrate over all of (E,p,R,z) space together with the cylindrical coordinate Jacobian, i.e. 
#
# N_FI = ∫ f(E,p,R,z)*2*pi*R*dE*dp*dR*dz
#
# where N_FI is the total number of fast ions.

#### The inputs are as follows:
# btipsign - Should be set to 1 if sign(dot(J,B))>0. Else, set to -1.
# e_range - The energy (E) grid boundaries, specified as a Tuple (E_min, E_max) in keV. If left as nothing, 
#           the energy grid extrema already in the TRANSP-NUBEAM data will be used - Tuple
# filepath_distr - The path to the TRANSP fast-ion distribution in .cdf file format that you want to use - String
# folderpath_OWCF - The path to the OWCF-folder. Remember to finish with '/' since it's a folder and not a file. - String
# folderpath_out - The path to the folder in which you want your output. Remember to finish with '/'. - String
# nR - The number of major radius position grid points that you would like your (E,p,R,z) distribution to have - Int64
# nz - The number of vertical position grid points that you would like your (E,p,R,z) distribution to have - Int64
# p_range - The pitch (p=v_||/v) boundaries, specified as a Tuple (p_min,p_max). If left as nothing,
#           the pitch grid extrema already in the TRANSP-NUBEAM data will be used - Tuple
# plotting - If set to true, the TRANSP-NUBEAM fast-ion data will be plotted as it is loaded - Bool
# Rmin - The minimum value (lower boundary) for your major radius position grid points. In centimeters. - Float64
# Rmax - The maximum value (lower boundary) for your major radius position grid points. In centimeters. - Float64
# save_plots - If set to true, the plots plotted by 'plotting' input variable will be saved upon script completion - Bool
# species - The TRANSP particle species index. Usually set to 1. But if several fast-ion species are present in the 
#           TRANSP-NUBEAM data, other integers >1 might need to be used - Int64
# verbose - If set to true, the script will talk a lot! - Bool
# vverbose - If set to true, the script will talk too much! - Bool
# zmin - The minimum value (lower boundary) for your vertical position grid points. In centimeters. - Float64
# zmax - The maximum value (lower boundary) for your vertical position grid points. In centimeters. - Float64

#### Other
# 

# Script written by Henrik Järleblad. Last maintained 2025-05-09.
######################################################################################################

## First you have to set the system specifications
using Distributed # Needed to be loaded, even though multi-core computations are not needed for createCustomLOS.jl.
folderpath_OWCF = "" # OWCF folder path. Finish with '/'

## Navigate to the OWCF folder and activate the OWCF environment
cd(folderpath_OWCF)
using Pkg
Pkg.activate(".")

## -----------------------------------------------------------------------------
@everywhere begin
    btipsign=1 # '1' means that sign(J ⋅ B) > 0. '-1' means that sign(J ⋅ B) < 0
    e_range = nothing # Specified as (a,b). keV. If all TRANSP energies are desired, leave as 'nothing'
    filepath_distr = "" # for example '/Users/anna/TRANSP/JET/99500/V05/99500V05_fi_1.cdf'
    folderpath_OWCF = $folderpath_OWCF # Copy the 'folderpath_OWCF' input variable to all CPU processes
    folderpath_out = ""
    nR = 0 # The number of major radius R grid points
    nz = 0 # The number of vertical coordinate z grid points
    p_range = nothing # Specified as (a,b). If all TRANSP pitch values are desired, leave as 'nothing'
    plotting=true # If true, the script will plot results along the way. For safety checks.
    Rmin = 0.0 # cm. 183.5 suitable for JET
    Rmax = 0.0 # cm. 389.1 suitable for JET
    save_plots = false
    species=1 # '1' means the first species of the TRANSP file (usually the main fast-ion species)
    verbose=true # If true, the script will talk a lot!
    vverbose=false # If true, the script will talk a lot more! WARNING! You don't want this... It's just too much.
    zmin = 0.0 # cm. -174.5 suitable for JET
    zmax = 0.0 # cm. 198.4 suitable for JET
end
## -----------------------------------------------------------------------------
# Change directory to OWCF-folder on all external processes. Activate the environment there to ensure correct package versions
# as specified in the Project.toml and Manuscript.toml files.
@everywhere begin
    using Pkg
    cd(folderpath_OWCF)
    Pkg.activate(".")
end

## -----------------------------------------------------------------------------
# Then you execute the script
include(folderpath_OWCF*"extra/TRANSP_fastion_data_to_EpRz.jl")