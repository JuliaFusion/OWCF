################################ start_createCustomMagneticEquilibrium_template.jl #########################################
# This file contains all the inputs that the script OWCF/extra/createCustomMagneticEquilibrium.jl needs to create a custom
# toroidally symmetric magnetic equilibrium, modelled as a Solov'ev equilibrium. The model is based on the content in the 
# publication “One size fits all” analytic solutions to the Grad–Shafranov equation in the journal Phys. Plasmas 17, 032502 
# (2010), https://doi.org/10.1063/1.3328818. It is also based on the equations in the subchapter "Ideal MHD" by Freidberg, 
# Chapter 6.6.1. After the inputs are defined, this file executes the script OWCF/extra/createCustomMagneticEquilibrium.jl.
#
# OWCF/extra/createCustomMagneticEquilibrium.jl allows the user of the OWCF to create a custom magnetic equilibrium to be 
# used in the OWCF computations. The magnetic equilibrium will be saved as a .jld2 file that can be used by other OWCF scripts, 
# often as the 'filepath_equil' input variable.

#### The inputs are as follows:
# folderpath_OWCF - The path to where the OWCF folder is saved on your computed - String

# B0 - The magnetic field strength at the magnetic axis. In teslas - Float64
# R0 - The major radius position of the magnetic axis. In meters - Float64
# ϵ - The inverse aspect ratio a/R0 (a=minor radius) - Float64
# δ - The triangularity of the plasma - Float64
# κ - The elongation of the plasma - Float64
# α - The constant relating beta regime. However, if β is specified (see below), then the value of α does not matter. 
#     The algorithm will instead use the value of β to create the magnetic equilibrium - Float64
# β - The constant determining the beta regime. It is the ratio of the plasma pressure to the magnetic pressure. 
#     That is,
#               β = 2 * μ_0 * <p> / B^2
#     where μ_0 is the vacuum permeability, <p> is the mean plasma pressure and B is the mean total magnetic field strength.
#     For standard tokamaks, β usually has a value of around 0.04. The value of β can be left unspecified (:UNSPECIFIED). 
#     If so, the algorithm will instead use the value of α (see above) to create the magnetic equilibrium - Float64
# qstar - The kink safety factor (1.57 in ITER for example). It is defined as 
#               qstar = ϵB0/avg(B_pol)
#         where ϵ is the inverse aspect ratio, B0 is the magnetic field strength at the magnetic axis and avg(B_pol) is the 
#         mean poloidal magnetic field strength on the plasma surface. By default, qstar is pre-specified as ϵ*10, since in
#         many tokamaks, an approximation for the poloidal magnetic field strength is 1/10 of the toroidal magnetic field 
#         strength - Float64
# filepath_wall - The filepath to a .jld2 file containing (optional) data for the tokamak wall. The .jld2 file needs to have 
#                 the keys "R" (major radius coordinates in meters) and "z" (vertical coordinates in meters). The lengths of 
#                 R and z need to be equal. For every R_i element in R, the corresponding z value is assumed to be z_i in z.
#                 HOWEVER, some standard walls are available. If you specify filepath_wall to be either "JET" or "ITER",
#                 the (R,z) coordinates for those tokamaks will be used - String
# filename_out - The name of the magnetic equilibrium output data .jld2 file. By default, it is unspecified ("") and the default 
#                filename format is used to name the output data .jld2 file 
#                (see OWCF/extra/createCustomMagneticEquilibrium.jl for more info).
#                PLEASE NOTE! Do NOT include the .jld2 filename extension in the 'filename_out' input variable - String
# folderpath_out - The folderpath to the output folder where the results will be saved - String
# plot_equilibrium - If set to true, the magnetic equilibrium will be plotted and the plot will be saved in .png file format - Bool
# verbose - If true, the script will talk a lot! - Bool

### Optional input arguments are:
# B0_dir - The direction of the magnetic field. '1' means counter-clockwise (viewed from above). '-1' means the other way.
# Ip_dir - The direction of the plasma current. '1' means counter-clockwise (viewed from above). '-1' means the other way.
# diverted - If true, then the Solov'ev equilibrium will be diverted with an xpoint at (R0*(1-1.1*δ*ϵ),-R0*1.1*κ*ϵ) by default. This can also be specified manually by modifying the code below.
# debug - If true, the OWCF/extra/createCustomMagneticEquilibrium.jl script will run in debug mode

#### Other

# Script written by Henrik Järleblad. Last maintained 2026-09-24.
############################################################################################################################

## First you have to set the system specifications
using Distributed # Needed to be loaded, even though multi-core computations are not needed for createCustomLOS.jl.
folderpath_OWCF = "" # OWCF folder path. Finish with '/'

## -----------------------------------------------------------------------------
# Change directory to OWCF-folder on all external processes. Activate the environment there to ensure correct package versions
# as specified in the Project.toml and Manuscript.toml files.
@everywhere begin
    folderpath_OWCF = $folderpath_OWCF
    cd(folderpath_OWCF)
    using Pkg
    Pkg.activate(".")
end

## -----------------------------------------------------------------------------
@everywhere begin
    B0 = 0.0 # Magnetic field on-axis. Tesla
    R0 = 0.0 # Major radius position of magnetic axis. Meters
    ϵ = 0.0/R0 # Inverse aspect ratio, a/R0 where a is the radial (R) distance from the magnetic axis to the separatrix
    δ = 0.0 # Triangularity
    κ = 0.0 # Plasma elongation
    α = 0.0 # Constant relating beta regime. -0.155 works for ITER and results in a β-value of approx. 0.05
    β = :UNSPECIFIED
    qstar = ϵ*10 # Kink safety factor 
    filepath_wall = "JET" # Specify as "JET", "ITER", "/path/to/a/wall/file.jld2" or leave as ""
    filename_out = ""
    folderpath_out = ""
    plot_equilibrium = false # If set to true, the magnetic equilibrium will be plotted (poloidal cross section) after creation. For validation purposes.
    verbose = true

    # Optional inputs
    B0_dir = -1 # By default, assume the toroidal magnetic field points in the clockwise direction around the torus (viewed from above)
    Ip_dir = -1 # By default, assume that the plasma current runs clockwise around the torus (viewed from above)
    diverted = true # By default, assume that the user would like a divertor
    debug = false # Set to true for extra β(α) plot, if β is specified
end

## -----------------------------------------------------------------------------
# Then you execute the script
include(folderpath_OWCF*"extra/createCustomMagneticEquilibrium.jl")