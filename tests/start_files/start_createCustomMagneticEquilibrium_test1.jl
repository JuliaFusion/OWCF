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
#
# B0 - The magnetic field strength at the magnetic axis - Tesla
# R0 - The major radius position of the magnetic axis - Meters
# ϵ - The inverse aspect ratio a/R0 (a=minor radius)
# δ - The triangularity of the plasma
# κ - The elongation of the plasma
# α - The constant relating beta regime (This must be chosen manually to achieve desired total plasma β-value)
# qstar - The kink safety factor (1.57 in ITER for example)
# filepath_wall - The filepath to a .jld2 file containing (optional) data for the tokamak wall - String
# filename_out - ...
# folderpath_out - The folderpath to the output folder where the results will be saved - String
# plot_equilibrium - If set to true, the magnetic equilibrium will be plotted when created - Bool
# save_equil_plot - If set to true, the magnetic equilibrium plot will be saved in .png format when plotted - Bool
# verbose - If true, the script will talk a lot! - Bool

### Optional input arguments are:
# B0_dir - The direction of the magnetic field. '1' means counter-clockwise (viewed from above). '-1' means the other way.
# Ip_dir - The direction of the plasma current. '1' means counter-clockwise (viewed from above). '-1' means the other way.
# diverted - If true, then the Solov'ev equilibrium will be diverted with an xpoint at (R0*(1-1.1*δ*ϵ),-R0*1.1*κ*ϵ) by default. This can also be specified manually by modifying the code below.

#### Other
# s

# Script written by Henrik Järleblad. Last maintained 2025-04-14.
############################################################################################################################

## First you have to set the system specifications
using Distributed # Needed to be loaded, even though multi-core computations are not needed for createCustomLOS.jl.

## -----------------------------------------------------------------------------
@everywhere begin
    B0 = 3.0 # Magnetic field on-axis. Tesla
    R0 = 3.0 # Major radius position of magnetic axis. Meters
    ϵ = 1.0/R0 # Inverse aspect ratio, a/R0 where a is the radial (R) distance from the magnetic axis to the separatrix
    δ = 0.3 # Triangularity
    κ = 1.4 # Plasma elongation
    α = -0.155 # Constant relating beta regime. This must be chosen freely. -0.155 works for ITER and results in a β-value of 0.05
    qstar = 3.33 # Kink safety factor. qstar = ϵB0/avg(B_pol) where the average of the poloidal magnetic field B_pol is at the separatrix 
    filepath_wall = "" # Specify, or leave as ""
    folderpath_OWCF = $folderpath_OWCF
    filename_out = "createCustomMagneticEquilibrium_test1"
    folderpath_out = folderpath_OWCF*"tests/outputs/"

    plot_equilibrium = $plot_test_results # SET TO TRUE, VIA THE plot_test_results INPUT VARIABLE IN THE OWCF/tests/run_tests.jl SCRIPT
    plot_equilibrium && (save_equil_plot = true) # If set to true, the magnetic equilibrium plot will be saved in .png file format
    verbose = true

    # Optional inputs
    B0_dir = -1 # By default, assume the toroidal magnetic field points in the clockwise direction around the torus (viewed from above)
    Ip_dir = -1 # By default, assume that the plasma current runs clockwise around the torus (viewed from above)
    diverted = true # By default, assume that the user would like a divertor
end

## -----------------------------------------------------------------------------
# Then you execute the script
include(folderpath_OWCF*"extra/createCustomMagneticEquilibrium.jl")