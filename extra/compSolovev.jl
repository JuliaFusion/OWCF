################################ compSolovev.jl ######################################################
# This script will compute a Solov'ev magnetic equilibrium given the inputs. If no wall is specified
# (nothing), the script will assume the wall to be given by the Equilibrium.boundary() function (please 
# see the solovev.jl file of the Equilibrium.jl package, likely located at C:/Users/[your name]/.julia/
# packages/Equilibrium/) or type 'using Pkg; Pkg.add("https://github.com/JuliaFusion/Equilibrium.jl"); 
# using Equilibrium; ?Equilibrium.boundary()). Basically, the wall will then be right outside the plasma
# given by the Solov'ev equilibrium. If a .jld2 file is specified, the script will just try to load the 
# 'R' and 'z' data and use that as wall. It is assumed that 'R' and 'z' are vectors of Float64, where each 
# major radius coordinate in 'R' has a corresponding vertical coordinate in 'z'. The units are assumed to 
# be meters.

### The inputs are as follows:
# folderpath_OWCF - The folderpath to the OWCF folder - String
#
# B0 - The magnetic field strength at the magnetic axis - Tesla
# R0 - The major radius position of the magnetic axis - Meters
# ϵ - The inverse aspect ratio a/R0 (a=minor radius)
# δ - The triangularity of the plasma
# κ - The elongation of the plasma
# α - The constant relating beta regime (This must be chosen manually to achieve desired total plasma β-value)
# qstar - The kink safety factor (1.57 in ITER for example)
# filepath_wall - The filepath to a .jld2 file containing (optional) data for the tokamak wall - String
# folderpath_o - The folderpath to the output folder where the results will be saved - String
# verbose - If true, the script will talk a lot! - Bool

### Optional input arguments are:
# B0_dir - The direction of the magnetic field. '1' means counter-clockwise (viewed from above). '-1' means the other way.
# Ip_dir - The direction of the plasma current. '1' means counter-clockwise (viewed from above). '-1' means the other way.
# diverted - If true, then the Solov'ev equilibrium will be diverted with an xpoint at (R0*(1-1.1*δ*ϵ),-R0*1.1*κ*ϵ) by default. This can also be specified manually by modifying the code below.

### The output file will be named
# solovev_equilibrium_[DATE].jld2
# And it will have the fields
#   S - The Solov'ev Equilbrium.jl object
#   wall - The wall data

### Other
# To find a suitable total plasma β-value, you can try creating a Solov'ev equilibrium with a trial α value.
# You can then check the resulting total β-value by doing
# 
# using JLD2
# using Equilibrium
# myfile = jldopen("solovev_equilibrium_[DATE].jld2",false,false,false,IOStream)
# S = myfile["S"]
# close(myfile)
# βt = S.beta_t # Toroidal plasma β
# βp = S.beta_p # Poloidal plasma β
# β = βt*βp / (βp+βt)
#
# You can then increase/decrease your α value, until you find the desired total plasma β-value.

# Script written by Henrik Järleblad. Last maintained 2023-01-04.
######################################################################################################

## --------------------------------------------------------------------------
# Please specify the OWCF folder and let the script change directory to the 
# OWCF folder when compSolovev.jl is executed. This is to be able to load the
# correct versions of the Julia packages as specified in the Project.toml and 
# Manifest.toml files.
folderpath_OWCF = "" # Finish with '/'
cd(folderpath_OWCF)
using Pkg
Pkg.activate(".")

## --------------------------------------------------------------------------
# Required inputs
B0 = 1.6 # Magnetic field on-axis. Tesla
R0 = 0.894 # Major radius position of magnetic axis. Meters
ϵ = 0.274/R0 # Inverse aspect ratio
δ = 0.49 # Triangularity
κ = 1.8 # Plasma elongation
α = -2.4 # Constant relating beta regime. This must be chosen freely. -0.155 works for ITER and results in a β-value of 0.05
qstar = 3.5 # Kink safety factor
filepath_wall = ""
folderpath_o = ""
verbose = true

# Optional inputs
B0_dir = -1
Ip_dir = -1
diverted = true

## --------------------------------------------------------------------------
# Loading Julia packages
verbose && println("Loading Julia packages... ")
using Equilibrium
using JLD2
using Dates

## --------------------------------------------------------------------------
# Creating Solov'ev equilibrium from inputs
verbose && println("Creating Solov'ev equilibrium... ")
xpoint = (diverted ? (R0*(1-1.1*δ*ϵ),-R0*1.1*κ*ϵ) : nothing)
S = Equilibrium.solovev(B0, R0, ϵ, δ, κ, α, qstar; B0_dir = B0_dir, Ip_dir = Ip_dir, diverted = diverted, xpoint = xpoint, symmetric = (xpoint === nothing))

## --------------------------------------------------------------------------
# Checking wall file for data. If not available, create default wall for Solov'ev equilibrium
if !(isfile(filepath_wall))
    verbose && println("No wall data found in 'filepath_wall' (or could not be loaded). Creating default tokamak wall... ")
    wall = Equilibrium.boundary(S)
else
    verbose && println("Loadable file detected at "*filepath_wall*"... ")
    verbose && print("Attempting to load wall data... ")
    myfile = jldopen(filepath_wall,false,false,false,IOStream)
    R_array = myfile["R"]
    z_array = myfile["z"]
    close(myfile)
    verbose && println("Success!")
    verbose && println("Creating wall from data... ")
    wall = Equilibrium.boundary(R_array, z_array)
end

## --------------------------------------------------------------------------
# Save the data
verbose && println("Saving results... ")
date_and_time = split("$(Dates.now())","T")[1]
global filepath_output_orig = folderpath_o*"solovev_equilibrium_"*date_and_time
global filepath_output = deepcopy(filepath_output_orig)
global count = 1
while isfile(filepath_output*".jld2") # To take care of not overwriting files. Add _(1), _(2) etc
    global filepath_output = filepath_output_orig*"_($(Int64(count)))"
    global count += 1 # global scope, to surpress warnings
end
global filepath_output = filepath_output*".jld2"
myfile = jldopen(filepath_output,true,true,false,IOStream)
write(myfile,"S",S)
write(myfile,"wall",wall)
close(myfile)