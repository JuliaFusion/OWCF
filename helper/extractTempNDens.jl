######################################### extractTempNDens.jl #########################################

#### Description:
# This script allows you to extract the temperature and density profiles of a certain tokamak shot,
# provided that you have the corresponding TRANSP .cdf shot file, the fast-ion TRANSP .cdf file for the shot
# and the .eqdsk/.geqdsk magnetic equilibrium file.
# The bulk ion temperature and bulk ion density will be saved as function of normalized flux coordinate,
# that is T(ρ_pol) and n(ρ_pol). Where
#
# ρ_pol = sqrt((ψ - ψ_axis)  / (ψ_sep - ψ_axis) )
#
# where ψ is the flux function, ψ_axis is the flux function value at the magnetic axis and
# ψ_sep is the flux function value at the separatrix. When specifying the ion SPECIES ("D", "T" etc), please
# remember that its data must exist in the .cdf fast-ion file that you provide.

#### Inputs
# folderpath_OWCF - The path to the OWCF on your computer - String
# filepath_TRANSP - The path to the TRANSP .cdf file of the shot. E.g. 96100J01.cdf - String
# filepath_TRANSP_FI - The path to the fast-ion .cdf file of the TRANSP shot file - String
# filepath_equil - The path to the .eqdsk or .jld2 file - String
# SPECIES - The ion SPECIES whose bulk temperature and bulk density you want - String
# folderpath_o - The path to the desired output folder
# timepoint - The timepoint of the TRANSP NUBEAM fast-ion distribution simulation. Format XX,YYYY where XX are seconds and YYYY are decimals - String
# verbose - If set to true, the script will be very talkative! - Bool

#### Outputs
# -

#### Saved files
# [TRANSP_id]_at[timepoint]s_[SPECIES]_TempNDens.jld2
# This saved file will have the fields:
#   rho_p_array - The normalized flux coordinate grid points - Array{Float64,1}
#   T_i - The bulk ion temperature values for the grid points in rho_p_array - Array{Float64,1}
#   n_i - The bulk ion density values for the grid points in rho_p_array - Array{Float64,1}

### Other
# Before you run the script, don't forget to specify the path to the OWCF folder via the
# folderpath_OWCF variable. This is to ensure Python is called correctly by Julia. For example
#
# folderpath_OWCF = "C:/Users/bobjk/codes/OWCF/"
#
# Change C:/Users/bobjk/codes/OWCF/ to the corresponding correct path on your own computer.
# Last, and maybe least, the rho_p_array will be computed in the equitorial plane (same cylindrical
# coordinate z level as the magnetic axis).

# Script written by Henrik Järleblad. Last maintained 2022-08-26.
#######################################################################################################

## --------------------------------------------------------------------------
# Please specify the OWCF folder and let the app change directory to the 
# OWCF folder when extractTempNDens.jl is executed. This is to be able to load the
# correct versions of the Julia packages as specified in the Project.toml and 
# Manifest.toml files.
folderpath_OWCF = "/Volumes/GoogleDrive/My Drive/DTU/codes/OWCF/" # Finish with '/'
cd(folderpath_OWCF)
using Pkg
Pkg.activate(".")

## --------------------------------------------------------------------------
# Define the input variables
filepath_TRANSP = "" # for example /Users/henrikj/codes/OWCF/TRANSP/JET/94701V01/94701V01.cdf
filepath_TRANSP_FI = "" # for example /Users/henrikj/codes/OWCF/TRANSP/JET/96100J01/96100J01_fi_1.cdf
filepath_equil = "" # for example g96100_0-53.0012.eqdsk
SPECIES = "D" # for example D, T, he3
folderpath_o = "" # The path to the desired output folder. Remember to finish path with '/'
timepoint = nothing # Format "XX,YYYY" where XX are seconds and YYYY are decimals. If unknown, leave as nothing
verbose = true

## --------------------------------------------------------------------------
verbose && println("Loading Julia packages... ")
using PyCall # For using Python code in Julia
using JLD2 # To write/open .jld2-files (Julia files, basically)
using FileIO # To write/open files in general
using EFIT # To handle EFIT structures
using Equilibrium # To handle tokamak equilibria
using NetCDF # To easily load the shot timepoint of the TRANSP fast-ion distribution
pushfirst!(PyVector(pyimport("sys")."path"), "") # To import DRESS code Python files

## --------------------------------------------------------------------------
verbose && println("Loading Python packages... ")
py"""
import h5py
import numpy as np
from netCDF4 import Dataset

import forward
import transp_dists
import transp_output
import vcone
"""

## --------------------------------------------------------------------------
verbose && println("Loading TRANSP id information... ")
TRANSP_id = (split((split(filepath_TRANSP,"."))[1],"/"))[end] # Assume part between last '/' and '.' is the TRANSP id

## --------------------------------------------------------------------------
verbose && println("Loading tokamak equilibrium... ")
if ((split(filepath_equil,"."))[end] == "eqdsk") || ((split(filepath_equil,"."))[end] == "geqdsk") 
    M, wall = read_geqdsk(filepath_equil,clockwise_phi=false) # Assume counter-clockwise phi-direction
    jdotb = M.sigma # The sign of the dot product between the plasma current and the magnetic field

    # Extract timepoint information from .eqdsk/.geqdsk file
    eqdsk_array = split(filepath_equil,".")
    XX = (split(eqdsk_array[end-2],"-"))[end] # Assume format ...-XX.YYYY.eqdsk where XX are the seconds and YYYY are the decimals
    YYYY = eqdsk_array[end-1] # Assume format ...-XX.YYYY.eqdsk where XX are the seconds and YYYY are the decimals
    timepoint = (timepoint == nothing ? XX*","*YYYY : timepoint) # Format XX,YYYY to avoid "." when including in filename of saved output
else # Otherwise, assume magnetic equilibrium is a saved .jld2 file
    myfile = jldopen(filepath_equil,false,false,false,IOStream)
    M = myfile["S"]
    wall = myfile["wall"]
    close(myfile)
    jdotb = (M.sigma_B0)*(M.sigma_Ip)
    timepoint = (timepoint == nothing ? "00,0000" : timepoint)
end

## --------------------------------------------------------------------------
verbose && println("Defining R_array... ")
da = round((maximum(wall.r)-M.axis[1])/1000,sigdigits=1) # Use da to start ρ_pol grid not exactly on the magnetic axis
dw = round((maximum(wall.r)-M.axis[1])/20,sigdigits=1) # Use dw to end ρ_pol grid approximately at the separatrix. Increase 20 for closer to LFS wall, decrease 20 for further from the LFS wall.
R_array = collect(range(M.axis[1]+da,maximum(wall.r)-dw, length=50))

## --------------------------------------------------------------------------
verbose && println("Defining z_array... ")
z_array = (M.axis[2]) .* ones(length(R_array)) # Equitorial plane
verbose && println("Computing rho_p_array from R_array, z_array... ")
rho_p_array = zeros(length(R_array))
for i=1:length(R_array)
    rho_p_array[i] = Equilibrium.rho_p(M,R_array[i],z_array[i])
end

## --------------------------------------------------------------------------
py"""
# Load TRANSP simulation data
print("Extracting TRANSP data... ")
tr_out = transp_output.TranspOutput($TRANSP_id, step=1,
                                    out_file=$filepath_TRANSP,
                                    fbm_files=[$filepath_TRANSP_FI])
print("Setting bulk thermal distribution... ")
bulk_dist = transp_dists.Thermal(tr_out, ion=$SPECIES)
"""

## --------------------------------------------------------------------------
py"""
print("Calculating T_bulk from R_array, z_array... ")
T_bulk = bulk_dist.get_temperature($R_array, $z_array)
print("Calculating n_bulk from R_array, z_array... ")
n_bulk = bulk_dist.get_density($R_array, $z_array)
"""

## --------------------------------------------------------------------------
TIME = round((ncread(filepath_TRANSP_FI,"TIME"))[1],digits=4)
TIME_array = split("$(TIME)",".") # Will be on format XX.YYYY
XX = TIME_array[1]
YYYY = TIME_array[2]
timepoint = XX*","*YYYY # Format XX,YYYY to avoid "." when including in filename of saved output

## --------------------------------------------------------------------------
verbose && println("Saving rho_p_array, T_bulk and n_bulk... ")
global filepath_output_orig = folderpath_o*TRANSP_id*"_at"*timepoint*"s_"*SPECIES*"_TempNDens"
global filepath_output = deepcopy(filepath_output_orig)
global count = 1
while isfile(filepath_output*".jld2") # To take care of not overwriting files. Add _(1), _(2) etc
    global filepath_output = filepath_output_orig*"_($(Int64(count)))"
    global count += 1 # global scope, to surpress warnings
end
global filepath_output = filepath_output*".jld2"
myfile = jldopen(filepath_output,true,true,false,IOStream)
write(myfile,"rho_pol",rho_p_array)
write(myfile,"thermal_temp",(py"T_bulk")[:])
write(myfile,"thermal_dens",(py"n_bulk")[:])
close(myfile)

println("~~~~~~~Done!~~~~~~")
