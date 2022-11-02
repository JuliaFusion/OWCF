################################################## F_os_1Dto3D.jl ##################################################

#### Description:
# This script is essentially an unnecessary complicated way of perfoming map_orbits(og, F_os). But it was realized
# that such a script would be useful to an OWCF user, if nothing else but to showcase how it should be done correctly.
# If ps2os.jl was run without the input variable include1Dto3D set to true, the F_os_1Dto3D.jl script can be used
# to inflate the 1D compressed fast-ion distribution vector into its full 3D inflated form.

#### Inputs:
# folderpath_OWCF - The path to the OWCF folder on your computer - String
# 
# batch_job - If true, then script will assume you're executing via a high-performance computing (HPC) cluster job - Boolean
# distributed - If true, parallel computing will be used - Boolean
# numOcores - The number of cores that you have specified for your HPC cluster job. If batch_job==true - Int64
#
# filepath_F_os - The path to the file with the compressed 1D fast-ion orbit-space distribution data. Should be an output file from ps2os.jl - String
# filepath_equil - The path to a .eqdsk/.geqdsk/.jld2 file containing the magnetic equilibrium and wall data - String
# FI_species - The fast-ion species of the fast-ion distribution. E.g. "D", "T", "3he" etc - String
# folderpath_o - The path to the folder where you want your output. Always finish all paths to folders with '/' - String
# timepoint - The tokamak shot timepoint of the magnetic equilibrium. Format XX,YYYY where XX are seconds and YYYY are decimals - String
# verbose - If set to true, the F_os_1Dto3D.jl helper script will talk a lot! - Bool

#### Outputs:
# -

#### Other:
# -

# Script written by Henrik Järleblad. Last maintained 2022-08-27.
####################################################################################################################

## --------------------------------------------------------------------------
# Please specify the OWCF folder and let the app change directory to the 
# OWCF folder when F_os_1Dto3D.jl is executed. This is to be able to load the
# correct versions of the Julia packages as specified in the Project.toml and 
# Manifest.toml files.
folderpath_OWCF = "G:/My Drive/DTU/codes/OWCF/" # Finish with '/'
cd(folderpath_OWCF)
using Pkg
Pkg.activate(".")

## --------------------------------------------------------------------------
# Please specify system inputs variables. The F_os_1Dto3D script can be run as part
# of a batch job on a computational SLURM cluster by setting batch_job=true. If so,
# please also specify the correct amount of CPU cores with the numOcores input variable.
using Distributed # Needed, even though distributed might be set to false. This is to export all inputs to all workers right away, if needed.
batch_job = false
distributed = true
numOcores = 1 # When executing script via HPC cluster job, make sure you know how many cores you have requested for your batch job

if batch_job && distributed
    # Load the SLURM CPU cores
    using ClusterManagers
    addprocs(SlurmManager(numOcores))
    hosts = []
    pids = []
    for i in workers()
        host, pid = fetch(@spawnat i (gethostname(), getpid()))
        push!(hosts, host)
        push!(pids, pid)
    end
    @show hosts
end

if !batch_job && distributed # Assume you are executing the script on a local laptop (/computer)
    println("Adding processes... ")
    addprocs(abs(nprocs() - length(Sys.cpu_info()) )) # If you didn't execute this script as an HPC cluster job, then you need to add processors like this. Add all remaining available cores.
end

## --------------------------------------------------------------------------
# Must load JLD2 package first, to be able to check filepath_F_os for 'extra_kw_args'
using JLD2

## --------------------------------------------------------------------------
## Defining input variables
filepath_F_os = "" # The filepath to the 1D (raw) orbit-space distribution (could be a ps2os_results file)
filepath_equil = "" # The filepath to the magnetic equilibrium
FI_species = "" # D, p, T, 3he etc
folderpath_o = "" # The path to the folder where you want your output
timepoint = nothing # Format "XX,YYYY" where XX are seconds and YYYY are decimals. If unknown, left as nothing. The algorithm will try to figure it out automatically
verbose = false

# EXTRA KEYWORD ARGUMENTS BELOW (these will go into the orbit_grid() and get_orbit() functions from GuidingCenterOrbits.jl)
myfile = jldopen(filepath_F_os,false,false,false,IOStream)
if haskey(myfile,"extra_kw_args")
    extra_kw_args = myfile["extra_kw_args"]
else
    verbose && println("No extra keyword arguments for orbit integration found in fast-ion distribution file. Assuming :toa, :limit_phi and :maxiter=0")
    extra_kw_args = Dict(:toa => true, :limit_phi => true, :maxiter => 0)
    # toa is 'try only adaptive'
    # limits the number of toroidal turns for orbits
    # The orbit integration algorithm will try progressively smaller timesteps these number of times
end
close(myfile)

## --------------------------------------------------------------------------
## Loading packages
verbose && println("Loading JULIA packages and OWCF dependencies... ")
@everywhere begin
    using EFIT # For equilibrium purposes
    using JLD2 # For saving variables to file. For back-up.
    using FileIO # -||-
    using Equilibrium # For plasma equilibrium loading (read_geqdsk)
    using GuidingCenterOrbits # For orbit and Jacobian calculations
    using OrbitTomography # That's what this all is 'bout!
    include(folderpath_OWCF*"extra/dependencies.jl")
    include(folderpath_OWCF*"misc/species_func.jl") # To convert species labels to particle mass
end

## --------------------------------------------------------------------------
verbose && println("Loading fast-ion orbit-space distribution in compressed raw 1D form... ")
myfile_F_os = jldopen(filepath_F_os,false,false,false,IOStream)
F_os_raw = myfile_F_os["F_os_raw"]
nfast = myfile_F_os["nfast"]
E_array = myfile_F_os["E_array"]
pm_array = myfile_F_os["pm_array"]
Rm_array = myfile_F_os["Rm_array"]
close(myfile_F_os)

## --------------------------------------------------------------------------
verbose && println("Computing correctly scaled fast-ion distribution from total number of fast ions... ")
F_os = (nfast/sum(F_os_raw)) .*F_os_raw

## --------------------------------------------------------------------------
# Loading tokamak equilibrium
verbose && println("Loading tokamak equilibrium... ")
if ((split(filepath_equil,"."))[end] == "eqdsk") || ((split(filepath_equil,"."))[end] == "geqdsk") 
    M, wall = read_geqdsk(filepath_equil,clockwise_phi=false) # Assume counter-clockwise phi-direction
    jdotb = M.sigma # The sign of the dot product between the plasma current and the magnetic field

    # Extract timepoint information from .eqdsk/.geqdsk file
    eqdsk_array = split(filepath_equil,".")
    XX = (split(eqdsk_array[end-2],"-"))[end] # Assume format ...-XX.YYYY.eqdsk where XX are the seconds and YYYY are the decimals
    YYYY = eqdsk_array[end-1] # Assume format ...-XX.YYYY.eqdsk where XX are the seconds and YYYY are the decimals
    timepoint = XX*","*YYYY # Format XX,YYYY to avoid "." when including in filename of saved output
else # Otherwise, assume magnetic equilibrium is a saved .jld2 file
    myfile = jldopen(filepath_equil,false,false,false,IOStream)
    M = myfile["S"]
    wall = myfile["wall"]
    close(myfile)
    jdotb = (M.sigma_B0)*(M.sigma_Ip)

    timepoint = (timepoint == nothing ? "00,0000" : timepoint)
end

filepath_F_array = split(filepath_F_os,"_")
time_string = filepath_F_array[end-2] # Assume this to start with
if (time_string[1]=='a') && (time_string[2]=='t') && (time_string[end]=='s')
    # Correct format identified. Superseeds user input and/or timepoint found via equilibrium file
    timepoint = time_string[3:end-1] # First to elements are "at" and last element is "s"
else
    timepoint = (timepoint == nothing ? "00,0000" : timepoint)
end

## --------------------------------------------------------------------------
verbose && println("Calculating orbit grid... ")
og_orbs, og = OrbitTomography.orbit_grid(M, E_array, pm_array, Rm_array; amu=getSpeciesAmu(FI_species), q=getSpeciesEcu(FI_species), wall=wall, extra_kw_args... )
og_orbs = nothing # Minimize memory

## --------------------------------------------------------------------------
verbose && println("Inflating 1D distribution to full 3D form... ")
F_os_3D = map_orbits_OWCF(og, F_os, false) # Don't assume equidistant orbit-space grid

## --------------------------------------------------------------------------
verbose && println("Done! Saving... ")
myfile = jldopen(folderpath_o*"F_os_3D"*"_at"*timepoint*"s_"*FI_species*"_"*"$(length(E_array))x$(length(pm_array))x$(length(Rm_array)).jld2",true,true,false,IOStream)
write(myfile, "F_os_3D", F_os_3D)
write(myfile, "nfast_orig", nfast)
write(myfile, "E_array", E_array)
write(myfile, "pm_array", pm_array)
write(myfile, "Rm_array", Rm_array)
close(myfile)

## --------------------------------------------------------------------------
verbose && println("~~~~~~~~ F_os_1Dto3D completed successfully! ~~~~~~~~")
