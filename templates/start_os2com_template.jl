################################ start_os2com_template.jl ##################################################
# This file contains all the inputs that the script os2com.jl needs to transform a fast-ion distribution
# from (E,pm,Rm) coordinates to (E,Λ,Pϕ_n;σ) coordinates. This can also be 
# achieved by simply utilizing the os2COM() function in extra/dependencies.jl. However, this script
# let's you avoid having to think about it, and just does it for you. It also works with weight functions,
# topological maps or boundaries, as well as poloidal and toroidal transit times, and null orbits. Λ is the 
# normalized magnetic moment given by Λ=μ*B0/E where μ is the magnetic moment, B0=B(mag. axis) and E is the 
# energy. Pϕ_n is the normalized toroidal canonical momentum, given by Pϕ_n=Pϕ/(q*|Ψ_w|) where q is the charge 
# of the fast ion and Ψ_w is the poloidal magnetic flux at the last closed flux surfaces (LCFS). If Ψ(LCFS)==0 
# for some reason (e.g. due to some convention), Ψ_w=Ψ(mag. axis) is assumed instead. σ is a binary coordinate. 
# σ=-1 (Julia index 1) corresponds to counter-current orbits. σ=+1 (Julia index 2) corresponds to co-current orbits.
#
# batch_job_SLURM - If true, the script assumes that it will be executed as part of an high-performance computational 
#                   cluster batch job within the SLURM cluster environment. The script will act accordingly. - Bool
# distributed - If true, parallel computing with multiple CPU cores will be used. - Bool
# folderpath_OWCF - The path to where the OWCF-folder is saved on your computer - String
# numOcores - The number of CPU cores that will be used if distributed is set to true. - Int64
#
# filepath_equil - The path to the magnetic equilibrium file (.eqdsk/.geqdsk/.jld2 file format) - String
# filepath_Q - The path to the 3D or 4D (E,pm,Rm) orbit-space quantity. Should be an output file from 
#              calcOrbWeights.jl, calcTopoMap.jl, F_os_1Dto3D.jl etc. - String
# folderpath_o - The path to the folder in which you want the os2com.jl outputs - String
# nmu - The number of gridpoints in normalized magnetic moment - Int64
# nPphi - The number of gridpoints in normalized toroidal canonical angular momentum - Int64
# timepoint - The timepoint of the tokamak shot for the magnetic equilibrium. Format XX,YYYY where XX are seconds and YYYY are decimals - String
# tokamak - The name of the tokamak. If unknown, leave as "" - String
# verbose - If set to true, the os2com.jl script will talk a lot! - String
# FI_species - The fast-ion species. E.g. "D", "T", "3he", "alpha" etc - String
# 
# extra_kw_args - Here you can define extra keyword arguments to be used when integrating the equations-of-motion
# for the guiding-center orbits. Please consult the OWCF manual for further info on this.

# Script written by Henrik Järleblad. Last maintained 2025-10-07.
#############################################################################################################

## First you have to set the system specifications
using Distributed # Needed, even though distributed might be set to false. This is to export all inputs to all workers right away, if needed.
batch_job_SLURM = false
distributed = true
folderpath_OWCF = "" # OWCF folder path. Finish with '/'
numOcores = 4 # When executing the script as part of an HPC cluster batch job, the number of CPU cores will be detected and set automatically (i.e. the value of the numOcores variable does not matter)

## Navigate to the OWCF folder and activate the OWCF environment
cd(folderpath_OWCF)
using Pkg
Pkg.activate(".")

## If running as a batch job on a SLURM HPC cluster with CPUs
if batch_job_SLURM
    # Load the SLURM package and add the CPU cores
    using SlurmClusterManager
    addprocs(SlurmManager()) # The number of CPU cores are detected automatically
    @everywhere println("X-wing $(rand(["Red","Green","Blue"])) $(myid()) reporting for duty at $(gethostname())")
end

## If running locally and multi-threaded
if !batch_job_SLURM && distributed # Assume you are executing the script on a local laptop (/computer)
    println("Adding processes... ")
    addprocs(numOcores-(nprocs()-1)) # If you didn't execute this script as an HPC cluster job, then you need to add processors like this. Add all remaining available cores.
    # The '-(nprocs()-1)' part is simply to ensure to extra processes are added, in case script needs to be restarted on a local computer
end

## -----------------------------------------------------------------------------
# Inputs
@everywhere begin
    filepath_equil = ""
    filepath_Q = "" 
    folderpath_o = ""
    folderpath_OWCF = $folderpath_OWCF # Set path to OWCF folder to same as main process (hence the '$')
    nmu = 1
    nPphi = 1
    timepoint = nothing # If unknown, just leave as nothing. The algorithm will try to figure it out automatically.
    tokamak = ""

    verbose = false

    FI_species = "D"

    # EXTRA KEYWORD ARGUMENTS BELOW (these will go into the orbit_grid() and get_orbit() functions from GuidingCenterOrbits.jl)
    extra_kw_args = Dict(:limit_phi => true, :maxiter => 0)
    # limits the number of toroidal turns for orbits
    # The orbit integration algorithm will try progressively smaller timesteps these number of times
end

## -----------------------------------------------------------------------------
# Change directory to OWCF-folder on all external processes. Activate the environment there to ensure correct package versions
# as specified in the Project.toml and Manuscript.toml files. Do this for every 
# CPU processor (@everywhere)
@everywhere begin
    using Pkg
    cd(folderpath_OWCF)
    Pkg.activate(".")
end

## -----------------------------------------------------------------------------
# Then you execute the script
include("helper/os2com.jl")

## -----------------------------------------------------------------------------
# Then you clean up after yourself
if batch_job_SLURM || distributed
    for i in workers()
        rmprocs(i)
    end
end