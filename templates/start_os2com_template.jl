################################ start_os2com_template.jl ##################################################
# This file contains all the inputs that the script os2com.jl needs to transform an orbit weight
# function matrix from (E,pm,Rm) coordinates to (E,mu,Pphi;alpha) coordinates. This can also be 
# achieved by simply utilizing the os2COM() function in extra/dependencies.jl. However, this script
# let's you avoid having to think about it, and just does it for you. It also works with topological
# maps or boundarie, as well as poloidal and toroidal transit times, and null orbits.
#
# batch_job - If true, the script assumes that it will be executed via an HPC batch job. The script will act accordingly. - Bool
# distributed - If true, parallel computing with multiple CPU cores will be used. - Bool
# numOcores - The number of CPU cores that will be used if distributed is set to true. - Int64
#
# filepath_equil - The path to the magnetic equilibrium file (.eqdsk/.geqdsk/.jld2 file format) - String
# filepath_Q - The path to the 3D or 4D (E,pm,Rm) orbit-space quantity - String
# folderpath_o - The path to the folder in which you want the os2com.jl outputs - String
# folderpath_OWCF - The path to where the OWCF-folder is saved on your computer - String
# nmu - The number of gridpoints in magnetic moment - Int64
# nPphi - The number of gridpoints in toroidal canonical angular momentum - Int64
# timepoint - The timepoint of the tokamak shot for the magnetic equilibrium. Format XX,YYYY where XX are seconds and YYYY are decimals - String
# tokamak - The name of the tokamak. If unknown, leave as "" - String
# verbose - If set to true, the os2com.jl script will talk a lot! - String
# FI_species - The fast-ion species. E.g. "D", "T", "3he", "alpha" etc - String
# 
# extra_kw_args - Here you can define extra keyword arguments to be used when integrating the equations-of-motion
# for the guiding-center orbits. Please consult the OWCF manual for further info on this.

# Script written by Henrik Järleblad. Last maintained 2022-10-11.
#############################################################################################################

## First you have to set the system specifications
using Distributed # Needed, even though distributed might be set to false. This is to export all inputs to all workers right away, if needed.
batch_job = false
distributed = true
numOcores = 4 # When executing script via HPC cluster job, make sure you know how many cores you have requested for your batch job

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
    addprocs(numOcores-(nprocs()-1)) # If you didn't execute this script as an HPC cluster job, then you need to add processors like this. Add all remaining available cores.
    # The '-(nprocs()-1)' part is simply to ensure to extra processes are added, in case script needs to be restarted on a local computer
end

## -----------------------------------------------------------------------------
# Inputs
@everywhere begin
    filepath_equil = ""
    filepath_Q = "" 
    folderpath_o = ""
    folderpath_OWCF = ""
    nmu = 1
    nPphi = 1
    timepoint = nothing # If unknown, just leave as nothing. The algorithm will try to figure it out automatically.
    tokamak = ""

    verbose = false

    FI_species = "D"

    # EXTRA KEYWORD ARGUMENTS BELOW (these will go into the orbit_grid() and get_orbit() functions from GuidingCenterOrbits.jl)
    extra_kw_args = Dict(:toa => true, :limit_phi => true, :maxiter => 0)
    # toa is 'try only adaptive'
    # limits the number of toroidal turns for orbits
    # The orbit integration algorithm will try progressively smaller timesteps these number of times
end

## -----------------------------------------------------------------------------
# Change directory to OWCF-folder. Activate the environment there to ensure correct package versions
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
if batch_job && distributed
    for i in workers()
        rmprocs(i)
    end
end
