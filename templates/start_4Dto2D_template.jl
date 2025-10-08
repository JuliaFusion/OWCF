################################ start_4Dto2D_template.jl ###################################################
# This file contains all the inputs that the script orbWeights_4Dto2D.jl needs to compress the 4D inflated form
# of orbit weight functions into their compressed 2D form. This file also executes the script
# orbWeights_4Dto2D.jl after the inputs are defined. The inputs are as follows:
#
# batch_job_SLURM - If true, the script assumes that it will be executed as part of an high-performance computational 
#                   cluster batch job within the SLURM cluster environment. The script will act accordingly. - Bool
# distributed - If true, parallel computing with multiple CPU cores will be used. - Bool
# numOcores - The number of CPU cores that will be used if distributed is set to true. - Int64
#
# diagnostic - The diagnostic of the orbit weights - String
# filepath_equil - The path to the .eqdsk-file (or .geqdsk/.jld2-file) with the tokamak magnetic equilibrium and geometry - String
# filepath_W - The path to the .jld2 weights file, containing orbit weights to be compressed - String
# folderpath_OWCF - The path to the OWCF folder - String
# FI_species - The ion species for which the orbit weights are computed. E.g. "D", "T" etc. See misc/species_func.jl for further info - String
# timepoint - The timepoint of the tokamak shot for the magnetic equilibrium. Format XX,YYYY where XX are seconds and YYYY are decimals - String
# tokamak - The tokamak of the orbit weights - String
# TRANSP_id - The TRANSP id to extract tokamak/equilibrium data from - String
# verbose - If set to true, the script will talk a lot - Bool
#
# extra_kw_args - Here you can define extra keyword arguments to be used when integrating the equations-of-motion
# for the guiding-center orbits. Please consult the OWCF manual for further info on this.
#
# Script written by Henrik Järleblad. Last maintained 2025-10-07.
#############################################################################################################

## First you have to set the system specifications
using Distributed # Needed, even though distributed might be set to false. This is to export all inputs to all workers right away, if needed.
batch_job_SLURM = false
distributed = true
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
    diagnostic = "" # TOFOR, AB etc
    filepath_equil = "" # for example, g96100_0-53.0012.eqdsk
    filepath_W = ""
    folderpath_OWCF = ""

    FI_species = "" # D, T, he3 etc
    timepoint = nothing # If unknown, just leave as nothing. The algorithm will try to figure it out automatically.
    tokamak = "" # JET, DIIID, ITER etc
    TRANSP_id = "" # for example, 94701J01

    verbose = false

    # EXTRA KEYWORD ARGUMENTS BELOW (these will go into the orbit_grid() and get_orbit() functions from GuidingCenterOrbits.jl)
    extra_kw_args = Dict(:limit_phi => true, :toa => true)
    # limits the number of toroidal turns for orbits
    # toa forces the orbit integration algorithm to only use adaptive integration schemes
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
include("helper/orbWeights_4Dto2D.jl")

## -----------------------------------------------------------------------------
# Then you clean up after yourself
if batch_job_SLURM || distributed
    for i in workers()
        rmprocs(i)
    end
end