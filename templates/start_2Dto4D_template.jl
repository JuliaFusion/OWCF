################################ start_2Dto4D_template.jl ##################################################
# This file contains all the inputs that the script orbWeights_2Dto4D.jl needs to inflate the 2D compact form
# weight matrix into its inflated 4D form. This file also executes the script
# orbWeights_2Dto4D.jl after the inputs are defined. The inputs are as follows:
#
# batch_job - If true, the script assumes that it will be executed via an HPC batch job. The script will act accordingly. - Bool
# distributed - If true, parallel computing with multiple CPU cores will be used. - Bool
# folderpath_OWCF - The path to the OWCF folder - String
# numOcores - The number of CPU cores that will be used if distributed is set to true. - Int64
#
# diagnostic_name - The name of the diagnostic of the orbit weights (esthetic) - String
# Ed_energies - Can be type Nothing, Int64 or Array{Float64,1}. 
#   # Nothing  - If Ed_energies is specified as 'nothing', all orbit weight functions (all diagnostic energy bins)
#                will be transformed from 1D to 3D. They will be put together into a 4D matrix.
#   # Int64 - If Ed_energies is specified as an Int64 (for example '10'), the algorithm will transform that many
#                orbit weight functions from 1D to 3D. They will be (approximately) evenly spaced across all diagnostic energy bins.
#                If Ed_energies is specified as an Int64, it cannot be greater than half of the number of rows of the weight matrix.
#                For example, if the weight matrix has 80 rows (diagnostic energy bins), then Ed_energies must be <= 40.
#                Continuing the example, if Ed_energies is specified as '10', then the orbit weight functions corresponding to 
#                weight matrix rows 1, 9, 17, 25, 33, 41, 49, 57, 65, 73 will be transformed from 1D to 3D, and then put 
#                together into a 4D matrix.
#   # Array{Float64,1} - If Ed_energies is specified as an Array{Float64,1}, the algorithm will attempt to find the 
#                diagnostic energy bins with the values closest to the values in Ed_energies. If several values are close 
#                to the same diagnostic energy bin, no duplicates will be produced (the weight function will be transformed only once).  
# filepath_equil - The path to the .eqdsk-file/.jld2-file with the tokamak magnetic equilibrium and geometry - String
# filepath_W - The path to the .jld2 weights file, containing orbit weights to be enflated - String
# FI_species - The ion species for which the orbit weights are computed. E.g. "D", "T" etc. See misc/species_func.jl for further info - String
# timepoint - The timepoint of the tokamak shot for the magnetic equilibrium. Format XX,YYYY where XX are seconds and YYYY are decimals - String
# tokamak - The tokamak of the orbit weights - String
# TRANSP_id - The TRANSP id to extract tokamak/equilibrium data from - String
# verbose - If set to true, the script will talk a lot! - Bool
#
# extra_kw_args - Here you can define extra keyword arguments to be used when integrating the equations-of-motion
# for the guiding-center orbits. Please consult the OWCF manual for further info on this.
#
# Script written by Henrik Järleblad. Last maintained 2022-10-11.
############################################################################################################

## First you have to set the system specifications
using Distributed # Needed, even though distributed might be set to false. This is to export all inputs to all workers right away, if needed.
batch_job = false
distributed = true
folderpath_OWCF = "" # OWCF folder path. Finish with '/'
numOcores = 4 # When executing script via HPC cluster job, make sure you know how many cores you have requested for your batch job

## Navigate to the OWCF folder and activate the OWCF environment
cd(folderpath_OWCF)
using Pkg
Pkg.activate(".")

## If running as a batch job on a SLURM CPU cluster
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

## If running locally and multi-threaded
if !batch_job && distributed # Assume you are executing the script on a local laptop (/computer)
    println("Adding processes... ")
    addprocs(numOcores-(nprocs()-1)) # If you didn't execute this script as an HPC cluster job, then you need to add processors like this. Add all remaining available cores.
    # The '-(nprocs()-1)' part is simply to ensure to extra processes are added, in case script needs to be restarted on a local computer
end

## -----------------------------------------------------------------------------
# Inputs
@everywhere begin
    diagnostic_name = "" # TOFOR, AB etc
    Ed_energies = nothing # Specific diagnostic energy bin levels to try to extract. Can be type Nothing, Int64 or Array{Float64,1} (see script description above)
    filepath_equil = "" # for example, g96100_0-53.0012.eqdsk
    folderpath_OWCF = $folderpath_OWCF # Set path to OWCF folder to same as main process (hence the '$')
    filepath_W = ""
    FI_species = "" # D, T, he3 etc
    timepoint = nothing # If unknown, just leave as nothing. The algorithm will try to figure it out automatically.
    tokamak = "" # JET, DIIID, ITER etc. Can be left blank
    TRANSP_id = "" # for example, 94701J01. Can be left blank

    verbose = false

    # EXTRA KEYWORD ARGUMENTS BELOW (these will go into the orbit_grid() and get_orbit() functions from GuidingCenterOrbits.jl)
    extra_kw_args = Dict(:limit_phi => true, :max_tries => 0)
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
include("helper/orbWeights_2Dto4D.jl")

## -----------------------------------------------------------------------------
# Then you clean up after yourself
if batch_job && distributed
    for i in workers()
        rmprocs(i)
    end
end
