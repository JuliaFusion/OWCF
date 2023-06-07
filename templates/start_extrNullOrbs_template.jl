################################ start_extrNullOrbs_template.jl #########################################

#### Description:
# This script contains all the inputs for the extrNullOrbs.jl script. It defines the computational
# environment (batch job, cores etc), loads the inputs and then executes the extrNullOrbs.jl script.

#### Inputs
# batch_job - If true, the script assumes that it will be executed via an HPC batch job. The script will act accordingly. - Bool
# distributed - If true, parallel computing with multiple CPU cores will be used. - Bool
# folderpath_OWCF - The path the OWCF-folder on your computer - String
# numOcores - The number of CPU cores that will be used if distributed is set to true. - Int64
#
# signalNullThreshold - A fraction between 0.0 and 1.0. Below this fraction of maximum(S), the diagnostic signal
#                 is counted as zero and null measurement regions will be examined - Float64
# orbNullThreshold - A fraction between 0.0 and 1.0. Below this fraction of maximum(W(E_d)), the orbit is 
#                 counted as a null orbit - Float64
# filepath_W - The file path to the 2D orbit weight functions .jld2 file from which to extract the
#              null orbits. It's the main output file of calcOrbWeights.jl. The rows of the 
#              weight matrix must be equal to the diagnostic energy bins of the signal file. - String
# filepath_S - The file path to the signal which is to be examined for zero-valued elements. The
#              diagnostic energy bins used to compute the signal must be equal to the diagnostic
#              energy bins corresponding to the rows of the weight matrix loaded from filepath_W.
#              filepath_S can be the output file from calcSpec.jl or ps2WF.jl. - String
# FI_species - The fast-ion species of the orbit space being used - String
# folderpath_o - The name of the output folder to which to save the results - String
# include2Dto4D - If true, then orbit-space indices of null orbits will be computed and saved as well - Bool
# filepath_equil - If include2Dto4D, you will need to specify an .eqdsk (or .geqdsk/.jld2) file to load magnetic equilibrium from - String
# timepoint - The timepoint of the tokamak shot for the magnetic equilibrium. Format XX,YYYY where XX are seconds and YYYY are decimals - String
# verbose - If true, then the script will be very talkative! - Bool
#
# extra_kw_args - Here you can define extra keyword arguments to be used when integrating the equations-of-motion
# for the guiding-center orbits. Please consult the OWCF manual for further info on this.

#### Outputs
# -

#### Saved files
# nullOrbs_[tokamak]_[TRANSP ID]_[diagnostic]_[FI species]_[nE]x[npm]x[nRm].jld2
#   This saved file will have the fields:
#   nullOrbs - A 1D array where '1.0' means a null orbit, and '0.0' is not. Same ordering of orbits as that of F_os_raw from ps2os.jl - Array{Float64,1}
#   if include2Dto4D
#       nullOrbs_indices - A vector containing cartesian indices marking the orbit-space location of the null orbits - Vector{CartesianIndex{3}}
#   E_array - The fast-ion energy grid array used for orbit space - Array{Float64,1}
#   pm_array - The fast-ion pm grid array used for orbit space - Array{Float64,1}
#   Rm_array - The fast-ion Rm grid array used for orbit space - Array{Float64,1}
#   filepath_W - The filepath to the orbit weights file used - String
#   filepath_S - The filepath to the signal file used - String 

# Script written by Henrik Järleblad. Last maintained 2022-10-11.
########################################################################################

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
@everywhere begin
    signalNullThreshold = 0.0
    orbNullThreshold = 0.0
    filepath_W = ""
    filepath_S = ""
    folderpath_OWCF = $folderpath_OWCF # Set path to OWCF folder to same as main process (hence the '$')
    FI_species = "" # D, p, T, 3he etc
    folderpath_o = ""
    (include2Dto4D = false) && (filepath_equil = "")
    timepoint = nothing # If unknown, just leave as nothing. The algorithm will try to figure it out automatically.
    verbose = true

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
include("helper/extractNullOrbits.jl")

## -----------------------------------------------------------------------------
# Then you clean up after yourself
if batch_job && distributed
    for i in workers()
        rmprocs(i)
    end
end
