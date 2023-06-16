################################ start_ps2os_template.jl #########################################
# This file contains all the inputs needed by the script ps2os.jl. Copy this file, and change
# inputs accordingly. After reading the inputs, the script will execute ps2os.jl. As of this version
# of the OWCF, the ps2os.jl script will use Monte-Carlo sampling to transform from (E,p,R,z) to 
# (E,pm,Rm). In future version, other methods of transforming are planned to be selectable.
#
### The inputs are:
#
# batch_job - If true, the script assumes that it will be executed via an HPC batch job. The script will act accordingly. - Bool
# distributed - If true, parallel computing with multiple CPU cores will be used. - Bool
# folderpath_OWCF - The path to the OWCF folder - String
# numOcores - The number of CPU cores that will be used if distributed is set to true. - Int64
#
# include1Dto3D - If true, then the 3D OS FI distribution will be saved as well - Boolean
# filepath_distr - The path to the fast-ion distribution file to extract tokamak/fast-ion data from. Must be .h5 or .jld2 file format - String
# filepath_equil - The path to the .eqdsk-file (or .geqdsk/.jld2-file) with the tokamak magnetic equilibrium and geometry - String
# filepath_W - The path to the weights file to use as blueprint for the orbit grid. If desired. - String
# filepath_OG - The path to the output file from calcOrbGrid.jl. If specified, ps2os.jl will use that as the (E,pm,Rm) grid points - String
# flip_F_EpRz_pitch - If true, then the FI PS distribution will be mirrored in pitch - Boolean
# folderpath_o - Path to folder where you want your results. End with '/' (or '\' if Windows...) - String'
# numOsamples - The number of orbit samples you want. The more, the better the transform. But takes longer to compute! - Int64
# nbatch - The algorithm will save the progress every nbatch sample. Useful if the sampling gets cancelled unexpectedly - Int64
# os_equidistant - If true, script assumes equidistant orbit grid. False not supported yet - Boolean
# timepoint - The timepoint of the tokamak shot for the magnetic equilibrium. Format XX,YYYY where XX are seconds and YYYY are decimals - String
# tokamak - The tokamak for which the topological map is created - String
# TRANSP_id - The TRANSP id to extract tokamak/equilibrium data from - String
# useWeightsFile - If true, ps2os.jl will use filepath_W to load info to create orbit grid - Boolean
# useOrbGridFile - If true, ps2os.jl will use filepath_OG to load into to create orbit grid - Boolean
# verbose - If true, lots of outputs will be printed during script execution - Boolean
# visualizeProgress - If true, then progress bar will be shown. Good for local calcs. Bad for clusters. - Boolean
#
## Then there are input variables that you only need to specify is you set useWeightsFile/useOrbGridFile to false:
# Emin - The lower bound for the fast-ion energy in orbit space - Float64
# Emax - The upper bound for the fast-ion energy in orbit space - Float64
# FI_species - The particle species of the fast-ion distribution. E.g. "D", "T", "3he" etc - String
# inclPrideRockOrbs - If true, then pride rock (counter-stagnation) orbits will be included
# nE - The number of fast-ion energy grid points in orbit space - Int64
# npm - The number of fast-ion pm grid points in orbit space - Int64
# nRm - The number of fast-ion Rm grid points in orbit space - Int64
# pm_min - The lower bound for the pitch maximum in orbit space - Float64
# pm_max - The upper bound for the pitch maximum in orbit space - Float64
# Rm_min - The lower boundary for the orbit-grid Rm values - Float64
# Rm_max - The upper boundary for the orbit-grid Rm values - Float64
#
## Then, finally, you could specify to have the fast-ion interpolated onto a coarser/finer
# (E,p,R,z) grid before transforming it (via sampling) to orbit space. This is useful to 
# test robustness and the effect of grid resolution etc. The inputs are as follows
# interp - If set to true, the algorithm will interpolate the loaded fast-ion distribution - Bool
# nE_ps - The number of fast-ion energy grid points in particle space to interpolate onto - Int64
# np_ps - The number of fast-ion pitch grid points in particle space to interpolate onto - Int64
# nR_ps - The number of fast-ion R grid points in particle space to interpolate onto - Int64
# nz_ps - The number of fast-ion z grid points in particle space to interpolate onto - Int64

### Other
# By particle space (ps), we mean the guiding-center (E,p,R,z) four-dimensional space

# Script written by Henrik Järleblad. Last maintained 2022-10-11.
######################################################################################################

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
    include1Dto3D = false
    filepath_distr = "" # for example "c12_2_FI_distribution_800x50x43x49(+1).h5" or "c15_3_FI_distribution_101x52x45x47.jld2"
    filepath_equil = "" #JET/g94701/g94701_0-50.7932.eqdsk  #JET/g96100/g96100_0-53.0012.eqdsk
    filepath_W = ""
    filepath_OG = ""
    flip_F_EpRz_pitch = false
    folderpath_o = "../OWCF_results/template/" # Output folder path. Finish with '/'
    folderpath_OWCF = $folderpath_OWCF # Set path to OWCF folder to same as main process (hence the '$')
    numOsamples = 12_000_000 # The total number of monte-carlo samples for the fast-ion distribution
    nbatch = 100_000 # The algorithm will save the sampling progress every nbatch sample
    os_equidistant = true # False is not possible yet.
    timepoint = nothing # If unknown, just leave as nothing. The algorithm will try to figure it out automatically.
    tokamak = "" # If unknown, leave as ""
    TRANSP_id = "" # If unknown, leave as ""
    useWeightsFile = false # Set to true, if you want to load the orbit grid from the orbit weights file filepath_W
    useOrbGridFile = false # Set to true, if you want to load the orbit grid from the orbit grid file filepath_OG
    verbose = true
    visualizeProgress = false

    # \/ Only necessary if !useWeightsFile && !useOrbGridFile
    Emin = 0.0 # keV
    Emax = 0000.0 # keV
    FI_species = "D" # D,p,T,He3 etc
    inclPrideRockOrbs = true # If true, then pride rock orbits will be included. Otherwise, minimum(Rm) = magnetic axis.
    nE = 0
    npm = 0
    nRm = 0
    pm_min = -1.0
    pm_max = 1.0
    Rm_min = nothing # Automatically use magnetic axis or 4/5 from HFS wall to magnetic axis. Override by specifying together with Rm_max
    Rm_max = nothing # Automatically use LFS wall. Override by specifying together with Rm_min
    # /\ Only necessary if useWeightsFile set to false

    interp = true
    # \/ Only necessary if interp is set to true
    nE_ps = 201
    np_ps = 56
    nR_ps = 59
    nz_ps = 67
    # /\ Only necessary if interp is set to true

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
include("ps2os.jl")
## -----------------------------------------------------------------------------
# Then you clean up after yourself
if batch_job && distributed
    for i in workers()
        rmprocs(i)
    end
end
