################################ start_ps2os_template.jl #########################################
# This file contains all the inputs needed by the script ps2os.jl. Copy this file, and change
# inputs accordingly. After reading the inputs, the script will execute ps2os.jl. As of this version
# of the OWCF, the ps2os.jl script will use Monte-Carlo sampling to transform from (E,p,R,z) to 
# (E,pm,Rm). In future version, other methods of transforming are planned to be selectable.
#
### The inputs are:
#
# batch_job_SLURM - If true, the script assumes that it will be executed as part of an high-performance computational 
#                   cluster batch job within the SLURM cluster environment. The script will act accordingly. - Bool
# distributed - If true, parallel computing with multiple CPU cores will be used. - Bool
# folderpath_OWCF - The path to the OWCF folder - String
# numOcores - The number of CPU cores that will be used if distributed is set to true. - Int64
#
# include1Dto3D - If true, then the 3D (E,pm,Rm) quantity will be saved as well - Boolean
# filepath_EpRz - The path to the (E,p,R,z) quantity to be transformed to (E,pm,Rm). Must be .h5 or .jld2 file format - String
# filepath_equil - The path to the .eqdsk-file (or .geqdsk/.jld2-file) with the tokamak magnetic equilibrium and geometry - String
# filepath_W - The path to the weights file to use as blueprint for the orbit grid. If desired. - String
# filepath_OG - The path to the output file from calcOrbGrid.jl. If specified, ps2os.jl will use that as the (E,pm,Rm) grid points - String
# folderpath_o - Path to folder where you want your results. End with '/' (or '\' if Windows...) - String'
# h5_is_rowmajor - If the 'filepath_EpRz' file is an .h5/.hdf5 file and it was saved using a row-major programming language (e.g. C/C++ or 
#                  NumPy in Python, please see https://en.wikipedia.org/wiki/Row-_and_column-major_order for more info), this input variable 
#                  should be set to true - Bool
# numOsamples - The number of orbit samples you want. The more, the better the transform. But takes longer to compute! - Int64
# nbatch - The algorithm will save the progress every nbatch sample. Useful if the sampling gets cancelled unexpectedly - Int64
# os_equidistant - If true, script assumes equidistant orbit grid. False not supported yet - Boolean
# sign_o_pitch_wrt_B - If true, script assumes sign of pitch (v_||/v) of input quantity to be defined w.r.t the B-field, and not the plasma current (OWCF standard) - Boolean
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
# FI_species - The particle species of the (E,p,R,z) quantity. E.g. "D", "T", "3he" etc - String
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
# interp - If set to true, the algorithm will interpolate the loaded (E,p,R,z) quantity - Bool
# nE_ps - The number of fast-ion energy grid points in particle space to interpolate onto - Int64
# np_ps - The number of fast-ion pitch grid points in particle space to interpolate onto - Int64
# nR_ps - The number of fast-ion R grid points in particle space to interpolate onto - Int64
# nz_ps - The number of fast-ion z grid points in particle space to interpolate onto - Int64

### Other
# By particle space (ps), we mean the guiding-center (E,p,R,z) four-dimensional space

# Script written by Henrik Järleblad. Last maintained 2025-10-07.
######################################################################################################

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
@everywhere begin
    include1Dto3D = false
    filename_start = reduce(*,split(split(@__FILE__,"/")[end],".")[1:end-1]) # @__FILE__ returns the path and file name of the start file. Remove path and don't include the .jl file extension of the start file
    filepath_EpRz = "" # For example "c12_2_FI_distribution_800x50x43x49(+1).h5" or "c15_3_FI_distribution_101x52x45x47.jld2" or "EpRzTopoMap_JET_99965K73_at48,4058s_D_100x50x51x54_wLost_wJac.jld2"
    filepath_equil = "" #E.g. "equilibrium/JET/g94701/g94701_0-50.7932.eqdsk"  or "equilibrium/JET/g96100/g96100_0-53.0012.eqdsk"
    filepath_W = ""
    filepath_OG = ""
    folderpath_o = "../OWCF_results/template/" # Output folder path. Finish with '/'
    folderpath_OWCF = $folderpath_OWCF # Set path to OWCF folder to same as main process (hence the '$')
    h5_is_rowmajor = false # Set to true, if 'filepath_EpRz' is an .h5/.hdf5 and it was saved using a row-major programming language
    numOsamples = 12_000_000 # The total number of monte-carlo samples for the fast-ion distribution, if that is the (E,p,R,z) quantity
    nbatch = 100_000 # The algorithm will save the sampling progress every nbatch sample
    os_equidistant = true # False is not possible yet.
    sign_o_pitch_wrt_B = false # Set to true, if sign(dot(J,B))<0 and the sign of the pitch of the input quantity at filepath_EpRz is defined w.r.t. to the B-field (and not the plasma current)
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
include("ps2os.jl")
## -----------------------------------------------------------------------------
# Then you clean up after yourself
if batch_job_SLURM || distributed
    for i in workers()
        rmprocs(i)
    end
end