################################ start_calcOG_template.jl #########################################
# This file contains all the inputs that the script calcOrbGrid.jl needs to calculate an orbit
# grid and the pertaining valid orbits. This file also executes the script calcOrbGrid.jl after the 
# inputs are defined.

# The inputs are as follows:
#
# batch_job - If true, the script assumes that it will be executed via an HPC batch job. The script will act accordingly. - Bool
# distributed - If true, parallel computing with multiple CPU cores will be used. - Bool
# folderpath_OWCF - The path to where the OWCF folder is saved on your computer - String
# numOcores - The number of CPU cores that will be used if distributed is set to true. - Int64
#
# filepath_equil - The path to the .eqdsk-file (or .geqdsk/.jld2-file) with the tokamak magnetic equilibrium and geometry - String
# folderpath_o - The path to the folder in which you want the output results - String
# FI_species - The ion species for which the orbit weights are computed. E.g. "D", "T" etc. See misc/species_func.jl for further info - String
# timepoint - The timepoint of the tokamak shot for the magnetic equilibrium. Format XX,YYYY where XX are seconds and YYYY are decimals - String
# tokamak - The tokamak of the magnetic equilibrium. "JET", "ITER" etc. If not applicable, set to "" - String
# TRANSP_id - The TRANSP id to identify tokamak/equilibrium data from - String
# verbose - If set to true, the script will talk a lot! - Bool
# visualzeProgress - If set to true, the progress of the computations will be clearly visualized via a progress bar - Bool
# E_array - The fast-ion energy (keV) grid points of your (E,pm,Rm) grid. If set to 'nothing': nE, Emin and Emax must be specified - Vector
# pm_array - The fast-ion pm grid points of your (E,pm,Rm) grid. If set to 'nothing': npm, pm_min and pm_max must be specified - Vector
# Rm_array - The fast-ion Rm grid points of your (E,pm,Rm) grid. If set to 'nothing': nRm, Rm_min and Rm_max must be specified - Vector
# nE - If E_array set to 'nothing', please specify the number of fast-ion energy grid points - Int64
# npm - If pm_array set to 'nothing', please specify the number of fast-ion pm grid points - Int64
# nRm - If Rm_array set to 'nothing', please specify the number of fast-ion Rm grid points - Int64
# Emin - If E_array set to 'nothing', please specify the lower fast-ion energy (keV) boundary - Float64
# Emax - If E_array set to 'nothing', please specify the upper fast-ion energy (keV) boundary - Float64
# pm_min - If pm_array set to 'nothing', please specify the lower fast-ion pm boundary - Float64
# pm_max - If pm_array set to 'nothing', please specify the upper fast-ion pm boundary - Float64
# Rm_min - If Rm_array set to 'nothing', please either specify lower fast-ion Rm (meters) boundary or leave as 'nothing' - Float64
# Rm_max - If Rm_array set to 'nothing', please either specify upper fast-ion Rm (meters) boundary or leave as 'nothing' - Float64
# inclPrideRockOrbs - If Rm_array, Rm_min and Rm_max are all set to 'nothing', the OWCF will automatically include the Pride-rock orbits with this boolean - Bool
# 
# extra_kw_args - Here you can define extra keyword arguments to be used when integrating the equations-of-motion
# for the guiding-center orbits. Please consult the OWCF manual for further info on this.

### Other
# 

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
    debug = false # If you want to run the calcOrbGrid.jl script in debug-mode, please set to true. Otherwise, it should always be set to false.
    filepath_equil = "" #for example "equilibrium/JET/g96100/g96100_0-53.0012.eqdsk" Do NOT finish with '/', it's a file.
    folderpath_o = "../OWCF_results/template/" # Output folder path. Finish with '/', it's a folder.
    folderpath_OWCF = $folderpath_OWCF # Set path to OWCF folder to same as main process (hence the '$')
    FI_species = "D" # D, T, p, 3he etc
    # PLEASE NOTE! Specify alpha particles as '4he' or '4He' (NOT 'he4' or 'He4'). Same goes for helium-3 (specify as '3he', NOT 'he3')
    timepoint = nothing # If unknown, just leave as nothing. The algorithm will try to figure it out automatically.
    tokamak = "" # You could specify the tokamak, if you know it. Please note, it's only for esthetic purposes
    TRANSP_id = "" # You could specify the TRANSP ID, if you know it. Please note, it's only for esthetic purposes
    verbose = true # If true, then the program will be very talkative!
    visualizeProgress = false # If false, progress bar will not be displayed for computations

    E_array = nothing # keV. Array can be specified manually. Otherwise, leave as 'nothing'
    pm_array = nothing # Array can be specified manually. Otherwise, leave as 'nothing'
    Rm_array = nothing # Meter. Array can be specified manually. Otherwise, leave as 'nothing'

    # If E_array and/or pm_array and/or Rm_array has been set to 'nothing', please specify these
    nE = 0
    npm = 0
    nRm = 0
    Emin = 0.0 # keV
    Emax = 000.0 # keV
    pm_min = -1.0
    pm_max = 1.0

    Rm_min = nothing # Automatically use magnetic axis or 4/5 from HFS wall to magnetic axis. Override by specifying together with Rm_max. Meter
    Rm_max = nothing # Automatically use LFS wall. Override by specifying together with Rm_min. Meter
    inclPrideRockOrbs = true # If true, then pride rock orbits will be included. Otherwise, minimum(Rm) = magnetic axis.

    # EXTRA KEYWORD ARGUMENTS BELOW (these will go into the orbit_grid() and get_orbit() functions from GuidingCenterOrbits.jl)
    extra_kw_args = Dict(:limit_phi => true, :max_tries => 0, :debug => debug)
    # limit_phi - Limits the number of toroidal turns for orbits
    # maxiter - The orbit integration algorithm will try progressively smaller timesteps these number of times
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
include("helper/calcOrbGrid.jl")

## -----------------------------------------------------------------------------
# Then you clean up after yourself
if batch_job && distributed
    for i in workers()
        rmprocs(i)
    end
end
