################################ start_ps2WF_template.jl #########################################
# This file contains all the inputs that the script ps2WF.jl needs to calculate ultra-high
# resolution orbit weight function signals for a given diagnostic. If calcWFOs is set to true, the 
# ps2WF.jl script will also compute orbit splits of the diagnostic signal, fast-ion distribution,
# orbit weight functions and signal density. This can be visualized using e.g. signalWebApp.jl.
# This file also executes the script ps2WF.jl after the inputs are defined.

# The inputs are as follows:
#
# batch_job - If true, the script assumes that it will be executed via an HPC batch job. The script will act accordingly. - Bool
# distributed - If true, parallel computing with multiple CPU cores will be used. - Bool
# folderpath_OWCF - The path to the OWCF folder - String
# numOcores - The number of CPU cores that will be used if distributed is set to true. - Int64
#
# calcWFOs - If true, then the algorithm will compute the WF densities in orbit format. PLEASE NOTE! Takes extra RAM memory - Bool
# diagnostic_filepath - The path to the LINE21 data diagnostic line-of-sight file. Leave as "" for assumed sperical emission - String
# diagnostic_name - The name of the diagnostic. Purely for esthetic purposes - String
# Ed_min - The lower bound for the diagnostic energy grid (signal/measurement grid) - Float64
# Ed_max - The upper bound for the diagnostic energy grid (signal/measurement grid) - Float64
# Ed_diff - The width of the diagnostic energy bins (PLEASE NOTE! THIS WILL INDIRECTLY DEFINE THE NUMBER OF DIAGNOSTIC ENERGY GRID POINTS.) - Float64
# E_array - If specified, the ps2WF.jl will use the values in this array as the energy grid points in its (E,pm,Rm) grid. extrema(E_array) should be within extrema(E) where 'E' are the energy grid points of f(E,p,R,z) - Array{Float64,1}
# Emin - The lower bound for the fast-ion energy in orbit space, if E_array is left unspecified. Emin > minimum(E) should be respected, where 'E' are the energy grid points of f(E,p,R,z) - Float64
# Emax - The upper bound for the fast-ion energy in orbit space, if E_array is left unspecified. Emax < maximum(E) should be respected, where 'E' are the energy grid points of f(E,p,R,z) - Float64
# filepath_equil - The path to the .eqdsk-file (or .geqdsk/.jld2-file) with the tokamak magnetic equilibrium and geometry - String
# filepath_FI_distr - The path to the fast-ion distribution file to extract tokamak/fast-ion data from. Must be .h5 or .jld2 file format - String
# filepath_FI_TRANSP_shot - The path to the TRANSP .cdf fast-ion file, if filepath_thermal_distr is specified as a .cdf TRANSP file - String
# filepath_thermal_distr - The path to the thermal species distribution file. Must be .cdf (TRANSP) or .jld2 file format. Otherwise "" - String
# folderpath_o - The path to the folder where the results will be saved - String
# h5_is_rowmajor - If the 'filepath_FI_distr' file is an .h5/.hdf5 file and it was saved using a row-major programming language (e.g. C/C++ or 
#                  NumPy in Python, please see https://en.wikipedia.org/wiki/Row-_and_column-major_order for more info), this input variable 
#                  should be set to true - Bool
# inclPrideRockOrbs - If true, then pride-rock/pinch orbits will be included. Otherwise, minimum(Rm) = magnetic axis by default - Bool
# nE - The number of fast-ion energy grid points in orbit space, if E_array is left unspecified - Int64
# nEbatch - The batch size of your energy grid in orbit space. Must be <= nE - Int64
# npm - The number of fast-ion pm grid points in orbit space, if pm_array is left unspecified - Int64
# nRm - The number of fast-ion Rm grid points in orbit space, if Rm_array is left unspecified - Int64
# pm_array - If specified, the ps2WF.jl will use the values in this array as the pm grid points in its (E,pm,Rm) grid - Array{Float64,1}
# pm_min - The lower boundary for the orbit-grid pitch values - Float64
# pm_max - THe upper boundary for the orbit-grid pitch values - Float64
# reaction - The nuclear fusion reaction that you want to simulate. Please see OWCF/misc/availReacts.jl for available fusion reactions - String
# Rm_array - If specified, the ps2WF.jl will use the values in this array as the Rm grid points in its (E,pm,Rm) grid - Array{Float64,1}
# Rm_min - The lower boundary for the orbit-grid maximum major radius grid points. Only used if specified together with Rm_max - Float64
# Rm_max - The upper boundary for the orbit-grid maximum major radius grid points. Only used if specified together with Rm_min - Float64
# timepoint - The timepoint of the tokamak shot for the magnetic equilibrium. Format XX,YYYY where XX are seconds and YYYY are decimals - String
# thermal_dens_axis - The density of the thermal species distribution on axis, if filepath_thermal_distr is not specified - Float64
# thermal_temp_axis - The temperature of the thermal species distribution on axis, if filepath_thermal_distr is not specified - Float64
# visualizeProgress - If false, progress bar will not be displayed during computations - Boolean
# verbose - If true, lots of information will be printed during execution - Boolean
#
# You can also specify the tokamak if you want. It is purely for esthetic purposes, and might be overwritten anyway.
# tokamak - The tokamak for which the WF signal is computed - String
#
# extra_kw_args - Here you can define extra keyword arguments to be used when integrating the equations-of-motion
# for the guiding-center orbits. Please consult the OWCF manual for further info on this.

### Other
# Please note, when specifying files (distribution files etc.) finish the path without '/'.
# When specifying folders (folderpath_o etc.) finish the path with '/'.
#
# Also, please note, that if you specify a TRANSP .cdf file for the thermal distribution,
# you must also specify another TRANSP .cdf file with filepath_FI_TRANSP_shot.
# Further info below at the variable specifications. If no TRANSP .cdf file has been specified
# for filepath_thermal_distr, then filepath_FI_TRANSP_shot can remain unspecified.
# Why this complicated scheme of file specification? Because if a TRANSP .cdf
# file has been specified, the algorithm needs to know the fast-ion time window to extract data from.
# That is only stored in the TRANSP .cdf fast-ion file

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
    calcWFOs = false # Set to true if you want the orbit WF densities to be computed and saved as well. Will take extra computational time
    diagnostic_filepath = ""
    diagnostic_name = ""
    Ed_min = 0000.0 # keV
    Ed_max = 0000.0 # keV
    Ed_diff = 00.0 # keV
    E_array = nothing
    Emin = 0.0 # keV
    Emax = 0000.0 # keV
    filepath_equil = "" #g94701/g94701_0-50.7932.eqdsk" #g96100/g96100_0-53.0012.eqdsk"
    filepath_FI_distr = "" # for example "c12_2_FI_distribution_800x50x43x49(+1).h5" or "c15_3_FI_distribution_101x52x45x47.jld2"
    filepath_FI_TRANSP_shot = "" # As an example, if filepath_thermal_distr=="96100J01.cdf" then this variable should be "96100J01_fi_1.cdf".
    filepath_thermal_distr = "" # for example "96100J01.cdf", "c21_3_thermal_profiles.jld2" or ""
    folderpath_o = "../OWCF_results/template/" # Output folder path. Finish with '/'
    folderpath_OWCF = $folderpath_OWCF # Set path to OWCF folder to same as main process (hence the '$')
    h5_is_rowmajor = false # Set to true, if 'filepath_EpRz' is an .h5/.hdf5 and it was saved using a row-major programming language
    inclPrideRockOrbs = false # If true, then does what it says. Otherwise, magnetic axis will be used.
    nE = 0
    nEbatch = 2
    npm = 0
    nRm = 0
    og_filepath = nothing
    pm_array = nothing
    pm_min = -1.0
    pm_max = 1.0
    reaction = "D(d,n)3He" # Specified on the form a(b,c)d where a is thermal ion, b is fast ion, c is emitted particle and d is the product nucleus.
    # PLEASE NOTE! Specify alpha particles as '4he' or '4He' (NOT 'he4' or 'He4'). Same goes for helium-3 (specify as '3he', NOT 'he3')
    Rm_array = nothing
    Rm_min = nothing # Automatically use magnetic axis or 4/5 from HFS wall to magnetic axis. Override by specifying together with Rm_max
    Rm_max = nothing # Automatically use LFS wall. Override by specifying together with Rm_min
    timepoint = nothing # If unknown, just leave as nothing. The algorithm will try to figure it out automatically.
    thermal_dens_axis = nothing # m^-3. Please specify if filepath_thermal_distr is not specified.
    thermal_temp_axis = nothing # keV. Please specify if filepath_thermal_distr is not specified.
    visualizeProgress = false # If false, progress bar will not be displayed for computations
    verbose = true

    # You might also specify the tokamak for which the script is executed.
    # However, it is purely for esthetic reasons
    tokamak = "" # Might also be overwritten depending on the specified diagnostic

    # EXTRA KEYWORD ARGUMENTS BELOW (these will go into the orbit_grid() and get_orbit() functions from GuidingCenterOrbits.jl)
    extra_kw_args = Dict(:limit_phi => true, :max_tries => 0)
    # limits the number of toroidal turns for orbits
    # The orbit integration algorithm will try progressively smaller timesteps these number of times
    debugging = false # For de-bugging purposes
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
include("ps2WF.jl")

## -----------------------------------------------------------------------------
# Then you clean up after yourself
if batch_job && distributed
    for i in workers()
        rmprocs(i)
    end
end
