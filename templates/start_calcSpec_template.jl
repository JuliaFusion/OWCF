################################ start_calcSpec_template.jl #########################################
# This file contains all the inputs that the script calcSpec.jl needs to calculate an expected diagnostic
# signal for a specific diagnostic. It also executes the script calcSpec.jl after the inputs are defined.
# Please also read the description in calcSpec.jl for further input information.
# The inputs are as follows:
#
# batch_job - If true, the script assumes that it will be executed via an HPC batch job. The script will act accordingly. - Bool
# distributed - If true, parallel computing with multiple CPU cores will be used. - Bool
# numOcores - The number of CPU cores that will be used if distributed is set to true. - Int64
#
# addNoise - If set to true, the script will add noise to the synthetic signal, defined by input variables noiseLevel_b and _s. - Bool
# diagnostic_filepath - The filepath to the geometry-data for the diagnostic for which to calculate the synthetic signal. If not specified, spherical emission is assumed - String
# diagnostic_name - The name of the diagnostic for which to compute the synthetic signal. Purely for esthetic purposes- String
# Ed_min - The lower bound for the diagnostic energy grid (measurement grid) - Float64
# Ed_max - The upper bound for the diagnostic energy grid (measurement grid) - Float64
# Ed_diff - The width of the diagnostic energy bins (PLEASE NOTE! THIS WILL INDIRECTLY DEFINE THE NUMBER OF DIAGNOSTIC ENERGY GRID POINTS.) - Float64
# filepath_equil - The path to the .eqdsk-file (or .geqdsk/.jld2-file) with the tokamak magnetic equilibrium and geometry - String
# filepath_FI_distr - The path to the fast-ion distribution file to extract tokamak/fast-ion data from. Must be .cdf, .h5 or .jld2 file format - String
# filepath_FI_TRANSP_shot - The path to the TRANSP .cdf fast-ion file, if filepath_FI_distr or filepath_thermal_distr is specified as a .cdf TRANSP file. If so, should be same as filepath_FI_distr - String
# filepath_thermal_distr - The path to the thermal plasma distribution file. Must be .cdf (TRANSP) or .jld2 file format. Otherwise "" - String
# filepath_TRANSP_shot - The path to the TRANSP .cdf shot file, if filepath_FI_distr or filepath_thermal_distr is specified as a .cdf TRANSP file - String
# folderpath_o - The path to the folder where the results will be saved - String
# mc_samples - The number of Monte-Carlo samples to use for the signal computation (sample from the fast-ion distribution) - Int64
# mc_chunk - The number of Monte-Carlo samples in a single sample-partitioned chunk - Int64
# noiseLevel_b - The level of the background noise for the synthetic signal, if addNoise is set to true - Float64
# noiseLevel_s - The level of the signal noise for the synthetic signal, if addNoise is set to true - Float64
# folderpath_OWCF - The script needs to know where it is located. Yup, I know. A little bit awkward isn't it? - String
# reaction - The nuclear fusion reaction that you want to simulate. Please see OWCF/misc/availReacts.jl for available fusion reactions - String
# calcProjVel - If true, then the synthetic diagnostic spectrum will be the result of projected fast-ion velocities. Remember to set 'reaction' to 'proj-X' where 'X' is the fast-ion species which projected velocity you want to compute - Bool
# timepoint - The timepoint of the tokamak shot for the magnetic equilibrium. Format XX,YYYY where XX are seconds and YYYY are decimals - String
# thermal_temp_axis - The temperature of the thermal plasma distribution on axis, if filepath_thermal_distr is not specified - Float64
# thermal_dens_axis - The density of the thermal plasma distribution on axis, if filepath_thermal_distr is not specified - Float64
# verbose - If true, lots of information will be printed during execution. Fun! - Boolean
# visualizeProgress - If false, progress bar will not be displayed during computations - Boolean
#
# If you have specified a .h5 or .jld2 fast-ion distribution file, you can ask to have it interpolated
# so as to have the grid resolution of your choice. This can be good for investigating how grid resolution
# affects the resulting synthetic signal. The extra inputs are as follows:
# interp - If true, then the loaded fast-ion distribution will be interpolated - Boolean
# nE_ps - The energy grid resolution after interpolation - Int64
# np_ps - The pitch grid resolution after interpolation - Int64
# nR_ps - The R grid resolution after interpolation - Int64
# nz_ps - The z grid resolution after interpolation - Int64

# If you want, you can also manually specify a tokamak. Only for esthetic purposes. It is
# not necessary to specify and might be overwritten anyway. The input is as follows:
# tokamak - The tokamak in which the diagnostic and fast-ion distribution is located - String

### Other
# Please note, when specifying files (distribution files etc.) finish the path without '/'.
# When specifying folders (folderpath_o etc.) finish the path with '/'.
# In short, finish path to files without '/'. Finish path to folder with '/'.
#
# Also, please note, that if you specify a TRANSP .cdf file, regardless of it being thermal or fast-ion (FI),
# you must also specify the other TRANSP .cdf file with filepath_TRANSP_shot/filepath_FI_TRANSP_shot.
# Further info below at the variable specifications. If no TRANSP .cdf files have been specified both
# for filepath_FI_distr and filepath_thermal_distr, then filepath_TRANSP_shot and filepath_FI_TRANSP_shot
# can remain unspecified. Why this complicated scheme of file specification? Because if a TRANSP .cdf
# file has been specified, the algorithm needs to know the fast-ion time window to extract data from.
# That is only stored in the TRANSP .cdf fast-ion file

# Script written by Henrik Järleblad. Last maintained 2022-10-11.
######################################################################################################

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
# Then you load the script inputs
@everywhere begin
    addNoise = false # If true, then the algorithm will automatically add noiseLevel % noise to your synthetic diagnostic signal
    diagnostic_filepath = "" # Currently supported: "TOFOR", "AB" and ""
    diagnostic_name = ""
    Ed_min = 0000.0 # keV
    Ed_max = 0000.0 # keV
    Ed_diff = 00.0 # keV
    filepath_equil = "" # for example "g94701_0-50.7932.eqdsk" or "g96100_0-53.0012.eqdsk"
    filepath_FI_distr = "" # for example "96100J01_fi_1.cdf", "c12_2_FI_distribution_800x50x43x49(+1).h5" or "c15_3_FI_distribution_101x52x45x47.jld2"
    filepath_FI_TRANSP_shot = "" # As an example, if filepath_thermal_distr=="96100J01.cdf" then this variable should be "96100J01_fi_1.cdf". If filepath_FI_distr=="96100J01_fi_1.cdf", then filepath_FI_TRANSP_shot==filepath_FI_distr.
    filepath_thermal_distr = "" # for example "96100J01.cdf", "c21_3_thermal_profiles.jld2" or ""
    filepath_TRANSP_shot = "" # As an example, if filepath_FI_distr=="96100J01_fi_1.cdf" then this variable should be "96100J01.cdf". If filepath_thermal_distr=="96100J01.cdf", then filepath_TRANSP_shot==filepath_thermal_distr.
    folderpath_o = "../OWCF_results/template/" # Output folder path. Finish with '/'
    mc_samples = 1_000_000
    mc_chunk = 10_000 # Should be smaller than mc_samples. If unsure, do not alter from default value of 10_000
    noiseLevel_b = 0.05 # Background noise level. As a decimal fraction of 1.0. Could also be greater than 1.0 but why would you want that?
    noiseLevel_s = 0.05 # Signal noise level. As a decimal fraction of 1.0. Could also be greater than 1.0 but why would you want that?
    folderpath_OWCF = "" # The OWCF folder path. Finish with '/'
    reaction = "D(d,n)3He" # Specified on the form a(b,c)d where a is thermal ion, b is fast ion, c is emitted particle and d is the product nucleus.
    # PLEASE NOTE! Specify alpha particles as '4he' or '4He' (NOT 'he4' or 'He4'). Same goes for helium-3 (specify as '3he', NOT 'he3')
    calcProjVel = false # If true, the synthetic diagnostic spectrum will be computed from projected velocities. If true, remember to set the 'reaction' input variable to the format 'proj-X' where 'X' is the fast-ion species of interest. If true, remember to also leave the thermal plasma input variables unspecified.
    timepoint = nothing # If unknown, just leave as nothing. The algorithm will try to figure it out automatically.
    thermal_temp_axis = nothing # keV. Please specify if filepath_thermal_distr is not specified.
    thermal_dens_axis = nothing # m^-3. Please specify if filepath_thermal_distr is not specified.
    verbose = true # If set to true, the script will talk a lot! Yay!
    visualizeProgress = false # If set to true, a progress bar will be displayed during distributed computation

    # If you would like to have the fast-ion loaded from .h5/.jld2 file interpolated onto specific grid dimensions
    interp = true
    nE_ps = 201
    np_ps = 61
    nR_ps = 59
    nz_ps = 63

    tokamak = "" # Set to "" if unknown
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
include("calcSpec.jl")
## -----------------------------------------------------------------------------
# Then you clean up after yourself
if batch_job && distributed
    for i in workers()
        rmprocs(i)
    end
end
