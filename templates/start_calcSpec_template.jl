################################ start_calcSpec_template.jl #########################################
# This file contains all the inputs that the script calcSpec.jl needs to calculate an expected diagnostic
# signal for a specific diagnostic. It also executes the script calcSpec.jl after the inputs are defined.
# Please also read the description in calcSpec.jl for further input information.
# The inputs are as follows:
#
# batch_job - If true, the script assumes that it will be executed via an HPC batch job. The script will act accordingly. - Bool
# distributed - If true, parallel computing with multiple CPU cores will be used. - Bool
# folderpath_OWCF - The script needs to know where the OWCF folder is located. Yup, I know. A little bit awkward isn't it? - String
# numOcores - The number of CPU cores that will be used if distributed is set to true. - Int64
#
# addNoise - If set to true, the script will add noise to the synthetic signal, defined by input variables noiseLevel_b and _s. - Bool
# diagnostic_filepath - The filepath to the geometry-data for the diagnostic for which to calculate the synthetic signal. If not specified, spherical emission is assumed - String
# diagnostic_name - The name of the diagnostic for which to compute the synthetic signal. Purely for esthetic purposes- String
# instrumental_response_filepath - The path to three .txt-files or one .jld2-file, containing the necessary data for representing 
#                                the instrumental response of the diagnostic. If paths to three .txt-files are specified, they should 
#                                be specified together in a vector of three strings. That is,
#                                
#                                instrumental_response_filepath = ["/path/to/reponse/matrix.txt","/path/to/particle/inputs.txt","/path/to/diagnostic/outputs.txt"]
#
#                                The first string is the filepath to the response matrix. The size of the matrix is ni x no, where ni is the number of diagnostic input 
#                                grid points and no is the number of output (actually being measured) grid points. The second string is the filepath to the input grid 
#                                points (vector). The third string is the filepath to the output grid points (vector). So, for example, for a proton recoil diagnostic, 
#                                the input could be incoming neutron energies and the output could be proton impact positions. If, instead, the path to one .jld2 file 
#                                is specified, it should be specified as 
#
#                                instrumental_response_filepath = "/path/to/diagnostic/response/data.jld2"
#
#                                The keys of the .jld2-file should then be "response_matrix", "input" and "output". Input data should be a vector of length ni, output 
#                                data should be a vector of length no and response_matrix should be a matrix of size ni x no.
# Ed_min - The lower bound for the diagnostic energy grid (measurement grid) - Float64
# Ed_max - The upper bound for the diagnostic energy grid (measurement grid) - Float64
# Ed_diff - The width of the diagnostic energy bins (PLEASE NOTE! THIS WILL INDIRECTLY DEFINE THE NUMBER OF DIAGNOSTIC ENERGY GRID POINTS.) - Float64
# filepath_equil - The path to the .eqdsk-file (or .geqdsk/.jld2-file) with the tokamak magnetic equilibrium and geometry - String
# filepath_FI_distr - The path to the fast-ion distribution file to extract tokamak/fast-ion data from. Must be .cdf, .h5 or .jld2 file format - String
# filepath_FI_TRANSP_shot - The path to the TRANSP .cdf fast-ion file, if filepath_FI_distr or filepath_thermal_distr is specified as a .cdf TRANSP file. If so, should be same as filepath_FI_distr - String
# filepath_thermal_distr - The path to the thermal species distribution file. Must be .cdf (TRANSP) or .jld2 file format. Otherwise "" - String
# filepath_TRANSP_shot - The path to the TRANSP .cdf shot file, if filepath_FI_distr or filepath_thermal_distr is specified as a .cdf TRANSP file - String
# folderpath_o - The path to the folder where the results will be saved - String
# mc_samples - The number of Monte-Carlo samples to use for the signal computation (sample from the fast-ion distribution) - Int64
# mc_chunk - The number of Monte-Carlo samples in a single sample-partitioned chunk - Int64
# noiseLevel_b - The level of the background noise for the synthetic signal, if addNoise is set to true - Float64
# noiseLevel_s - The level of the signal noise for the synthetic signal, if addNoise is set to true - Float64
# phase_space_point_of_interest - Specific E-, p-, R- and/or z-values can specified. This is done by specifying a Float64 value in place of the corresponding 'nothing' vector element (see below).
#                                 To illustrate the function of this input variable, consider the following. If only an E value has been specified, the expected diagnostic signal for f(p,R,z) at E 
#                                 will be computed. If only a p value has been specified, the expected diagnostic signal for f(E,R,z) at p will be computed. And so on.
#                                 If both an E and a p value (only) have been specified, the expected diagnostic signal for f(R,z) at (E,p) will be computed.
#                                 If both an R and a z value (only) have been specified, the expected diagnostic signal for f(E,p) at (R,z) will be computed. And so on.
#                                 If an E-, a p- and an R-value (three coordinates) have been specified, the expected diagnstic signal for f(z) at (E,p,R) will be computed.
#                                 And so on. PLEASE NOTE! This feature ONLY works with .h5 or .jld2 file inputs for the filepath_FI_distr input variable.
#                                 Furthermore, PLEASE NOTE that the R and/or z grid values in the .h5 or .jld2 file might be given in centimeters or meters.
#                                 So please specify the R- and/or z-values accordingly. E- values should always be given in keV. - Vector{Union{Nothing,Float64}}
# reaction - The nuclear fusion reaction that you want to simulate. Please see OWCF/misc/availReacts.jl for available fusion reactions - String
# calcProjVel - If true, then the synthetic diagnostic spectrum will be the result of projected fast-ion velocities. Remember to set 'reaction' to 'proj-X' where 'X' is the fast-ion species which projected velocity you want to compute - Bool
# timepoint - The timepoint of the tokamak shot for the magnetic equilibrium. Format XX,YYYY where XX are seconds and YYYY are decimals - String
# thermal_temp_axis - The temperature of the thermal species distribution on axis, if filepath_thermal_distr is not specified - Float64
# thermal_dens_axis - The density of the thermal species distribution on axis, if filepath_thermal_distr is not specified - Float64
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
# Then you load the script inputs
@everywhere begin
    addNoise = false # If true, then the algorithm will automatically add noiseLevel % noise to your synthetic diagnostic signal
    diagnostic_filepath = "" # Currently supported: "TOFOR", "AB" and ""
    diagnostic_name = ""
    instrumental_response_filepath = "" # Should be the filepaths to three .txt-files or one .jld2-file. Otherwise, leave as ""
    Ed_min = 0000.0 # keV
    Ed_max = 0000.0 # keV
    Ed_diff = 00.0 # keV
    filepath_equil = "" # for example "g94701_0-50.7932.eqdsk" or "g96100_0-53.0012.eqdsk"
    filepath_FI_distr = "" # for example "96100J01_fi_1.cdf", "c12_2_FI_distribution_800x50x43x49(+1).h5" or "c15_3_FI_distribution_101x52x45x47.jld2"
    filepath_FI_TRANSP_shot = "" # As an example, if filepath_thermal_distr=="96100J01.cdf" then this variable should be "96100J01_fi_1.cdf". If filepath_FI_distr=="96100J01_fi_1.cdf", then filepath_FI_TRANSP_shot==filepath_FI_distr.
    filepath_thermal_distr = "" # for example "96100J01.cdf", "my_custom_thermal_profiles.jld2" or ""
    filepath_TRANSP_shot = "" # As an example, if filepath_FI_distr=="96100J01_fi_1.cdf" then this variable should be "96100J01.cdf". If filepath_thermal_distr=="96100J01.cdf", then filepath_TRANSP_shot==filepath_thermal_distr.
    folderpath_o = "../OWCF_results/template/" # Output folder path. Finish with '/'
    folderpath_OWCF = $folderpath_OWCF # Set path to OWCF folder to same as main process (hence the '$')
    include_flr_effects = false
    mc_samples = 1_000_000
    mc_chunk = 10_000 # Should be smaller than mc_samples. If unsure, do not alter from default value of 10_000
    noiseLevel_b = 0.05 # Background noise level. As a decimal fraction of 1.0. Could also be greater than 1.0 but why would you want that?
    noiseLevel_s = 0.05 # Signal noise level. As a decimal fraction of 1.0. Could also be greater than 1.0 but why would you want that?
    phase_space_point_of_interest = [nothing,nothing,nothing,nothing] # To specify a specific E-, p-, R- and/or z-value, change the corresponding 'nothing' element to a Float64 value. Check the .h5 or .jld2 fast-ion file to decide on meters or cm for the R- and/or z-values [keV, -, meters or centimeters, meters or centimeters]
    reaction = "D(d,n)3He" # Specified on the form a(b,c)d where a is thermal ion, b is fast ion, c is emitted particle and d is the product nucleus.
    # PLEASE NOTE! Specify alpha particles as '4he' or '4He' (NOT 'he4' or 'He4'). Same goes for helium-3 (specify as '3he', NOT 'he3')
    calcProjVel = false # If true, the synthetic diagnostic spectrum will be computed from projected velocities. If true, remember to set the 'reaction' input variable to the format 'proj-X' where 'X' is the fast-ion species of interest. If true, remember to also leave the thermal species input variables unspecified.
    timepoint = nothing # If unknown, just leave as nothing. The algorithm will try to figure it out automatically.
    thermal_temp_axis = nothing # keV. Please specify if filepath_thermal_distr is not specified.
    thermal_dens_axis = nothing # m^-3. Please specify if filepath_thermal_distr is not specified.
    verbose = true # If set to true, the script will talk a lot! Yay!
    visualizeProgress = false # If set to true, a progress bar will be displayed during distributed computation

    # If you would like to have the fast-ion loaded from .h5/.jld2 file interpolated onto specific grid dimensions
    # If any specific E-, p-, R- or z-values have been specified in the input variable phase_space_point_of_interest,
    # the interpolation will be performed before using a specific f(p,R,z), f(E,p), f(z) etc
    interp = true
    nE_ps = 201
    np_ps = 61
    nR_ps = 59
    nz_ps = 63

    tokamak = "" # Set to "" if unknown
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
include("calcSpec.jl")
## -----------------------------------------------------------------------------
# Then you clean up after yourself
if batch_job && distributed
    for i in workers()
        rmprocs(i)
    end
end
