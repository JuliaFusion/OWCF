################################ start_calcOW_template.jl #########################################
# This file contains all the inputs that the script calcOrbWeights.jl needs to calculate orbit
# weight functions. This file also executes the script calcOrbWeights.jl after the inputs are defined.
# The computed orbit weight functions will be two-dimensional and have dimensions (channels, valid orbits)
# where 'channels' is ((Ed_max-Ed_min)/Ed_diff)-1 and 'valid orbits' is the number of valid orbits for
# the orbit grid for which the orbit weight functions have been computed.

# The inputs are as follows:
#
# batch_job - If true, the script assumes that it will be executed via an HPC batch job. The script will act accordingly. - Bool
# distributed - If true, parallel computing with multiple CPU cores will be used. - Bool
# folderpath_OWCF - The path to the OWCF folder - String
# numOcores - The number of CPU cores that will be used if distributed is set to true. - Int64
#
# analyticalOWs - If set to true, projected velocities will be used to compute the orbit weight functions. In that case, no thermal data is needed - Bool
# debug - If true, then the script will run in debug-mode. Should almost always be set to false - Bool
# diagnostic_filepath - The path to the LINE21 data diagnostic line-of-sight file. Leave as "" for assumed sperical emission - String
# diagnostic_name - The name of the diagnostic. Purely for esthetic purposes - String
# instrumental_response_filepath - The path to three .txt-files or one .jld2-file, containing the necessary data for representing 
#                                the instrumental response of the diagnostic. If paths to three .txt-files are specified, they should 
#                                be specified together in a vector of three strings. That is,
#                                
#                                instrumental_response_filepath = ["/path/to/response/matrix.txt","/path/to/particle/inputs.txt","/path/to/diagnostic/outputs.txt"]
#
#                                The first string is the filepath to the response matrix. The size of the matrix is 'ni x no', where 'ni' is the number of diagnostic input 
#                                grid points and 'no' is the number of output (actually being measured) grid points. The second string is the filepath to the input grid 
#                                points (vector). The third string is the filepath to the output grid points (vector). So, for example, for a proton recoil diagnostic, 
#                                the input should be incoming neutron energies (keV) and the output could be proton impact positions (e.g. cm). If, instead, the path to one 
#                                .jld2 file is specified, it should be specified as 
#
#                                instrumental_response_filepath = "/path/to/diagnostic/response/data.jld2"
#
#                                The keys of the .jld2-file should then be "response_matrix", "input" and "output". Input data should be a vector of length 'ni', output 
#                                data should be a vector of length 'no' and response_matrix should be a matrix of size 'ni x no'.
# instrumental_response_output_units - The units of measurement of the output grid points specified via the output in the 'instrumental_response_filepath' input variable.
#                                      Please see the OWCF/misc/convert_units.jl script for a list of all acceptable units in the OWCF - String
# Ed_min - The lower boundary for the diagnostic energy grid (measurement grid) - Float64
# Ed_max - The upper boundary for the diagnostic energy grid (measurement grid) - Float64
# Ed_diff - The width of the diagnostic energy bins (PLEASE NOTE! THIS WILL INDIRECTLY DEFINE THE NUMBER OF DIAGNOSTIC ENERGY GRID POINTS.) - Float64
# E_array - The fast-ion energy (keV) grid points of your (E,pm,Rm) grid. If set to 'nothing': nE, Emin and Emax must be specified - Vector
# Emin - The lower boundary for the fast-ion energy in orbit space - Float64
# Emax - The upper boundary for the fast-ion energy in orbit space - Float64
# filepath_equil - The path to the file with the tokamak magnetic equilibrium and geometry - String
# filepath_FI_cdf - Can be specified, if filepath_thermal_distr is a TRANSP .cdf shot file. See below for specifications - String
# filepath_thermal_distr - The path to the thermal distribution file to extract thermal species data from. Must be TRANSP .cdf, .jld2 file format or "" - String
# folderpath_o - The path to the folder where the results will be saved - String
# iiimax - If specified to be greater than 1, several copies of weight functions will be calculated. For comparison. - Int64
# inclPrideRockOrbs - If true, then pride-rock/pinch orbits will be included. Otherwise, minimum(Rm) = magnetic axis by default - Bool
# include2Dto4D - If true, then (in addition) the 2D orbit weight functions will be enflated to their 4D form (channels,nE,npm,nRm) and saved in the folderpath_o folder as well - Bool
# nE - The number of fast-ion energy grid points in orbit space - Int64
# npm - The number of fast-ion pm grid points in orbit space - Int64
# nRm - The number of fast-ion Rm grid points in orbit space - Int64
# og_filepath - The path to a .jld2 file containing the output of the calcOrbGrid.jl script. If specified, the values of E_array, pm_array, Rm_array, nE, npm, nRm, Emin, pm_min, Rm_min, Emax, pm_max and Rm_max won't matter. Leave as 'nothing' to use them instead - String
# pm_array - The fast-ion pm grid points of your (E,pm,Rm) grid. If set to 'nothing': npm, pm_min and pm_max must be specified - Vector
# pm_min - The lower boundary for the orbit-grid pm values - Float64
# pm_max - THe upper boundary for the orbit-grid pm values - Float64
# reaction - The nuclear fusion reaction that you want to simulate. Please see OWCF/misc/availReacts.jl for available fusion reactions - String
# Rm_array - The fast-ion Rm grid points of your (E,pm,Rm) grid. If set to 'nothing': nRm, Rm_min and Rm_max must be specified - Vector
# Rm_min - The lower boundary for the orbit-grid Rm values - Float64
# Rm_max - The upper boundary for the orbit-grid Rm values - Float64
# timepoint - The timepoint of the tokamak shot for the magnetic equilibrium. Format XX,YYYY where XX are seconds and YYYY are decimals - String
# thermal_temp_axis - The temperature of the thermal species distribution on axis, if filepath_thermal_distr is not specified - Float64
# thermal_dens_axis - The density of the thermal species distribution on axis, if filepath_thermal_distr is not specified - Float64
# verbose - If true, lots of information will be printed during execution - Bool
# visualizeProgress - If false, progress bar will not be displayed during computations - Bool
#
# extra_kw_args - Here you can define extra keyword arguments to be used when integrating the equations-of-motion
# for the guiding-center orbits. Please consult the OWCF manual for further info on this.

### Other
# If filepath_thermal_distr is not specified, then an interpolation object will be
# used as 'analytical' thermal temperature and thermal density profiles, respectively. The thermal_temp_axis
# and thermal_dens_axis variables will be used to scale the polynomial profiles to match the specified
# thermal temperature and thermal density at the magnetic axis. Please see the /misc/temp_n_dens.jl script for info.

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
    analyticalOWs = false # If true, then no thermal species data is needed. The weight functions will be computed solely from the projected velocity of the ion as it crosses the diagnostic sightline.
    debug = false
    diagnostic_filepath = "" # The filepath to the LINE21 file containing viewing cone data. Or, automatic sightlines include "TOFOR", "AB" and "". If diagnostic_angles wish to be used, set as "" or nothing
    diagnostic_name = ""
    instrumental_response_filepath = "" # Should be the filepaths to three .txt-files or one .jld2-file. Otherwise, leave as ""
    instrumental_response_output_units = "" # Should be specified as described in OWCF/misc/convert_units.jl. If instrumental_response_filepath=="", leave as ""
    Ed_min = 0000.0 # keV (or m/s if analyticalOWs===true)
    Ed_max = 0000.0 # keV (or m/s if analyticalOWs===true)
    Ed_diff = 00.0 # keV (or m/s if analyticalOWs===true)
    E_array = nothing # keV. Array can be specified manually. Otherwise, leave as 'nothing'
    Emin = 0.0 # keV
    Emax = 000.0 # keV
    filepath_equil = "" # for example "equilibrium/JET/g96100/g96100_0-53.0012.eqdsk" or "myOwnSolovev.jld2"
    filepath_FI_cdf = "" # If filepath_thermal_distr=="96100J01.cdf", then filepath_FI_cdf could be "96100J01_fi_1.cdf" for example. To auto-extract timepoint
    filepath_thermal_distr = "" # for example "96100J01.cdf", "myOwnThermalDistr.jld2" or ""
    folderpath_o = "../OWCF_results/template/" # Output folder path. Finish with '/'
    folderpath_OWCF = $folderpath_OWCF # Set path to OWCF folder to same as main process (hence the '$')
    iiimax = 1 # The script will calculate iiimax number of weight functions. They can then be examined in terms of similarity (to determine MC noise influence etc).
    inclPrideRockOrbs = true # If true, then pride rock orbits will be included. Otherwise, minimum(Rm) = magnetic axis.
    include2Dto4D = true # If true, then the 2D orbit weight functions will be enflated to their 4D form
    nE = 0
    npm = 0
    nRm = 0
    og_filepath = nothing # nothing, by default. Specify a string with the path to a .jld2 file computed with the calcOrbGrid.jl to use that orbit grid
    pm_array = nothing # Array can be specified manually. Otherwise, leave as 'nothing'
    pm_min = -1.0
    pm_max = 1.0
    product_state = "1L" #!# Can be GS, 1L or 2L. GS means Ground State. 1L and 2L are the first and second excited nucleus levels, respectively. 1L and 2L are only relevant if the 'reaction' input variable is set to a two-step gamma-ray reaction
    reaction = "9Be(4He,12C)n" # Specified in the form a(b,c)d where a is thermal ion, b is fast ion, c is emitted particle and d is the product nucleus. In the case of a 2-step fusion reaction with the goal of measuring gamma-rays, the c in a(b,c)d refers to the product particle that will de-excite and emit the gamma-ray. If analyticalOWs==true then 'reaction' should be provided in the format 'proj-X' where 'X' is the fast ion species ('D', 'T' etc)
    # PLEASE NOTE! Specify alpha particles as '4he' or '4He' (NOT 'he4' or 'He4'). Same goes for helium-3 (specify as '3he', NOT 'he3')
    Rm_array = nothing # Meter. Array can be specified manually. Otherwise, leave as 'nothing'
    Rm_min = nothing # Automatically use magnetic axis or 4/5 from HFS wall to magnetic axis. Override by specifying together with Rm_max. Meter
    Rm_max = nothing # Automatically use LFS wall. Override by specifying together with Rm_min. Meter
    timepoint = nothing # If unknown, just leave as nothing. The algorithm will try to figure it out automatically. Should be the absolute time. I.e. in JET, the 40 second start-up time should be included.
    thermal_temp_axis = 0.0 # keV. Please specify this if filepath_thermal_distr and filepath_FI_cdf are not specified
    thermal_dens_axis = 0.0e20 # m^-3. Please specify this if filepath_thermal_distr and filepath_FI_cdf are not specified

    verbose = true # If true, then the program will be very talkative!
    visualizeProgress = false # If false, progress bar will not be displayed for computations

    # EXTRA KEYWORD ARGUMENTS BELOW (these will go into the orbit_grid() and get_orbit() functions from GuidingCenterOrbits.jl)
    extra_kw_args = Dict(:limit_phi => true, :max_tries => 0)
    # limits the number of toroidal turns for orbits
    # The orbit integration algorithm will try progressively smaller timesteps these number of times

    # You could specify the tokamak as well, if you know it. Please note, it's only for esthetic purposes
    tokamak = "JET"
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
include("calcOrbWeights.jl")

## -----------------------------------------------------------------------------
# Then you clean up after yourself
if batch_job && distributed
    for i in workers()
        rmprocs(i)
    end
end
