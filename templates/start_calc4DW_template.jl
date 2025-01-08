################################ start_calc4DW_template.jl #########################################
# This file contains all the inputs that the script calc4DWeights.jl needs to calculate 4D (E,p,R,z)
# weight functions. This file also executes the script calc4DWeights.jl after the inputs are defined.
# The computed 4D weight functions will be five-dimensional and have dimensions (nEd, nE, np, nR, nz)
# where nEd is the number of diagnostic measurement bins, nE is the number of energy grid points,
# np is the number of pitch grid points, nR is the number of major radius grid points and nz is the 
# number of vertical grid points.

# The inputs are as follows:
#
# batch_job - If true, the script assumes that it will be executed via an HPC batch job. The script will act accordingly. - Bool
# distributed - If true, parallel computing with multiple CPU cores will be used. - Bool
# folderpath_OWCF - The path to the OWCF folder - String
# numOcores - The number of CPU cores that will be used if distributed is set to true. - Int64
#
# analytical4DWs - If set to true, projected velocities will be used to compute the weight functions. In that case, no thermal data is needed - Bool
# debug - If true, then the script will run in debug-mode. Should almost always be set to false - Bool
# diagnostic_filepath - The path to the diagnostic line-of-sight file- (either from LINE21 or extra/createCustomLOS.jl). 
#                       Leave as "" for assumed sperical emission - String
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
# Ed_min - The lower boundary for the diagnostic measurement bins (measurement grid) - Float64
# Ed_max - The upper boundary for the diagnostic measurement bins (measurement grid) - Float64
# Ed_diff - The width of the diagnostic energy bins (PLEASE NOTE! THIS WILL INDIRECTLY DEFINE THE NUMBER OF DIAGNOSTIC ENERGY GRID POINTS.) - Float64
# E_array - The fast-ion energy (keV) grid points of your (E,p,R,z) grid. If set to 'nothing': nE, Emin and Emax must be specified - Vector
# Emin - The lower boundary for the fast-ion energy in orbit space - Float64
# Emax - The upper boundary for the fast-ion energy in orbit space - Float64
# filepath_equil - The path to the file with the tokamak magnetic equilibrium and geometry - String
# filepath_FI_cdf - To be specified, if filepath_thermal_distr is a TRANSP .cdf shot file. See below for specifications - String
# filepath_thermal_distr - The path to the thermal distribution file to extract thermal species data from. Must be TRANSP .cdf, .jld2 file format or "" - String
# folderpath_o - The path to the folder where the results will be saved - String
# gyro_samples - The number of points to discretize the gyro-motion for computation of synthetic spectra - Int64
# iiimax - If specified to be greater than 1, several copies of weight functions will be calculated. For comparison. - Int64
# iii_average - If set to true, an average of all the copies of weight functions will be computed. Saved without the "_i" suffix. - Bool
# nE - The number of fast-ion energy grid points in (E,p,R,z) space - Int64
# np - The number of fast-ion pitch grid points in (E,p,R,z) space - Int64
# nR - The number of fast-ion major radius grid points in (E,p,R,z) space - Int64
# nz - The number of fast-ion vertical grid points in (E,p,R,z) space - Int64
# p_array - The fast-ion pitch grid points of your (E,p,R,z) grid. If set to 'nothing': np, p_min and p_max must be specified - Vector
# p_min - The lower boundary for the (E,p,R,z) grid p values - Float64
# p_max - THe upper boundary for the (E,p,R,z) grid p values - Float64
# R_array - The major radius grid points of your (E,p,R,z) grid. If set to 'nothing': nR, R_min and R_max must be specified - Vector
# R_min - The lower boundary for the (E,p,R,z) grid R values - Float64
# R_max - The upper boundary for the (E,p,R,z) grid R values - Float64
# reaction - The nuclear fusion reaction that you want to simulate. Please see OWCF/misc/availReacts.jl for available fusion reactions - String
# saveOnlyNonZeroWeights - If set to true, only the non-zero weight will be saved, as values and (E,p,R,z) indices - Bool
# saveVparaVperpWeights - If set to true, the weight functions will be saved on a (vpara,vperp) grid, in addition to (E,p) - Bool
# timepoint - The timepoint of the tokamak shot for the magnetic equilibrium. Format XX,YYYY where XX are seconds and YYYY are decimals - String
# thermal_temp - The temperature of the thermal species distribution at the (R,z) point of interest - Float64
# thermal_temp_axis - The temperature of the thermal species distribution on axis, if filepath_thermal_distr is not specified - Float64
# thermal_dens - The density of the thermal species distribution at the (R,z) point of interest - Float64
# thermal_dens_axis - The density of the thermal species distribution on axis, if filepath_thermal_distr is not specified - Float64
# verbose - If true, lots of information will be printed during execution - Bool
# visualizeProgress - If false, progress bar will not be displayed during computations - Bool
# z_array - The vertical grid points of your (E,p,R,z) grid. If set to 'nothing': nz, z_min and z_max must be specified - Vector
# z_min - The lower boundary for the (E,p,R,z) grid z values - Float64
# z_max - The upper boundary for the (E,p,R,z) grid z values - Float64

### Other
# If filepath_thermal_distr is not specified, then an interpolation object will be
# used as 'analytical' thermal temperature and thermal density profiles, respectively. The thermal_temp_axis
# and thermal_dens_axis variables will be used to scale the polynomial profiles to match the specified
# thermal temperature and thermal density at the magnetic axis. Please see the /misc/temp_n_dens.jl script for info.

# Script written by Henrik Järleblad. Last maintained 2023-12-18.
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
    analytical4DWs = false # If true, then no thermal species data is needed. The weight functions will be computed solely from the projected velocity of the ion onto the diagnostic line-of-sight.
    debug = false
    diagnostic_filepath = "" # Currently supported: "TOFOR", "AB" and ""
    diagnostic_name = ""
    instrumental_response_filepath = "" # Should be the filepath to three .txt-files or one .jld2-file. Otherwise, leave as ""
    instrumental_response_output_units = "" # Should be specified as described in OWCF/misc/convert_units.jl. If instrumental_response_filepath=="", leave as ""
    Ed_min = 0000.0 # keV (or m/s if analyticalOWs===true)
    Ed_max = 0000.0 # keV (or m/s if analyticalOWs===true)
    Ed_diff = 00.0 # keV (or m/s if analyticalOWs===true)
    E_array = nothing # keV. Array can be specified manually. Otherwise, leave as 'nothing'
    Emin = 0.0 # keV
    Emax = 000.0 # keV
    filepath_equil = "" # for example "equilibrium/JET/g96100/g96100_0-53.0012.eqdsk" or "myOwnSolovev.jld2"
    filepath_FI_cdf = "" # If filepath_thermal_distr=="96100J01.cdf", then filepath_FI_cdf should be "96100J01_fi_1.cdf" for example
    filepath_thermal_distr = "" # for example "96100J01.cdf", "myOwnThermalDistr.jld2" or ""
    folderpath_o = "../OWCF_results/template/" # Output folder path. Finish with '/'
    folderpath_OWCF = $folderpath_OWCF # Set path to OWCF folder to same as main process (hence the '$')
    gyro_samples = 50 # 50 is the default discretization number for the gyro-motion
    iiimax = 1 # The script will calculate iiimax number of weight matrices. They can then be examined in terms of similarity (to determine MC noise influence etc).
    (iiimax > 1) && (iii_average = false) # If true, the average of all weight matrices will be computed and saved. Without the "_i" suffix.
    nE = 0 # Please specify, if E_array is set to nothing
    np = 0 # Please specify, if p_array is set to nothing
    nR = 0 # Please specify, if R_array is set to nothing
    nz = 0 # Please specify, if z_array is set to nothing
    p_array = nothing # Array can be specified manually. Otherwise, leave as 'nothing'
    p_min = -1.0
    p_max = 1.0
    R_array = nothing # Major radius grid points in meters
    R_min = :r_mag # The major radius lower boundary, if R_array is set to nothing. Specify in meters e.g. "3.0", "3.4" etc. Can also be specified as a symbol :r_mag, then the major radius coordinate of the magnetic axis will automatically be used
    R_max = :r_mag # The major radius upper boundary, if R_array is set to nothing. Specify in meters e.g. "3.0", "3.4" etc. Can also be specified as a symbol :r_mag, then the major radius coordinate of the magnetic axis will automatically be used
    saveVparaVperpWeights = false # Set to true, and the weight functions will be saved in (vpara,vperp), in addition to (E,p)
    reaction = "D(d,n)3He" # Specified on the form a(b,c)d where a is thermal ion, b is fast ion, c is emitted particle and d is the product nucleus. However, if analyticalOWs==true then 'reaction' should be provided in the format 'proj-X' where 'X' is the fast ion species ('D', 'T' etc)
    # PLEASE NOTE! Specify alpha particles as '4he' or '4He' (NOT 'he4' or 'He4'). Same goes for helium-3 (specify as '3he', NOT 'he3')
    timepoint = nothing # If unknown, just leave as nothing. The algorithm will try to figure it out automatically.
    thermal_temp = nothing # The thermal species temperature for the (R,z) point of interest. If filepath_thermal_distr is provided, leave as nothing
    thermal_temp_axis = 0.0 # keV. Please specify this if filepath_thermal_distr, filepath_FI_cdf and thermal_temp are not specified
    thermal_dens = nothing # The thermal species density for the (R,z) point of interest. If filepath_thermal_distr is provided, leave as nothing
    thermal_dens_axis = 0.0e20 # m^-3. Please specify this if filepath_thermal_distr and filepath_FI_cdf are not specified

    verbose = true # If true, then the program will be very talkative!
    visualizeProgress = false # If false, progress bar will not be displayed for computations

    # You could specify the tokamak as well, if you know it. Please note, it's only for esthetic purposes
    tokamak = "JET"

    z_array = nothing # Vertical grid points in meters
    z_min = :z_mag # The vertical lower boundary. Specify in meters e.g. "0.3", "0.4" etc. Can also be specified as a symbol :z_mag, then the vertical coordinate of the magnetic axis will automatically be used
    z_max = :z_mag # The vertical upper boundary. Specify in meters e.g. "0.3", "0.4" etc. Can also be specified as a symbol :z_mag, then the vertical coordinate of the magnetic axis will automatically be used
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
include("calc4DWeights.jl")

## -----------------------------------------------------------------------------
# Then you clean up after yourself
if batch_job && distributed
    for i in workers()
        rmprocs(i)
    end
end
