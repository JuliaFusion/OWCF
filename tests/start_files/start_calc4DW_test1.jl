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
# filepath_thermal_distr - The path to the thermal distribution file to extract thermal species data from. Must be a valid TRANSP .cdf file, a valid .jld2 file or "".
#                          Please see the file description of the OWCF/calc4DWeights.jl script for more info - String
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
# plot_results - If true, the weight functions will be plotted upon completion of computation. The figure(s) will be saved as .png files in the 
#                folderpath_o folder - Bool
# R_array - The major radius grid points of your (E,p,R,z) grid. If set to 'nothing': nR, R_min and R_max must be specified - Vector
# R_min - The lower boundary for the (E,p,R,z) grid R values - Float64
# R_max - The upper boundary for the (E,p,R,z) grid R values - Float64
# reaction - The nuclear fusion reaction that you want to simulate. Please see OWCF/misc/availReacts.jl for available fusion reactions - String
# saveOnlyNonZeroWeights - If set to true, only the non-zero weight will be saved, as values and (E,p,R,z) indices - Bool
# saveVparaVperpWeights - If set to true, the weight functions will be saved on a (vpara,vperp) grid, in addition to (E,p) - Bool
# timepoint - The timepoint of the tokamak shot for the magnetic equilibrium. Format XX,YYYY where XX are seconds and YYYY are decimals - String
# thermal_temp_axis - The temperature of the thermal species distribution on axis, if filepath_thermal_distr is not specified - Float64
# thermal_temp_profile_type - Specifies the thermal ion temperature profile type, if filepath_thermal_distr is not specified. 
#                             Currently supported options include :DEFAULT and :FLAT. Please see the OWCF/misc/temp_n_dens.jl file 
#                             for more info on the :DEFAULT profile - Symbol
# thermal_dens_axis - The density of the thermal species distribution on axis, if filepath_thermal_distr is not specified - Float64
# thermal_dens_profile_type - Specifies the thermal ion density profile type, if filepath_thermal_distr is not specified. 
#                             Currently supported options include :DEFAULT and :FLAT. Please see the OWCF/misc/temp_n_dens.jl file 
#                             for more info on the :DEFAULT profile - Symbol
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

# Script written by Henrik Järleblad. Last maintained 2025-06-11.
######################################################################################################

## First you have to set the system specifications
using Distributed # Needed, even though distributed might be set to false. This is to export all inputs to all workers right away, if needed.
#batch_job = false
distributed = false

############---------------------------------------------------------------------------------------###
# We need to thoroughly deduce if the user wants the 'plot_test_results' input variable to be true or false

# First, check if the length of the Julia input arguments list is greater than 1
if length(ARGS)>1
    if "plot_test_results" in lowercase.(ARGS) # If the argument list contains the 'plot_test_results' input variable
        # Assume that the boolean value for the 'plot_test_results' input variable is provided as the input argument directly after the 'plot_test_results' input argument
        i_bool = findfirst(x-> x=="plot_test_results", lowercase.(ARGS))+1
        try 
            # Declare global scope. To be able to use the 'plot_test_results' input variable outside of this try-catch statement
            global plot_test_results = parse(Bool, ARGS[i_bool])
        catch
            # If anything goes wrong, assume that the 'plot_test_results' input variable should be set to false
            global plot_test_results = false
        end
    end
elseif @isdefined plot_test_results # If not, check if the 'plot_test_results' variable has already been defined
    plot_test_results = plot_test_results # Use that value (might have been set in a super script, with this script run via 'include("OWCF/start_files/start_..._test....jl")')
else # If nothing else, assume that the 'plot_test_results' input variable should be set to false
    plot_test_results = false
end
###------------------------------------------------------------------------------------------------###

############---------------------------------------------------------------------------------------###
# Define the folderpath_OWCF variable, if not already defined in a super script
if !(@isdefined folderpath_OWCF)
    folderpath_OWCF = reduce(*,map(x-> "/"*x,split(@__DIR__,"/")[2:end-2]))*"/" # We know that the test start file is located in the OWCF/tests/start_files/ folder. Deduce the full OWCF folder path from that information
end
# Create the OWCF/tests/outputs/ folder, if it does not already exist
if !isdir(folderpath_OWCF*"tests/outputs/")
    print("The folder $(folderpath_OWCF)tests/outputs/ does not exist. Creating... ")
    mkdir(folderpath_OWCF*"tests/outputs")
    println("ok!")
end
# Change the working directory to the OWCF/ folder, and activate the OWCF Julia environment
@everywhere begin
    using Pkg
    cd(folderpath_OWCF)
    Pkg.activate(".")
end
###------------------------------------------------------------------------------------------------###

#numOcores = 4 # When executing script via HPC cluster job, make sure you know how many cores you have requested for your batch job

## Navigate to the OWCF folder and activate the OWCF environment
#cd(folderpath_OWCF)
#using Pkg
#Pkg.activate(".")

## If running as a batch job on a SLURM CPU cluster
#if batch_job && distributed
#    # Load the SLURM CPU cores
#    using ClusterManagers
#    addprocs(SlurmManager(numOcores))
#    hosts = []
#    pids = []
#    for i in workers()
#        host, pid = fetch(@spawnat i (gethostname(), getpid()))
#        push!(hosts, host)
#        push!(pids, pid)
#    end
#    @show hosts
#end

## If running locally and multi-threaded
#if !batch_job && distributed # Assume you are executing the script on a local laptop (/computer)
#    println("Adding processes... ")
#    addprocs(numOcores-(nprocs()-1)) # If you didn't execute this script as an HPC cluster job, then you need to add processors like this. Add all remaining available cores.
#    # The '-(nprocs()-1)' part is simply to ensure to extra processes are added, in case script needs to be restarted on a local computer
#end

## -----------------------------------------------------------------------------
@everywhere begin
    debug = false
    diagnostic_filepath = folderpath_OWCF*"vc_data/TOFOR/TOFOR.vc" # In addition to a file path, the following Strings also supported: "TOFOR", "AB" and ""
    diagnostic_name = "calc4DW_test1"
    instrumental_response_filepath = "" # Should be the filepath to three .txt-files or one .jld2-file. Otherwise, leave as ""
    instrumental_response_output_units = "" # Should be specified as described in OWCF/misc/convert_units.jl. If instrumental_response_filepath=="", leave as ""
    Ed_min = 19600.0 # keV (or m/s if 'reaction' input variable is on form (3) (please see OWCF/misc/availReacts.jl for further explanation))
    Ed_max = 20200.0 # keV (or m/s if 'reaction' input variable is on form (3) (please see OWCF/misc/availReacts.jl for further explanation))
    Ed_diff = 10.0 # keV (or m/s if 'reaction' input variable is on form (3) (please see OWCF/misc/availReacts.jl for further explanation))
    E_array = nothing # keV. Array can be specified manually. Otherwise, leave as 'nothing'
    Emin = 20.0 # keV
    Emax = 500.0 # keV
    filename_o = "calc4DW_test1" # For example "myCustomName". Otherwise, leave as nothing for default output file naming by the OWCF
    filepath_equil = folderpath_OWCF*"equilibrium/JET/g94701/g94701_0-50.7932.eqdsk" # for example "equilibrium/JET/g96100/g96100_0-53.0012.eqdsk" or "myOwnSolovev.jld2"
    filepath_FI_cdf = "" # If filepath_thermal_distr=="96100J01.cdf", then filepath_FI_cdf should be "96100J01_fi_1.cdf" for example
    filepath_thermal_distr = "" # for example "96100J01.cdf", "myOwnThermalDistr.jld2" or ""
    folderpath_o = ($folderpath_OWCF)*"tests/outputs/" # Output folder path. Finish with '/'
    folderpath_OWCF = $folderpath_OWCF # Set path to OWCF folder to same as main process (hence the '$')
    gyro_samples = 50 # 50 is the default discretization number for the gyro-motion
    iiimax = 2 # The script will calculate iiimax number of weight matrices. They can then be examined in terms of similarity (to determine MC noise influence etc).
    (iiimax > 1) && (iii_average = true) # If true, the average of all weight matrices will be computed and saved. Without the "_i" suffix.
    nE = 11 # Please specify, if E_array is set to nothing
    np = 12 # Please specify, if p_array is set to nothing
    nR = 5 # Please specify, if R_array is set to nothing
    nz = 6 # Please specify, if z_array is set to nothing
    p_array = nothing # Array can be specified manually. Otherwise, leave as 'nothing'
    p_min = -1.0
    p_max = 1.0
    plot_results = $plot_test_results # If true, results will be plotted and figures saved as .png files
    R_array = nothing # Major radius grid points in meters
    R_min = :HFSwall # The major radius lower boundary, if R_array is set to nothing. Specify in meters e.g. "3.0", "3.4" etc. Can also be specified as the symbols :r_mag or :HFSwall, then the major radius coordinate of the magnetic axis or (the minimum of) the HFS wall will automatically be used, respectively
    R_max = :LFSwall # The major radius upper boundary, if R_array is set to nothing. Specify in meters e.g. "3.0", "3.4" etc. Can also be specified as the symbols :r_mag or :LFSwall, then the major radius coordinate of the magnetic axis or (the maximum of) the LFS wall will automatically be used, respectively
    saveVparaVperpWeights = true # Set to true, and the weight functions will be saved in (vpara,vperp), in addition to (E,p)
    ################################################################################
    # The 'reaction' input variable below should be specified using one of the following forms:
    # (1) "a(b,c)d" 
    # (2) "a(b,c)d-l" 
    # (3) "b" 
    # where a is thermal ion, b is fast ion, c is fusion product particle of interest, d is fusion product particle of disinterest and l is the nuclear energy state of c. 
    # l can be GS, 1L or 2L, corresponding to Ground State (GS), 1st excited energy level (1L) and 2nd excited energy level (2L).
    # For lists of available fusion reactions and particle species, please see OWCF/misc/availReacts.jl and OWCF/misc/species_func.jl. The reaction forms imply:
    # (1) The standard fusion reaction form. The nuclear energy level 'l' of the fusion product particle of interest 'c' is automatically assumed to be GS (if relevant).
    # (2) The advanced fusion reaction form. The nuclear energy level 'l' of the fusion product particle of interest 'c' is taken to be 'l'.
    # (3) The projected velocity reaction form. No fusion reaction is computed. Instead, the 4D weight functions are computed from the velocity vectors of the ion 'b', projected onto the diagnostic line-of-sight (LOS), using the (E,p,R,z) points that are inside the LOS 
    # PLEASE NOTE! Specify alpha particles as '4he' or '4He' (NOT 'he4' or 'He4'). Same goes for helium-3 (specify as '3he', NOT 'he3'). Etc.
    reaction = "T(h,g)4He-GS"
    ################################################################################
    # PLEASE NOTE! Specify alpha particles as '4he' or '4He' (NOT 'he4' or 'He4'). Same goes for helium-3 (specify as '3he', NOT 'he3')
    timepoint = nothing # If unknown, just leave as nothing. The algorithm will try to figure it out automatically.
    thermal_temp_axis = 5.0 # keV. Please specify this if filepath_thermal_distr, filepath_FI_cdf are not specified
    thermal_temp_profile_type = :DEFAULT # Currently supported values are :DEFAULT and :FLAT
    thermal_dens_axis = 1.0e20 # m^-3. Please specify this if filepath_thermal_distr and filepath_FI_cdf are not specified
    thermal_dens_profile_type = :DEFAULT # Currently supported values are :DEFAULT and :FLAT

    verbose = true # If true, then the program will be very talkative!
    visualizeProgress = false # If false, progress bar will not be displayed for computations

    # You could specify the tokamak as well, if you know it. Please note, it's only for esthetic purposes
    tokamak = "JET"

    z_array = nothing # Vertical grid points in meters
    z_min = :floor # The vertical lower boundary. Specify in meters e.g. "0.3", "0.4" etc. Can also be specified as the symbols :z_mag or :floor, then the vertical coordinate of the magnetic axis or the vertical minimum of the first wall will automatically be used, respectively
    z_max = :ceil # The vertical upper boundary. Specify in meters e.g. "0.3", "0.4" etc. Can also be specified as the symbols :z_mag or :ceil, then the vertical coordinate of the magnetic axis or the vertical minimum of the first wall will automatically be used, respectively
end

## -----------------------------------------------------------------------------
# Then you execute the script
include(folderpath_OWCF*"calc4DWeights.jl")

## -----------------------------------------------------------------------------
# Then you clean up after yourself
#if batch_job && distributed
#    for i in workers()
#        rmprocs(i)
#    end
#end
