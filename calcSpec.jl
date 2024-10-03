################################ calcSpec.jl #########################################

#### Description:
# This script calculates the expected, synthetic diagnostic signal (S) from a
# specific fast-ion distribution and thermal species distribution.

### The fast-ion distribution can be provided via a TRANSP (fast-ion) .cdf file,
### via a .h5 file or via a .jld2 file.
## If you choose to provide a TRANSP fast-ion .cdf file (for example '94701V01_fi_2.cdf')
# you must also make sure that you have the matching TRANSP shot .cdf file (following
# the example, '94701V01.cdf' matches '94701V01_fi_2.cdf'). This is to be able to
# process the fast-ion information together with the TRANSP shot data in the Python framework.
#
## The .h5 file must have the keys:
# 'f' - The 4D matrix containing the fast-ion distribution.
# 'energy' - The 1D array containing the energy grid points
# 'pitch ' - The 1D array containing the pitch grid points
# 'R' - The 1D array containing the R grid points
# 'Z' or 'z' - The 1D array containing the z grid points
#
## The .jld2 file must have the keys:
# 'F_ps', 'f' or 'F_EpRz' - The 4D matrix containing the fast-ion distribution.
# 'energy' - The 1D array containing the energy grid points
# 'pitch ' - The 1D array containing the pitch grid points
# 'R' - The 1D array containing the R grid points
# 'Z' or 'z' - The 1D array containing the z grid points
#
### For the thermal species distribution, you have three options.
## Firstly, you can specify a TRANSP shot .cdf file (for example '94701V01.cdf'). The thermal
# species temperature and density profiles will then be extracted from the file.
#
## Secondly, you could instead specify a .jld2 file that contains information about the
# thermal species temperature and density profiles as a function of normalized flux coordinate
# ρ_pol. The .jld2 file must have the keys:
# 'thermal_temp' - The 1D array containing the thermal temperature for the ρ_pol grid points
# 'thermal_dens' - The 1D array containing the thermal density for the ρ_pol grid points
# 'rho_pol' - The 1D array containing the ρ_pol grid points
#
## Finally, you could also choose not to specify a thermal species distribution. The thermal
# species temperature and density profiles will then be set to default profiles with specified
# thermal_temp_axis and thermal_dens_axis temperature and density on-axis values, respectively.
#
# Please see the start_calcSpec_template.jl file for further input information.

#### Inputs (Units given when defined in script)
# Given via input file start_calcSpec_template.jl, for example. Copy the template file and replace
# 'template' by  a run-id. (you decide upon a run-id yourself. To be able to keep track of your work).

#### Outputs
# -

#### Saved files
# spec_[tokamak]_[TRANSP_id]_at[timepoint]s_[diagnostic_name]_[reaction].jld2
# The saved file will have the fields:
#   S - The computed expected (clean) diagnostic signal - Array{Float64,1}
#   Ed_array - The diagnostic measurement grid pertaining to S - Array{Float64,1}
#   nfast - The total number of fast ions in the fast-ion distribution function - Union{Float64,Int64}
# If addNoise input was set to true, there will also be the keys
#   S_noise - The computed expected (noisy) diagnostic signal - Array{Float64,1}
#   noise - The signal noise level for every measurement. length(S_noise)==length(S)==length(noise) - Array{Float64,1}
# If the input calcProjVel was set to true (i.e. the measurements are binned projected velocities), there will also be the key
#   calcProjVel - True, if calcProjVel was set to true in the calcSpec.jl inputs - Bool

### Other
# Please note that the diagnostic energy grid will be created as bin centers.
# That is, the first diagnostic energy grid value will be (Ed_min+Ed_diff/2) and so on.
#
# Please note that if you input a .h5 or a .jld2 fast-ion distribution file, the same sampling
# method will be used. But (!), if you input a TRANSP .cdf fast-ion distribution file, a different
# sampling method will be used. In theory this should not matter. In practice, it's not always so.
# The .h5/.jld2 Monte-Carlo sampling is done with Julia. The .cdf Monte-Carlo sampling is done
# with Python.
#
# Finally, this script has many 'if' statements and I would expect it to be able to be
# optimized beyond its current state. This will be done in future version of the OWCF.

# Script written by Henrik Järleblad. Last maintained 2023-11-07.
################################################################################################

## ---------------------------------------------------------------------------------------------
# Determine the file extensions for the fast-ion and thermal species distribution files
# This is the very first thing that needs to happen, because it determines which packages and dependencies
# files will be loaded.
fileext_FI = (split(filepath_FI_distr,"."))[end] # Assume last part after final '.' is the file extension
fileext_FI = lowercase(fileext_FI)
fileext_thermal = (split(filepath_thermal_distr,"."))[end] # Assume last part after final '.' is the file extension
fileext_thermal = lowercase(fileext_thermal)

if !(fileext_FI=="h5" || fileext_FI=="jld2" || fileext_FI=="cdf")
    println("Fast-ion distribution file format: ."*fileext_FI)
    error("Unknown fast-ion distribution file format. Please re-specify file and re-try.")
end

if !(fileext_thermal=="cdf" || fileext_thermal=="jld2" || fileext_thermal=="")
    println("Thermal distribution file format: ."*fileext_thermal)
    error("Unknown thermal distribution file format. Please re-specify file and re-try.")
end

if fileext_FI=="cdf"
    filepath_FI_TRANSP_shot = filepath_FI_distr # Safety measure. Assume TRANSP .cdf fast-ion file has been specified via an .cdf file
end

if fileext_thermal=="cdf"
    filepath_TRANSP_shot = filepath_thermal_distr # Safety measure. Assume TRANSP .cdf file has been specified via an .cdf file
end

if fileext_FI=="cdf" || fileext_thermal=="cdf"
    fileext_TRANSP_shot = (split(filepath_TRANSP_shot,"."))[end] # Assume last part after final '.' is the file extension
    fileext_TRANSP_shot = lowercase(fileext_TRANSP_shot)
    fileext_FI_TRANSP_shot = (split(filepath_FI_TRANSP_shot,"."))[end] # Assume last part after final '.' is the file extension
    fileext_FI_TRANSP_shot = lowercase(fileext_FI_TRANSP_shot)
end

if fileext_FI=="cdf" && !(fileext_TRANSP_shot=="cdf")
    error("filepath_FI_distr was specified as a TRANSP .cdf fast-ion file. But no corresponding TRANSP shot file was specified for filepath_TRANSP_shot. Please correct and re-try.")
end

if fileext_thermal=="cdf" && !(fileext_FI_TRANSP_shot=="cdf")
    error("filepath_thermal_distr was specified as a TRANSP .cdf file. But no corresponding TRANSP .cdf fast-ion file was specified for filepath_FI_TRANSP_shot. Please correct and re-try.")
end
## ---------------------------------------------------------------------------------------------
# Ensure that noise level is above 0.0
if addNoise
    if (noiseLevel_b < 0.0) || (noiseLevel_s < 0.0)
        error("Specified noise level is negative. Please correct and re-try.")
    end
end
## ---------------------------------------------------------------------------------------------
verbose && println("Loading Julia packages... ")
@everywhere begin
    cd(folderpath_OWCF) # Necessary to move all the workers to the correct folder
    using Distributions # For the add_noise() function
    using PyCall # For using Python code in Julia
    using EFIT # For calculating magn½etic equilibrium quantities
    using Equilibrium # For loading flux function data, tokamak geometry data etc.
    using JLD2 # To write/open .jld2-files (Julia files, basically)
    using FileIO # To write/open files in general
    using ProgressMeter # To visualize parallel computational process
    using NetCDF # To write/open .cdf files
    include("misc/diag2tokamak.jl") # To deduce tokamak machines from diagnostic's specifications
    include("misc/availReacts.jl") # To examine fusion reaction and extract thermal and fast-ion species
    include("misc/rewriteReacts.jl") # To rewrite a fusion reaction from the A(b,c)D format to the A-b=c-D format
    pushfirst!(PyVector(pyimport("sys")."path"), "") # To add the forward, transp_dists, transp_output and vcone modules (scripts in current path)
end
## ---------------------------------------------------------------------------------------------
if addNoise
    verbose && println("Defining the add_noise() function... ")
    """
        add_noise(S,0.05)
        add_noise(-||-, k=0.15)

    This function is also defined in extra/dependencies.jl. It is only defined here to avoid having to load the whole of dependencies.jl.
    """
    function add_noise(s, b; k=0.1)
        sn = max.(s,0.0) .+ k.*(mean(sqrt.(abs.(s)))).*rand.(Normal.(0.0, max.(sqrt.(max.(s,0)), sqrt.(b))))
        err = k.*mean(sqrt.(abs.(s))).*max.(sqrt.(max.(s,0)), sqrt.(b))
        return sn, err
    end
end
## ---------------------------------------------------------------------------------------------
verbose && println("Checking fusion reaction... ")
reaction_full = deepcopy(reaction) # Make a fully independent copy of the fusion reaction variable
#@everywhere reaction_full = $reaction_full # Not yet necessary. Might be necessary when two-step fusion reactions can be computed with the OWCF via the DRESS code
reaction = full2reactsOnly(reaction; projVelocity=calcProjVel) # Converts from 'a(b,c)d' format to 'a-b' format (reactants only)
@everywhere reaction = $reaction # Transfer to all external processes
emittedParticleHasCharge = false # By default, assume that the emitted particle 'c' in a(b,c)d does NOT have charge (is neutral)
RHEPWC = ["D-3He", "3He-D"] # RHEPWC means 'reaction has emitted particle with charge'
if reaction in RHEPWC # However, there are some fusion reactions which WILL produce an emitted particle with non-zero charge
    emittedParticleHasCharge = true
end

if emittedParticleHasCharge
    verbose && println("")
    verbose && println("The emitted "*getEmittedParticle(reaction_full)*" particle of the "*reaction_full*" reaction has non-zero charge!")
    verbose && println("For emitted particles with non-zero charge, the OWCF currently only supports computing the expected energy spectrum from the plasma as a whole (4*pi emission).")
    verbose && println("Therefore, the 'diagnostic_name' and 'diagnostic_filepath' input variables will be forcibly set to ''.")
    verbose && println("")
    diagnostic_name = ""
    diagnostic_filepath = ""
end

## ---------------------------------------------------------------------------------------------
# Safety check for projected velocity spectrum computations
if calcProjVel
    !(split(reaction,"-")[1] == "proj") && error("Projected velocity diagnostic spectrum computation was specified, but input variable 'reaction' was not correctly specified (it should be specified as 'proj-X' where 'X' is the fast-ion species). Please correct and re-try.")
end
if !calcProjVel
    (split(reaction,"-")[1] == "proj") && error("Normal synthetic diagnostic spectrum computation was specified, but input variable 'reaction' was not correctly specified (it should be specified as 'a(b,c)d' where a is thermal ion, b is fast ion, c is emitted particle and d is the product nucleus). Please correct and re-try.")
end
## ---------------------------------------------------------------------------------------------
verbose && println("Determining sampling variables... ")
if !(typeof(mc_samples)==Int)
    try 
        global mc_samples = Int(mc_samples)
    catch
        error("Input variable 'mc_samples' not specified as an integer. Please correct and re-try.")
    end
end
if !(typeof(mc_chunk)==Int)
    try
        global mc_chunk = Int(mc_chunk)
    catch
        error("Input variable 'mc_chunk' not specified as an integer. Please correct and re-try.")
    end
end
if mc_chunk > 10_000
    @warn "mc_chunk set to greater than 10 000. Depending on your memory settings, the computation might run out of RAM."
end
if mc_samples <= mc_chunk
    @warn "Warning! mc_samples <= mc_chunk. mc_samples will not be partitioned into smaller chunks."
    mc_sample_chunks = [mc_samples]
else
    mc_sample_chunks = zeros(Int(ceil(mc_samples/mc_chunk))) # Create an array to accomodate all the sample chunks
    mc_sample_chunks[1:end-1] .= mc_chunk # Almost all the chunks will be of size mc_chunk...
    mc_sample_chunks[end] = round(((mc_samples/mc_chunk)-floor(mc_samples/mc_chunk))*mc_chunk) #... Except for the last chunk, which will be the remainder
    mc_sample_chunks[end] = mc_sample_chunks[end] == 0.0 ? mc_chunk : mc_sample_chunks[end]
    mc_sample_chunks = Int.(mc_sample_chunks) # Make sure all the elements of the chunk array are integers
    if !(sum(mc_sample_chunks)==mc_samples)
        error("Something went wrong! Please try different values for mc_samples and mc_chunk, and re-try.")
    end
end
## ---------------------------------------------------------------------------------------------
verbose && println("If necessary, loading dependencies... ")

if fileext_FI=="h5" || fileext_FI=="jld2" || isfile(instrumental_response_filepath) || (typeof(instrumental_response_filepath)==Vector{String})
    verbose && println("OWCF/extra/dependencies.jl is needed. Loading (might take a while)... ")
    @everywhere begin
        include("extra/dependencies.jl") # Need the dependencies for sampling processes
    end
end

## ---------------------------------------------------------------------------------------------
verbose && println("Loading thermal .jld2 file data or TRANSP info... ")

if fileext_thermal=="jld2"
    myfile = jldopen(filepath_thermal_distr,false,false,false,IOStream)
    thermal_temp_array = myfile["thermal_temp"]
    @everywhere thermal_temp_array = $thermal_temp_array # These are sent to external processes here, for efficiency
    thermal_dens_array = myfile["thermal_dens"]
    @everywhere thermal_dens_array = $thermal_dens_array # These are sent to external processes here, for efficiency
    ρ_pol_array = myfile["rho_pol"]
    @everywhere ρ_pol_array = $ρ_pol_array # These are sent to external processes here, for efficiency
    close(myfile)

    @everywhere begin
        verbose && println("Loading the Interpolations package (needed because .jld2 thermal file was specified)... ")
        using Interpolations
    end
end

# Loading tokamak information and TRANSP RUN-ID from thermal .cdf file
if fileext_thermal=="cdf"
    verbose && println("Loading tokamak and TRANSP information... ")
    TRANSP_id = (split((split(filepath_thermal_distr,"."))[1],"/"))[end] # Assume part between last '/' and '.' is the TRANSP id
else
    TRANSP_id = ""
end

if tokamak==""
    verbose && println("Trying to deduce tokamak from diagnostic_name... ")
    if !(diagnostic_name=="")
        tokamak = diag2tokamak(diagnostic_name)
    end
end

if filepath_thermal_distr=="" && (!(typeof(thermal_temp_axis)==Float64) && !(typeof(thermal_temp_axis)==Int64))
    @everywhere thermal_temp_axis = 3.0
    @warn "filepath_thermal_distr was not specified, and thermal_temp_axis was not specified correctly. thermal_temp_axis will be set to default value of 3.0 keV."
end

if filepath_thermal_distr=="" && !(typeof(thermal_dens_axis)==Float64)
    @everywhere thermal_dens_axis = 1.0e19
    @warn "filepath_thermal_distr was not specified, and thermal_dens_axis was not specified correctly. thermal_dens_axis will be set to default value of 1.0e19 m^-3."
end

## ---------------------------------------------------------------------------------------------
# Define all variables so far only defined on the main process, also on the external processes
# Send the variables from the main process to all remote workers. The '$' indicates the main process
@everywhere fileext_thermal = $fileext_thermal
@everywhere fileext_FI = $fileext_FI
@everywhere TRANSP_id = $TRANSP_id
@everywhere tokamak = $tokamak
@everywhere mc_sample_chunks = $mc_sample_chunks
@everywhere calcProjVel = $calcProjVel

## ---------------------------------------------------------------------------------------------
# Loading tokamak equilibrium
verbose && println("Loading tokamak equilibrium... ")
if ((split(filepath_equil,"."))[end] == "eqdsk") || ((split(filepath_equil,"."))[end] == "geqdsk")
    M, wall = read_geqdsk(filepath_equil,clockwise_phi=false) # Assume counter-clockwise phi-direction
    jdotb = M.sigma # The sign of the dot product between the plasma current and the magnetic field

    # Extract timepoint information from .eqdsk/.geqdsk file
    eqdsk_array = split(filepath_equil,".")
    XX = (split(eqdsk_array[end-2],"-"))[end] # Assume format ...-XX.YYYY.eqdsk where XX are the seconds and YYYY are the decimals
    YYYY = eqdsk_array[end-1] # Assume format ...-XX.YYYY.eqdsk where XX are the seconds and YYYY are the decimals
    timepoint = XX*","*YYYY # Format XX,YYYY to avoid "." when including in filename of saved output
else # Otherwise, assume magnetic equilibrium is a saved .jld2 file
    myfile = jldopen(filepath_equil,false,false,false,IOStream)
    M = myfile["S"]
    wall = myfile["wall"]
    close(myfile)
    jdotb = (M.sigma_B0)*(M.sigma_Ip)

    if typeof(timepoint)==String && length(split(timepoint,","))==2
        timepoint = timepoint
    else
        timepoint = "00,0000" # Unknown timepoint for magnetic equilibrium
    end
end

if fileext_FI=="cdf"
    # If the user has specified a TRANSP .cdf file with pertaining NUBEAM fast-ion distribution data...
    # Load the time, and overwrite timepoint. TRANSP time data superseeds .eqdsk time data
    TIME = round((ncread(filepath_FI_distr,"TIME"))[1],digits=4)
    TIME_array = split("$(TIME)",".") # Will be on format XX.YYYY
    XX = TIME_array[1]
    YYYY = TIME_array[2]
    timepoint = XX*","*YYYY # Format XX,YYYY to avoid "." when including in filename of saved output
end
fileext_FI=="cdf" || fileext_thermal=="cdf"
    if fileext_FI_TRANSP_shot=="cdf"
        # If the user has specified a TRANSP .cdf file with pertaining NUBEAM fast-ion distribution data...
        # Load the time, and overwrite timepoint. TRANSP time data superseeds .eqdsk time data
        TIME = round((ncread(filepath_FI_TRANSP_shot,"TIME"))[1],digits=4)
        TIME_array = split("$(TIME)",".") # Will be on format XX.YYYY
        XX = TIME_array[1]
        YYYY = TIME_array[2]
        timepoint = XX*","*YYYY # Format XX,YYYY to avoid "." when including in filename of saved output
    end
end

## ---------------------------------------------------------------------------------------------
# If available, load instrumental response and process it accordingly
instrumental_response = false
if isfile(instrumental_response_filepath) # Returns true both for Strings and Vector{String}
    if typeof(instrumental_response_filepath)==String # If it is a single file, it must be a .jld2 file (specified in start file)
        myfile = jldopen(instrumental_response_filepath,false,false,false,IOStream)
        instrumental_response_matrix = myfile["response_matrix"]
        instrumental_response_input = myfile["input"]
        instrumental_response_output = myfile["output"]
        close(myfile)
    else # Otherwise, it must be a Vector{String}
        matrix_filepath = instrumental_response_filepath[1]
        input_filepath = instrumental_response_filepath[2]
        output_filepath = instrumental_response_filepath[3]

        py"""
        import numpy as np
        instrumental_response_matrix = np.loadtxt($matrix_filepath)
        instrumental_response_input = np.loadtxt($input_filepath)
        instrumental_response_output = np.loadtxt($output_filepath)
        """
        instrumental_response_matrix = py"instrumental_response_matrix"
        instrumental_response_input = py"instrumental_response_input"
        instrumental_response_output = py"instrumental_response_output"
    end
    instrumental_response = true
end

## ---------------------------------------------------------------------------------------------
# Printing script info and inputs
println("")
println("--------------------------------------------------------calcSpec.jl--------------------------------------------------------")
if !(tokamak=="")
    print("Tokamak specified: "*tokamak*"      ")
else
    print("Tokamak not specified.         ")
end
if !(TRANSP_id=="")
    print("TRANSP ID specified: "*TRANSP_id*"         ")
else
    print("TRANSP ID not specified."*"         ")
end
if !(diagnostic_filepath=="")
    println("Diagnostic specified: "*diagnostic_name)
else
    println("Diagnostic not specified. Spherical emission will be assumed.")
end
if instrumental_response
    println("instrumental_response_filepath specified. Diagnostic response included.")
else
    println("instrumental_response_filepath not specified. Diagnostic response not included.")
end
println("")
if include_flr_effects
    println("---> Finite Larmor radius (FLR) effects included? Yes!")
else
    println("---> Finite Larmor radius (FLR) effects included? No.")
end
println("")
println("Fast-ion distribution file specified: "*filepath_FI_distr)
println("Timepoint: "*timepoint*" seconds")
if interp
    println("Fast-ion distribution will be interpolated onto a ($(nE_ps),$(np_ps),$(nR_ps),$(nz_ps)) grid in (E,p,R,z) phase space.")
end
if !isempty(findall(x-> x!==nothing,phase_space_point_of_interest))
    const_inds = findall(x-> x!==nothing,phase_space_point_of_interest)
    println("It will be sampled while keeping the following dimensions and values constant: ")
    if 1 in const_inds
        println("--- E: $(phase_space_point_of_interest[1]) keV")
    end
    if 2 in const_inds
        println("--- p: $(phase_space_point_of_interest[2])")
    end
    if 3 in const_inds
        println("--- R: $(phase_space_point_of_interest[3])")
    end
    if 4 in const_inds
        println("--- z: $(phase_space_point_of_interest[4])")
    end
end
println("")
if !(filepath_thermal_distr=="")
    println("Thermal species profiles file specified: "*filepath_thermal_distr)
else
    println("No thermal species file specified. Using default temperature/density profiles with: ")
    println("Thermal species temperature on axis: $(thermal_temp_axis) keV")
    println("Thermal species density on axis: $(thermal_dens_axis) m^-3")
end
if fileext_FI=="cdf" || fileext_thermal=="cdf"
    println("TRANSP shot file to be used: "*filepath_TRANSP_shot)
end
if fileext_FI=="cdf" || fileext_thermal=="cdf"
    println("Fast-ion TRANSP file to be used: "*filepath_FI_TRANSP_shot)
end
println("Equilibrium file specified: "*filepath_equil)
println("")
if distributed
    println("Parallel computing will be used with $(nprocs()) processes (1 main + $(nprocs()-1) workers).")
else
    println("Single-threaded computing.")
end
println("")
println("Fusion reaction specified: "*reaction_full)
println("Fast-ion species specified: "*split(reaction,"-")[2])
if emittedParticleHasCharge
    println("The emitted "*getEmittedParticle(reaction_full)*" particle of the "*reaction_full*" reaction has non-zero charge!")
    println("The resulting energy distribution for "*getEmittedParticle(reaction_full)*" from the plasma as a whole will be computed.")
end
println("")
println("")
println("There will be $(length(range(Ed_min,stop=Ed_max,step=Ed_diff))-1) diagnostic energy bins with")
println("Lower diagnostic energy bound: $(Ed_min) keV")
println("Upper diagnostic energy bound: $(Ed_max) keV")
println("")
println("$(mc_samples) Monte-Carlo samples will be used, and they will be partitioned into")
println("chunks of size $(mc_chunk) samples.")
println("")
if calcProjVel
    println("Projected fast-ion velocities will be used to create the synthetic signal spectrum.")
    println("")
end
if verbose && addNoise
    println("Background noise will be added to the synthetic signal. It will have a magnitude of $(noiseLevel_b*100) % of the peak signal level.")
    println("Noise will also be added to the synthetical signal itself. It will have a magnitude of $(noiseLevel_s*100) % of the peak signal level.")
    println("")
end
println("Results will be saved to: ")
println(folderpath_o*"spec_"*tokamak*"_"*TRANSP_id*"_"*diagnostic_name*"_"*pretty2scpok(reaction_full)*".jld2")
println("")
println("Please remove previously saved files with the same file name (if any) prior to script completion. Quickly!")
println("")
println("If you would like to change any settings, please edit the start_calcSpec_template.jl file or similar.")
println("")
println("Written by Henrik Järleblad. Last maintained 2023-11-07.")
println("------------------------------------------------------------------------------------------------------------------------------")
println("")
## ---------------------------------------------------------------------------------------------
verbose && println("Loading Python packages...")
@everywhere begin
    py"""
    import h5py
    import numpy as np

    import forward
    import transp_dists
    import transp_output
    import vcone
    """
end

## If no thermal distribution has been specified, we are going to need the default temp. and dens. profiles
if filepath_thermal_distr==""
    @everywhere begin
        include("misc/temp_n_dens.jl")
    end
end

verbose && println("Defining diagnostic energy grid vector... ")
@everywhere begin
    py"""
    # Define the diagnostic energy bins to which samples will be binned
    Ed_bins = np.arange($Ed_min,$Ed_max,$Ed_diff) # diagnostic spectrum bin edges (e.g. keV)
    Ed_vals = 0.5*(Ed_bins[1:] + Ed_bins[:-1])  # bin centers (e.g. keV)
    """
end
## ---------------------------------------------------------------------------------------------
verbose && println("Defining the thermal species distribution (and TRANSP output)... ")
@everywhere begin
    py"""
    reaction = $reaction
    forwardmodel = forward.Forward($diagnostic_filepath) # Pre-initialize the forward model

    thermal_species = (reaction.split("-"))[0] # Assume first species specified in reaction to be the thermal species. For example, in 'p-t' the 'p' will be assumed thermal. [0] is the first element in Python...

    # Load thermal and/or fast-ion TRANSP data
    if $fileext_thermal=="cdf" and (not $fileext_FI=="cdf"):
        tr_out = transp_output.TranspOutput($TRANSP_id, step=1, out_file=$filepath_thermal_distr,fbm_files=[$filepath_FI_TRANSP_shot]) # Load the TRANSP shot file
        thermal_dist = transp_dists.Thermal(tr_out, ion=thermal_species) # Then load the thermal ion distribution from that .cdf file
    elif $fileext_thermal=="cdf" and $fileext_FI=="cdf":
        tr_out = transp_output.TranspOutput($TRANSP_id, step=1, # The TRANSP_id and the step number (always 1 for 1 fbm_file)
                                            out_file=$filepath_thermal_distr, # The thermal distribution file
                                            fbm_files=[$filepath_FI_distr]) # The fast-ion distribution file
        thermal_dist = transp_dists.Thermal(tr_out, ion=thermal_species) # Then load the thermal ion distribution from that .cdf file
    elif (not $fileext_thermal=="cdf") and $fileext_FI=="cdf":
        tr_out = transp_output.TranspOutput($TRANSP_id, step=1, # The TRANSP_id and the step number (always 1 for 1 fbm_file)
                                            out_file=$filepath_TRANSP_shot, # The TRANSP shot file
                                            fbm_files=[$filepath_FI_distr]) # The fast-ion distribution file...
        thermal_dist = "" # thermal distribution will have been specified in some other way
    else:
        thermal_dist = "" # Otherwise, just let the thermal_dist variable be the empty string
        # This doesn't mean that the thermal distribution will be non-existent. It just means the Python framework won't use
        # it's internal processes to sample from the thermal distribution. That is up to Julia instead.
    """
end

if fileext_thermal=="jld2"
    verbose && println("Creating thermal temperature and density interpolations objects... ")
    thermal_temp_itp = Interpolations.interpolate((ρ_pol_array,), thermal_temp_array, Gridded(Linear()))
    thermal_dens_itp = Interpolations.interpolate((ρ_pol_array,), thermal_dens_array, Gridded(Linear()))
    thermal_temp_etp = Interpolations.extrapolate(thermal_temp_itp,0) # Add what happens in extrapolation scenario. 0 means we assume vacuum scrape-off layer (SOL)
    thermal_dens_etp = Interpolations.extrapolate(thermal_dens_itp,0) # Add what happens in extrapolation scenario. 0 means we assume vacuum scrape-off layer (SOL)

    @everywhere thermal_temp_etp = $thermal_temp_etp
    @everywhere thermal_dens_etp = $thermal_dens_etp
end

## ---------------------------------------------------------------------------------------------
verbose && println("Defining the fast-ion distribution... ")
if fileext_FI=="h5"
    frdvols_cumsum_vector, subs, E_array, p_array, R_array, z_array, dE_vector, dp_vector, dR_vector, dz_vector, nfast = h5toSampleReady(filepath_FI_distr, interp=interp, nE_ps=nE_ps, np_ps=np_ps, nR_ps=nR_ps, nz_ps=nz_ps, slices_of_interest=phase_space_point_of_interest, verbose=verbose)
elseif fileext_FI=="jld2"
    frdvols_cumsum_vector, subs, E_array, p_array, R_array, z_array, dE_vector, dp_vector, dR_vector, dz_vector, nfast = jld2toSampleReady(filepath_FI_distr, interp=interp, nE_ps=nE_ps, np_ps=np_ps, nR_ps=nR_ps, nz_ps=nz_ps, slices_of_interest=phase_space_point_of_interest, verbose=verbose)
elseif fileext_FI=="cdf"
    py"""
    fastion_species = (reaction.split("-"))[1] # Assume second species specified in reaction to be the fast-ion species. For example, in 'p-t' the 't' will be assumed fast-ion. [1] is the second element in Python...
    if $verbose:print("Setting fast-ion distribution from TRANSP .cdf file... ")
    fbm_dist = transp_dists.FBM(tr_out, ion=fastion_species, jdotb=$jdotb)
    fbm_dist.EP_sampling_method = 'inv_cdf' # Setting sampling method for Python sampling
    """
    fbm_dist = py"fbm_dist"
else
    println("Fast-ion distribution file format: ."*fileext_FI)
    error("Unknown fast-ion distribution file format. Please re-specify file and re-try.")
end
if fileext_FI=="h5" || fileext_FI=="jld2"
    verbose && println("Transferring defined fast-ion quantities to external processes... ")
    @everywhere frdvols_cumsum_vector = $frdvols_cumsum_vector
    @everywhere subs = $subs
    @everywhere E_array = $E_array
    @everywhere p_array = $p_array
    @everywhere R_array = $R_array
    @everywhere z_array = $z_array
    @everywhere dE_vector = $dE_vector
    @everywhere dp_vector = $dp_vector
    @everywhere dR_vector = $dR_vector
    @everywhere dz_vector = $dz_vector
    @everywhere nfast = $nfast
end

if fileext_FI=="cdf"
    verbose && println("Transferring defined fast-ion quantities to external processes... ")
    @everywhere fbm_dist = $fbm_dist # From main to externals in Julia
    @everywhere begin
        py"""
        fbm_dist = $fbm_dist # From main to externals in Python
        """
    end
end

verbose && println("Going into the distributed part... ")
if distributed # If parallel computating is desired...
    if visualizeProgress # if you want the progress to be visualized...
        prog = Progress(length(mc_sample_chunks)) # Create a progress bar that is length(mc_sample_chunks) long
        channel = RemoteChannel(()->Channel{Bool}(mc_samples), 1) # Utilize a channel
        spec_tot = fetch(@sync begin
            @async while take!(channel)
                ProgressMeter.next!(prog)
            end
            @async begin
                spec = @distributed (+) for mc_sample_chunk in mc_sample_chunks # Sampling all mc_samples in chunks. To balance RAM memory usage with computational speed

                    E = Array{Float64}(undef,(mc_sample_chunk,)) # Pre-allocate an empty vector of length mc_sample_chunk for energy samples
                    p = Array{Float64}(undef,(mc_sample_chunk,)) # Pre-allocate an empty vector of length mc_sample_chunk for pitch samples
                    R = Array{Float64}(undef,(mc_sample_chunk,)) # Pre-allocate an empty vector of length mc_sample_chunk for R samples
                    z = Array{Float64}(undef,(mc_sample_chunk,)) # Pre-allocate an empty vector of length mc_sample_chunk for z samples
                    B = Array{Float64}(undef,(3, mc_sample_chunk)) # Pre-allocate an empty matrix of size (3,mc_sample_chunk) for the magnetic field

                    if fileext_FI=="h5" || fileext_FI=="jld2"
                        rand_cumsum_values = rand(mc_sample_chunk)*frdvols_cumsum_vector[end] # Get mc_sample_chunk number of random samples from the cumulative-sum fr vector
                        for i=1:mc_sample_chunk
                            rand_cumsum_value = rand_cumsum_values[i]
                            j = searchsortedfirst(frdvols_cumsum_vector,rand_cumsum_value,Base.Order.Forward)
                            inds = collect(Tuple(subs[j])) # Cumsum sample index
                            rand_hyperdiff = rand(4) .- 0.5 # 4 stands for the number of dimensions. 0.5 to sample within a hypercube
                            E[i] = max(E_array[inds[1]] + rand_hyperdiff[1]*dE_vector[inds[1]], 0.0)
                            p[i] = p_array[inds[2]] + rand_hyperdiff[2]*dp_vector[inds[2]]
                            R[i] = R_array[inds[3]] + rand_hyperdiff[3]*dR_vector[inds[3]]
                            z[i] = z_array[inds[4]] + rand_hyperdiff[4]*dz_vector[inds[4]]
                        end
                        w = (nfast/mc_samples) .* ones(size(E)) # Each weight has to be equal to the total number of fast ions, divided by the number of Monte-Carlo samples
                    else
                        py"""
                        # Sample from the fast ion distribution
                        R, z, E, p = fbm_dist.sample(n_samples=$mc_sample_chunk) # Acquire mc_sample_chunk samples from the fast-ion distribution
                        w = fbm_dist.N * np.ones($mc_sample_chunk) # Each weight has to be equal to the total number of fast ions, since the forward model expects sum(w)=total number of fast ions.
                        """
                        E, p, R, z, w = py"E", py"p", py"R", py"z", (py"w")/mc_samples # Convert from Python to Julia, and normalize weights by number of Monte-Carlo samples
                    end

                    if typeof(py"thermal_dist")<:AbstractString # This boolean operator '<:' means 'a subtype of'
                        thermal_temp = Array{Float64}(undef,(mc_sample_chunk,)) # Pre-allocate an empty vector of length mc_sample_chunk for thermal temperatures
                        thermal_dens = Array{Float64}(undef,(mc_sample_chunk,)) # Pre-allocate an empty vector of length mc_sample_chunk for thermal densities
                        for i=1:mc_sample_chunk
                            B[:,i] = collect(reshape(Bfield(M,R[i],z[i]),3,1)) # Calculating magnetic field for R,z samples. Need to do this in Julia, not Python, because of need to include poloidal magnetic field
                            ψ_rz = M.psi_rz(R[i],z[i]) # Get the interpolated poloidal flux function at the R,z point
                            ρ_pol_rz = sqrt(ψ_rz-M.psi[1])/(M.psi[end]-M.psi[1]) # The formula for the normalized flux coordinate ρ_pol = (ψ-ψ_axis)/(ψ_edge-ψ_axis)

                            if fileext_thermal=="jld2"
                                thermal_temp[i] = thermal_temp_etp(ρ_pol_rz) # Interpolate onto ρ_pol_rz using the data from the .jld2 file
                                thermal_dens[i] = thermal_dens_etp(ρ_pol_rz) # -||-
                            elseif filepath_thermal_distr==""
                                thermal_temp[i] = getAnalyticalTemp(thermal_temp_axis, ρ_pol_rz)
                                thermal_dens[i] = getAnalyticalDens(thermal_dens_axis, ρ_pol_rz)
                            else
                                error("This should not be possible. Please double-check your inputs and re-try.")
                            end
                        end
                    else
                        for i=1:mc_sample_chunk
                            B[:,i] = collect(reshape(Bfield(M,R[i],z[i]),3,1)) # Calculating magnetic field for R,z samples. Need to do this in Julia, not Python, because of need to include poloidal magnetic field
                        end
                        thermal_temp = "" # If thermal_dist is defined, then we don't need to manually define thermal_temp (forward.py will take care of it automatically)
                        thermal_dens = "" # If thermal_dist is defined, then we don't need to manually define thermal_dens (forward.py will take care of it automatically)
                    end

                    # Compute the spectra from the E,p,R,z MC sample
                    py"""
                    spec_i = forwardmodel.calc($E, $p, $R, $z, $w, thermal_dist, Ed_bins, $B, n_repeat=50, reaction=reaction, bulk_temp=$thermal_temp, bulk_dens=$thermal_dens, flr=$include_flr_effects)
                    """

                    put!(channel, true) # Update the progress bar
                    spec_i = vec(py"spec_i") # Convert from Python to Julia (and vectorize, just in case)
                    spec_i # Declare this spectrum to be added to the parallel computing reduction (@distributed (+))
                end
                put!(channel, false) # Update the progress bar
                spec # Declare the total spectrum as ready for fetching
            end
        end)
    else
        spec_tot = @distributed (+) for mc_sample_chunk in mc_sample_chunks # Sampling all mc_samples in chunks. To balance RAM memory usage with computational speed

            E = Array{Float64}(undef,(mc_sample_chunk,)) # Pre-allocate an empty vector of length mc_sample_chunk for energy samples
            p = Array{Float64}(undef,(mc_sample_chunk,)) # Pre-allocate an empty vector of length mc_sample_chunk for pitch samples
            R = Array{Float64}(undef,(mc_sample_chunk,)) # Pre-allocate an empty vector of length mc_sample_chunk for R samples
            z = Array{Float64}(undef,(mc_sample_chunk,)) # Pre-allocate an empty vector of length mc_sample_chunk for z samples
            B = Array{Float64}(undef,(3, mc_sample_chunk)) # Pre-allocate an empty matrix of size (3,mc_sample_chunk) for the magnetic field

            if fileext_FI=="h5" || fileext_FI=="jld2"
                rand_cumsum_values = rand(mc_sample_chunk)*frdvols_cumsum_vector[end] # Get mc_sample_chunk number of random samples from the cumulative-sum fr vector
                for i=1:mc_sample_chunk
                    rand_cumsum_value = rand_cumsum_values[i]
                    j = searchsortedfirst(frdvols_cumsum_vector,rand_cumsum_value,Base.Order.Forward)
                    inds = collect(Tuple(subs[j])) # Cumsum sample index
                    rand_hyperdiff = rand(4) .- 0.5 # 4 stands for the number of dimensions. 0.5 to sample within a hypercube
                    E[i] = max(E_array[inds[1]] + rand_hyperdiff[1]*dE_vector[inds[1]], 0.0)
                    p[i] = p_array[inds[2]] + rand_hyperdiff[2]*dp_vector[inds[2]]
                    R[i] = R_array[inds[3]] + rand_hyperdiff[3]*dR_vector[inds[3]]
                    z[i] = z_array[inds[4]] + rand_hyperdiff[4]*dz_vector[inds[4]]
                end
                w = (nfast/mc_samples) .* ones(size(E)) # Each weight has to be equal to the total number of fast ions, divided by the number of Monte-Carlo samples
            else
                py"""
                # Sample from the fast ion distribution
                R, z, E, p = fbm_dist.sample(n_samples=$mc_sample_chunk) # Acquire mc_sample_chunk samples from the fast-ion distribution
                w = fbm_dist.N * np.ones($mc_sample_chunk) # Each weight has to be equal to the total number of fast ions, since the forward model expects sum(w)=total number of fast ions.
                """
                E, p, R, z, w = py"E", py"p", py"R", py"z", (py"w")/mc_samples # Convert from Python to Julia, and normalize weights by number of Monte-Carlo samples
            end

            if typeof(py"thermal_dist")<:AbstractString # This boolean operator '<:' means 'a subtype of'
                thermal_temp = Array{Float64}(undef,(mc_sample_chunk,)) # Pre-allocate an empty vector of length mc_sample_chunk for thermal temperatures
                thermal_dens = Array{Float64}(undef,(mc_sample_chunk,)) # Pre-allocate an empty vector of length mc_sample_chunk for thermal densities
                for i=1:mc_sample_chunk
                    B[:,i] = collect(reshape(Bfield(M,R[i],z[i]),3,1)) # Calculating magnetic field for R,z samples. Need to do this in Julia, not Python, because of need to include poloidal magnetic field
                    ψ_rz = M.psi_rz(R[i],z[i]) # Get the interpolated poloidal flux function at the R,z point
                    ρ_pol_rz = sqrt(ψ_rz-M.psi[1])/(M.psi[end]-M.psi[1]) # The formula for the normalized flux coordinate ρ_pol = (ψ-ψ_axis)/(ψ_edge-ψ_axis)

                    if fileext_thermal=="jld2"
                        thermal_temp[i] = thermal_temp_etp(ρ_pol_rz) # Interpolate onto ρ_pol_rz using the data from the .jld2 file
                        thermal_dens[i] = thermal_dens_etp(ρ_pol_rz) # -||-
                    elseif filepath_thermal_distr==""
                        thermal_temp[i] = getAnalyticalTemp(thermal_temp_axis, ρ_pol_rz)
                        thermal_dens[i] = getAnalyticalDens(thermal_dens_axis, ρ_pol_rz)
                    else
                        error("This should not be possible. Please double-check your inputs and re-try.")
                    end
                end
            else
                for i=1:mc_sample_chunk
                    B[:,i] = collect(reshape(Bfield(M,R[i],z[i]),3,1)) # Calculating magnetic field for R,z samples. Need to do this in Julia, not Python, because of need to include poloidal magnetic field
                end
                thermal_temp = "" # If thermal_dist is defined, then we don't need to manually define thermal_temp (forward.py will take care of it automatically)
                thermal_dens = "" # If thermal_dist is defined, then we don't need to manually define thermal_dens (forward.py will take care of it automatically)
            end

            # Compute the spectra from the E,p,R,z MC sample
            py"""
            spec_i = forwardmodel.calc($E, $p, $R, $z, $w, thermal_dist, Ed_bins, $B, n_repeat=50, reaction=reaction, bulk_temp=$thermal_temp, bulk_dens=$thermal_dens, flr=$include_flr_effects)
            """
            spec_i = vec(py"spec_i") # Convert from Python to Julia (and vectorize, just in case)
            spec_i # Declare this spectrum to be added to the parallel computing reduction (@distributed (+))
        end
    end
else # ... if you do not use multiple cores, good luck!
    spec_tot = zeros(length(py"Ed_vals"))
    for (mci, mc_sample_chunk) in enumerate(mc_sample_chunks) # Sampling all mc_samples in chunks. To balance RAM memory usage with computational speed
        verbose && println("Sampling chunk $(mci) of $(length(mc_sample_chunks))... ")

        E = Array{Float64}(undef,(mc_sample_chunk,)) # Pre-allocate an empty vector of length mc_sample_chunk for energy samples
        p = Array{Float64}(undef,(mc_sample_chunk,)) # Pre-allocate an empty vector of length mc_sample_chunk for pitch samples
        R = Array{Float64}(undef,(mc_sample_chunk,)) # Pre-allocate an empty vector of length mc_sample_chunk for R samples
        z = Array{Float64}(undef,(mc_sample_chunk,)) # Pre-allocate an empty vector of length mc_sample_chunk for z samples
        B = Array{Float64}(undef,(3, mc_sample_chunk)) # Pre-allocate an empty matrix of size (3,mc_sample_chunk) for the magnetic field

        if fileext_FI=="h5" || fileext_FI=="jld2"
            rand_cumsum_values = rand(mc_sample_chunk)*frdvols_cumsum_vector[end] # Get mc_sample_chunk number of random samples from the cumulative-sum fr vector
            for i=1:mc_sample_chunk
                rand_cumsum_value = rand_cumsum_values[i]
                j = searchsortedfirst(frdvols_cumsum_vector,rand_cumsum_value,Base.Order.Forward)
                inds = collect(Tuple(subs[j])) # Cumsum sample index
                rand_hyperdiff = rand(4) .- 0.5 # 4 stands for the number of dimensions. 0.5 to sample within a hypercube
                E[i] = max(E_array[inds[1]] + rand_hyperdiff[1]*dE_vector[inds[1]], 0.0)
                p[i] = p_array[inds[2]] + rand_hyperdiff[2]*dp_vector[inds[2]]
                R[i] = R_array[inds[3]] + rand_hyperdiff[3]*dR_vector[inds[3]]
                z[i] = z_array[inds[4]] + rand_hyperdiff[4]*dz_vector[inds[4]]
            end
            w = (nfast/mc_samples) .* ones(size(E)) # Each weight has to be equal to the total number of fast ions, divided by the number of Monte-Carlo samples
        else
            py"""
            # Sample from the fast ion distribution
            R, z, E, p = fbm_dist.sample(n_samples=$mc_sample_chunk) # Acquire mc_sample_chunk samples from the fast-ion distribution
            w = fbm_dist.N * np.ones($mc_sample_chunk) # Each weight has to be equal to the total number of fast ions, since the forward model expects sum(w)=total number of fast ions.
            """
            E, p, R, z, w = py"E", py"p", py"R", py"z", (py"w")/mc_samples # Convert from Python to Julia, and normalize weights by number of Monte-Carlo samples
        end

        if typeof(py"thermal_dist")<:AbstractString # This boolean operator '<:' means 'a subtype of'
            thermal_temp = Array{Float64}(undef,(mc_sample_chunk,)) # Pre-allocate an empty vector of length mc_sample_chunk for thermal temperatures
            thermal_dens = Array{Float64}(undef,(mc_sample_chunk,)) # Pre-allocate an empty vector of length mc_sample_chunk for thermal densities
            for i=1:mc_sample_chunk
                B[:,i] = collect(reshape(Bfield(M,R[i],z[i]),3,1)) # Calculating magnetic field for R,z samples. Need to do this in Julia, not Python, because of need to include poloidal magnetic field
                ψ_rz = M.psi_rz(R[i],z[i]) # Get the interpolated poloidal flux function at the R,z point
                ρ_pol_rz = sqrt(ψ_rz-M.psi[1])/(M.psi[end]-M.psi[1]) # The formula for the normalized flux coordinate ρ_pol = (ψ-ψ_axis)/(ψ_edge-ψ_axis)

                if fileext_thermal=="jld2"
                    thermal_temp[i] = thermal_temp_etp(ρ_pol_rz) # Interpolate onto ρ_pol_rz using the data from the .jld2 file
                    thermal_dens[i] = thermal_dens_etp(ρ_pol_rz) # -||-
                elseif filepath_thermal_distr==""
                    thermal_temp[i] = getAnalyticalTemp(thermal_temp_axis, ρ_pol_rz)
                    thermal_dens[i] = getAnalyticalDens(thermal_dens_axis, ρ_pol_rz)
                else
                    error("This should not be possible. Please double-check your inputs and re-try.")
                end
            end
        else
            for i=1:mc_sample_chunk
                B[:,i] = collect(reshape(Bfield(M,R[i],z[i]),3,1)) # Calculating magnetic field for R,z samples. Need to do this in Julia, not Python, because of need to include poloidal magnetic field
            end
            thermal_temp = "" # If thermal_dist is defined, then we don't need to manually define thermal_temp (forward.py will take care of it automatically)
            thermal_dens = "" # If thermal_dist is defined, then we don't need to manually define thermal_dens (forward.py will take care of it automatically)
        end

        # Compute the spectra from the E,p,R,z MC sample
        py"""
        spec_tot = forwardmodel.calc($E, $p, $R, $z, $w, thermal_dist, Ed_bins, $B, n_repeat=50, reaction=reaction, bulk_temp=$thermal_temp, bulk_dens=$thermal_dens, flr=$include_flr_effects)
        """

        global spec_tot .+= vec(py"spec_tot") # Convert from Python to Julia, vectorize and add spectrum sample to total spectrum
    end
end
println("Done with sampling!")
Ed_array = py"Ed_vals" # Extract the measurement bin centers from Python
## ---------------------------------------------------------------------------------------------
# Apply instrumental response to the computed signal
if instrumental_response
    spec_raw = deepcopy(spec_tot)
    Ed_array_raw = deepcopy(Ed_array)
    spec_tot = apply_instrumental_response(spec_raw,Ed_array_raw,instrumental_response_input,instrumental_response_output,instrumental_response_matrix)
    Ed_array = instrumental_response_output
end

## ---------------------------------------------------------------------------------------------
# Add noise to the signal, if that was specified
# Make sure there are no negative measurements
if addNoise
    S_clean = deepcopy(spec_tot)
    spec_tot, noise = add_noise(spec_tot,noiseLevel_b; k=noiseLevel_s)
    spec_tot = map(x-> x<0.0 ? 0.0 : x, spec_tot)
end
## ---------------------------------------------------------------------------------------------
verbose && println("Saving results to file... ")
global filepath_output_orig = folderpath_o*"spec_"*tokamak*"_"*TRANSP_id*"_at"*timepoint*"s_"
if !isempty(findall(x-> x!==nothing,phase_space_point_of_interest))
    const_inds = findall(x-> x!==nothing,phase_space_point_of_interest)
    if 1 in const_inds
        filepath_output_orig *= "E"*replace("$(phase_space_point_of_interest[1])","."=>",")*"keV_"
    end
    if 2 in const_inds
        filepath_output_orig *= "p"*replace("$(phase_space_point_of_interest[2])","."=>",")*"_"
    end
    if 3 in const_inds
        filepath_output_orig *= "R"*replace("$(phase_space_point_of_interest[3])","."=>",")*"m_"
    end
    if 4 in const_inds
        filepath_output_orig *= "z"*replace("$(phase_space_point_of_interest[4])","."=>",")*"m_"
    end
end
filepath_output_orig *= diagnostic_name*"_"*pretty2scpok(reaction_full)
global filepath_output = deepcopy(filepath_output_orig)
global count = 1
while isfile(filepath_output*".jld2") # To take care of not overwriting files. Add _(1), _(2) etc
    global filepath_output = filepath_output_orig*"_($(Int64(count)))"
    global count += 1 # global scope, to surpress warnings
end
global filepath_output = filepath_output*".jld2"
myfile = jldopen(filepath_output,true,true,false,IOStream)
write(myfile,"S",spec_tot)
write(myfile,"Ed_array",Ed_array)
if instrumental_response
    write(myfile,"S_raw",spec_raw)
    write(myfile,"Ed_array_raw",Ed_array_raw)
end
if addNoise
    write(myfile,"S_clean",S_clean)
    write(myfile,"noise",noise)
end
if fileext_FI=="h5" || fileext_FI=="jld2"
    write(myfile,"nfast",nfast)
else # Must have been .cdf file format for FI distribution
    write(myfile,"nfast",py"fbm_dist.N")
end
if calcProjVel
    write(myfile,"calcProjVel",calcProjVel)
end
if include_flr_effects
    write(myfile,"include_flr_effects",include_flr_effects)
end
close(myfile)

println("------ Done! ------")
