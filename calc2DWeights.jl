#########################################  calc2DWeights.jl #########################################

#### Description:
# This script computes 2D weight functions, it's as simple as that. In the current version of the OWCF, 
# the calc2DWeights.jl script is taylored to utilize the DRESS code (J. Eriksson et al, CPC, 199, 40-46, 2016) 
# to compute the weight functions. In future versions, it can be easily modified to instead save the weighted 
# (E,p) points in a file readable by e.g. FIDASIM. The calc2DWeights.jl script computes the weight functions
# on a grid in (E,p) space. The weights on the corresponding (v_para, v_perp) grid will be computed and
# saved if the user has set the saveVparaVperpWeights input variable to true.
# 
# The DRESS code is written in Python and calc2DWeights.jl utilizes the DRESS code via the Julia-Python
# interface package PyCall.jl
#
# The calc2DWeights.jl script computes weight functions for an equidistant, rectangular grid 
# in (E,p) space. Irregular/non-equidistant grid points are currently not supported.
#
# The (R,z) point of interest for the weights is chosen in the input file as R_of_interest and 
# z_of_interest. If set to :r_mag and :z_mag, the magnetic axis will be chosen as the point of 
# interest.
#
# An instrumental response function can be specified in the start file. The weight functions will
# then be saved both with and without instrumental response. Please see list of outputs below for further info.
#
### For the thermal species distribution (filepath_thermal_distr), you have three options.
## First, you can specify a TRANSP shot .cdf file (for example '94701V01.cdf'). The thermal
# species temperature and density profiles will then be extracted from the file. Please note, if you
# specify a .cdf TRANSP output file as filepath_thermal_distr, you have to specify either a .cdf TRANSP 
# fast-ion (NUBEAM) output file for the 'filepath_FI_cdf' input variable, or a valid value for the 
# 'timepoint' input variable. This is to let the script know the correct time windows to extract data from
# the TRANSP shot file. If none of these are specified when a .cdf TRANSP output file is specified 
# as the 'filepath_thermal_distr' input variable, the script will try to deduce the timepoint on its own.
#
## Secondly, you could instead specify a .jld2 file that contains information about the
# thermal species temperature and density profiles as a function of normalized flux coordinate
# ρ_pol. The .jld2 file must have the keys:
# 'thermal_temp' - The 1D array containing the thermal temperature for the ρ_pol grid points
# 'thermal_dens' - The 1D array containing the thermal density for the ρ_pol grid points
# 'rho_pol' - The 1D array containing the ρ_pol grid points
#
## Finally, you could also choose not to specify a thermal species distribution. The thermal
# species temperature and density profiles will then be set to default/temp profiles with specified
# thermal_temp_axis and thermal_dens_axis temperature and density on-axis values, respectively.
#
### The 2D weight functions can also be computed using the analytic equations implemented in the OWCF/forward.jl 
# script. The analytic equations are derived and described in the following publications:
#   - https://doi.org/10.1088/1741-4326/adc1df (two-step fusion reactions),  
#   - https://doi.org/10.1088/1741-4326/ad9bc8 (one-step fusion reactions). 
# If so, no thermal temperature data is needed, since the analytic equations assume zero temperature
# for the thermal species.
#
# Please see the start_calc2DW_template.jl file for further input information.

#### Inputs (Units given when defined in script)
# Given via input file start_calc2DW_template.jl, for example. 'template' should be replaced by whatever.

#### Outputs
# -

#### Saved files
# velWeights_[FLR]_[tokamak]_[TRANSP_id]_at[timepoint]s_[diagnostic_name]_[nEd]x[nE]x[np].jld2 - If iiimax == 1
# velWeights_[FLR]_[tokamak]_[TRANSP_id]_at[timepoint]s_[diagnostic_name]_$(i).jld2 - If iiimax != 1
# Regardless of saved file name, this saved file will have the fields:
#   W - The computed (E,p) weights. Dimensions are (nEd, nE, np) - Array{Float64,3}
#   E_array - The fast-ion energy grid array used for (E,p) space - Array{Float64,1}
#   p_array - The fast-ion p grid array used for (E,p) space - Array{Float64,1}
#   Ed_array - The diagnostic energy bin centers - Array{Float64,1}
#   reaction - The nuclear fusion reaction for which the weights are computed - String
#   filepath_thermal_distr - The filepath of the thermal species distribution. For reference - String
#   filepath_start - The filepath to the start file used to execute calc2DWeights.jl. For reference - String
# If plasma_rot was set to true, the output file will contain
#   plasma_rot_at_Rz - An (R,phi,z) vector with the plasma rotation at the (R,z) point of interest - Array{Float,1}
# If saveVparaVperpWeights was set to true, the output file will also contain
#   W_vel - The 2D weights but in (v_para,v_perp). Dimensions are (nEd, nvpara, nvperp) - Array{Float64,3}
#   vpara_array - The fast-ion vpara grid array used for (v_para, v_perp) space - Array{Float64,1}
#   vperp_array - The fast-ion vpara grid array used for (v_para, v_perp) space - Array{Float64,1}
# If an instrumental response function has been specified in the start file, the output file will also contain
#   W_raw - The computed (E,p) weights, without instrumental response. Dimensions are same as W - Array{Float64,3}
#   Ed_array_raw - The diagnostic bin centers, if there was no instrumental response - Array{Float64,1}
#   instrumental_response_input - The input grid points for the instrumental response model - Array{Float64,1}
#   instrumental_response_output - The output grid points for the instrumental response model - Array{Float64,1}
#   instrumental_response_matrix - The instrumental response matrix, i.e. the model - Matrix{Float64}
# If, in addition, saveVparaVperpWeights was also set to true, the output file will also contain 
#   W_vel_raw - The 2D weights but in (vpara,vperp), without instrumental response. Dimensions are same as W_vel - Array{Float64,3}

### Other
# Please note that the diagnostic energy grid will be created as bin centers.
# That is, the first diagnostic energy grid value will be (Ed_min+Ed_diff/2) and so on.

# Script written by Henrik Järleblad. Last maintained 2025-09-01.
################################################################################################

## ---------------------------------------------------------------------------------------------
verbose && println("Loading Julia packages... ")
@everywhere begin
    cd(folderpath_OWCF) # Necessary to move all the workers to the correct folder
    using EFIT # For calculating magn½etic equilibrium quantities
    using Equilibrium # For loading flux function data, tokamak geometry data etc.
    using FileIO # To write/open files in general
    using Interpolations # To be able to interpolate, if no thermal distribution is specified
    using JLD2 # To write/open .jld2 files (Julia files, basically)
    using NetCDF # To enable write/open .cdf files
    using Printf # To be able to print specific formats
    using ProgressMeter # To display computational progress during parallel computations
    using PyCall # For using Python code in Julia
    using SparseArrays # To enable utilization of sparse matrices/vectors
    plot_results && (using Plots) # To be able to plot results, if requested in start file
    include("misc/availReacts.jl") # To examine fusion reaction and extract thermal and fast-ion species
    include("misc/convert_units.jl") # To be able to work with units of measurement
    include("extra/dependencies.jl") # To utilize the OWCF dependencies library
    include("misc/rewriteReacts.jl") # To rewrite a fusion reaction from the A(b,c)D format to the A-b=c-D format
    include("misc/species_func.jl") # To convert species labels to particle mass
    include("misc/temp_n_dens.jl") # To use default OWCF temperature/density profiles, as well as load temp/dens from only TRANSP output file (without NUBEAM output file)
    pushfirst!(PyVector(pyimport("sys")."path"), "") # To add the forward, transp_dists, transp_output and vcone modules (scripts in current path)
end

if analytic
    verbose && println("---> 'analytic' input variable set to true. Loading analytic equations for forward modelling in forward.jl... ")
    thermal_temp_axis = 0.0 # keV. Just to be sure, since analytic equations require T_i = 0 keV
    @everywhere begin
        include("forward.jl")
        include("vcone.jl")
    end
end

## ---------------------------------------------------------------------------------------------
# Fusion reaction, diagnostic and particle-related checks
if getReactionForm(reaction)==3 && analytic 
    error("The 'analytic' input variable was set to true, but the 'reaction' input variable was specified on form (3)(fast-ion species only). This is not allowed. Please correct and re-try.")
end
if getReactionForm(reaction)==1 # If no excited energy state for the emitted particle (in the case it is an atomic nucleus) has been specified...
    verbose && println("No energy state specified for the emitted particle $(getEmittedParticle(reaction)). Assuming ground state (GS), if relevant... ")
    reaction *= "-GS"
end
if !reactionIsAvailable(reaction)
    error("Fusion reaction $(reaction) is not yet available in the OWCF. The following reactions are available: $(OWCF_AVAILABLE_FUSION_REACTIONS). For projected-velocity computations, the following particle species are available: $(OWCF_SPECIES). Please correct and re-try.")
end
if analytic && !reactionIsAvailableAnalytically(reaction)
    error("Expected spectra from fusion reaction $(reaction) is currently not available for computation via analytic equations. Currently analytically available fusion reactions include: $(OWCF_AVAILABLE_FUSION_REACTIONS_FOR_ANALYTIC_COMPUTATION). Please correct and re-try.")
end
@everywhere reaction = $reaction # Copy the reaction variable to all external (CPU) processes

emittedParticleHasCharge = false
if !(getSpeciesCharge(getEmittedParticle(reaction))==0) && lowercase(getEmittedParticleEnergyLevel(reaction))=="gs" && !(getReactionForm(reaction)==3)
    verbose && println("")
    verbose && println("The fusion product particle of interest ($(getEmittedParticle(reaction))) of the $(reaction) reaction has non-zero charge!")
    verbose && println("- A 1-step or 2-step gamma-ray reaction is NOT assumed, since the energy level of the $(getEmittedParticle(reaction)) particle is specified to be in ground state (GS).")
    verbose && println("- Computation of 2D weight functions from projected velocities is also NOT assumed, since the 'reaction' input variable is NOT on form (3)($(reaction))")
    verbose && println("---> For emitted particles with non-zero charge, the OWCF currently only supports computing the expected energy spectrum from the plasma as a whole (4*pi emission).")
    verbose && println("---> Therefore, the 'diagnostic_name' and 'diagnostic_filepath' input variables will be forcibly set to \"\".")
    verbose && println("")
    diagnostic_name = ""
    diagnostic_filepath = ""
    emittedParticleHasCharge = true
end

if !(lowercase(getEmittedParticleEnergyLevel(reaction))=="gs") && (diagnostic_filepath=="")
    error("The fusion product particle of interest ($(getEmittedParticle(reaction))) of the $(reaction) reaction is on an excited state, but no diagnostic line-of-sight was specified (the diagnostic_filepath input variable was left unspecified). This is not yet allowed. Please correct and re-try.")
end
if !(lowercase(getEmittedParticleEnergyLevel(reaction))=="gs")
    verbose && println("")
    verbose && println("The fusion product particle of interest ($(getEmittedParticle(reaction))) of the $(reaction) reaction is on an excited state!")
    verbose && println("The OWCF will calculate the spectrum of gamma-rays emitted from the de-excitation of the $(getEmittedParticle(reaction)) particle, towards the specified detector.")
    verbose && println("")
end


## ---------------------------------------------------------------------------------------------
# Determine filepath_thermal_distr file extension
fileext_thermal = lowercase((split(filepath_thermal_distr,"."))[end]) # Assume last part after final '.' is the file extension
fileext_FI_cdf = lowercase((split(filepath_FI_cdf,"."))[end]) # Assume last part after final '.' is the file extension
@everywhere fileext_thermal = $fileext_thermal
@everywhere fileext_FI_cdf = $fileext_FI_cdf

# Work with output filename
if !(typeof(filename_o)==String)
    error("The 'filename_o' input variable was specified as $(filename_o). This is not a String type. Please correct and re-try.")
end

## ---------------------------------------------------------------------------------------------
# Error checks
if (nprocs()>1) && !distributed
    error(ErrorException("Number of processes greater than 1, but single-threaded computation specified. Please set distributed to true."))
end

if !(fileext_thermal=="cdf" || fileext_thermal=="jld2" || fileext_thermal=="")
    error("Unknown thermal distribution file format ($(filepath_thermal_distr)). Please re-specify file and re-try.")
end

## ---------------------------------------------------------------------------------------------
# Determine fast-ion and thermal (thermal) species from inputs in start file
thermal_reactant, fast_reactant = String.(getFusionReactants(reaction)) # Check the specified fusion reaction, and extract thermal and fast-ion species
@everywhere thermal_reactant = $thermal_reactant # Transfer variable to all external processes
@everywhere fast_reactant = $fast_reactant # Transfer variable to all external processes

# projVel variable. To clarify when 2D weight functions are computed from projected velocities, in code below
projVel = false
if getReactionForm(reaction)==3 # If fusion reaction is specified as a single particle species..
    projVel = true # ...2D weight functions will be computed using projected velocities!
end
@everywhere projVel = $projVel # Export to all external CPU processes

## ---------------------------------------------------------------------------------------------
# Loading tokamak information and TRANSP RUN-ID
verbose && println("Loading thermal .jld2 file data or TRANSP info... ")

if fileext_thermal=="jld2"
    myfile = jldopen(filepath_thermal_distr,false,false,false,IOStream)
    thermal_temp_array = myfile["thermal_temp"]
    thermal_dens_array = myfile["thermal_dens"]
    ρ_pol_array = myfile["rho_pol"]
    @everywhere thermal_temp_array = $thermal_temp_array # These are sent to external processes here, for efficiency
    @everywhere thermal_dens_array = $thermal_dens_array # These are sent to external processes here, for efficiency
    @everywhere ρ_pol_array = $ρ_pol_array # These are sent to external processes here, for efficiency
    close(myfile)
end
# Loading tokamak information and TRANSP RUN-ID from thermal .cdf file
if fileext_thermal=="cdf"
    verbose && println("Loading TRANSP id information... ")
    TRANSP_id = (split((split(filepath_thermal_distr,"."))[1],"/"))[end] # Assume part between last '/' and '.' is the TRANSP id
else
    TRANSP_id = ""
end
@everywhere TRANSP_id = $TRANSP_id
timepoint_source = "UNKNOWN" # Initially assume we initially know nothing about where we got the timepoint data from

if filepath_thermal_distr=="" && (!(typeof(thermal_temp_axis)==Float64) && !(typeof(thermal_temp_axis)==Int64)) && !projVel
    @everywhere thermal_temp_axis = 3.0
    @warn "filepath_thermal_distr was not specified, and thermal_temp_axis was not specified correctly. thermal_temp_axis will be set to default value of 3.0 keV."
end
if filepath_thermal_distr=="" && !(typeof(thermal_dens_axis)==Float64) && !projVel
    @everywhere thermal_dens_axis = 1.0e20
    @warn "filepath_thermal_distr was not specified, and thermal_dens_axis was not specified correctly. thermal_dens_axis will be set to default value of 1.0e20 m^-3."
end

@everywhere TRANSP_id = $TRANSP_id
@everywhere tokamak = $tokamak

## ---------------------------------------------------------------------------------------------
# Loading tokamak equilibrium and timepoint
timepoint_start = deepcopy(timepoint) # Keep timepoint value set in start file

verbose && println("Loading magnetic equilibrium... ")
M, wall, jdotb = nothing, nothing, nothing # Initialize global magnetic equilibrium variables
try
    global M; global wall; global jdotb # Declare global scope
    M, wall = read_geqdsk(filepath_equil,clockwise_phi=false) # Assume counter-clockwise phi-direction
    jdotb = M.sigma # The sign of the dot product between the plasma current and the magnetic field

    # Extract timepoint information from .eqdsk/.geqdsk file
    eqdsk_array = split(filepath_equil,".")
    try
        global timepoint; global timepoint_source # Declare global scope
        XX = (split(eqdsk_array[end-2],"-"))[end] # Assume format ...-XX.YYYY.eqdsk where XX are the seconds and YYYY are the decimals
        YYYY = eqdsk_array[end-1] # Assume format ...-XX.YYYY.eqdsk where XX are the seconds and YYYY are the decimals
        timepoint = "$(XX*","*YYYY)" # Format XX,YYYY to avoid "." when including in filename of saved output
        timepoint_source = "EQDSK"
        verbose && println("---> Found timepoint data in magnetic equilibrium file! Loading... ")
    catch
        global timepoint; global timepoint_source # Declare global scope
        timepoint_source, timepoint = "UNKNOWN", "00,0000" # (SOURCE, VALUE). Unknown timepoint for magnetic equilibrium
    end
catch # Otherwise, assume magnetic equilibrium is a saved .jld2 file
    global M; global wall; global jdotb; global timepoint_source; global timepoint; local myfile
    myfile = jldopen(filepath_equil,false,false,false,IOStream)
    M = myfile["S"]
    wall = myfile["wall"]
    close(myfile)
    jdotb = (M.sigma_B0)*(M.sigma_Ip)

    if typeof(timepoint)==String && length(split(timepoint,","))==2
        timepoint_source, timepoint = "STARTFILE", timepoint # (SOURCE, VALUE)
    else
        timepoint_source, timepoint = "UNKNOWN", "00,0000" # (SOURCE, VALUE). Unknown timepoint for magnetic equilibrium
    end
end
psi_axis, psi_bdry = psi_limits(M) # The limits of the flux function
@everywhere M = $M
@everywhere wall = $wall
@everywhere jdotb = $jdotb
@everywhere psi_axis = $psi_axis
@everywhere psi_bdry = $psi_bdry

if tokamak=="JET" && parse(Float64,replace(timepoint,","=>"."))>=40.0
    # Standard JET pulse needed 40 seconds to prepare.
    # This is sometimes included.
    # The OWCF (and TRANSP) does not count this preparation time.
    # Therefore, deduct it.
    timepoint = replace("$(round(parse(Float64,replace(timepoint,","=>"."))-40.0,digits=6))", "." => ",")
end

if (fileext_thermal=="cdf") && (fileext_FI_cdf=="cdf")
    verbose && println("Checking TRANSP .cdf files for timepoint data... ")
    # If the user has specified a TRANSP .cdf file with pertaining NUBEAM fast-ion distribution data...
    # Load the time, and overwrite timepoint. TRANSP time data superseeds .eqdsk time data
    TIME = round((ncread(filepath_FI_cdf,"TIME"))[1],digits=4)
    TIME_array = split("$(TIME)",".") # Will be on format XX.YYYY
    XX = TIME_array[1]
    YYYY = TIME_array[2]
    timepoint_source, timepoint = "TRANSP", XX*","*YYYY # (SOURCE, VALUE). Format XX,YYYY to avoid "." when including in filename of saved output
    if timepoint_source=="UNKNOWN"
        verbose && println("---> Timepoint data found! Loading... ")
    else
        verbose && println("---> Timepoint data found! Overwriting previously found timepoint data... ")
    end
end

if !isnothing(timepoint_start)
    verbose && println("Checking start file for timepoint data... ")
    if timepoint_source=="UNKNOWN"
        verbose && println("---> Timepoint data found! Loading... ")
    else
        verbose && println("---> Timepoint data found! Overwriting previously found timepoint data... ")
    end
    timepoint_source, timepoint = "STARTFILE", timepoint_start # (SOURCE, VALUE)
end

## ---------------------------------------------------------------------------------------------
# Defining energy-pitch grid vectors
verbose && println("Defining energy-pitch grid vectors... ")
if !(E_array == nothing)
    Emin = minimum(E_array)
    Emax = maximum(E_array)
    nE = length(E_array)
else
    E_array = collect(range(Emin,stop=Emax,length=nE))
end
if !(p_array == nothing)
    p_min = minimum(p_array)
    p_max = maximum(p_array)
    np = length(p_array)
else
    p_array = collect(range(p_min, stop=p_max, length=np))
end

## ---------------------------------------------------------------------------------------------
# Python code (load Python packages for DRESS code)
# This is how you write Python code in Julia: py""" [Python code] """
verbose && println("Loading Python modules... ")
@everywhere begin
    py"""
    import os.path
    import numpy as np
    import forward
    import spec
    import vcone
    """
end

## ---------------------------------------------------------------------------------------------
# If available, load instrumental response and process it accordingly
global instrumental_response # Declare global scope
instrumental_response = false
if isfile(instrumental_response_filepath) # Returns true both for Strings and Vector{String}
    verbose && println("Loading instrumental response... ")
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
# Pre-processing (R,z) of interest input data
if (R_of_interest==:r_mag)
    R_of_interest = magnetic_axis(M)[1]
end
if (z_of_interest==:z_mag)
    z_of_interest = magnetic_axis(M)[2]
end
@everywhere R_of_interest = $R_of_interest
@everywhere z_of_interest = $z_of_interest

## ---------------------------------------------------------------------------------------------
# Load/compute the B-field at the (R,z) point of interest
B_at_Rz = Equilibrium.Bfield(M,R_of_interest,z_of_interest) # Calculate the B-field vector at the (R,z) point
B_at_Rz = [B_at_Rz[1], B_at_Rz[2], B_at_Rz[3]] # Re-structure for Python
B_at_Rz = reshape(B_at_Rz,(length(B_at_Rz),1)) # Re-structure for Python
@everywhere B_at_Rz = $B_at_Rz # Export to all external CPU processes

## ---------------------------------------------------------------------------------------------
# If specified, include plasma rotation in the weight computations
# If TRANSP data available and desired, load plasma rotation data from TRANSP file
if plasma_rot
    if plasma_rot_speed_data_source_type == :TRANSP
        include("misc/load_TRANSP_interp_object.jl") # Need to load plasma rotation interpolation object from TRANSP
        omega_itp = getTRANSPInterpObject("OMEGA",parse(Float64,replace(timepoint,","=>".")),plasma_rot_speed_data_source)
        ψ_rz = M(R_of_interest,z_of_interest)
        psi_on_axis, psi_at_bdry = psi_limits(M)
        ρ_pol_rz = sqrt((ψ_rz-psi_on_axis)/(psi_at_bdry-psi_on_axis)) # The formula for the normalized flux coordinate ρ_pol = sqrt((ψ-ψ_axis)/(ψ_edge-ψ_axis))
        omega = omega_itp(ρ_pol_rz) # Get omega at (R,z) of interest
        plasma_rot_speed = (magnetic_axis(M)[1])*omega # v = R * Ω (in m/s)
    elseif plasma_rot_speed_data_source_type == :MANUAL
        plasma_rot_speed = plasma_rot_speed
    else
        error("Invalid data source type for the plasma rotation. Please correct and re-try.")
    end
    if plasma_rot_dir == :BFIELD
        plasma_rot_unit_vec = B_at_Rz ./norm(B_at_Rz) # Plasma is assumed to rotate in the direction of the B-field
    elseif plasma_rot_dir == :TOROIDAL
        plasma_rot_unit_vec = [0.0, jdotb*sign(B_at_Rz[2]), 0.0] # Plasma is assumed to rotate in the toroidal direction. Clockwise/counter-clockwise is determined by the plasma current
    else
        error("Invalid direction for plasma rotation. Please correct and re-try.")
    end
    plasma_rot_at_Rz = plasma_rot_speed .* plasma_rot_unit_vec
else
    plasma_rot_at_Rz = [0.0, 0.0, 0.0]
end
plasma_rot_at_Rz = reshape(plasma_rot_at_Rz, (length(plasma_rot_at_Rz),1)) # Reshape for Python
@everywhere plasma_rot_at_Rz = $plasma_rot_at_Rz # Export to all external CPU processes

## ---------------------------------------------------------------------------------------------
# Printing script info and inputs
println("")
println("-------------------------------------------------calc2DWeights.jl-------------------------------------------------")
if !(tokamak=="")
    print("Tokamak: "*tokamak*"      ")
else
    print("Tokamak: N/A"*"      ")
end
if !(TRANSP_id=="")
    print("TRANSP ID: "*TRANSP_id*"      ")
else
    print("TRANSP ID: N/A."*"      ")
end
if !(diagnostic_name=="")
    print("Diagnostic (name): "*diagnostic_name*"      ")
else
    print("Diagnostic(name): N/A"*"      ")
end
println("Timepoint: $(timepoint) seconds (from $(timepoint_source))")
println("")
if projVel
    println("Forward model: Projected velocities")
elseif analytic
    println("Forward model: OWCF/forward.jl")
else
    println("Forward model: DRESS")
end
println("")
if !projVel
    println("Fusion reaction: "*reaction)
else
    println("Fusion reaction: {projected velocities}")
end
println("---> Fast-ion species: "*fast_reactant)
if emittedParticleHasCharge && !projVel
    println("---> The emitted "*getEmittedParticle(reaction)*" particle of the "*reaction*" reaction has non-zero charge!")
    println("------> The resulting energy distribution for "*getEmittedParticle(reaction)*" from the plasma as a whole will be computed.")
end
println("---> Velocity samples per (E,p) point: $(n_gyro)")
print("---> Finite Larmor radius (FLR) effects included? ")
if include_FLR_effects
    println("Yes!")
else
    println("No.")
end
print("---> Plasma rotation included? ")
if plasma_rot
    println("Yes!")
    print("------> Plasma rotation vector: ")
    println("$(round.(reshape(plasma_rot_at_Rz,(1,3)),sigdigits=4)) m/s")
    analytic && (@warn "Plasma rotation is currently not supported for computations using analytic equations (i.e. expressions in OWCF/forward.jl). Plasma rotation will be ignored")
else
    println("No.")
end
println("")
if !(diagnostic_filepath=="")
    println("Diagnostic (model): "*diagnostic_filepath)
else
    println("Diagnostic (model): N/A (4*pi spherical emission)" )
end
println("---> There will be $(length(range(Ed_min,stop=Ed_max-Ed_diff,step=Ed_diff))-2) diagnostic measurement bins with")
su = projVel ? "m/s" : "keV"
println("------> Lower diagnostic measurement bound: $(Ed_min) $(su)")
println("------> Upper diagnostic measurement bound: $(Ed_max) $(su)")
if instrumental_response
    println("---> Instrumental response: " .*instrumental_response_filepath)
else
    println("---> Instrumental response: No.")
end
println("")
if isfile(filepath_thermal_distr)
    println("Bulk (thermal) plasma profile (data): "*filepath_thermal_distr)
    println("---> Bulk (thermal) plasma species: "*thermal_reactant)
else
    if projVel
        println("No bulk (thermal) plasma distribution specified. Projected velocities of $(fast_reactant) ions to be computed.")
    else
        thermal_profile_type_available_options = [:DEFAULT, :FLAT]
        !(thermal_temp_profile_type in thermal_profile_type_available_options) && error("'thermal_temp_profile_type' input variable was not correctly specified. Currently available options include $(thermal_profile_type_available_options). Please correct and re-try.")
        !(thermal_dens_profile_type in thermal_profile_type_available_options) && error("'thermal_dens_profile_type' input variable was not correctly specified. Currently available options include $(thermal_profile_type_available_options). Please correct and re-try.")
        println("No bulk (thermal) plasma distribution file specified. Using default OWCF temperature+density profiles with: ")
        println("---> Thermal ($(thermal_reactant)) temperature on-axis: $(thermal_temp_axis) keV")
        println("------> Thermal ion ($(thermal_reactant)) temperature profile: $(thermal_temp_profile_type)")
        println("---> Thermal ($(thermal_reactant)) density on-axis: $(thermal_dens_axis) m^-3")
        println("------> Thermal ion ($(thermal_reactant)) density profile: $(thermal_dens_profile_type)")
    end
end
println("")
println("Magnetic (B-field) equilibrium: "*filepath_equil)
println("---> |B| on-axis: $(round(norm(Equilibrium.Bfield(M,magnetic_axis(M)...)),digits=2)) T")
println("")
sw = iiimax > 1 ? "s" : ""
println("2D weight function$(sw) will be computed for a $(length(E_array))x$(length(p_array)) (E,p)-space grid with")
println("Energy: [$(minimum(E_array)),$(maximum(E_array))] keV")
println("Pitch: [$(minimum(p_array)),$(maximum(p_array))]")
println("The (R,z) point of interest is: ($(R_of_interest),$(z_of_interest)) m")
println("---> $(iiimax) weight matrix/matrices will be computed.")
if saveVparaVperpWeights
    println("---> Weight functions will be saved in both (E,p) and (vpara,vperp) coordinates.")
end
println("")
println("Results will be saved to: ")
sFLR = include_FLR_effects ? "FLR_" : ""; @everywhere sFLR = $sFLR
filepath_o_s_fp = !(filename_o=="") ? filename_o : ("velWeights_$(sFLR)"*tokamak*"_"*TRANSP_id*"_at"*timepoint*"s_"*diagnostic_name*"_"*pretty2scpok(reaction; projVel = projVel))
if iiimax == 1
    println(folderpath_o*filepath_o_s_fp*"_$(length(range(Ed_min,stop=Ed_max,step=Ed_diff))-1)x$(nE)x$(np).jld2")
else
    println(folderpath_o*filepath_o_s_fp*"_1.jld2")
    println("... ")
    println(folderpath_o*filepath_o_s_fp*"_$(iiimax).jld2")
    if iii_average
        println("---> Average of all files will be computed and saved to: ")
        println("---> "*folderpath_o*filepath_o_s_fp*".jld2")
    end
end
println("")
if distributed
    println("Parallel computing will be used with $(nprocs()) processes (1 main + $(nprocs()-1) workers).")
else
    println("Single-threaded computing the weights... Good luck!")
end
if debug
    println("")
    println("!!!!!! DEBUGGING SPECIFIED. ALGORITHM WILL DEBUG. !!!!!!")
    println("")
end
println("")
println("If you would like to change any settings, please edit the start_calc2DW_template.jl file or similar.")
println("")
println("Written by Henrik Järleblad. Last maintained 2025-07-11.")
println("--------------------------------------------------------------------------------------------------------------------")
println("")

## ---------------------------------------------------------------------------------------------
# Loading TRANSP data (if any) and initializing the forward model of choice
verbose && println("Loading TRANSP data (if any) and initializing forward model on all distributed workers/processes... ")
if analytic
    verbose && println("---> Loading diagnostic viewing cone model data... ")
    @everywhere begin
        viewing_cone_model = ViewingCone(diagnostic_filepath)
    end
end
@everywhere begin
    py"""
    # Start by assuming that the bulk (thermal) plasma distribution is "", i.e. to be provided by custom non-TRANSP data
    thermal_dist = ""
    """
end
@everywhere begin
    py"""
    # The '$' in front of many Python variables means that the variable is defined in Julia, not in Python.
    test_thermal_particle = spec.Particle($thermal_reactant) # Check so that thermal species is available in DRESS code
    projVel = $projVel

    # Load TRANSP simulation data
    if (($fileext_thermal=="cdf") and ($fileext_FI_cdf=="cdf")) and (not projVel): # If there is a pair of TRANSP-bulk+TRANSP-FI files specified and we do not want to simply compute projected velocities...
        import transp_output
        import transp_dists
        $verbose and print("---> Loading TRANSP output from TRANSP files... ")
        tr_out = transp_output.TranspOutput($TRANSP_id, step=1, out_file=$filepath_thermal_distr,fbm_files=[$filepath_FI_cdf]) # Load the TRANSP shot file. Assume first step. This is likely to be extended in the future.
        $verbose and print("---> Setting bulk (thermal) plasma distribution... ")
        thermal_dist = transp_dists.Thermal(tr_out, ion=$thermal_reactant) # Then change from "" to the thermal ion distribution from that .cdf file

    $verbose and print("---> Initializing Python forward model object with diagnostic viewing cone (if any), fusion reaction (if any) and bulk (thermal) plasma distribution (if any)... ") 
    forwardmodel = forward.Forward($diagnostic_filepath, $reaction, thermal_dist) # Initialize the forward model

    $verbose and print("---> Initializing diagnostic measurement bins... ")
    Ed_bin_edges = np.arange($Ed_min,$Ed_max,$Ed_diff) # diagnostic spectrum bin edges (keV or m/s)
    if len(Ed_bin_edges)==1: # Make sure that there are at least one lower and one upper bin edge
        dEd = (($Ed_max)-($Ed_min))/10
        Ed_bin_edges = np.arange($Ed_min,($Ed_max)+dEd,$Ed_diff)
    Ed_vals = 0.5*(Ed_bin_edges[1:] + Ed_bin_edges[:-1]) # bin centers (keV or m/s)
    nEd = len(Ed_vals)
    """
    nEd = py"nEd"
    Ed_array = vec(py"Ed_vals")
end

## ---------------------------------------------------------------------------------------------
# Pre-processing bulk (thermal) plasma density and temperature data

verbose && println("Preparing bulk (thermal) plasma temperature and density data... ")
(verbose && projVel) && println("---{projected velocities}---> (not used)")
if !(py"thermal_dist"=="")
    # If !(py"thermal_dist"==""), thermal ion density and temperature data will be loaded from TRANSP file
    (verbose && !(py"thermal_dist"=="")) && println("---> Thermal ion temperature and density data loaded from TRANSP file ($(filepath_thermal_distr))")
    thermal_dens_etp = nothing
    thermal_temp_etp = nothing
else
    if isfile(filepath_thermal_distr)
        # If a thermal ion data file was specified
        if fileext_thermal=="jld2"
            # And the thermal ion data file was a file with a .jld2 file extension
            thermal_temp_itp = Interpolations.interpolate((ρ_pol_array,), thermal_temp_array, Gridded(Linear()))
            thermal_dens_itp = Interpolations.interpolate((ρ_pol_array,), thermal_dens_array, Gridded(Linear()))
            thermal_temp_etp = Interpolations.extrapolate(thermal_temp_itp,0) # Add what happens in extrapolation scenario. 0 means we assume vacuum scrape-off layer (SOL)
            thermal_dens_etp = Interpolations.extrapolate(thermal_dens_itp,0) # Add what happens in extrapolation scenario. 0 means we assume vacuum scrape-off layer (SOL)
            verbose && println("---> Thermal ion temperature and density data loaded from $(filepath_thermal_distr)")
            (verbose && projVel) && println("---{projected velocities}---> (not used)")
        elseif fileext_thermal=="cdf"
            if !(typeof(timepoint)==String && length(split(timepoint,","))==2) && !(typeof(timepoint)==Float64) && !(typeof(timepoint)==Int64)
                error("$(filepath_thermal_distr) specified as bulk (thermal) plasma distribution, but no valid timepoint (in seconds) could be inferred. Possible solutions: (1) Manually specify timepoint in start file (2) Use an .eqdsk file with timepoint data included in file name")
            end
            if typeof(timepoint)==String
                tp = parse(Float64,replace(timepoint, "," => "."))
            end
            thermal_temp_etp = getTempProfileFromTRANSP(tp, filepath_thermal_distr, thermal_reactant) # Get the temperature from TRANSP as an interpolation object
            thermal_dens_etp = getDensProfileFromTRANSP(tp, filepath_thermal_distr, thermal_reactant) # Get the density from TRANSP as an interpolation object
        else
            error("This error should be impossible to reach. Please post an issue at www.github.com/JuliaFusion/OWCF")
        end
    else # No bulk (thermal) plasma distribution file specified in start file
        if thermal_temp_profile_type==:DEFAULT
            thermal_temp_etp = x-> getAnalyticalTemp(thermal_temp_axis,x)
            verbose && println("---> Default OWCF temperature profile with T_axis=$(thermal_temp_axis) keV will be used as thermal ion temperature data (more info in OWCF/misc/ folder)")
            (verbose && projVel) && println("---{projected velocities}---> (not used)")
        else # Otherwise, just use a flat bulk (thermal) plasma temperature profile
            thermal_temp_etp = x-> thermal_temp_axis
            # Don't print if projected velocities are to be computed. Then, bulk (thermal) temp+dens are irrelevant
            (verbose && !projVel) && println("---> Flat thermal ion temperature ($(thermal_temp_axis) keV) profile assumed")
        end
        if thermal_dens_profile_type==:DEFAULT
            thermal_dens_etp = x-> getAnalyticalDens(thermal_dens_axis,x)
            verbose && println("---> Default OWCF density profile with n_axis=$(thermal_dens_axis) m^-3 will be used as thermal ion density data (more info in OWCF/misc/ folder)")
        else # Otherwise, just use a flat bulk (thermal) plasma density profile
            thermal_dens_etp = x-> thermal_dens_axis
            # Don't print if projected velocities are to be computed. Then, bulk (thermal) temp+dens are irrelevant
            (verbose && !projVel) && println("---> Flat thermal ion density ($(thermal_dens_axis) m^-3) profile assumed")
        end
    end
end
@everywhere thermal_dens_etp = $thermal_dens_etp
@everywhere thermal_temp_etp = $thermal_temp_etp

verbose && println("Defining bulk (thermal) plasma temperature and density functions... ")
(verbose && projVel) && println("---{projected velocities}---> (not used)")
@everywhere begin
    if py"thermal_dist"=="" # Custom data
        function getThermalSpeciesTemperature(R,z)
            ψ_rz = M.(R,z) # Get the interpolated poloidal flux function at the R,z point
            ρ_pol_rz = sqrt.(max.(0.0, (ψ_rz .- psi_axis) ./(psi_bdry-psi_axis))) # The formula for the normalized flux coordinate ρ_pol = (ψ-ψ_axis)/(ψ_edge-ψ_axis)
            return thermal_temp_etp.(ρ_pol_rz) # Interpolate onto ρ_pol_rz using the data from the .jld2 file
        end
        function getThermalSpeciesDensity(R,z)
            ψ_rz = M.(R,z) # Get the interpolated poloidal flux function at the R,z point
            ρ_pol_rz = sqrt.(max.(0.0, (ψ_rz .- psi_axis) ./(psi_bdry-psi_axis))) # The formula for the normalized flux coordinate ρ_pol = (ψ-ψ_axis)/(ψ_edge-ψ_axis)
            return thermal_dens_etp.(ρ_pol_rz) # Interpolate onto ρ_pol_rz using the data from the .jld2 file
        end
    else # TRANSP data
        function getThermalSpeciesTemperature(R,z)
            temp = py"thermal_dist.get_temperature"(R,z)
            if length(R)==1
                return temp[1]
            end
            return temp
        end
        function getThermalSpeciesDensity(R,z)
            dens = py"thermal_dist.get_density"(R,z)
            if length(R)==1
                return dens[1]
            end
            return dens
        end
    end
end

verbose && println("---> Pre-computing bulk (thermal) plasma density and temperature at (R,z)=($(R_of_interest) m, $(z_of_interest) m)... ")
(verbose && projVel) && println("---{projected velocities}---> (not used)")
R_of_interest = vcat(R_of_interest); z_of_interest = vcat(z_of_interest) # Vectorize scalars
bulk_temp_at_Rz = vcat(getThermalSpeciesTemperature(R_of_interest, z_of_interest))
bulk_dens_at_Rz = vcat(getThermalSpeciesDensity(R_of_interest, z_of_interest))
@everywhere R_of_interest = $R_of_interest
@everywhere z_of_interest = $z_of_interest
@everywhere bulk_temp_at_Rz = $bulk_temp_at_Rz
@everywhere bulk_dens_at_Rz = $bulk_dens_at_Rz

## ---------------------------------------------------------------------------------------------
# Prepare functions for computing synthetic measurements, based on inputs in start file
verbose && println("Creating synthetic measurements computation function... ")
@everywhere begin
    function helper_func(E, p, w)
        E = vcat(E); p = vcat(p); w = vcat(w); N = length(E) # Assume length(E)==length(p)==length(w)
        R = repeat(R_of_interest, N); z = repeat(z_of_interest, N)
        B_vecs = repeat(B_at_Rz, inner=(1,N))
        return E, p, R, z, w, B_vecs
    end
    if include_FLR_effects # Need to take (R,z) variation into account
        if analytic # Analytic equations in OWCF/forward.jl as forward model
            function compute_measurements(E, p, w)
                # Before gyro sampling
                E, p, R, z, w, B_vecs = helper_func(E, p, w)

                # Gyro sampling
                py"""
                _, x = forwardmodel.add_gyration($E, $p, $R, $z, $B_vecs, n_gyro=$n_gyro) # No need to specify particle species, already specified via Forward.__init__()
                """

                # After gyro sampling
                E_gyro = repeat(E, inner=n_gyro) # Vector (M,) where M = length(E)*n_gyro
                p_gyro = repeat(p, inner=n_gyro) # Vector (M,)
                R_gyro = py"x[0,:]" # Vector (M,)
                z_gyro = py"x[2,:]" # Vector (M,)
                w_gyro = inv(n_gyro) .*repeat(w, inner=n_gyro) # Vector (N*n_gyro,), weights rescaled by number of gyro-orbit samples
                B_gyro = reduce(hcat,map(i-> Equilibrium.Bfield(M,R_gyro[i],z_gyro[i]), eachindex(R_gyro))) # Create a (3,M) array with the magnetic field vector at each (R,z) point
                bulk_dens = getThermalSpeciesDensity(R_gyro, z_gyro) # Bulk plasma densities in m^-3

                # Call analytic function
                # Please note! Ed_array are the diagnostic measurement bin centers, NOT the edges as in py"Ed_bin_edges" for the DRESS code
                return analytic_calc(viewing_cone_model, E_gyro, p_gyro, R_gyro, z_gyro, w_gyro, Ed_array, B_gyro; reaction=reaction, bulk_dens=bulk_dens)
            end
        else # DRESS code as forward model
            if py"thermal_dist"=="" # Custom bulk (thermal) plasma distribution data
                function compute_measurements(E, p, w)
                    # Before gyro sampling
                    E, p, R, z, w, B_vecs = helper_func(E, p, w)

                    # Gyro sampling
                    py"""
                    v, x = forwardmodel.add_gyration($E, $p, $R, $z, $B_vecs, n_gyro=$n_gyro) # No need to specify particle species, already specified via Forward.__init__()
                    """

                    # After gyro sampling
                    R_gyro = py"x[0,:]" # Vector (M,) where M = length(E)*n_gyro
                    z_gyro = py"x[2,:]" # Vector (M,)
                    v = py"v" .+ plasma_rot_at_Rz # Add plasma rotation
                    w_gyro = inv(n_gyro) .*repeat(w, inner=n_gyro) # Vector (N*n_gyro,), weights rescaled by number of gyro-orbit samples
                    B_gyro = reduce(hcat,map(i-> Equilibrium.Bfield(M,R_gyro[i],z_gyro[i]), eachindex(R_gyro))) # Create a (3,M) array with the magnetic field vector at each (R,z) point
                    bulk_temp = getThermalSpeciesTemperature(R_gyro, z_gyro) # Bulk plasma temperature in keV
                    bulk_dens = getThermalSpeciesDensity(R_gyro, z_gyro) # Bulk plasma density in m^-3

                    # Call DRESS
                    py"""
                    spec = forwardmodel.compute_spectrum(Ed_bin_edges, $R_gyro, $z_gyro, $v, $w_gyro, $B_gyro, bulk_temp=$bulk_temp, bulk_dens=$bulk_dens)
                    """

                    return vcat(py"spec")
                end
            else # TRANSP bulk (thermal) plasma distribution data
                function compute_measurements(E, p, w)
                    # Before gyro sampling
                    E, p, R, z, w, B_vecs = helper_func(E, p, w)
                    
                    # Gyro sampling
                    py"""
                    v, x = forwardmodel.add_gyration($E, $p, $R, $z, $B_vecs, n_gyro=$n_gyro) # No need to specify particle species, already specified via Forward.__init__()
                    """

                    # After gyro sampling
                    R_gyro = py"x[0,:]" # Vector (M,) where M = length(E)*n_gyro
                    z_gyro = py"x[2,:]" # Vector (M,)
                    v = py"v" .+ plasma_rot_at_Rz # Add plasma rotation
                    w_gyro = inv(n_gyro) .*repeat(w, inner=n_gyro)
                    B_gyro = reduce(hcat,map(i-> Equilibrium.Bfield(M,R_gyro[i],z_gyro[i]), eachindex(R_gyro))) # Create a (3,M) array with the magnetic field vector at each (R,z) point
                    # bulk_temp <--- Already in forwardmodel via instantiation
                    # bulk_dens <--- -||-

                    # Call DRESS
                    py"""
                    spec = forwardmodel.compute_spectrum(Ed_bin_edges, $R_gyro, $z_gyro, $v, $w_gyro, $B_gyro)
                    """

                    return vcat(py"spec")
                end
            end
        end
    else # (R,z) is fixed, since 2D weight functions
        if analytic # Analytic equations in OWCF/forward.jl as forward model
            function compute_measurements(E, p, w)
                # Before gyro sampling
                E, p, R, z, w, B_vecs = helper_func(E, p, w)

                # No sampling needed because of analytic equations and no FLR effects

                # After gyro sampling
                bulk_dens = repeat(bulk_dens_at_Rz, length(E)) # Vector (N,)

                # Call analytic function
                # Please note! Ed_array are the diagnostic measurement bin centers, NOT the edges as in py"Ed_bin_edges" for the DRESS code
                return analytic_calc(viewing_cone_model, E, p, R, z, w, Ed_array, B_vecs; reaction=reaction, bulk_dens=bulk_dens)
            end
        else # DRESS code as forward model
            if py"thermal_dist"=="" # Custom bulk (thermal) plasma distribution data
                function compute_measurements(E, p, w)
                    # Before gyro sampling
                    E, p, R, z, w, B_vecs = helper_func(E, p, w)

                    # Gyro sampling, to get velocity vectors
                    py"""
                    v, _ = forwardmodel.add_gyration($E, $p, $R, $z, $B_vecs, n_gyro=$n_gyro) # No need to specify particle species, already specified via Forward.__init__()
                    """

                    # After gyro sampling
                    R_gc = repeat(R, inner=n_gyro) # Ignore FLR effects
                    z_gc = repeat(z, inner=n_gyro) # Ignore FLR effects
                    v = py"v" .+ plasma_rot_at_Rz # Add plasma rotation
                    w_gc = inv(n_gyro) .*repeat(w, inner=n_gyro) # Vector (N*n_gyro,), weights rescaled by number of gyro-orbit samples
                    B_gc = repeat(B_vecs, inner=(1,n_gyro)) # Ignore FLR effects
                    bulk_temp = repeat(bulk_temp_at_Rz, inner=length(E)*n_gyro) # Ignore FLR effects
                    bulk_dens = repeat(bulk_dens_at_Rz, inner=length(E)*n_gyro) # Ignore FLR effects

                    # Call DRESS
                    py"""
                    spec = forwardmodel.compute_spectrum(Ed_bin_edges, $R_gc, $z_gc, $v, $w_gc, $B_gc, bulk_temp=$bulk_temp, bulk_dens=$bulk_dens)
                    """

                    return vcat(py"spec")
                end
            else # TRANSP bulk (thermal) plasma distribution data
                function compute_measurements(E, p, w)
                    # Before gyro sampling
                    E, p, R, z, w, B_vecs = helper_func(E, p, w)

                    # Gyro sampling, to get velocity vectors
                    py"""
                    v, _ = forwardmodel.add_gyration($E, $p, $R, $z, $B_vecs, n_gyro=$n_gyro) # No need to specify particle species, already specified via Forward.__init__()
                    """

                    # After gyro sampling
                    R_gc = repeat(R, inner=n_gyro) # Ignore FLR effects
                    z_gc = repeat(z, inner=n_gyro) # Ignore FLR effects
                    v = py"v" .+ plasma_rot_at_Rz # Add plasma rotation
                    w_gc = inv(n_gyro) .*repeat(w, inner=n_gyro) # Vector (N*n_gyro,), weights rescaled by number of gyro-orbit samples
                    B_gc = repeat(B_vecs, inner=(1,n_gyro)) # Ignore FLR effects
                    # bulk_temp <--- Already in forwardmodel via instantiation
                    # bulk_dens <--- -||-

                    # Call DRESS
                    py"""
                    spec = forwardmodel.compute_spectrum(Ed_bin_edges, $R_gc, $z_gc, $v, $w_gc, $B_gc)
                    """

                    return vcat(py"spec")
                end
            end
        end
    end
end

## ---------------------------------------------------------------------------------------------
E_iE_p_ip_array = zip(repeat(E_array,outer=length(p_array)),repeat(eachindex(E_array),outer=length(p_array)),repeat(p_array,inner=length(E_array)),repeat(eachindex(p_array),inner=length(E_array)))
npoints = length(E_iE_p_ip_array)

list_o_saved_filepaths = Vector{String}(undef,iiimax) # Create a list to keep track of all saved files (needed for averaging)

# Calculating the weights
verbose && println("Starting the "*diagnostic_name*" weights calculations... ")
for iii=1:iiimax
    verbose && println("iii: $(iii)")
    if distributed && !debug # If parallel computating is desired (and !debug)...
        if visualizeProgress # if you want the progress to be visualized...
            prog = Progress(npoints) # Create a progress bar the is E_iE_p_ip_array long
            channel = RemoteChannel(()->Channel{Bool}(npoints), 1) # Utilize a channel
            Wtot = fetch(@sync begin
                @async while take!(channel)
                    ProgressMeter.next!(prog)
                end
                @async begin
                    W = @distributed (+) for E_iE_p_ip in collect(E_iE_p_ip_array)
                        E = E_iE_p_ip[1]
                        iE = E_iE_p_ip[2]
                        p = E_iE_p_ip[3]
                        ip = E_iE_p_ip[4]
                        
                        spec = compute_measurements(E, p, 1.0)

                        W_iEip = zeros(length(spec),nE,np)
                        W_iEip[:,iE,ip] = spec
                        put!(channel, true) # Update the progress bar
                        W_iEip # Declare this orbit weight to be added to the parallel computing reduction (@distributed (+))
                    end
                    put!(channel, false) # Update the progress bar
                    W # Declare the total weight function as ready for fetching
                end
            end)
        else
            Wtot = @distributed (+) for E_iE_p_ip in collect(E_iE_p_ip_array)
                E = E_iE_p_ip[1]
                iE = E_iE_p_ip[2]
                p = E_iE_p_ip[3]
                ip = E_iE_p_ip[4]

                spec = compute_measurements(E, p, 1.0)

                W_iEip = zeros(length(spec),nE,np)
                W_iEip[:,iE,ip] = spec
                W_iEip # Declare this orbit weight to be added to the parallel computing reduction (@distributed (+))
            end
        end
    else # Single-threaded computing
        if debug
            verbose && println("Debugging specified. Only single-threaded debugging is possible. Will debug.")
            # WRITE CODE TO DEBUG QUANTITIES OF INTEREST

        else
            E_iE_p_ip = collect(E_iE_p_ip_array)[1]
            E = E_iE_p_ip[1]
            iE = E_iE_p_ip[2]
            p = E_iE_p_ip[3]
            ip = E_iE_p_ip[4]

            spec = compute_measurements(E, p, 1.0)

            Wtot = zeros(length(spec),nE,np)
            Wtot[:,iE,ip] = spec
        end

        for i=2:length(E_iE_p_ip_array)
            if debug
                verbose && println("Debugging (E,p) point $(i) of $(length(E_iE_p_ip_array))... ")
                # WRITE CODE TO DEBUG QUANTITIES OF INTEREST

            else
                verbose && println("Calculating spectra for (E,p) point $(i) of $(length(E_iE_p_ip_array))... ")
                E_iE_p_ip = collect(E_iE_p_ip_array)[i]
                E = E_iE_p_ip[1]
                iE = E_iE_p_ip[2]
                p = E_iE_p_ip[3]
                ip = E_iE_p_ip[4]

                spec = compute_measurements(E, p, 1.0)

                Wtot[:,iE,ip] = spec
            end
        end
    end

    if saveVparaVperpWeights
        verbose && println("Transforming (E,p) weights to (vpara,vperp)... ")
        nvpara = Int64(round(sqrt(nE*np))) # We know the default values for Ep2VparaVperp()
        nvperp = nvpara # We know the default values for Ep2VparaVperp()
        W_vel = zeros(nEd,nvpara,nvperp)
        global vpara_array; global vperp_array # Declare global scope, to avoid annoying warnings when executing OWCF/tests/run_tests.jl
        W_vel[1,:,:], vpara_array, vperp_array = Ep2VparaVperp(E_array, p_array, Wtot[1,:,:], returnAbscissas=true)
        for iEd=2:nEd
            W_vel[iEd,:,:] = Ep2VparaVperp(E_array, p_array, Wtot[iEd,:,:], returnAbscissas=false)
        end
    end

    if instrumental_response
        verbose && println("Applying diagnostic response to weight functions... ")
        Wtot_raw = deepcopy(Wtot) # To be able to save the weight functions without instrumental response as well
        Wtot = apply_instrumental_response(Wtot, Ed_array, instrumental_response_input, instrumental_response_output, instrumental_response_matrix)
        if saveVparaVperpWeights
            W_vel_raw = deepcopy(W_vel)
            W_vel = apply_instrumental_response(W_vel, Ed_array, instrumental_response_input, instrumental_response_output, instrumental_response_matrix)
        end
    end

    if !debug
        # Set the output file name 'filepath_out'
        # First, if the user specified a custom output file name (filename_o), use that instead of the default OWCF output file name
        filepath_o_s = !(filename_o=="") ? filename_o : ("velWeights_$(sFLR)"*tokamak*"_"*TRANSP_id*"_at$(replace("$(timepoint)", "," => "."))s_"*diagnostic_name*"_"*pretty2scpok(reaction; projVel = projVel))
        if iiimax==1 # If you intend to calculate only one weight matrix
            global filepath_output_orig = folderpath_o*filepath_o_s*"_$(length(Ed_array))x$(nE)x$(np)"
        else # If you intend to calculate several (identical) weight matrices
            global filepath_output_orig = folderpath_o*filepath_o_s*"_$(iii)"
        end
        global filepath_output = deepcopy(filepath_output_orig)
        global count = 1
        while isfile(filepath_output*".jld2") # To take care of not overwriting files. Add _(1), _(2) etc
            global filepath_output = filepath_output_orig*"_($(Int64(count)))"
            global count += 1 # global scope, to surpress warnings
        end

        if plot_results
            plot_font = "Computer Modern"
            Plots.default(fontfamily=plot_font)
            verbose && println("Plotting weight function data... ")

            if instrumental_response
                W_raw_plt = Wtot_raw # Bad variable names...
                W_plt = Wtot # Bad variable names...
                if saveVparaVperpWeights
                    W_vel_raw_plt = W_vel_raw # Bad variable names...
                    W_vel_plt = W_vel # Bad variable names...
                end
            else
                W_raw_plt = Wtot # Bad variable names...
                if saveVparaVperpWeights
                    W_vel_raw_plt = W_vel # Bad variable names...
                end
            end

            # Without instrumental response (raw)
            N_bins = length(Ed_array)
            if N_bins>=5
                # Algorithm to find suitable indices of Ed to visualize
                W_gross = dropdims(sum(W_raw_plt , dims=(2,3)), dims=(2,3))
                iEd_mid = argmax(W_gross)
                iEd_half = argmin(abs.(W_gross .- W_gross[iEd_mid]/2))
                if iEd_mid == N_bins
                    iEd_hig = iEd_mid
                    iEd_low = iEd_half
                    iEd_mid = iEd_mid - 1
                elseif iEd_mid == 1
                    iEd_hig = iEd_half
                    iEd_low = iEd_mid
                    iEd_mid = iEd_mid + 1
                else
                    iEd_diff = iEd_mid - iEd_half
                    iEd_low = iEd_mid - abs(iEd_diff)
                    iEd_hig = iEd_mid + abs(iEd_diff)
                end
                plt_raw_inds = [iEd_low, iEd_mid, iEd_hig]
            elseif N_bins==4
                plt_raw_inds = [2,3,4]
            elseif N_bins==3
                plt_raw_inds = [1,2,3]
            elseif N_bins==2
                plt_raw_inds = [1,1,2]
            else # N_bins==1
                plt_raw_inds = [1,1,1]
            end

            Ed_low = @sprintf "%.2E" Ed_array[plt_raw_inds[1]]; W_raw_max_low = @sprintf "%.2E" maximum(W_raw_plt[plt_raw_inds[1],:,:])
            Ed_mid = @sprintf "%.2E" Ed_array[plt_raw_inds[2]]; W_raw_max_mid = @sprintf "%.2E" maximum(W_raw_plt[plt_raw_inds[2],:,:])
            Ed_hi = @sprintf "%.2E" Ed_array[plt_raw_inds[3]]; W_raw_max_hi = @sprintf "%.2E" maximum(W_raw_plt[plt_raw_inds[3],:,:])
            plt_Ep_raw_low = Plots.heatmap(E_array, p_array, transpose(W_raw_plt[plt_raw_inds[1],:,:]),title="w($(Ed_low),E,p), max(w): $(W_raw_max_low) keV^-1")
            plt_Ep_raw_mid = Plots.heatmap(E_array, p_array, transpose(W_raw_plt[plt_raw_inds[2],:,:]),title="w($(Ed_mid),E,p), max(w): $(W_raw_max_mid) keV^-1")
            plt_Ep_raw_hi = Plots.heatmap(E_array, p_array, transpose(W_raw_plt[plt_raw_inds[3],:,:]),title="w($(Ed_hi),E,p), max(w): $(W_raw_max_hi) keV^-1")
            plt_Ep_raw = Plots.plot(plt_Ep_raw_low, plt_Ep_raw_mid, plt_Ep_raw_hi, layout=(1,3), size=(1200,400), dpi=200, xlabel="Energy [keV]", ylabel="Pitch [-]", 
                                fillcolor=cgrad([:white, :yellow, :orange, :red, :black]), titlefontsize=10, colorbar=false, bottom_margin=6Plots.mm,
                                left_margin=6Plots.mm)
            png(plt_Ep_raw,filepath_output*"_Ep_raw")

            if saveVparaVperpWeights
                W_vel_raw_max_low = @sprintf "%.2E" maximum(W_vel_raw_plt[plt_raw_inds[1],:,:])
                W_vel_raw_max_mid = @sprintf "%.2E" maximum(W_vel_raw_plt[plt_raw_inds[2],:,:])
                W_vel_raw_max_hi = @sprintf "%.2E" maximum(W_vel_raw_plt[plt_raw_inds[3],:,:])
                plt_VpaVpe_raw_low = Plots.heatmap(vpara_array, vperp_array, transpose(W_vel_raw_plt[plt_raw_inds[1],:,:]),title="w($(Ed_low),E,p), max(w): $(W_vel_raw_max_low) keV^-1")
                plt_VpaVpe_raw_mid = Plots.heatmap(vpara_array, vperp_array, transpose(W_vel_raw_plt[plt_raw_inds[2],:,:]),title="w($(Ed_mid),E,p), max(w): $(W_vel_raw_max_mid) keV^-1")
                plt_VpaVpe_raw_hi = Plots.heatmap(vpara_array, vperp_array, transpose(W_vel_raw_plt[plt_raw_inds[3],:,:]),title="w($(Ed_hi),E,p), max(w): $(W_vel_raw_max_hi) keV^-1")
                plt_VpaVpe_raw = Plots.plot(plt_VpaVpe_raw_low, plt_VpaVpe_raw_mid, plt_VpaVpe_raw_hi, layout=(1,3), size=(1200,400), dpi=200, xlabel="vpara [m/s]", ylabel="vperp [m/s]", 
                                    fillcolor=cgrad([:white, :yellow, :orange, :red, :black]), titlefontsize=10, colorbar=false, bottom_margin=6Plots.mm,
                                    left_margin=6Plots.mm)
                png(plt_VpaVpe_raw,filepath_output*"_VpaVpe_raw")
            end

            if instrumental_response
                # With instrumental response
                N_bins = length(instrumental_response_output)
                if N_bins>=5
                    # Algorithm to find suitable indices of Ed to visualize
                    W_gross = dropdims(sum(W_plt , dims=(2,3)), dims=(2,3))
                    iEd_mid = argmax(W_gross)
                    iEd_half = argmin(abs.(W_gross .- W_gross[iEd_mid]/2))
                    if iEd_mid == N_bins
                        iEd_hig = iEd_mid
                        iEd_low = iEd_half
                        iEd_mid = iEd_mid - 1
                    elseif iEd_mid == 1
                        iEd_hig = iEd_half
                        iEd_low = iEd_mid
                        iEd_mid = iEd_mid + 1
                    else
                        iEd_diff = iEd_mid - iEd_half
                        iEd_low = iEd_mid - abs(iEd_diff)
                        iEd_hig = iEd_mid + abs(iEd_diff)
                    end
                    plt_inds = [iEd_low, iEd_mid, iEd_hig]
                elseif N_bins==4
                    plt_inds = [2,3,4]
                elseif N_bins==3
                    plt_inds = [1,2,3]
                elseif N_bins==2
                    plt_inds = [1,1,2]
                else # N_bins==1
                    plt_inds = [1,1,1]
                end

                Ed_low = @sprintf "%.2E" instrumental_response_output[plt_inds[1]]; W_max_low = @sprintf "%.2E" maximum(W_plt[plt_inds[1],:,:])
                Ed_mid = @sprintf "%.2E" instrumental_response_output[plt_inds[2]]; W_max_mid = @sprintf "%.2E" maximum(W_plt[plt_inds[2],:,:])
                Ed_hi = @sprintf "%.2E" instrumental_response_output[plt_inds[3]]; W_max_hi = @sprintf "%.2E" maximum(W_plt[plt_inds[3],:,:])
                plt_Ep_low = Plots.heatmap(E_array, p_array, transpose(W_plt[plt_inds[1],:,:]),title="w($(Ed_low),E,p), max(w): $(W_max_low) $(units_inverse(instrumental_response_output_units))")
                plt_Ep_mid = Plots.heatmap(E_array, p_array, transpose(W_plt[plt_inds[2],:,:]),title="w($(Ed_mid),E,p), max(w): $(W_max_mid) $(units_inverse(instrumental_response_output_units))")
                plt_Ep_hi = Plots.heatmap(E_array, p_array, transpose(W_plt[plt_inds[3],:,:]),title="w($(Ed_hi),E,p), max(w): $(W_max_hi) $(units_inverse(instrumental_response_output_units))")
                plt_Ep = Plots.plot(plt_Ep_low, plt_Ep_mid, plt_Ep_hi, layout=(1,3), size=(1200,400), dpi=200, xlabel="Energy [keV]", ylabel="Pitch [-]", 
                                    fillcolor=cgrad([:white, :yellow, :orange, :red, :black]), titlefontsize=10, colorbar=false, bottom_margin=6Plots.mm,
                                    left_margin=6Plots.mm)
                png(plt_Ep,filepath_output*"_Ep")

                if saveVparaVperpWeights
                    W_vel_max_low = @sprintf "%.2E" maximum(W_vel_plt[plt_inds[1],:,:])
                    W_vel_max_mid = @sprintf "%.2E" maximum(W_vel_plt[plt_inds[2],:,:])
                    W_vel_max_hi = @sprintf "%.2E" maximum(W_vel_plt[plt_inds[3],:,:])
                    plt_VpaVpe_low = Plots.heatmap(vpara_array, vperp_array, transpose(W_vel_plt[plt_inds[1],:,:]),title="w($(Ed_low),E,p), max(w): $(W_vel_max_low) $(units_inverse(instrumental_response_output_units))")
                    plt_VpaVpe_mid = Plots.heatmap(vpara_array, vperp_array, transpose(W_vel_plt[plt_inds[2],:,:]),title="w($(Ed_mid),E,p), max(w): $(W_vel_max_mid) $(units_inverse(instrumental_response_output_units))")
                    plt_VpaVpe_hi = Plots.heatmap(vpara_array, vperp_array, transpose(W_vel_plt[plt_inds[3],:,:]),title="w($(Ed_hi),E,p), max(w): $(W_vel_max_hi) $(units_inverse(instrumental_response_output_units))")
                    plt_VpaVpe = Plots.plot(plt_VpaVpe_low, plt_VpaVpe_mid, plt_VpaVpe_hi, layout=(1,3), size=(1200,400), dpi=200, xlabel="vpara [m/s]", ylabel="vperp [m/s]", 
                                        fillcolor=cgrad([:white, :yellow, :orange, :red, :black]), titlefontsize=10, colorbar=false, bottom_margin=6Plots.mm,
                                        left_margin=6Plots.mm)
                    png(plt_VpaVpe,filepath_output*"_VpaVpe")
                end
            end
        end

        verbose && println("Saving weight function data to file... ")
        global filepath_output = filepath_output*".jld2"
        list_o_saved_filepaths[iii] = filepath_output # Store a record of the saved file
        myfile_s = jldopen(filepath_output,true,true,false,IOStream)
        write(myfile_s,"W",Wtot)
        write(myfile_s,"E_array",vec(E_array))
        write(myfile_s,"p_array",vec(p_array))
        write(myfile_s,"reaction",reaction)
        write(myfile_s,"R_of_interest",R_of_interest)
        write(myfile_s,"z_of_interest",z_of_interest)
        if projVel
            write(myfile_s,"projVel",projVel)
        end
        if saveVparaVperpWeights
            write(myfile_s,"W_vel",W_vel)
            write(myfile_s,"vpara_array",vpara_array)
            write(myfile_s,"vperp_array",vperp_array)
        end
        if instrumental_response
            write(myfile_s,"W_raw",Wtot_raw)
            write(myfile_s,"Ed_array",instrumental_response_output)
            write(myfile_s,"Ed_array_units",instrumental_response_output_units)
            write(myfile_s,"Ed_array_raw",Ed_array)
            write(myfile_s,"Ed_array_raw_units",projVel ? "m_s^-1" : "keV") # The raw output abscissa of calc2DWeights.jl is always in m/s or keV
            write(myfile_s,"instrumental_response_input",instrumental_response_input)
            write(myfile_s,"instrumental_response_output",instrumental_response_output)
            write(myfile_s,"instrumental_response_matrix",instrumental_response_matrix)
            if saveVparaVperpWeights
                write(myfile_s,"W_vel_raw",W_vel_raw)
            end
        else
            write(myfile_s,"Ed_array",Ed_array)
            write(myfile_s,"Ed_array_units",projVel ? "m_s^-1" : "keV") # Otherwise, the output abscissa of calcSpec.jl is always in m/s or keV
        end
        if analytic
            write(myfile_s,"analytic","Weight functions computed using the equations in e.g. A. Valentini et al, 2025, Nucl. Fusion 65, 046031")
        end
        if plasma_rot
            write(myfile_s,"plasma_rot_at_Rz",plasma_rot_at_Rz)
        end
        write(myfile_s,"filepath_thermal_distr",filepath_thermal_distr)
        write(myfile_s,"filepath_start",filepath_start)
        close(myfile_s)
    else
        verbose && println("Saving debugged quantities... ")
        # WRITE WHATEVER CODE TO SAVE THE DEBUGGED QUANTITIES
    end
end # The end of the iii=1:iiimax for-loop. (to enable computation of several identical orbit weight matrices, to post-analyze MC noise influence for example)
if iiimax>1 # If we were supposed to compute more than one weight matrix...
    if iii_average # If an average of all saved files were supposed to be computed... compute it!
        verbose && println("Computing average of all weight matrices... ")
        filepath_first_W = list_o_saved_filepaths[1]
        myfile = jldopen(filepath_first_W,false,false,false,IOStream)
        W_total = zeros(size(myfile["W"]))
        E_array = myfile["E_array"]
        p_array = myfile["p_array"]
        Ed_array = myfile["Ed_array"]
        reaction = myfile["reaction"]
        R_of_interest = myfile["R_of_interest"]
        z_of_interest = myfile["z_of_interest"]
        if projVel
            projVel = myfile["projVel"]
        end
        Ed_array_units = myfile["Ed_array_units"]
        if saveVparaVperpWeights
            W_vel_total = zeros(size(myfile["W_vel"]))
            vpara_array = myfile["vpara_array"]
            vperp_array = myfile["vperp_array"]
        end
        if instrumental_response
            W_raw_total = zeros(size(myfile["W_raw"]))
            Ed_array_raw = myfile["Ed_array_raw"]
            Ed_array_raw_units = myfile["Ed_array_raw_units"]
            instrumental_response_input = myfile["instrumental_response_input"]
            instrumental_response_output = myfile["instrumental_response_output"]
            instrumental_response_matrix = myfile["instrumental_response_matrix"]
            if saveVparaVperpWeights
                W_vel_raw_total = zeros(size(myfile["W_vel_raw"]))
            end
        end
        if analytic
            analytic = myfile["analytic"]
        end
        if plasma_rot
            plasma_rot_at_Rz = myfile["plasma_rot_at_Rz"]
        end
        filepath_thermal_distr = myfile["filepath_thermal_distr"]
        filepath_start = myfile["filepath_start"]
        close(myfile)

        for filepath in list_o_saved_filepaths
            local myfile = jldopen(filepath,false,false,false,IOStream)
            W_total[:,:,:] += myfile["W"]
            if saveVparaVperpWeights
                W_vel_total[:,:,:] += myfile["W_vel"]
            end
            if instrumental_response
                W_raw_total[:,:,:] += myfile["W_raw"]
                if saveVparaVperpWeights
                    W_vel_raw_total[:,:,:] += myfile["W_vel_raw"]
                end
            end
            close(myfile)
        end

        W_total ./= length(list_o_saved_filepaths)
        if saveVparaVperpWeights
            W_vel_total ./= length(list_o_saved_filepaths)
        end
        if instrumental_response
            W_raw_total ./= length(list_o_saved_filepaths)
            if saveVparaVperpWeights
                W_vel_raw_total ./= length(list_o_saved_filepaths)
            end
        end

        verbose && println("Success! Saving average in separate file... ")
        filepath_o_avg = reduce(*,map(x -> x*"_",split(filepath_first_W,"_")[1:end-1]))[1:end-1]*".jld2"
        myfile = jldopen(filepath_o_avg,true,true,false,IOStream)
        write(myfile,"W",W_total)
        write(myfile,"E_array",E_array)
        write(myfile,"p_array",p_array)
        write(myfile,"Ed_array",Ed_array)
        write(myfile,"reaction",reaction)
        write(myfile,"R_of_interest",R_of_interest)
        write(myfile,"z_of_interest",z_of_interest)
        if projVel
            write(myfile,"projVel",projVel)
        end
        write(myfile,"Ed_array_units",Ed_array_units)
        if saveVparaVperpWeights
            write(myfile,"W_vel",W_vel_total)
            write(myfile,"vpara_array",vpara_array)
            write(myfile,"vperp_array",vperp_array)
        end
        if instrumental_response
            write(myfile,"W_raw",W_raw_total)
            write(myfile,"Ed_array_raw",Ed_array_raw)
            write(myfile,"Ed_array_raw_units",Ed_array_raw_units)
            write(myfile,"instrumental_response_input",instrumental_response_input)
            write(myfile,"instrumental_response_output",instrumental_response_output)
            write(myfile,"instrumental_response_matrix",instrumental_response_matrix)
            if saveVparaVperpWeights
                write(myfile,"W_vel_raw",W_vel_raw_total)
            end
        end
        if analytic
            write(myfile,"analytic",analytic)
        end
        if plasma_rot
            write(myfile,"plasma_rot_at_Rz",plasma_rot_at_Rz)
        end
        write(myfile,"filepath_thermal_distr",filepath_thermal_distr)
        write(myfile,"filepath_start",filepath_start)
        close(myfile)
    end
end
println("------ Done! ------")
