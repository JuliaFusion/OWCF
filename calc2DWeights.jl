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
# specify a .cdf TRANSP shot file as filepath_thermal_distr, you have to specify a .cdf TRANSP fast-ion file
# for filepath_FI_cdf as well. This is to let the script know the correct time windows to extract from
# the TRANSP shot file.
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
# If thermal_profiles_type==:FLAT, flat thermal plasma temperature and density profiles will be used.
#
### The 2D weight functions can also be computed using the analytic equations implemented in the OWCF/forward.jl 
# script. https://doi.org/10.1088/1741-4326/adc1df (two-step) and https://doi.org/10.1088/1741-4326/ad9bc8 (one-step). 
# If so, no thermal temperature data is needed, since the analytic equations assume zero temperature
# for the thermal species.
#
# Please see the start_calc2DW_template.jl file for further input information.

#### Inputs (Units given when defined in script)
# Given via input file start_calc2DW_template.jl, for example. 'template' should be replaced by whatever.

#### Outputs
# -

#### Saved files
# velWeights_[tokamak]_[TRANSP_id]_at[timepoint]s_[diagnostic_name]_[nEd]x[nE]x[np].jld2 - If iiimax == 1
# velWeights_[tokamak]_[TRANSP_id]_at[timepoint]s_[diagnostic_name]_$(i).jld2 - If iiimax != 1
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

# Script written by Henrik Järleblad. Last maintained 2025-03-17.
################################################################################################

## ---------------------------------------------------------------------------------------------
verbose && println("Loading Julia packages... ")
@everywhere begin
    cd(folderpath_OWCF) # Necessary to move all the workers to the correct folder
    using PyCall # For using Python code in Julia
    using EFIT # For calculating magn½etic equilibrium quantities
    using Equilibrium # For loading flux function data, tokamak geometry data etc.
    using ProgressMeter # To display computational progress during parallel computations
    using JLD2 # To write/open .jld2 files (Julia files, basically)
    using FileIO # To write/open files in general
    using SparseArrays # To enable utilization of sparse matrices/vectors
    using NetCDF # To enable write/open .cdf files
    using Interpolations # To be able to interpolate, if no thermal distribution is specified
    include("misc/species_func.jl") # To convert species labels to particle mass
    include("misc/availReacts.jl") # To examine fusion reaction and extract thermal and fast-ion species
    include("misc/rewriteReacts.jl") # To rewrite a fusion reaction from the A(b,c)D format to the A-b=c-D format
    include("extra/dependencies.jl")
    pushfirst!(PyVector(pyimport("sys")."path"), "") # To add the forward, transp_dists, transp_output and vcone modules (scripts in current path)
end

if analytic
    verbose && println("'analytic' input variable set to true. Loading analytic equations for forward modelling in forward.jl... ")
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
fileext_thermal = (split(filepath_thermal_distr,"."))[end] # Assume last part after final '.' is the file extension
fileext_thermal = lowercase(fileext_thermal)
fileext_FI_cdf = (split(filepath_FI_cdf,"."))[end] # Assume last part after final '.' is the file extension
fileext_FI_cdf = lowercase(fileext_FI_cdf)
@everywhere fileext_thermal = $fileext_thermal
@everywhere fileext_FI_cdf = $fileext_FI_cdf

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
thermal_species, FI_species = getFusionReactants(reaction) # Check the specified fusion reaction, and extract thermal and fast-ion species
@everywhere thermal_species = $thermal_species # Transfer variable to all external processes
@everywhere FI_species = $FI_species # Transfer variable to all external processes

# projVel variable. To clarify when 2D weight functions are computed from projected velocities, in code below
projVel = false
if getReactionForm(reaction)==3 # If fusion reaction is specified as a single particle species..
    projVel = true # ...2D weight functions will be computed using projected velocities!
end

## ---------------------------------------------------------------------------------------------
# Loading tokamak information and TRANSP RUN-ID
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
end
# Loading tokamak information and TRANSP RUN-ID from thermal .cdf file
if fileext_thermal=="cdf"
    verbose && println("Loading TRANSP id information... ")
    TRANSP_id = (split((split(filepath_thermal_distr,"."))[1],"/"))[end] # Assume part between last '/' and '.' is the TRANSP id
else
    TRANSP_id = ""
end


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

if tokamak=="JET" && parse(Float64,replace(timepoint,","=>"."))>=40.0
    # Standard JET pulse needed 40 seconds to prepare.
    # This is sometimes included.
    # The OWCF (and TRANSP) does not count this preparation time.
    # Therefore, deduct it.
    timepoint = replace("$(parse(Float64,replace(timepoint,","=>"."))-40.0)", "." => ",")
end

if (fileext_thermal=="cdf") && (fileext_FI_cdf=="cdf")
    # If the user has specified a TRANSP .cdf file with pertaining NUBEAM fast-ion distribution data...
    # Load the time, and overwrite timepoint. TRANSP time data superseeds .eqdsk time data
    TIME = round((ncread(filepath_FI_cdf,"TIME"))[1],digits=4)
    TIME_array = split("$(TIME)",".") # Will be on format XX.YYYY
    XX = TIME_array[1]
    YYYY = TIME_array[2]
    timepoint = XX*","*YYYY # Format XX,YYYY to avoid "." when including in filename of saved output
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
println("----------------------------------------calc2DWeights.jl----------------------------------------")
if !(tokamak=="")
    print("Tokamak specified: "*tokamak*"      ")
else
    print("Tokamak not specified.         ")
end
if !(TRANSP_id=="")
    print("TRANSP ID specified: "*TRANSP_id*"      ")
else
    print("TRANSP ID not specified."*"      ")
end
println("Timepoint: "*timepoint*" seconds")
println("")
if !projVel
    println("Fusion reaction specified: "*reaction)
else
    println("2D weight functions will be computed using the projected velocity (u) of fast $(getFastParticleSpecies(reaction)) ions for all (E,p) points.")
end
println("Fast-ion species specified: "*FI_species)
if emittedParticleHasCharge && !projVel
    println("The emitted "*getEmittedParticle(reaction)*" particle of the "*reaction*" reaction has non-zero charge!")
    println("The resulting energy distribution for "*getEmittedParticle(reaction)*" from the plasma as a whole will be computed.")
end
println("")
print("---> Finite Larmor radius (FLR) effects included? ")
if flr_effects
    println("Yes!")
else
    println("No.")
end
println("")
if !(diagnostic_name=="")
    println("Diagnostic name specified: "*diagnostic_name)
end
if !(diagnostic_filepath=="")
    println("Diagnostic filepath specified: "*diagnostic_filepath)
else
    println("diagnostic_filepath not specified. Spherical emission will be assumed.")
end
if instrumental_response
    println("Instrumental response filepath: "*instrumental_response_filepath)
else
    println("instrumental_response_filepath not specified. Diagnostic response not included.")
end
println("")
print("---> Plasma rotation included? ")
if plasma_rot
    println("Yes!")
    print("Plasma rotation vector at (R,z) is (in m/s): ")
    println(round.(reshape(plasma_rot_at_Rz,(1,3)),sigdigits=4))
    println("")
else
    println("No.")
end
println("")
if isfile(filepath_thermal_distr)
    println("Thermal species $(getSlowParticleSpecies(reaction)) distribution file specified: "*filepath_thermal_distr)
else
    if projVel
        println("2D weight functions will be computed using the projected velocity (u) of the fast $(getFastParticleSpecies(reaction)) ions.")
    else
        if thermal_profiles_type==:FLAT
            println("Flat thermal plasma profiles will be assumed.")
            println("Thermal ion ($(thermal_species)) temperature will be set to $(thermal_temp_axis) keV.")
            println("Thermal ion ($(thermal_species)) density will be set to $(thermal_dens_axis) m^-3.")
        elseif thermal_profiles_type==:DEFAULT
            println("Default OWCF thermal plasma profiles will be used (OWCF/misc/default_temp_n_dens.png).")
            println("Thermal ion ($(thermal_species)) temperature on-axis will be set to $(thermal_temp_axis) keV.")
            println("Thermal ion ($(thermal_species)) density on axis will be set to $(thermal_dens_axis) m^-3.")
        else
            error("The 'thermal_profiles_type' input variable was not correctly specified (available options are :FLAT and :DEFAULT). Please correct and re-try.")
        end
    end
end
println("")
println("Magnetic equilibrium file specified: "*filepath_equil)
println("")
println("$(iiimax) weight matrix/matrices will be computed.")
println("")
println("2D weight function(s) will be computed for a $(length(E_array))x$(length(p_array)) (E,p)-space grid with")
println("Energy: [$(minimum(E_array)),$(maximum(E_array))] keV")
println("Pitch: [$(minimum(p_array)),$(maximum(p_array))]")
println("The (R,z) point of interest is: ($(R_of_interest),$(z_of_interest)) m")
println("")
if !projVel
    println("There will be $(length(range(Ed_min,stop=Ed_max,step=Ed_diff))-1) diagnostic energy bin(s) with")
    println("Lower diagnostic energy bound: $(Ed_min) keV")
    println("Upper diagnostic energy bound: $(Ed_max) keV")
else
    println("There will be $(length(range(Ed_min,stop=Ed_max,step=Ed_diff))-1) projected velocity (u) bin(s) with")
    println("Lower u bound: $(Ed_min) m/s")
    println("Upper u bound: $(Ed_max) m/s")
end
println("")
if saveVparaVperpWeights
    println("Weight functions will be saved in both (E,p) and (vpara,vperp) coordinates.")
    println("")
end
println("Results will be saved to: ")
if iiimax == 1
    println(folderpath_o*"velWeights_"*tokamak*"_"*TRANSP_id*"_at"*timepoint*"s_"*diagnostic_name*"_"*pretty2scpok(reaction; projVel = projVel)*"_$(length(range(Ed_min,stop=Ed_max,step=Ed_diff))-1)x$(nE)x$(np).jld2")
else
    println(folderpath_o*"velWeights_"*tokamak*"_"*TRANSP_id*"_at"*timepoint*"s_"*diagnostic_name*"_"*pretty2scpok(reaction; projVel = projVel)*"_1.jld2")
    println("... ")
    println(folderpath_o*"velWeights_"*tokamak*"_"*TRANSP_id*"_at"*timepoint*"s_"*diagnostic_name*"_"*pretty2scpok(reaction; projVel = projVel)*"_$(iiimax).jld2")
    if iii_average
        println("---> Average of all files will be computed and saved to: ")
        println("---> "*folderpath_o*"velWeights_"*tokamak*"_"*TRANSP_id*"_at"*timepoint*"s_"*diagnostic_name*"_"*pretty2scpok(reaction; projVel = projVel)*".jld2")
    end
end
println("")
if distributed
    println("Parallel computing will be used with $(nprocs()) processes (1 main + $(nprocs()-1) workers).")
else
    println("Single-threaded computing the weights... Good luck!")
end
println("")
if debug
    println("")
    println("!!!!!! DEBUGGING SPECIFIED. ALGORITHM WILL DEBUG. !!!!!!")
    println("")
end
println("Please remove previously saved files with the same file name (if any) prior to script completion. Quickly!")
println("")
println("If you would like to change any settings, please edit the start_calc2DW_template.jl file or similar.")
println("")
println("Written by Henrik Järleblad. Last maintained 2025-03-17.")
println("--------------------------------------------------------------------------------------------------")
println("")

## ---------------------------------------------------------------------------------------------
# Python code (essentially calling the forward model black box)
# This is how you write Python code in Julia: py""" [Python code] """
verbose && println("Loading Python modules... ")
@everywhere begin
    py"""
    import os.path
    import h5py
    import numpy as np
    from netCDF4 import Dataset

    import forward
    import transp_dists
    import transp_output
    import vcone
    """
end

## ---------------------------------------------------------------------------------------------
# If no thermal distribution file was specified, and default temp. and dens. profiles were specified,
# we need to load the OWCF/misc/temp_n_dens.jl script 
if !isfile(filepath_thermal_distr) && (thermal_profiles_type==:DEFAULT)
    @everywhere begin
        include("misc/temp_n_dens.jl")
    end
end

## ---------------------------------------------------------------------------------------------
# Specifying forward model
if analytic
    # If analytic equations are to be used as forward model to compute weight functions, we only need the diagnostic viewing cone
    myVC = ViewingCone(diagnostic_filepath)
else
    @everywhere begin
        py"""
        # The '$' in front of many Python variables means that the variable is defined in Julia, not in Python.
        reaction = $reaction
        test_thermal_particle = Particle($thermal_species) # Check so that thermal species is available in DRESS code
        thermal_species = $thermal_species
        projVel = $projVel
        if $verbose:
            print("From Python: Loading forward model with diagnostic... ") 
        forwardmodel = forward.Forward($diagnostic_filepath) # Pre-initialize the forward model
        """
    end
end
@everywhere begin
    py"""
    Ed_bin_edges = np.arange($Ed_min,$Ed_max,$Ed_diff) # diagnostic spectrum bin edges (keV or m/s)
    if len(Ed_bin_edges)==1: # Make sure that there are at least one lower and one upper bin edge
        dEd = (($Ed_max)-($Ed_min))/10
        Ed_bin_edges = np.arange($Ed_min,($Ed_max)+dEd,$Ed_diff)
    Ed_vals = 0.5*(Ed_bin_edges[1:] + Ed_bin_edges[:-1]) # bin centers (keV or m/s)
    nEd = len(Ed_vals)
    """
    Ed_array = Vector(py"Ed_vals"); nEd = length(Ed_array)
end

## ---------------------------------------------------------------------------------------------
# Pre-processing thermal density and temperature data
ψ_rz = M(R_of_interest,z_of_interest)
psi_on_axis, psi_at_bdry = psi_limits(M)
ρ_pol_rz = sqrt((ψ_rz-psi_on_axis)/(psi_at_bdry-psi_on_axis)) # The formula for the normalized flux coordinate ρ_pol = sqrt((ψ-ψ_axis)/(ψ_edge-ψ_axis))

if lowercase(fileext_thermal)=="jld2" && !projVel
    verbose && println("Setting thermal particle species temperature and density using data in $(filepath_thermal_distr)... ")
    thermal_temp_itp = Interpolations.interpolate((ρ_pol_array,), thermal_temp_array, Gridded(Linear()))
    thermal_dens_itp = Interpolations.interpolate((ρ_pol_array,), thermal_dens_array, Gridded(Linear()))
    thermal_temp_etp = Interpolations.extrapolate(thermal_temp_itp,0) # Add what happens in extrapolation scenario. 0 means we assume vacuum scrape-off layer (SOL)
    thermal_dens_etp = Interpolations.extrapolate(thermal_dens_itp,0) # Add what happens in extrapolation scenario. 0 means we assume vacuum scrape-off layer (SOL)
    
    thermal_temp = [thermal_temp_etp(ρ_pol_rz)]
    thermal_dens = [thermal_dens_etp(ρ_pol_rz)]
    thermal_dist = "" # No distribution. Temp and dens modelled by thermal_temp and thermal_dens
end

# If a thermal species distribution has not been specified, we simply need the thermal species temperature and density on axis
if !isfile(filepath_thermal_distr) && !projVel
    verbose && println("The 'filepath_thermal_distr' was not specified. Checking the 'thermal_profiles_type' input variable... ")
    if thermal_profiles_type==:FLAT
        verbose && println("Found :FLAT! Thermal profiles will be assumed constant throughout the plasma.")
        thermal_temp = [thermal_temp_axis]
        thermal_dens = [thermal_dens_axis]
    elseif thermal_profiles_type==:DEFAULT
        verbose && println("Found :DEFAULT! Thermal profiles will be as returned by getAnalyticalTemp() and getAnalyticalDens() functions in OWCF/misc/temp_n_dens.jl.")
        thermal_temp = [getAnalyticalTemp(thermal_temp_axis, ρ_pol_rz)]
        thermal_dens = [getAnalyticalDens(thermal_dens_axis, ρ_pol_rz)]
    else
        error("The 'thermal_profiles_type' input variable was not specified correctly. Currently available options are :FLAT and :DEFAULT. Please correct and re-try.")
    end
    thermal_dist = "" # No distribution. Temp and dens modelled by thermal_temp and thermal_dens
end

# If a TRANSP .cdf file and a TRANSP FI .cdf file have been specified for the thermal species distribution, we have everything we need in the Python thermal_dist variable
if lowercase(fileext_thermal)=="cdf" && lowercase(fileext_FI_cdf)=="cdf" && !projVel
    if !analytic # If not analytic, we know that the 'thermal_dist' variable (defined above) will take care of everything
        verbose && println("Setting all Python variables and structures on all distributed workers/processes... ")
        @everywhere begin
            py"""
            if $verbose:
                print("From Python: Loading TRANSP output from TRANSP files... ")
            tr_out = transp_output.TranspOutput($TRANSP_id, step=1, out_file=$filepath_thermal_distr,fbm_files=[$filepath_FI_cdf]) # Load the TRANSP shot file. Assume first step. This is likely to be patched in the future.
            if $verbose:
                print("From Python: Setting thermal distribution... ")
            thermal_dist = transp_dists.Thermal(tr_out, ion=$thermal_species) # Then load the thermal ion distribution from that .cdf file
            """
        end
        thermal_temp = ""
        thermal_dens = ""
    else
        verbose && println("Input variable 'analytic' set to true, 'filepath_thermal_distr' specified as TRANSP output file ($(filepath_thermal_distr)) and input variable 'filepath_FI_cdf' specified ($(filepath_FI_cdf)). Loading thermal plasma density from TRANSP files...  ")
        # If analytic, we cannot use the 'thermal_dist' variable (want to avoid Python, to optimize speed)
        # Therefore, load the density directly from TRANSP .cdf files
        thermal_temp = [0.0] # Analytic equations assume zero temperature for thermal particle species
        thermal_dens = [getDensProfileFromTRANSP(filepath_FI_cdf, filepath_thermal_distr, thermal_species; verbose=verbose)(ρ_pol_rz)]
        thermal_dist = "" # No distribution. Temp and dens modelled by thermal_temp and thermal_dens
    end
end

# If there is a specified valid timepoint,
# and if a TRANSP .cdf file has been specified, but NOT a TRANSP FI .cdf file, 
if lowercase(fileext_thermal)=="cdf" && !(lowercase(fileext_FI_cdf)=="cdf") && !projVel
    if typeof(timepoint)==String && length(split(timepoint,","))==2
        verbose && println("Input variable 'filepath_thermal_distr' is assumed to be a TRANSP output file ($(filepath_thermal_distr)) but input variable 'filepath_FI_cdf' was not specified. Timepoint ($(timepoint) s) inferred. Loading thermal temperature and density using 'filepath_thermal_distr' and timepoint...")
        thermal_temp = [getTempProfileFromTRANSP(timepoint, filepath_thermal_distr, thermal_species; verbose=verbose)(ρ_pol_rz)]
        thermal_dens = [getDensProfileFromTRANSP(timepoint, filepath_thermal_distr, thermal_species; verbose=verbose)(ρ_pol_rz)]
    else
        @warn "Input variable 'filepath_thermal_distr' specified ($(filepath_thermal_distr)) while input variable 'filepath_FI_cdf' left unspecified. No valid timepoint could be inferred ($(timepoint)). Default OWCF temperature profile (with 5 keV on-axis) and density profile (with 1.0e20 m^-3 on-axis) used instead."
        thermal_temp = [getAnalyticalTemp(5.0, ρ_pol_rz)]
        thermal_dens = [getAnalyticalDens(1.0e20, ρ_pol_rz)]
    end
    thermal_dist = "" # No distribution. Temp and dens modelled by thermal_temp and thermal_dens
end

if analytic
    # Just ensure that thermal temperature is 0, when using analytic equations to computed weight functions 
    # This is purely for aesthetic purposes, since the thermal temperature won't be given as input to the analytic equations forward model
    thermal_temp = [0.0]
end

if !projVel # 'projVel' and 'analytic' cannot both be true at the same time (see initial script safety checks). No risk of conflicts.
    if (thermal_dist=="")
        verbose && println("---> Thermal plasma temperature to be used in weight function computation): $(thermal_temp[1]) keV")
        verbose && println("---> Thermal plasma density to be used in weight function computation): $(thermal_dens[1]) m^-3")
    else
        verbose && println("---> Thermal plasma temperature and density modelled via data in specified TRANSP output files specified as input... ")
    end
else
    verbose && println("---> No thermal plasma distribution/profiles specified. Weight functions to be computed via projected fast-ion velocities.")
    # If we are computing projected velocities, we don't need any thermal plasma distribution. At all.
    thermal_dist = ""
    thermal_temp = ""
    thermal_dens = ""
end

# Transfer the thermal_temp and thermal_dens variables to external processes
@everywhere thermal_temp = $thermal_temp
@everywhere thermal_dens = $thermal_dens
@everywhere thermal_dist = $thermal_dist

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
                        
                        if analytic
                            spec = forward_calc(myVC, [E], clamp.([p],-1,1), [R_of_interest], [z_of_interest], [1.0], Ed_array, B_at_Rz; reaction=reaction, bulk_dens=thermal_dens) # Please see the forward.jl script, for further specifications
                        else
                            py"""
                            # Please note, the '$' symbol is used below to convert objects from Julia to Python. Even Julia PyObjects
                            spec = forwardmodel.calc($([E]), $([p]), $([R_of_interest]), $([z_of_interest]), $([1.0]), thermal_dist, Ed_bin_edges, $B_at_Rz, n_repeat=$gyro_samples, reaction=$reaction, bulk_temp=$thermal_temp, bulk_dens=$thermal_dens, flr=$flr_effects, v_rot=$plasma_rot_at_Rz) # Please see the forward.py script, for further specifications
                            """
                            spec = Vector(py"spec")
                        end

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

                if analytic
                    spec = forward_calc(myVC, [E], clamp.([p],-1,1), [R_of_interest], [z_of_interest], [1.0], Ed_array, B_at_Rz; reaction=reaction, bulk_dens=thermal_dens) # Please see the forward.jl script, for further specifications
                else
                    py"""
                    # Please note, the '$' symbol is used below to convert objects from Julia to Python. Even Julia PyObjects
                    spec = forwardmodel.calc($([E]), $([p]), $([R_of_interest]), $([z_of_interest]), $([1.0]), thermal_dist, Ed_bin_edges, $B_at_Rz, n_repeat=$gyro_samples, reaction=$reaction, bulk_temp=$thermal_temp, bulk_dens=$thermal_dens, flr=$flr_effects, v_rot=$plasma_rot_at_Rz) # Please see the forward.py script, for further specifications
                    """
                    spec = Vector(py"spec")
                end

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

            if analytic
                spec = forward_calc(myVC, [E], clamp.([p],-1,1), [R_of_interest], [z_of_interest], [1.0], Ed_array, B_at_Rz; reaction=reaction, bulk_dens=thermal_dens) # Please see the forward.jl script, for further specifications
            else
                py"""
                # Please note, the '$' symbol is used below to convert objects from Julia to Python. Even Julia PyObjects
                spec = forwardmodel.calc($([E]), $([p]), $([R_of_interest]), $([z_of_interest]), $([1.0]), thermal_dist, Ed_bin_edges, $B_at_Rz, n_repeat=$gyro_samples, reaction=$reaction, bulk_temp=$thermal_temp, bulk_dens=$thermal_dens, flr=$flr_effects, v_rot=$plasma_rot_at_Rz) # Please see the forward.py script, for further specifications
                """
                spec = Vector(py"spec")
            end

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

                if analytic
                    spec = forward_calc(myVC, [E], clamp.([p],-1,1), [R_of_interest], [z_of_interest], [1.0], Ed_array, B_at_Rz; reaction=reaction, bulk_dens=thermal_dens) # Please see the forward.jl script, for further specifications
                else
                    py"""
                    # Please note, the '$' symbol is used below to convert objects from Julia to Python. Even Julia PyObjects
                    spec = forwardmodel.calc($([E]), $([p]), $([R_of_interest]), $([z_of_interest]), $([1.0]), thermal_dist, Ed_bin_edges, $B_at_Rz, n_repeat=$gyro_samples, reaction=$reaction, bulk_temp=$thermal_temp, bulk_dens=$thermal_dens, flr=$flr_effects, v_rot=$plasma_rot_at_Rz) # Please see the forward.py script, for further specifications
                    """
                    spec = Vector(py"spec")
                end

                Wtot[:,iE,ip] = spec
            end
        end
    end

    if saveVparaVperpWeights
        verbose && println("Transforming (E,p) weights to (vpara,vperp)... ")
        nvpara = Int64(round(sqrt(nE*np))) # We know the default values for Ep2VparaVperp()
        nvperp = nvpara # We know the default values for Ep2VparaVperp()
        W_vel = zeros(nEd,nvpara,nvperp)
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
        verbose && println("Saving weight function matrix... ")
        if iiimax==1 # If you intend to calculate only one weight matrix
            global filepath_output_orig = folderpath_o*"velWeights_"*tokamak*"_"*TRANSP_id*"_at"*timepoint*"s_"*diagnostic_name*"_"*pretty2scpok(reaction; projVel = projVel)*"_$(length(range(Ed_min,stop=Ed_max,step=Ed_diff))-1)x$(nE)x$(np)"
        else # If you intend to calculate several (identical) weight matrices
            global filepath_output_orig = folderpath_o*"velWeights_"*tokamak*"_"*TRANSP_id*"_at"*timepoint*"s_"*diagnostic_name*"_"*pretty2scpok(reaction; projVel = projVel)*"_$(iii)"
        end
        global filepath_output = deepcopy(filepath_output_orig)
        global count = 1
        while isfile(filepath_output*".jld2") # To take care of not overwriting files. Add _(1), _(2) etc
            global filepath_output = filepath_output_orig*"_($(Int64(count)))"
            global count += 1 # global scope, to surpress warnings
        end
        global filepath_output = filepath_output*".jld2"
        list_o_saved_filepaths[iii] = filepath_output # Store a record of the saved file
        myfile_s = jldopen(filepath_output,true,true,false,IOStream)
        write(myfile_s, "W", Wtot)
        write(myfile_s, "E_array", vec(E_array))
        write(myfile_s, "p_array", vec(p_array))
        if instrumental_response
            write(myfile_s, "Ed_array", instrumental_response_output)
        else
            write(myfile_s, "Ed_array", Ed_array)
        end
        write(myfile_s, "reaction", reaction)
        write(myfile_s, "R", R_of_interest)
        write(myfile_s, "z", z_of_interest)
        if projVel
            write(myfile_s, "projVel", projVel)
        end
        if saveVparaVperpWeights
            write(myfile_s, "W_vel", W_vel)
            write(myfile_s, "vpara_array", vpara_array)
            write(myfile_s, "vperp_array", vperp_array)
        end
        if instrumental_response
            write(myfile_s, "W_raw", Wtot_raw)
            write(myfile_s, "Ed_array_raw",Ed_array)
            write(myfile_s, "instrumental_response_input", instrumental_response_input)
            write(myfile_s, "instrumental_response_output", instrumental_response_output)
            write(myfile_s, "instrumental_response_matrix", instrumental_response_matrix)
            if saveVparaVperpWeights
                write(myfile_s, "W_vel_raw", W_vel_raw)
            end
        end
        if plasma_rot
            write(myfile_s, "plasma_rot_at_Rz", plasma_rot_at_Rz)
        end
        write(myfile_s, "filepath_thermal_distr", filepath_thermal_distr)
        write(myfile_s, "filepath_start", filepath_start)
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
        R = myfile["R"]
        z = myfile["z"]
        if projVel
            projVel = myfile["projVel"]
        end
        if saveVparaVperpWeights
            W_vel_total = zeros(size(myfile["W_vel"]))
            vpara_array = myfile["vpara_array"]
            vperp_array = myfile["vperp_array"]
        end
        if instrumental_response
            W_raw_total = zeros(size(myfile["W_raw"]))
            Ed_array_raw = myfile["Ed_array_raw"]
            instrumental_response_input = myfile["instrumental_response_input"]
            instrumental_response_output = myfile["instrumental_response_output"]
            instrumental_response_matrix = myfile["instrumental_response_matrix"]
            if saveVparaVperpWeights
                W_vel_raw_total = zeros(size(myfile["W_vel_raw"]))
            end
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
        myfile = jldopen(folderpath_o*"velWeights_"*tokamak*"_"*TRANSP_id*"_at"*timepoint*"s_"*diagnostic_name*"_"*pretty2scpok(reaction; projVel = projVel)*".jld2",true,true,false,IOStream)
        write(myfile,"W",W_total)
        write(myfile,"E_array",E_array)
        write(myfile,"p_array",p_array)
        write(myfile,"Ed_array",Ed_array)
        write(myfile,"reaction",reaction)
        write(myfile,"R",R)
        write(myfile,"z",z)
        if projVel
            write(myfile,"projVel",projVel)
        end
        if saveVparaVperpWeights
            write(myfile,"W_vel",W_vel_total)
            write(myfile,"vpara_array",vpara_array)
            write(myfile,"vperp_array",vperp_array)
        end
        if instrumental_response
            write(myfile,"W_raw",W_raw_total)
            write(myfile,"Ed_array_raw",Ed_array_raw)
            write(myfile,"instrumental_response_input",instrumental_response_input)
            write(myfile,"instrumental_response_output",instrumental_response_output)
            write(myfile,"instrumental_response_matrix",instrumental_response_matrix)
            if saveVparaVperpWeights
                write(myfile,"W_vel_raw",W_vel_raw_total)
            end
        end
        if plasma_rot
            write(myfile, "plasma_rot_at_Rz", plasma_rot_at_Rz)
        end
        write(myfile,"filepath_thermal_distr",filepath_thermal_distr)
        write(myfile,"filepath_start",filepath_start)
        close(myfile)
    end
end
println("------ Done! ------")
