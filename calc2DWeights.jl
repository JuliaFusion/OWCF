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
#   reaction_full - The nuclear fusion reaction for which the weights are computed - String
#   filepath_thermal_distr - The filepath of the thermal species distribution. For reference - String
# If analytical weight functions (with projected velocities) are computed, the saved file will also have the key
#   analytical2DWs - True if analytical weight functions were computed. False otherwise - Bool
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

# Script written by Henrik Järleblad. Last maintained 2023-12-18.
################################################################################################

## ---------------------------------------------------------------------------------------------
verbose && println("Loading Julia packages... ")
@everywhere begin
    cd(folderpath_OWCF) # Necessary to move all the workers to the correct folder
    using PyCall # For using Python code in Julia
    using EFIT # For calculating magn½etic equilibrium quantities
    using Equilibrium # For loading flux function data, tokamak geometry data etc.
    using GuidingCenterOrbits # For calculating guiding-center orbits
    using OrbitTomography # This is what all this is about!
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

## ---------------------------------------------------------------------------------------------
if !analytical2DWs
    reaction_full = deepcopy(reaction) # Make a fully independent copy of the fusion reaction variable
    reaction = full2reactsOnly(reaction; verbose=verbose, projVelocity=analytical2DWs) # Converts from 'a(b,c)d' format to 'a-b' format (reactants only)
else
    reaction_full = deepcopy(reaction) # Make a fully independent copy of the fusion reaction variable
end
#@everywhere reaction_full = $reaction_full # Not yet necessary. Might be necessary when starting to compute two-step fusion reactions with the OWCF via the DRESS code
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
    println("Thermal distribution file format: ."*fileext_thermal)
    error("Unknown thermal distribution file format. Please re-specify file and re-try.")
end

# REWRITE THESE ERROR CHECKS TO TAKE THE DIFFERENT FILE EXTENSIONS INTO ACCOUNT
if fileext_thermal=="cdf" && !(fileext_FI_cdf=="cdf")
    error("filepath_thermal_distr specified as TRANSP .cdf file, but filepath_FI_cdf was wrongly specified. Please specify and re-try.")
end
## ---------------------------------------------------------------------------------------------
# Safety check for analytical orbit weight function computations
if analytical2DWs
    !(split(reaction,"-")[1] == "proj") && error("Analytical orbit weight function computation was specified, but input variable 'reaction' was not correctly specified (it should be specified as 'proj-X' where 'X' is the fast-ion species). Please correct and re-try.")
end
if !analytical2DWs
    (split(reaction,"-")[1] == "proj") && error("Normal orbit weight function computation was specified, but input variable 'reaction' was not correctly specified (it should be specified as 'a(b,c)d' where a is thermal ion, b is fast ion, c is emitted particle and d is the product nucleus. Please correct and re-try.")
end

## ---------------------------------------------------------------------------------------------
# Determine fast-ion and thermal (thermal) species from inputs in start file
if !analytical2DWs
    thermal_species, FI_species = checkReaction(reaction_full) # Check the specified fusion reaction, and extract thermal and fast-ion species
else
    FI_species = split(reaction,"-")[2] # Assumed 'proj-X' format
    thermal_species = split(reaction,"-")[1] # Will just be 'proj', assumed 'proj-X' format
end
@everywhere thermal_species = $thermal_species # Transfer variable to all external processes
@everywhere FI_species = $FI_species # Transfer variable to all external processes

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

    @everywhere begin
        verbose && println("Loading the Interpolations package (needed because .jld2 thermal file was specified)... ")
        using Interpolations
    end
end
# Loading tokamak information and TRANSP RUN-ID from thermal .cdf file
if fileext_thermal=="cdf"
    verbose && println("Loading TRANSP id information... ")
    TRANSP_id = (split((split(filepath_thermal_distr,"."))[1],"/"))[end] # Assume part between last '/' and '.' is the TRANSP id
else
    TRANSP_id = ""
end


if filepath_thermal_distr=="" && (!(typeof(thermal_temp_axis)==Float64) && !(typeof(thermal_temp_axis)==Int64)) && !analytical2DWs
    @everywhere thermal_temp_axis = 3.0
    @warn "filepath_thermal_distr was not specified, and thermal_temp_axis was not specified correctly. thermal_temp_axis will be set to default value of 3.0 keV."
end
if filepath_thermal_distr=="" && !(typeof(thermal_dens_axis)==Float64) && !analytical2DWs
    @everywhere thermal_dens_axis = 1.0e20
    @warn "filepath_thermal_distr was not specified, and thermal_dens_axis was not specified correctly. thermal_dens_axis will be set to default value of 1.0e20 m^-3."
end

@everywhere TRANSP_id = $TRANSP_id
@everywhere tokamak = $tokamak

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
# Pre-processing (R,z) of interest input data
if (R_of_interest==:r_mag)
    R_of_interest, zdummy = magnetic_axis(M)
end
if (z_of_interest==:z_mag)
    rdummy, z_of_interest = magnetic_axis(M)
end

## ---------------------------------------------------------------------------------------------
# If available, load instrumental response and process it accordingly
global instrumental_response # Declare global scope
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
println("----------------------------------------calc2DWeights.jl----------------------------------------")
if !(tokamak=="")
    print("Tokamak specified: "*tokamak*"      ")
else
    print("Tokamak not specified.         ")
end
if !(TRANSP_id=="")
    println("TRANSP ID specified: "*TRANSP_id)
else
    println("TRANSP ID not specified.")
end
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
if !(filepath_thermal_distr=="")
    println("Thermal distribution file specified: "*filepath_thermal_distr)
elseif (filepath_thermal_distr=="") && !analytical2DWs
    println("Thermal distribution file not specified.")
    println("Thermal ion ($((split(reaction,"-"))[1])) temperature on axis will be set to $(thermal_temp_axis) keV.")
    println("Thermal ion ($((split(reaction,"-"))[1])) density on axis will be set to $(thermal_dens_axis) m^-3.")
else
    println("Analytical 2D weight functions will be computed using the projected velocity (u) of the fast ions.")
end
println("Magnetic equilibrium file specified: "*filepath_equil)
println("")
if !analytical2DWs
    println("Fusion reaction specified: "*reaction_full)
else
    println("Projected velocity (u) will be used as weights for the weight functions.")
end
println("Fast-ion species specified: "*FI_species)
if emittedParticleHasCharge && !analytical2DWs
    println("The emitted "*getEmittedParticle(reaction_full)*" particle of the "*reaction_full*" reaction has non-zero charge!")
    println("The resulting energy distribution for "*getEmittedParticle(reaction_full)*" from the plasma as a whole will be computed.")
end
println("")
if distributed
    println("Parallel computing will be used with $(nprocs()) processes (1 main + $(nprocs()-1) workers).")
else
    println("Single-threaded computing the weights... Good luck!")
end
println("")
println("$(iiimax) weight matrix/matrices will be computed.")
println("")
println("2D weight function(s) will be computed for a $(length(E_array))x$(length(p_array)) (E,p)-space grid with")
println("Energy: [$(minimum(E_array)),$(maximum(E_array))] keV")
println("Pitch: [$(minimum(p_array)),$(maximum(p_array))]")
println("The (R,z) point of interest is: ($(R_of_interest),$(z_of_interest)) m")
println("")
if !analytical2DWs
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
    println(folderpath_o*"velWeights_"*tokamak*"_"*TRANSP_id*"_at"*timepoint*"s_"*diagnostic_name*"_"*pretty2scpok(reaction_full; projVel = analytical2DWs)*"_$(length(range(Ed_min,stop=Ed_max,step=Ed_diff))-1)x$(nE)x$(np).jld2")
else
    println(folderpath_o*"velWeights_"*tokamak*"_"*TRANSP_id*"_at"*timepoint*"s_"*diagnostic_name*"_"*pretty2scpok(reaction_full; projVel = analytical2DWs)*"_1.jld2")
    println("... ")
    println(folderpath_o*"velWeights_"*tokamak*"_"*TRANSP_id*"_at"*timepoint*"s_"*diagnostic_name*"_"*pretty2scpok(reaction_full; projVel = analytical2DWs)*"_$(iiimax).jld2")
    if iii_average
        println("---> Average of all files will be computed and saved to: ")
        println("---> "*folderpath_o*"velWeights_"*tokamak*"_"*TRANSP_id*"_at"*timepoint*"s_"*diagnostic_name*"_"*pretty2scpok(reaction_full; projVel = analytical2DWs)*".jld2")
    end
end
if debug
    println("")
    println("!!!!!! DEBUGGING SPECIFIED. ALGORITHM WILL DEBUG. !!!!!!")
    println("")
end
println("Please remove previously saved files with the same file name (if any) prior to script completion. Quickly!")
println("")
println("If you would like to change any settings, please edit the start_calc2DW_template.jl file or similar.")
println("")
println("Written by Henrik Järleblad. Last maintained 2023-11-18.")
println("--------------------------------------------------------------------------------------------------")
println("")

## ---------------------------------------------------------------------------------------------
# Python code (essentially calling the forward model black box)
# This is how you write Python code in Julia: py""" [Python code] """
verbose && println("Loading Python modules... ")
@everywhere begin
    py"""
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
## If no thermal distribution has been specified, we are going to need the OWCF default temp. and dens. profiles
if filepath_thermal_distr==""
    @everywhere begin
        include("misc/temp_n_dens.jl")
    end
end

## ---------------------------------------------------------------------------------------------
# Setting Python variables and structures on all distributed workers/processes...
verbose && println("Setting all Python variables and structures on all distributed workers/processes... ")
@everywhere begin
    py"""
    # The '$' in front of many Python variables means that the variable is defined in Julia, not in Python.
    reaction = $reaction
    thermal_species = $thermal_species
    analytical2DWs = $analytical2DWs
    if $verbose:
        print("From Python: Loading forward model with diagnostic... ") 
    forwardmodel = forward.Forward($diagnostic_filepath) # Pre-initialize the forward model

    # Load TRANSP simulation data
    if (not $TRANSP_id=="") and (not analytical2DWs): # If there is some TRANSP_id specified and we do not want to simply compute projected velocities...
        if ($fileext_thermal).lower()=="cdf": # If there is some TRANSP .cdf output file specified...
            if $verbose:
                print("From Python: Loading TRANSP output from TRANSP files... ")
            tr_out = transp_output.TranspOutput($TRANSP_id, step=1, out_file=$filepath_thermal_distr,fbm_files=[$filepath_FI_cdf]) # Load the TRANSP shot file. Assume first step. This is likely to be patched in the future.
            if $verbose:
                print("From Python: Setting thermal distribution... ")
            thermal_dist = transp_dists.Thermal(tr_out, ion=thermal_species) # Then load the thermal ion distribution from that .cdf file
        else:
            raise ValueError('From Python: TRANSP_id was specified, but filepath_thermal_distr was not (this should be impossible). Please correct and re-try.')
    else:
        thermal_dist = "" # Otherwise, just let the thermal_dist variable be the empty string

    Ed_bin_edges = np.arange($Ed_min,$Ed_max,$Ed_diff) # diagnostic spectrum bin edges (keV or m/s)
    if len(Ed_bin_edges)==1: # Make sure that there are at least one lower and one upper bin edge
        dEd = (($Ed_max)-($Ed_min))/10
        Ed_bin_edges = np.arange($Ed_min,($Ed_max)+dEd,$Ed_diff)
    Ed_vals = 0.5*(Ed_bin_edges[1:] + Ed_bin_edges[:-1]) # bin centers (keV or m/s)
    nEd = len(Ed_vals)
    """
end
nEd = py"nEd"
global Ed_array = vec(py"Ed_vals")

## ---------------------------------------------------------------------------------------------
# Pre-processing thermal density and temperature data
if !isnothing(thermal_dens) && (typeof(thermal_dens)==Float64 || typeof(thermal_dens)==Int64)
    thermal_dens = [thermal_dens]
else
    if fileext_thermal=="jld2"
        verbose && println("Creating thermal density interpolations object... ")
        thermal_dens_itp = Interpolations.interpolate((ρ_pol_array,), thermal_dens_array, Gridded(Linear()))
        thermal_dens_etp = Interpolations.extrapolate(thermal_dens_itp,0) # Add what happens in extrapolation scenario. 0 means we assume vacuum scrape-off layer (SOL)
        ψ_rz = M(R_of_interest,z_of_interest)
        psi_on_axis, psi_at_bdry = psi_limits(M)
        ρ_pol_rz = sqrt((ψ_rz-psi_on_axis)/(psi_at_bdry-psi_on_axis)) # The formula for the normalized flux coordinate ρ_pol = sqrt((ψ-ψ_axis)/(ψ_edge-ψ_axis))
        thermal_dens = [thermal_dens_etp(ρ_pol_rz)]
    end

    # If a thermal species distribution has not been specified, we simply need the thermal species temperature on axis
    if !isfile(filepath_thermal_distr)
        ψ_rz = M(R_of_interest,z_of_interest)
        psi_on_axis, psi_at_bdry = psi_limits(M)
        ρ_pol_rz = sqrt((ψ_rz-psi_on_axis)/(psi_at_bdry-psi_on_axis)) # The formula for the normalized flux coordinate ρ_pol = sqrt((ψ-ψ_axis)/(ψ_edge-ψ_axis))
        thermal_dens = [getAnalyticalDens(thermal_dens_axis,ρ_pol_rz)] # Use the default OWCF density profile
    end

    # If a TRANSP .cdf file has been specified for the thermal species distribution, we have everything we need in the Python thermal_dist variable
    if lowercase(fileext_thermal)=="cdf"
        thermal_dens = ""
    end

    # Transfer the thermal_dens variables to external processes
    @everywhere thermal_dens = $thermal_dens
end

if !isnothing(thermal_temp) && (typeof(thermal_temp)==Float64 || typeof(thermal_temp)==Int64)
    thermal_temp = [thermal_temp]
else
    if fileext_thermal=="jld2"
        verbose && println("Creating thermal temperature interpolations object... ")
        thermal_temp_itp = Interpolations.interpolate((ρ_pol_array,), thermal_temp_array, Gridded(Linear()))
        thermal_temp_etp = Interpolations.extrapolate(thermal_temp_itp,0) # Add what happens in extrapolation scenario. 0 means we assume vacuum scrape-off layer (SOL)
        ψ_rz = M(R_of_interest,z_of_interest)
        psi_on_axis, psi_at_bdry = psi_limits(M)
        ρ_pol_rz = sqrt((ψ_rz-psi_on_axis)/(psi_at_bdry-psi_on_axis)) # The formula for the normalized flux coordinate ρ_pol = sqrt((ψ-ψ_axis)/(ψ_edge-ψ_axis))
        thermal_temp = [thermal_temp_etp(ρ_pol_rz)]
    end

    # If a thermal species distribution has not been specified, we simply need the thermal species temperature on axis
    if !isfile(filepath_thermal_distr)
        ψ_rz = M(R_of_interest,z_of_interest)
        psi_on_axis, psi_at_bdry = psi_limits(M)
        ρ_pol_rz = sqrt((ψ_rz-psi_on_axis)/(psi_at_bdry-psi_on_axis)) # The formula for the normalized flux coordinate ρ_pol = sqrt((ψ-ψ_axis)/(ψ_edge-ψ_axis))
        thermal_temp = [getAnalyticalDens(thermal_temp_axis,ρ_pol_rz)] # Use the default OWCF temperature profile
    end

    # If a TRANSP .cdf file has been specified for the thermal species distribution, we have everything we need in the Python thermal_dist variable
    if lowercase(fileext_thermal)=="cdf"
        thermal_temp = ""
    end

    # Transfer the thermal_temp variables to external processes
    @everywhere thermal_temp = $thermal_temp
end

## ---------------------------------------------------------------------------------------------
E_iE_p_ip_array = zip(repeat(E_array,outer=length(p_array)),repeat(eachindex(E_array),outer=length(p_array)),repeat(p_array,inner=length(E_array)),repeat(eachindex(p_array),inner=length(E_array)))
B_at_Rz = Equilibrium.Bfield(M,R_of_interest,z_of_interest) # Calculate the B-field vector at the (R,z) point
B_at_Rz = [B_at_Rz[1], B_at_Rz[2], B_at_Rz[3]] # Re-structure for Python
B_at_Rz = reshape(B_at_Rz,(length(B_at_Rz),1)) # Re-structure for Python
@everywhere B_at_Rz = $B_at_Rz
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

                        py"""
                        #forwardmodel = $forwardmodel # Convert the Forward Python object from Julia to Python.
                        # Please note, the '$' symbol is used below to convert objects from Julia to Python. Even Julia PyObjects
                        spec = forwardmodel.calc($([E]), $([p]), $([R_of_interest]), $([z_of_interest]), $([1.0]), thermal_dist, Ed_bin_edges, $B_at_Rz, n_repeat=$gyro_samples, reaction=$reaction, bulk_temp=$thermal_temp, bulk_dens=$thermal_dens) # Please see the forward.py script, for further specifications
                        """

                        spec = vec(py"spec")
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

                py"""
                #forwardmodel = $forwardmodel # Convert the Forward Python object from Julia to Python.
                # Please note, the '$' symbol is used below to convert objects from Julia to Python. Even Julia PyObjects
                spec = forwardmodel.calc($([E]), $([p]), $([R_of_interest]), $([z_of_interest]), $([1.0]), thermal_dist, Ed_bin_edges, $B_at_Rz, n_repeat=$gyro_samples, reaction=$reaction, bulk_temp=$thermal_temp, bulk_dens=$thermal_dens) # Please see the forward.py script, for further specifications
                """

                spec = vec(py"spec")
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

            py"""
            #forwardmodel = $forwardmodel # Convert the Forward Python object from Julia to Python.
            # Please note, the '$' symbol is used below to convert objects from Julia to Python. Even Julia PyObjects
            spec = forwardmodel.calc($([E]), $([p]), $([R_of_interest]), $([z_of_interest]), $([1.0]), thermal_dist, Ed_bin_edges, $B_at_Rz, n_repeat=$gyro_samples, reaction=$reaction, bulk_temp=$thermal_temp, bulk_dens=$thermal_dens) # Please see the forward.py script, for further specifications
            """

            spec = vec(py"spec")
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

                py"""
                #forwardmodel = $forwardmodel # Convert the Forward Python object from Julia to Python.
                # Please note, the '$' symbol is used below to convert objects from Julia to Python. Even Julia PyObjects
                spec = forwardmodel.calc($([E]), $([p]), $([R_of_interest]), $([z_of_interest]), $([1.0]), thermal_dist, Ed_bin_edges, $B_at_Rz, n_repeat=$gyro_samples, reaction=$reaction, bulk_temp=$thermal_temp, bulk_dens=$thermal_dens) # Please see the forward.py script, for further specifications
                """

                spec = vec(py"spec")
                Wtot[:,iE,ip] = spec
            end
        end
    end

    if saveVparaVperpWeights
        verbose && println("Transforming (E,p) weights to (vpara,vperp)... ")
        nvpara = Int64(round(sqrt(nE*np)))
        nvperp = nvpara
        W_vel = zeros(nEd,nvpara,nvperp)
        vpara_array, vperp_array, W_vel[1,:,:] = Ep2VparaVperp(E_array, p_array, Wtot[1,:,:])
        for iEd=2:nEd
            vpara_array, vperp_array, W_vel[iEd,:,:] = Ep2VparaVperp(E_array, p_array, Wtot[iEd,:,:])
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
            global filepath_output_orig = folderpath_o*"velWeights_"*tokamak*"_"*TRANSP_id*"_at"*timepoint*"s_"*diagnostic_name*"_"*pretty2scpok(reaction_full; projVel = analytical2DWs)*"_$(length(range(Ed_min,stop=Ed_max,step=Ed_diff))-1)x$(nE)x$(np)"
        else # If you intend to calculate several (identical) weight matrices
            global filepath_output_orig = folderpath_o*"velWeights_"*tokamak*"_"*TRANSP_id*"_at"*timepoint*"s_"*diagnostic_name*"_"*pretty2scpok(reaction_full; projVel = analytical2DWs)*"_$(iii)"
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
        write(myfile_s, "reaction_full", reaction_full)
        write(myfile_s, "R", R_of_interest)
        write(myfile_s, "z", z_of_interest)
        if analytical2DWs
            write(myfile_s, "analytical2DWs", analytical2DWs)
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
        write(myfile_s, "filepath_thermal_distr", filepath_thermal_distr)
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
        reaction_full = myfile["reaction_full"]
        R = myfile["R"]
        z = myfile["z"]
        if analytical2DWs
            analytical2DWs = myfile["analytical2DWs"]
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
        filepath_thermal_distr = myfile["filepath_thermal_distr"]
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
        myfile = jldopen(folderpath_o*"velWeights_"*tokamak*"_"*TRANSP_id*"_at"*timepoint*"s_"*diagnostic_name*"_"*pretty2scpok(reaction_full; projVel = analytical2DWs)*".jld2",true,true,false,IOStream)
        write(myfile,"W",W_total)
        write(myfile,"E_array",E_array)
        write(myfile,"p_array",p_array)
        write(myfile,"Ed_array",Ed_array)
        write(myfile,"reaction_full",reaction_full)
        write(myfile,"R",R)
        write(myfile,"z",z)
        if analytical2DWs
            write(myfile,"analytical2DWs",analytical2DWs)
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
        write(myfile,"filepath_thermal_distr",filepath_thermal_distr)
        close(myfile)
    end
end
println("------ Done! ------")
