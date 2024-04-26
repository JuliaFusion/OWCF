#########################################  calcOrbWeights.jl #########################################

#### Description:
# This script computes orbit weight functions, it's as simple as that. It can be considered the flagship
# of the OWCF. As the rest of the OWCF, it utilizes Luke Stagner's orbit computation codes to compute
# the guiding-center orbits of fast ions in toroidally symmetric fusion devices. The orbit computation
# codes are (as of 2022-08-26) packaged into the Julia packages GuidingCenterOrbits.jl and OrbitTomography.jl.
# Furthermore, in the current version of the OWCF, the calcOrbWeights.jl script is taylored to utilize
# the DRESS code (J. Eriksson et al, CPC, 199, 40-46, 2016) to compute the orbit weight functions.
# The calcOrbWeights.jl script does this by computing weighted (E,p,R,z) points for each orbit, and then
# sends them into the DRESS code which then returns an expected signal for that (E,pm,Rm) grid point.
# In future versions, it can be easily modified to instead save the weighted (E,p,R,z) points of the
# orbits in a file readable by e.g. FIDASIM. A direct call to FIDASIM or other codes can also be easily
# implemented by modifying the calcOrbSpec.jl helper script, on which calcOrbWeights.jl depends.
# 
# The DRESS code is written in Python and calcOrbWeights.jl utilizes the DRESS code via the Julia-Python
# interface package PyCall.jl
#
# The calcOrbWeights.jl script computes orbit weight functions for an equidistant, rectangular grid 
# in (E,pm,Rm) orbit space. Irregular/non-equidistant grid points are currently not supported.
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
# Please see the start_calcOW_template.jl file for further input information.

#### Inputs (Units given when defined in script)
# Given via input file start_calcOW_template.jl, for example. 'template' should be replaced by whatever.

#### Outputs
# -

#### Saved files
# orbWeights_[tokamak]_[TRANSP_id]_at[timepoint]s_[diagnostic_name]_[nE]x[npm]x[nRm].jld2 - If iiimax == 1
# orbWeights_[tokamak]_[TRANSP_id]_at[timepoint]s_[diagnostic_name]_$(i).jld2 - If iiimax != 1
# Regardless of saved file name, this saved file will have the fields:
#   W - The computed orbit weights. Dimensions are (channels, valid orbits) - Array{Float64,2}
#   E_array - The fast-ion energy grid array used for orbit space - Array{Float64,1}
#   pm_array - The fast-ion pm grid array used for orbit space - Array{Float64,1}
#   Rm_array - The fast-ion Rm grid array used for orbit space - Array{Float64,1}
#   Ed_array - The diagnostic energy bin centers - Array{Float64,1}
#   reaction_full - The nuclear fusion reaction for which the orbit weights are computed - String
#   filepath_thermal_distr - The filepath of the thermal species distribution. For reference - String
# If analytical orbit weight functions (with projected velocities) are computed, the saved file will also have the key
#   analyticalOWs - True if analytical orbit weight functions were computed. False otherwise - Bool
# If an orbit-space grid file was used to define the orbit grid for the orbit weight functions, the saved file will also have the key
#   og_filepath - The path to the .jld2-file (should be an output of calcOrbGrid.jl) used as orbit-space grid input - String

### Other
# Please note that the diagnostic energy grid will be created as bin centers.
# That is, the first diagnostic energy grid value will be (Ed_min+Ed_diff/2) and so on.

# Script written by Henrik Järleblad. Last maintained 2022-10-05.
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
    include("misc/temp_n_dens.jl")
    include("misc/availReacts.jl") # To examine fusion reaction and extract thermal and fast-ion species
    include("misc/rewriteReacts.jl") # To rewrite a fusion reaction from the A(b,c)D format to the A-b=c-D format
    pushfirst!(PyVector(pyimport("sys")."path"), "") # To add the forward, transp_dists, transp_output and vcone modules (scripts in current path)
end

## ---------------------------------------------------------------------------------------------
if !analyticalOWs
    reaction_full = deepcopy(reaction) # Make a fully independent copy of the fusion reaction variable
    reaction = full2reactsOnly(reaction; verbose=verbose, projVelocity=analyticalOWs) # Converts from 'a(b,c)d' format to 'a-b' format (reactants only)
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
if analyticalOWs
    !(split(reaction,"-")[1] == "proj") && error("Analytical orbit weight function computation was specified, but input variable 'reaction' was not correctly specified (it should be specified as 'proj-X' where 'X' is the fast-ion species). Please correct and re-try.")
end
if !analyticalOWs
    (split(reaction,"-")[1] == "proj") && error("Normal orbit weight function computation was specified, but input variable 'reaction' was not correctly specified (it should be specified as 'a(b,c)d' where a is thermal ion, b is fast ion, c is emitted particle and d is the product nucleus. Please correct and re-try.")
end

## ---------------------------------------------------------------------------------------------
# Determine fast-ion and thermal (thermal) species from inputs in start file
if !analyticalOWs
    thermal_species, FI_species = checkReaction(reaction_full) # Check the specified fusion reaction, and extract thermal and fast-ion species
else
    FI_species = split(reaction,"-")[2] # Assumed 'proj-X' format
    thermal_species = split(reaction,"-")[1] # Will just be 'proj', assumed 'proj-X' format
end
@everywhere thermal_species = $thermal_species # Transfer variable to all external processes
@everywhere FI_species = $FI_species # Transfer variable to all external processes

## ---------------------------------------------------------------------------------------------
# Loading tokamak information and TRANSP RUN-ID
if fileext_thermal=="jld2"
    verbose && println("Loading thermal .jld2 file data... ")
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


if filepath_thermal_distr=="" && (!(typeof(thermal_temp_axis)==Float64) && !(typeof(thermal_temp_axis)==Int64)) && !analyticalOWs
    @everywhere thermal_temp_axis = 3.0
    @warn "filepath_thermal_distr was not specified, and thermal_temp_axis was not specified correctly. thermal_temp_axis will be set to default value of 3.0 keV."
end
if filepath_thermal_distr=="" && !(typeof(thermal_dens_axis)==Float64) && !analyticalOWs
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
        timepoint = "TIMELESS" # Unknown timepoint for magnetic equilibrium
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
# Defining orbit grid vectors
if !(og_filepath===nothing)
    verbose && println("Filepath to .jld2 file containing orbit grid was specified. Loading orbit grid... ")
    myfile = jldopen(og_filepath,false,false,false,IOStream)
    og = myfile["og"]
    og_orbs = myfile["og_orbs"]
    FI_species_loaded = myfile["FI_species"]
    extra_kw_args = myfile["extra_kw_args"]
    close(myfile)
    E_array = og.energy
    pm_array = og.pitch
    Rm_array = og.r
    og = nothing # Clear memory, to minimize memory usage
    if !(lowercase(FI_species)==lowercase(FI_species_loaded))
        @warn "Fast-ion species ("*FI_species_loaded*") loaded from orbit-grid .jld2 file does not match fast-ion species in specified fusion reaction ("*FI_species*"). Did you confuse the thermal species ("*thermal_species*") for the fast-ion species ("*FI_species*")?"
    end
    if !(lowercase(FI_species)==lowercase(FI_species_loaded)) && !(lowercase(thermal_species)==lowercase(FI_species_loaded))
        error("Fast-ion species ("*FI_species_loaded*") loaded from orbit-grid .jld2 file does not match any reactants in specified fusion reaction ("*reaction*"). Please correct and re-try.")
    end
else
    verbose && println("Defining orbit grid vectors... ")
    if !(E_array == nothing)
        Emin = minimum(E_array)
        Emax = maximum(E_array)
        nE = length(E_array)
    else
        E_array = vec(range(Emin,stop=Emax,length=nE))
    end
    if !(pm_array == nothing)
        pm_min = minimum(pm_array)
        pm_max = maximum(pm_array)
        npm = length(pm_array)
    else
        pm_array = vec(range(pm_min, stop=pm_max, length=npm))
    end
    if !(Rm_array == nothing)
        Rm_min = minimum(Rm_array)
        Rm_max = maximum(Rm_array)
        nRm = length(Rm_array)
    else
        if (Rm_min==nothing) || (Rm_max==nothing)
            if inclPrideRockOrbs
                # 4/5 of the distance from the HFS wall to the magnetic axis is usually enough to capture all the Pride Rock orbits
                Rm_array = vec(range((4*magnetic_axis(M)[1]+minimum(wall.r))/5, stop=maximum(wall.r), length=nRm))
            else
                Rm_array = vec(range(magnetic_axis(M)[1], stop=maximum(wall.r), length=nRm))
            end
        else
            Rm_array = vec(range(Rm_min, stop=Rm_max, length=nRm))
        end
    end
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
@everywhere instrumental_response = $instrumental_response

## ---------------------------------------------------------------------------------------------
# Printing script info and inputs
println("")
println("----------------------------------------calcOrbWeights.jl----------------------------------------")
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
    println("instrumental_response_filepath specified. Diagnostic response included.")
else
    println("instrumental_response_filepath not specified. Diagnostic response not included.")
end
if !(filepath_thermal_distr=="")
    println("Thermal distribution file specified: "*filepath_thermal_distr)
elseif (filepath_thermal_distr=="") && !analyticalOWs
    println("Thermal distribution file not specified.")
    println("Thermal ion ($((split(reaction,"-"))[1])) temperature on axis will be set to $(thermal_temp_axis) keV.")
    println("Thermal ion ($((split(reaction,"-"))[1])) density on axis will be set to $(thermal_dens_axis) m^-3.")
else
    println("Analytical orbit weight functions will be computed using the projected velocity (u) of the fast ions.")
end
println("Magnetic equilibrium file specified: "*filepath_equil)
println("")
if !analyticalOWs
    println("Fusion reaction specified: "*reaction_full)
else
    println("Projected velocity (u) will be used as weights for the weight functions.")
end
println("Fast-ion species specified: "*FI_species)
if emittedParticleHasCharge && !analyticalOWs
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
println("Orbit weight function(s) will be computed for a $(length(E_array))x$(length(pm_array))x$(length(Rm_array)) orbit-space grid with")
println("Energy: [$(minimum(E_array)),$(maximum(E_array))] keV")
println("Pitch maximum: [$(minimum(pm_array)),$(maximum(pm_array))]")
println("Radius maximum: [$(minimum(Rm_array)),$(maximum(Rm_array))] m")
println("")
if include2Dto4D
    println("Orbit weight matrix will be inflated and saved in its full 4D format upon completion of computations.")
end
if !analyticalOWs
    println("There will be $(length(range(Ed_min,stop=Ed_max,step=Ed_diff))-1) diagnostic energy bin(s) with")
    println("Lower diagnostic energy bound: $(Ed_min) keV")
    println("Upper diagnostic energy bound: $(Ed_max) keV")
else
    println("There will be $(length(range(Ed_min,stop=Ed_max,step=Ed_diff))-1) projected velocity (u) bin(s) with")
    println("Lower u bound: $(Ed_min) m/s")
    println("Upper u bound: $(Ed_max) m/s")
end
println("")
println("The orbit integration algorithm will used the extra keyword arguments: ")
println(extra_kw_args)
println("")
println("Results will be saved to: ")
if iiimax == 1
    println(folderpath_o*"orbWeights_"*tokamak*"_"*TRANSP_id*"_at"*timepoint*"s_"*diagnostic_name*"_"*pretty2scpok(reaction_full; projVel = analyticalOWs)*"_$(length(range(Ed_min,stop=Ed_max,step=Ed_diff))-1)x[NUMBER OF VALID ORBITS].jld2")
else
    println(folderpath_o*"orbWeights_"*tokamak*"_"*TRANSP_id*"_at"*timepoint*"s_"*diagnostic_name*"_"*pretty2scpok(reaction_full; projVel = analyticalOWs)*"_1.jld2")
    println("... ")
    println(folderpath_o*"orbWeights_"*tokamak*"_"*TRANSP_id*"_at"*timepoint*"s_"*diagnostic_name*"_"*pretty2scpok(reaction_full; projVel = analyticalOWs)*"_$(iiimax).jld2")
end
if debug
    println("")
    println("!!!!!! DEBUGGING SPECIFIED. ALGORITHM WILL DEBUG. !!!!!!")
    println("")
end
println("Please remove previously saved files with the same file name (if any) prior to script completion. Quickly!")
println("")
println("If you would like to change any settings, please edit the start_calcOW_template.jl file or similar.")
println("")
println("Written by Henrik Järleblad. Last maintained 2022-10-05.")
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
# Helper function (to make the code more transparent)
verbose && println("Loading helper functions... ")
@everywhere begin
    include("helper/calcOrbSpec.jl") # Helper function for computing diagnostic spectra from orbits
end

## If no thermal distribution has been specified, we are going to need the default temp. and dens. profiles
if filepath_thermal_distr==""
    @everywhere begin
        include("misc/temp_n_dens.jl")
    end
end

## ---------------------------------------------------------------------------------------------
# Calculating orbit grid
verbose && debug && println("")
if !(og_filepath===nothing)
    verbose && println("Orbit grid and pertaining valid orbits were found in "*og_filepath*"... ")
else
    verbose && println("Calculating the orbit grid... ")
    og_orbs, og = OrbitTomography.orbit_grid(M, E_array, pm_array, Rm_array; amu=getSpeciesAmu(FI_species), q=getSpeciesEcu(FI_species), wall=wall, extra_kw_args...)
    og = nothing # Memory minimization
    # og_orbs is a vector of Orbit Julia-structs (found in L. Stagner's GuidingCenterOrbits.jl/orbit.jl package file)
# og is an OrbitGrid Julia-struct (found in L.Stagner's OrbitTomography.jl/orbits.jl package file)
end
norbs = length(og_orbs) # The number of valid orbits for the orbit grid og

# For de-bugging purposes
if debug
    println("Debugging specified. Saving norbs file... ")
    myfile = jldopen(folderpath_o*"norbs_"*diagnostic_name*"_$(nE)x$(npm)x$(nRm).jld2",true,true,false,IOStream)
    write(myfile,"norbs",norbs)
    close(myfile)
end
## ---------------------------------------------------------------------------------------------
F_os = 1.0 .*ones(size(og_orbs)) # Assume one ion for every orbit. This is the default and should not be changed lightly.

## ---------------------------------------------------------------------------------------------
# Setting Python variables and structures on all distributed workers/processes...
verbose && println("Setting all Python variables and structures on all distributed workers/processes... ")
@everywhere begin
    py"""
    # The '$' in front of many Python variables means that the variable is defined in Julia, not in Python.
    reaction = $reaction
    thermal_species = $thermal_species
    analyticalOWs = $analyticalOWs
    if $verbose:
        print("From Python: Loading forward model with diagnostic... ") 
    forwardmodel = forward.Forward($diagnostic_filepath) # Pre-initialize the forward model

    # Load TRANSP simulation data
    if (not $filepath_FI_cdf=="") and (not analyticalOWs): # If there is some TRANSP_id specified and we do not want to simply compute projected velocities...
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
    """
end
Ed_array = vec(py"Ed_vals")
@everywhere Ed_array = $Ed_array

## ---------------------------------------------------------------------------------------------
# If a .jld2 file has been specified for the thermal species distribution, we will need interpolation objects
if fileext_thermal=="jld2"
    verbose && println("Creating thermal temperature and density interpolations objects... ")
    thermal_temp_itp = Interpolations.interpolate((ρ_pol_array,), thermal_temp_array, Gridded(Linear()))
    thermal_dens_itp = Interpolations.interpolate((ρ_pol_array,), thermal_dens_array, Gridded(Linear()))
    thermal_temp_etp = Interpolations.extrapolate(thermal_temp_itp,0) # Add what happens in extrapolation scenario. 0 means we assume vacuum scrape-off layer (SOL)
    thermal_dens_etp = Interpolations.extrapolate(thermal_dens_itp,0) # Add what happens in extrapolation scenario. 0 means we assume vacuum scrape-off layer (SOL)
    thermal_temp = thermal_temp_etp
    thermal_dens = thermal_dens_etp
end

# If a thermal species distribution has not been specified, we simply need the thermal species temperature and density on axis
if filepath_thermal_distr==""
    thermal_temp = thermal_temp_axis
    thermal_dens = thermal_dens_axis
end

# If a TRANSP .cdf file and a TRANSP FI .cdf file have been specified for the thermal species distribution, we have everything we need in the Python thermal_dist variable
if lowercase(fileext_thermal)=="cdf" && isfile(filepath_FI_cdf)
    thermal_temp = nothing
    thermal_dens = nothing
end

# If there is a specified valid timepoint,
# and if a TRANSP .cdf file has been specified, but NOT a TRANSP FI .cdf file, 
# and we are NOT computing analytical orbit weight functions
if typeof(timepoint)==String && length(split(timepoint,","))==2 && lowercase(fileext_thermal)=="cdf" && !isfile(filepath_FI_cdf) && !analyticalOWs
    thermal_temp = getTempProfileFromTRANSP(timepoint, filepath_thermal_distr; verbose=verbose)
    thermal_dens = getDensProfileFromTRANSP(timepoint, filepath_thermal_distr, thermal_species; verbose=verbose)
end

# Transfer the thermal_temp and thermal_dens variables to external processes
@everywhere thermal_temp = $thermal_temp
@everywhere thermal_dens = $thermal_dens

## ---------------------------------------------------------------------------------------------
# Calculating the orbit weights
verbose && println("Starting the "*diagnostic_name*" weights calculations... ")
for iii=1:iiimax
    global instrumental_response
    global instrumental_response_matrix
    global Ed_array
    verbose && println("iii: $(iii)")
    if distributed && !debug # If parallel computating is desired (and !debug)...
        if visualizeProgress # if you want the progress to be visualized...
            prog = Progress(norbs) # Create a progress bar the is norbs long
            channel = RemoteChannel(()->Channel{Bool}(norbs), 1) # Utilize a channel
            Wtot = fetch(@sync begin
                @async while take!(channel)
                    ProgressMeter.next!(prog)
                end
                @async begin
                    W = @distributed (+) for i=1:norbs
                        spec = calcOrbSpec(M, og_orbs[i], F_os[i], py"forwardmodel", py"thermal_dist", py"Ed_bin_edges", reaction; thermal_temp=thermal_temp, thermal_dens=thermal_dens) # Calculate the diagnostic energy spectrum for the orbit
                        rows = append!(collect(1:length(spec)),length(spec)) # To be able to tell the sparse framework about the real size of the weight matrix
                        cols = append!(i .*ones(Int64, length(spec)), norbs) # To be able to tell the sparse framework about the real size of the weight matrix

                        # NOTE! Important to have this append!(spec,0.0) line AFTER the rows and cols assignment lines
                        append!(spec,0.0) # To be able to tell the sparse framework about the real size of the weight matrix
                        Wi = dropzeros(sparse(rows,cols,spec)) # Create a sparse weight matrix, with the current (i) column non-zero. Remove all redundant zeros.

                        put!(channel, true) # Update the progress bar
                        Wi # Declare this orbit weight to be added to the parallel computing reduction (@distributed (+))
                    end
                    put!(channel, false) # Update the progress bar
                    W # Declare the total weight function as ready for fetching
                end
            end)
        else
            Wtot = @distributed (+) for i=1:norbs
                spec = calcOrbSpec(M, og_orbs[i], F_os[i], py"forwardmodel", py"thermal_dist", py"Ed_bin_edges", reaction; thermal_temp=thermal_temp, thermal_dens=thermal_dens) # Calculate the diagnostic energy spectrum for the orbit
                rows = append!(collect(1:length(spec)),length(spec)) # Please see similar line earlier in the script
                cols = append!(i .*ones(Int64, length(spec)), norbs) # Please see similar line earlier in the script

                # NOTE! Important to have this append!(spec,0.0) line AFTER the rows and cols assignment lines
                append!(spec,0.0) # Please see similar line earlier in the script
                Wi = dropzeros(sparse(rows,cols,spec)) # Please see similar line earlier in the script
                Wi # Please see similar line earlier in the script
            end
        end
    else # ... if you do not use multiple cores, good luck!
        if debug
            verbose && println("Debugging specified. Only single-threaded debugging is possible. Will debug.")
            # WRITE CODE TO DEBUG QUANTITIES OF INTEREST

        else
            spec = calcOrbSpec(M, og_orbs[1], F_os[1], py"forwardmodel", py"thermal_dist", py"Ed_bin_edges", reaction; thermal_temp=thermal_temp, thermal_dens=thermal_dens) # Calculate the diagnostic energy spectrum for the orbit
            rows = append!(collect(1:length(spec)),length(spec)) # # Please see similar line earlier in the script
            cols = append!(1 .*ones(Int64, length(spec)), norbs) # # Please see similar line earlier in the script

            # NOTE! Important to have this append!(spec,0.0) line AFTER the rows and cols assignment lines
            append!(spec,0.0)
            Wtot = dropzeros(sparse(rows,cols,spec)) # Pre-allocate a sparse matrix
        end

        for i=2:norbs
            if debug
                verbose && println("Debugging orbit $(i) of $(norbs)... ")
                # WRITE CODE TO DEBUG QUANTITIES OF INTEREST

            else
                verbose && println("Calculating spectra for orbit $(i) of $(norbs)... ")
                local spec = calcOrbSpec(M, og_orbs[i], F_os[i], py"forwardmodel", py"thermal_dist", py"Ed_bin_edges", reaction; thermal_temp=thermal_temp, thermal_dens=thermal_dens) # Calculate the diagnostic energy spectrum for the orbit
                local rows = append!(collect(1:length(spec)),length(spec)) # Please see similar line earlier in the script
                local cols = append!(i .*ones(Int64, length(spec)), norbs) # Please see similar line earlier in the script

                # NOTE! Important to have this append!(spec,0.0) line AFTER the rows and cols assignment lines
                append!(spec,0.0) # Please see similar line earlier in the script
                Wtot += dropzeros(sparse(rows,cols,spec)) # Please see similar line earlier in the script
            end
        end
    end

    if instrumental_response
        verbose && println("Applying diagnostic response to weight functions... ")
        lo = findfirst(x-> x>minimum(Ed_array),instrumental_response_input)
        hi = findlast(x-> x<maximum(Ed_array),instrumental_response_input)
        if (typeof(lo)==Nothing) || (typeof(hi)==Nothing)
            @warn "Instrumental response matrix input completely outside weight function measurement bin range. No diagnostic response will be applied."
            instrumental_response = false
        else
            if lo==1
                @warn "Lower bound of instrumental response matrix input might not be low enough to cover weight function measurement bin range. Diagnostic response representation might be inaccurate."
            end
            if hi==length(instrumental_response_input)
                @warn "Upper bound of instrumental response matrix input might not be high enough to cover weight function measurement bin range. Diagnostic response representation might be inaccurate."
            end
            Wtot_withDiagResp = zeros(length(instrumental_response_output),size(Wtot,2))
            instrumental_response_matrix = (instrumental_response_matrix[lo:hi,:])'
            for io=1:size(Wtot,2)
                W = Wtot[:,io]
                itp = LinearInterpolation(Ed_array,W)
                W_itp = itp.(instrumental_response_input[lo:hi])
                W_out = instrumental_response_matrix * W_itp # The diagnostic response
                Wtot_withDiagResp[:,io] = W_out
            end
            Wtot = Wtot_withDiagResp # Update the outputs of calcOrbWeights.jl with the diagnostic response
            Ed_array = instrumental_response_output # Update the outputs of calcOrbWeights.jl with the diagnostic response
        end
    end

    if !debug
        verbose && println("Saving orbit weight function matrix in its 2D form... ")
        if iiimax==1 # If you intend to calculate only one weight function
            global filepath_output_orig = folderpath_o*"orbWeights_"*tokamak*"_"*TRANSP_id*"_at"*timepoint*"s_"*diagnostic_name*"_"*pretty2scpok(reaction_full; projVel = analyticalOWs)*"_$(length(vec(py"Ed_vals")))x$(norbs)"
        else # If you intend to calculate several (identical) weight functions
            global filepath_output_orig = folderpath_o*"orbWeights_"*tokamak*"_"*TRANSP_id*"_at"*timepoint*"s_"*diagnostic_name*"_"*pretty2scpok(reaction_full; projVel = analyticalOWs)*"_$(iii)"
        end
        global filepath_output = deepcopy(filepath_output_orig)
        global count = 1
        while isfile(filepath_output*".jld2") # To take care of not overwriting files. Add _(1), _(2) etc
            global filepath_output = filepath_output_orig*"_($(Int64(count)))"
            global count += 1 # global scope, to surpress warnings
        end
        global filepath_output = filepath_output*".jld2"
        myfile_s = jldopen(filepath_output,true,true,false,IOStream)
        write(myfile_s, "W", Wtot)
        write(myfile_s, "E_array", vec(E_array))
        write(myfile_s, "pm_array", vec(pm_array))
        write(myfile_s, "Rm_array", vec(Rm_array))
        write(myfile_s, "Ed_array", Ed_array)
        if instrumental_response
            write(myfile_s, "instrumental_response_input", instrumental_response_input)
            write(myfile_s, "instrumental_response_output", instrumental_response_output)
            write(myfile_s, "instrumental_response_matrix", instrumental_response_matrix)
        end
        write(myfile_s, "reaction_full", reaction_full)
        if analyticalOWs
            write(myfile_s, "analyticalOWs", analyticalOWs)
        end
        write(myfile_s, "filepath_thermal_distr", filepath_thermal_distr)
        write(myfile_s, "extra_kw_args", extra_kw_args)
        if !(og_filepath===nothing)
            write(myfile_s, "og_filepath", og_filepath)
        end
        close(myfile_s)
        if include2Dto4D
            (iiimax==1) && verbose && println("Enflating orbit weight functions to 4D form... ")
            (iiimax>1) && verbose && println("Enflating orbit weight functions $(iii) (of $(iiimax)) to 4D form... ")
            global filepath_W = filepath_output
            include("helper/orbWeights_2Dto4D.jl")
        end
    else
        verbose && println("Saving debugged quantities... ")
        # WRITE WHATEVER CODE TO SAVE THE DEBUGGED QUANTITIES
    end
end # The end of the iii=1:iiimax for-loop. (to enable computation of several identical orbit weight functions, to post-analyze MC noise influence for example)
println("------ Done! ------")
