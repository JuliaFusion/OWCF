#########################################  calcOrbWeights.jl #########################################

#### Description:
# This script computes orbit weight functions, it's as simple as that. It can be considered the flagship
# of the OWCF. As the rest of the OWCF, it utilizes Luke Stagner's orbit computation codes to compute
# the guiding-center (drift) orbits of fast ions in toroidally symmetric fusion devices. The orbit computation
# codes are (as of 2025-07-11) packaged into the Julia packages GuidingCenterOrbits.jl and OrbitTomography.jl.
# In the current version of the OWCF, the calcOrbWeights.jl script can compute orbit weight functions in 
# the following ways:
#
#   - By sending (E,p,R,z) points into the DRESS code (J. Eriksson et al, CPC, 199, 40-46, 2016)
#   - By using (E,p,R,z) points to compute projected fast-ion velocities
#   - By using (E,p,R,z) points together with the equations in A. Valentini et al, Nucl. Fusion (2025)
#
# The (E,p,R,z) points are computed for each orbit, using the GuidingCenterOrbits.jl and OrbitTomography.jl 
# packages. The orbit can be uniquely identified using an (E,pm,Rm) coordinate where E is the energy (in keV)
# pm is the pitch (v/v_||) at the maximum major radius point Rm of the orbit. Regardless of which approach is 
# chosen from the list above, an expected 1D signal for the (E,pm,Rm) coordinate is returned and placed into a matrix.
# Once signals for all (E,pm,Rm) coordinates of interest have been computed, the orbit weight function matrix 
# is complete and saved into an output file. In future versions of the OWCF, more codes/models can be added 
# to the list e.g. FIDASIM.
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
# specify a .cdf TRANSP shot file as filepath_thermal_distr, you have to (usually) specify a .cdf TRANSP fast-ion file
# for filepath_FI_cdf as well. This is to let the script know the correct time windows to extract from
# the TRANSP shot file. Otherwise, the script will try to load the TRANSP bulk (thermal) plasma distribution 
# data using either the .eqdsk-inferred timepoint, or (prioritized) user-specified timepoint (in seconds).
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
# orbWeights_[FLR]_[tokamak]_[TRANSP_id]_at[timepoint]s_[diagnostic_name]_[nE]x[npm]x[nRm].jld2 - If iiimax == 1
# orbWeights_[FLR]_[tokamak]_[TRANSP_id]_at[timepoint]s_[diagnostic_name]_$(i).jld2 - If iiimax != 1
# Regardless of saved file name, this saved file will have the fields:
#   W - The computed orbit weights. Dimensions are (channels, valid orbits) - Array{Float64,2}
#   E_array - The fast-ion energy grid array used for orbit space - Array{Float64,1}
#   pm_array - The fast-ion pm grid array used for orbit space - Array{Float64,1}
#   Rm_array - The fast-ion Rm grid array used for orbit space - Array{Float64,1}
#   Ed_array - The diagnostic energy bin centers - Array{Float64,1}
#   reaction - The nuclear fusion reaction for which the orbit weights are computed - String
#   filepath_thermal_distr - The filepath of the thermal species distribution. For reference - String
# If an orbit-space grid file was used to define the orbit grid for the orbit weight functions, the saved file will also have the key
#   og_filepath - The path to the .jld2-file (should be an output of calcOrbGrid.jl) used as orbit-space grid input - String

### Other
# Please note that the diagnostic energy grid will be created as bin centers.
# That is, the first diagnostic energy grid value will be (Ed_min+Ed_diff/2) and so on.

# Script written by Henrik Järleblad. Last maintained 2025-07-11.
################################################################################################

## ---------------------------------------------------------------------------------------------
verbose && println("Loading Julia packages... ")
@everywhere begin
    cd(folderpath_OWCF) # Necessary to move all the workers to the correct folder
    using EFIT # For calculating magn½etic equilibrium quantities
    using Equilibrium # For loading flux function data, tokamak geometry data etc.
    using FileIO # To write/open files in general
    using JLD2 # To write/open .jld2 files (Julia files, basically)
    using GuidingCenterOrbits # For calculating guiding-center orbits
    using Interpolations # To be able to interpolate, if no thermal distribution is specified
    using LinearAlgebra # To be able to easily take the norm of arrays
    using NetCDF # To enable write/open .cdf files
    using OrbitTomography # This is what all this is about!
    using ProgressMeter # To display computational progress during parallel computations
    using Printf # To be able to use specific String prints
    using PyCall # For using Python code in Julia
    using SparseArrays # To enable utilization of sparse matrices/vectors
    plot_results && (using LaTeXStrings) # If results are to be plotted, we need LaTeXStrings.jl for text formatting
    plot_results && (using Plots) # If results are to be plotted, we need the Plots.jl package
    include("misc/availReacts.jl") # To examine fusion reaction and extract thermal and fast-ion species
    include("misc/rewriteReacts.jl") # To rewrite a fusion reaction from the A(b,c)D format to the A-b=c-D format
    include("misc/species_func.jl") # To convert species labels to particle mass and charge
    include("misc/temp_n_dens.jl")
    (include2Dto4D || plot_results) && (include("extra/dependencies.jl")) # Load 2D-to-4D tool from the OWCF dependencies
    pushfirst!(PyVector(pyimport("sys")."path"), "") # To add the forward, transp_dists, transp_output and vcone modules (scripts in current path)
end

## ---------------------------------------------------------------------------------------------
## Fusion reaction, diagnostic and particle-related checks
# Reaction
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

# Charge of fusion product particle of interest
emittedParticleHasCharge = false
if !(getSpeciesCharge(getEmittedParticle(reaction))==0) && lowercase(getEmittedParticleEnergyLevel(reaction))=="gs" && !(getReactionForm(reaction)==3)
    verbose && println("")
    verbose && println("The fusion product particle of interest ($(getEmittedParticle(reaction))) of the $(reaction) reaction has non-zero charge!")
    verbose && println("- A 1-step or 2-step gamma-ray reaction is NOT assumed, since the energy level of the $(getEmittedParticle(reaction)) particle is specified to be in ground state (GS).")
    verbose && println("- Computation of orbit weight function from projected velocities is also NOT assumed, since the 'reaction' input variable is NOT on form (3)($(reaction))")
    verbose && println("---> For emitted particles with non-zero charge, the OWCF currently only supports computing the expected energy spectrum from the plasma as a whole (4*pi emission).")
    verbose && println("---> Therefore, the 'diagnostic_name' and 'diagnostic_filepath' input variables will be forcibly set to \"\".")
    verbose && println("")
    diagnostic_name = ""
    diagnostic_filepath = ""
    emittedParticleHasCharge = true
end

# Energy state of the nucleus of the fusion product particle of interest
if !(lowercase(getEmittedParticleEnergyLevel(reaction))=="gs") && (diagnostic_filepath=="")
    error("The fusion product particle of interest ($(getEmittedParticle(reaction))) of the $(reaction) reaction is on an excited state, but no diagnostic line-of-sight was specified (the diagnostic_filepath input variable was left unspecified). This is not allowed. Please correct and re-try.")
end
if !(lowercase(getEmittedParticleEnergyLevel(reaction))=="gs")
    verbose && println("")
    verbose && println("The fusion product particle of interest ($(getEmittedParticle(reaction))) of the $(reaction) reaction is on an excited state!")
    verbose && println("The OWCF will calculate the spectrum of gamma-rays emitted from the de-excitation of the $(getEmittedParticle(reaction)) particle, towards the specified detector.")
    verbose && println("")
end

# Currently, the analytic equations in OWCF/forward.jl do not support 4*pi spherical emission. I.e. a synthetic diagnostic viewing cone model file 
# must be specified. Otherwise, an error is thrown
if analytic && !isfile(diagnostic_filepath)
    error("Input variable 'analytic' was set to true, but the 'diagnostic_filepath' was not specified as a String with the file path to a valid diagnostic viewing cone model file (i.e. an output file of either the OWCF/extra/createCustomLOS.jl tool or the LINE21 code). 4*pi (spherical) emission is currently not supported for analytic==true. Please correct and re-try.")
end

## ---------------------------------------------------------------------------------------------
# Determine filepath_thermal_distr file extension
fileext_thermal = lowercase((split(filepath_thermal_distr,"."))[end]) # Assume last part after final '.' is the file extension
fileext_FI_cdf = lowercase((split(filepath_FI_cdf,"."))[end]) # Assume last part after final '.' is the file extension
@everywhere fileext_thermal = $fileext_thermal
@everywhere fileext_FI_cdf = $fileext_FI_cdf

## ---------------------------------------------------------------------------------------------
# Misc error checks
if (nprocs()>1) && !distributed
    error(ErrorException("Number of processes greater than 1, but single-threaded computation specified. Please set distributed to true."))
end

if !(fileext_thermal=="cdf" || fileext_thermal=="jld2" || fileext_thermal=="")
    println("Thermal distribution file format: ."*fileext_thermal)
    error("Unknown thermal distribution file format. Please re-specify file and re-try.")
end

## ---------------------------------------------------------------------------------------------
# Determine fast-ion and thermal species from inputs in start file
thermal_reactant, fast_reactant = getFusionReactants(reaction) # Check the specified fusion reaction, and extract thermal and fast-ion species
@everywhere thermal_reactant = $thermal_reactant # Transfer variable to all external processes
@everywhere fast_reactant = $fast_reactant # Transfer variable to all external processes

# projVel variable. To clarify when orbit weight functions are computed from projected velocities, in code below
projVel = getReactionForm(reaction)==3 ? true : false # If fusion reaction is specified as a single particle species..
@everywhere projVel = $projVel # Export 'projVel' variable to all external CPU processors

## ---------------------------------------------------------------------------------------------
# Loading thermal information and TRANSP RUN-ID, and perform checks
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
# Loading magnetic equilibrium and deduce timepoint, if possible
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
psi_axis, psi_bdry = psi_limits(M) # The limits of the flux function
@everywhere M = $M
@everywhere wall = $wall
@everywhere jdotb = $jdotb
@everywhere psi_axis = $psi_axis
@everywhere psi_bdry = $psi_bdry

if (fileext_thermal=="cdf") && (fileext_FI_cdf=="cdf")
    # If the user has specified a TRANSP .cdf file with pertaining NUBEAM fast-ion distribution data...
    # Load the time, and overwrite timepoint. TRANSP time data superseeds .eqdsk time data
    TIME = round((ncread(filepath_FI_cdf,"TIME"))[1],digits=4)
    TIME_array = split("$(TIME)",".") # Will be on format XX.YYYY
    XX = TIME_array[1]
    YYYY = TIME_array[2]
    timepoint = XX*","*YYYY # Format XX,YYYY to avoid "." when including in filename of saved output
end
@everywhere timepoint = $timepoint

## ---------------------------------------------------------------------------------------------
# Defining orbit grid vectors
if isfile(og_filepath)
    verbose && println("Filepath to .jld2 file containing orbit grid was specified. Loading orbit grid... ")
    myfile = jldopen(og_filepath,false,false,false,IOStream)
    og = myfile["og"]
    og_orbs = myfile["og_orbs"]
    fast_reactant_loaded = myfile["fast_reactant"]
    extra_kw_args = myfile["extra_kw_args"]
    close(myfile)
    E_array = og.energy
    pm_array = og.pitch
    Rm_array = og.r
    if !(lowercase(fast_reactant)==lowercase(fast_reactant_loaded))
        @warn "Fast-ion species ("*fast_reactant_loaded*") loaded from orbit-grid .jld2 file does not match fast-ion species in specified fusion reaction ("*fast_reactant*"). Did you confuse the thermal species ("*thermal_reactant*") for the fast-ion species ("*fast_reactant*")?"
    end
    if !(lowercase(fast_reactant)==lowercase(fast_reactant_loaded)) && !(lowercase(thermal_reactant)==lowercase(fast_reactant_loaded))
        error("Fast-ion species ("*fast_reactant_loaded*") loaded from orbit-grid .jld2 file does not match any reactants in specified fusion reaction ("*reaction*"). Please correct and re-try.")
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

# Convert grid vector elements to Float64, just in case of any user-defined Int64 elements
E_array = Float64.(E_array)
pm_array = Float64.(pm_array)
Rm_array = Float64.(Rm_array)

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
println("---------------------------------------------------- calcOrbWeights.jl ----------------------------------------------------")
if !(tokamak=="")
    print("Tokamak: "*tokamak*"      ")
else
    print("Tokamak: N/A"*"      ")
end
if !(TRANSP_id=="")
    print("TRANSP ID: "*TRANSP_id*"      ")
else
    print("TRANSP ID: N/A"*"      ")
end
if !(diagnostic_name=="")
    println("Diagnostic (name): "*diagnostic_name)
else
    println("Diagnostic (name): N/A")
end
println("")
if !(diagnostic_filepath=="")
    println("Diagnostic (file): "*diagnostic_filepath)
else
    println("Diagnostic (file): N/A")
end
if instrumental_response
    println("Detector instrumental response included? Yes!")
    println("---> Instrumental response: $(instrumental_response_filepath)")
else
    println("---> Detector instrumental response included? No.")
end
println("")
if !(filepath_thermal_distr=="")
    println("Bulk (thermal) plasma distribution: "*filepath_thermal_distr)
    println("---> Thermal species specified: "*thermal_reactant)
elseif (filepath_thermal_distr=="") && !projVel
    println("Thermal distribution file not specified.")
    println("---> Thermal ion ($(thermal_reactant)) temperature on-axis will be set to $(thermal_temp_axis) keV.")
    println("------> Thermal ion ($(thermal_reactant)) temperature profile: $(thermal_temp_profile_type)")
    println("---> Thermal ion ($(thermal_reactant)) density on-axis will be set to $(thermal_dens_axis) m^-3.")
    println("------> Thermal ion ($(thermal_reactant)) density profile: $(thermal_dens_profile_type)")
else
    println("Bulk (thermal) plasma distribution: N/A (projected velocities of $(getFastParticleSpecies(reaction)) ions will be used as weights)")
end
println("")
println("Magnetic (B-field) equilibrium: "*filepath_equil)
println("---> |B| on-axis: $(round(norm(Equilibrium.Bfield(M,magnetic_axis(M)...)),digits=2)) T")
println("")
if !projVel
    println("Fusion reaction specified: "*reaction)
else
    println("Fusion reaction specified: N/A (projected velocities computation)")
end
println("---> Fast-ion species specified: "*fast_reactant)
if emittedParticleHasCharge && !projVel
    println("PLEASE NOTE! The emitted "*getEmittedParticle(reaction)*" particle of the "*reaction*" reaction has non-zero charge!")
    println("PLEASE NOTE! The resulting energy distribution for "*getEmittedParticle(reaction)*" from the plasma as a whole will be computed.")
end
println("")
if !analytic
    println("To compute the up-/down-shift of the nominal birth energy of the $(getEmittedParticle(reaction)) particle, for every guiding-center point")
    println("along the (drift) orbit trajectory, the number of velocity vector samples will be: $(n_gyro)")
else
    println("Analytic equations in OWCF/forward.jl will be used to compute the synthetic diagnostic measurements.")
end
print("---> Finite Larmor radius (FLR) (spatial) effects included? ")
if include_FLR_effects
    println("Yes!")
else
    println("No.")
end
println("")
if distributed
    println("Parallel computing will be used with $(nprocs()) processes (1 main + $(nprocs()-1) workers).")
else
    println("Single-threaded computing the weights... Good luck!")
end
println("")
sWM = iiimax > 1 ? "matrices" : "matrix"
println("$(iiimax) $(sWM) will be computed.")
println("")
println("Orbit weight functions will be computed for a $(length(E_array))x$(length(pm_array))x$(length(Rm_array)) orbit-space grid with")
println("Energy: [$(minimum(E_array)),$(maximum(E_array))] keV")
println("Pitch maximum: [$(minimum(pm_array)),$(maximum(pm_array))]")
println("Radius maximum: [$(minimum(Rm_array)),$(maximum(Rm_array))] m")
if include2Dto4D
    println("---> Orbit weight matrix will be inflated and saved in its full 4D format upon completion of computations.")
end
println("")
if !projVel
    println("There will be $(length(range(Ed_min,stop=Ed_max,step=Ed_diff))-2) diagnostic energy bin(s) with")
    println("Lower diagnostic energy bound: $(Ed_min) keV")
    println("Upper diagnostic energy bound: $(Ed_max) keV")
else
    println("There will be $(length(range(Ed_min,stop=Ed_max,step=Ed_diff))-2) projected velocity (u) bin(s) with")
    println("Lower u bound: $(Ed_min) m/s")
    println("Upper u bound: $(Ed_max) m/s")
end
println("")
println("The equations-of-motion integration algorithm will use the extra keyword arguments: ")
println(extra_kw_args)
println("")
println("Results will be saved to: ")
sFLR = include_FLR_effects ? "FLR_" : ""
if iiimax == 1
    println(folderpath_o*"orbWeights_$(sFLR)"*tokamak*"_"*TRANSP_id*"_at"*timepoint*"s_"*diagnostic_name*"_"*pretty2scpok(reaction; projVel = projVel)*"_$(length(range(Ed_min,stop=Ed_max,step=Ed_diff))-2)x[NUMBER OF VALID ORBITS].jld2")
else
    println(folderpath_o*"orbWeights_$(sFLR)"*tokamak*"_"*TRANSP_id*"_at"*timepoint*"s_"*diagnostic_name*"_"*pretty2scpok(reaction; projVel = projVel)*"_1.jld2")
    println("... ")
    println(folderpath_o*"orbWeights_$(sFLR)"*tokamak*"_"*TRANSP_id*"_at"*timepoint*"s_"*diagnostic_name*"_"*pretty2scpok(reaction; projVel = projVel)*"_$(iiimax).jld2")
end
println("")
if plot_results
    println("--- Plots included ---> A figure with plots of (a subset of) the results will be saved in .png file format in the output folder.")
end
if debug
    println("")
    println("!!!!!! DEBUGGING SPECIFIED. ALGORITHM WILL DEBUG. !!!!!!")
    println("")
end
println("")
println("If you would like to change any settings, please edit the start_calcOW_template.jl file or similar.")
println("")
println("Written by Henrik Järleblad. Last maintained 2025-07-11.")
println("--------------------------------------------------------------------------------------------------------------------------")
println("")

## ---------------------------------------------------------------------------------------------
# Loading necessary Python and/or packages, depending on the forward model of choice.
# Currently, the three available forward models are
#   - Analytic equations (A. Valentini et al, Nucl. Fusion, 2025)
#   - Projected velocities
#   - The DRESS code (J. Eriksson et al, Comp. Phys. Comm., 2016)

verbose && println("Loading Python packages for forward model object... ")
@everywhere begin
    py"""
    import numpy as np
    import forward
    import spec
    import vcone
    """
end
if analytic
    verbose && println("Loading OWCF/forward.jl and OWCF/vcone.jl for analytic computations... ")
    @everywhere begin
        include("forward.jl") # Load analytic equations in the OWCF/forward.jl file
        include("vcone.jl") # Load the Julia version of vcone.py, to work with viewing cone data files purely in Julia
    end
end

## ---------------------------------------------------------------------------------------------
# Helper function (to make the code more transparent)
verbose && println("Loading helper functions... ")
@everywhere begin
    include("helper/calcOrbSpec.jl") # Helper function for computing diagnostic spectra from orbits
end

## If no thermal distribution has been specified, we are going to need the default temp. and dens. profiles
if !isfile(filepath_thermal_distr)
    @everywhere begin
        include("misc/temp_n_dens.jl")
    end
end

## ---------------------------------------------------------------------------------------------
# Calculating orbit grid
verbose && debug && println("")
if isfile(og_filepath)
    verbose && println("Orbit grid and pertaining valid orbits were found in "*og_filepath*"... ")
else
    verbose && println("Calculating the orbit grid... ")
    og_orbs, og = OrbitTomography.orbit_grid(M, E_array, pm_array, Rm_array; amu=getSpeciesAmu(fast_reactant), q=getSpeciesEcu(fast_reactant), wall=wall, extra_kw_args...)
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
        function getThermalSpeciesTemperature(R,z) # This function is not actually used (an equivalent function is used automatically in the forward.Forward object)
            temp = py"thermal_dist.get_temperature"(R,z)
            if length(R)==1
                return temp[1]
            end
            return temp
        end
        function getThermalSpeciesDensity(R,z) # This function is not actually used (an equivalent function is used automatically in the forward.Forward object)
            dens = py"thermal_dist.get_density"(R,z)
            if length(R)==1
                return dens[1]
            end
            return dens
        end
    end
end

## ---------------------------------------------------------------------------------------------
# Prepare functions for computing synthetic measurements, based on inputs in start file
verbose && println("Creating synthetic measurements computation function... ")
@everywhere begin
    if analytic # Analytic equations in OWCF/forward.jl as forward model
        function analytic_comp_helper(E, p, R, z, w)
            N = length(R) # Assume length(E)==length(p)==length(R)==length(z)==length(w)
            E = vcat(E); p = vcat(p); R = vcat(R); z = vcat(z); w = vcat(w) # If not Vectors, make into 1-element Vectors. If already Vectors, keep as is
            B_gc = zeros(3,N)
            for i in eachindex(R)
                B_gc[:,i] = collect(reshape(Equilibrium.Bfield(M,R[i],z[i]),3,1))
            end
            return E, p, R, z, w, B_gc
        end
        if include_FLR_effects
            function compute_measurements(E, p, R, z, w)
                E, p, R, z, w_gc, B_gc = analytic_comp_helper(E, p, R, z, w)
                py"""
                v, x = forwardmodel.add_gyration($E, $p, $R, $z, $B_gc, n_gyro=$n_gyro) # No need to specify particle species, already specified via Forward.__init__()
                """
                E_gyro = repeat(E, inner=n_gyro) # Vector (M,) where M = length(E)*n_gyro
                p_gyro = repeat(p, inner=n_gyro) # Vector (M,)
                R_gyro = py"x[0,:]" # Vector (M,)
                z_gyro = py"x[2,:]" # Vector (M,)
                w_gyro = inv(n_gyro) .*repeat(w_gc, inner=n_gyro) # Vector (N*n_gyro,), weights rescaled by number of gyro-orbit samples
                B_gyro = reduce(hcat,map(i-> Equilibrium.Bfield(M,R_gyro[i],z_gyro[i]), eachindex(R_gyro))) # Create a (3,M) array with the magnetic field vector at each (R,z) point
                bulk_dens = getThermalSpeciesDensity(R_gyro, z_gyro) # Bulk plasma densities in m^-3
                # Please note! Ed_array are the diagnostic measurement bin centers, NOT the edges as in py"Ed_bin_edges" for the DRESS code
                return analytic_calc(viewing_cone_model, E_gyro, p_gyro, R_gyro, z_gyro, w_gyro, Ed_array, B_gyro; reaction=reaction, bulk_dens=bulk_dens)
            end
        else
            function compute_measurements(E, p, R, z, w)
                E, p, R_gc, z_gc, w_gc, B_gc = analytic_comp_helper(E, p, R, z, w)
                bulk_dens = getThermalSpeciesDensity(R_gc, z_gc) # Bulk plasma density in m^-3
                # Please note! Ed_array are the diagnostic measurement bin centers, NOT the edges as in py"Ed_bin_edges" for the DRESS code
                return analytic_calc(viewing_cone_model, E, p, R_gc, z_gc, w_gc, Ed_array, B_gc; reaction=reaction, bulk_dens=bulk_dens)
            end
        end
    else # DRESS code as forward model
        function DRESS_comp_helper(E, p, R, z, w)
            N = length(R) # Assume length(E)==length(p)==length(R)==length(z)==length(w)
            E = vcat(E); p = vcat(p); R = vcat(R); z = vcat(z); w = vcat(w) # If not Vectors, make into 1-element Vectors. If already Vectors, keep as is
            B_vecs = zeros(3,N)
            for i in eachindex(R)
                B_vecs[:,i] = collect(reshape(Equilibrium.Bfield(M,R[i],z[i]),3,1))
            end
            py"""
            v, x = forwardmodel.add_gyration($E, $p, $R, $z, $B_vecs, n_gyro=$n_gyro) # No need to specify particle species, already specified via Forward.__init__()
            """
            weights = inv(n_gyro) .*repeat(w, inner=n_gyro) # Vector (N*n_gyro,), weights rescaled by number of gyro-orbit samples
            return R, z, weights, B_vecs, py"v", py"x" 
        end
        if include_FLR_effects
            if py"thermal_dist"=="" # Custom bulk (thermal) plasma distribution data
                function compute_measurements(E, p, R, z, w)
                    R_gc, z_gc, w_gyro, B_gc, v, x_gyro = DRESS_comp_helper(E, p, R, z, w)
                    # PLEASE NOTE! comp_helper() is only used to avoid repeating code. Note: length(E)=N but length(weights)=N*n_gyro
                    R_gyro = x_gyro[1,:] # Vector (M,) where M = length(E)*n_gyro
                    z_gyro = x_gyro[3,:] # Vector (M,)
                    B_gyro = reduce(hcat,map(i-> Equilibrium.Bfield(M,R_gyro[i],z_gyro[i]), eachindex(R_gyro))) # Create a (3,M) array with the magnetic field vector at each (R,z) point
                    bulk_temp = getThermalSpeciesTemperature(R_gyro, z_gyro) # Bulk plasma temperature in keV
                    bulk_dens = getThermalSpeciesDensity(R_gyro, z_gyro) # Bulk plasma density in m^-3
                    py"""
                    spec = forwardmodel.compute_spectrum(Ed_bin_edges, $R_gyro, $z_gyro, $v, $w_gyro, $B_gyro, bulk_temp=$bulk_temp, bulk_dens=$bulk_dens)
                    """
                    return vec(vcat(py"spec"))
                end
            else # TRANSP bulk (thermal) plasma distribution data
                function compute_measurements(E, p, R, z, w)
                    R_gc, z_gc, w_gyro, B_gc, v, x_gyro = DRESS_comp_helper(E, p, R, z, w)
                    # PLEASE NOTE! comp_helper() is only used to avoid repeating code. Note: length(E)=N but length(weights)=N*n_gyro
                    R_gyro = x_gyro[1,:] # Vector (M,) where M = length(E)*n_gyro
                    z_gyro = x_gyro[3,:] # Vector (M,)
                    B_gyro = reduce(hcat,map(i-> Equilibrium.Bfield(M,R_gyro[i],z_gyro[i]), eachindex(R_gyro))) # Create a (3,M) array with the magnetic field vector at each (R,z) point
                    # bulk_temp <--- Already in forwardmodel via instantiation
                    # bulk_dens <--- -||-
                    py"""
                    spec = forwardmodel.compute_spectrum(Ed_bin_edges, $R_gyro, $z_gyro, $v, $w_gyro, $B_gyro)
                    """
                    return vec(vcat(py"spec"))
                end
            end
        else
            if py"thermal_dist"=="" # Custom bulk (thermal) plasma distribution data
                function compute_measurements(E, p, R, z, w)
                    R_gc, z_gc, w_gyro, B_gc, v, x_gyro = DRESS_comp_helper(E, p, R, z, w)
                    # PLEASE NOTE! comp_helper() is only used to avoid repeating code. Note: length(E)=N but length(weights)=N*n_gyro
                    R_gc = repeat(R_gc, inner=n_gyro) # Ignore FLR effects. 'inner' keyword argument included to match with Python's default repeat() functionality
                    z_gc = repeat(z_gc, inner=n_gyro) # Ignore FLR effects
                    B_gc = repeat(B_gc, inner=(1,n_gyro)) # Array (3,M)
                    bulk_temp = getThermalSpeciesTemperature(R_gc, z_gc) # Bulk plasma temperature in keV
                    bulk_dens = getThermalSpeciesDensity(R_gc, z_gc) # Bulk plasma density in m^-3
                    py"""
                    spec = forwardmodel.compute_spectrum(Ed_bin_edges, $R_gc, $z_gc, $v, $w_gyro, $B_gc, bulk_temp=$bulk_temp, bulk_dens=$bulk_dens)
                    """
                    return vec(vcat(py"spec"))
                end
            else # TRANSP bulk (thermal) plasma distribution data
                function compute_measurements(E, p, R, z, w)
                    R_gc, z_gc, w_gyro, B_gc, v, x_gyro = DRESS_comp_helper(E, p, R, z, w)
                    # PLEASE NOTE! comp_helper() is only used to avoid repeating code. Note: length(E)=N but length(weights)=N*n_gyro
                    R_gc = repeat(R_gc, inner=n_gyro) # Ignore FLR effects. 'inner' keyword argument included to match with Python's default repeat() functionality
                    z_gc = repeat(z_gc, inner=n_gyro) # Ignore FLR effects
                    B_gc = repeat(B_gc, inner=(1,n_gyro)) # Array (3,M)
                    # bulk_temp <--- Already in forwardmodel via instantiation
                    # bulk_dens <--- -||-
                    py"""
                    spec = forwardmodel.compute_spectrum(Ed_bin_edges, $R_gc, $z_gc, $v, $w_gyro, $B_gc)
                    """
                    return vec(vcat(py"spec"))
                end
            end
        end
    end
end

## ---------------------------------------------------------------------------------------------
# Take care of local/global scope issues, if plotting and/or 2D->4D
if include2Dto4D && plot_results
    W4D = nothing # Pre-define this variable identifier, so that it has global scope
end

## ---------------------------------------------------------------------------------------------
# Calculating the orbit weights
verbose && println("Starting the "*diagnostic_name*" weights calculations... ")
for iii=1:iiimax
    global instrumental_response # Use the instrumental_response variable from the global scope
    global instrumental_response_matrix # -||- instrumental_response_matrix -||-
    global Ed_array # -|| Ed_array -||-
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
                        spec = calcOrbSpec(og_orbs[i], F_os[i], compute_measurements) # Calculate the expected diagnostic energy spectrum for the orbit, using the function compute_measurements()
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
                spec = calcOrbSpec(og_orbs[i], F_os[i], compute_measurements) # Calculate the expected diagnostic energy spectrum for the orbit, using the function compute_measurements()
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
            verbose && println("Calculating spectra for orbit 1 of $(norbs)... ")
            spec = calcOrbSpec(og_orbs[1], F_os[1], compute_measurements) # Calculate the expected diagnostic energy spectrum for the orbit, using the function compute_measurements()
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
                local spec = calcOrbSpec(og_orbs[i], F_os[i], compute_measurements) # Calculate the expected diagnostic energy spectrum for the orbit, using the function compute_measurements()
                local rows = append!(collect(1:length(spec)),length(spec)) # Please see similar line earlier in the script
                local cols = append!(i .*ones(Int64, length(spec)), norbs) # Please see similar line earlier in the script

                # NOTE! Important to have this append!(spec,0.0) line AFTER the rows and cols assignment lines
                append!(spec,0.0) # Please see similar line earlier in the script
                Wtot += dropzeros(sparse(rows,cols,spec)) # Please see similar line earlier in the script
            end
        end
    end

    Wtot[findall(x-> isnan(x), Wtot)] .= 0.0 # If NaNs, set them to 0.0

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
            Wtot_withInstrResp = zeros(length(instrumental_response_output),size(Wtot,2))
            instrumental_response_matrix = (instrumental_response_matrix[lo:hi,:])'
            for io=1:size(Wtot,2)
                W = Wtot[:,io]
                itp = LinearInterpolation(Ed_array,W)
                W_itp = itp.(instrumental_response_input[lo:hi])
                W_out = instrumental_response_matrix * W_itp # The diagnostic response
                Wtot_withInstrResp[:,io] = W_out
            end
            Wtot = Wtot_withInstrResp # Update the outputs of calcOrbWeights.jl with the diagnostic response
            Ed_array = instrumental_response_output # Update the outputs of calcOrbWeights.jl with the diagnostic response
        end
    end

    if !debug
        verbose && println("Saving orbit weight function matrix in its 2D form... ")
        if iiimax==1 # If you intend to calculate only one weight function
            global filepath_output_orig = folderpath_o*"orbWeights_$(sFLR)"*tokamak*"_"*TRANSP_id*"_at"*timepoint*"s_"*diagnostic_name*"_"*pretty2scpok(reaction; projVel = projVel)*"_$(length(vec(py"Ed_vals")))x$(norbs)"
        else # If you intend to calculate several (identical) weight functions
            global filepath_output_orig = folderpath_o*"orbWeights_$(sFLR)"*tokamak*"_"*TRANSP_id*"_at"*timepoint*"s_"*diagnostic_name*"_"*pretty2scpok(reaction; projVel = projVel)*"_$(iii)"
        end
        global filepath_output = deepcopy(filepath_output_orig)
        global count = 1
        while isfile(filepath_output*".jld2") # To take care of not overwriting files. Add _(1), _(2) etc
            global filepath_output = filepath_output_orig*"_($(Int64(count)))"
            global count += 1 # global scope, to surpress warnings
        end

        myfile_s = jldopen(filepath_output*".jld2",true,true,false,IOStream)
        write(myfile_s,"W2D",Wtot)
        write(myfile_s,"E_array",vec(E_array))
        write(myfile_s,"pm_array",vec(pm_array))
        write(myfile_s,"Rm_array",vec(Rm_array))
        write(myfile_s,"Ed_array",Ed_array)
        if instrumental_response
            write(myfile_s,"Ed_array_units",instrumental_response_output_units)
            write(myfile_s,"instrumental_response_input",instrumental_response_input)
            write(myfile_s,"instrumental_response_output",instrumental_response_output)
            write(myfile_s,"instrumental_response_matrix",instrumental_response_matrix)
        else
            write(myfile_s,"Ed_array_units",projVel ? "m_s^-1" : "keV") # Otherwise, the output abscissa of calcSpec.jl is always in m/s or keV
        end
        write(myfile_s, "reaction", reaction)
        if projVel
            write(myfile_s, "projVel", projVel)
        end
        write(myfile_s,"filepath_thermal_distr",filepath_thermal_distr)
        write(myfile_s,"extra_kw_args",extra_kw_args)
        if !(og_filepath===nothing)
            write(myfile_s,"og_filepath",og_filepath)
        end
        close(myfile_s)
        if include2Dto4D || plot_results
            splt = plot_results ? " (for plotting)" : ""
            verbose && println("Inflating orbit weight function matrix ($(iii) of $(iiimax)) to 4D form$(splt)... ")
            global W4D # Get the global scope variable identifier
            W4D = zeros(length(Ed_array),length(E_array),length(pm_array),length(Rm_array))
            for iEd in eachindex(Ed_array)
                W3D = OWCF_map_orbits(og, Vector(Wtot[iEd,:]), true; weights=true) # Assume equidistant orbit-space grid. For now.
                W4D[iEd,:,:,:] .= W3D
            end
        end
        if include2Dto4D
            verbose && println("Saving orbit weight function matrix in its inflated 4D form... ")
            fpath, fname = split(filepath_output, "/")[1:end-1], split(filepath_output, "/")[end] # Split filepath_output into path and filename
            fname_calcOW, fname_rest = split(fname, "_")[1], split(fname, "_")[2:end] # Split filename into 'orbWeights' and the rest
            fname_new = String(reduce(*, map(x-> x*"_", vcat(fname_calcOW*"4D", fname_rest)))[1:end-1]) # Change 'orbWeights' to 'orbWeights4D', then re-create the filename
            fpath = String(reduce(*,map(x-> x*"/", fpath))) # Re-create the path
            filepath_4D_output = fpath*fname_new # Create the new 4D full file path
            myfile = jldopen(filepath_4D_output*".jld2",true,true,false,IOStream)
            write(myfile,"W", W4D)
            write(myfile,"E_array", E_array)
            write(myfile,"pm_array", pm_array)
            write(myfile,"Rm_array", Rm_array)
            write(myfile,"Ed_array", Ed_array) # Save only diagnostic energy bins that have been used
            if (@isdefined Ed_array_units)
                write(myfile,"Ed_array_units", Ed_array_units)
            end
            if (@isdefined reaction)
                write(myfile,"reaction", reaction)
            end
            if (@isdefined reaction_full)
                write(myfile,"reaction_full",reaction_full)
            end
            if projVel
                write(myfile, "projVel", projVel)
            end
            if (@isdefined filepath_thermal_distr)
                write(myfile,"filepath_thermal_distr",filepath_thermal_distr)
            end
            close(myfile)
        end
        if plot_results

            plot_font = "Computer Modern"
            Plots.default(fontfamily=plot_font)
            verbose && println("Plotting (a subset of) the results of the orbit weight functions computations... ")

            # The indices of the measurement bin centers
            N_bins = length(Ed_array)
            if N_bins==1
                plt_Ed_inds = [1,1,1]
            elseif N_bins==2
                plt_Ed_inds = [1,2,2]
            elseif N_bins==3
                plt_Ed_inds = [1,2,3]
            else
                # Algorithm to find suitable indices of Ed to visualize
                W_gross = dropdims(sum(W4D, dims=(2,3,4)), dims=(2,3,4))
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
                plt_Ed_inds = [iEd_low, iEd_mid, iEd_hig]
            end
            iEd_low, iEd_mid, iEd_hig = plt_Ed_inds

            # The indices of the energy slices to be visualized
            global nE
            nE = length(E_array)
            if nE>=5
                plt_E_inds = Int64.(round.(collect(range(1, stop=nE, length=5))[2:4]))
            elseif nE==4
                plt_E_inds = [2,3,4]
            elseif nE==3
                plt_E_inds = [1,2,3]
            elseif nE==2
                plt_E_inds = [1,1,2]
            else # nE==1
                plt_E_inds = [1,1,1]
            end
            iE_low, iE_mid, iE_hig = plt_E_inds

            # The values of the measurement bin centers to be visualized, as Strings
            Ed_low = @sprintf "%.3E" Ed_array[iEd_low]
            Ed_mid = @sprintf "%.3E" Ed_array[iEd_mid]
            Ed_hig = @sprintf "%.3E" Ed_array[iEd_hig]
            units_plt = projVel ? "m/s" : "keV"
        
            # The values of the energy slices to be visualized
            E_low = @sprintf "%.2E" E_array[iE_low]
            E_mid = @sprintf "%.2E" E_array[iE_mid]
            E_hig = @sprintf "%.2E" E_array[iE_hig]

            my_color_map = cgrad([:white, :darkblue, :green, :yellow, :orange, :red])
            plt_lowlow = Plots.heatmap(Rm_array, pm_array, W4D[iEd_low,iE_low,:,:], xlabel = L"R_{m} \, [m]"*"\n"*L"\textbf{Ed:} \, \,"*"$(Ed_low) [$(units_plt)]", ylabel=L"p_{m} \, [-]", fillcolor=my_color_map, colorbar=false, title="E: $(E_low) keV")
            plt_lowmid = Plots.heatmap(Rm_array, pm_array, W4D[iEd_low,iE_mid,:,:], xlabel = L"R_{m} \, [m]",   ylabel=L"p_{m} \, [-]", fillcolor=my_color_map, colorbar=false, title="E: $(E_mid) keV")
            plt_lowhig = Plots.heatmap(Rm_array, pm_array, W4D[iEd_low,iE_hig,:,:], xlabel = L"R_{m} \, [m]",   ylabel=L"p_{m} \, [-]", fillcolor=my_color_map, colorbar=false, title="E: $(E_hig) keV")

            plt_midlow = Plots.heatmap(Rm_array, pm_array, W4D[iEd_mid,iE_low,:,:], xlabel = L"R_{m} \, [m]"*"\n"*L"\textbf{Ed:} \, \,"*"$(Ed_mid) [$(units_plt)]", ylabel=L"p_{m} \, [-]", fillcolor=my_color_map, colorbar=false, title="E: $(E_low) keV")
            plt_midmid = Plots.heatmap(Rm_array, pm_array, W4D[iEd_mid,iE_mid,:,:], xlabel = L"R_{m} \, [m]",   ylabel=L"p_{m} \, [-]"  , fillcolor=my_color_map, colorbar=false, title="E: $(E_mid) keV")
            plt_midhig = Plots.heatmap(Rm_array, pm_array, W4D[iEd_mid,iE_hig,:,:], xlabel = L"R_{m} \, [m]",   ylabel=L"p_{m} \, [-]"  , fillcolor=my_color_map, colorbar=false, title="E: $(E_hig) keV")

            plt_higlow = Plots.heatmap(Rm_array, pm_array, W4D[iEd_hig,iE_low,:,:], xlabel = L"R_{m} \, [m]"*"\n"*L"\textbf{Ed:} \, \,"*"$(Ed_hig) [$(units_plt)]", ylabel=L"p_{m} \, [-]", fillcolor=my_color_map, colorbar=false, title="E: $(E_low) keV")
            plt_higmid = Plots.heatmap(Rm_array, pm_array, W4D[iEd_hig,iE_mid,:,:], xlabel = L"R_{m} \, [m]",   ylabel=L"p_{m} \, [-]"  , fillcolor=my_color_map, colorbar=false, title="E: $(E_mid) keV")
            plt_highig = Plots.heatmap(Rm_array, pm_array, W4D[iEd_hig,iE_hig,:,:], xlabel = L"R_{m} \, [m]",   ylabel=L"p_{m} \, [-]"  , fillcolor=my_color_map, colorbar=false, title="E: $(E_hig) keV")

            plt_tot = Plots.plot(plt_lowhig, plt_midhig, plt_highig,
                                    plt_lowmid, plt_midmid, plt_higmid,
                                    plt_lowlow, plt_midlow, plt_higlow,
                                    layout=(3,3), size=(1200,1200), 
                                    top_margin=4Plots.mm, dpi=200)
            display(plt_tot)

            verbose && println("---> Saving figure at $(filepath_output).png... ")
            png(plt_tot, filepath_output)
        end
    else
        verbose && println("Saving debugged quantities... ")
        # WRITE WHATEVER CODE TO SAVE THE DEBUGGED QUANTITIES
    end
end # The end of the iii=1:iiimax for-loop. (to enable computation of several identical orbit weight functions, to post-analyze MC noise influence for example)
println("------ calcOrbWeights.jl completed successfully! ------")
