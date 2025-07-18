######################################## ps2WF.jl ##############################################

### Description:
# The need for this script arose after realizing that there is a need for a really fine energy resolution
# in orbit space when working with ICRF shots. The tail of the fast-ion distribution extends all
# the way up into the several MeVs, while still having the thermal of the population down below
# a hundred keVs. By using clever partitioning of both particle space and orbit space, it is possible
# to calculate an immensely fine-resolved computation of the WF signal. Basically,
# the calculation is divided into many smaller parts. This can be written as
#
# S = WF = (W_1)*(F_1) + (W_2)*(F_2) + ... + (W_n)*(F_n)
#
# Where '1' means the distribution and weight function spans E=1-10 keV (for example) and '2'
# means the distribution and weight function spans E=11-20 keV (for example, and so on).
# This partitioning can be made arbitrarily fine-resolved.
#
# If calcWFOs is set to true, the ps2WF.jl will also compute so-called FOs, WOs and WFOs. What are those?
# They are the fast-ion distribution f, the orbit weight functions w and signal densities wf split into their
# orbit-space constituents. These splits are computed for each dimensional direction in orbit space. That is:
#
# f(E) = f_{co-passing}(E) + f_{trapped}(E) + ...
# f(pm) = f_{co-passing}(pm) + f_{trapped}(pm) + ...
# f(Rm) = f_{co-passing}(Rm) + f_{trapped}(Rm) + ...
# w(E_{d},E) = w_{co-passing}(E_{d},E) + w_{trapped}(E_{d},E) + ...
# w(E_{d},pm) = w_{co-passing}(E_{d},pm) + w_{trapped}(E_{d},pm) + ...
# w(E_{d},Rm) = w_{co-passing}(E_{d},Rm) + w_{trapped}(E_{d},Rm) + ...
# wf(E_{d},E) = wf_{co-passing}(E_{d},E) + wf_{trapped}(E_{d},E) + ...
# wf(E_{d},pm) = wf_{co-passing}(E_{d},pm) + wf_{trapped}(E_{d},pm) + ...
# wf(E_{d},Rm) = wf_{co-passing}(E_{d},Rm) + wf_{trapped}(E_{d},Rm) + ...
#
# The signal densities are point-wise multiplied together for every point in 3D orbit space before
# being categorized into a co-passing/trapped/... contribution. The integration is performed via binning. 
#
# The orbit weight functions are normalized by the number of grid points of a specific orbit type, for each orbit type quantity.
# For example, for the co-passing energy dependence:
#
# w_{co-passing}(E_{d}, E) = w_{co-passing}(E_{d},E) / N_{co-passing}(E)
#
# where N_{co-passing}(E) is the number of co-passing grid points with a fast-ion energy E. This is to ensure consistency. For example,
# the metric is different in (E,pm,Rm) and (E,mu,Pphi;sigma) where E, mu, Pphi (and sigma) are the constants-of-motion
# coordinates. That is, the different topological regions might appear smaller of larger in the different coordinate systems.
# To ensure e.g. w_{co-passing}(E_{d},E) is the same regardlss of coordinate system, we normalize by the number of grid points for 
# each orbit type.

### The fast-ion distribution can be provided via a .h5 file or via a .jld2 file.
## The .h5 file must have the keys:
# 'f' - The 4D matrix containing the fast-ion distribution.
# 'energy' - The 1D array containing the energy grid points
# 'pitch ' - The 1D array containing the pitch grid points
# 'R' - The 1D array containing the R grid points
# 'Z' or 'z' - The 1D array containing the z grid points
# It can be that the order of the fast-ion data 'f' is reversed. That is, it has the dimensional order
# [z, R, pitch, energy]. If so, the OWCF will detect it and permute it accordingly.
#
## The .jld2 file must have the keys:
# 'F_ps' or 'F_EpRz' - The 4D matrix containing the fast-ion distribution. Order [energy,pitch,R,z] (forwards).
# 'energy' - The 1D array containing the energy grid points
# 'pitch ' - The 1D array containing the pitch grid points
# 'R' - The 1D array containing the R grid points
# 'Z' or 'z' - The 1D array containing the z grid points

#### Inputs (Units given when defined in script)
# Given easily via input file start_ps2WF_template.jl. 'template' should be renamed with whatever.

#### Outputs
# -

#### Saved files
# ps2WF_results_[FLR]_[tokamak]_[TRANSP_id]_at[timepoint]s_[diagnostic_name]_[nE]x[npm]x[nRm].jld2
# Regardless of saved file name, this saved file will have the fields:
#   S_WF - The computed WF signal - Array{Float64,2}
#   E_array - The total fast-ion energy grid array used for orbit space - Array{Float64,1}
#   pm_array - The fast-ion pm grid array used for orbit space - Array{Float64,1}
#   Rm_array - The fast-ion Rm grid array used for orbit space - Array{Float64,1}
#   Ed_array - The diagnostic energy bin centers - Array{Float64,1}
#   filepath_FI_distr - The path to the fast-ion distribution provided as input by the user - String
#   filepath_thermal_distr - The path to the thermal distribution/profiles provided as input by the user - String
#   nfast_orig - The total number of fast ions in the fast-ion distribution provided as input by the user - Float64 (or Int64)
#   nfast_chunks - The total number of fast ions split into chunks for every WF_n chunk. sum(nfast_chunks)==nfast_orig at least approximately. Otherwise, something went wrong. - Array{Float64,1}
#   extra_kw_args - The extra keyword arguments provided as input by the user. Used as keyword arguments in the orbit integration - Dictionary
# If calcWFOs is set to true, it will also have the fields
#   WO_E - The orbit weights as functions of diagnostic measurement bin, fast-ion energy and orbit type. w(E_d, E, o) where 'o' is the orbit type - Array{Float64,3}
#   WO_pm - The orbit weights as functions of diagnostic measurement bin, pitch maximum and orbit type. w(E_d, pm, o) where 'o' is the orbit type - Array{Float64,3}
#   WO_Rm - The orbit weights as functions of diagnostic measurement bin, radius maximum and orbit type. w(E_d, Rm, o) where 'o' is the orbit type - Array{Float64,3}
#   FO_E - The fast-ion distribution function as function of fast-ion energy and orbit type. f(E,o) where 'o' is the orbit type - Array{Float64,2}
#   FO_pm - The fast-ion distribution function as function of pitch maximum and orbit type. f(pm,o) where 'o' is the orbit type - Array{Float64,2}
#   FO_Rm - The fast-ion distribution function as function of radius maximum and orbit type. f(Rm,o) where 'o' is the orbit type - Array{Float64,2}
#   WFO_E - The synthetic diagnostic signal densities as functions of diagnostic measurement bin, fast-ion energy and orbit type. WFO_E(E_d, E, o) where 'o' is the orbit type. That is, sum(WFO_E, dims=(2,3))=WF - Array{Float64,3}
#   WFO_pm - The synthetic diagnostic signal densities as functions of diagnostic measurement bin, pitch maximum and orbit type. WFO_pm(E_d, pm, o) where 'o' is the orbit type. That is, sum(WFO_pm, dims=(2,3))=WF - Array{Float64,3}
#   WFO_E - The synthetic diagnostic signal densities as functions of diagnostic measurement bin, radius maximum and orbit type. WFO_Rm(E_d, Rm, o) where 'o' is the orbit type. That is, sum(WFO_Rm, dims=(2,3))=WF - Array{Float64,3}

### Other
# Please note that the diagnostic energy grid will be created up to and NOT including Ed_max.
# That is, the resulting diagnostic energy grid array will have the boundaries [Ed_min,Ed_max).
#
# Also, this script assumes that the orbit grid is equidistant.
#
# Lastly, if the script is terminated mid-execution without having finished, a progress file will
# have been created (in this/the execution folder). Just restart the script via
# 'julia start_ps2WF_template.jl' (for example) and the program will take care of the rest. The
# script will then resume at the last checkpoint, and continue until the script is finished.
# The progress file is then deleted once ps2WF.jl completes successfully.

# Script written by Henrik Järleblad. Last maintained 2025-07-11.
################################################################################################


## ---------------------------------------------------------------------------------------------
# Determine the file extensions for the fast-ion and thermal species distribution files
# This is the very first thing that needs to happen, because it determines which packages and dependencies
# files will be loaded.
@everywhere begin
    fileext_FI = (split(filepath_FI_distr,"."))[end] # Assume last part after final '.' is the file extension
    fileext_FI = lowercase(fileext_FI)
    fileext_thermal = (split(filepath_thermal_distr,"."))[end] # Assume last part after final '.' is the file extension
    fileext_thermal = lowercase(fileext_thermal)
    fileext_FI_TRANSP_shot = (split(filepath_FI_TRANSP_shot,"."))[end] # Assume last part after final '.' is the file extension
    fileext_FI_TRANSP_shot = lowercase(fileext_FI_TRANSP_shot)
end
if !(fileext_FI=="h5" || fileext_FI=="jld2")
    println("Fast-ion distribution file format: ."*fileext_FI)
    error("Unknown fast-ion distribution file format. Please re-specify file and re-try.")
end

if !(fileext_thermal=="cdf" || fileext_thermal=="jld2" || fileext_thermal=="")
    println("Thermal distribution file format: ."*fileext_thermal)
    error("Unknown thermal distribution file format. Please re-specify file and re-try.")
end

## ---------------------------------------------------------------------------------------------
verbose && println("Loading the Julia packages and OWCF functions... ")
@everywhere begin
    cd(folderpath_OWCF) # Necessary to move all the workers to the correct folder
    using EFIT # For calculating magn½etic equilibrium quantities
    using Equilibrium # For loading flux function data, tokamak geometry data etc.
    using FileIO # To write/open files in general
    using GuidingCenterOrbits # For calculating guiding-center orbits
    using HDF5 # To write/open .hdf5 and .h5 files
    using Interpolations # To be able to create interpolation objects
    using JLD2 # To write/open .jld2-files (Julia files, basically)
    using NetCDF # To write/open .cdf files
    using OrbitTomography # This is what all this is about!
    using Printf # For print formatting
    using ProgressMeter # To display computational progress during parallel computations
    using PyCall # For using Python code in Julia
    using SparseArrays # To enable utilization of sparse matrices/vectors
    pushfirst!(PyVector(pyimport("sys")."path"), "") # To add the forward, transp_dists, transp_output and vcone Python modules (scripts in OWCF folder)
    include("misc/temp_n_dens.jl") # To be able to work with custom temperature and density profiles
    include("misc/species_func.jl") # To enable utilization of particle mass selection functions
    include("misc/availReacts.jl") # To examine fusion reaction and extract thermal and fast-ion species
    include("misc/diag2tokamak.jl") # To be able to determine a tokamak from a diagnostic label
    include("misc/rewriteReacts.jl") # To rewrite a fusion reaction from the A(b,c)D format to the A-b=c-D format
end

## ---------------------------------------------------------------------------------------------
verbose && println("Checking fusion reaction... ")
reaction_full = deepcopy(reaction) # Make a fully independent copy of the fusion reaction variable
@everywhere reaction_full = $reaction_full
reactants = full2reactsOnly(reaction) # Converts from 'a(b,c)d' format to 'a-b' format (reactants only)
@everywhere reaction = $reaction # Transfer to all external processes
emittedParticleHasCharge = false # By default, assume that the emitted particle 'c' in a(b,c)d does NOT have charge (is neutral)
RHEPWC = ["D-3He", "3He-D"] # RHEPWC means 'reaction has emitted particle with charge'
if reactants in RHEPWC # However, there are some fusion reactions which WILL produce an emitted particle with non-zero charge
    emittedParticleHasCharge = true
end

## ---------------------------------------------------------------------------------------------
verbose && println("Loading OWCF/extra/dependencies.jl... ")
@everywhere begin
    include("extra/dependencies.jl") # Load the functions in dependencies.jl.
end

## ---------------------------------------------------------------------------------------------
# Determine fast-ion and thermal species from inputs in start file
thermal_reactant, fast_reactant = getFusionReactants(reaction_full)
@everywhere thermal_reactant = $thermal_reactant
@everywhere fast_reactant = $fast_reactant

## ---------------------------------------------------------------------------------------------
# Error checks
if (nprocs()>1) && !distributed
    error(ErrorException("Number of processes greater than 1, but single-threaded computation specified. Please set distributed to true."))
end

try Int64(nE/nEbatch)
catch
    error("nE not divisible by nEbatch. Please correct and re-try.")
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
# Loading tokamak information and TRANSP RUN-ID from thermal .cdf file
if fileext_thermal=="cdf"
    verbose && println("Loading tokamak and TRANSP information... ")
    TRANSP_id = (split((split(filepath_thermal_distr,"."))[1],"/"))[end] # Assume part between last '/' and '.' is the TRANSP id
else
    TRANSP_id = ""
    if !(diagnostic_name=="")
        tokamak = diag2tokamak(diagnostic_name)
    end
end
@everywhere TRANSP_id = $TRANSP_id
timepoint_source = "UNKNOWN" # Initially assume we know nothing about where we got the timepoint data

## ---------------------------------------------------------------------------------------------
# Load the particle-space fast-ion distribution and corresponding grid arrays
if fileext_FI=="h5" || fileext_FI=="hdf5"
    verbose && println("Loading fast-ion distribution from .h5/.hdf5 file... ")
    F_ps, energy, pitch, R, z = h5to4D(filepath_FI_distr; rowmajor=h5_is_rowmajor, verbose = verbose)
elseif fileext_FI=="jld2"
    verbose && println("Loading fast-ion distribution from .jld2 file... ")
    myfile = jldopen(filepath_FI_distr,false,false,false,IOStream)
    if haskey(myfile,"F_ps")
        F_ps = myfile["F_ps"]
    elseif haskey(myfile,"F_EpRz")
        F_ps = myfile["F_EpRz"]
    else
        error("Fast-ion distribution .jld2 file did not have 'F_ps' nor 'F_EpRz'. Please correct and re-try.")
    end
    energy = myfile["energy"]
    pitch = myfile["pitch"]
    R = myfile["R"]
    if haskey(myfile,"z")
        z = myfile["z"]
    else
        z = myfile["Z"]
    end
    close(myfile)
else
    println("Fast-ion distribution file format: ."*fileext_FI)
    error("Unknown fast-ion distribution file format. Please re-specify file and re-try.")
end

## ---------------------------------------------------------------------------------------------
# Loading tokamak equilibrium
verbose && println("Loading tokamak equilibrium... ")
try
    global M; global wall; global jdotb; global timepoint; global timepoint_source
    M, wall = read_geqdsk(filepath_equil,clockwise_phi=false) # Assume counter-clockwise phi-direction
    jdotb = M.sigma # The sign of the dot product between the plasma current and the magnetic field

    # Extract timepoint information from .eqdsk/.geqdsk file
    eqdsk_array = split(filepath_equil,".")
    XX = (split(eqdsk_array[end-2],"-"))[end] # Assume format ...-XX.YYYY.eqdsk where XX are the seconds and YYYY are the decimals
    YYYY = eqdsk_array[end-1] # Assume format ...-XX.YYYY.eqdsk where XX are the seconds and YYYY are the decimals
    timepoint_source = "EQDSK"
    timepoint = XX*","*YYYY # (SOURCE, VALUE). Format XX,YYYY to avoid "." when including in filename of saved output
catch # Otherwise, assume magnetic equilibrium is a saved .jld2 file
    global M; global wall; global jdotb; global timepoint; local myfile; global timepoint_source
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

if fileext_FI_TRANSP_shot=="cdf" 
    # If the user has specified a TRANSP .cdf file with pertaining NUBEAM fast-ion distribution data...
    # Load the time, and overwrite timepoint. TRANSP time data superseeds .eqdsk time data
    TIME = round((ncread(filepath_FI_TRANSP_shot,"TIME"))[1],digits=4)
    TIME_array = split("$(TIME)",".") # Will be on format XX.YYYY
    XX = TIME_array[1]
    YYYY = TIME_array[2]
    timepoint_source, timepoint = "TRANSP", XX*","*YYYY # (SOURCE, VALUE). Format XX,YYYY to avoid "." when including in filename of saved output
end

@everywhere timepoint = $timepoint
@everywhere timepoint_source = $timepoint_source

## ---------------------------------------------------------------------------------------------
# Defining orbit grid vectors
if isfile(og_filepath)
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
    if !(fast_reactant===FI_species_loaded)
        @warn "Fast-ion species ("*FI_species_loaded*") loaded from orbit-grid .jld2 file does not match fast-ion species in specified fusion reaction ("*fast_reactant*"). Did you confuse the thermal species ("*thermal_reactant*") for the fast-ion species ("*fast_reactant*")?"
    end
    if !(fast_reactant===FI_species_loaded) && !(thermal_reactant===FI_species_loaded)
        error("Fast-ion species ("*FI_species_loaded*") loaded from orbit-grid .jld2 file does not match any reactants in specified fusion reaction ("*reaction*"). Please correct and re-try.")
    end
else
    verbose && println("Defining orbit grid vectors... ")
    E_array = range(Emin,stop=Emax,length=nE)
    pm_array = range(pm_min, stop=pm_max, length=npm)
    if isnothing(Rm_min) || isnothing(Rm_max)
        if inclPrideRockOrbs
            # 4/5 of the distance from the HFS wall to the magnetic axis is usually enough to capture all the Pride Rock orbits
            Rm_array = range((4*M.axis[1]+minimum(wall.r))/5, stop=maximum(wall.r), length=nRm)
        else
            Rm_array = range(M.axis[1], stop=maximum(wall.r), length=nRm)
        end
    else
        Rm_array = range(Rm_min, stop=Rm_max, length=nRm)
    end
end

if maximum(energy)<maximum(E_array)
    @warn "Specified orbit-space fast-ion energy range extends beyond the energy range of the loaded fast-ion distribution. Results may be inaccurate."
end
if minimum(energy)>minimum(E_array)
    @warn "Specified orbit-space fast-ion energy range extends below the energy range of the loaded fast-ion distribution. Results may be inaccurate."
end

## ---------------------------------------------------------------------------------------------
# Printing script info and inputs
println("")
println("-------------------------------------------------------ps2WF.jl---------------------------------------------------")
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
println("Diagnostic (model): "*diagnostic_filepath)
println("---> There will be $(length(range(Ed_min,stop=Ed_max-Ed_diff,step=Ed_diff))-2) diagnostic energy bins with")
println("------> Lower diagnostic energy bound: $(Ed_min) keV")
println("------> Upper diagnostic energy bound: $(Ed_max) keV")
println("")
println("Fusion reaction: "*reaction_full)
println("---> Fast-ion species: "*fast_reactant)
println("---> Fast-ion distribution (file): "*filepath_FI_distr)
if fileext_FI=="cdf"
    @warn ".cdf TRANSP file specified for fast-ion distribution. Sampling will take a relatively long time, compared to .h5 or .jld2 format."
end
if emittedParticleHasCharge
    println("---> The emitted "*getEmittedParticle(reaction_full)*" particle of the "*reaction_full*" reaction has non-zero charge!")
    println("---> The resulting energy distribution for "*getEmittedParticle(reaction_full)*" from the plasma as a whole will be computed.")
end
println("")
if !(filepath_thermal_distr=="")
    println("Bulk (thermal) plasma profile (data): "*filepath_thermal_distr)
    println("---> Bulk (thermal) plasma species: "*thermal_reactant)
else
    println("No bulk (thermal) plasma distribution file specified. Using default OWCF temperature+density profiles with: ")
    println("---> Thermal ($(thermal_reactant)) temperature on-axis: $(thermal_temp_axis) keV")
    println("------> Thermal ion ($(thermal_reactant)) temperature profile: $(thermal_temp_profile_type)")
    println("---> Thermal ($(thermal_reactant)) density on-axis: $(thermal_dens_axis) m^-3")
    println("------> Thermal ion ($(thermal_reactant)) density profile: $(thermal_dens_profile_type)")
end
println("")
println("Magnetic (B-field) equilibrium: "*filepath_equil)
println("---> |B| on-axis: $(round(norm(Equilibrium.Bfield(M,magnetic_axis(M)...)),digits=2)) T")
println("")
if distributed
    println("Parallel computing will be used with $(nprocs()) processes (1 main + $(nprocs()-1) workers).")
else
    println("Single-threaded computing... Good luck! (I hope you don't have too big an orbit grid)")
end
if !(diagnostic_filepath=="")
    println(""*diagnostic_name*" WF signal will be computed for a $(nE)x$(npm)x$(nRm) orbit-space grid with")
else
    println("No diagnostic_filepath specified. Spherical (4*pi sterradians) emission will be assumed.")
end
println("---> Energy: [$(minimum(E_array)),$(maximum(E_array))] keV")
println("---> Pitch maximum: [$(minimum(pm_array)),$(maximum(pm_array))]")
println("---> Radius maximum: [$(minimum(Rm_array)),$(maximum(Rm_array))] m")
println("")
println("The WF computation will be divided into $(Int64(nE/nEbatch)) sub-WF-computations.")
println("")
if calcWFOs
    println("WF densities, W and F (split into orbit classes) as functions of E, pm and Rm will be computed.")
    println("")
end
println("Results will be saved to: ")
sFLR = include_FLR_effects ? "FLR_" : ""
println("---> "*folderpath_o*"ps2WF_results_$(sFLR)"*tokamak*"_"*TRANSP_id*"_"*diagnostic_name*"_"*pretty2scpok(reaction_full)*"_$(nE)x$(npm)x$(nRm).jld2")
println("")
println("Please remove previously saved files with the same file name (if any) prior to script completion!")
println("")
println("If you would like to change any settings, please edit the ps2WF start file.")
println("Written by Henrik Järleblad. Last maintained 2025-07-11.")
println("------------------------------------------------------------------------------------------------------------------")
println("")

## ---------------------------------------------------------------------------------------------
# Python code (essentially calling the NES forward model black box)
println("Loading Python modules... ")
@everywhere begin
    py"""
    import numpy as np
    import forward
    import spec
    import vcone
    """
end

## ---------------------------------------------------------------------------------------------
# Helper function (to make the code more transparent)
println("Loading helper functions... ")
@everywhere begin
    include("helper/calcOrbSpec.jl")
end

## ---------------------------------------------------------------------------------------------
verbose && println("Initializing forward model and loading bulk (thermal) plasma distribution data (if specified)... ")
@everywhere begin
    py"""
    # The '$' in front of many Python variables means that the variable is defined in Julia, not in Python.
    test_thermal_particle = spec.Particle($thermal_reactant) # Check so that thermal species is available in DRESS code
    reaction_full = $reaction_full
    thermal_reactant = $thermal_reactant
    thermal_dist = ""

    # Load thermal and/or fast-ion TRANSP data
    if ($fileext_thermal=="cdf") and (not $filepath_FI_TRANSP_shot==""):
        import transp_output
        import transp_dists
        tr_out = transp_output.TranspOutput($TRANSP_id, step=1, # The TRANSP_id and the step number (always 1 for 1 fbm_file)
                                            out_file=$filepath_thermal_distr, # The thermal distribution file
                                            fbm_files=[$filepath_FI_TRANSP_shot]) # Load the TRANSP shot file
        thermal_dist = transp_dists.Thermal(tr_out, ion=thermal_reactant) # Then load the thermal ion distribution from that TRANSP object

    forwardmodel = forward.Forward($diagnostic_filepath, reaction_full, thermal_dist) # Pre-initialize the forward model
    $verbose and print("---> Initializing diagnostic measurement bins... ")
    Ed_bin_edges = np.arange($Ed_min,$Ed_max,$Ed_diff) # diagnostic spectrum bin edges (keV or m/s)
    if len(Ed_bin_edges)==1: # Make sure that there are at least one lower and one upper bin edge
        dEd = (($Ed_max)-($Ed_min))/10
        Ed_bin_edges = np.arange($Ed_min,($Ed_max)+dEd,$Ed_diff)
    Ed_vals = 0.5*(Ed_bin_edges[1:] + Ed_bin_edges[:-1]) # bin centers (keV or m/s)
    nEd = len(Ed_vals)
    """
    nEd = py"nEd"; Ed_array = vec(py"Ed_vals")
end

if fileext_thermal=="jld2"
    verbose && println("Loading thermal temperature and density profiles from .jld2 file... ")
    myfile = jldopen(filepath_thermal_distr,false,false,false,IOStream)
    thermal_temp_array = myfile["thermal_temp"]
    thermal_dens_array = myfile["thermal_dens"]
    ρ_pol_array = myfile["rho_pol"]
    close(myfile)
    @everywhere thermal_temp_array = $thermal_temp_array
    @everywhere thermal_dens_array = $thermal_dens_array
    @everywhere ρ_pol_array = $ρ_pol_array
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
            if timepoint_source=="UNKNOWN"
                error("$(filepath_thermal_distr) specified as bulk (thermal) plasma distribution, but no valid timepoint could be inferred ($(timepoint) seconds, $(timepoint_source)). Possible solutions: (1) Manually specify timepoint in start file (2) Use an .eqdsk file with timepoint data included in file name")
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
            verbose && println("---> Flat thermal ion temperature ($(thermal_temp_axis) keV) profile assumed")
        end
        if thermal_dens_profile_type==:DEFAULT
            thermal_dens_etp = x-> getAnalyticalDens(thermal_dens_axis,x)
            verbose && println("---> Default OWCF density profile with n_axis=$(thermal_dens_axis) m^-3 will be used as thermal ion density data (more info in OWCF/misc/ folder)")
        else # Otherwise, just use a flat bulk (thermal) plasma density profile
            thermal_dens_etp = x-> thermal_dens_axis
            # Don't print if projected velocities are to be computed. Then, bulk (thermal) temp+dens are irrelevant
            verbose && println("---> Flat thermal ion density ($(thermal_dens_axis) m^-3) profile assumed")
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

## ---------------------------------------------------------------------------------------------
# The partitioning of the energy array
verbose && println("Partitioning the orbit-space energy grid... ")
E_chunks = reshape(E_array, (nEbatch,div(length(E_array),nEbatch)))

## ---------------------------------------------------------------------------------------------
# Load already computed progress from file
# Load i, S_WF
if isfile("ps2WF_progress_$(sFLR)"*tokamak*"_"*TRANSP_id*"_at"*timepoint*"s_"*diagnostic_name*"_$(nE)x$(npm)x$(nRm).jld2")
    verbose && println("Found ps2WF progress file! Loading... ")
    myfile = jldopen("ps2WF_progress_$(sFLR)"*tokamak*"_"*TRANSP_id*"_at"*timepoint*"s_"*diagnostic_name*"_$(nE)x$(npm)x$(nRm).jld2",false,false,false,IOStream)
    i0 = myfile["i"]
    S_WF = myfile["S_WF"]
    nfast_chunks = myfile["nfast_chunks"]
    if calcWFOs
        WO_E = myfile["WO_E"]
        WO_pm = myfile["WO_pm"]
        WO_Rm = myfile["WO_Rm"]

        FO_E = myfile["FO_E"]
        FO_pm = myfile["FO_pm"]
        FO_Rm = myfile["FO_Rm"]

        WFO_E = myfile["WFO_E"]
        WFO_pm = myfile["WFO_pm"]
        WFO_Rm = myfile["WFO_Rm"]

        NO_E = myfile["NO_E"]
        NO_pm = myfile["NO_pm"]
        NO_Rm = myfile["NO_Rm"]
    end
    close(myfile)
    # Assuming E_chunks and everything else is the same.
else
    S_WF = zeros(length(Ed_array))
    nfast_chunks = zeros(length(E_chunks))
    if calcWFOs
        WO_E = zeros(length(Ed_array),length(E_array),6) # There are 6 types of valid orbits. 1 = stagnation, 2 = trapped, 3 = co-passing, 4 = counter-passing, 5 = potato, 6 = counter-stagnation
        WO_pm = zeros(length(Ed_array),length(pm_array),6) # There are 6 types of valid orbits
        WO_Rm = zeros(length(Ed_array),length(Rm_array),6) # There are 6 types of valid orbits

        FO_E = zeros(length(E_array),6) # There are 6 types of valid orbits. 1 = stagnation, 2 = trapped, 3 = co-passing, 4 = counter-passing, 5 = potato, 6 = counter-stagnation
        FO_pm = zeros(length(pm_array),6) # There are 6 types of valid orbits
        FO_Rm = zeros(length(Rm_array),6) # There are 6 types of valid orbits

        WFO_E = zeros(length(Ed_array),length(E_array),6) # There are 6 types of valid orbits. 1 = stagnation, 2 = trapped, 3 = co-passing, 4 = counter-passing, 5 = potato, 6 = counter-stagnation
        WFO_pm = zeros(length(Ed_array),length(pm_array),6) # There are 6 types of valid orbits
        WFO_Rm = zeros(length(Ed_array),length(Rm_array),6) # There are 6 types of valid orbits

        NO_E = zeros(length(E_array),6) # There are 6 types of valid orbits. 1 = stagnation, 2 = trapped, 3 = co-passing, 4 = counter-passing, 5 = potato, 6 = counter-stagnation
        NO_pm = zeros(length(pm_array),6) # There are 6 types of valid orbits
        NO_Rm = zeros(length(Rm_array),6) # There are 6 types of valid orbits
    end
    i0 = 1
end

## ---------------------------------------------------------------------------------------------
# Send all necessary variables to external processes (like E_array etc)
# Also, send nE instead of length(E_array) etc, to minimize memory usage
# using @everywhere x = $x
nE = length(E_array)
npm = length(pm_array)
nRm = length(Rm_array)
nEd = length(Ed_array)
dE = diff(E_array)[1] # Assume equidistant
dpm = diff(pm_array)[1] # Assume equidistant
dRm = diff(Rm_array)[1] # Assume equidistant
@everywhere nE = $nE
@everywhere npm = $npm
@everywhere nRm = $nRm
@everywhere nEd = $nEd
@everywhere dE = $dE
@everywhere dpm = $dpm
@everywhere dRm = $dRm

## ---------------------------------------------------------------------------------------------
# Define WFO function (in case calcWFO===true)
# This makes cluster jobs as efficient as possible.
"""

    calcWFO()

Just a helper function to the big for-loop below. E_array should be total array, not E_chunk.
"""
function calcWFO(M::AbstractEquilibrium, og::OrbitGrid, i::Int, F_os_chunk::AbstractVector, E_array::AbstractVector, Ed_array::AbstractVector, nEbatch::Int, og_orbs::Vector{Orbit{Float64, EPRCoordinate{Float64}}}, W_chunk::AbstractMatrix;verbose=false)

    F_os_3D_chunk = OWCF_map_orbits(og, Vector(F_os_chunk)) # Transform F_os_chunk from a 1D Vector with units [ion] to a 3D array with units [ion/(orbit-space volume)]
    subs = CartesianIndices(size(F_os_3D_chunk))

    non_zero_orbInds = findall(x -> x != 0.0, og.orbit_index)

    verbose && println("calcWFO(): Computing the @distributed loop... ")
    WFO_res = @distributed (+) for oii in non_zero_orbInds
        # There are 6 possible valid basic orbit types
        WFO_res_sub = [zeros(nE,6),zeros(npm,6),zeros(nRm,6),zeros(nEd,nE,6),zeros(nEd,npm,6),zeros(nEd,nRm,6),zeros(nEd,nE,6),zeros(nEd,npm,6),zeros(nEd,nRm,6),zeros(nE,6),zeros(npm,6),zeros(nRm,6)]
        oi = (og.orbit_index)[oii]
        Ei = (i-1)*nEbatch+subs[oii][1]
        pmi = subs[oii][2]
        Rmi = subs[oii][3]
        oint = class2int(M, og_orbs[oi];plot=false)
        (WFO_res_sub[1])[Ei,oint] += F_os_3D_chunk[subs[oii]] * dpm * dRm
        (WFO_res_sub[2])[pmi,oint] += F_os_3D_chunk[subs[oii]] * dE * dRm
        (WFO_res_sub[3])[Rmi,oint] += F_os_3D_chunk[subs[oii]] * dE * dpm

        for Edi=1:nEd
            W_chunk_3D = OWCF_map_orbits(og, Vector(W_chunk[Edi,:]); weights=true) # weights=true to not have inv(dE*dpm*Rm) automatically included
            WF_chunk_3D = W_chunk_3D .* F_os_3D_chunk

            (WFO_res_sub[4])[Edi, Ei, oint] += W_chunk_3D[subs[oii]]
            (WFO_res_sub[5])[Edi, pmi, oint] += W_chunk_3D[subs[oii]]
            (WFO_res_sub[6])[Edi, Rmi, oint] += W_chunk_3D[subs[oii]]

            (WFO_res_sub[7])[Edi, Ei, oint] += WF_chunk_3D[subs[oii]] * dpm * dRm
            (WFO_res_sub[8])[Edi, pmi, oint] += WF_chunk_3D[subs[oii]] * dE * dRm
            (WFO_res_sub[9])[Edi, Rmi, oint] += WF_chunk_3D[subs[oii]] * dE * dpm
        end

        (WFO_res_sub[10])[Ei,oint] += 1.0
        (WFO_res_sub[11])[pmi,oint] += 1.0
        (WFO_res_sub[12])[Rmi,oint] += 1.0

        WFO_res_sub # Declare for @distributed (+) reduction
    end

    return WFO_res
end

## ---------------------------------------------------------------------------------------------
# Save the standard output channel to be able to switch back to it, after surpressing outputs from orbit_grid()
oldstd = stdout

## ---------------------------------------------------------------------------------------------
# Computing the total number of fast ions (nfast_orig) in the fast-ion distribution data, for reference
verbose && println("Computing total number of fast ions (N_FI) in the plasma, assuming N_FI = ∫ f(E,p,R,z) 2*pi*R dEdpdRdz... ")
dE4D, dp4D, dR4D, dz4D = get4DDiffs(energy, pitch, R, z)
fr = F_ps .* reshape(R,(1,1,length(R),1))
nfast_orig = sum((2*pi) .* fr .* dE4D .*dp4D .*dR4D .*dz4D)

## ---------------------------------------------------------------------------------------------
# The actual, very long, loop over energies... !
for i = i0:size(E_chunks,2)
    verbose && println("---------- Computing WF-chunk $(i) of $(size(E_chunks,2)) ----------")
    E_chunk = E_chunks[:,i]
    verbose && println("Computing orbit grid $(i) of $(size(E_chunks,2))... ")
    redirect_stdout(devnull) # Surpress progress bar here. Would mess up the log file
    og_orbs, og = OrbitTomography.orbit_grid(M, E_chunk, pm_array, Rm_array; wall=wall, amu=getSpeciesAmu(fast_reactant), q=getSpeciesEcu(fast_reactant), extra_kw_args...)
    redirect_stdout(oldstd) # recover original stdout. Print to stdout again
    verbose && println(og)

    # Pick out the correct F_ps chunk
    verbose && println("Interpolating F_ps $(i) of $(size(E_chunks,2))... ")
    F_ps_chunk = interpFps(F_ps, energy, pitch, R, z, E_chunk, pitch, R, z; debug=debugging)

    verbose && println("Extrema(energy): $(extrema(energy))")
    verbose && println("Extrema(E_chunk): $(extrema(E_chunk))")


    verbose && println("Performing ps2os $(i) of $(size(E_chunks,2))... ")
    F_os_raw_chunk, class_distr, nfast_chunk = ps2os(M, wall, F_ps_chunk, E_chunk, pitch, R, z, og; 
                                        numOsamples = Int64(12*length(og_orbs)), distributed = distributed, nbatch = 100_000, 
                                        saveProgress=false, FI_species=fast_reactant, extra_kw_args...)

    wasteOtime = false
    if sum(F_os_raw_chunk)==0 || nfast_chunk==0
        wasteOtime = true
    end

    verbose && print("sum(F_os_raw_chunk): $(Int64(sum(F_os_raw_chunk)))         "); nfast_chunk_s = @sprintf "%.2E" nfast_chunk
    verbose && println("nfast_chunk: $(nfast_chunk_s) ($(round(100*nfast_chunk/nfast_orig,sigdigits=2))% of total number of fast ions in the plasma)")
    (verbose && wasteOtime && !calcWFOs) && println("Chunk is a waste of time. Skipping orbit weights computation... ")

    if !wasteOtime || calcWFOs
        F_os_chunk = (nfast_chunk/(sum(F_os_raw_chunk)==0.0 ? 1.0 : sum(F_os_raw_chunk))) .*F_os_raw_chunk # If sum(F_os_raw_chunk) is zero, divide by 1.0 instead (to avoid divide by zero).
    end

    # The OW computation part
    if !wasteOtime || calcWFOs
        verbose && println("Calculating W-chunk $(i) of $(size(E_chunks,2))... ")
        W_chunk = calcOrbSpecs(og_orbs, 1.0 .*ones(size(og_orbs)), compute_measurements; distributed=distributed)
    end
    
    if sum(isnan.(W_chunk))!=0.0 && debugging
        verbose && println("Encountered NaNs in weight function! Saving debugging info... ")
        local myfile = jldopen(folderpath_o*"ps2WF_debug_$(nE)x$(npm)x$(nRm).jld2",true,true,false,IOStream)
        sumOnans = sum(isnan.(W_chunk))
        write(myfile,"W_chunk",W_chunk)
        write(myfile,"i",i)
        write(myfile,"E_chunk",E_chunk)
        write(myfile,"sumOnans",sumOnans)
        close(myfile)
        error("ENCOUNTERED NaNs IN WEIGHT FUNCTION CHUNK. THROWING ERROR AND TERMINATING COMPUTATION TO STOP WASTING CPU RESOURCES. PLEASE EXAMINE THE DEBUG FILE SAVED IN THE folderpath_o FOLDER.")
    end

    if calcWFOs
        verbose && println("Computing WF densities $(i) of $(size(E_chunks,2))...  ")
        WFO_res = calcWFO(M, og, i, F_os_chunk, E_array, Ed_array, nEbatch, og_orbs, W_chunk;verbose=verbose)
    end

    # Save to actual result variables
    if calcWFOs
        verbose && println("Adding computed WFO quantities to global quantities...")
        global FO_E += WFO_res[1]
        global FO_pm += WFO_res[2]
        global FO_Rm += WFO_res[3]
        global WO_E += WFO_res[4]
        global WO_pm += WFO_res[5]
        global WO_Rm += WFO_res[6]
        global WFO_E += WFO_res[7]
        global WFO_pm += WFO_res[8]
        global WFO_Rm += WFO_res[9]
        global NO_E += WFO_res[10]
        global NO_pm += WFO_res[11]
        global NO_Rm += WFO_res[12]
    end

    # The S_WF computation part
    verbose && println("Calculating S_WF signal $(i) of $(size(E_chunks,2))... ")
    if !wasteOtime
        S_WF .+= W_chunk*F_os_chunk # Please note! F_os_chunk is in units of [ion]
    else
        S_WF .+= zeros(size(Ed_array))
    end
    nfast_chunks[i] = nfast_chunk

    # Save current progress, including i, F_os_chunk, WF_chunk and S_WF
    verbose && println("Force-removing ps2WF_progress_$(sFLR)"*tokamak*"_"*TRANSP_id*"_at"*timepoint*"s_"*diagnostic_name*"_$(nE)x$(npm)x$(nRm).jld2 file... ")
    rm("ps2WF_progress_$(sFLR)"*tokamak*"_"*TRANSP_id*"_at"*timepoint*"s_"*diagnostic_name*"_$(nE)x$(npm)x$(nRm).jld2",force=true)
    verbose && println("Creating new ps2WF_progress_$(sFLR)"*tokamak*"_"*TRANSP_id*"_at"*timepoint*"s_"*diagnostic_name*"_$(nE)x$(npm)x$(nRm).jld2 file... ")
    myfile = jldopen("ps2WF_progress_$(sFLR)"*tokamak*"_"*TRANSP_id*"_at"*timepoint*"s_"*diagnostic_name*"_$(nE)x$(npm)x$(nRm).jld2",true,true,false,IOStream)

    verbose && println("Saving progress... ")
    write(myfile,"i",i)
    write(myfile,"progress",round(i/size(E_chunks,2),sigdigits=3)) # Progress. 0 means no progress, 1 means fully completed, 0.5 means halfway completed etc
    write(myfile,"F_os_chunk", F_os_chunk)
    write(myfile,"S_WF",S_WF)
    write(myfile,"nfast_chunks",nfast_chunks)
    if calcWFOs
        write(myfile,"WO_E",WO_E)
        write(myfile,"WO_pm",WO_pm)
        write(myfile,"WO_Rm",WO_Rm)

        write(myfile,"FO_E",FO_E)
        write(myfile,"FO_pm",FO_pm)
        write(myfile,"FO_Rm",FO_Rm)

        write(myfile,"WFO_E",WFO_E)
        write(myfile,"WFO_pm",WFO_pm)
        write(myfile,"WFO_Rm",WFO_Rm)

        write(myfile,"NO_E",NO_E)
        write(myfile,"NO_pm",NO_pm)
        write(myfile,"NO_Rm",NO_Rm)
    end
    close(myfile)
end
verbose && println("Done!")

if calcWFOs
    # Normalize the WOs to account for orbit-space measure
    verbose && println("Normalizing orbit weight function splits by orbit type... ")
    WNO_E = zeros(size(WO_E))
    WNO_pm = zeros(size(WO_pm))
    WNO_Rm = zeros(size(WO_Rm))
    for iEd=1:nEd
        WNO_E[iEd,:,:] = WO_E[iEd,:,:] ./ NO_E
        WNO_pm[iEd,:,:] = WO_pm[iEd,:,:] ./ NO_pm
        WNO_Rm[iEd,:,:] = WO_Rm[iEd,:,:] ./ NO_Rm
    end
    WNO_E = map(x-> isnan(x) ? 0.0 : x, WNO_E) # We know that NaNs would correspond to zero weight anyway (e.g. w(stagnation) / N(stagnation) can be NaN only if both w(stagnation) and N(stagnation) were 0.0)
    WNO_pm = map(x-> isnan(x) ? 0.0 : x, WNO_pm) # We know that NaNs would correspond to zero weight anyway (e.g. w(stagnation) / N(stagnation) can be NaN only if both w(stagnation) and N(stagnation) were 0.0)
    WNO_Rm = map(x-> isnan(x) ? 0.0 : x, WNO_Rm) # We know that NaNs would correspond to zero weight anyway (e.g. w(stagnation) / N(stagnation) can be NaN only if both w(stagnation) and N(stagnation) were 0.0)
end

# Removing progress files and saving results... 
verbose && println("Force-removing ps2WF_progress_$(sFLR)"*tokamak*"_"*TRANSP_id*"_at"*timepoint*"s_"*diagnostic_name*"_$(nE)x$(npm)x$(nRm).jld2 file... ")
rm("ps2WF_progress_$(sFLR)"*tokamak*"_"*TRANSP_id*"_at"*timepoint*"s_"*diagnostic_name*"_$(nE)x$(npm)x$(nRm).jld2",force=true)
verbose && println("Saving to results file... ")

global filepath_output_orig = folderpath_o*"ps2WF_results_$(sFLR)"*tokamak*"_"*TRANSP_id*"_at"*timepoint*"s_"*diagnostic_name*"_"*pretty2scpok(reaction_full)*"_$(nE)x$(npm)x$(nRm)"
global filepath_output = deepcopy(filepath_output_orig)
global count = 1
while isfile(filepath_output*".jld2") # To take care of not overwriting files. Add _(1), _(2) etc
    global filepath_output = filepath_output_orig*"_($(Int64(count)))"
    global count += 1 # global scope, to surpress warnings
end
global filepath_output = filepath_output*".jld2"
myfile = jldopen(filepath_output,true,true,false,IOStream)
write(myfile,"S_WF",S_WF)
write(myfile,"E_array",E_array)
write(myfile,"pm_array",pm_array)
write(myfile,"Rm_array",Rm_array)
write(myfile,"Ed_array",Ed_array)
if calcWFOs
    write(myfile,"FO_E",FO_E)
    write(myfile,"FO_pm",FO_pm)
    write(myfile,"FO_Rm",FO_Rm)

    write(myfile,"WFO_E",WFO_E)
    write(myfile,"WFO_pm",WFO_pm)
    write(myfile,"WFO_Rm",WFO_Rm)

    write(myfile,"WO_E",WO_E)
    write(myfile,"WO_pm",WO_pm)
    write(myfile,"WO_Rm",WO_Rm)

    write(myfile,"WNO_E",WNO_E)
    write(myfile,"WNO_pm",WNO_pm)
    write(myfile,"WNO_Rm",WNO_Rm)

    write(myfile,"NO_E",NO_E)
    write(myfile,"NO_pm",NO_pm)
    write(myfile,"NO_Rm",NO_Rm)
end
write(myfile,"filepath_FI_distr",filepath_FI_distr)
write(myfile,"filepath_thermal_distr",filepath_thermal_distr)
write(myfile,"nfast_orig",nfast_orig)
write(myfile,"nfast_chunks",nfast_chunks)
write(myfile,"extra_kw_args",extra_kw_args)
close(myfile)

println("~~~~ps2WF.jl completed!~~~~")
