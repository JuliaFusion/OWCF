#########################################  calc4DWeights.jl #########################################

#### Description:
# This script computes 4D weight functions in (E,p,R,z). They can also be saved in (vpara, vperp, R, z)
# should the user wish. In the current version of the OWCF, the calc4DWeights.jl script is taylored 
# to utilize the DRESS code (J. Eriksson et al, CPC, 199, 40-46, 2016) to compute the weight functions. 
# In future versions, it can be easily modified to instead save the weighted (E,p,R,z) points in a file 
# readable by e.g. FIDASIM. The calc4DWeights.jl script computes the weight functions on a grid in 
# (E,p,R,z) space. The weights on the corresponding (vpara, vperp, R, z) grid will be computed and
# saved if the user has set the 'saveVparaVperpWeights' input variable to true.
# 
# The DRESS code is written in Python and calc4DWeights.jl utilizes the DRESS code via the Julia-Python
# interface package PyCall.jl
#
# The calc4DWeights.jl script computes weight functions for an equidistant, rectangular grid 
# in (E,p,R,z) space. Irregular/non-equidistant grid points are currently not supported.
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
# The 'thermal_temp_profile_type' and 'thermal_dens_profile_type' input variables can be used 
# to specify a flat or OWCF default profile. Please see OWCF/templates/start_calc4DW_template.jl and 
# OWCF/misc/temp_n_dens.jl for more info.
#
# Please see the start_calc4DW_template.jl file for further input information.

#### Inputs (Units given when defined in script)
# Given via input file start_calc4DW_template.jl, for example. 'template' should be replaced by whatever.

#### Outputs
# -

#### Saved files
# EpRzWeights_[FLR]_[tokamak]_[TRANSP_id]_at[timepoint]s_[diagnostic_name]_[nEd]x[nE]x[np]x[nR]x[nz].jld2 - If iiimax == 1
# EpRzWeights_[FLR]_[tokamak]_[TRANSP_id]_at[timepoint]s_[diagnostic_name]_$(i).jld2 - If iiimax != 1
# EpRzWeights_[FLR]_[tokamak]_[TRANSP_id]_at[timepoint]s_[diagnostic_name].jld2 - If iiimax != 1 && iii_average
# Regardless of saved file name, the saved file will have the fields:
#   VAL - The non-zero values of the (E,p,R,z) weight matrix (actually a 5D quantity, but arranged as a 2D matrix). - Vector{Float64}
#   ROW - The row indices of the VAL elements. The ROW elements point to a specific element in Ed_array, since the number of rows of the weight matrix equals length(Ed_array) - Vector{Int64}
#   COL - The column indices of the VAL elements. The COL elements point to a specific index in EpRz_coords (see below), since the number of columns of the weight matrix equals the number of phase-space points - Vector{Int64}
#   m_W - The total number of rows of the (E,p,R,z) weight matrix. Saved to ensure the SparseArrays.jl API can correctly re-create the weight matrix, since ROW might not contain the index of the last row of the weight matrix - Int64
#   n_W - The total number of columns of the (E,p,R,z) weight matrix. -||-, since COL might not contain the index of the last column of the weight matrix - Int64
#   EpRz_coords - The (E,p,R,z) indices and grid points for each column of the weight matrix. Element with index i in EpRz_coords[:] corresponds to column i in the weight matrix.
#                Each element in EpRz_coords is a tuple with four elements, and each element in that tuple is a tuple with two elements. An example of an element in EpRz_coords looks like:
#                ((9,45.3),(4,-0.35),(6,3.23),(3,0.12)) 
#                where 9 denotes that 45.3 keV is the 9th element in E_array, 4 denotes that -0.35 is the 4th element in p_array, 6 denotes that 3.23 meters is the 6th element in R_array
#                and 3 denotes that 0.12 meters is the 3rd element in z_array. In this way, each column of the weight matrix can be mapped to a specific (E,p,R,z) point, both in terms of indices and values,
#                via EpRz_coords[:]. - Array{NTuple{4, Tuple{Int64, Float64}}, 4}
#   E_array - The fast-ion energy grid array used for (E,p,R,z) space - Vector{Float64}
#   p_array - The fast-ion p grid array used for (E,p,R,z) space - Vector{Float64}
#   R_array - The fast-ion R grid array used for (E,p,R,z) space - Vector{Float64}
#   z_array - The fast-ion z grid array used for (E,p,R,z) space - Vector{Float64}
#   Ed_array - The diagnostic energy bin centers - Vector{Float64}
#   reaction - The nuclear fusion reaction for which the weights are computed - String
#   filepath_thermal_distr - The filepath of the thermal species distribution. For reference - String
# If weight functions are computed using projected velocities, the saved file will also have the key
#   projVel - True if weight functions were computed using projected velocities. - Bool
# If saveVparaVperpWeights was set to true, the output file will also contain
#   VAL_vel - The non-zero values of the (vpara,vperp,R,z) weight matrix (actually a 5D quantity, but arranged as a 2D matrix). - Vector{Float64}
#   ROW_vel - The row indices of the VAL_vel elements. The ROW_vel elements point to a specific element in Ed_array, since the number of rows of the weight matrix equals length(Ed_array) - Vector{Int64}
#   COL_vel - The column indices of the VAL_vel elements. The COL_vel elements point to a specific index in VparVperRz_coords (see below), since the number of columns of the weight matrix equals the number of phase-space points - Vector{Int64}
#   m_W_vel - The total number of rows of the (vpara,vperp,R,z) weight matrix. Saved to ensure the SparseArrays.jl API can correctly re-create the weight matrix, since ROW_vel might not contain the index of the last row of the weight matrix - Int64
#   n_W_vel - The total number of columns of the (vpara,vperp,R,z) weight matrix. -||-, since COL_vel might not contain the index of the last column of the weight matrix - Int64
#   VparVperRz_coords - The (vpara,vperp,R,z) indices and grid points for each column of the weight matrix. Element with index i in VparVperRz_coords[:] corresponds to column i in the weight matrix.
#                Each element in VparVperRz_coords is a tuple with four elements, and each element in that tuple is a tuple with two elements. An example of an element in VparVperRz_coords looks like:
#                ((7,-123000.27),(3,68000.56),(4,3.63),(5,0.27)) 
#                where 7 denotes that -123000.27 m/s is the 7th element in vpara_array, 3 denotes that 68000.56 m/s is the 3rd element in vperp_array, 4 denotes that 3.63 meters is the 4th element in R_array
#                and 5 denotes that 0.27 meters is the 5th element in z_array. In this way, each column of the weight matrix can be mapped to a specific (vpara,vperp,R,z) point, both in terms of indices and values,
#                via VparVperRz_coords[:]. - Array{NTuple{4, Tuple{Int64, Float64}}, 4}
#   vpara_array - The fast-ion vpara grid array used for (vpara, vperp) space - Vector{Float64}
#   vperp_array - The fast-ion vpara grid array used for (vpara, vperp) space - Vector{Float64}
# If an instrumental response function has been specified in the start file, the output file will also contain
#   VAL_raw - The non-zero values of the raw (without instrumental response) (E,p,R,z) weight matrix (actually a 5D quantity, but arranged as a 2D matrix). - Vector{Float64}
#   ROW_raw - The row indices of the VAL_raw elements. The ROW_raw elements point to a specific element in Ed_array_raw, since the number of rows of the weight matrix equals length(Ed_array_raw) - Vector{Int64}
#   COL_raw - The column indices of the VAL_raw elements. The COL_raw elements point to a specific index in EpRz_coords (see above), since the number of columns of the weight matrix equals the number of phase-space points - Vector{Int64}
#   m_W_raw - The total number of rows of the raw (E,p,R,z) weight matrix. Saved to ensure the SparseArrays.jl API can correctly re-create the weight matrix, since ROW_raw might not contain the index of the last row of the weight matrix - Int64
#   n_W_raw - The total number of columns of the raw (E,p,R,z) weight matrix. -||-, since COL_raw might not contain the index of the last column of the weight matrix - Int64
#   Ed_array_raw - The diagnostic bin centers, if there was no instrumental response - Vector{Float64}
#   instrumental_response_input - The input grid points for the instrumental response model - Vector{Float64}
#   instrumental_response_output - The output grid points for the instrumental response model - Vector{Float64}
#   instrumental_response_matrix - The instrumental response matrix, i.e. the model - Matrix{Float64}
# If, in addition, saveVparaVperpWeights was also set to true, the output file will also contain 
#   VAL_vel_raw - The non-zero values of the raw (without instrumental response) (vpara,vperp,R,z) weight matrix (actually a 5D quantity, but arranged as a 2D matrix). - Vector{Float64}
#   ROW_vel_raw - The row indices of the VAL_vel_raw elements. The ROW_vel_raw elements point to a specific element in Ed_array, since the number of rows of the weight matrix equals length(Ed_array) - Vector{Int64}
#   COL_vel_raw - The column indices of the VAL_vel_raw elements. The COL_vel_raw elements point to a specific index in VparVperRz_coords (see above), since the number of columns of the weight matrix equals the number of phase-space points - Vector{Int64}
#   m_W_vel_raw - The total number of rows of the raw (vpara,vperp,R,z) weight matrix. Saved to ensure the SparseArrays.jl API can correctly re-create the weight matrix, since ROW_vel_raw might not contain the index of the last row of the weight matrix - Int64
#   n_W_vel_raw - The total number of columns of the raw (vpara,vperp,R,z) weight matrix. -||-, since COL_vel_raw might not contain the index of the last column of the weight matrix - Int64

# To rebuild the weight matrices from the VAL, ROW and COL vectors, one does the following:
# *Load the SparseArrays.jl package
# *Load the calc4DWeights.jl output file
# W = dropzeros(sparse(append!(ROW,m_W),append!(COL,n_W),append!(VAL,0.0)))
# Woi = W[134,:]
# Woi_4D = zeros(length(E_array),length(p_array),length(R_array),length(z_array))
# for (i,c) in enumerate(EpRz_coords[:])
#     Woi_4D[c[1][1],c[2][1],c[3][1],c[4][1]] = Woi[i]
# end
# W will then be your 2D weight matrix, Woi will be the 134th row of W and Woi_4D will be the inflated 4D (E,p,R,z) array of Woi.
# The same loading process can be repeated with the _vel, _raw and _vel_raw data.

### Other
# Please note that the diagnostic energy grid will be created as bin centers.
# That is, the first diagnostic energy grid value will be (Ed_min+Ed_diff/2) and so on.
#
# WARNING! Please note that the output files of calc4DWeights.jl will be LARGE. This is due to the relatively high dimensionality
# of the weight functions.

# Script written by Henrik Järleblad. Last maintained 2025-09-01.
################################################################################################

## ---------------------------------------------------------------------------------------------
verbose && println("Loading Julia packages... ")
@everywhere begin
    cd(folderpath_OWCF) # Necessary to move all the workers to the correct folder
    using PyCall # For using Python code in Julia
    using Printf # To be able to print specific formats
    plot_results && (using Plots)
    include("misc/availReacts.jl") # To examine fusion reaction and extract thermal and fast-ion species
    include("extra/dependencies.jl") # To enable usage of OWCF functions library
    include("misc/rewriteReacts.jl") # To rewrite a fusion reaction from e.g. the a(b,c)d format to the a-b=c-d format
    include("misc/species_func.jl") # To map particle species labels to particle mass etc
    pushfirst!(PyVector(pyimport("sys")."path"), "") # To add the forward, transp_dists, transp_output and vcone modules (scripts in current path)
end

## ---------------------------------------------------------------------------------------------
# Fusion reaction, diagnostic and particle-related checks
if getReactionForm(reaction)==1 # If no excited energy state for the emitted particle (in the case it is an atomic nucleus) has been specified...
    verbose && println("No energy state specified for the emitted particle $(getEmittedParticle(reaction)). Assuming ground state (GS), if relevant... ")
    reaction *= "-GS"
end
if !reactionIsAvailable(reaction)
    error("Fusion reaction $(reaction) is not yet available in the OWCF. The following reactions are available: $(OWCF_AVAILABLE_FUSION_REACTIONS). For projected-velocity computations, the following particle species are available: $(OWCF_SPECIES). Please correct and re-try.")
end
@everywhere reaction = $reaction # Copy the reaction variable to all external (CPU) processes

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

if !(lowercase(getEmittedParticleEnergyLevel(reaction))=="gs") && (diagnostic_filepath=="")
    error("The fusion product particle of interest ($(getEmittedParticle(reaction))) of the $(reaction) reaction is on an excited state, but no diagnostic line-of-sight was specified (the diagnostic_filepath input variable was left unspecified). This is not allowed. Please correct and re-try.")
end
if !(lowercase(getEmittedParticleEnergyLevel(reaction))=="gs")
    verbose && println("")
    verbose && println("The fusion product particle of interest ($(getEmittedParticle(reaction))) of the $(reaction) reaction is on an excited state!")
    verbose && println("The OWCF will calculate the spectrum of gamma-rays emitted from the de-excitation of the $(getEmittedParticle(reaction)) particle, towards the specified detector.")
    verbose && println("")
end

## ---------------------------------------------------------------------------------------------
# Extract filepath_thermal_distr file extension
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
    println("Thermal distribution file format: ."*fileext_thermal)
    error("Unknown thermal distribution file format. Please re-specify file and re-try.")
end

if fileext_thermal=="cdf" && !(fileext_FI_cdf=="cdf")
    error("filepath_thermal_distr specified as TRANSP .cdf file, but filepath_FI_cdf was wrongly specified. Please specify and re-try.")
end

## ---------------------------------------------------------------------------------------------
# Determine fast-ion and thermal (thermal) species from inputs in start file
thermal_reactant, fast_reactant = getFusionReactants(reaction) # Check the specified fusion reaction, and extract thermal and fast-ion species
@everywhere thermal_reactant = $thermal_reactant # Transfer variable to all external processes
@everywhere fast_reactant = $fast_reactant # Transfer variable to all external processes

# projVel variable. To clarify when weight functions are computed from projected velocities, in code below
projVel = getReactionForm(reaction)==3 ? true : false
@everywhere projVel = $projVel

## ---------------------------------------------------------------------------------------------
# Loading thermal data and/or TRANSP RUN-ID
verbose && println("Loading thermal .jld2 file data and/or TRANSP info... ")

if fileext_thermal=="jld2"
    myfile = jldopen(filepath_thermal_distr,false,false,false,IOStream)
    thermal_temp_array = myfile["thermal_temp"]
    thermal_dens_array = myfile["thermal_dens"]
    ρ_pol_array = myfile["rho_pol"]
    @everywhere thermal_temp_array = $thermal_temp_array # These are sent to external processes here, for efficiency
    @everywhere thermal_dens_array = $thermal_dens_array # These are sent to external processes here, for efficiency
    @everywhere ρ_pol_array = $ρ_pol_array # These are sent to external processes here, for efficiency
    close(myfile)

    @everywhere begin
        verbose && println("Loading the Interpolations package (needed because .jld2 thermal file was specified)... ")
        using Interpolations
    end
end
# Loading TRANSP RUN-ID from thermal .cdf file
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
# Loading magnetic equilibrium
verbose && println("Loading magnetic equilibrium... ")
M, wall, jdotb = nothing, nothing, nothing # Initialize global magnetic equilibrium variables
try
    global M; global wall; global jdotb # Declare global scope
    M, wall = read_geqdsk(filepath_equil,clockwise_phi=false) # Assume counter-clockwise phi-direction
    jdotb = M.sigma # The sign of the dot product between the plasma current and the magnetic field

    # Extract timepoint information from .eqdsk/.geqdsk file
    eqdsk_array = split(filepath_equil,".")
    try
        global timepoint#; global timepoint_source # Declare global scope
        XX = (split(eqdsk_array[end-2],"-"))[end] # Assume format ...-XX.YYYY.eqdsk where XX are the seconds and YYYY are the decimals
        YYYY = eqdsk_array[end-1] # Assume format ...-XX.YYYY.eqdsk where XX are the seconds and YYYY are the decimals
        timepoint = "$(XX*","*YYYY)" # Format XX,YYYY to avoid "." when including in filename of saved output
        #timepoint_source = "EQDSK"
        verbose && println("---> Found timepoint data in magnetic equilibrium file! Loading... ")
    catch
        global timepoint#; global timepoint_source # Declare global scope
        timepoint = "00,0000" # (SOURCE, VALUE). Unknown timepoint for magnetic equilibrium
        #timepoint_source = "UNKNOWN"
    end
catch # Otherwise, assume magnetic equilibrium is a saved .jld2 file
    global M; global wall; global jdotb; global timepoint#; global timepoint_source
    myfile = jldopen(filepath_equil,false,false,false,IOStream)
    M = myfile["S"]
    wall = myfile["wall"]
    close(myfile)
    jdotb = (M.sigma_B0)*(M.sigma_Ip)

    if typeof(timepoint)==String && length(split(timepoint,","))==2
        timepoint = timepoint # (SOURCE, VALUE)
        #timepoint_source = "STARTFILE"
    else
        timepoint = "00,0000" # (SOURCE, VALUE). Unknown timepoint for magnetic equilibrium
        #timepoint_source = "UNKNOWN"
    end
end
psi_axis, psi_bdry = psi_limits(M)
@everywhere psi_axis = $psi_axis # Export to workers
@everywhere psi_bdry = $psi_bdry # Export to workers
@everywhere M = $M # Export to workers
@everywhere wall = $wall # Export to workers

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
# Defining (E,p,R,z) grid vectors
verbose && println("Defining (E,p,R,z) grid vectors... ")
if !isnothing(E_array)
    Emin = minimum(E_array)
    Emax = maximum(E_array)
    nE = length(E_array)
else
    E_array = collect(range(Emin,stop=Emax,length=nE))
end
if !isnothing(p_array)
    p_min = minimum(p_array)
    p_max = maximum(p_array)
    np = length(p_array)
else
    p_array = collect(range(p_min, stop=p_max, length=np))
end
if !isnothing(R_array)
    R_min = minimum(R_array)
    R_max = maximum(R_array)
    nR = length(R_array)
else
    if R_min==:r_mag
        R_min = magnetic_axis(M)[1]
    elseif R_min==:HFSwall
        R_min = minimum(wall.r)
    elseif typeof(R_min) <: Real
    else
        error("R_min wrongly specified. Should be either :r_mag, :HFSwall or an Int64/Float64. Please correct and re-try.")
    end
    if R_max==:r_mag
        R_max = magnetic_axis(M)[1]
    elseif R_max==:LFSwall
        R_max = maximum(wall.r)
    elseif typeof(R_max) <: Real
    else
        error("R_max wrongly specified. Should be either :r_mag, :LFSwall or an Int64/Float64. Please correct and re-try.")
    end
    R_array = collect(range(R_min, stop=R_max, length=nR))
end
if !isnothing(z_array)
    z_min = minimum(z_array)
    z_max = maximum(z_array)
    nz = length(z_array)
else
    if z_min==:z_mag
        z_min = magnetic_axis(M)[2]
    elseif z_min==:floor
        z_min = minimum(wall.z)
    elseif typeof(z_min) <: Real
    else
        error("z_min wrongly specified. Should be either :z_mag, :floor or an Int64/Float64. Please correct and re-try.")
    end
    if z_max==:z_mag
        z_max = magnetic_axis(M)[1]
    elseif z_max==:ceil
        z_max = maximum(wall.z)
    elseif typeof(z_max) <: Real
    else
        error("z_max wrongly specified. Should be either :z_mag, :ceil or an Int64/Float64. Please correct and re-try.")
    end
    z_array = collect(range(z_min, stop=z_max, length=nz))
end
@everywhere E_array = $E_array 
@everywhere p_array = $p_array
@everywhere R_array = $R_array
@everywhere z_array = $z_array

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
println("----------------------------------------calc4DWeights.jl----------------------------------------")
if !(tokamak=="")
    print("Tokamak: "*tokamak*"    ")
else
    print("Tokamak not specified.    ")
end
if !(TRANSP_id=="")
    print("TRANSP ID: "*TRANSP_id*"    ")
else
    print("TRANSP ID not specified.     ")
end
if !(diagnostic_name=="")
    println("Diagnostic name: "*diagnostic_name)
end
println("")
if !(diagnostic_filepath=="")
    println("Diagnostic line-of-sight file: "*diagnostic_filepath)
else
    println("Diagnostic line-of-sight file not specified. Spherical emission will be assumed.")
end
if instrumental_response
    println("Instrumental response file: "*instrumental_response_filepath)
else
    println("Instrumental response file not specified. Instrumental response will not be included.")
end
println("")
if include_FLR_effects
    println("---> Finite Larmor radius (FLR) effects included? Yes!")
else
    println("---> Finite Larmor radius (FLR) effects included? No.")
end
println("")
if !(filepath_thermal_distr=="")
    println("Bulk (thermal) plasma distribution file: "*filepath_thermal_distr)
    println("---> Bulk (thermal) plasma species: $(thermal_reactant)")
elseif (filepath_thermal_distr=="") && !projVel
    println("Bulk (thermal) plasma distribution file not specified.")
    println("Thermal ion ($(getThermalParticle(reaction))) temperature on-axis will be set to $(thermal_temp_axis) keV.")
    (thermal_temp_profile_type==:DEFAULT) && println("---> Default OWCF temperature profile will be used (see OWCF/misc/default_temp_n_dens.png)") 
    (thermal_temp_profile_type==:FLAT) && println("---> Flat temperature profile will be assumed.")
    println("Thermal ion ($(getThermalParticle(reaction))) density on-axis will be set to $(thermal_dens_axis) m^-3.")
    (thermal_dens_profile_type==:DEFAULT) && println("---> Default OWCF density profile will be used (see OWCF/misc/default_temp_n_dens.png)") 
    (thermal_dens_profile_type==:FLAT) && println("---> Flat density profile will be assumed.")
else
    println("4D weight functions will be computed using the projected velocity (u) of the fast ions.")
end
println("")
println("Magnetic equilibrium file: "*filepath_equil)
println("Magnetic field strength on-axis is $(round(norm(Equilibrium.Bfield(M,magnetic_axis(M)...)),digits=2)) T")
println("")
if !projVel
    println("Fusion reaction: "*reaction)
else
    println("4D weight functions will be computed using the projected velocity (u) of fast $(fast_reactant) ions for all relevant (E,p,R,z) points.")
end
println("---> Fast-ion species: "*fast_reactant)
if emittedParticleHasCharge && !projVel
    println("The emitted "*getEmittedParticle(reaction)*" particle of the "*reaction*" reaction has non-zero charge!")
    println("The resulting energy distribution for "*getEmittedParticle(reaction)*" from the plasma as a whole will be computed.")
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
println("4D weight function(s) will be computed for a $(length(E_array))x$(length(p_array))x$(length(R_array))x$(length(z_array)) (E,p,R,z)-space grid with")
println("Energy: [$(minimum(E_array)),$(maximum(E_array))] keV")
println("Pitch: [$(minimum(p_array)),$(maximum(p_array))]")
println("R: [$(minimum(R_array)),$(maximum(R_array))] meters")
println("z: [$(minimum(z_array)),$(maximum(z_array))] meters")
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
    println("Weight functions will be saved in both (E,p,R,z) and (vpara,vperp,R,z) coordinates.")
    println("")
end
println("Results will be saved to: ")
sFLR = !include_FLR_effects ? "" : "FLR_"
t = isnothing(tokamak) ? "" : tokamak
T = isnothing(TRANSP_id) ? "" : TRANSP_id
ti = isnothing(timepoint) ? "" : timepoint
d = isnothing(diagnostic_name) ? "" : diagnostic_name
filepath_o_s_fp = !(filename_o=="") ? filename_o : "EpRzWeights_$(sFLR)"*t*"_"*T*"_at"*ti*"s_"*d*"_"*pretty2scpok(reaction; projVel = projVel)
if iiimax == 1
    println(folderpath_o*filepath_o_s_fp*"_$(length(range(Ed_min,stop=Ed_max,step=Ed_diff))-2)x$(nE)x$(np)x$(nR)x$(nz).jld2")
else
    println(folderpath_o*filepath_o_s_fp*"_1.jld2")
    println("... ")
    println(folderpath_o*filepath_o_s_fp*"_$(iiimax).jld2")
    if iii_average
        println("---> Average of all files will be computed and saved to: ")
        println("---> "*folderpath_o*filepath_o_s_fp*".jld2")
    end
end
if debug
    println("")
    println("!!!!!! DEBUGGING SPECIFIED. ALGORITHM WILL DEBUG. !!!!!!")
    println("")
end
println("")
println("If you would like to change any settings, please edit the start_calc4DW_template.jl file or similar.")
println("")
println("Written by Henrik Järleblad. Last maintained 2025-07-30.")
println("--------------------------------------------------------------------------------------------------")
println("")

## ---------------------------------------------------------------------------------------------
# Python code (essentially calling the forward model black box)
# This is how you write Python code in Julia: py""" [Python code] """
verbose && println("Loading Python modules... ")
@everywhere begin
    py"""
    import numpy as np
    import forward
    import vcone
    import spec
    """
end

## ---------------------------------------------------------------------------------------------
## If no thermal distribution has been specified, we are going to need the OWCF default temp. and dens. profiles
if !isfile(filepath_thermal_distr)
    @everywhere begin
        include("misc/temp_n_dens.jl")
    end
end

## ---------------------------------------------------------------------------------------------
# Loading TRANSP data (if any) and initializing forward model
verbose && println("Loading TRANSP data (if any) and initializing forward model on all distributed workers/processes... ")
@everywhere begin
    py"""
    # The '$' in front of many Python variables means that the variable is defined in Julia, not in Python.
    test_thermal_particle = spec.Particle($thermal_reactant) # Check so that thermal species is available in DRESS code
    projVel = $projVel

    # Load TRANSP simulation data
    if (not $TRANSP_id=="") and (not projVel): # If there is some TRANSP_id specified and we do not want to simply compute projected velocities...
        if ($fileext_thermal).lower()=="cdf": # If there is some TRANSP .cdf output file specified...
            import transp_output
            import transp_dists
            $verbose and print("Loading TRANSP output from TRANSP files... ")
            tr_out = transp_output.TranspOutput($TRANSP_id, step=1, out_file=$filepath_thermal_distr,fbm_files=[$filepath_FI_cdf]) # Load the TRANSP shot file. Assume first step. This is likely to be patched in the future.
            $verbose and print("Setting bulk (thermal) plasma distribution... ")
            thermal_dist = transp_dists.Thermal(tr_out, ion=$thermal_reactant) # Then load the thermal ion distribution from that .cdf file
        else:
            raise ValueError('TRANSP_id was specified, but filepath_thermal_distr was not (this should be impossible). Please correct and re-try.')
    else:
        thermal_dist = "" # Otherwise, just let the thermal_dist variable be the empty string

    $verbose and print("Initializing forward model with diagnostic viewing cone, fusion reaction and bulk (thermal) plasma distribution... ") 
    forwardmodel = forward.Forward($diagnostic_filepath, $reaction, thermal_dist) # Initialize the forward model

    Ed_bin_edges = np.arange($Ed_min,$Ed_max,$Ed_diff) # diagnostic spectrum bin edges (keV or m/s)
    if len(Ed_bin_edges)==1: # Make sure that there are at least one lower and one upper bin edge
        dEd = (($Ed_max)-($Ed_min))/10
        Ed_bin_edges = np.arange($Ed_min,($Ed_max)+dEd,$Ed_diff)
    Ed_vals = 0.5*(Ed_bin_edges[1:] + Ed_bin_edges[:-1]) # bin centers (keV or m/s)
    nEd = len(Ed_vals)
    """
end
nEd = py"nEd"
Ed_array = vec(py"Ed_vals")
@everywhere nEd = $nEd
@everywhere Ed_array = $Ed_array

## ---------------------------------------------------------------------------------------------
# Pre-processing thermal ion density and temperature data

verbose && println("Preparing bulk (thermal) plasma temperature and density data... ")
@everywhere begin
    if !(py"thermal_dist"=="")
        # If !(py"thermal_dist"==""), thermal ion density and temperature data will be loaded from TRANSP file
        (verbose && !(py"thermal_dist"=="")) && println("---> Thermal ion temperature and density data loaded from TRANSP file ($(filepath_thermal_distr))")
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
end

verbose && println("Defining bulk (thermal) plasma temperature and density functions... ")
@everywhere begin
    if !(fileext_thermal=="cdf") # Custom data
        function getThermalSpeciesTemperature(R,z)
            ψ_rz = M(R,z) # Get the interpolated poloidal flux function at the R,z point
            ρ_pol_rz = sqrt(max(0.0, (ψ_rz-psi_axis)/(psi_bdry-psi_axis))) # The formula for the normalized flux coordinate ρ_pol = (ψ-ψ_axis)/(ψ_edge-ψ_axis)
            return thermal_temp_etp(ρ_pol_rz) # Interpolate onto ρ_pol_rz using the data from the .jld2 file
        end
        function getThermalSpeciesDensity(R,z)
            ψ_rz = M(R,z) # Get the interpolated poloidal flux function at the R,z point
            ρ_pol_rz = sqrt(max(0.0, (ψ_rz-psi_axis)/(psi_bdry-psi_axis))) # The formula for the normalized flux coordinate ρ_pol = (ψ-ψ_axis)/(ψ_edge-ψ_axis)
            return thermal_dens_etp(ρ_pol_rz) # Interpolate onto ρ_pol_rz using the data from the .jld2 file
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
# Preparing (E,p,R,z) points
verbose && println("Preparing (E,p,R,z) grid... ")
@everywhere begin
    E_zip = zip(eachindex(E_array),E_array)
    p_zip = zip(eachindex(p_array),p_array)
    R_zip = zip(eachindex(R_array),R_array)
    z_zip = zip(eachindex(z_array),z_array)
    EpRz_zip_itr = Iterators.product(E_zip,p_zip,R_zip,z_zip) # An iterable with all possible combinations of (E,p,R,z) coordinate-index pairs
    EpRz_zip_array = collect(EpRz_zip_itr); npoints = length(EpRz_zip_itr)
end
## ---------------------------------------------------------------------------------------------
# Prepare functions for computing synthetic measurements, based on inputs in start file
verbose && println("Creating synthetic measurements computation function... ")
@everywhere begin
    if include_FLR_effects && py"forwardmodel.bulk_dist_name"=="TRANSP"
        function compute_measurements(E, p, R, z, w)
            N = length(R) # Assume length(E)==length(p)==length(R)==length(z)==length(w)
            E = vcat(E); p = vcat(p); R = vcat(R); z = vcat(z); w = vcat(w) # If not Vectors, make into 1-element Vectors. If already Vectors, keep as is
            B_vecs = zeros(3,N)
            for i in eachindex(R)
                B_vecs[:,i] = collect(reshape(Equilibrium.Bfield(M,R[i],z[i]),3,1))
            end
            py"""
            v, x = forwardmodel.add_gyration($E, $p, $R, $z, $B_vecs, n_gyro=$n_gyro) # No need to specify particle species, already specified via Forward.__init__()
            """

            R_gyro = py"x[0,:]" # Vector (M,) where M >= N
            z_gyro = py"x[2,:]" # Vector (M,)

            weights = inv(n_gyro) .*repeat(w, inner=n_gyro) # Vector (M,), weights rescaled by number of gyro-orbit samples

            B_vecs = reduce(hcat,map(i-> Equilibrium.Bfield(M,R_gyro[i],z_gyro[i]), eachindex(R_gyro))) # Create a (3,M) array with the magnetic field vector at each (R,z) point

            py"""
            spec = forwardmodel.compute_spectrum(Ed_bin_edges, x[0,:], x[2,:], v, $weights, $B_vecs) # Bulk temp and dens already included then forwardmodel was initialized
            """

            return vec(vcat(py"spec"))
        end
    elseif include_FLR_effects && !(py"forwardmodel.bulk_dist_name"=="TRANSP")
        function compute_measurements(E, p, R, z, w)
            N = length(R) # Assume length(E)==length(p)==length(R)==length(z)==length(w)
            E = vcat(E); p = vcat(p); R = vcat(R); z = vcat(z); w = vcat(w) # If not Vectors, make into 1-element Vectors. If already Vectors, keep as is
            B_vecs = zeros(3,N)
            for i in eachindex(R)
                B_vecs[:,i] = collect(reshape(Equilibrium.Bfield(M,R[i],z[i]),3,1))
            end
            py"""
            v, x = forwardmodel.add_gyration($E, $p, $R, $z, $B_vecs, n_gyro=$n_gyro) # No need to specify particle species, already specified via Forward.__init__()
            """

            R_gyro = py"x[0,:]" # Vector (M,) where M >= N
            z_gyro = py"x[2,:]" # Vector (M,)

            weights = inv(n_gyro) .*repeat(w, inner=n_gyro) # Vector (M,), weights rescaled by number of gyro-orbit samples

            B_vecs = reduce(hcat,map(i-> Equilibrium.Bfield(M,R_gyro[i],z_gyro[i]), eachindex(R_gyro))) # Create a (3,M) array with the magnetic field vector at each (R,z) point

            bulk_temp = map(i-> getThermalSpeciesTemperature(R_gyro[i], z_gyro[i]), eachindex(R_gyro)) # Bulk plasma temperature in keV
            bulk_dens = map(i-> getThermalSpeciesDensity(R_gyro[i], z_gyro[i]), eachindex(R_gyro)) # Bulk plasma density in m^-3

            py"""
            spec = forwardmodel.compute_spectrum(Ed_bin_edges, x[0,:], x[2,:], v, $weights, $B_vecs, bulk_temp=$bulk_temp, bulk_dens=$bulk_dens)
            """

            return vec(vcat(py"spec"))
        end
    elseif !include_FLR_effects && py"forwardmodel.bulk_dist_name"=="TRANSP"
        function compute_measurements(E, p, R, z, w)
            N = length(R) # Assume length(E)==length(p)==length(R)==length(z)==length(w)
            E = vcat(E); p = vcat(p); R = vcat(R); z = vcat(z); w = vcat(w) # If not Vectors, make into 1-element Vectors. If already Vectors, keep as is
            B_vecs = zeros(3,N)
            for i in eachindex(R)
                B_vecs[:,i] = collect(reshape(Equilibrium.Bfield(M,R[i],z[i]),3,1))
            end
            py"""
            v, x = forwardmodel.add_gyration($E, $p, $R, $z, $B_vecs, n_gyro=$n_gyro) # No need to specify particle species, already specified via Forward.__init__()
            """

            R = repeat(R, inner=n_gyro) # Ignore FLR effects. 'inner' keyword argument included to match with Python's default repeat() functionality
            z = repeat(z, inner=n_gyro) # Ignore FLR effects
            weights = inv(n_gyro) .*repeat(w, inner=n_gyro) # Vector (M,), weights rescaled by number of gyro-orbit samples
            B_vecs = repeat(B_vecs, inner=(1,n_gyro)) # Array (3,M)

            py"""
            spec = forwardmodel.compute_spectrum(Ed_bin_edges, $R, $z, v, $weights, $B_vecs) # Bulk temp and dens already included then forwardmodel was initialized
            """

            return vec(vcat(py"spec"))
        end
    else # !include_FLR_effects && !(py"forwardmodel.bulk_dist_name"=="TRANSP")
        function compute_measurements(E, p, R, z, w)
            N = length(R) # Assume length(E)==length(p)==length(R)==length(z)==length(w)
            E = vcat(E); p = vcat(p); R = vcat(R); z = vcat(z); w = vcat(w) # If not Vectors, make into 1-element Vectors. If already Vectors, keep as is
            B_vecs = zeros(3,N)
            for i in eachindex(R)
                B_vecs[:,i] = collect(reshape(Equilibrium.Bfield(M,R[i],z[i]),3,1))
            end
            py"""
            v, x = forwardmodel.add_gyration($E, $p, $R, $z, $B_vecs, n_gyro=$n_gyro) # No need to specify particle species, already specified via Forward.__init__()
            """

            R = repeat(R, inner=n_gyro) # Ignore FLR effects. 'inner' keyword argument included to match with Python's default repeat() functionality
            z = repeat(z, inner=n_gyro) # Ignore FLR effects
            weights = inv(n_gyro) .*repeat(w, inner=n_gyro) # Vector (M,), weights rescaled by number of gyro-orbit samples
            B_vecs = repeat(B_vecs, inner=(1,n_gyro)) # Array (3,M)

            bulk_temp = map(i-> getThermalSpeciesTemperature(R[i], z[i]), eachindex(R)) # Bulk plasma temperature in keV
            bulk_dens = map(i-> getThermalSpeciesDensity(R[i], z[i]), eachindex(R)) # Bulk plasma density in m^-3

            py"""
            spec = forwardmodel.compute_spectrum(Ed_bin_edges, $R, $z, v, $weights, $B_vecs, bulk_temp=$bulk_temp, bulk_dens=$bulk_dens)
            """

            return vec(vcat(py"spec"))
        end
    end
end

## ---------------------------------------------------------------------------------------------
# Create a list to keep track of all saved files (needed for averaging)
list_o_saved_filepaths = Vector{String}(undef,iiimax)

## ---------------------------------------------------------------------------------------------
# Calculating the weights
verbose && println("Starting the "*diagnostic_name*" weights calculations... ")
for iii=1:iiimax
    verbose && println("iii: $(iii)")
    if distributed && !debug # If parallel computating is desired (and !debug)...
        if visualizeProgress # if you want the progress to be visualized...
            prog = Progress(npoints) # Create a progress bar, the EpRz_zip_itr is very long
            channel = RemoteChannel(()->Channel{Bool}(npoints), 1) # Utilize a channel
            Wtot = fetch(@sync begin
                @async while take!(channel)
                    ProgressMeter.next!(prog)
                end
                @async begin
                    W = @distributed (+) for i=1:npoints
                        EpRz_zip = EpRz_zip_array[i]
                        E = EpRz_zip[1][2]; p = EpRz_zip[2][2]; R = EpRz_zip[3][2]; z = EpRz_zip[4][2]

                        spec = compute_measurements([E], [p], [R], [z], [1.0]) # Compute the synthetic measurement for the (E, p, R, z) point, weight 1.0
                        rows = append!(collect(1:length(spec)),length(spec)) # To be able to tell the sparse framework about the real size of the weight matrix
                        cols = append!(i .*ones(Int64, length(spec)), npoints) # To be able to tell the sparse framework about the real size of the weight matrix

                        # NOTE! Important to have this append!(spec,0.0) line AFTER the rows and cols assignment lines
                        append!(spec,0.0) # To be able to tell the sparse framework about the real size of the weight matrix
                        Wi = dropzeros(sparse(rows,cols,spec)) # Create a sparse weight matrix, with the current (i) column non-zero. Remove all redundant zeros.

                        put!(channel, true) # Update the progress bar
                        Wi # Declare this weight matrix to be added to the parallel computing reduction (@distributed (+))
                    end
                    put!(channel, false) # Update the progress bar
                    W # Declare the total weight function as ready for fetching
                end
            end)
        else
            Wtot = @distributed (+) for i=1:npoints
                EpRz_zip = EpRz_zip_array[i]
                E = EpRz_zip[1][2]; p = EpRz_zip[2][2]; R = EpRz_zip[3][2]; z = EpRz_zip[4][2]

                spec = compute_measurements([E], [p], [R], [z], [1.0]) # Compute the synthetic measurement for the (E, p, R, z) point, weight 1.0
                rows = append!(collect(1:length(spec)),length(spec)) # Please see similar line earlier in the script
                cols = append!(i .*ones(Int64, length(spec)), npoints) # Please see similar line earlier in the script

                # NOTE! Important to have this append!(spec,0.0) line AFTER the rows and cols assignment lines
                append!(spec,0.0) # Please see similar line earlier in the script
                Wi = dropzeros(sparse(rows,cols,spec)) # Please see similar line earlier in the script
                Wi # Please see similar line earlier in the script
            end
        end
    else # Single-threaded computing... Good luck!
        if debug
            verbose && println("Debugging specified. Only single-threaded debugging is possible. Will debug.")
            # WRITE CODE TO DEBUG QUANTITIES OF INTEREST

        else
            verbose && println("Calculating spectra for (E,p,R,z) point 1 of $(length(EpRz_zip_array))... ")
            EpRz_zip = EpRz_zip_array[1]
            E = EpRz_zip[1][2]; p = EpRz_zip[2][2]; R = EpRz_zip[3][2]; z = EpRz_zip[4][2]

            spec = compute_measurements([E], [p], [R], [z], [1.0]) # Compute the synthetic measurement for the (E, p, R, z) point, weight 1.0
            rows = append!(collect(1:length(spec)),length(spec)) # # Please see similar line earlier in the script
            cols = append!(1 .*ones(Int64, length(spec)), npoints) # # Please see similar line earlier in the script

            # NOTE! Important to have this append!(spec,0.0) line AFTER the rows and cols assignment lines
            append!(spec,0.0)
            Wtot = dropzeros(sparse(rows,cols,spec)) # Pre-allocate a sparse matrix
        end

        for i in eachindex(EpRz_zip_array)[2:end]
            if debug
                verbose && println("Debugging (E,p) point $(i) of $(length(E_iE_p_ip_array))... ")
                # WRITE CODE TO DEBUG QUANTITIES OF INTEREST

            else
                verbose && println("Calculating spectra for (E,p,R,z) point $(i) of $(length(EpRz_zip_itr))... ")
                EpRz_zip = EpRz_zip_array[i]
                E = EpRz_zip[1][2]; p = EpRz_zip[2][2]; R = EpRz_zip[3][2]; z = EpRz_zip[4][2]

                spec = compute_measurements([E], [p], [R], [z], [1.0]) # Compute the synthetic measurement for the (E, p, R, z) point, weight 1.0
                local rows = append!(collect(1:length(spec)),length(spec)) # Please see similar line earlier in the script
                local cols = append!(i .*ones(Int64, length(spec)), npoints) # Please see similar line earlier in the script

                # NOTE! Important to have this append!(spec,0.0) line AFTER the rows and cols assignment lines
                append!(spec,0.0) # Please see similar line earlier in the script
                Wtot += dropzeros(sparse(rows,cols,spec)) # Please see similar line earlier in the script
            end
        end
    end

    verbose && println("Computation of weight matrix $(iii) of $(iiimax): Done!")
    verbose && println("Size of weight matrix: $(size(Wtot))")

    if saveVparaVperpWeights
        # Please note, this algorithm is designed so as to avoid ludicrously large arrays (e.g. (E,p,R,z) weight matrices, essentially huge 5D arrays)
        # There is always a trade-off. Either you have an algorithm that only needs to run through EpRz_zip_itr once, that requires a LOT (!) of RAM.
        # Or, you have an algorithm that needs to run through EpRz_zip_itr many times, and takes longer. But, in return, it does not require as much RAM.
        # Here, to avoid huge 5D arrays (or e.g. equivalent hashmaps), we have chosen the second option.
        verbose && println("Transforming (E,p,R,z) weight matrix to (vpara,vperp,R,z) for row 1 of $(size(Wtot,1))... ")
        nvpara = Int64(round(sqrt(nE*np))) # We know in advance how the Ep2VparaVperp algorithm works by default, so we use sqrt(nE*np)
        nvperp = nvpara
        rows = [1,size(Wtot,1)] # The first and last row of Wtot (and also Wtot_vel)
        cols = [1,nvpara*nvperp*nR*nz] # The first and last column of Wtot_vel (but non Wtot, since (vpara,vperp) will result in another set of weight matrix columns)
        Wtot_vel = dropzeros(sparse(rows,cols,[0.0,0.0])) # Instantiate a sparse matrix for the (vpara,vperp,R,z) weight matrix
        VparVperRz_coords = CartesianIndices((nvpara,nvperp,nR,nz)) # (vpara,vperp,R,z) points corresponding to the columns of the Wtot_vel sparse weight matrix

        Rz_inds = collect(Iterators.product(1:nR,1:nz))
        # Do the first case separately, first weight matrix row and first (R,z) point, to avoid having to re-assign vpara_array, vperp_array all the time
        verbose && println("---> Transforming w(E,p) for (R,z) point 1 of $(length(Rz_inds)) from (E,p) to (vpara,vperp)... ")
        ### DEBUG PRINT ### println(""); println(Wtot); println("")
        iR = Rz_inds[1][1]
        iz = Rz_inds[1][2]
        W_Ep = zeros(nE,np)
        for (i,EpRz_zip) in enumerate(EpRz_zip_itr) # Go through all (E,p,R,z) points
            if (EpRz_zip[3][1]==iR && EpRz_zip[4][1]==iz) # If the (R,z) indices match the indices of interest
                W_Ep[EpRz_zip[1][1], EpRz_zip[2][1]] = Wtot[1, i] # Build w(E,p) from the first row of the weight matrix Wtot
            end
        end
        global vpara_array; global vperp_array
        ### DEBUG ### println(W_Ep); println("")
        W_vel, vpara_array, vperp_array = Ep2VparaVperp(E_array, p_array, W_Ep; returnAbscissas=true) # Transform w(E,p) to w(vpara,vperp)
        ### DEBUG ### println(W_vel); println("")
        ### Now we have to take the elements of W_vel (w(vpara,vperp)) and put them in the right place in the sparse weight matrix Wtot_vel
        cinds_wvel = CartesianIndices(W_vel)
        ROW = append!(1 .*ones(Int64,length(cinds_wvel)),size(Wtot,1)) # Here, we are only doing the first row of the sparse weight matrix Wtot_vel. Append dummy element, to create correct size
        COL = append!(zeros(Int64,length(cinds_wvel)),nvpara*nvperp*nR*nz) # But many different columns. Append dummy element, to create correct size
        VAL = append!(zeros(length(cinds_wvel)),0.0) # And values. Append dummy element, to create correct size
        for (iv,vel_inds) in enumerate(cinds_wvel) # Go through all (vpara,vperp) points
            ivpar = vel_inds[1]
            ivper = vel_inds[2]
            # Here, we have to find the matching column in Wtot_vel using ivpar, ivper, ir and iz.
            COL[iv] = findfirst(c -> (c[1]==ivpar && c[2]==ivper && c[3]==iR && c[4]==iz), VparVperRz_coords[:]) # Should be only one matching element
            VAL[iv] = W_vel[ivpar,ivper] # The weight value
        end
        Wtot_vel += dropzeros(sparse(ROW,COL,VAL)) # Build upon the (vpara,vperp,R,z) sparse matrix Wtot_vel, instantiated above

        # NOW, REPEAT FOR ALL (R,z) POINTS FOR THE FIRST WEIGHT MATRIX ROW
        global prog_proc = []
        for iRz in 2:length(Rz_inds)
            if !(floor(100*iRz/length(Rz_inds)) in prog_proc)
                append!(prog_proc,floor(100*iRz/length(Rz_inds)))
                verbose && println("---> Transforming w(E,p) to w(vpara,vperp) for the rest of the (R,z) points for row 1 of $(size(Wtot,1)) ($(prog_proc[end]) %)... ")
            end
            iR = Rz_inds[iRz][1]
            iz = Rz_inds[iRz][2]
            W_Ep = zeros(nE,np)
            for (i,EpRz_zip) in enumerate(EpRz_zip_itr) # Go through all (E,p,R,z) points
                if (EpRz_zip[3][1]==iR && EpRz_zip[4][1]==iz) # If the (R,z) indices match the indices of interest
                    W_Ep[EpRz_zip[1][1], EpRz_zip[2][1]] = Wtot[1, i] # Build w(E,p) from the first row of the weight matrix Wtot
                end
            end
            W_vel = Ep2VparaVperp(E_array, p_array, W_Ep; returnAbscissas=false) # Transform w(E,p) to w(vpara,vperp)
            # Now we have to take the elements of W_vel (w(vpara,vperp)) and put them in the right place in the sparse weight matrix Wtot_vel
            cinds_wvel = CartesianIndices(W_vel)
            ROW = append!(1 .*ones(Int64,length(cinds_wvel)),size(Wtot,1)) # Here, we are only doing the first row of the sparse weight matrix Wtot_vel. Append dummy element, to create correct size
            COL = append!(zeros(Int64,length(cinds_wvel)),nvpara*nvperp*nR*nz) # But many different columns. Append dummy element, to create correct size
            VAL = append!(zeros(length(cinds_wvel)),0.0) # And values. Append dummy element, to create correct size
            for (iv,vel_inds) in enumerate(cinds_wvel) # Go through all (vpara,vperp) points
                ivpar = vel_inds[1]
                ivper = vel_inds[2]
                # Here, we have to find the matching column in Wtot_vel using ivpar, ivper, ir and iz.
                COL[iv] = findfirst(c -> (c[1]==ivpar && c[2]==ivper && c[3]==iR && c[4]==iz), VparVperRz_coords[:]) # Should be only one matching element
                VAL[iv] = W_vel[ivpar,ivper] # The weight value
            end
            Wtot_vel += dropzeros(sparse(ROW,COL,VAL)) # Build upon the (vpara,vperp,R,z) sparse matrix Wtot_vel, instantiated way above
        end

        # NOW, REPEAT FOR ALL WEIGHT MATRIX ROWS 
        verbose && println("Transforming (E,p,R,z) weight matrix to (vpara,vperp,R,z) for all rows >1.. ")
        if distributed && !debug # If parallel computating is desired (and !debug)...
            if visualizeProgress
                prog = Progress(size(Wtot,1)-1) # Create a progress bar
                channel = RemoteChannel(()->Channel{Bool}(size(Wtot,1)-1), 1) # Utilize a channel
                Wtot_vel_exceptFirstRow = fetch(@sync begin
                    @async while take!(channel)
                        ProgressMeter.next!(prog)
                    end
                    @async begin
                        Wtot_vel_exceptFirstRow_preFetch = @distributed (+) for iRow in collect(2:size(Wtot,1))
                            rows = [1,size(Wtot,1)] # The first and last row of Wtot (and also Wtot_vel)
                            cols = [1,nvpara*nvperp*nR*nz] # The first and last column of Wtot_vel (but non Wtot, since (vpara,vperp) will result in another set of weight matrix columns)
                            Wtot_vel_row = dropzeros(sparse(rows,cols,[0.0,0.0])) # Instantiate a sparse matrix for the (vpara,vperp,R,z) weight matrix (but only the iRow:th row)
                            for iRz in 1:length(Rz_inds)
                                iR = Rz_inds[iRz][1]
                                iz = Rz_inds[iRz][2]
                                W_Ep = zeros(nE,np)
                                for (i,EpRz_zip) in enumerate(EpRz_zip_itr) # Go through all (E,p,R,z) points
                                    if (EpRz_zip[3][1]==iR && EpRz_zip[4][1]==iz) # If the (R,z) indices match the indices of interest
                                        W_Ep[EpRz_zip[1][1], EpRz_zip[2][1]] = Wtot[iRow, i] # Build w(E,p) from the iRow:th row of the weight matrix Wtot
                                    end
                                end
                                W_vel = Ep2VparaVperp(E_array, p_array, W_Ep; returnAbscissas=false) # Transform w(E,p) to w(vpara,vperp)
                                # Now we have to take the elements of W_vel (w(vpara,vperp)) and put them in the right place in the sparse weight matrix Wtot_vel
                                cinds_wvel = CartesianIndices(W_vel)
                                ROW = append!(iRow .*ones(Int64,length(cinds_wvel)),size(Wtot,1)) # Here, we are only doing the iRow:th row of the sparse weight matrix Wtot_vel. Append dummy element, to create correct size
                                COL = append!(zeros(Int64,length(cinds_wvel)),nvpara*nvperp*nR*nz) # But many different columns. Append dummy element, to create correct size
                                VAL = append!(zeros(length(cinds_wvel)),0.0) # And values. Append dummy element, to create correct size
                                for (iv,vel_inds) in enumerate(cinds_wvel) # Go through all (vpara,vperp) points
                                    ivpar = vel_inds[1]
                                    ivper = vel_inds[2]
                                    # Here, we have to find the matching column in Wtot_vel using ivpar, ivper, ir and iz.
                                    COL[iv] = findfirst(c -> (c[1]==ivpar && c[2]==ivper && c[3]==iR && c[4]==iz), VparVperRz_coords[:]) # Should be only one matching element
                                    VAL[iv] = W_vel[ivpar,ivper] # The weight value
                                end
                                Wtot_vel_row += dropzeros(sparse(ROW,COL,VAL)) # Build upon the (vpara,vperp,R,z) sparse matrix Wtot_vel_row, instantiated above
                            end
                        
                            put!(channel, true) # Update the progress bar
                            Wtot_vel_row # Declare this row of the weight matrix to be added to the parallel computing reduction (@distributed (+))
                        end
                        put!(channel, false) # Update the progress bar
                        Wtot_vel_exceptFirstRow_preFetch # Declare the total (except first row) weight matrix as ready for fetching
                    end
                end)
            else
                Wtot_vel_exceptFirstRow = @distributed (+) for iRow in collect(2:size(Wtot,1))
                    rows = [1,size(Wtot,1)] # The first and last row of Wtot (and also Wtot_vel)
                    cols = [1,nvpara*nvperp*nR*nz] # The first and last column of Wtot_vel (but non Wtot, since (vpara,vperp) will result in another set of weight matrix columns)
                    Wtot_vel_row = dropzeros(sparse(rows,cols,[0.0,0.0])) # Instantiate a sparse matrix for the (vpara,vperp,R,z) weight matrix (but only the iRow:th row)
                    for iRz in 1:length(Rz_inds)
                        iR = Rz_inds[iRz][1]
                        iz = Rz_inds[iRz][2]
                        W_Ep = zeros(nE,np)
                        for (i,EpRz_zip) in enumerate(EpRz_zip_itr) # Go through all (E,p,R,z) points
                            if (EpRz_zip[3][1]==iR && EpRz_zip[4][1]==iz) # If the (R,z) indices match the indices of interest
                                W_Ep[EpRz_zip[1][1], EpRz_zip[2][1]] = Wtot[iRow, i] # Build w(E,p) from the iRow:th row of the weight matrix Wtot
                            end
                        end
                        W_vel = Ep2VparaVperp(E_array, p_array, W_Ep; returnAbscissas=false) # Transform w(E,p) to w(vpara,vperp)
                        # Now we have to take the elements of W_vel (w(vpara,vperp)) and put them in the right place in the sparse weight matrix Wtot_vel
                        cinds_wvel = CartesianIndices(W_vel)
                        ROW = append!(iRow .*ones(Int64,length(cinds_wvel)),size(Wtot,1)) # Here, we are only doing the iRow:th row of the sparse weight matrix Wtot_vel. Append dummy element, to create correct size
                        COL = append!(zeros(Int64,length(cinds_wvel)),nvpara*nvperp*nR*nz) # But many different columns. Append dummy element, to create correct size
                        VAL = append!(zeros(length(cinds_wvel)),0.0) # And values. Append dummy element, to create correct size
                        for (iv,vel_inds) in enumerate(cinds_wvel) # Go through all (vpara,vperp) points
                            ivpar = vel_inds[1]
                            ivper = vel_inds[2]
                            # Here, we have to find the matching column in Wtot_vel using ivpar, ivper, ir and iz.
                            COL[iv] = findfirst(c -> (c[1]==ivpar && c[2]==ivper && c[3]==iR && c[4]==iz), VparVperRz_coords[:]) # Should be only one matching element
                            VAL[iv] = W_vel[ivpar,ivper] # The weight value
                        end
                        Wtot_vel_row += dropzeros(sparse(ROW,COL,VAL)) # Build upon the (vpara,vperp,R,z) sparse matrix Wtot_vel_row, instantiated above
                    end
                    Wtot_vel_row # Declare this row of the weight matrix to be added to the parallel computing reduction (@distributed (+))
                end
            end
        else # Single-threaded computing of the rest of the weight matrix transformation from (E,p,R,z) to (vpara,vperp,R,z)
            rows = [1,size(Wtot,1)] # The first and last row of Wtot (and also Wtot_vel)
            cols = [1,nvpara*nvperp*nR*nz] # The first and last column of Wtot_vel (but non Wtot, since (vpara,vperp) will result in another set of weight matrix columns)
            Wtot_vel_exceptFirstRow = dropzeros(sparse(rows,cols,[0.0,0.0]))
            for iRow in collect(2:size(Wtot,1))
                verbose && println("Transforming row $(iRow) of $(size(Wtot,1)) from (E,p) to (vpara,vperp)... ")
                rows = [1,size(Wtot,1)] # The first and last row of Wtot (and also Wtot_vel)
                cols = [1,nvpara*nvperp*nR*nz] # The first and last column of Wtot_vel (but non Wtot, since (vpara,vperp) will result in another set of weight matrix columns)
                Wtot_vel_row = dropzeros(sparse(rows,cols,[0.0,0.0])) # Instantiate a sparse matrix for the (vpara,vperp,R,z) weight matrix (but only the iRow:th row)
                global prog_proc = []
                for iRz in 1:length(Rz_inds)
                    if !(floor(100*iRz/length(Rz_inds)) in prog_proc)
                        append!(prog_proc,floor(100*iRz/length(Rz_inds)))
                        verbose && println("---> Transforming w(E,p) to w(vpara,vperp) for all (R,z) points for row $(iRow) of $(size(Wtot,1)) ($(prog_proc[end]) %)... ")
                    end
                    iR = Rz_inds[iRz][1]
                    iz = Rz_inds[iRz][2]
                    W_Ep = zeros(nE,np)
                    for (i,EpRz_zip) in enumerate(EpRz_zip_itr) # Go through all (E,p,R,z) points
                        if (EpRz_zip[3][1]==iR && EpRz_zip[4][1]==iz) # If the (R,z) indices match the indices of interest
                            W_Ep[EpRz_zip[1][1], EpRz_zip[2][1]] = Wtot[iRow, i] # Build w(E,p) from the iRow:th row of the weight matrix Wtot
                        end
                    end
                    W_vel = Ep2VparaVperp(E_array, p_array, W_Ep; returnAbscissas=false) # Transform w(E,p) to w(vpara,vperp)
                    # Now we have to take the elements of W_vel (w(vpara,vperp)) and put them in the right place in the sparse weight matrix Wtot_vel
                    cinds_wvel = CartesianIndices(W_vel)
                    ROW = append!(iRow .*ones(Int64,length(cinds_wvel)),size(Wtot,1)) # Here, we are only doing the iRow:th row of the sparse weight matrix Wtot_vel. Append dummy element, to create correct size
                    COL = append!(zeros(Int64,length(cinds_wvel)),nvpara*nvperp*nR*nz) # But many different columns. Append dummy element, to create correct size
                    VAL = append!(zeros(length(cinds_wvel)),0.0) # And values. Append dummy element, to create correct size
                    for (iv,vel_inds) in enumerate(cinds_wvel) # Go through all (vpara,vperp) points
                        ivpar = vel_inds[1]
                        ivper = vel_inds[2]
                        # Here, we have to find the matching column in Wtot_vel using ivpar, ivper, ir and iz.
                        COL[iv] = findfirst(c -> (c[1]==ivpar && c[2]==ivper && c[3]==iR && c[4]==iz), VparVperRz_coords[:]) # Should be only one matching element
                        VAL[iv] = W_vel[ivpar,ivper] # The weight value
                    end
                    Wtot_vel_row += dropzeros(sparse(ROW,COL,VAL)) # Build upon the (vpara,vperp,R,z) sparse matrix Wtot_vel_row, instantiated above
                end
                Wtot_vel_exceptFirstRow += Wtot_vel_row
            end
        end
        Wtot_vel += Wtot_vel_exceptFirstRow
        VparVperRz_coords = map(c -> ((c[1],vpara_array[c[1]]),(c[2],vperp_array[c[2]]),(c[3],R_array[c[3]]),(c[4],z_array[c[4]])),VparVperRz_coords) # Include values with indices
    end

    if instrumental_response
        verbose && println("Applying diagnostic response to weight functions... ")
        Wtot_raw = deepcopy(Wtot) # To be able to save the weight functions without instrumental response as well
        Wtot = zeros(length(instrumental_response_output),size(Wtot,2))
        for iCol=1:size(Wtot,2)
            Wtot[:,iCol] = apply_instrumental_response(Wtot_raw[:,iCol], Ed_array, instrumental_response_input, instrumental_response_output, instrumental_response_matrix)
        end
        Wtot = dropzeros(sparse(Wtot)) # Remake into a sparse matrix representation
        if saveVparaVperpWeights
            Wtot_vel_raw = deepcopy(Wtot_vel)
            Wtot_vel = zeros(length(instrumental_response_output),size(Wtot_vel,2))
            for iCol=1:size(Wtot_vel,2)
                Wtot_vel[:,iCol] = apply_instrumental_response(Wtot_vel_raw[:,iCol], Ed_array, instrumental_response_input, instrumental_response_output, instrumental_response_matrix)
            end
            Wtot_vel = dropzeros(sparse(Wtot_vel))
        end
    end

    if !debug
        # Set the output file name 'filepath_out'
        # First, if the user specified a custom output file name (filename_o), use that instead of the default OWCF output file name
        filepath_o_s = !(filename_o=="") ? filename_o : "EpRzWeights_$(sFLR)"*tokamak*"_"*TRANSP_id*"_at"*timepoint*"s_"*diagnostic_name*"_"*pretty2scpok(reaction; projVel = projVel)
        if iiimax==1 # If you intend to calculate only one weight matrix
            global filepath_output_orig = folderpath_o*filepath_o_s*"_$(length(range(Ed_min,stop=Ed_max,step=Ed_diff))-2)x$(nE)x$(np)x$(nR)x$(nz)"
        else # If you intend to calculate several (identical) weight matrices
            global filepath_output_orig = folderpath_o*filepath_o_s*"_$(iii)"
        end
        global filepath_output = deepcopy(filepath_output_orig)
        global count = 1
        while isfile(filepath_output*".jld2") # To take care of not overwriting files. Add _(1), _(2) etc
            global filepath_output; global count
            filepath_output = filepath_output_orig*"_($(Int64(count)))"
            count += 1 # global scope, to surpress warnings
        end

        if plot_results
            # DEVELOPER COMMENT: Hopefully, sometime in the future, these plot routines are going to be re-written as one, single function in OWCF/extra/gui.jl. But for now, let the ugly code be...
            # DEVELOPER COMMENT: Hopefully, sometime in the future, these plot routines are going to be re-written as one, single function in OWCF/extra/gui.jl. But for now, let the ugly code be...
            # DEVELOPER COMMENT: Hopefully, sometime in the future, these plot routines are going to be re-written as one, single function in OWCF/extra/gui.jl. But for now, let the ugly code be...
            plot_font = "Computer Modern"
            Plots.default(fontfamily=plot_font)
            verbose && println("Plotting weight function data... ")

            if instrumental_response
                W_raw_plt = Wtot_raw # Bad variable names...
                W_plt = Wtot # Bad variable names...
                if saveVparaVperpWeights
                    W_vel_raw_plt = Wtot_vel_raw # Bad variable names...
                    W_vel_plt = Wtot_vel # Bad variable names...
                end
            else
                W_raw_plt = Wtot # Bad variable names...
                if saveVparaVperpWeights
                    W_vel_raw_plt = Wtot_vel # Bad variable names...
                end
            end

            # Without instrumental response (raw), (E,p)
            N_bins = length(Ed_array)
            if N_bins==1
                plt_raw_inds = [1,1,1]
            elseif N_bins==2
                plt_raw_inds = [1,2,2]
            elseif N_bins==3
                plt_raw_inds = [1,2,3]
            else
                # Algorithm to find suitable indices of Ed to visualize
                W_gross = dropdims(sum(W_raw_plt, dims=2), dims=2)
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
                    iEd_low = max(1, iEd_mid - abs(iEd_diff))
                    iEd_hig = min(N_bins, iEd_mid + abs(iEd_diff))
                end
                plt_raw_inds = [iEd_low, iEd_mid, iEd_hig]
            end

            Ed_low = @sprintf "%.2E" Ed_array[plt_raw_inds[1]]
            Ed_mid = @sprintf "%.2E" Ed_array[plt_raw_inds[2]]
            Ed_hi = @sprintf "%.2E" Ed_array[plt_raw_inds[3]]
            W_raw_plt_low = W_raw_plt[plt_raw_inds[1],:] # w_low(E,p,R,z) but vectorized
            W_raw_plt_mid = W_raw_plt[plt_raw_inds[2],:] # w_mid(E,p,R,z) but vectorized
            W_raw_plt_hi = W_raw_plt[plt_raw_inds[3],:] # w_hi(E,p,R,z) but vectorized

            W_raw_plt_low_Ep = zeros(length(E_array),length(p_array)) # w_low(E,p) = inv(dRdz) ∫ w_low(E,p,R,z) dRdz 
            W_raw_plt_mid_Ep = zeros(length(E_array),length(p_array)) # w_mid(E,p) = inv(dRdz) ∫ w_mid(E,p,R,z) dRdz 
            W_raw_plt_hi_Ep = zeros(length(E_array),length(p_array)) # w_hi(E,p) = inv(dRdz) ∫ w_hi(E,p,R,z) dRdz 
            W_raw_plt_low_Rz = zeros(length(R_array),length(z_array)) # w_low(R,z) = inv(dEdp) ∫ w_low(E,p,R,z) dEdp 
            W_raw_plt_mid_Rz = zeros(length(R_array),length(z_array)) # w_mid(R,z) = inv(dEdp) ∫ w_mid(E,p,R,z) dEdp 
            W_raw_plt_hi_Rz = zeros(length(R_array),length(z_array)) # w_hi(R,z) = inv(dEdp) ∫ w_hi(E,p,R,z) dEdp 

            for (i,elem) in enumerate(EpRz_zip_array)
                iE = elem[1][1]; ip = elem[2][1]; iR = elem[3][1]; iz = elem[4][1]

                W_raw_plt_low_Ep[iE,ip] += W_raw_plt_low[i]
                W_raw_plt_mid_Ep[iE,ip] += W_raw_plt_mid[i]
                W_raw_plt_hi_Ep[iE,ip] += W_raw_plt_hi[i]
                W_raw_plt_low_Rz[iR,iz] += W_raw_plt_low[i]
                W_raw_plt_mid_Rz[iR,iz] += W_raw_plt_mid[i]
                W_raw_plt_hi_Rz[iR,iz] += W_raw_plt_hi[i]
            end

            W_raw_plt_low_Ep_max = @sprintf "%.2E" maximum(W_raw_plt_low_Ep)
            W_raw_plt_mid_Ep_max = @sprintf "%.2E" maximum(W_raw_plt_mid_Ep)
            W_raw_plt_hi_Ep_max = @sprintf "%.2E" maximum(W_raw_plt_hi_Ep)
            W_raw_plt_low_Rz_max = @sprintf "%.2E" maximum(W_raw_plt_low_Rz)
            W_raw_plt_mid_Rz_max = @sprintf "%.2E" maximum(W_raw_plt_mid_Rz)
            W_raw_plt_hi_Rz_max = @sprintf "%.2E" maximum(W_raw_plt_hi_Rz)

            plt_Ep_raw_low = Plots.heatmap(E_array, p_array, transpose(W_raw_plt_low_Ep),title="∫dRdz w($(Ed_low),E,p,R,z), max(w): $(W_raw_plt_low_Ep_max) keV^-1", colorbar=false, xlims=extrema(E_array), ylims=extrema(p_array))
            plt_Ep_raw_mid = Plots.heatmap(E_array, p_array, transpose(W_raw_plt_mid_Ep),title="∫dRdz w($(Ed_mid),E,p,R,z), max(w): $(W_raw_plt_mid_Ep_max) keV^-1",colorbar=false, xlims=extrema(E_array), ylims=extrema(p_array))
            plt_Ep_raw_hi = Plots.heatmap(E_array, p_array, transpose(W_raw_plt_hi_Ep),title="∫dRdz w($(Ed_hi),E,p,R,z), max(w): $(W_raw_plt_hi_Ep_max) keV^-1",colorbar=false, xlims=extrema(E_array), ylims=extrema(p_array))
            plt_Rz_raw_low = Plots.heatmap(R_array, z_array, transpose(W_raw_plt_low_Rz),title="∫dEdp w($(Ed_low),E,p,R,z), max(w): $(W_raw_plt_low_Rz_max) keV^-1",colorbar=false)
            plt_Rz_raw_mid = Plots.heatmap(R_array, z_array, transpose(W_raw_plt_mid_Rz),title="∫dEdp w($(Ed_mid),E,p,R,z), max(w): $(W_raw_plt_mid_Rz_max) keV^-1",colorbar=false)
            plt_Rz_raw_hi = Plots.heatmap(R_array, z_array, transpose(W_raw_plt_hi_Rz),title="∫dEdp w($(Ed_hi),E,p,R,z), max(w): $(W_raw_plt_hi_Rz_max) keV^-1",colorbar=false)
            
            # Add labels etc to (E,p) plots
            plt_Ep_raw_low = Plots.plot!(plt_Ep_raw_low, xlabel="Energy [keV]", ylabel="Pitch [-]")
            plt_Ep_raw_mid = Plots.plot!(plt_Ep_raw_mid, xlabel="Energy [keV]", ylabel="Pitch [-]")
            plt_Ep_raw_hi = Plots.plot!(plt_Ep_raw_hi, xlabel="Energy [keV]", ylabel="Pitch [-]")

            # Add tokamak wall, labels etc to (R,z) plots
            plt_Rz_raw_low = Plots.plot!(plt_Rz_raw_low, wall.r, wall.z, linewidth=2.0, label="", xlabel="R [m]", ylabel="z [m]", aspect_ratio=:equal)
            plt_Rz_raw_mid = Plots.plot!(plt_Rz_raw_mid, wall.r, wall.z, linewidth=2.0, label="", xlabel="R [m]", ylabel="z [m]", aspect_ratio=:equal)
            plt_Rz_raw_hi = Plots.plot!(plt_Rz_raw_hi, wall.r, wall.z, linewidth=2.0, label="", xlabel="R [m]", ylabel="z [m]", aspect_ratio=:equal)

            plt_EpRz_raw = Plots.plot(plt_Ep_raw_low, plt_Ep_raw_mid, plt_Ep_raw_hi,
                                      plt_Rz_raw_low, plt_Rz_raw_mid, plt_Rz_raw_hi,
                                      layout=(2,3), size=(1200,800), dpi=200, 
                                fillcolor=cgrad([:white, :yellow, :orange, :red, :black]), titlefontsize=10, colorbar=false, bottom_margin=6Plots.mm,
                                left_margin=6Plots.mm, right_margin=5Plots.mm)
            png(plt_EpRz_raw,filepath_output*"_EpRz_raw")

            if saveVparaVperpWeights
                # Without instrumental response (raw), (vpara,vperp)
                W_vel_raw_plt_low = W_vel_raw_plt[plt_raw_inds[1],:] # w_low(vpara,vperp,R,z) but vectorized
                W_vel_raw_plt_mid = W_vel_raw_plt[plt_raw_inds[2],:] # w_mid(vpara,vperp,R,z) but vectorized
                W_vel_raw_plt_hi  = W_vel_raw_plt[plt_raw_inds[3],:] # w_hi(vpara,vperp,R,z)  but vectorized

                W_vel_raw_plt_low_VpaVpe = zeros(length(vpara_array),length(vperp_array)) # w_low(vpara,vperp) = inv(dRdz) ∫ w_low(vpara,vperp,R,z) dRdz 
                W_vel_raw_plt_mid_VpaVpe = zeros(length(vpara_array),length(vperp_array)) # w_mid(vpara,vperp) = inv(dRdz) ∫ w_mid(vpara,vperp,R,z) dRdz 
                W_vel_raw_plt_hi_VpaVpe = zeros(length(vpara_array),length(vperp_array)) # w_hi(vpara,vperp) = inv(dRdz) ∫ w_hi(vpara,vperp,R,z) dRdz 
                W_vel_raw_plt_low_Rz = zeros(length(R_array),length(z_array)) # w_low(R,z) = inv(dVpa*dVpe) ∫ w_low(vpara,vperp,R,z) dVpa*dVpe 
                W_vel_raw_plt_mid_Rz = zeros(length(R_array),length(z_array)) # w_mid(R,z) = inv(dVpa*dVpe) ∫ w_mid(vpara,vperp,R,z) dVpa*dVpe 
                W_vel_raw_plt_hi_Rz = zeros(length(R_array),length(z_array)) # w_hi(R,z) = inv(dVpa*dVpe) ∫ w_hi(vpara,vperp,R,z) dVpa*dVpe

                for (i,elem) in enumerate(VparVperRz_coords)
                    iVpa = elem[1][1]; iVpe = elem[2][1]; iR = elem[3][1]; iz = elem[4][1]

                    W_vel_raw_plt_low_VpaVpe[iVpa,iVpe] += W_vel_raw_plt_low[i]
                    W_vel_raw_plt_mid_VpaVpe[iVpa,iVpe] += W_vel_raw_plt_mid[i]
                    W_vel_raw_plt_hi_VpaVpe[iVpa,iVpe] += W_vel_raw_plt_hi[i]
                    W_vel_raw_plt_low_Rz[iR,iz] += W_vel_raw_plt_low[i]
                    W_vel_raw_plt_mid_Rz[iR,iz] += W_vel_raw_plt_mid[i]
                    W_vel_raw_plt_hi_Rz[iR,iz] += W_vel_raw_plt_hi[i]
                end

                W_vel_raw_plt_low_VpaVpe_max = @sprintf "%.2E" maximum(W_vel_raw_plt_low_VpaVpe)
                W_vel_raw_plt_mid_VpaVpe_max = @sprintf "%.2E" maximum(W_vel_raw_plt_mid_VpaVpe)
                W_vel_raw_plt_hi_VpaVpe_max = @sprintf "%.2E" maximum(W_vel_raw_plt_hi_VpaVpe)
                W_vel_raw_plt_low_Rz_max = @sprintf "%.2E" maximum(W_vel_raw_plt_low_Rz)
                W_vel_raw_plt_mid_Rz_max = @sprintf "%.2E" maximum(W_vel_raw_plt_mid_Rz)
                W_vel_raw_plt_hi_Rz_max = @sprintf "%.2E" maximum(W_vel_raw_plt_hi_Rz)

                plt_vel_VpaVpe_raw_low = Plots.heatmap(vpara_array, vperp_array, transpose(W_vel_raw_plt_low_VpaVpe),title="∫dRdz w($(Ed_low),vpa,vpe,R,z), max(w): $(W_vel_raw_plt_low_VpaVpe_max) keV^-1", colorbar=false, xlims=extrema(vpara_array), ylims=extrema(vperp_array))
                plt_vel_VpaVpe_raw_mid = Plots.heatmap(vpara_array, vperp_array, transpose(W_vel_raw_plt_mid_VpaVpe),title="∫dRdz w($(Ed_mid),vpa,vpe,R,z), max(w): $(W_vel_raw_plt_mid_VpaVpe_max) keV^-1",colorbar=false, xlims=extrema(vpara_array), ylims=extrema(vperp_array))
                plt_vel_VpaVpe_raw_hi = Plots.heatmap(vpara_array, vperp_array, transpose(W_vel_raw_plt_hi_VpaVpe),title="∫dRdz w($(Ed_hi),vpa,vpe,R,z), max(w): $(W_vel_raw_plt_hi_VpaVpe_max) keV^-1",colorbar=false, xlims=extrema(vpara_array), ylims=extrema(vperp_array))
                plt_vel_Rz_raw_low = Plots.heatmap(R_array, z_array, transpose(W_vel_raw_plt_low_Rz),title="∫dVpa*dVpe w($(Ed_low),vpa,vpe,R,z), max(w): $(W_vel_raw_plt_low_Rz_max) keV^-1",colorbar=false)
                plt_vel_Rz_raw_mid = Plots.heatmap(R_array, z_array, transpose(W_vel_raw_plt_mid_Rz),title="∫dVpa*dVpe w($(Ed_mid),vpa,vpe,R,z), max(w): $(W_vel_raw_plt_mid_Rz_max) keV^-1",colorbar=false)
                plt_vel_Rz_raw_hi = Plots.heatmap(R_array, z_array, transpose(W_vel_raw_plt_hi_Rz),title="∫dVpa*dVpe w($(Ed_hi),vpa,vpe,R,z), max(w): $(W_vel_raw_plt_hi_Rz_max) keV^-1",colorbar=false)

                # Add labels etc to (vpa,vpe) plots
                plt_vel_VpaVpe_raw_low = Plots.plot!(plt_vel_VpaVpe_raw_low, xlabel="Vpa [m/s]", ylabel="Vpe [m/s]")
                plt_vel_VpaVpe_raw_mid = Plots.plot!(plt_vel_VpaVpe_raw_mid, xlabel="Vpa [m/s]", ylabel="Vpe [m/s]")
                plt_vel_VpaVpe_raw_hi = Plots.plot!(plt_vel_VpaVpe_raw_hi, xlabel="Vpa [m/s]", ylabel="Vpe [m/s]")

                # Add tokamak wall, labels etc to (R,z) plots
                plt_vel_Rz_raw_low = Plots.plot!(plt_vel_Rz_raw_low, wall.r, wall.z, linewidth=2.0, label="", xlabel="R [m]", ylabel="z [m]", aspect_ratio=:equal)
                plt_vel_Rz_raw_mid = Plots.plot!(plt_vel_Rz_raw_mid, wall.r, wall.z, linewidth=2.0, label="", xlabel="R [m]", ylabel="z [m]", aspect_ratio=:equal)
                plt_vel_Rz_raw_hi = Plots.plot!(plt_vel_Rz_raw_hi, wall.r, wall.z, linewidth=2.0, label="", xlabel="R [m]", ylabel="z [m]", aspect_ratio=:equal)

                plt_VpaVpeRz_raw = Plots.plot(plt_vel_VpaVpe_raw_low, plt_vel_VpaVpe_raw_mid, plt_vel_VpaVpe_raw_hi,
                                      plt_vel_Rz_raw_low, plt_vel_Rz_raw_mid, plt_vel_Rz_raw_hi,
                                      layout=(2,3), size=(1200,800), dpi=200, 
                                fillcolor=cgrad([:white, :yellow, :orange, :red, :black]), titlefontsize=10, colorbar=false, bottom_margin=6Plots.mm,
                                left_margin=6Plots.mm, right_margin=5Plots.mm)
                png(plt_VpaVpeRz_raw,filepath_output*"_VpaVpeRz_raw")
            end

            if instrumental_response
                # With instrumental response, (E,p)
                N_bins = length(instrumental_response_output)
                if N_bins==1
                    plt_inds = [1,1,1]
                elseif N_bins==2
                    plt_inds = [1,2,2]
                elseif N_bins==3
                    plt_inds = [1,2,3]
                else
                    # Algorithm to find suitable indices of Ed to visualize
                    W_gross = dropdims(sum(W_plt, dims=2), dims=2)
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
                        iEd_low = max(1, iEd_mid - abs(iEd_diff))
                        iEd_hig = min(N_bins, iEd_mid + abs(iEd_diff))
                    end
                    plt_inds = [iEd_low, iEd_mid, iEd_hig]
                end

                Ed_low = @sprintf "%.2E" instrumental_response_output[plt_inds[1]]
                Ed_mid = @sprintf "%.2E" instrumental_response_output[plt_inds[2]]
                Ed_hi = @sprintf "%.2E" instrumental_response_output[plt_inds[3]]
                W_plt_low = W_plt[plt_inds[1],:] # w_low(E,p,R,z) but vectorized
                W_plt_mid = W_plt[plt_inds[2],:] # w_mid(E,p,R,z) but vectorized
                W_plt_hi = W_plt[plt_inds[3],:] # w_hi(E,p,R,z) but vectorized

                W_plt_low_Ep = zeros(length(E_array),length(p_array)) # w_low(E,p) = inv(dRdz) ∫ w_low(E,p,R,z) dRdz 
                W_plt_mid_Ep = zeros(length(E_array),length(p_array)) # w_mid(E,p) = inv(dRdz) ∫ w_mid(E,p,R,z) dRdz 
                W_plt_hi_Ep = zeros(length(E_array),length(p_array)) # w_hi(E,p) = inv(dRdz) ∫ w_hi(E,p,R,z) dRdz 
                W_plt_low_Rz = zeros(length(R_array),length(z_array)) # w_low(R,z) = inv(dEdp) ∫ w_low(E,p,R,z) dEdp 
                W_plt_mid_Rz = zeros(length(R_array),length(z_array)) # w_mid(R,z) = inv(dEdp) ∫ w_mid(E,p,R,z) dEdp 
                W_plt_hi_Rz = zeros(length(R_array),length(z_array)) # w_hi(R,z) = inv(dEdp) ∫ w_hi(E,p,R,z) dEdp

                for (i,elem) in enumerate(EpRz_zip_array)
                    iE = elem[1][1]; ip = elem[2][1]; iR = elem[3][1]; iz = elem[4][1]

                    W_plt_low_Ep[iE,ip] += W_plt_low[i]
                    W_plt_mid_Ep[iE,ip] += W_plt_mid[i]
                    W_plt_hi_Ep[iE,ip] += W_plt_hi[i]
                    W_plt_low_Rz[iR,iz] += W_plt_low[i]
                    W_plt_mid_Rz[iR,iz] += W_plt_mid[i]
                    W_plt_hi_Rz[iR,iz] += W_plt_hi[i]
                end

                W_plt_low_Ep_max = @sprintf "%.2E" maximum(W_plt_low_Ep)
                W_plt_mid_Ep_max = @sprintf "%.2E" maximum(W_plt_mid_Ep)
                W_plt_hi_Ep_max = @sprintf "%.2E" maximum(W_plt_hi_Ep)
                W_plt_low_Rz_max = @sprintf "%.2E" maximum(W_plt_low_Rz)
                W_plt_mid_Rz_max = @sprintf "%.2E" maximum(W_plt_mid_Rz)
                W_plt_hi_Rz_max = @sprintf "%.2E" maximum(W_plt_hi_Rz)

                plt_Ep_low = Plots.heatmap(E_array, p_array, transpose(W_plt_low_Ep),title="∫dRdz w($(Ed_low),E,p,R,z), max(w): $(W_plt_low_Ep_max) $(units_inverse(instrumental_response_output_units))", colorbar=false, xlims=extrema(E_array), ylims=extrema(p_array))
                plt_Ep_mid = Plots.heatmap(E_array, p_array, transpose(W_plt_mid_Ep),title="∫dRdz w($(Ed_mid),E,p,R,z), max(w): $(W_plt_mid_Ep_max) $(units_inverse(instrumental_response_output_units))",colorbar=false, xlims=extrema(E_array), ylims=extrema(p_array))
                plt_Ep_hi = Plots.heatmap(E_array, p_array, transpose(W_plt_hi_Ep),title="∫dRdz w($(Ed_hi),E,p,R,z), max(w): $(W_plt_hi_Ep_max) $(units_inverse(instrumental_response_output_units))",colorbar=false, xlims=extrema(E_array), ylims=extrema(p_array))
                plt_Rz_low = Plots.heatmap(R_array, z_array, transpose(W_plt_low_Rz),title="∫dEdp w($(Ed_low),E,p,R,z), max(w): $(W_plt_low_Rz_max) $(units_inverse(instrumental_response_output_units))",colorbar=false)
                plt_Rz_mid = Plots.heatmap(R_array, z_array, transpose(W_plt_mid_Rz),title="∫dEdp w($(Ed_mid),E,p,R,z), max(w): $(W_plt_mid_Rz_max) $(units_inverse(instrumental_response_output_units))",colorbar=false)
                plt_Rz_hi = Plots.heatmap(R_array, z_array, transpose(W_plt_hi_Rz),title="∫dEdp w($(Ed_hi),E,p,R,z), max(w): $(W_plt_hi_Rz_max) $(units_inverse(instrumental_response_output_units))",colorbar=false)

                # Add labels etc to (E,p) plots
                plt_Ep_low = Plots.plot!(plt_Ep_low, xlabel="Energy [keV]", ylabel="Pitch [-]")
                plt_Ep_mid = Plots.plot!(plt_Ep_mid, xlabel="Energy [keV]", ylabel="Pitch [-]")
                plt_Ep_hi = Plots.plot!(plt_Ep_hi, xlabel="Energy [keV]", ylabel="Pitch [-]")

                # Add tokamak wall, labels etc to (R,z) plots
                plt_Rz_low = Plots.plot!(plt_Rz_low, wall.r, wall.z, linewidth=2.0, label="", xlabel="R [m]", ylabel="z [m]", aspect_ratio=:equal)
                plt_Rz_mid = Plots.plot!(plt_Rz_mid, wall.r, wall.z, linewidth=2.0, label="", xlabel="R [m]", ylabel="z [m]", aspect_ratio=:equal)
                plt_Rz_hi = Plots.plot!(plt_Rz_hi, wall.r, wall.z, linewidth=2.0, label="", xlabel="R [m]", ylabel="z [m]", aspect_ratio=:equal)

                plt_EpRz = Plots.plot(plt_Ep_low, plt_Ep_mid, plt_Ep_hi,
                                      plt_Rz_low, plt_Rz_mid, plt_Rz_hi,
                                      layout=(2,3), size=(1200,800), dpi=200, 
                                      fillcolor=cgrad([:white, :yellow, :orange, :red, :black]), titlefontsize=10, colorbar=false, bottom_margin=6Plots.mm,
                                      left_margin=6Plots.mm, right_margin=5Plots.mm)
                png(plt_EpRz,filepath_output*"_EpRz")

                if saveVparaVperpWeights
                    # With instrumental response, (vpara,vperp)
                    W_vel_plt_low = W_vel_plt[plt_inds[1],:] # w_low(vpara,vperp,R,z) but vectorized
                    W_vel_plt_mid = W_vel_plt[plt_inds[2],:] # w_mid(vpara,vperp,R,z) but vectorized
                    W_vel_plt_hi  = W_vel_plt[plt_inds[3],:] # w_hi(vpara,vperp,R,z)  but vectorized

                    W_vel_plt_low_VpaVpe = zeros(length(vpara_array),length(vperp_array)) # w_low(vpara,vperp) = inv(dRdz) ∫ w_low(vpara,vperp,R,z) dRdz 
                    W_vel_plt_mid_VpaVpe = zeros(length(vpara_array),length(vperp_array)) # w_mid(vpara,vperp) = inv(dRdz) ∫ w_mid(vpara,vperp,R,z) dRdz 
                    W_vel_plt_hi_VpaVpe = zeros(length(vpara_array),length(vperp_array)) # w_hi(vpara,vperp) = inv(dRdz) ∫ w_hi(vpara,vperp,R,z) dRdz 
                    W_vel_plt_low_Rz = zeros(length(R_array),length(z_array)) # w_low(R,z) = inv(dVpa*dVpe) ∫ w_low(vpara,vperp,R,z) dVpa*dVpe 
                    W_vel_plt_mid_Rz = zeros(length(R_array),length(z_array)) # w_mid(R,z) = inv(dVpa*dVpe) ∫ w_mid(vpara,vperp,R,z) dVpa*dVpe 
                    W_vel_plt_hi_Rz = zeros(length(R_array),length(z_array)) # w_hi(R,z) = inv(dVpa*dVpe) ∫ w_hi(vpara,vperp,R,z) dVpa*dVpe

                    for (i,elem) in enumerate(VparVperRz_coords)
                        iVpa = elem[1][1]; iVpe = elem[2][1]; iR = elem[3][1]; iz = elem[4][1]

                        W_vel_plt_low_VpaVpe[iVpa,iVpe] += W_vel_plt_low[i]
                        W_vel_plt_mid_VpaVpe[iVpa,iVpe] += W_vel_plt_mid[i]
                        W_vel_plt_hi_VpaVpe[iVpa,iVpe] += W_vel_plt_hi[i]
                        W_vel_plt_low_Rz[iR,iz] += W_vel_plt_low[i]
                        W_vel_plt_mid_Rz[iR,iz] += W_vel_plt_mid[i]
                        W_vel_plt_hi_Rz[iR,iz] += W_vel_plt_hi[i]
                    end

                    W_vel_plt_low_VpaVpe_max = @sprintf "%.2E" maximum(W_vel_plt_low_VpaVpe)
                    W_vel_plt_mid_VpaVpe_max = @sprintf "%.2E" maximum(W_vel_plt_mid_VpaVpe)
                    W_vel_plt_hi_VpaVpe_max = @sprintf "%.2E" maximum(W_vel_plt_hi_VpaVpe)
                    W_vel_plt_low_Rz_max = @sprintf "%.2E" maximum(W_vel_plt_low_Rz)
                    W_vel_plt_mid_Rz_max = @sprintf "%.2E" maximum(W_vel_plt_mid_Rz)
                    W_vel_plt_hi_Rz_max = @sprintf "%.2E" maximum(W_vel_plt_hi_Rz)

                    plt_vel_VpaVpe_low = Plots.heatmap(vpara_array, vperp_array, transpose(W_vel_plt_low_VpaVpe),title="∫dRdz w($(Ed_low),vpa,vpe,R,z), max(w): $(W_vel_plt_low_VpaVpe_max) $(units_inverse(instrumental_response_output_units))", colorbar=false, xlims=extrema(vpara_array), ylims=extrema(vperp_array))
                    plt_vel_VpaVpe_mid = Plots.heatmap(vpara_array, vperp_array, transpose(W_vel_plt_mid_VpaVpe),title="∫dRdz w($(Ed_mid),vpa,vpe,R,z), max(w): $(W_vel_plt_mid_VpaVpe_max) $(units_inverse(instrumental_response_output_units))",colorbar=false, xlims=extrema(vpara_array), ylims=extrema(vperp_array))
                    plt_vel_VpaVpe_hi = Plots.heatmap(vpara_array, vperp_array, transpose(W_vel_plt_hi_VpaVpe),title="∫dRdz w($(Ed_hi),vpa,vpe,R,z), max(w): $(W_vel_plt_hi_VpaVpe_max) $(units_inverse(instrumental_response_output_units))",colorbar=false, xlims=extrema(vpara_array), ylims=extrema(vperp_array))
                    plt_vel_Rz_low = Plots.heatmap(R_array, z_array, transpose(W_vel_plt_low_Rz),title="∫dVpa*dVpe w($(Ed_low),vpa,vpe,R,z), max(w): $(W_vel_plt_low_Rz_max) $(units_inverse(instrumental_response_output_units))",colorbar=false)
                    plt_vel_Rz_mid = Plots.heatmap(R_array, z_array, transpose(W_vel_plt_mid_Rz),title="∫dVpa*dVpe w($(Ed_mid),vpa,vpe,R,z), max(w): $(W_vel_plt_mid_Rz_max) $(units_inverse(instrumental_response_output_units))",colorbar=false)
                    plt_vel_Rz_hi = Plots.heatmap(R_array, z_array, transpose(W_vel_plt_hi_Rz),title="∫dVpa*dVpe w($(Ed_hi),vpa,vpe,R,z), max(w): $(W_vel_plt_hi_Rz_max) $(units_inverse(instrumental_response_output_units))",colorbar=false)

                    # Add labels etc to (vpa,vpe) plots
                    plt_vel_VpaVpe_low = Plots.plot!(plt_vel_VpaVpe_low, xlabel="Vpa [m/s]", ylabel="Vpe [m/s]")
                    plt_vel_VpaVpe_mid = Plots.plot!(plt_vel_VpaVpe_mid, xlabel="Vpa [m/s]", ylabel="Vpe [m/s]")
                    plt_vel_VpaVpe_hi = Plots.plot!(plt_vel_VpaVpe_hi, xlabel="Vpa [m/s]", ylabel="Vpe [m/s]")

                    # Add tokamak wall, labels etc to (R,z) plots
                    plt_vel_Rz_low = Plots.plot!(plt_vel_Rz_low, wall.r, wall.z, linewidth=2.0, label="", xlabel="R [m]", ylabel="z [m]", aspect_ratio=:equal)
                    plt_vel_Rz_mid = Plots.plot!(plt_vel_Rz_mid, wall.r, wall.z, linewidth=2.0, label="", xlabel="R [m]", ylabel="z [m]", aspect_ratio=:equal)
                    plt_vel_Rz_hi = Plots.plot!(plt_vel_Rz_hi, wall.r, wall.z, linewidth=2.0, label="", xlabel="R [m]", ylabel="z [m]", aspect_ratio=:equal)

                    plt_VpaVpeRz = Plots.plot(plt_vel_VpaVpe_low, plt_vel_VpaVpe_mid, plt_vel_VpaVpe_hi,
                                        plt_vel_Rz_low, plt_vel_Rz_mid, plt_vel_Rz_hi,
                                        layout=(2,3), size=(1200,800), dpi=200, 
                                    fillcolor=cgrad([:white, :yellow, :orange, :red, :black]), titlefontsize=10, colorbar=false, bottom_margin=6Plots.mm,
                                    left_margin=6Plots.mm, right_margin=5Plots.mm)
                    png(plt_VpaVpeRz,filepath_output*"_VpaVpeRz")
                end
            end
        end

        verbose && println("Saving weight function matrix... ")
        filepath_output = filepath_output*".jld2"
        list_o_saved_filepaths[iii] = filepath_output # Store a record of the saved file
        myfile_s = jldopen(filepath_output,true,true,false,IOStream)
        ROW, COL, VAL = findnz(Wtot) # Extract rows, columns and values of non-zero indices of (E,p,R,z) weight matrix
        write(myfile_s, "VAL", VAL)
        write(myfile_s, "ROW", ROW)
        write(myfile_s, "COL", COL)
        write(myfile_s, "m_W", Wtot.m)
        write(myfile_s, "n_W", Wtot.n)
        write(myfile_s, "EpRz_coords", EpRz_zip_array) # All the (E,p,R,z) points and indices, corresponding to the columns of the weight matrix
        write(myfile_s, "E_array", vec(E_array))
        write(myfile_s, "p_array", vec(p_array))
        write(myfile_s, "R_array", vec(R_array))
        write(myfile_s, "z_array", vec(z_array))
        write(myfile_s, "reaction", reaction)
        if projVel
            write(myfile_s, "projVel", projVel)
        end
        if saveVparaVperpWeights
            ROW_vel, COL_vel, VAL_vel = findnz(Wtot_vel) # Extract rows, columns and values of non-zero indices of (vpara,vperp,R,z) weight matrix
            write(myfile_s, "VAL_vel", VAL_vel)
            write(myfile_s, "ROW_vel", ROW_vel)
            write(myfile_s, "COL_vel", COL_vel)
            write(myfile_s, "m_W_vel", Wtot_vel.m)
            write(myfile_s, "n_W_vel", Wtot_vel.n)
            write(myfile_s, "VparVperRz_coords", VparVperRz_coords) # All the (vpara,vperp,R,z) points and indices, corresponding to the columns of the weight matrix
            write(myfile_s, "vpara_array", vpara_array)
            write(myfile_s, "vperp_array", vperp_array)
        end
        if instrumental_response
            ROW_raw, COL_raw, VAL_raw = findnz(Wtot_raw) # Extract rows, columns and values of non-zero indices of raw (E,p,R,z) weight matrix
            write(myfile_s, "VAL_raw", VAL_raw)
            write(myfile_s, "ROW_raw", ROW_raw)
            write(myfile_s, "COL_raw", COL_raw)
            write(myfile_s, "m_W_raw", Wtot_raw.m)
            write(myfile_s, "n_W_raw", Wtot_raw.n)
            write(myfile_s, "Ed_array", instrumental_response_output)
            write(myfile_s, "Ed_array_units", instrumental_response_output_units)
            write(myfile_s, "Ed_array_raw", Ed_array)
            write(myfile_s, "Ed_array_raw_units", projVel ? "m_s^-1" : "keV")
            write(myfile_s, "instrumental_response_input", instrumental_response_input)
            write(myfile_s, "instrumental_response_output", instrumental_response_output)
            write(myfile_s, "instrumental_response_matrix", instrumental_response_matrix)
            if saveVparaVperpWeights
                ROW_vel_raw, COL_vel_raw, VAL_vel_raw = findnz(Wtot_vel_raw) # Extract rows, columns and values of non-zero indices of raw (vpara,vperp,R,z) weight matrix
                write(myfile_s, "VAL_vel_raw", VAL_vel_raw)
                write(myfile_s, "ROW_vel_raw", ROW_vel_raw)
                write(myfile_s, "COL_vel_raw", COL_vel_raw)
                write(myfile_s, "m_W_vel_raw", Wtot_vel_raw.m)
                write(myfile_s, "n_W_vel_raw", Wtot_vel_raw.n)
            end
        else
            write(myfile_s, "Ed_array", Ed_array)
            write(myfile_s, "Ed_array_units", projVel ? "m_s^-1" : "keV") # Otherwise, the output abscissa of calc4DWeights.jl is always in m/s or keV
        end
        write(myfile_s, "filepath_thermal_distr", filepath_thermal_distr)
        close(myfile_s)
    else
        verbose && println("Saving debugged quantities... ")
        # WRITE WHATEVER CODE TO SAVE THE DEBUGGED QUANTITIES
    end
end # The end of the iii=1:iiimax for-loop. (to enable computation of several identical orbit weight matrices, to post-analyze MC noise influence for example)

if iiimax>1 && iii_average # If we were supposed to compute more than one weight matrix and an average of all saved files were supposed to be computed... compute it!
    verbose && println("Computing average of all weight matrices... ")
    filepath_first_W = list_o_saved_filepaths[1]
    myfile = jldopen(filepath_first_W,false,false,false,IOStream)
    VAL_total = myfile["VAL"] # Vector containing the non-zero values of the (E,p,R,z) weight matrix
    ROW = myfile["ROW"] # Vector containing the corresponding row elements
    COL = myfile["COL"] # Vector containing the corresponding column elements
    m_W = myfile["m_W"] # Total number of rows (including zero elements not included in R and C) of the (E,p,R,z) weight matrix
    n_W = myfile["n_W"] # Total number of columns (including zero elements not included in R and C) of the (E,p,R,z) weight matrix
    global W_total; W_total = dropzeros(sparse(append!(ROW,m_W),append!(COL,n_W),append!(VAL_total,0.0))) # Make the corresponding sparse matrix (to be able to use addition on several matrices effectively). Add one dummy element to re-create correct size. Don't store it though, only the total matrix size.
    E_array = myfile["E_array"]
    p_array = myfile["p_array"]
    R_array = myfile["R_array"]
    z_array = myfile["z_array"]
    Ed_array = myfile["Ed_array"]
    Ed_array_units = myfile["Ed_array_units"]
    reaction = myfile["reaction"]
    EpRz_coords = myfile["EpRz_coords"] # The indices and (E,p,R,z) coordinates for all the columns of the (E,p,R,z) weight matrix
    if projVel
        projVel = myfile["projVel"]
    end
    if saveVparaVperpWeights
        VAL_vel_total = myfile["VAL_vel"] # Vector containing the non-zero values of the (vpara,vperp,R,z) weight matrix
        ROW_vel = myfile["ROW_vel"] # Vector containing the corresponding row elements
        COL_vel = myfile["COL_vel"] # Vector containing the corresponding column elements
        m_W_vel = myfile["m_W_vel"] # Total number of rows (including zero elements not included in R and C) of the (vpara,vperp,R,z) weight matrix
        n_W_vel = myfile["n_W_vel"] # Total number of columns (including zero elements not included in R and C) of the (vpara,vperp,R,z) weight matrix
        global W_vel_total; W_vel_total = dropzeros(sparse(append!(ROW_vel,m_W_vel),append!(COL_vel,n_W_vel),append!(VAL_vel_total,0.0))) # Make the corresponding sparse matrix (to be able to use addition on several matrices effectively). Add one dummy element to re-create correct size. Don't store it though, only the total matrix size.
        global vpara_array; vpara_array = myfile["vpara_array"]
        global vperp_array; vperp_array = myfile["vperp_array"]
        VparVperRz_coords = myfile["VparVperRz_coords"] # The indices and (vpara,vperp,R,z) coordinates for all the columns of the (vpara,vperp,R,z) weight matrix
    end
    if instrumental_response
        VAL_raw_total = myfile["VAL_raw"] # Vector containing the non-zero raw (without instrumental response) values of the (E,p,R,z) weight matrix
        ROW_raw = myfile["ROW_raw"] # Vector containing the corresponding row elements
        COL_raw = myfile["COL_raw"] # Vector containing the corresponding column elements
        m_W_raw = myfile["m_W_raw"] # Total number of rows (including zero elements not included in R and C) of the raw (without instrumental response) (E,p,R,z) weight matrix
        n_W_raw = myfile["n_W_raw"] # Total number of columns (including zero elements not included in R and C) of the raw (without instrumental response) (E,p,R,z) weight matrix
        global W_raw_total; W_raw_total = dropzeros(sparse(append!(ROW_raw,m_W_raw),append!(COL_raw,n_W_raw),append!(VAL_raw_total,0.0))) # Make the corresponding sparse matrix (to be able to use addition on several matrices effectively). Add one dummy element to re-create correct size. Don't store it though, only the total matrix size.
        Ed_array_raw = myfile["Ed_array_raw"]
        Ed_array_raw_units = myfile["Ed_array_raw_units"]
        instrumental_response_input = myfile["instrumental_response_input"]
        instrumental_response_output = myfile["instrumental_response_output"]
        instrumental_response_matrix = myfile["instrumental_response_matrix"]
        if saveVparaVperpWeights
            VAL_vel_raw_total = myfile["VAL_vel_raw"] # Vector containing the non-zero raw (without instrumental response) values of the (vpara,vperp,R,z) weight matrix
            ROW_vel_raw = myfile["ROW_vel_raw"] # Vector containing the corresponding row elements
            COL_vel_raw = myfile["COL_vel_raw"] # Vector containing the corresponding column elements
            m_W_vel_raw = myfile["m_W_vel_raw"] # Total number of rows (including zero elements not included in R and C) of the raw (without instrumental response) (vpara,vperp,R,z) weight matrix
            n_W_vel_raw = myfile["n_W_vel_raw"] # Total number of columns (including zero elements not included in R and C) of the raw (without instrumental response) (vpara,vperp,R,z) weight matrix
            global W_vel_raw_total; W_vel_raw_total = dropzeros(sparse(append!(ROW_vel_raw,m_W_vel_raw),append!(COL_vel_raw,n_W_vel_raw),append!(VAL_vel_raw_total,0.0))) # Make the corresponding sparse matrix (to be able to use addition on several matrices effectively). Add one dummy element to re-create correct size. Don't store it though, only the total matrix size.
        end
    end
    filepath_thermal_distr = myfile["filepath_thermal_distr"]
    close(myfile)

    for il in 2:length(list_o_saved_filepaths)
        local myfile = jldopen(list_o_saved_filepaths[il],false,false,false,IOStream)
        VAL_il = myfile["VAL"]
        ROW_il = myfile["ROW"]
        COL_il = myfile["COL"]
        m_W_il = myfile["m_W"]
        n_W_il = myfile["n_W"]
        global W_total += dropzeros(sparse(append!(ROW_il,m_W_il),append!(COL_il,n_W_il),append!(VAL_il,0.0))) # Add the il sparse weight matrix to the total sparse weight matrix
        if saveVparaVperpWeights
            VAL_vel_il = myfile["VAL_vel"]
            ROW_vel_il = myfile["ROW_vel"]
            COL_vel_il = myfile["COL_vel"]
            m_W_vel_il = myfile["m_W_vel"]
            n_W_vel_il = myfile["n_W_vel"]
            global W_vel_total += dropzeros(sparse(append!(ROW_vel_il,m_W_vel_il),append!(COL_vel_il,n_W_vel_il),append!(VAL_vel_il,0.0))) # Add the il sparse weight matrix to the total sparse weight matrix
        end
        if instrumental_response
            VAL_raw_il = myfile["VAL_raw"]
            ROW_raw_il = myfile["ROW_raw"]
            COL_raw_il = myfile["COL_raw"]
            m_W_raw_il = myfile["m_W_raw"]
            n_W_raw_il = myfile["n_W_raw"]
            global W_raw_total += dropzeros(sparse(append!(ROW_raw_il,m_W_raw_il),append!(COL_raw_il,n_W_raw_il),append!(VAL_raw_il,0.0))) # Add the il sparse weight matrix to the total sparse weight matrix
            if saveVparaVperpWeights
                VAL_vel_raw_il = myfile["VAL_vel_raw"]
                ROW_vel_raw_il = myfile["ROW_vel_raw"]
                COL_vel_raw_il = myfile["COL_vel_raw"]
                m_W_vel_raw_il = myfile["m_W_vel_raw"]
                n_W_vel_raw_il = myfile["n_W_vel_raw"]
                global W_vel_raw_total += dropzeros(sparse(append!(ROW_vel_raw_il,m_W_vel_raw_il),append!(COL_vel_raw_il,n_W_vel_raw_il),append!(VAL_vel_raw_il,0.0))) # Add the il sparse weight matrix to the total sparse weight matrix
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
    ROW, COL, VAL = findnz(W_total)
    write(myfile,"VAL",VAL)
    write(myfile,"ROW",ROW)
    write(myfile,"COL",COL)
    write(myfile,"m_W", W_total.m)
    write(myfile,"n_W", W_total.n)
    write(myfile,"E_array",E_array)
    write(myfile,"p_array",p_array)
    write(myfile,"R_array",R_array)
    write(myfile,"z_array",z_array)
    write(myfile,"EpRz_coords",EpRz_coords)
    write(myfile,"Ed_array",Ed_array)
    write(myfile,"Ed_array_units",Ed_array_units)
    write(myfile,"reaction",reaction)
    if projVel
        write(myfile,"projVel",projVel)
    end
    if saveVparaVperpWeights
        ROW_vel, COL_vel, VAL_vel = findnz(W_vel_total)
        write(myfile,"VAL_vel",VAL_vel)
        write(myfile,"ROW_vel",ROW_vel)
        write(myfile,"COL_vel",COL_vel)
        write(myfile,"m_W_vel",W_vel_total.m)
        write(myfile,"n_W_vel",W_vel_total.n)
        write(myfile,"vpara_array",vpara_array)
        write(myfile,"vperp_array",vperp_array)
        write(myfile,"VparVperRz_coords",VparVperRz_coords)
    end
    if instrumental_response
        ROW_raw, COL_raw, VAL_raw = findnz(W_raw_total)
        write(myfile,"VAL_raw",VAL_raw)
        write(myfile,"ROW_raw",ROW_raw)
        write(myfile,"COL_raw",COL_raw)
        write(myfile,"m_W_raw",W_raw_total.m)
        write(myfile,"n_W_raw",W_raw_total.n)
        write(myfile,"Ed_array_raw",Ed_array_raw)
        write(myfile,"Ed_array_raw_units",Ed_array_raw_units)
        write(myfile,"instrumental_response_input",instrumental_response_input)
        write(myfile,"instrumental_response_output",instrumental_response_output)
        write(myfile,"instrumental_response_matrix",instrumental_response_matrix)
        if saveVparaVperpWeights
            ROW_vel_raw, COL_vel_raw, VAL_vel_raw = findnz(W_vel_raw_total)
            write(myfile,"VAL_vel_raw",VAL_vel_raw)
            write(myfile,"ROW_vel_raw",ROW_vel_raw)
            write(myfile,"COL_vel_raw",COL_vel_raw)
            write(myfile,"m_W_vel_raw",W_vel_raw_total.m)
            write(myfile,"n_W_vel_raw",W_vel_raw_total.n)
        end
    end
    write(myfile,"filepath_thermal_distr",filepath_thermal_distr)
    close(myfile)
end
println("------ Done! ------")