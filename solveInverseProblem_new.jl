#########################################  solveInverseProblem.jl #########################################

#### Description:
# This script solves an inverse problem on the form s=W*f where 
# 


#
# Please see the start_solveInverseProblem_template.jl file for further input information.

#### Inputs (Units given when defined in script)
# Given via input file start_solveInverseProblem_template.jl, for example. 'template' should be replaced by whatever.

#### Outputs
# -

#### Saved files
# tomographic_reconstruction_[tokamak]_[date and time].jld2
# This saved file will have the keys:
#   F - 

### Other
# 

# Script written by Henrik Järleblad. Last maintained 2024-10-15.
###########################################################################################################

# SECTION 0
verbose && println("Loading Julia packages... ")
using Base.Iterators
using FileIO
using HDF5
using Interpolations
using JLD2
using LinearAlgebra
using SparseArrays
using Statistics
using SCS, Convex
include(folderpath_OWCF*"misc/convert_units.jl")
if plot_solutions || gif_solutions
    using Plots
end
if "collisions" in lowercase.(String.(regularization))
    verbose && println("Collision physics specified as included regularization. OWCF dependencies needed.")
    include(folderpath_OWCF*"extra/dependencies.jl")
    include(folderpath_OWCF*"misc/temp_n_dens.jl")
end
if rescale_W && lowercase(String(rescale_W_F_ref_source))=="file"
    verbose && println("rescale_W set to true and F_ref specified to be loaded from file. OWCF dependencies needed.")
    include(folderpath_OWCF*"extra/dependencies.jl")
end

###########################################################################################################
# SECTION 1: PERFORM SAFETY CHECKS
if !(length(filepaths_S)==length(filepaths_W))
    error("Number of measurement files: $(length(filepaths_S)). Number of weight matrix files: $(length(filepaths_W)). Number of measurement files and number of weight matrix files must match. Please correct and re-try.")
end

###########################################################################################################
# SECTION 2: LOAD MEASUREMENTS

verbose && print("Loading measurements... ")
global S; S = Vector{Vector{Float64}}(undef,length(filepaths_S)) # Pre-allocate measurements for all diagnostics
global S_units; S_units = Vector{String}(undef,length(filepaths_S)) # Pre-allocate units of measurements for all diagnostics
global S_errs; S_errs = Vector{Vector{Float64}}(undef,length(filepaths_S)) # Pre-allocate errors (uncertainties) for all diagnostics
global S_errs_units; S_errs_units = Vector{String}(undef,length(filepaths_S)) # Pre-allocate units of errors (uncertainties) for all diagnostics
global S_abscissas; S_abscissas = Vector{Vector{Float64}}(undef,length(filepaths_S)) # Pre-allocate measurement bin centers for all diagnostics
global S_abscissas_units; S_abscissas_units = Vector{String}(undef,length(filepaths_S)) # Pre-allocate units of measurement bin centers for all diagnostics
for (i,f) in enumerate(filepaths_S)
    f_suffix = lowercase(split(f,".")[end])
    if f_suffix=="jld2"
        myfile = jldopen(f,false,false,false,IOStream)
        read_func = x -> x # Identity function. For .jld2 files
    elseif f_suffix=="h5" || f_suffix=="hdf5"
        myfile = h5open(f,"r")
        read_func = read # HDF5 read function. For .hdf5 (.h5) files
    else
        error("Measurement files should be in either .jld2 or .hdf5 (.h5) file format. Measurement file "*f*"is not. Please correct and re-try.")
    end
    S[i] = read_func(myfile["S"])
    S_units[i] = read_func(myfile["S_unit"])
    S_errs[i] = read_func(myfile["err"])
    S_errs_units[i] = read_func(myfile["err_unit"])
    S_abscissas[i] = read_func(myfile["Ed_array"])
    S_abscissas_units[i] = read_func(myfile["Ed_array_unit"])
    close(myfile)
end
verbose && println("Done!")


global ok # Declare global scope
ok = true # Simplifying print bool dummy variable
verbose && print("Ensuring measurements units consistency between S and err... ")
for (i,S_unit) in enumerate(S_units)
    if !(S_unit==S_errs_units[i])
        ok = false
        verbose && println("")
        verbose && print("Unit of signal $(i): $(S_unit). Unit of err $(i): $(err_units[i]). Attempting to convert err to match S unit... ")
        S_errs[i] = units_conversion_factor(S_errs_units[i],S_unit) .*S_errs[i]
        S_errs_units[i] = S_unit
        verbose && println("Success!")
    end
end
ok && verbose && println("Done!")

###########################################################################################################
# SECTION 3: LOAD WEIGHT FUNCTIONS

verbose && println("Loading weight functions... ")
if isempty(min_array) || isempty(max_array) || isempty(n_array) # If any of them are empty...
    # Use filepaths_W[1] to determine dimensionality
    verbose && println("min_array, max_array or n_array was not specified (empty). Determining problem dimensionality via the first weight function matrix in filepaths_W... ")
    if !isfile(filepaths_W[1])
        error("$(filepaths_W[1]) is not a valid filepath. Please correct and re-try.")
    end
    file_ext = lowercase(split(filepaths_W[1],".")[end]) # Assume final part of filepath after "." is the file extension
    if file_ext=="jld2"
        myfile = jldopen(filepaths_W[1],false,false,false,IOStream)
        if !(lowercase(scriptSources_W[1])=="calc4dweights" || lowercase(scriptSources_W[1])=="calc4dweights.jl")
            W = myfile["W"] # The W file key should be present in all types of weight function files, except...
            DIM = length(size(W))
        else # When we are loading a calc4DWeights.jl output file. Then, we know the dimension to be 5 (one measurement dimension, plus 4 reconstruction space dimensions)
            DIM = 5
        end
        close(myfile)
    elseif file_ext=="hdf5" || file_ext=="h5"
        myfile = h5open(filepaths_W[1],"r")
        W = read(myfile["W"])
        DIM = length(size(W))
        close(myfile)
    else
        error("$(filepaths_W[1]) was not in .jld2, .hdf5 or .h5 file format. Please correct and re-try.")
    end
else # If none of them are empty...
    DIM = 1+length(n_array) # Simply use n_array. One measurement dimension, plus all reconstruction space dimensions
end

global W_inflated; W_inflated = Vector{Array{Float64,DIM}}(undef,length(filepaths_W)) # The (inflated) weight functions, for each diagnostic. PLEASE NOTE! All weight functions need to have the same number of dimensions!
global W_abscissas; W_abscissas = Vector{Vector{Vector{Real}}}(undef,length(filepaths_W)) # The abscissas (measurement bin centers + phase-space grid points), for each set of weight functions
global W_abscissas_units; W_abscissas_units = Vector{Vector{String}}(undef,length(filepaths_W)) # The units of measurement, for each abscissa, for each set of weight functions
for (i,f) in enumerate(filepaths_W)
    f_suffix = lowercase(split(f,".")[end]) # Check file extension for .jld2 or .hdf5 (.h5)
    if f_suffix=="jld2"
        myfile = jldopen(f,false,false,false,IOStream)
        read_func = x -> x # Identity function. For .jld2 files
    elseif f_suffix=="h5" || f_suffix=="hdf5"
        myfile = h5open(f,"r")
        read_func = read # HDF5 read function. For .hdf5 (.h5) files
    else
        error("Weight function files should be in either .jld2 or .hdf5 (.h5) file format. Weight function file "*f*"is not. Please correct and re-try.")
    end
    sswi = lowercase(scriptSources_W[i])
    if sswi=="calc2dweights" || sswi=="calc2dweights.jl" # Account for misinterpretation by user (include .jl file extension)
        if !haskey(myfile,"W")
            error("In the filepaths_W input array, the $(f) file does not have a 'W' file key, even though '$(scriptSources_W[i])' was specified in the scriptSources_W input array. Please correct and re-try.")
        end
        W_inflated[i] = read_func(myfile["W"])
        Ed_vector = read_func(myfile["Ed_array"])
        Ed_vector_unit = read_func(myfile["Ed_array_unit"])
        if !(Ed_vector_unit==S_abscissas_units[i])
            verbose && print("Unit of measurement bin centers of weight function $(i): $(Ed_vector_unit). Unit of measurement bin centers of signal $(i): $(S_abscissas_units[i]). Attempting to convert measurement bin centers of weight function $(i) to match unit of signal $(i)... ")
            Ed_vector = units_conversion_factor(Ed_vector_unit,S_abscissas_units[i]) .*Ed_vector
            Ed_vector_unit = S_abscissas_units[i]
            verbose && println("Success!")
        end
        E_vector = read_func(myfile["E_array"])
        E_vector_unit = "keV" # from calc2DWeights.jl, the energy grid points will always be in keV
        p_vector = read_func(myfile["p_array"])
        p_vector_unit = "-" # from calc2DWeights.jl, the pitch grid points will always be dimensionless
        close(myfile)
        W_abscissas[i] = [Ed_vector, E_vector, p_vector]
        W_abscissas_units[i] = [Ed_vector_unit, E_vector_unit, p_vector_unit]
    elseif sswi=="orbweights_2dto4d" || sswi=="orbweights_2dto4d.jl" || sswi=="calcorbweights" || sswi=="calcorbweights.jl" # Account for misinterpretation by user (include .jl file extension). PLEASE NOTE! THESE FILES SHOULD REALLY BE OUTPUT FILES FROM THE orbWeights_2Dto4D.jl SCRIPT. THE CODE BELOW JUST TAKES INTO ACCOUNT THAT USERS MIGHT BE CONFUSED
        if haskey(myfile,"W2D")
            error("In the filepaths_W input array, the $(f) file is an output file of the calcOrbWeights.jl script. This is not a file that you can specify as input to the solveInverseProblem.jl script. Even though one might think so. Please specify an output file of the OWCF/helper/orbWeights_2Dto4D.jl script instead.")
        end
        if !haskey(myfile,"W")
            error("In the filepaths_W input array, the $(f) file does not have a 'W' file key, even though '$(scriptSources_W[i])' was specified in the scriptSources_W input array. Please correct and re-try.")
        end
        W_inflated[i] = read_func(myfile["W"])
        Ed_vector = read_func(myfile["Ed_array"])
        Ed_vector_unit = read_func(myfile["Ed_array_unit"])
        if !(Ed_vector_unit==S_abscissas_units[i])
            verbose && print("Unit of measurement bin centers of weight function $(i): $(Ed_vector_unit). Unit of measurement bin centers of signal $(i): $(S_abscissas_units[i]). Attempting to convert measurement bin centers of weight function $(i) to match unit of signal $(i)... ")
            Ed_vector = units_conversion_factor(Ed_vector_unit,S_abscissas_units[i]) .*Ed_vector
            Ed_vector_unit = S_abscissas_units[i]
            verbose && println("Success!")
        end
        E_vector = read_func(myfile["E_array"])
        E_vector_unit = "keV" # from calcOrbWeights.jl, the energy grid points will always be in keV
        pm_vector = read_func(myfile["pm_array"])
        pm_vector_unit = "-" # from calcOrbWeights.jl, the pitch maximum grid points will always be dimensionless
        Rm_vector = read_func(myfile["Rm_array"])
        Rm_vector_unit = "m" # from calcOrbWeights.jl, the major radius maximum grid points will always be in meters
        close(myfile)
        W_abscissas[i] = [Ed_vector, E_vector, pm_vector, Rm_vector]
        W_abscissas_units[i] = [Ed_vector_unit, E_vector_unit, pm_vector_unit, Rm_vector_unit]
    elseif sswi=="calc4dweights" || sswi=="calc4dweights.jl"
        # Due to its size, the calc4DWeights.jl case needs special treatment
        # This code can most likely only be run on clusters with a HUGE amount of RAM
        VAL = read_func(myfile["VAL"]) # Vector containing the non-zero values of the (E,p,R,z) weight matrix
        ROW = read_func(myfile["ROW"]) # Vector containing the corresponding row elements
        COL = read_func(myfile["COL"]) # Vector containing the corresponding column elements
        m_W = read_func(myfile["m_W"]) # Total number of rows (including zero elements not included in R and C) of the (E,p,R,z) weight matrix
        n_W = read_func(myfile["n_W"]) # Total number of columns (including zero elements not included in R and C) of the (E,p,R,z) weight matrix
        Ed_vector = read_func(myfile["Ed_array"])
        Ed_vector_unit = read_func(myfile["Ed_array_unit"])
        if !(Ed_vector_unit==S_abscissas_units[i])
            verbose && print("Unit of measurement bin centers of weight function $(i): $(Ed_vector_unit). Unit of measurement bin centers of signal $(i): $(S_abscissas_units[i]). Attempting to convert measurement bin centers of weight function $(i) to match unit of signal $(i)... ")
            Ed_vector = units_conversion_factor(Ed_vector_unit,S_abscissas_units[i]) .*Ed_vector
            Ed_vector_unit = S_abscissas_units[i]
            verbose && println("Success!")
        end
        E_vector = read_func(myfile["E_array"])
        E_vector_unit = "keV" # from calc4DWeights.jl, the energy grid points will always be in keV
        p_vector = read_func(myfile["p_array"])
        p_vector_unit = "-" # from calc4DWeights.jl, the pitch grid points will always be dimensionless
        R_vector = read_func(myfile["R_array"])
        R_vector_unit = "m" # from calc4DWeights.jl, the major radius grid points will always be in meters
        z_vector = read_func(myfile["z_array"])
        z_vector_unit = "m" # from calc4DWeights.jl, the vertical coordinate grid points will always be in meters
        EpRz_coords = read_func(myfile["EpRz_coords"]) # The indices and (E,p,R,z) coordinates for all the columns of the (E,p,R,z) weight matrix (W_2D, see below)
        close(myfile)
        W_2D = dropzeros(sparse(append!(ROW,m_W),append!(COL,n_W),append!(VAL,0.0)))
        W_5D = zeros(length(Ed_vector),length(E_array),length(p_array),length(R_array),length(z_array))
        verbose && println("Re-creating the 5D (Ed,E,p,R,z) weight function matrix... ")
        for iEd in 1:size(W_2D,1)
            for (i,c) in enumerate(EpRz_coords[:])
                W_5D[iEd,c[1][1],c[2][1],c[3][1],c[4][1]] = W_2D[iEd,i]
            end
        end
        W_inflated[i] = W_5D
        W_abscissas[i] = [Ed_vector, E_vector, p_vector, R_vector, z_vector]
        W_abscissas_units[i] = [Ed_vector_unit, E_vector_unit, p_vector_unit, R_vector_unit, z_vector_unit]
    else
        # The general weight functions loading case
        W_inflated[i] = read_func(myfile["W"])
        Ed_vector = read_func(myfile["Ed_array"])
        Ed_vector_unit = read_func(myfile["Ed_array_unit"])
        if !(Ed_vector_unit==S_abscissas_units[i])
            verbose && print("Unit of measurement bin centers of weight function $(i): $(Ed_vector_unit). Unit of measurement bin centers of signal $(i): $(S_abscissas_units[i]). Attempting to convert measurement bin centers of weight function $(i) to match unit of signal $(i)... ")
            Ed_vector = units_conversion_factor(Ed_vector_unit,S_abscissas_units[i]) .*Ed_vector
            Ed_vector_unit = S_abscissas_units[i]
            verbose && println("Success!")
        end
        abscissas = []
        units = []
        if haskey(f,"D1_array")
            append!(abscissas,read_func(myfile["D1_array"]))
            append!(units,read_func(myfile["D1_array_unit"]))
        end
        if haskey(f,"D2_array")
            append!(abscissas,read_func(myfile["D2_array"]))
            append!(units,read_func(myfile["D2_array_unit"]))
        end
        if haskey(f,"D3_array")
            append!(abscissas,read_func(myfile["D3_array"]))
            append!(units,read_func(myfile["D3_array_unit"]))
        end
        if haskey(f,"D4_array")
            append!(abscissas,read_func(myfile["D4_array"]))
            append!(units,read_func(myfile["D4_array_unit"]))
        end
        if haskey(f,"D5_array")
            append!(abscissas,read_func(myfile["D5_array"]))
            append!(units,read_func(myfile["D5_array_unit"]))
        end
        if haskey(f,"D6_array")
            append!(abscissas,read_func(myfile["D6_array"]))
            append!(units,read_func(myfile["D6_array_unit"]))
        end
        close(myfile)
        W_abscissas[i] = Float64.(abscissas) # Transform to Vector{Float64} from Vector{any}
        W_abscissas_units[i] = map(x-> "$(x)",units) # Transform to Vector{String} from Vector{any}
    end
end

###########################################################################################################
# SECTION 4: CHECK THAT ALL ABSCISSAS MATCH IN TERMS OF DIMENSIONS AND UNITS.
# OTHERWISE, CORRECT THEM SO THAT ALL MATCH.

verbose && println("Performing dimension and unit checks for the weight functions... ")
for i in eachindex(filepaths_W)
    if !(length(W_abscissas[i])==length(W_abscissas[1]))
        error("Number of abscissas found in $(filepaths_W[1]): $(length(W_abscissas[1])). Number of abscissas found in $(filepaths_W[i]): $(length(W_abscissas[i])). Number of abscissas must match for all weight functions files. Please correct and re-try.")
    end

    # Check if the units of the measurement bin centers of the weight functions match the units of the measurement bin centers of the signal
    # If not, convert both the abscissa units and weight functions units (since weight functions have units signal/ion)
    if !(W_abscissas_units[i][1]==S_abscissas_units[i])
        verbose && print("Unit of abscissa (measurement bin centers) of $(filepaths_S[i]): $(S_abscissas_units[i]). Unit of FIRST abscissa (measurement bin centers) of $(filepaths_W[i]): $(W_abscissas_units[i][1]). Converting weight function and the first abscissa to match unit of signal abscissa... ")
        ucf = units_conversion_factor(W_abscissas_units[i][1],S_abscissas_units[i])
        W_inflated[i] = (1/ucf) .*W_inflated[i] # Since weight functions have units signal/ion and the signal will have 1/units of the abscissa, multiply by the inverse of the unit conversion factor
        W_abscissas[i][1] = ucf .*W_abscissas[i][1]
        W_abscissas_units[i][1] = S_abscissas_units[i]
    end

    # Check that the units of the reconstruction space abscissas of all weight functions match
    if i>1 # Only compare with previous file if i>1. Otherwise, will cause out-of-bounds error
        for j in eachindex(W_abscissas[i]) # For each reconstruction space dimension
            if (j>1) && !(W_abscissas_units[i][j]==W_abscissas_units[i-1][j]) # If the units of the same dimension for different weight function files don't match (except for the measurement bin center units of course (they will almost always differ), hence the j>1 condition)...
                verbose && println("Units of abscissa $(j) in $(filepaths_W[i]) and $(filepaths_W[i-1]) do not match. Converting abscissa $(j) in $(filepaths_W[i]) to match the unit of abscissa $(j) in $(filepaths_W[i-1])... ")
                ucf = units_conversion_factor(W_abscissas_units[i][j],W_abscissas_units[i-1][j])
                W_abscissas[i][j] = ucf .*W_abscissas[i-1][j]
                W_abscissas_units[i][j] = W_abscissas_units[i-1][j]
            end
        end
    end
end

###########################################################################################################
# SECTION 5: INTERPOLATE ONTO THE GRID SPECIFIED BY min_array, max_array and n_array 
# IF THEY ARE NOT SPECIFIED, USE THE ABSCISSAS OF THE FIRST WEIGHT MATRIX IN filepaths_W

verbose && println("Attemping to interpolate all weight functions onto a common reconstruction space grid... ")
verbose && println("---> Creating interpolation grid (query points)")
if isempty(min_array) || isempty(max_array) || isempty(n_array) # If any of them are empty...
    verbose && println("------> min_array, max_array and/or n_array were not specified. Data grid of $(filepaths_W[1]) will be used")
    query_vecs_n_inds = () # A tuple to hold all query point vectors and their indices. Structure: ((vector,indices),(vector,indices),...) with length equal to reconstruction grid dimension
    for abscissa in W_abscissas[1][2:end] # Use all reconstruction space abscissas of the first weight function file as reconstruction grid...
        query_vecs_n_inds = tuple(query_vecs_n_inds[:]...,collect(zip(abscissa,1:length(abscissa)))) # Add the (vector,indices) pairs one by one  into a big tuple (tuples are immutable, hence the cumbersome code)
    end
else
    verbose && println("------> min_array, max_array and n_array were specified. Utilizing... ")
    query_vecs_n_inds = () # A tuple to hold all query points and their indices. Structure: ((vector,indices),(vector,indices),...)
    for i in eachindex(n_array) # For all reconstruction grid dimensions... 
        query_vecs_n_inds = tuple(query_vecs_n_inds[:]...,collect(zip(collect(range(min_array[i],stop=max_array[i],length=n_array[i])),1:n_array[i]))) # Add the (vector,indices) pairs one by one  into a big tuple (tuples are immutable, hence the cumbersome code)
    end
end
query_points_n_coords = Iterators.product(query_vecs_n_inds...) # Create a long list of all reconstruction space grid points and their coordinates by computing a product between all query point-index vectors. Example structure (if 3 dimensions): [((x1_1,1),(x2_1,1),(x3_1,1)),((x1_2,2),(x2_1,1),(x3_1,1)),...]

for (i,w_inflated) in enumerate(W_inflated)
    w_inflated_interpolated = zeros(tuple(size(w_inflated,1),(length.(query_vecs_n_inds))...)) # Pre-allocate a new inflated weight function, to store the interpolated values
    nodes = () # The nodes of the interpolation object
    for abscissa in W_abscissas[i][2:end] # Use the reconstruction space abscissas of the weight functions as the nodes of the interpolation object...
        nodes = tuple(nodes[:]...,abscissa) # Put them all in a tuple (tuples are immutable, hence the cumbersome code)
    end
    node_coords = CartesianIndices(length.(nodes)) # A trick to specify all node coordinates of an array of general size
    verbose && println("---> Interpolating weight matrix $(i) of $(length(W_inflated))... ")
    for j in eachindex(1:size(w_inflated,1))
        itp = Interpolations.interpolate(nodes,w_inflated[j,node_coords],Gridded(Linear()))
        etp = Interpolations.extrapolate(itp,Flat()) # If outside of interpolation region, use edge values to extrapolate
        for query_point_n_coord in query_points_n_coords
            point = map(x-> x[1],query_point_n_coord) # The point to interpolate at. E.g. (100.0,0.3) in energy (keV),pitch
            coord = map(x-> x[2],query_point_n_coord) # The coordinate of that point. E.g. (53,14)
            w_inflated_interpolated[j,coord] = etp(point...)
        end
    end
    W_inflated[i] = w_inflated_interpolated # Replace the non-interpolated with the interpolated
    W_abscissas[i] = map(x-> map(xx-> xx[1],x), query_vecs_n_inds) # Replace the original grid points with the query points
end

###########################################################################################################
# SECTION 6: RESHAPE ALL INFLATED WEIGHT MATRICES INTO THEIR 2D SHAPE, TO BE USED IN INVERSE PROBLEMS

global W; W = Vector{Array{Float64,2}}(undef,length(filepaths_W)) # The weight matrices, for each diagnostic
for (i,w_inflated) in enumerate(W_inflated)
    ws = size(w_inflated)
    W[i] = reshape(w_inflated,(ws[1],reduce(*,ws[2:end])))
end
W_inflated = nothing # Clear memory. To minimize memory usage

###########################################################################################################
# SECTION 7: DEFINE FUNCTIONS THAT MIGHT NEED TO BE USED SEVERAL TIMES IN LATER SECTIONS, E.G. IF
# - WEIGHT FUNCTIONS ARE TO BE RESCALED AND THE REFERENCE FAST-ION DISTRIBUTION IS TO BE LOADED FROM FILE
# - COLLISIONAL PHYSICS IS TO BE USED AS REGULARIZATION WHEN SOLVING THE INVERSE PROBLEM
# IN FUTURE VERSIONS, ADD MORE FUNCTIONAL CHECKS HERE IF NECESSARY.

if (rescale_W && lowercase(String(rescale_W_F_ref_source))=="file") || ("collisions" in lowercase.(String.(regularization)))
    verbose && println("Defining necessary coordinate space deduction functions... ")
    function is_energy_pitch(w_abscissas::Vector{Vector{T}} where {T<:Real}, w_abscissas_units::Vector{String}; verbose=false, returnExtra=false)
        if length(w_abscissas_units)!=3
            verbose && println("is_energy_pitch(): Number of reconstruction space abscissas is not equal to 2. Returning false... ")
            return (returnExtra ? (false, [], []) : false)
        end

        units_1 = w_abscissas_units[2] # The second dimension of the weight matrix is the first dimension of the reconstruction space
        units_2 = w_abscissas_units[3] # The third dimension of the weight matrix is the second dimension of the reconstruction space
        units_tot = vcat(units_1, units_2)
        
        w_energy_ind = findall(x-> x in ENERGY_UNITS || x in ENERGY_UNITS_LONG, units_tot)
        w_pitch_ind = findall(x-> x in DIMENSIONLESS_UNITS || x in DIMENSIONLESS_UNITS_LONG, units_tot)
        if !(length(w_energy_ind)==1 && length(w_pitch_ind)==1)
            verbose && println("is_energy_pitch(): (E,p) coordinates not found. Returning false... ")
            return (returnExtra ? (false, [], []) : false)
        else
            verbose && println("is_energy_pitch(): (E,p) coordinates confirmed! Returning true... ")
            w_energy_ind .+= 1 # To align it with the original w_abscissas and w_abscissas_units indices
            w_pitch_ind .+= 1 # To align it with the original w_abscissas and w_abscissas_units indices
            return (returnExtra ? (true, w_energy_ind, w_pitch_ind) : true)
        end
    end

    function is_vpara_vperp(w_abscissas::Vector{Vector{T}} where {T<:Real}, w_abscissas_units::Vector{String}; verbose=false, returnExtra=false)
        if length(w_abscissas_units)!=3
            verbose && println("is_vpara_vperp(): Number of reconstruction space abscissas is not equal to 2. Returning false... ")
            return (returnExtra ? (false, [], []) : false)
        end

        units_1 = w_abscissas_units[2] # The second dimension of the weight matrix is the first dimension of the reconstruction space
        units_2 = w_abscissas_units[3] # The third dimension of the weight matrix is the second dimension of the reconstruction space
        units_tot = vcat(units_1, units_2)
        w_speed_inds = findall(x-> units_are_speed(x), units_tot)

        if !(length(w_speed_inds)==2)
            verbose && println("is_vpara_vperp(): (vpara,vperp) coordinates not found. Returning false... ")
            return (returnExtra ? (false, [], []) : false)
        else
            w_speed_inds .+= 1 # To align it with the original w_abscissas and w_abscissas_units indices
            verbose && print("is_vpara_vperp(): (vpara,vperp) coordinates found! Distinguishing (vpara, vperp) arrays...")
            w_vel_arrays = w_abscissas[w_speed_inds]
            if minimum(w_vel_arrays[1])<0 && minimum(w_vel_arrays[2])>0 # vpara can be negative. vperp cannot
                verbose && println("ok!")
                w_vpara_ind = [w_speed_inds[1]] # Order of abscissas in w_vel_arrays is the same as the order in w_speed_inds
                w_vperp_ind = [w_speed_inds[2]] # Order of abscissas in w_vel_arrays is the same as the order in w_speed_inds
            elseif minimum(w_vel_arrays[2]<0 && minimum(w_vel_arrays[1]>0)) # vpara can be negative. vperp cannot
                verbose && println("ok!")
                w_vpara_ind = [w_speed_inds[2]] # Order of abscissas in w_vel_arrays is the same as the order in w_speed_inds
                w_vperp_ind = [w_speed_inds[1]] # Order of abscissas in w_vel_arrays is the same as the order in w_speed_inds
            else
                verbose && println("")
                @warn "Could not distinguish (vpara,vperp) arrays from weight function abscissas. Assuming abscissa with index $(w_speed_inds[1]) to be vpara and abscissa with index $(w_speed_inds[2]) to be vperp."
                w_vpara_ind = [w_speed_inds[1]] # Order of abscissas in w_vel_arrays is the same as the order in w_speed_inds
                w_vperp_ind = [w_speed_inds[2]] # Order of abscissas in w_vel_arrays is the same as the order in w_speed_inds
            end
            verbose && println("is_vpara_vperp(): Returning true... ")
            return (returnExtra ? (true, w_vpara_ind, w_vperp_ind) : true)
        end
    end

    function is_EpRz(w_abscissas::Vector{Vector{T}} where {T<:Real}, w_abscissas_units::Vector{String}; verbose=false, returnExtra=false)
        if length(w_abscissas_units)!=5
            verbose && println("is_EpRz(): Number of reconstruction space abscissas is not equal to 4. Returning false... ")
            return (returnExtra ? (false, [], [], [], [], [], []) : false)
        end
        units_1 = w_abscissas_units[2] # The second dimension of the weight matrix is the first dimension of the reconstruction space
        units_2 = w_abscissas_units[3] # The third dimension of the weight matrix is the second dimension of the reconstruction space
        units_3 = w_abscissas_units[4] # The fourth dimension of the weight matrix is the third dimension of the reconstruction space
        units_4 = w_abscissas_units[5] # The fifth dimension of the weight matrix is the fourth dimension of the reconstruction space

        units_tot = vcat(units_1, units_2, units_3, units_4)
        w_energy_ind = findall(x-> x in ENERGY_UNITS || x in ENERGY_UNITS_LONG, units_tot)
        w_pitch_ind = findall(x-> x in DIMENSIONLESS_UNITS || x in DIMENSIONLESS_UNITS_LONG, units_tot)
        w_Rz_inds = findall(x-> x in LENGTH_UNITS || x in LENGTH_UNITS_LONG, units_tot)

        if !(length(w_energy_ind)==1 && length(w_pitch_ind)==1 && length(w_Rz_inds)==2)
            verbose && println("is_EpRz(): (E,p,R,z) coordinates not found. Returning false... ")
            return (returnExtra ? (false, [], [], [], [], [], []) : false)
        else
            w_energy_ind .+= 1 # To align it with the original w_abscissas and w_abscissas_units indices
            w_pitch_ind .+= 1 # To align it with the original w_abscissas and w_abscissas_units indices
            w_Rz_inds .+= 1 # To align it with the original w_abscissas and w_abscissas_units indices

            verbose && print("is_EpRz(): (E,p,R,z) coordinates found! Distinguishing (R,z) arrays... ")
            w_RnZ_arrays = w_abscissas[w_Rz_inds]
            if minimum(w_RnZ_arrays[2])<0 && minimum(w_RnZ_arrays[1])>0 # If the second LENGTH_UNITS abscissa has negative elements, and the first one does not..
                verbose && println("ok!")
                w_R_ind = w_Rz_inds[1] # Order of abscissas in w_RnZ_arrays is the same as the order in w_Rz_inds
                w_z_ind = w_Rz_inds[2] # Order of abscissas in w_RnZ_arrays is the same as the order in w_Rz_inds
                R_of_interests = w_RnZ_arrays[1] # The first LENGTH_UNITS abscissa is very likely to be the R grid points...
                z_of_interests = w_RnZ_arrays[2] # ...and the second LENGTH_UNITS abscissa is very likely to be the z grid points
            elseif minimum(w_RnZ_arrays[1])<0 && minimum(w_RnZ_arrays[2])>0 # If it's the other way around...
                verbose && println("ok!")
                w_R_ind = w_Rz_inds[2] # Order of abscissas in w_RnZ_arrays is the same as the order in w_Rz_inds
                w_z_ind = w_Rz_inds[1] # Order of abscissas in w_RnZ_arrays is the same as the order in w_Rz_inds
                R_of_interests = w_RnZ_arrays[2] # ...it's very likely to be the other way around.
                z_of_interests = w_RnZ_arrays[1] # ...it's very likely to be the other way around.
            else
                verbose && println("")
                @warn "Could not deduce (R,z) arrays from weight function abscissas. Assuming abscissa with index $(w_Rz_inds[1]) to be R and abscissa with index $(w_Rz_inds[2]) to be z."
                w_R_ind = w_Rz_inds[1] # Order of abscissas in w_RnZ_arrays is the same as the order in w_Rz_inds
                w_z_ind = w_Rz_inds[2] # Order of abscissas in w_RnZ_arrays is the same as the order in w_Rz_inds
                R_of_interests = w_RnZ_arrays[1]
                z_of_interests = w_RnZ_arrays[2]
            end
            verbose && println("is_EpRz(): Returning true... ")
            return (returnExtra ? (true, w_energy_ind, w_pitch_ind, w_R_ind, w_z_ind, R_of_interests, z_of_interests) : true)
        end
    end

    function get_energy_abscissa(w_abscissas::Vector{Vector{T}} where {T<:Real}, w_abscissas_units::Vector{String}; verbose=false, returnExtra=false)
        # Assume there is only one energy abscissa
        w_energy_ind = findfirst(x-> x in ENERGY_UNITS || x in ENERGY_UNITS_LONG, w_abscissas_units)
        if isnothing(w_energy_ind)
            verbose && println("No energy abscissa found in input abscissas! Returning nothing... ")
            return returnExtra ? (nothing,nothing) : nothing
        end

        return returnExtra ? (w_abscissas[w_energy_ind], w_energy_ind) : w_abscissas[w_energy_ind]
    end

    function get_pitch_abscissa(w_abscissas::Vector{Vector{T}} where {T<:Real}, w_abscissas_units::Vector{String}; verbose=false, returnExtra=false)
        w_dimensionless_inds = findall(x-> x in DIMENSIONLESS_UNITS || x in DIMENSIONLESS_UNITS_LONG, w_abscissas_units)
    
        if isempty(w_dimensionless_inds)
            verbose && println("No dimensionless abscissa found in input abscissas! Returning nothing...")
            return returnExtra ? (nothing,nothing) : nothing
        end
    
        w_pitch_ind = 0
        for w_dimensionless_ind in w_dimensionless_inds
            w_abscissa = w_abscissas[w_dimensionless_ind]
            if maximum(w_abscissa)<=1 && minimum(w_abscissa)>=-1
                w_pitch_ind = w_dimensionless_ind
                break # Assume there is only one pitch-like coordinate
            end
        end
    
        if w_pitch_ind==0
            verbose && println("No pitch-like dimensionless abscissa found in input abscissas! Returning nothing...")
            return returnExtra ? (nothing,nothing) : nothing
        end
    
        return returnExtra ? (w_abscissas[w_pitch_ind], w_pitch_ind) : w_abscissas[w_pitch_ind]
    end

    function get_vpara_abscissa(w_abscissas::Vector{Vector{T}} where {T<:Real}, w_abscissas_units::Vector{String}; verbose=false, returnExtra=false)
        w_speed_inds = findall(x-> units_are_speed(x), w_abscissas_units)
        if isempty(w_speed_inds)
            verbose && println("No speed abscissa found in input abscissas! Returning nothing...")
            return returnExtra ? (nothing,nothing) : nothing
        end
        if length(w_speed_inds)>2
            verbose && println("Impossible to determine vpara since more than two abscissas with unit of measurement 'speed' was found. Returning nothing... ")
            return returnExtra ? (nothing,nothing) : nothing
        end

        if length(w_speed_inds)==1
            verbose && println("One abscissa with unit of measurement 'speed' was found. Assuming it is vpara. Returning... ")
            return returnExtra ? (w_abscissas[w_speed_inds], w_speed_inds[1]) : w_abscissas[w_speed_inds]
        end

        # w_speed_inds must be length 2
        w_speed_abscissa_1 = w_speed_inds[1]
        w_speed_abscissa_2 = w_speed_inds[2]
        if minimum(w_speed_abscissa_1)<0 # Must be vpara since vperp>= always holds
            return returnExtra ? (w_speed_abscissa_1, w_speed_inds[1]) : w_speed_abscissa_1
        elseif minimum(w_speed_abscissa_2)<0 # -||-
            return returnExtra ? (w_speed_abscissa_2, w_speed_inds[2]) : w_speed_abscissa_2
        else # Ambiguous
            @warn "Cannot deduce vpara from w_abscissas and w_abscissas_units. Returning abscissas with smallest index as vpara."
            return returnExtra ? (w_speed_abscissa_1, w_speed_inds[1]) : w_speed_abscissa_1
        end
    end

    function get_vperp_abscissa(w_abscissas::Vector{Vector{T}} where {T<:Real}, w_abscissas_units::Vector{String}; verbose=false, returnExtra=false)
        w_speed_inds = findall(x-> units_are_speed(x), w_abscissas_units)
        # DO SOMETHING HERE TO DEDUCE vperp!!
        w_vpara_abscissa, w_vpara_ind = get_vpara_abscissa(w_abscissas, w_abscissas_units; verbose=verbose, returnExtra=returnExtra)
        # DO SOMETHING HERE TO DEDUCE vperp!!!
        
    end
end
###########################################################################################################
# SECTION 8: IF rescale_W WAS SET TO true IN THE START FILE, LOAD OR COMPUTE FAST-ION DISTRIBUTION(S).
# USE THIS/THESE DISTRIBUTION(S) TOGETHER WITH THE WEIGHT MATRICES TO COMPUTE REFERENCE MEASUREMENTS.
# COMPUTE RESCALING FACTORS TO BE ABLE TO RESCALE WEIGHT FUNCTIONS TO HAVE THE REFERENCE MEASUREMENTS MATCH 
# THE EXPERIMENTAL MEASUREMENTS.

if rescale_W
    verbose && println("Rescaling weight functions... ")
    if lowercase(String(rescale_W_type))=="mean"
        verbose && println("---> mean(S)/mean(W*F_ref) will be used to rescale weight functions... ")
        rescale_func = Statistics.mean # Set the mean() function as the rescale_func() function
    elseif lowercase(String(rescale_W_type))=="maximum" || lowercase(String(rescale_W_type))=="max" 
        verbose && println("---> max(S)/max(W*F_ref) will be used to rescale weight functions... ")
        rescale_func = Base.maximum # Set the maximum() function as the rescale_func() function
    else
        error("rescale_W_type=$(rescale_W_type). Currently supported options include :MEAN and :MAXIMUM. Please correct and re-try.")
    end
    if lowercase(String(rescale_W_F_ref_source))=="file"
        currently_supported_abscissas = "(E,p), (vpara, vperp), (E,p,R,z)" # UPDATE THIS WHEN NEEDED!
        standard_rescale_W_error = "rescale_W_F_ref_source=$(rescale_W_F_ref_source) was specified, but weight function abscissas did not match any currently supported groups. Currently supported abscissas for weight functions include $(currently_supported_abscissas). Please change rescale_W_F_ref_source to :GAUSSIAN, or specify filepaths_W to be filepath(s) to file(s) containing weight functions with supported abscissas."
        w_abscissas_units = W_abscissas_units[1] # Units of abscissas of first weight function matrix are the units of the abscissas of all weight functions matrices, because of section 5
        w_abscissas = W_abscissas[1] # Abscissas of first weight function matrix are the abscissas of all weight functions matrices, because of section 5

        verbose && println("---> rescale_W_F_ref_source=$(rescale_W_F_ref_source) was specified... ")
        verbose && println("------> Checking weight function reconstruction space... ")
        Ep_bool, w_energy_ind, w_pitch_ind = is_energy_pitch(w_abscissas, w_abscissas_units; verbose=verbose, returnExtra=true)
        VEL_bool, w_vpara_ind, w_vperp_ind = is_vpara_vperp(w_abscissas, w_abscissas_units; verbose=verbose, returnExtra=true)
        EpRz_bool, w_energy_ind, w_pitch_ind, w_R_ind, w_z_ind, R_of_interests, z_of_interests = is_EpRz(w_abscissas, w_abscissas_units; verbose=verbose, returnExtra=true)
        # ADD MORE CHECKS HERE IN FUTURE VERSIONS, IF NEEDED <---------------------------------------------------------------------
        if Ep_bool
            w_rec_space = "(E,p)"
            R_of_interests = vcat(units_conversion_factor(R_of_interest_unit,"m")*R_of_interest) # R_of_interests will be the same, regardless of 2D coordinates
            z_of_interests = vcat(units_conversion_factor(z_of_interest_unit,"m")*z_of_interest) # z_of_interests will be the same, regardless of 2D coordinates
        elseif VEL_bool
            w_rec_space = "(vpara,vperp)"
            R_of_interests = vcat(units_conversion_factor(R_of_interest_unit,"m")*R_of_interest) # R_of_interests will be the same, regardless of 2D coordinates
            z_of_interests = vcat(units_conversion_factor(z_of_interest_unit,"m")*z_of_interest) # z_of_interests will be the same, regardless of 2D coordinates
        elseif EpRz_bool
            w_rec_space = "(E,p,R,z)"
            R_of_interests = units_conversion_factor(w_abscissas_units[w_R_ind],"m") .*R_of_interests # Need unit of measurement meter for CDFto4D function (if that will be used)
            z_of_interests = units_conversion_factor(w_abscissas_units[w_z_ind],"m") .*z_of_interests # Need unit of measurement meter for CDFto4D function (if that will be used)
        else
            error(standard_rescale_W_error)
        end

        verbose && println("------> Loading F_ref from file... ")
        file_ext = lowercase(split(rescale_W_F_file_path,".")[end]) # Assume final part of filepath after "." is the file extension
        if file_ext=="jld2"
            F_ref, E_ref, p_ref, R_ref, z_ref = JLD2to4D(rescale_W_F_file_path)
        elseif file_ext=="hdf5" || file_ext=="h5"
            F_ref, E_ref, p_ref, R_ref, z_ref = h5to4D(rescale_W_F_file_path;backwards=h5file_of_nonJulia_origin)
        elseif file_ext=="cdf"
            F_ref, E_ref, p_ref, R_ref, z_ref = CDFto4D(rescale_W_F_file_path, R_of_interests, z_of_interests; btipsign=btipsign)
        else
            error("Input variable 'rescale_W_F_file_path' has unknown file extension ($(file_ext)). Accepted file extensions are .jld2, .hdf5, .h5 and .cdf. Please correct and re-try.")
        end

        if w_rec_space=="(E,p)"
            verbose && println("------> Computing W*F_ref assuming weight functions are given in (E,p) coordinates (or similar)... ")
            w_energy_unit = w_abscissas_units[w_energy_ind]
            w_pitch_unit = w_abscissas_units[w_pitch_ind]
            w_energy_abscissa = w_abscissas[w_energy_ind]; nwE = length(w_energy_abscissa)
            w_pitch_abscissa = w_abscissas[w_pitch_ind]; nwp = length(w_pitch_abscissa)

            query_points_n_coords = Iterators.product(zip(w_energy_abscissa,1:nwE),zip(w_pitch_abscissa,1:nwp),zip(R_of_interests,1:1),zip(z_of_interests,1:1)) # R and z should already be 1-element Vectors with unit of measurement 'm'
            F_ref_interpolated = zeros(nwE,nwp,1,1) # Pre-allocate the interpolated f(E,p) distribution

            E_ref = units_conversion_factor("keV",w_energy_unit) .*E_ref # Convert the energy units from keV (which we know the JLD2to4D, h5to4D and CDFto4D functions output) to w_energy_unit
            F_ref = units_conversion_factor("keV^-1",units_inverse(w_energy_unit)) .*F_ref # Also, convert the fast-ion distribution using the inverse of the w_energy_unit

            nodes = (E_ref,p_ref,R_ref,z_ref)
            itp = Interpolations.interpolate(nodes,F_ref,Gridded(Linear()))
            etp = Interpolations.extrapolate(itp,Flat()) # If outside of interpolation region, use edge values to extrapolate
            for query_point_n_coord in query_points_n_coords
                point = map(x-> x[1],query_point_n_coord) # The point to interpolate at. E.g. (100.0,0.3) in energy (keV),pitch
                coord = map(x-> x[2],query_point_n_coord) # The coordinate of that point. E.g. (53,14)
                F_ref_interpolated[coord] = etp(point...)
            end
            F_ref_interpolated = dropdims(F_ref_interpolated,dims=(3,4)) # From shape (nwE,nwp,1,1) to (nwE,nwp)
            F_ref_interp_1D = reshape(F_ref_interpolated,(nwE*nwp,))

            rescale_W_factors = ones(length(filepaths_W)) # Pre-allocate weight re-scale factors
            for i in eachindex(filepaths_W)
                WF = W[i]*F_ref_interp_1D # The weight function matrix (2D) multiplied with the (1D/vectorized) reference fast-ion distribution
                rescale_W_factors[i] = rescale_func(S[i])/rescale_func(WF) # S[i] is the experimental data vector of (fast-ion) diagnostic 'i'
            end
        elseif w_rec_space=="(vpara,vperp)"
            verbose && println("------> Using reference fast-ion distribution assuming weight functions are given in (vpara,vperp) coordinates (or similar)... ")
            w_vpara_unit = w_abscissas_units[w_vpara_ind]
            w_vperp_unit = w_abscissas_units[w_vperp_ind]
            w_vpara_abscissa = w_abscissas[w_vpara_ind]; nwvpa = length(w_vpara_abscissa)
            w_vperp_abscissa = w_abscissas[w_vperp_ind]; nwvpe = length(w_vpara_abscissa)
            nfE = length(E_ref); nfp = length(p_ref)

            nodes = (E_ref,p_ref,R_ref,z_ref)
            itp = Interpolations.interpolate(nodes,F_ref,Gridded(Linear()))
            etp = Interpolations.extrapolate(itp,Flat()) # If outside of interpolation region, use edge values to extrapolate
            query_points_n_coords = Iterators.product(zip(E_ref,1:nfE),zip(p_ref,1:nfp),zip(R_of_interests,1:1),zip(z_of_interests,1:1))
            F_ref_interpolated = zeros(nfE,nfp,1,1) # Pre-allocate the interpolated f(E,p) distribution
            for query_point_n_coord in query_points_n_coords
                point = map(x-> x[1],query_point_n_coord) # The point to interpolate at. E.g. (100.0,0.3) in energy (keV),pitch
                coord = map(x-> x[2],query_point_n_coord) # The coordinate of that point. E.g. (53,14)
                F_ref_interpolated[coord] = etp(point...)
            end
            F_ref_interpolated = dropdims(F_ref_interpolated,dims=(3,4)) # From shape (nfE,nfp,1,1) to (nfE,nfp)

            F_ref_VEL, vpara_ref, vperp_ref = Ep2VparaVperp(E_ref, p_ref, F_ref_interpolated; my_gcp=FI_species, needJac=true, returnAbscissas=true)
            vpara_ref = units_conversion_factor("m_s^-1",w_vpara_unit) .*vpara_ref # Match the units of w_abscissas
            vperp_ref = units_conversion_factor("m_s^-1",w_vperp_unit) .*vperp_ref # Match the units of w_abscissas
            F_ref_VEL = units_conversion_factor("m^-2_s^2",units_inverse(w_vpara_unit)*"_"*units_inverse(w_vperp_unit)) .* F_ref_VEL

            nodes = (vpara_ref,vperp_ref)
            itp = Interpolations.interpolate(nodes,F_ref_VEL,Gridded(Linear()))
            etp = Interpolations.extrapolate(itp,Flat())
            query_points_n_coords = Iterators.product(zip(w_vpara_abscissa,1:nwvpa),zip(w_vperp_abscissa,1:nwvpe))
            F_ref_VEL_interp = zeros(nwvpa,nwvpe)
            for query_point_n_coord in query_points_n_coords
                point = map(x-> x[1],query_point_n_coord) # The point to interpolate at. E.g. (-0.6e6,0.2e6) in vpara (m/s), vperp (m/s)
                coord = map(x-> x[2],query_point_n_coord) # The coordinate of that point. E.g. (53,14)
                F_ref_VEL_interp[coord] = etp(point...)
            end
            F_ref_VEL_interp_1D = reshape(F_ref_VEL_interp,(nwvpa*nwvpe,))

            rescale_W_factors = ones(length(filepaths_W))
            for i in eachindex(filepaths_W)
                WF = W[i]*F_ref_VEL_interp_1D # The weight function matrix (2D) multiplied with the (1D/vectorized) reference fast-ion distribution
                rescale_W_factors[i] = rescale_func(S[i])/rescale_func(WF) # S[i] is the experimental data vector of (fast-ion) diagnostic 'i'
            end
        elseif w_rec_space=="(E,p,R,z)"
            verbose && println("------> Using reference fast-ion distribution assuming weight functions are given in (E,p,R,z) coordinates (or similar)... ")
            w_energy_unit = w_abscissas_units[w_energy_ind]
            w_pitch_unit = w_abscissas_units[w_pitch_ind]
            w_R_unit = w_abscissas_units[w_R_ind]
            w_z_unit = w_abscissas_units[w_z_ind]
            w_energy_abscissa = w_abscissas[w_energy_ind]; nwE = length(w_energy_abscissa)
            w_pitch_abscissa = w_abscissas[w_pitch_ind]; nwp = length(w_pitch_abscissa)
            w_R_abscissa = w_abscissas[w_R_ind]; nwR = length(w_R_abscissa)
            w_z_abscissa = w_abscissas[w_z_ind]; nwz = length(w_z_abscissa)

            E_ref = units_conversion_factor("keV",w_energy_unit) .*E_ref # Convert the energy units from keV (which we know the JLD2to4D, h5to4D and CDFto4D functions output) to w_energy_unit
            R_ref = units_conversion_factor("m",w_R_unit) .*R_ref # Convert the R units from m (which we know the JLD2to4D, h5to4D and CDFto4D functions output) to w_R_unit
            z_ref = units_conversion_factor("m",w_z_unit) .*z_ref # Convert the z units from m (which we know the JLD2to4D, h5to4D and CDFto4D functions output) to w_z_unit
            F_ref = units_conversion_factor("keV^-1_m^-3",units_inverse(w_energy_unit)*"_"*units_inverse(w_R_unit)*"_"*units_inverse(w_z_unit)) .*F_ref # Also, convert the fast-ion distribution from keV^-1_m^-3 to the inverse of w_energy_unit, w_R_unit and w_z_unit 

            query_points_n_coords = Iterators.product(zip(w_energy_abscissa,1:nwE),zip(w_pitch_abscissa,1:nwp),zip(w_R_abscissa,1:nwR),zip(w_z_abscissa,1:nwz))
            F_ref_interpolated = zeros(nwE,nwp,nwR,nwz) # Pre-allocate the interpolated f(E,p,R,z) distribution

            nodes=(E_ref,p_ref,R_ref,z_ref)
            itp = Interpolations.interpolate(nodes,F_ref,Gridded(Linear()))
            etp = Interpolations.extrapolate(itp,Flat())
            for query_point_n_coord in query_points_n_coords
                point = map(x-> x[1],query_point_n_coord) # The point to interpolate at. E.g. (120.0,0.64,3.1,0.2) in E [keV], p [-], R [m], z [m] 
                coord = map(x-> x[2],query_point_n_coord) # The coordinate of that point. E.g. (53,34,15,16)
                F_ref_interpolated[coord] = etp(point...)
            end
            F_ref_interp_1D = reshape(F_ref_interpolated,(nwE*nwp*nwR*nwz,))

            rescale_W_factors = zeros(length(filepaths_W))
            for i in eachindex(filepaths_W)
                WF = W[i]*F_ref_interp_1D # The weight function matrix (2D) multiplied with the (1D/vectorized) reference fast-ion distribution
                rescale_W_factors[i] = rescale_func(S[i])/rescale_func(WF) # S[i] is the experimental data vector of (fast-ion) diagnostic 'i'
            end
        elseif false # ADD MORE CHECKS HERE IN FUTURE VERSIONS, IF NEEDED <---------------------------------------------------------------------
        else
            error(standard_rescale_W_error)
        end
    elseif lowercase(String(rescale_W_F_ref_source))=="gaussian"
        w_abscissas_units = W_abscissas_units[1] # Units of abscissas of first weight function matrix are the units of the abscissas of all weight functions matrices, because of section 5
        w_abscissas = W_abscissas[1] # Abscissas of first weight function matrix are the abscissas of all weight functions matrices, because of section 5

        verbose && println("---> rescale_W_F_ref_source=$(rescale_W_F_ref_source) was specified... ")
        verbose && println("---> The average of the WF-products for $(5*10^(length(w_abscissas)-1)) Gaussian reference FI distributions will be used to rescale the weight functions...")
        
        rec_space_lengths = length.(w_abscissas[2:end])
        rec_space_diffs = [x[1] for x in diff.(w_abscissas[2:end])] # Assume equidistant grids in reconstruction space
        rec_space_extrema = [[x[1],x[2]] for x in extrema.(w_abscissas[2:end])]
        F_ref_abscissas = [collect(range(x[1],stop=x[2],length=10)) for x in rec_space_extrema] # Same abscissas as reconstruction space, but a 10x10x...x10 grid
        std_abscissas = [collect(range(0.05*x[1],stop=0.25*x[1],length=5)) for x in diff.(rec_space_extrema)] # For each reconstruction space dimension, a set of standard deviations from 0.05 to 0.25 times the width of the dimension
        gauss_abscissas = vcat(F_ref_abscissas,std_abscissas) # Put together rec space abscissas and their standard deviations
        gauss_coords = Iterators.product(gauss_abscissas...) # Create an iterator of the product between all elements
        rec_space_minima = [x[1] for x in rec_space_extrema]
        rec_space_maxima = [x[2] for x in rec_space_extrema]
        WF_avg = [zeros(size(W[i],1)) for i in 1:length(filepaths_W)] # Pre-allocate WF_avg. The number of diagnostic measurement bins for each diagnostic
        for gauss_coord in gauss_coords # For all the Gaussians...
            l = length(gauss_coord)
            mu = gauss_coord[1:l/2] # The first half will be from F_ref_abscissas
            std = gauss_coord[(l/2)+1:end] # The second half will be from std_abscissas
            F_gauss = gaussian(mu, std; mn=rec_space_minima, mx=rec_space_maxima, n=rec_space_lengths, floor_level=0.001) # Set all valus below 0.001*maximum(F_gauss) to zero
            F_gauss = (rescale_W_F_gaus_N_FI_tot/sum(reduce(*,rec_space_diffs) .*F_gauss)) .*F_gauss # Rescale F_gauss to integrate to rescale_W_F_gaus_N_FI_tot
            F_gauss_1D = reshape(F_gauss,(reduce(*,rec_space_lengths),))
            for i in 1:length(filepaths_W) # For all diagnostics...
                WF_avg[i] += W[i]*F_gauss_1D # The weight function matrix (2D) multiplied with the (1D/vectorized) reference fast-ion distribution
            end
        end
        WF_avg ./= length(gauss_coords)
        rescale_W_factors = [rescale_func(S[i])/rescale_func(wf) for (i,wf) in enumerate(WF_avg)] # WF_avg will have same length as S and W (and filepaths_S and filepaths_W)
    else
        error("'rescale_W' was set to true but 'rescale_W_F_ref_source' was not specified correctly. Currently supported options include :FILE and :GAUSSIAN. Please correct and re-try.")
    end
else
    verbose && println("Weight functions will NOT be rescaled... ")
    rescale_W_factors = ones(length(filepaths_W))
end

###########################################################################################################
# SECTION 9: IF :COLLISIONS WAS INCLUDED IN THE regularization INPUT VARIABLE, AND THE RECONSTRUCTION 
# SPACE IS A SPACE IN WHICH COLLISIONS REGULARIZATION IS SUPPORTED, COMPUTE SLOWING-DOWN BASIS FUNCTIONS

if "collisions" in lowercase.(String.(regularization))
    COMPATIBLE_OPTIONS = [is_energy_pitch, is_vpara_vperp] # ADD MORE CHECKS HERE IN FUTURE VERSIONS, IF NEEDED <---------------------------------------------------------------------
    verbose && print(":COLLISIONS included in 'regularization' input variable. Performing safety checks... ")
    if sum([f(W_abscissas[1],W_abscissas_units[1]) for f in COMPATIBLE_OPTIONS])==0 # Is the reconstruction space NOT in any of the COMPATIBLE_OPTIONS?
        error(":COLLISIONS was specified in 'regularization' input variable, but reconstruction space is not compatible. Current compatible options include (E,p) and (vpara,vperp). Please correct and re-try.")
    end
    if !(length(regularization_thermal_ion_temp)==length(regularization_thermal_ion_dens)==length(regularization_thermal_ion_species))
        error("The lengths of the 'regularization_thermal_ion_temp', 'regularization_thermal_ion_dens' and 'regularization_thermal_ion_species' input variables were not equal. Please correct and re-try.")
    end
    verbose && println("ok!")
    verbose && print("---> Loading magnetic equilibrium from $(regularization_equil_filepath)... ")
    M, wall = read_geqdsk(regularization_equil_filepath,clockwise_phi=false) # Assume the phi-direction is pointing counter-clockwise when tokamak is viewed from above. This is true for almost all coordinate systems used in the field of plasma physics
    psi_axis, psi_bdry = psi_limits(M)
    verbose && println("ok!")

    verbose && println("---> Computing normalized poloidal flux coordinate (rho_pol) for (R_of_interest, z_of_interest) input variables... ")
    if length(W_abscissas_units[1])==3
        R_of_interests = vcat(units_conversion_factor(R_of_interest_unit,"m") .* R_of_interest) # Converting to meters to match Equilibrium.jl package standard
        z_of_interests = vcat(units_conversion_factor(z_of_interest_unit,"m") .* z_of_interest) # Converting to meters to match Equilibrium.jl package standard
        rho_of_interests = vcat(sqrt((M(R_of_interest,z_of_interest)-psi_axis)/(psi_bdry-psi_axis)))
    else
        # CODE FUTURE CASES HERE (3D, 4D ETC)
    end

    verbose && println("---> Preparing thermal ion data... ")
    thermal_ion_temperatures = zeros(length(regularization_thermal_ion_species),length(rho_of_interests)) # Ready for future versions with many rho_of_interest values
    thermal_ion_densities = zeros(length(regularization_thermal_ion_species),length(rho_of_interests)) # Ready for future versions with many rho_of_interest values
    
    for (i,ion) in regularization_thermal_ion_species
        verbose && println("------> Preparing thermal ion temperature data for species $(ion)... ")
        if (typeof(regularization_thermal_ion_temp[i]) <: Real)
            verbose && println("---------> Found an ion temperature of $(regularization_thermal_ion_temp[i]) keV. Incorporating... ")
            thermal_ion_temperatures[i,:] = repeat(vcat(regularization_thermal_ion_temp[i]),length(rho_of_interests))
        elseif (typeof(regularization_thermal_ion_temp[i]) <: String)
            if !isfile(regularization_thermal_ion_temp[i])
                error("regularization_thermal_ion_temp element with index $(i) was not specified correctly. $(regularization_thermal_ion_temp[i]) is not a valid filepath. Please correct and re-try.")
            end
            file_ext = split(regularization_thermal_ion_temp[i],".")[end] # Assume last part of String after "." is the file extension
            if lowercase(file_ext)=="jld2"
                myfile = jldopen(regularization_thermal_ion_temp[i],false,false,false,IOStream)
                rho_array = myfile["rho_pol"]
                temp_array = myfile["thermal_temp"]
                close(myfile)
                itp = Interpolations.interpolate((rho_array,),temp_array,Gridded(Linear()))
                etp = Interpolations.extrapolate(itp,0.0) # Assume vacuum scrape-off layer
                T_is = etp.(rho_of_interests) # T_i values at the rho_pol values of interest
                verbose && println("---------> Computed ion temperature(s) of $(T_is) keV")
                verbose && println("---------> at ρ_{pol}=$(round.(rho_of_interests,sigdigits=4))")
                verbose && println("---------> i.e. R = $(round.(R_of_interests,sigdigits=4)) meters")
                verbose && println("---------> i.e. z = $(round.(z_of_interests,sigdigits=4)) meters")
                verbose && println("---------> by loading $(regularization_thermal_ion_temp[i]) . Incorporating temperature(s)... ")
                thermal_ion_temperatures[i,:] = T_is
            elseif lowercase(file_ext)=="cdf"
                if (typeof(regularization_timepoint) <: Real)
                    t_mid = Float64.(regularization_timepoint)
                elseif (typeof(regularization_timepoint) <: String)
                    if !isfile(regularization_timepoint)
                        t_mid = regularization_timepoint # Assume format "XX.XXXX" or "XX,XXXX"
                    else
                        t_mid = (read_ncdf(regularization_timepoint,wanted_keys=["TIME"]))["TIME"] # Assume TRANSP-NUBEAM output file
                    end
                else
                    error("regularization_timepoint input variable has invalid type. Expected Float64, Int64 or String. regularization_timepoint input variable needed to load ion temperature data from $(regularization_thermal_ion_temp[i]). Please correct and re-try.")
                end
                etp = getTempProfileFromTRANSP(t_mid,regularization_thermal_ion_temp[i],regularization_thermal_ion_species[i])
                T_is = etp.(rho_of_interests)
                verbose && println("---------> Computed ion temperature(s) of $(T_is) keV")
                verbose && println("---------> at ρ_{pol}=$(round.(rho_of_interests,sigdigits=4))")
                verbose && println("---------> i.e. R = $(round.(R_of_interests,sigdigits=4)) meters")
                verbose && println("---------> i.e. z = $(round.(z_of_interests,sigdigits=4)) meters")
                verbose && println("---------> by loading $(regularization_thermal_ion_temp[i]) . Incorporating temperature(s)... ")
                thermal_ion_temperatures[i,:] = T_is
            else
                error("regularization_thermal_ion_temp element with index $(i) was not specified correctly. $(regularization_thermal_ion_temp[i]) has an invalid file extension. Expected .jld2 or .cdf. Please correct and re-try.")
            end
        else
            error("regularization_thermal_ion_temp element with index $(i) was not specified correctly. It is an invalid type $(regularization_thermal_ion_temp[i]). Please correct and re-try.")
        end

        verbose && println("------> Preparing thermal ion density data for species $(ion)... ")
        if (typeof(regularization_thermal_ion_dens[i]) <: Real)
            verbose && println("---------> Found an ion density of $(regularization_thermal_ion_dens[i]) m^-3. Incorporating... ")
            thermal_ion_densities[i,:] = repeat(vcat(regularization_thermal_ion_dens[i]),length(rho_of_interests))
        elseif (typeof(regularization_thermal_ion_dens[i]) <: String)
            if !isfile(regularization_thermal_ion_dens[i])
                error("regularization_thermal_ion_dens element with index $(i) was not specified correctly. $(regularization_thermal_ion_dens[i]) is not a valid filepath. Please correct and re-try.")
            end
            file_ext = split(regularization_thermal_ion_dens[i],".")[end] # Assume last part of String after "." is the file extension
            if lowercase(file_ext)=="jld2"
                myfile = jldopen(regularization_thermal_ion_dens[i],false,false,false,IOStream)
                rho_array = myfile["rho_pol"]
                dens_array = myfile["thermal_dens"]
                close(myfile)
                itp = Interpolations.interpolate((rho_array,),dens_array,Gridded(Linear()))
                etp = Interpolations.extrapolate(itp,0.0) # Assume vacuum scrape-off layer
                n_is = etp.(rho_of_interests) # n_i values at the rho_pol values of interest
                verbose && println("---------> Computed ion densitie(s) of $(n_is) m^-3")
                verbose && println("---------> at ρ_{pol}=$(round.(rho_of_interests,sigdigits=4))")
                verbose && println("---------> i.e. R = $(round.(R_of_interests,sigdigits=4)) meters")
                verbose && println("---------> i.e. z = $(round.(z_of_interests,sigdigits=4)) meters")
                verbose && println("---------> by loading $(regularization_thermal_ion_dens[i]) . Incorporating densitie(s)... ")
                thermal_ion_densities[i,:] = n_is
            elseif lowercase(file_ext)=="cdf"
                if (typeof(regularization_timepoint) <: Real)
                    t_mid = Float64.(regularization_timepoint)
                elseif (typeof(regularization_timepoint) <: String)
                    if !isfile(regularization_timepoint)
                        t_mid = regularization_timepoint # Assume format "XX.XXXX" or "XX,XXXX"
                    else
                        t_mid = (read_ncdf(regularization_timepoint,wanted_keys=["TIME"]))["TIME"] # Assume TRANSP-NUBEAM output file
                    end
                else
                    error("regularization_timepoint input variable has invalid type. Expected Float64, Int64 or String. regularization_timepoint input variable needed to load ion density data from $(regularization_thermal_ion_dens[i]). Please correct and re-try.")
                end
                etp = getDensProfileFromTRANSP(t_mid,regularization_thermal_ion_dens[i],regularization_thermal_ion_species[i])
                n_is = etp.(rho_of_interests)
                verbose && println("---------> Computed ion densities(s) of $(n_is) m^-3")
                verbose && println("---------> at ρ_{pol}=$(round.(rho_of_interests,sigdigits=4))")
                verbose && println("---------> i.e. R = $(round.(R_of_interests,sigdigits=4)) meters")
                verbose && println("---------> i.e. z = $(round.(z_of_interests,sigdigits=4)) meters")
                verbose && println("---------> by loading $(regularization_thermal_ion_dens[i]) . Incorporating densitie(s)... ")
                thermal_ion_densities[i,:] = n_is
            else
                error("regularization_thermal_ion_dens element with index $(i) was not specified correctly. $(regularization_thermal_ion_dens[i]) has an invalid file extension. Expected .jld2 or .cdf. Please correct and re-try.")
            end
        else
            error("regularization_thermal_ion_dens element with index $(i) was not specified correctly. It is an invalid type $(regularization_thermal_ion_dens[i]). Please correct and re-try.")
        end
    end
    # CODE ELECTRON TEMPERATURE AND DENSITY
    verbose && println("---> Preparing thermal electron data... ")
    thermal_electron_temperatures = zeros(length(rho_of_interests)) # Ready for future versions with many rho_of_interest values
    thermal_electron_densities = zeros(length(rho_of_interests)) # Ready for future versions with many rho_of_interest values
    
    verbose && println("------> Preparing thermal electron temperature data... ")
    if (typeof(regularization_thermal_electron_temp) <: Real)
        verbose && println("---------> Found an electron temperature of $(regularization_thermal_electron_temp) keV. Incorporating... ")
        thermal_electron_temperatures = repeat(vcat(regularization_thermal_electron_temp),length(rho_of_interests))
    elseif (typeof(regularization_thermal_electron_temp) <: String)
        if !isfile(regularization_thermal_electron_temp)
            error("regularization_thermal_electron_temp was not specified correctly. $(regularization_thermal_electron_temp) is not a valid filepath. Please correct and re-try.")
        end
        file_ext = split(regularization_thermal_electron_temp,".")[end] # Assume last part of String after "." is the file extension
        if lowercase(file_ext)=="jld2"
            myfile = jldopen(regularization_thermal_electron_temp,false,false,false,IOStream)
            rho_array = myfile["rho_pol"]
            temp_array = myfile["thermal_temp"]
            close(myfile)
            itp = Interpolations.interpolate((rho_array,),temp_array,Gridded(Linear()))
            etp = Interpolations.extrapolate(itp,0.0) # Assume vacuum scrape-off layer
            T_es = etp.(rho_of_interests) # T_i values at the rho_pol values of interest
            verbose && println("---------> Computed electron temperature(s) of $(T_es) keV")
            verbose && println("---------> at ρ_{pol}=$(round.(rho_of_interests,sigdigits=4))")
            verbose && println("---------> i.e. R = $(round.(R_of_interests,sigdigits=4)) meters")
            verbose && println("---------> i.e. z = $(round.(z_of_interests,sigdigits=4)) meters")
            verbose && println("---------> by loading $(regularization_thermal_electron_temp) . Incorporating temperature(s)... ")
            thermal_electron_temperatures = T_es
        elseif lowercase(file_ext)=="cdf"
            if (typeof(regularization_timepoint) <: Real)
                t_mid = Float64.(regularization_timepoint)
            elseif (typeof(regularization_timepoint) <: String)
                if !isfile(regularization_timepoint)
                    t_mid = regularization_timepoint # Assume format "XX.XXXX" or "XX,XXXX"
                else
                    t_mid = (read_ncdf(regularization_timepoint,wanted_keys=["TIME"]))["TIME"] # Assume TRANSP-NUBEAM output file
                end
            else
                error("regularization_timepoint input variable has invalid type. Expected Float64, Int64 or String. regularization_timepoint input variable needed to load electron temperature data from $(regularization_thermal_electron_temp). Please correct and re-try.")
            end
            etp = getTempProfileFromTRANSP(t_mid,regularization_thermal_electron_temp,"e")
            T_es = etp.(rho_of_interests)
            verbose && println("---------> Computed electron temperature(s) of $(T_es) keV")
            verbose && println("---------> at ρ_{pol}=$(round.(rho_of_interests,sigdigits=4))")
            verbose && println("---------> i.e. R = $(round.(R_of_interests,sigdigits=4)) meters")
            verbose && println("---------> i.e. z = $(round.(z_of_interests,sigdigits=4)) meters")
            verbose && println("---------> by loading $(regularization_thermal_electron_temp) . Incorporating temperature(s)... ")
            thermal_electron_temperatures = T_es
        else
            error("regularization_thermal_electron_temp was not specified correctly. $(regularization_thermal_electron_temp) has an invalid file extension. Expected .jld2 or .cdf. Please correct and re-try.")
        end
    else
        error("regularization_thermal_electron_temp was not specified correctly. It is an invalid type $(regularization_thermal_electron_temp). Please correct and re-try.")
    end
    
    verbose && println("------> Preparing thermal electron density data... ")
    if (typeof(regularization_thermal_electron_dens) <: Real)
        verbose && println("---------> Found an electron density of $(regularization_thermal_electron_dens) m^-3. Incorporating... ")
        thermal_electron_densities = repeat(vcat(regularization_thermal_electron_dens),length(rho_of_interests))
    elseif (typeof(regularization_thermal_electron_dens) <: String)
        if !isfile(regularization_thermal_electron_dens)
            error("regularization_thermal_electron_dens was not specified correctly. $(regularization_thermal_electron_dens) is not a valid filepath. Please correct and re-try.")
        end
        file_ext = split(regularization_thermal_electron_dens,".")[end] # Assume last part of String after "." is the file extension
        if lowercase(file_ext)=="jld2"
            myfile = jldopen(regularization_thermal_electron_dens,false,false,false,IOStream)
            rho_array = myfile["rho_pol"]
            dens_array = myfile["thermal_dens"]
            close(myfile)
            itp = Interpolations.interpolate((rho_array,),dens_array,Gridded(Linear()))
            etp = Interpolations.extrapolate(itp,0.0) # Assume vacuum scrape-off layer
            n_es = etp.(rho_of_interests) # n_e values at the rho_pol values of interest
            verbose && println("---------> Computed electron densitie(s) of $(n_es) m^-3")
            verbose && println("---------> at ρ_{pol}=$(round.(rho_of_interests,sigdigits=4))")
            verbose && println("---------> i.e. R = $(round.(R_of_interests,sigdigits=4)) meters")
            verbose && println("---------> i.e. z = $(round.(z_of_interests,sigdigits=4)) meters")
            verbose && println("---------> by loading $(regularization_thermal_electron_dens) . Incorporating densitie(s)... ")
            thermal_electron_densities = n_es
        elseif lowercase(file_ext)=="cdf"
            if (typeof(regularization_timepoint) <: Real)
                t_mid = Float64.(regularization_timepoint)
            elseif (typeof(regularization_timepoint) <: String)
                if !isfile(regularization_timepoint)
                    t_mid = regularization_timepoint # Assume format "XX.XXXX" or "XX,XXXX"
                else
                    t_mid = (read_ncdf(regularization_timepoint,wanted_keys=["TIME"]))["TIME"] # Assume TRANSP-NUBEAM output file
                end
            else
                error("regularization_timepoint input variable has invalid type. Expected Float64, Int64 or String. regularization_timepoint input variable needed to load electron density data from $(regularization_thermal_electron_dens). Please correct and re-try.")
            end
            etp = getDensProfileFromTRANSP(t_mid,regularization_thermal_electron_dens,"e")
            n_es = etp.(rho_of_interests)
            verbose && println("---------> Computed electron densities(s) of $(n_es) m^-3")
            verbose && println("---------> at ρ_{pol}=$(round.(rho_of_interests,sigdigits=4))")
            verbose && println("---------> i.e. R = $(round.(R_of_interests,sigdigits=4)) meters")
            verbose && println("---------> i.e. z = $(round.(z_of_interests,sigdigits=4)) meters")
            verbose && println("---------> by loading $(regularization_thermal_electron_dens) . Incorporating densitie(s)... ")
            thermal_electron_densities = n_es
        else
            error("regularization_thermal_electron_dens was not specified correctly. $(regularization_thermal_electron_dens) has an invalid file extension. Expected .jld2 or .cdf. Please correct and re-try.")
        end
    else
        error("regularization_thermal_electron_dens was not specified correctly. It is an invalid type $(regularization_thermal_electron_dens). Please correct and re-try.")
    end

    verbose && print("---> Computing slowing-down (collision physics) basis functions...")
    if is_energy_pitch(W_abscissas[1],W_abscissas_units[1])
        n_e = thermal_electron_densities[1] # We know it must be a 1-element Vector because reconstruction space is (E,p)
        T_e = thermal_electron_temperatures[1] # -||-
        thermal_ion_densities_at_rho = thermal_ion_densities[:,1] # -||-
        thermal_ion_temperatures_at_rho = thermal_ion_temperatures[:,1] # -||-

        verbose && println("for (E,p) space with grid point extrema:")
        w_E_abscissa = get_energy_abscissa(W_abscissas[1],W_abscissas_units[1])
        verbose && println("------> (E_min,E_max)=$(round.(extrema(w_E_abscissa),sigdigits=4))")
        w_p_abscissa = get_pitch_abscissa(W_abscissas[1],W_abscissas_units[1])
        verbose && println("------> (p_min,p_max)=$(round.(extrema(w_p_abscissa),sigdigits=4))")
        dE = diff(w_E_abscissa)[1]; dp = diff(w_p_abscissa)[1] # Assume equidistant (E,p) grid

        E_inj_array = w_E_abscissa[1:1:end-1] .+ diff(w_E_abscissa)[1]/2 # Assume uniform energy grid
        p_inj_array = w_p_abscissa[1:1:end-1] .+ diff(w_p_abscissa)[1]/2 # Assume uniform pitch grid

        verbose && println("------> Computing... ")
        F_SD = zeros(length(w_E_abscissa)*length(w_p_abscissa),length(E_inj_array)*length(p_inj_array))
        F_SD_coords = CartesianIndices((length(E_inj_array),length(p_inj_array)))
        for (ic,c) in enumerate(F_SD_coords)
            iE = c[1]; ip = c[2]; E_inj = E_inj_array[iE]; p_inj = p_inj_array[ip]
            f_SD = slowing_down_function(E_inj, p_inj, w_E_abscissa, w_p_abscissa, n_e, T_e, FI_species, regularization_thermal_ion_species, thermal_ion_densities_at_rho, thermal_ion_temperatures_at_rho; dampen=true, damp_type=:erfc, sigma=2.0)
            f_SD = f_SD ./sum((dE*dp) .*f_SD) # Normalize the basis function so they integrate to 1.0
            F_SD[:,ic] .= reshape(f_SD,(length(w_E_abscissa)*length(w_p_abscissa),1))
        end
    elseif is_vpara_vperp(W_abscissas[1],W_abscissas_units[1])
        v2E_rel = (v-> (FI_species*(GuidingCenterOrbits.c0)^2)*(sqrt(1/(1-(v/(GuidingCenterOrbits.c0))^(2)))-1)) # A one-line function to transform from relativistic speed (m/s) to energy (keV)

        n_e = thermal_electron_densities[1] # We know it must be a 1-element Vector because reconstruction space is (E,p)
        T_e = thermal_electron_temperatures[1] # -||-
        thermal_ion_densities_at_rho = thermal_ion_densities[:,1] # -||-
        thermal_ion_temperatures_at_rho = thermal_ion_temperatures[:,1] # -||-

        verbose && println("for (vpara,vperp) space with grid point extrema:")
        w_vpara_abscissa = get_vpara_abscissa(W_abscissas[1],W_abscissas_units[1]); nwvpa = length(w_vpara_abscissa)
        verbose && println("------> (vpara_min,vpara_max)=$(round.(extrema(w_vpara_abscissa),sigdigits=4))")
        w_vperp_abscissa = get_vperp_abscissa(W_abscissas[1],W_abscissas_units[1]); nwvpe = length(w_vperp_abscissa)
        verbose && println("------> (vperp_min,vperp_max)=$(round.(extrema(w_vperp_abscissa),sigdigits=4))")
        dvpa = diff(w_vpara_abscissa)[1]; dvpe = diff(w_vperp_abscissa)[1] # Assume equidistant (vpara,vperp) grid

        # To find minimum and maximum of (E,p) grid for
        VEL_points = Iterators.product(w_vpara_abscissa, w_vperp_abscissa)
        min_E, max_E = Inf, -Inf
        min_p, max_p = Inf, -Inf
        for VEL_point in VEL_points
            v = sqrt(VEL_point[1]^2 + VEL_point[2]^2)
            p = VEL_point[1]/v
            E = v2E_rel(v)
            if E<min_E
                min_E = E
            end
            if E>max_E
                max_E = E
            end
            if p<min_p
                min_p = p
            end
            if p>max_p
                max_p = p
            end
        end
        E_SD_array = collect(range(min_E,stop=max_E,length=sqrt(nwvpa*nwvpe)))
        verbose && println("------> Inferred (E_min,E_max)=$(round.(extrema(E_SD_array),sigdigits=4))")
        p_SD_array = collect(range(min_p,stop=max_p,length=sqrt(nwvpa*nwvpe)))
        dE = diff(E_SD_array)[1]; dp = diff(p_SD_array)[1]
        verbose && println("------> Inferred (p_min,p_max)=$(round.(extrema(p_SD_array),sigdigits=4))")

        E_inj_array = E_SD_array[1:1:end-1] .+ diff(E_SD_array)[1]/2
        p_inj_array = p_SD_array[1:1:end-1] .+ diff(p_SD_array)[1]/2

        verbose && println("------> Computing... ")
        F_SD = zeros(length(w_vpara_abscissa)*length(w_vperp_abscissa),length(E_inj_array)*length(p_inj_array))
        F_SD_coords = CartesianIndices((length(E_inj_array),length(p_inj_array)))
        for (ic,c) in enumerate(F_SD_coords)
            iE = c[1]; ip = c[2]; E_inj = E_inj_array[iE]; p_inj = p_inj_array[ip]
            f_SD_Ep = slowing_down_function(E_inj, p_inj, E_SD_array, p_SD_array, n_e, T_e, FI_species, regularization_thermal_ion_species, thermal_ion_densities_at_rho, thermal_ion_temperatures_at_rho; dampen=true, damp_type=:erfc, sigma=2.0)
            f_SD_VEL, vpara_SD_array, vperp_SD_array = Ep2VparaVperp(E_SD_array, p_SD_array, f_SD_Ep; my_gcp=getGCP(FI_species)(0.0,0.0,0.0,0.0), needJac=true, returnAbscissas=true)
            nodes = (vpara_SD_array, vperp_SD_array)
            itp = Interpolations.interpolate(nodes,f_SD_VEL,Gridded(Linear()))
            etp = Interpolations.extrapolate(itp,Flat())
            f_SD = zeros(length(w_vpara_abscissa),length(w_vperp_abscissa))
            for (ivpa,vpa) in enumerate(w_vpara_abscissa)
                for (ivpe,vpe) in enumerate(w_vperp_abscissa)
                    f_SD[ivpa,ivpe] = etp(vpa,vpe)
                end
            end
            f_SD = f_SD ./sum((dvpa*dvpe) .*f_SD) # Normalize the basis function so they integrate to 1.0
            F_SD[:,ic] .= reshape(f_SD,(length(w_vpara_abscissa)*length(w_vperp_abscissa),1))
        end
    elseif false # <----- ADD MORE OPTIONS HERE IN FUTURE VERSIONS!!!
    else
        error("This error should be impossible to reach. How did you do that?")
    end
    F_SD_safe = reduce(hcat,filter(col-> sum(col)!=0.0 && sum(isnan.(col))==0, eachcol(F_SD))) # Filter out all basis functions that are NOT identically zero and do NOT contain NaNs, and hcat them back into a matrix
end

###########################################################################################################
# SECTION 10: ENFORCE THE NOISE FLOOR

verbose && println("Enforcing a noise floor of $(noise_floor_factor). All error values smaller than $(noise_floor_factor)*maximum(S) will be lifted up to $(noise_floor_factor)*maximum(S)... ")
for (i,s) in enumerate(S)
    noise_floor_i = maximum(s)*noise_floor_factor
    verbose && println("---> Diagnostic $(i) has a noise floor of $(noise_floor_i)... ")
    below_noise_floor_inds = findall(x-> x<noise_floor_i, S_errs[i])
    verbose && println("---> Found $(length(below_noise_floor_inds)) error values smaller than $(noise_floor_i). Lifting them up... ")
    S_errs[below_noise_floor_inds] .= noise_floor_i
end

###########################################################################################################
# SECTION 11: EXCLUDE UNWANTED MEASUREMENT BINS

if !(length(excluded_measurement_intervals)==length(excluded_measurement_units)==length(filepaths_S))
    @warn "The lengths of input variables 'excluded_measurement_intervals' and 'excluded_measurement_units' were not equal to 'filepaths_S'. Cannot execute exclusion of unwanted measurement bins. All measurement bins will be included for all diagnostics."
else
    verbose && println("Performing exclusion of unwanted measurement bins... ")
    for i in eachindex(excluded_measurement_units)
        verbose && println("---> For diagnostic $(i)... ")
        excluded_measurement_intervals_i = excluded_measurement_intervals[i]
        excluded_measurement_units_i = excluded_measurement_units[i]
        W_i = W[i]
        S_i = S[i]
        S_errs_i = S_errs[i]
        S_abscissas_units_i = S_abscissas_units[i]
        S_abscissas_i = S_abscissas[i]
        ucf_i = units_conversion_factor(excluded_measurement_units_i,S_abscissas_units_i)
        bad_inds_i = []
        for interval in excluded_measurement_intervals_i
            interval_min = ucf_i*interval[1]
            interval_max = ucf_i*interval[2]
            verbose && println("------> Identified excluded measurement region (min,max)=($(interval_min),$(interval_max)) $(S_abscissas_units_i). Excluding... ") 
            append!(bad_inds, findall(xx-> (xx>=interval_min && xx<=interval_max),S_abscissas_i))
        end
        if exclude_zero_measurements
            zero_inds = findall(x-> x==0.0, S_i)
            verbose && println("---> For diagnostic $(i), found $(length(zero_inds)) measurements equal to 0. Excluding... ")
            append!(bad_inds, zero_inds)
        end
        bad_inds_i = unique(bad_inds_i) # The indices of the excluded measurement intervals
        good_inds_i = filter(xx -> !(xx in bad_inds_i), collect(1:length(S_abscissas_i))) # The indices of all measurement bins to keep
        verbose && println("---> A total of $(length(bad_inds_i)) measurement bin centers will be excluded for diagnostic $(i)...")
        S[i] = S_i[good_inds_i] # Update measurements
        S_errs[i] = S_errs_i[good_inds_i] # Update errors/uncertainties
        S_abscissas[i] = S_abscissas_i[good_inds_i] # Update measurement abscissa
        W[i] = W_i[good_inds_i,:] # Update weight function matrix
    end
end

###########################################################################################################
# SECTION 12: 
# CREATE 1st ORDER TIKHONOV MATRIX BEFORE EXCLUDING BAD RECONSTRUCTION SPACE POINTS
# CREATE 1st ORDER TIKHONOV MATRIX BEFORE EXCLUDING BAD RECONSTRUCTION SPACE POINTS
# CREATE 1st ORDER TIKHONOV MATRIX BEFORE EXCLUDING BAD RECONSTRUCTION SPACE POINTS

###########################################################################################################
# SECTION 13: EXCLUDE ALL RECONSTRUCTION SPACE POINTS FOR WHICH THE WEIGHT MATRICES OF ALL DIAGNOSTICS ARE ZERO
# PLEASE NOTE! THIS ALSO TAKES CARE OF INVALID ORBITS, SHOULD THE RECONSTRUCTION SPACE BE E.G. ORBIT SPACE 
# OR CONSTANTS-OF-MOTION SPACE.

verbose && println("Finding reconstruction space points for which the weight matrix (column) is zero (no sensitivity), for all diagnostics... ")
rec_space_coords = CartesianIndices(tuple(length.(W_abscissas[1][2:end]))) # All reconstruction space coordinates. All w abscissas are the same (because of sections 4 and 5), hence [1]. First w abscissa is the diagnostic measurement bins, hence [2:end].
rec_space_zero_dict = Dict() # A Dict of CartesianIndex to keep track of all reconstruction space coordinates for which the weight matrix is zero (the whole column is zero), and for how many diagnostics
for w in W # For all weight matrices (all diagnostics)
    for (j,c) in enumerate(eachcol(w)) # For each column in the weight matrix (and the column index j)
        if sum(c)==0 # If the sum of the column is zero
            if haskey(rec_space_zero_dict,rec_space_coords[j]) # Increase the corresponding reconstruction space coordinate dictionary element by 1
                rec_space_zero_dict[rec_space_coords[j]] += 1
            else # To avoid errors, separate the cases when the dictionary key is created, and when it's not
                rec_space_zero_dict[rec_space_coords[j]] = 1
            end
        end
    end
end
bad_inds = []
if !isempty(rec_space_zero_dict)
    for key in keys(rec_space_zero_dict)
        if rec_space_zero_dict[key]==length(W) # If the reconstruction space point has zero sensitivity, for all diagnostics
            verbose && println("---> The reconstruction space point (bad point) $([W_abscissas[1][i+1][key[i]] for i in 1:length(key)]) with units $([W_abscissas_units[1][i+1] for i in 1:length(key)]) has no sensitivity for all diagnostics. Excluding from reconstruction... ")
            append!(bad_inds, findfirst(x-> x==key, rec_space_coords))
        end
    end
    bad_inds = unique(bad_inds) # Should be redundant, but just in case
    good_inds = filter(x-> !(x in bad_inds), collect(1:length(rec_space_coords)))
    
    rec_space_coords = rec_space_coords[good_inds] # Update the list of reconstruction space coordinates, to be able to correctly inflate tomographic reconstruction vector to full N-dimensional array
    for (i,w) in enumerate(W)
        verbose && println("------> Excluding bad points from the weight matrix of diagnostic $(i)... ")
        W[i] = w[:,good_inds]
    end
    if "collisions" in lowercase.(String.(regularization)) # Exclude bad points from slowing-down (collision physics) basis functions
        verbose && println("------> Excluding bad points from the slowing-down (collision physics) basis functions... ")
        F_SD_safe = F_SD_safe[good_inds,:]
    end

    if "firsttikhonov" in lowercase.(String.(regularization)) # Exclude bad points from 1st order Tikhonov matrix

    end
else
    verbose && println("---> No reconstruction space points found.")
end

###########################################################################################################
# SECTION 14: INCLUDE COLLISIONAL PHYSICS BY MULTIPLYING THE WEIGHT MATRICES WITH THE SLOWING-DOWN BASIS FUNCTIONS

if "collisions" in lowercase.(String.(regularization)) # Don't need all the safety checks here, it was already taken care of in section 9
    verbose && println("Multiplying weight matrices with slowing-down (collision physics) matrix to include collision physics as prior information/regularization... ")
    for (i,w) in enumerate(W)
        verbose && print("---> For diagnostic $(i)... ")
        W[i] = w*F_SD_safe
        verbose && println("ok!")
    end
end

###########################################################################################################
# SECTION 15: RESCALE WEIGHT FUNCTIONS

if rescale_W
    verbose && println("Rescaling weight matrices using computed rescale factors... ")
    for (i,r) in enumerate(rescale_W_factors)
        verbose && print("---> Rescaling weight matrix for diagnostic $(i) with a factor of $(round(r,sigdigits=4))... ")
        W[i] = r .*W[i]
        verbose && println("ok!")
    end
end

###########################################################################################################
# SECTION 16: CONCATENATE ALL MEASUREMENTS AND BUILD TOTAL WEIGHT MATRIX

# CONTINUE CODING HERE!!!
# CONTINUE CODING HERE!!!
# CONTINUE CODING HERE!!!