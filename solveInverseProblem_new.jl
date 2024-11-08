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
if :COLLISIONS in regularization
    include(folderpath_OWCF*"extra/dependencies.jl")
    include(folderpath_OWCF*"misc/temp_n_dens.jl")
end
if lowercase(String(rescale_W_F_ref_source))=="file"
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
S = Vector{Vector{Float64}}(undef,length(filepaths_S)) # Pre-allocate measurements for all diagnostics
S_units = Vector{String}(undef,length(filepaths_S)) # Pre-allocate units of measurements for all diagnostics
S_errs = Vector{Vector{Float64}}(undef,length(filepaths_S)) # Pre-allocate errors (uncertainties) for all diagnostics
S_errs_units = Vector{String}(undef,length(filepaths_S)) # Pre-allocate units of errors (uncertainties) for all diagnostics
S_abscissas = Vector{Vector{Float64}}(undef,length(filepaths_S)) # Pre-allocate measurement bin centers for all diagnostics
S_abscissas_units = Vector{String}(undef,length(filepaths_S)) # Pre-allocate units of measurement bin centers for all diagnostics
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

W_inflated = Vector{Array{Float64,DIM}}(undef,length(filepaths_W)) # The (inflated) weight functions, for each diagnostic. PLEASE NOTE! All weight functions need to have the same number of dimensions!
W_abscissas = Vector{Vector{Vector{Real}}}(undef,length(filepaths_W)) # The abscissas (measurement bin centers + phase-space grid points), for each set of weight functions
W_abscissas_units = Vector{Vector{String}}(undef,length(filepaths_W)) # The units of measurement, for each abscissa, for each set of weight functions
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
# SECTION 4: CHECK THAT ALL ABSCISSAS MATCH IN TERMS OF DIMENSIONS AND UNITS

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
# SECTION 6: WRITE SOMETHING TO TAKE CARE OF ONLY INCLUDING VALID ORBITS 
# IF RECONSTRUCTING IN ORBIT SPACE

###########################################################################################################
# SECTION 7: RESHAPE ALL INFLATED WEIGHT MATRICES INTO THEIR 2D SHAPE, TO BE USED IN INVERSE PROBLEMS

W = Vector{Array{Float64,2}}(undef,length(filepaths_W)) # The weight matrices, for each diagnostic
for (i,w_inflated) in enumerate(W_inflated)
    ws = size(w_inflated)
    W[i] = reshape(w_inflated,(ws[1],reduce(*,ws[2:end])))
end
W_inflated = nothing # Clear memory. To minimize memory usage

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
        if length(w_abscissas_units)==3 # Deduce (E,p) or (vpara,vperp) from w_abscissas_units (length of 3 means reconstruction space dim of 2)
            verbose && print("------> Weight functions reconstruction space is 2D. Deducing coordinates type... ")
            units_1 = w_abscissas_units[2] # The second dimension of the weight matrix is the first dimension of the reconstruction space
            units_2 = w_abscissas_units[3] # The third dimension of the weight matrix is the second dimension of the reconstruction space

            units_tot = vcat(units_1, units_2)
            w_energy_ind = findall(x-> x in ENERGY_UNITS || x in ENERGY_UNITS_LONG, units_tot)
            w_pitch_ind = findall(x-> x in DIMENSIONLESS_UNITS || x in DIMENSIONLESS_UNITS_LONG, units_tot)
            w_speed_inds = findall(x-> units_are_speed(x), units_tot)
            # IN THE FUTURE, ADD MORE UNIT CHECKS HERE IF NEEDED
            if length(w_energy_ind)==1 && length(w_pitch_ind)==1
                w_energy_ind .+= 1 # To align it with the original w_abscissas and w_abscissas_units indices
                w_pitch_ind .+= 1 # To align it with the original w_abscissas and w_abscissas_units indices
                w_rec_space = "(E,p)"
                verbose && println(w_rec_space)
            elseif length(w_speed_inds)==2
                w_speed_inds .+= 1 # To align it with the original w_abscissas and w_abscissas_units indices
                w_rec_space = "(vpara,vperp)"
                verbose && print(w_rec_space*". Distinguishing (vpara, vperp) arrays...")
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
            elseif false # ADD MORE 2D CASES HERE IN FUTURE UPDATES
            else
                error(standard_rescale_W_error)
            end
            R_of_interests = vcat(units_conversion_factor(R_of_interest_unit,"m")*R_of_interest) # R_of_interests will be the same, regardless of 2D coordinates
            z_of_interests = vcat(units_conversion_factor(z_of_interest_unit,"m")*z_of_interest) # z_of_interests will be the same, regardless of 2D coordinates
        elseif length(w_abscissas_units)==5 # DEDUCE (E,p,R,z) FROM W_abscissas_units
            verbose && print("------> Weight functions reconstruction space is 4D. Deducing coordinates type... ")
            units_1 = w_abscissas_units[2] # The second dimension of the weight matrix is the first dimension of the reconstruction space
            units_2 = w_abscissas_units[3] # The third dimension of the weight matrix is the second dimension of the reconstruction space
            units_3 = w_abscissas_units[4] # The fourth dimension of the weight matrix is the third dimension of the reconstruction space
            units_4 = w_abscissas_units[5] # The fifth dimension of the weight matrix is the fourth dimension of the reconstruction space

            units_tot = vcat(units_1, units_2, units_3, units_4)
            w_energy_ind = findall(x-> x in ENERGY_UNITS || x in ENERGY_UNITS_LONG, units_tot)
            w_pitch_ind = findall(x-> x in DIMENSIONLESS_UNITS || x in DIMENSIONLESS_UNITS_LONG, units_tot)
            w_Rz_inds = findall(x-> x in LENGTH_UNITS || x in LENGTH_UNITS_LONG, units_tot)
            # IN THE FUTURE, ADD MORE UNIT CHECKS HERE IF NEEDED
            if length(w_energy_ind)==1 && length(w_pitch_ind)==1 && length(w_Rz_inds)==2
                w_energy_ind .+= 1 # To align it with the original w_abscissas and w_abscissas_units indices
                w_pitch_ind .+= 1 # To align it with the original w_abscissas and w_abscissas_units indices
                w_Rz_inds .+= 1 # To align it with the original w_abscissas and w_abscissas_units indices
                w_rec_space = "(E,p,R,z)"
                verbose && print(w_rec_space*". Distinguishing (R,z) arrays... ")

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
                R_of_interests = units_conversion_factor(w_abscissas_units[w_R_ind],"m") .*R_of_interests # Need unit of measurement meter for CDFto4D function (if that will be used)
                z_of_interests = units_conversion_factor(w_abscissas_units[w_z_ind],"m") .*z_of_interests # Need unit of measurement meter for CDFto4D function (if that will be used)
            elseif false # ADD MORE 4D CASES HERE IN FUTURE UPDATES
            else
                error(standard_rescale_W_error)
            end
        elseif false # ADD MORE (1D-6D) CASES HERE IN FUTURE UPDATES. SUCH AS (E,pm, Rm), (E,mu,Pphi), (E,mu,Pphi;sigma) etc
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
            F_ref = units_conversion_factor("keV^-1",units_inverse(w_energy_unit)) .*F_ref # simultaneously, convert the fast-ion distribution via the inverse of the energy units

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
            # CONTINUE CODING HERE!!!
            # CONTINUE CODING HERE!!!
            # CONTINUE CODING HERE!!!

        elseif w_rec_space=="(E,p,R,z)"
            verbose && println("------> Using reference fast-ion distribution assuming weight functions are given in (E,p,R,z) coordinates (or similar)... ")
            # CONTINUE CODING HERE!!!
            # CONTINUE CODING HERE!!!
            # CONTINUE CODING HERE!!!
            
        elseif false
        else
            error(standard_rescale_W_error)
        end
    elseif lowercase(String(rescale_W_F_ref_source))=="gaussian"
        # COMPUTE A BUNCH OF RANDOM GAUSSIANS FOR THE RECONSTRUCTION SPACE
        # RESHAPE AND MULTIPLY WITH WEIGHT FUNCTIONS
        # ADD TO AVERAGE 
        # LET THE AVERAGE BE THE RESCALE FACTORS FOR EACH DIAGNOSTIC
    else
        error("'rescale_W' was set to true but 'rescale_W_F_ref_source' was not specified correctly. Currently supported options include :FILE and :GAUSSIAN. Please correct and re-try.")
    end
else
    verbose && println("Weight functions will NOT be rescaled... ")
    rescale_W_factors = ones(length(filepaths_W))
end
# Plan suggestion:
# Check if file or gaussians
# If file, check the scriptSources_W
# If any of the known OWCF scripts, try jld2to4d, h5to4d and cdfto4d
# If not, raise error (not supported yet?)
# If gaussian, proceed quite straight forward