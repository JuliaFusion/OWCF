#########################################  solveInverseProblem.jl #########################################

#### Description:
# This script solves an inverse problem on the form S=W*F where 'S' is a vector of length m containing all the 
# measurements/data points, 'F' is a vector of length n containing the fast-ion distribution in vectorized form 
# and 'W' is an m x n matrix containing all the weight functions (as rows) relating 'S' to 'F'. The script can 
# solve standard fast-ion tomography inverse problems in any dimension D, where 1 <= D <= 6, including orbit 
# tomography problems. Several forms of prior information, constraints and regularization can be included when 
# solving the inverse problem. These are described in the template file OWCF/templates/start_solveInverseProblem_template.jl.
# Any number of diagnostics can be included in the inverse problem, as long as they are included together with 
# an equal number of weight matrices.
#
# Please see the OWCF/templates/start_solveInverseProblem_template.jl file for further input information.

#### Inputs (Units given when defined in script)
# Given via input file start_solveInverseProblem_template.jl, for example. 'template' should be replaced by whatever.

#### Outputs
# -

#### Saved files
# The results output file will have the filename format:
#   solInvProb_[date and time]_[tokamak]_S_[numODiag]_[numOMeasure]_W_[numOW]_[recSpaceGridSize]_[hyperParamGridSize].jld2
# (or instead with a .hdf5 file extension if the 'saveInHDF5format' input variable was set to true in the start file)
# where 
#   - 'date and time' is the date and time at which the inverse problem was solved
#   - 'tokamak' is the tokamak in which the inverse problem is solved. Specified in input file
#   - 'numODiag' is the number of diagnostics used to solve the inverse problem
#   - 'numOMeasure' is the total number of measurements used to solve the inverse problem
#   - 'numOW' is the number of weight matrices used to solve the inverse problem (equal to numODiag)
#   - 'recSpaceGridSize' is the size of the reconstruction space grid on which the fast-ion distribution is reconstructed
#   - 'hyperParamGridSize' is the size of the hyper-parameter value(s) grid used to regularize the inverse problem
# This results output file will have the keys:
#   sols - The solutions to the inverse problem in inflated format, i.e. recSpaceGridSize x hyperParamGridSize - Array{Float64}
#   sols_abscissa_i - where i is an integer >=1. The reconstruction space abscissas of the sols output data. E.g. sols_abscissa_i
#                     contains the grid points corresponding to the i:th dimension of sols - Vector{Float64}
#   sols_abscissa_units_i - where i is an integer >=1. The units of measurement of sols_abscissa_i - String
#   sols_hyper_abscissa_j - where j is an integer >=1. The hyper-parameter values of the j:th regularization dimension of sols.
#                           For example, if the reconstruction space is 2D and two hyper-parameters where needed to solved the inverse 
#                           problem (e.g. the regularization strength for 0th order Tikhonov and 1st order Tikhonov, respectively) 
#                           sols[:,:,j1,j2] corresponds to the solution with hyper-parameter values sols_hyper_abscissa_1[j1] and 
#                           sols_hyper_abscissa_2[j2]. - Vector{Float64}
#   FI_species - The particle species of the reconstructed fast-ion distribution - String
#   S_k - where k is an integer >=1. The measurements of the k:th diagnostic used to solve the inverse problem - Vector{Float64}
#   S_units_k - where k is an integer >=1. The units of S_k - String 
#   S_err_k - where k is an integer >=1. The uncertainties of the S_k measurements. S_err_k[l] corresponds to the uncertainty 
#             of measurement S_k[l] - Vector{Float64}
#   S_err_units_k - where k is an integer >=1. The units of S_err_k - String
#   S_abscissa_k - where k is an integer >=1. The measurement bin centers of S_k. S_abscissa_k[l] corresponds to the measurement 
#                  bin center of S_k[l] - Vector{Float64}
#   S_abscissa_units_k - where k is an integer >=1. The units of S_abscissa_k - String
#   l_curve_x - The x-axis values of the L-curve of the solutions to the inverse problem, corresponding to ||S-WF|| - Vector{Float64}
#   l_curve_y - The y-axis-values of the L-curve of the solutions to the inverse problem, corresponding to ||F|| - Vector{Float64}
#   l_curve_opt_index - The index of the L-curve point with the largest curvature. (l_curve_x[l_curve_opt_index], l_curve_y[l_curve_opt_index])
#                       will then correpond to the L-curve point with the largest curvature - Int64
#   l_curve_opt_hyper_point - The hyper-parameter point corresponding to the L-curve point with the largest curvature. I.e. the values of 
#                             λ_h where h is an integer >=1 and λ_h is the regularization strength for regularization type h - Vector{Float64}
#   l_curve_opt_sol - The solution to the inverse problem obtained using the l_curve_opt_hyper_point regularization strength(s). l_curve_opt_sol 
#                     will be an array of size recSpaceGridSize. l_curve_opt_sol is equivalent to sols[:,...,:,j1,j2,...,jH] where :,...,: denotes 
#                     all reconstruction space dimensions and jh = findfirst(x-> x==l_curve_opt_hyper_point[h],sols_hyper_abscissa_h) and h>=1 - Array{Float64}
#   filepath_start - The filepath to the start file used to execute the solveInverseProblem.jl script - String
#   If rescale_W was set to true
#       rescale_W_factors - The factors used to rescale the weight functions to have WF_synth match the experimental data - Vector{Float64}
#   If :COLLISIONS and/or :ICRF was included in the 'regularization' input variable
#       regularization_equil_filepath - The file path to the magnetic equilibrium data used to compute :COLLISIONS and/or :ICRF prior - String
#   If :COLLISIONS was included in the 'regularization' input variable
#       coll_phys_basis - The collision physics basis functions used to solve the inverse problem. The basis functions are saved as a multi-dimensional 
#                         array with size recSpaceGridSize x reduce(*, recSpaceGridSize) - Array{Float64}
#       coll_phys_thermal_species - The list of thermal ion species used to compute the collision physics basis functions - Vector{String}
#       coll_phys_Ti - The list of thermal ion temperatures used to compute the collision physics basis functions - Vector{Float64}
#       coll_phys_ni - The list of thermal ion densities used to compute the collision physics basis functions - Vector{Float64}
#       coll_phys_Te - The electron temperature used to compute the collision physics basis functions - Float64
#       coll_phys_ne - The electron density used to compute the collision physics basis functions - Float64
#       coll_phys_rho - The normalized poloidal flux coordinate(s) used to compute the collision physics basis functions - Vector{Float64}
#   If :ICRF was included in the 'regularization' input variable
#       L_ICRF - The matrix used to regularize the inverse problem using the physics of electromagnetic wave heating in the ion 
#                cyclotron range of frequencies (ICRF) - Matrix{Float64}
# If the 'plot_solutions' input variable was set to true in the start file, .png files will be saved showing plots of the solutions.
# The filename format of those .png files will be:
#   solInvProb_[date and time]_resultsPlot_[p]_of_[P].png
# where
#   - 'date and time' is the date and time at which the inverse problem was solved
#   - 'p' is an integer 1<=p<=P 
#   - 'P' is an integer equal to the reduce(*,hyperParamGridSize) (the total number of hyper-parameter points for all regularization combinations)
# If the 'gif_solutions' input variable was set to true in the start file, a .gif file will be saved showing all .png files together as an animation.
# The filename format of that .gif file will be:
#   solInvProb_[date and time].gif
# where 
#   - 'date and time' is the date and time at which the inverse problem was solved
# ALL OUTPUT FILES WILL BES SAVED IN THE 'folderpath_out' INPUT VARIABLE FOLDER, SPECIFIED IN THE START FILE.

### Other
# 

# Script written by Henrik Järleblad. Last maintained 2025-09-01.
###########################################################################################################

# A dictionary to keep a record of the names of all sections and their respective execution times
dictionary_of_sections = Dict()

# println section number tracker
prnt = 0

# Timestamp of script execution start
timestamps = [time()]

# SECTION 0
println("------------------------------------------------ SECTION $(prnt) ------------------------------------------------"); prnt += 1

verbose && println("Loading Julia packages... ")
using Base.Iterators
using Dates
using FileIO
using HDF5
using Interpolations
using JLD2
using LinearAlgebra
using Plots # If macros (@animate) did not expand at parse time, this would only need to be loaded if plot_solutions || gif_solutions
using SparseArrays
using Statistics
using SCS, Convex
include("extra/constants.jl")
include("misc/convert_units.jl")

if "collisions" in lowercase.(String.(regularization))
    verbose && println("Collision physics included in list of regularization!")
    verbose && print("---> ")
    include("extra/dependencies.jl")
    include("misc/temp_n_dens.jl")
    include("misc/species_func.jl")
end
if ("icrf" in lowercase.(String.(regularization)))
    verbose && println("ICRF physics included in list of regularization!")
    verbose && print("---> ")
    include("extra/dependencies.jl")
end
if ("firsttikhonov" in lowercase.(String.(regularization)))
    verbose && println("1st order Tikhonov in list of regularization!")
    verbose && print("---> ")
    include("extra/dependencies.jl")
end
if rescale_W
    verbose && println("rescale_W set to true!")
    verbose && print("---> ")
    include("extra/dependencies.jl")
end

append!(timestamps,time()) # The timestamp when all necessary packages have been loaded
dictionary_of_sections[prnt-1] = ("Loading Julia packages",diff(timestamps)[end])
###########################################################################################################
# SECTION: PERFORM SAFETY CHECKS
println("------------------------------------------------ SECTION $(prnt) ------------------------------------------------"); prnt += 1

verbose && println("Performing safety checks... ")
if !(length(filepaths_S)==length(filepaths_W))
    error("Number of measurement files: $(length(filepaths_S)). Number of weight matrix files: $(length(filepaths_W)). Number of measurement files and number of weight matrix files must match. Please correct and re-try.")
end

append!(timestamps,time()) # The timestamp when the foundational safety check has been performed
dictionary_of_sections[prnt-1] = ("Performing safety checks",diff(timestamps)[end])
###########################################################################################################
# SECTION: LOAD MEASUREMENTS
println("------------------------------------------------ SECTION $(prnt) ------------------------------------------------"); prnt += 1

verbose && print("Loading measurements... ")
global S; S = Vector{Vector{Float64}}(undef,length(filepaths_S)) # Pre-allocate measurements for all diagnostics
global S_units; S_units = Vector{String}(undef,length(filepaths_S)) # Pre-allocate units of measurements for all diagnostics
global S_errs; S_errs = Vector{Vector{Float64}}(undef,length(filepaths_S)) # Pre-allocate errors (uncertainties) for all diagnostics
global S_errs_units; S_errs_units = Vector{String}(undef,length(filepaths_S)) # Pre-allocate units of errors (uncertainties) for all diagnostics
global S_abscissas; S_abscissas = Vector{Vector{Float64}}(undef,length(filepaths_S)) # Pre-allocate measurement bin centers for all diagnostics
global S_abscissas_units; S_abscissas_units = Vector{String}(undef,length(filepaths_S)) # Pre-allocate units of measurement bin centers for all diagnostics
for (i,f) in enumerate(filepaths_S)
    local myfile # Declare local scope. Just for clarity
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
    S_units[i] = read_func(myfile["S_units"])
    S_errs[i] = read_func(myfile["err"])
    S_errs_units[i] = read_func(myfile["err_units"])
    S_abscissas[i] = read_func(myfile["Ed_array"])
    S_abscissas_units[i] = read_func(myfile["Ed_array_units"])
    close(myfile)
end
verbose && println("Done!")


global ok # Declare global scope
ok = true # Simplifying print bool dummy variable
verbose && print("Ensuring measurements units consistency between S and err... ")
for (i,S_unit) in enumerate(S_units)
    if !(S_unit==S_errs_units[i])
        global ok; ok = false
        verbose && println("")
        verbose && print("Units of signal $(i): $(S_unit). Units of err $(i): $(S_errs_units[i]). Attempting to convert err to match S units... ")
        S_errs[i] = units_conversion_factor(S_errs_units[i],S_unit) .*S_errs[i]
        S_errs_units[i] = S_unit
        verbose && println("Success!")
    end
end
ok && verbose && println("Done!")

append!(timestamps,time()) # The timestamp when the measurements have been loaded completely
dictionary_of_sections[prnt-1] = ("Loading measurements",diff(timestamps)[end])
###########################################################################################################
# SECTION: LOAD WEIGHT FUNCTIONS
println("------------------------------------------------ SECTION $(prnt) ------------------------------------------------"); prnt += 1

verbose && print("Loading weight functions... ")
if isempty(min_array) || isempty(max_array) || isempty(n_array) # If any of them are empty...
    # Use filepaths_W[1] to determine dimensionality
    verbose && println("")
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
    local myfile # Declare local scope. Just for clarity
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
        Ed_array = read_func(myfile["Ed_array"])
        Ed_array_units = read_func(myfile["Ed_array_units"])

        if lowercase(coordinate_system)=="(vpara,vperp)"
            W_inflated[i] = read_func(myfile["W_vel"])
            D1_array = read_func(myfile["vpara_array"])
            D1_array_units = "m_s^-1" # from calc2DWeights.jl, the vpara grid points will always be in m_s^-1
            D2_array = read_func(myfile["vperp_array"])
            D2_array_units = "m_s^-1" # from calc2DWeights.jl, the vperp grid points will always be in m_s^-1
        elseif lowercase(coordinate_system)=="(e,p)"
            W_inflated[i] = read_func(myfile["W"])
            D1_array = read_func(myfile["E_array"])
            D1_array_units = "keV" # from calc2DWeights.jl, the energy grid points will always be in keV
            D2_array = read_func(myfile["p_array"])
            D2_array_units = "-" # from calc2DWeights.jl, the pitch grid points will always be dimensionless
        else
            @warn "Element number $(i) of input variable Vector 'scriptSources_W' was specified as $(scriptSources_W[i]), but input variable 'coordinate_system' was incorrectly specified ($(coordinate_system)). Therefore, a coordinate system of (E,p) will be automatically assumed."
            W_inflated[i] = read_func(myfile["W"])
            D1_array = read_func(myfile["E_array"])
            D1_array_units = "keV" # from calc2DWeights.jl, the energy grid points will always be in keV
            D2_array = read_func(myfile["p_array"])
            D2_array_units = "-" # from calc2DWeights.jl, the pitch grid points will always be dimensionless
        end
        close(myfile)
        W_abscissas[i] = [Ed_array, D1_array, D2_array]
        W_abscissas_units[i] = [Ed_array_units, D1_array_units, D2_array_units]
    elseif sswi=="orbweights_2dto4d" || sswi=="orbweights_2dto4d.jl" || sswi=="calcorbweights" || sswi=="calcorbweights.jl" # Account for misinterpretation by user (include .jl file extension). PLEASE NOTE! THESE FILES SHOULD REALLY BE OUTPUT FILES FROM THE orbWeights_2Dto4D.jl SCRIPT. THE CODE BELOW JUST TAKES INTO ACCOUNT THAT USERS MIGHT BE CONFUSED
        if haskey(myfile,"W2D")
            error("In the filepaths_W input array, the $(f) file is an output file of the calcOrbWeights.jl script. This is not a file that you can specify as input to the solveInverseProblem.jl script. Even though one might think so. Please specify an output file of the OWCF/helper/orbWeights_2Dto4D.jl script instead.")
        end
        if !haskey(myfile,"W")
            error("In the filepaths_W input array, the $(f) file does not have a 'W' file key, even though '$(scriptSources_W[i])' was specified in the scriptSources_W input array. Please correct and re-try.")
        end
        W_inflated[i] = read_func(myfile["W"])
        Ed_array = read_func(myfile["Ed_array"])
        Ed_array_units = read_func(myfile["Ed_array_units"])
        E_array = read_func(myfile["E_array"])
        E_array_units = "keV" # from calcOrbWeights.jl, the energy grid points will always be in keV
        pm_array = read_func(myfile["pm_array"])
        pm_array_units = "-" # from calcOrbWeights.jl, the pitch maximum grid points will always be dimensionless
        Rm_array = read_func(myfile["Rm_array"])
        Rm_array_units = "m" # from calcOrbWeights.jl, the major radius maximum grid points will always be in meters
        close(myfile)
        W_abscissas[i] = [Ed_array, E_array, pm_array, Rm_array]
        W_abscissas_units[i] = [Ed_array_units, E_array_units, pm_array_units, Rm_array_units]
    elseif sswi=="calc4dweights" || sswi=="calc4dweights.jl"
        # Due to its size, the calc4DWeights.jl case needs special treatment
        # This code can most likely only be run on clusters with a HUGE amount of RAM
        VAL = read_func(myfile["VAL"]) # Vector containing the non-zero values of the (E,p,R,z) weight matrix
        ROW = read_func(myfile["ROW"]) # Vector containing the corresponding row elements
        COL = read_func(myfile["COL"]) # Vector containing the corresponding column elements
        m_W = read_func(myfile["m_W"]) # Total number of rows (including zero elements not included in R and C) of the (E,p,R,z) weight matrix
        n_W = read_func(myfile["n_W"]) # Total number of columns (including zero elements not included in R and C) of the (E,p,R,z) weight matrix
        Ed_array = read_func(myfile["Ed_array"])
        Ed_array_units = read_func(myfile["Ed_array_units"])
        E_array = read_func(myfile["E_array"])
        E_array_units = "keV" # from calc4DWeights.jl, the energy grid points will always be in keV
        p_array = read_func(myfile["p_array"])
        p_array_units = "-" # from calc4DWeights.jl, the pitch grid points will always be dimensionless
        R_array = read_func(myfile["R_array"])
        R_array_units = "m" # from calc4DWeights.jl, the major radius grid points will always be in meters
        z_array = read_func(myfile["z_array"])
        z_array_units = "m" # from calc4DWeights.jl, the vertical coordinate grid points will always be in meters
        EpRz_coords = read_func(myfile["EpRz_coords"]) # The indices and (E,p,R,z) coordinates for all the columns of the (E,p,R,z) weight matrix (W_2D, see below)
        close(myfile)
        W_2D = dropzeros(sparse(append!(ROW,m_W),append!(COL,n_W),append!(VAL,0.0)))
        W_5D = zeros(length(Ed_array),length(E_array),length(p_array),length(R_array),length(z_array))
        verbose && println("Re-creating the 5D (Ed,E,p,R,z) weight function matrix... ")
        for iEd in 1:size(W_2D,1)
            for (i,c) in enumerate(EpRz_coords[:])
                W_5D[iEd,c[1][1],c[2][1],c[3][1],c[4][1]] = W_2D[iEd,i]
            end
        end
        W_inflated[i] = W_5D
        W_abscissas[i] = [Ed_array, E_array, p_array, R_array, z_array]
        W_abscissas_units[i] = [Ed_array_units, E_array_units, p_array_units, R_array_units, z_array_units]
    else
        # The general weight functions loading case
        W_inflated[i] = read_func(myfile["W"])
        abscissas = []
        units = []
        append!(abscissas,[read_func(myfile["Ed_array"])]) # Must have 'Ed_array' file key
        append!(units,[read_func(myfile["Ed_array_units"])]) # Must have 'Ed_array_units' file key
        append!(abscissas,[read_func(myfile["D1_array"])]) # Must have at least 'D1_array' file key
        append!(units,[read_func(myfile["D1_array_units"])]) # Must have at least 'D1_array_units' file key
        if haskey(myfile,"D2_array")
            append!(abscissas,[read_func(myfile["D2_array"])])
            append!(units,[read_func(myfile["D2_array_units"])])
        end
        if haskey(myfile,"D3_array")
            append!(abscissas,[read_func(myfile["D3_array"])])
            append!(units,[read_func(myfile["D3_array_units"])])
        end
        if haskey(myfile,"D4_array")
            append!(abscissas,[read_func(myfile["D4_array"])])
            append!(units,[read_func(myfile["D4_array_units"])])
        end
        if haskey(myfile,"D5_array")
            append!(abscissas,[read_func(myfile["D5_array"])])
            append!(units,[read_func(myfile["D5_array_units"])])
        end
        if haskey(myfile,"D6_array")
            append!(abscissas,[read_func(myfile["D6_array"])])
            append!(units,[read_func(myfile["D6_array_units"])])
        end
        close(myfile)
        W_abscissas[i] = [Float64.(ab) for ab in abscissas] # Transform to Vector{Float64} from Vector{Any}
        W_abscissas_units[i] = map(x-> "$(x)",units) # Transform to Vector{String} from Vector{Any}
    end
end
verbose && println("Done!")

append!(timestamps,time()) # The timestamp when the weight functions have been loaded completely
dictionary_of_sections[prnt-1] = ("Loading weight functions",diff(timestamps)[end])
###########################################################################################################
# SECTION: CHECK THAT ALL ABSCISSAS MATCH IN TERMS OF DIMENSIONS AND UNITS.
# OTHERWISE, CORRECT THEM SO THAT ALL MATCH.
println("------------------------------------------------ SECTION $(prnt) ------------------------------------------------"); prnt += 1

verbose && println("Performing dimension and units checks for the weight functions... ")
for i in eachindex(filepaths_W)
    if !(length(W_abscissas[i])==length(W_abscissas[1]))
        error("Number of abscissas found in $(filepaths_W[1]): $(length(W_abscissas[1])). Number of abscissas found in $(filepaths_W[i]): $(length(W_abscissas[i])). Number of abscissas must match for all weight functions files. Please correct and re-try.")
    end

    # Check if the units of the measurement bin centers of the weight functions match the units of the measurement bin centers of the signal
    # If not, convert both the abscissa units and weight functions units (since weight functions have units signal/ion)
    if !(W_abscissas_units[i][1]==S_abscissas_units[i])
        verbose && println("---> Units of abscissa (measurement bin centers) of $(filepaths_S[i]) is $(S_abscissas_units[i]).")
        verbose && println("---> Units of FIRST abscissa (measurement bin centers) of $(filepaths_W[i]) is $(W_abscissas_units[i][1]).") 
        verbose && println("------> Converting weight function and its first abscissa to match units of signal abscissa... ")
        ucf = units_conversion_factor(W_abscissas_units[i][1],S_abscissas_units[i])
        W_inflated[i] = (1/ucf) .*W_inflated[i] # Since weight functions have units signal/ion and the signal will have 1/units of the abscissa, multiply by the inverse of the units conversion factor
        W_abscissas[i][1] = ucf .*W_abscissas[i][1]
        W_abscissas_units[i][1] = S_abscissas_units[i]
    end

    # Check that the units of the reconstruction space abscissas of all weight functions match
    if i>1 # Only compare with previous file if i>1. Otherwise, will cause out-of-bounds error
        for j in eachindex(W_abscissas[i]) # For each reconstruction space dimension (j>1)
            if (j>1) && !(W_abscissas_units[i][j]==W_abscissas_units[i-1][j]) # If the units of the same dimension for different weight function files don't match (except for the measurement bin center units of course (they will almost always differ), hence the j>1 condition)...
                verbose && println("Units of abscissa $(j) in $(filepaths_W[i]) and $(filepaths_W[i-1]) do not match. Converting abscissa $(j) in $(filepaths_W[i]) to match the units of abscissa $(j) in $(filepaths_W[i-1])... ")
                ucf = units_conversion_factor(W_abscissas_units[i][j],W_abscissas_units[i-1][j])
                W_abscissas[i][j] = ucf .*W_abscissas[i][j]
                W_abscissas_units[i][j] = W_abscissas_units[i][j]
            end
        end
    end
end

append!(timestamps,time()) # The timestamp when the abscissas units checks have been performed
dictionary_of_sections[prnt-1] = ("Performing abscissa units checks",diff(timestamps)[end])
###########################################################################################################
# SECTION: INTERPOLATE ONTO THE GRID SPECIFIED BY min_array, max_array and n_array 
# IF THEY ARE NOT SPECIFIED, USE THE ABSCISSAS OF THE FIRST WEIGHT MATRIX IN filepaths_W
println("------------------------------------------------ SECTION $(prnt) ------------------------------------------------"); prnt += 1

verbose && println("Attempting to interpolate all weight functions onto a common reconstruction space grid... ")
verbose && println("---> Creating interpolation grid (query points)")
global query_vecs_n_inds; query_vecs_n_inds = () # A tuple to hold all query point vectors and their indices. Structure: ((vector,indices),(vector,indices),...) with length equal to reconstruction grid dimension
if isempty(min_array) || isempty(max_array) || isempty(n_array) # If any of them are empty...
    verbose && println("------> min_array, max_array and/or n_array were not specified. Data grid of $(filepaths_W[1]) will be used")
    for abscissa in W_abscissas[1][2:end] # Use all reconstruction space abscissas of the first weight function file as reconstruction grid...
        global query_vecs_n_inds # Use global scope variable
        query_vecs_n_inds = tuple(query_vecs_n_inds[:]...,collect(zip(abscissa,1:length(abscissa)))) # Add the (vector,indices) pairs one by one  into a big tuple (tuples are immutable, hence the cumbersome code)
    end
else
    verbose && println("------> min_array, max_array and n_array were specified. Utilizing... ")
    for i in eachindex(n_array) # For all reconstruction grid dimensions... 
        global query_vecs_n_inds # Use global scope variable
        query_vecs_n_inds = tuple(query_vecs_n_inds[:]...,collect(zip(collect(range(min_array[i],stop=max_array[i],length=n_array[i])),1:n_array[i]))) # Add the (vector,indices) pairs one by one  into a big tuple (tuples are immutable, hence the cumbersome code)
    end
end
query_points_n_coords = Iterators.product(query_vecs_n_inds...) # Create a long list of all reconstruction space grid points and their coordinates by computing a product between all query point-index vectors. Example structure (if 3 dimensions): [((x1_1,1),(x2_1,1),(x3_1,1)),((x1_2,2),(x2_1,1),(x3_1,1)),...]

for (i,w_inflated) in enumerate(W_inflated)
    w_inflated_interpolated = zeros(tuple(size(w_inflated,1),(length.(query_vecs_n_inds))...)) # Pre-allocate a new inflated weight function, to store the interpolated values
    nodess = () # The nodes of the interpolation object (nodess name only due to Scope warnings)
    for abscissa in W_abscissas[i][2:end] # Use the reconstruction space abscissas of the weight functions as the nodes of the interpolation object...
        nodess = tuple(nodess[:]...,abscissa) # Put them all in a tuple (tuples are immutable, hence the cumbersome code)
    end
    node_coords = CartesianIndices(length.(nodess)) # A trick to specify all node coordinates of an array of general size
    verbose && println("---> Interpolating weight matrix $(i) of $(length(W_inflated))... ")
    for j=1:size(w_inflated,1)
        local itp; itp = Interpolations.interpolate(nodess,w_inflated[j,node_coords],Gridded(Linear()))
        local etp; etp = Interpolations.extrapolate(itp,Interpolations.Flat()) # If outside of interpolation region, use edge values to extrapolate
        for query_point_n_coord in query_points_n_coords
            point = map(x-> x[1],query_point_n_coord) # The point to interpolate at. E.g. (100.0,0.3) in energy (keV),pitch
            coord = map(x-> x[2],query_point_n_coord) # The coordinate of that point. E.g. (53,14)
            w_inflated_interpolated[j,coord...] = etp(point...)
        end
    end
    W_inflated[i] = w_inflated_interpolated # Replace the non-interpolated with the interpolated
    W_abscissas[i] = [W_abscissas[i][1], map(x-> map(xx-> xx[1],x), query_vecs_n_inds)...] # Replace the original grid points with the query points
end

append!(timestamps,time()) # The timestamp when the weight functions have been interpolated onto a common reconstruction space grid
dictionary_of_sections[prnt-1] = ("Interpolating weight functions onto common reconstruction space grid",diff(timestamps)[end])
###########################################################################################################
# SECTION: RESHAPE ALL INFLATED WEIGHT MATRICES INTO THEIR 2D SHAPE, TO BE USED IN INVERSE PROBLEMS
# THEN, INTERPOLATE EACH WEIGHT MATRIX W[i] TO MATCH THE MEASUREMENT BIN CENTER GRID OF THE CORRESPONDING SIGNAL S[i]
println("------------------------------------------------ SECTION $(prnt) ------------------------------------------------"); prnt += 1

verbose && println("For each diagnostic, restructuring weight functions into 2D matrix... ")
global W; W = Vector{Array{Float64,2}}(undef,length(filepaths_W)) # The weight matrices, for each diagnostic
for (i,w_inflated) in enumerate(W_inflated)
    verbose && println("---> Restructuring weight functions, creating weight matrix $(i) of $(length(W_inflated))... ")
    ws = size(w_inflated)
    W[i] = reshape(w_inflated,(ws[1],reduce(*,ws[2:end])))
end
W_inflated = nothing # Clear memory. To minimize memory usage

verbose && println("For each diagnostic, interpolating the weight matrix to match measurement data grid... ")
for (i,w) in enumerate(W)
    verbose && println("---> Interpolating weight matrix $(i) of $(length(W))... ")
    verbose && println("------> extrema(S abscissa grid): $(extrema(S_abscissas[i])) $(S_abscissas_units[i])")
    verbose && println("------> extrema(W abscissa grid): $(extrema(W_abscissas[i][1])) $(W_abscissas_units[i][1])")
    if !(W_abscissas_units[i][1]==S_abscissas_units[i])
        @warn "Units for measurement data grid ($(S_abscissas_units[i])) do not match units of weight matrix measurement bin centers ($(W_abscissas_units[i][1])). This should be impossible. Please contact the admins of the OWCF and attach the start file, to have the error investigated."
    end
    query_points = S_abscissas[i]
    w_new = zeros(length(query_points),size(w,2))
    for (ic,col) in enumerate(eachcol(w))
        local nodes; nodes = (W_abscissas[i][1],)
        local itp; itp = Interpolations.interpolate(nodes,col,Gridded(Linear()))
        local etp; etp = Interpolations.extrapolate(itp,0.0) # Assume zero sensitivity (no knowledge) outside of known spectrum range
        w_new[:,ic] = [etp(qp) for qp in query_points]
    end
    W[i] = w_new
    W_abscissas[i][1] = query_points
    W_abscissas_units[i][1] = S_abscissas_units[i] # This should already be the case
end

append!(timestamps,time()) # The timestamp when the weight functions have been reshaped into 2D and interpolated onto the corresponding measurements grid
dictionary_of_sections[prnt-1] = ("Constructing 2D weight matrices from weight functions, and interpolating to match diagnostics measurement bin centers",diff(timestamps)[end])
###########################################################################################################
# SECTION: IF rescale_W WAS SET TO true IN THE START FILE, LOAD OR COMPUTE FAST-ION DISTRIBUTION(S).
# USE THIS/THESE DISTRIBUTION(S) TOGETHER WITH THE WEIGHT MATRICES TO COMPUTE REFERENCE MEASUREMENTS.
# COMPUTE RESCALING FACTORS TO BE ABLE TO RESCALE WEIGHT FUNCTIONS TO HAVE THE REFERENCE MEASUREMENTS MATCH 
# THE EXPERIMENTAL MEASUREMENTS.
println("------------------------------------------------ SECTION $(prnt) ------------------------------------------------"); prnt += 1

if rescale_W
    verbose && println("Computing rescale factors for weight functions... ")
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
        w_abscissas_units = W_abscissas_units[1][2:end] # Units of reconstruction space (2:end) abscissas of first weight function matrix are the units of the abscissas of all weight functions matrices, because of section 5
        w_abscissas = W_abscissas[1][2:end] # Reconstruction space (2:end) abscissas of first weight function matrix are the abscissas of all weight functions matrices, because of section 5

        verbose && println("---> rescale_W_F_ref_source=$(rescale_W_F_ref_source) was specified... ")
        verbose && println("------> Checking weight function reconstruction space compatibility. Currently supported coordinate systems are $(currently_supported_abscissas)... ")
        Ep_bool = is_energy_pitch(w_abscissas, w_abscissas_units)
        VEL_bool = is_vpara_vperp(w_abscissas, w_abscissas_units)
        EpRz_bool = is_EpRz(w_abscissas, w_abscissas_units)
        # ADD MORE CHECKS HERE IN FUTURE VERSIONS, IF NEEDED <---------------------------------------------------------------------
        if Ep_bool
            w_rec_space = "(E,p)"
            dummy_bool, w_energy_ind, w_pitch_ind = is_energy_pitch(w_abscissas, w_abscissas_units; returnExtra=true)
            R_of_interests = vcat(units_conversion_factor(R_of_interest_units,"m")*R_of_interest) # R_of_interests will be the same, regardless of 2D coordinates
            z_of_interests = vcat(units_conversion_factor(z_of_interest_units,"m")*z_of_interest) # z_of_interests will be the same, regardless of 2D coordinates
        elseif VEL_bool
            w_rec_space = "(vpara,vperp)"
            dummy_bool, w_vpara_ind, w_vperp_ind = is_vpara_vperp(w_abscissas, w_abscissas_units; returnExtra=true)
            R_of_interests = vcat(units_conversion_factor(R_of_interest_units,"m")*R_of_interest) # R_of_interests will be the same, regardless of 2D coordinates
            z_of_interests = vcat(units_conversion_factor(z_of_interest_units,"m")*z_of_interest) # z_of_interests will be the same, regardless of 2D coordinates
        elseif EpRz_bool
            w_rec_space = "(E,p,R,z)"
            dummy_bool, w_energy_ind, w_pitch_ind, w_R_ind, w_z_ind, R_of_interests, z_of_interests = is_EpRz(w_abscissas, w_abscissas_units; returnExtra=true)
            R_of_interests = units_conversion_factor(w_abscissas_units[w_R_ind],"m") .*R_of_interests # Need unit of measurement meter for CDFto4D function (if that will be used)
            z_of_interests = units_conversion_factor(w_abscissas_units[w_z_ind],"m") .*z_of_interests # Need unit of measurement meter for CDFto4D function (if that will be used)
        else
            error(standard_rescale_W_error)
        end

        verbose && println("------> Loading F_ref from file $(rescale_W_F_file_path)... ")
        file_ext = lowercase(split(rescale_W_F_file_path,".")[end]) # Assume final part of filepath after "." is the file extension
        if file_ext=="jld2"
            F_ref, E_ref, p_ref, R_ref, z_ref = JLD2to4D(rescale_W_F_file_path)
        elseif file_ext=="hdf5" || file_ext=="h5"
            F_ref, E_ref, p_ref, R_ref, z_ref = h5to4D(rescale_W_F_file_path;rowmajor=h5_is_rowmajor)
        elseif file_ext=="cdf"
            F_ref, E_ref, p_ref, R_ref, z_ref = CDFto4D(rescale_W_F_file_path, R_of_interests, z_of_interests; btipsign=btipsign)
        else
            error("Input variable 'rescale_W_F_file_path' has unknown file extension ($(file_ext)). Accepted file extensions are .jld2, .hdf5, .h5 and .cdf. Please correct and re-try.")
        end

        if maximum(R_ref)>100.0 # Assume no tokamak has a major radius as large as 100 m...
            verbose && println("------> Loaded F_ref from file $(rescale_W_F_file_path) is assumed to be given in cm^-3. Converting to m^-3... ")
            R_ref = units_conversion_factor("cm","m") .* R_ref
            z_ref = units_conversion_factor("cm","m") .* z_ref
            F_ref = units_conversion_factor("cm^-3","m^-3") .* F_ref
        end

        if maximum(E_ref)>1.0e6 # Assume no tokamak has a fast-ion distribution with a maximum energy as great as 1.0e6 keV
            verbose && println("------> Loaded F_ref from file $(rescale_W_F_file_path) is assumed to be given in eV^-1. Converting to keV^-1... ")
            E_ref = units_conversion_factor("eV","keV") .* E_ref
            F_ref = units_conversion_factor("eV^-1","keV^-1") .* F_ref
        end

        if w_rec_space=="(E,p)"
            verbose && println("------> Computing W*F_ref assuming weight functions are given in (E,p) coordinates (or equivalent)... ")
            w_energy_units = w_abscissas_units[w_energy_ind]
            w_pitch_units = w_abscissas_units[w_pitch_ind]
            w_energy_abscissa = w_abscissas[w_energy_ind]; nwE = length(w_energy_abscissa)
            w_pitch_abscissa = w_abscissas[w_pitch_ind]; nwp = length(w_pitch_abscissa)

            query_points_n_coords = Iterators.product(zip(w_energy_abscissa,1:nwE),zip(w_pitch_abscissa,1:nwp),zip(R_of_interests,1:1),zip(z_of_interests,1:1)) # R and z should already be 1-element Vectors with unit of measurement 'm'
            F_ref_interpolated = zeros(nwE,nwp,1,1) # Pre-allocate the interpolated f(E,p) distribution

            E_ref = units_conversion_factor("keV",w_energy_units) .*E_ref # Convert the energy units from keV (which we know the JLD2to4D, h5to4D and CDFto4D functions output) to w_energy_units
            F_ref = units_conversion_factor("keV^-1",units_inverse(w_energy_units)) .*F_ref # Also, convert the fast-ion distribution using the inverse of the w_energy_units

            nodes = (E_ref,p_ref,R_ref,z_ref)
            itp = Interpolations.interpolate(nodes,F_ref,Gridded(Linear()))
            etp = Interpolations.extrapolate(itp,Interpolations.Flat()) # If outside of interpolation region, use edge values to extrapolate
            for query_point_n_coord in query_points_n_coords
                point = map(x-> x[1],query_point_n_coord) # The point to interpolate at. E.g. (100.0,0.3) in energy (keV),pitch
                coord = map(x-> x[2],query_point_n_coord) # The coordinate of that point. E.g. (53,14)
                F_ref_interpolated[coord...] = etp(point...)
            end
            F_ref_interpolated = dropdims(F_ref_interpolated,dims=(3,4)) # From shape (nwE,nwp,1,1) to (nwE,nwp)
            F_ref_interp_1D = reshape(F_ref_interpolated,(nwE*nwp,))

            rescale_W_factors = ones(length(filepaths_W)) # Pre-allocate weight re-scale factors
            for i in eachindex(filepaths_W)
                WF = W[i]*F_ref_interp_1D # The weight function matrix (2D) multiplied with the (1D/vectorized) reference fast-ion distribution
                rescale_W_factors[i] = rescale_func(S[i])/rescale_func(WF) # S[i] is the experimental data vector of (fast-ion) diagnostic 'i'
            end
        elseif w_rec_space=="(vpara,vperp)"
            verbose && println("------> Using reference fast-ion distribution assuming weight functions are given in (vpara,vperp) coordinates (or equivalent)... ")
            w_vpara_unit = w_abscissas_units[w_vpara_ind]
            w_vperp_unit = w_abscissas_units[w_vperp_ind]
            w_vpara_abscissa = w_abscissas[w_vpara_ind]; nwvpa = length(w_vpara_abscissa)
            w_vperp_abscissa = w_abscissas[w_vperp_ind]; nwvpe = length(w_vpara_abscissa)
            nfE = length(E_ref); nfp = length(p_ref)

            nodes = (E_ref,p_ref,R_ref,z_ref)
            itp = Interpolations.interpolate(nodes,F_ref,Gridded(Linear()))
            etp = Interpolations.extrapolate(itp,Interpolations.Flat()) # If outside of interpolation region, use edge values to extrapolate
            query_points_n_coords = Iterators.product(zip(E_ref,1:nfE),zip(p_ref,1:nfp),zip(R_of_interests,1:1),zip(z_of_interests,1:1))
            F_ref_interpolated = zeros(nfE,nfp,1,1) # Pre-allocate the interpolated f(E,p) distribution
            for query_point_n_coord in query_points_n_coords
                point = map(x-> x[1],query_point_n_coord) # The point to interpolate at. E.g. (100.0,0.3) in energy (keV),pitch
                coord = map(x-> x[2],query_point_n_coord) # The coordinate of that point. E.g. (53,14)
                F_ref_interpolated[coord...] = etp(point...)
            end
            F_ref_interpolated = dropdims(F_ref_interpolated,dims=(3,4)) # From shape (nfE,nfp,1,1) to (nfE,nfp)

            F_ref_VEL, vpara_ref, vperp_ref = Ep2VparaVperp(E_ref, p_ref, F_ref_interpolated; my_gcp=FI_species, needJac=true, returnAbscissas=true)
            vpara_ref = units_conversion_factor("m_s^-1",w_vpara_unit) .*vpara_ref # Match the units of w_abscissas
            vperp_ref = units_conversion_factor("m_s^-1",w_vperp_unit) .*vperp_ref # Match the units of w_abscissas
            F_ref_VEL = units_conversion_factor("m^-2_s^2",units_inverse(w_vpara_unit)*"_"*units_inverse(w_vperp_unit)) .* F_ref_VEL

            nodes = (vpara_ref,vperp_ref)
            itp = Interpolations.interpolate(nodes,F_ref_VEL,Gridded(Linear()))
            etp = Interpolations.extrapolate(itp,Interpolations.Flat())
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
            verbose && println("------> Using reference fast-ion distribution assuming weight functions are given in (E,p,R,z) coordinates (or equivalent)... ")
            w_energy_units = w_abscissas_units[w_energy_ind]
            w_pitch_units = w_abscissas_units[w_pitch_ind]
            w_R_units = w_abscissas_units[w_R_ind]
            w_z_units = w_abscissas_units[w_z_ind]
            w_energy_abscissa = w_abscissas[w_energy_ind]; nwE = length(w_energy_abscissa)
            w_pitch_abscissa = w_abscissas[w_pitch_ind]; nwp = length(w_pitch_abscissa)
            w_R_abscissa = w_abscissas[w_R_ind]; nwR = length(w_R_abscissa)
            w_z_abscissa = w_abscissas[w_z_ind]; nwz = length(w_z_abscissa)

            E_ref = units_conversion_factor("keV",w_energy_units) .*E_ref # Convert the energy units from keV (which we know the JLD2to4D, h5to4D and CDFto4D functions output) to w_energy_units
            R_ref = units_conversion_factor("m",w_R_units) .*R_ref # Convert the R units from m (which we know the JLD2to4D, h5to4D and CDFto4D functions output) to w_R_units
            z_ref = units_conversion_factor("m",w_z_units) .*z_ref # Convert the z units from m (which we know the JLD2to4D, h5to4D and CDFto4D functions output) to w_z_units
            F_ref = units_conversion_factor("keV^-1_m^-3",units_inverse(w_energy_units)*"_"*units_inverse(w_R_units)*"_"*units_inverse(w_z_units)) .*F_ref # Also, convert the fast-ion distribution from keV^-1_m^-3 to the inverse of w_energy_units, w_R_units and w_z_units 

            query_points_n_coords = Iterators.product(zip(w_energy_abscissa,1:nwE),zip(w_pitch_abscissa,1:nwp),zip(w_R_abscissa,1:nwR),zip(w_z_abscissa,1:nwz))
            F_ref_interpolated = zeros(nwE,nwp,nwR,nwz) # Pre-allocate the interpolated f(E,p,R,z) distribution

            nodes=(E_ref,p_ref,R_ref,z_ref)
            itp = Interpolations.interpolate(nodes,F_ref,Gridded(Linear()))
            etp = Interpolations.extrapolate(itp,Interpolations.Flat())
            for query_point_n_coord in query_points_n_coords
                point = map(x-> x[1],query_point_n_coord) # The point to interpolate at. E.g. (120.0,0.64,3.1,0.2) in E [keV], p [-], R [m], z [m] 
                coord = map(x-> x[2],query_point_n_coord) # The coordinate of that point. E.g. (53,34,15,16)
                F_ref_interpolated[coord...] = etp(point...)
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
            mu = [m for m in gauss_coord[1:l/2]] # The first half will be from F_ref_abscissas (Tuple to Vector)
            std = [s for s in gauss_coord[(l/2)+1:end]] # The second half will be from std_abscissas (Tuple to Vector)
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

append!(timestamps,time()) # The timestamp when the weight function rescale factors have been computed
dictionary_of_sections[prnt-1] = ("Computing weight functions rescale factors (if requested)",diff(timestamps)[end])
###########################################################################################################
# SECTION: IF :COLLISIONS WAS INCLUDED IN THE regularization INPUT VARIABLE, AND THE RECONSTRUCTION 
# SPACE IS A SPACE IN WHICH COLLISIONS REGULARIZATION IS SUPPORTED, COMPUTE SLOWING-DOWN BASIS FUNCTIONS
println("------------------------------------------------ SECTION $(prnt) ------------------------------------------------"); prnt += 1

if "collisions" in lowercase.(String.(regularization))
    COMPATIBLE_OPTIONS = [is_energy_pitch, is_vpara_vperp] # ADD MORE CHECKS HERE IN FUTURE VERSIONS, IF NEEDED <---------------------------------------------------------------------
    verbose && print(":COLLISIONS included in 'regularization' input variable. Performing safety checks... ")
    if sum([f(W_abscissas[1][2:end], W_abscissas_units[1][2:end]) for f in COMPATIBLE_OPTIONS])==0 # Is the reconstruction space NOT in any of the COMPATIBLE_OPTIONS?
        error(":COLLISIONS was specified in 'regularization' input variable, but reconstruction space is not compatible. Current compatible options include (E,p) and (vpara,vperp). Please correct and re-try.")
    end
    if !(length(regularization_thermal_ion_temp)==length(regularization_thermal_ion_dens)==length(regularization_thermal_ion_species))
        error("The lengths of the 'regularization_thermal_ion_temp', 'regularization_thermal_ion_dens' and 'regularization_thermal_ion_species' input variables were not equal. Please correct and re-try.")
    end
    verbose && println("ok!")
    verbose && print("---> Loading magnetic equilibrium from $(regularization_equil_filepath)... ")
    M, wall = nothing, nothing # Initialize global magnetic equilibrium variables
    try
        global M; global wall # Declare global scope
        M, wall = read_geqdsk(regularization_equil_filepath,clockwise_phi=false) # Assume the phi-direction is pointing counter-clockwise when tokamak is viewed from above. This is true for almost all coordinate systems used in the field of plasma physics
    catch # Otherwise, assume magnetic equilibrium is a saved .jld2 file
        global M; global wall; local myfile
        myfile = jldopen(regularization_equil_filepath,false,false,false,IOStream)
        M = myfile["S"]
        wall = myfile["wall"]
        close(myfile)
    end
    psi_axis, psi_bdry = psi_limits(M)
    verbose && println("ok!")

    verbose && println("---> Computing normalized poloidal flux coordinate (rho_pol) for (R_of_interest, z_of_interest) input variables... ")
    if length(W_abscissas_units[1])==3
        R_of_interests = vcat(units_conversion_factor(R_of_interest_units,"m") .* R_of_interest) # Converting to meters to match Equilibrium.jl package standard
        z_of_interests = vcat(units_conversion_factor(z_of_interest_units,"m") .* z_of_interest) # Converting to meters to match Equilibrium.jl package standard
        rho_of_interests = vcat(sqrt((M(R_of_interest,z_of_interest)-psi_axis)/(psi_bdry-psi_axis)))
    else
        # CODE FUTURE CASES HERE (3D, 4D ETC)
    end

    verbose && println("---> Preparing thermal ion data... ")
    thermal_ion_temperatures = zeros(length(regularization_thermal_ion_species),length(rho_of_interests)) # Ready for future versions with many rho_of_interest values
    thermal_ion_densities = zeros(length(regularization_thermal_ion_species),length(rho_of_interests)) # Ready for future versions with many rho_of_interest values
    
    for (i,ion) in enumerate(regularization_thermal_ion_species)
        verbose && println("------> Preparing thermal ion temperature data for species $(ion)... ")
        if (typeof(regularization_thermal_ion_temp[i]) <: Real)
            verbose && println("---------> Found an ion temperature of $(regularization_thermal_ion_temp[i]) keV. Incorporating... ")
            thermal_ion_temperatures[i,:] = repeat(vcat(regularization_thermal_ion_temp[i]),length(rho_of_interests))
        elseif (typeof(regularization_thermal_ion_temp[i]) <: String)
            if !isfile(regularization_thermal_ion_temp[i])
                error("regularization_thermal_ion_temp element with index $(i) was not specified correctly. $(regularization_thermal_ion_temp[i]) is not a valid filepath. Please correct and re-try.")
            end
            local file_ext; file_ext = split(regularization_thermal_ion_temp[i],".")[end] # Assume last part of String after "." is the file extension
            if lowercase(file_ext)=="jld2"
                local myfile; myfile = jldopen(regularization_thermal_ion_temp[i],false,false,false,IOStream)
                local rho_array; rho_array = myfile["rho_pol"]
                local temp_array; temp_array = myfile["thermal_temp"]
                close(myfile)
                local itp; itp = Interpolations.interpolate((rho_array,),temp_array,Gridded(Linear()))
                local etp; etp = Interpolations.extrapolate(itp,0.0) # Assume vacuum scrape-off layer
                T_is = etp.(rho_of_interests) # T_i values at the rho_pol values of interest
                verbose && println("---------> Computed ion temperature(s) of $(T_is) keV")
                verbose && println("---------> at ρ_{pol}=$(round.(rho_of_interests,sigdigits=4))")
                verbose && println("---------> i.e. R = $(round.(R_of_interests,sigdigits=4)) meters")
                verbose && println("---------> i.e. z = $(round.(z_of_interests,sigdigits=4)) meters")
                verbose && println("---------> by loading $(regularization_thermal_ion_temp[i]) . Incorporating temperature(s)... ")
                thermal_ion_temperatures[i,:] = T_is
            elseif lowercase(file_ext)=="cdf"
                if (typeof(regularization_timepoint) <: Real)
                    local t_mid; t_mid = Float64.(regularization_timepoint)
                elseif (typeof(regularization_timepoint) <: String)
                    if !isfile(regularization_timepoint)
                        local t_mid; t_mid = regularization_timepoint # Assume format "XX.XXXX" or "XX,XXXX"
                    else
                        local t_mid; t_mid = (read_ncdf(regularization_timepoint,wanted_keys=["TIME"]))["TIME"] # Assume TRANSP-NUBEAM output file
                    end
                else
                    error("regularization_timepoint input variable has invalid type. Expected Float64, Int64 or String. regularization_timepoint input variable needed to load ion temperature data from $(regularization_thermal_ion_temp[i]). Please correct and re-try.")
                end
                if typeof(t_mid)==String
                    t_mid = parse(Float64,replace(t_mid, "," => "."))
                end
                local etp; etp = getTempProfileFromTRANSP(t_mid,regularization_thermal_ion_temp[i],regularization_thermal_ion_species[i])
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
            local file_ext; file_ext = split(regularization_thermal_ion_dens[i],".")[end] # Assume last part of String after "." is the file extension
            if lowercase(file_ext)=="jld2"
                local myfile; myfile = jldopen(regularization_thermal_ion_dens[i],false,false,false,IOStream)
                local rho_array; rho_array = myfile["rho_pol"]
                local dens_array; dens_array = myfile["thermal_dens"]
                close(myfile)
                local itp; itp = Interpolations.interpolate((rho_array,),dens_array,Gridded(Linear()))
                local etp; etp = Interpolations.extrapolate(itp,0.0) # Assume vacuum scrape-off layer
                n_is = etp.(rho_of_interests) # n_i values at the rho_pol values of interest
                verbose && println("---------> Computed ion densitie(s) of $(n_is) m^-3")
                verbose && println("---------> at ρ_{pol}=$(round.(rho_of_interests,sigdigits=4))")
                verbose && println("---------> i.e. R = $(round.(R_of_interests,sigdigits=4)) meters")
                verbose && println("---------> i.e. z = $(round.(z_of_interests,sigdigits=4)) meters")
                verbose && println("---------> by loading $(regularization_thermal_ion_dens[i]) . Incorporating densitie(s)... ")
                thermal_ion_densities[i,:] = n_is
            elseif lowercase(file_ext)=="cdf"
                if (typeof(regularization_timepoint) <: Real)
                    local t_mid; t_mid = Float64.(regularization_timepoint)
                elseif (typeof(regularization_timepoint) <: String)
                    if !isfile(regularization_timepoint)
                        local t_mid; t_mid = regularization_timepoint # Assume format "XX.XXXX" or "XX,XXXX"
                    else
                        local t_mid; t_mid = (read_ncdf(regularization_timepoint,wanted_keys=["TIME"]))["TIME"] # Assume TRANSP-NUBEAM output file
                    end
                else
                    error("regularization_timepoint input variable has invalid type. Expected Float64, Int64 or String. regularization_timepoint input variable needed to load ion density data from $(regularization_thermal_ion_dens[i]). Please correct and re-try.")
                end
                if typeof(t_mid)==String
                    t_mid = parse(Float64,replace(t_mid, "," => "."))
                end
                local etp; etp = getDensProfileFromTRANSP(t_mid,regularization_thermal_ion_dens[i],regularization_thermal_ion_species[i])
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
            if typeof(t_mid)==String
                t_mid = parse(Float64,replace(t_mid, "," => "."))
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
    if is_energy_pitch(W_abscissas[1][2:end],W_abscissas_units[1][2:end])
        n_e = thermal_electron_densities[1] # We know it must be a 1-element Vector because reconstruction space is (E,p)
        T_e = thermal_electron_temperatures[1] # -||-
        thermal_ion_densities_at_rho = thermal_ion_densities[:,1] # -||-
        thermal_ion_temperatures_at_rho = thermal_ion_temperatures[:,1] # -||-

        verbose && println("for (E,p) space with grid point extrema:")
        w_E_abscissa, w_E_ind = get_energy_abscissa(W_abscissas[1][2:end],W_abscissas_units[1][2:end]; returnExtra=true) # First dimension is always the measurement bin centers
        verbose && println("------> (E_min,E_max)=$(round.(extrema(w_E_abscissa),sigdigits=4)) $(W_abscissas_units[1][w_E_ind+1])") # +1 because get_energy_abscissa() returns index shifted by 1, due to not including measurement bin centers abscissa as function input
        w_p_abscissa = get_pitch_abscissa(W_abscissas[1][2:end],W_abscissas_units[1][2:end]) # First dimension is always the measurement bin centers
        verbose && println("------> (p_min,p_max)=$(round.(extrema(w_p_abscissa),sigdigits=4))")
        dE = diff(w_E_abscissa)[1]; dp = diff(w_p_abscissa)[1] # Assume equidistant (E,p) grid

        E_inj_array = w_E_abscissa[1:1:end-1] .+ diff(w_E_abscissa)[1]/2 # Assume uniform energy grid
        p_inj_array = w_p_abscissa[1:1:end-1] .+ diff(w_p_abscissa)[1]/2 # Assume uniform pitch grid

        verbose && println("------> Computing... ")
        F_SD = zeros(length(w_E_abscissa)*length(w_p_abscissa),length(E_inj_array)*length(p_inj_array))
        F_SD_coords = CartesianIndices((length(E_inj_array),length(p_inj_array)))
        for (ic,c) in enumerate(F_SD_coords)
            iE = c[1]; ip = c[2]; E_inj = E_inj_array[iE]; p_inj = p_inj_array[ip]
            f_SD = slowing_down_function(E_inj, p_inj, w_E_abscissa, w_p_abscissa, n_e, T_e, FI_species, regularization_thermal_ion_species, thermal_ion_densities_at_rho, thermal_ion_temperatures_at_rho; dampen=true, damp_type=:erfc, sigma=5.0)
            f_SD = f_SD ./sum((dE*dp) .*f_SD) # Normalize the basis function so they integrate to 1.0
            F_SD[:,ic] .= reshape(f_SD,(length(w_E_abscissa)*length(w_p_abscissa),1))
        end
    elseif is_vpara_vperp(W_abscissas[1][2:end],W_abscissas_units[1][2:end])
        v2E_rel = (v-> inv(1000*GuidingCenterOrbits.e0)*getSpeciesMass(FI_species)*((GuidingCenterOrbits.c0)^2)*(sqrt(1/(1-(v/(GuidingCenterOrbits.c0))^(2)))-1)) # A one-line function to transform from relativistic speed (m/s) to energy (keV)

        n_e = thermal_electron_densities[1] # We know it must be a 1-element Vector because reconstruction space is (vpara,vperp)
        T_e = thermal_electron_temperatures[1] # -||-
        thermal_ion_densities_at_rho = thermal_ion_densities[:,1] # -||-
        thermal_ion_temperatures_at_rho = thermal_ion_temperatures[:,1] # -||-

        verbose && println("for (vpara,vperp) space with grid point extrema:")
        w_vpara_abscissa, w_vpara_ind = get_vpara_abscissa(W_abscissas[1][2:end],W_abscissas_units[1][2:end]; returnExtra=true); nwvpa = length(w_vpara_abscissa)
        verbose && println("------> (vpara_min,vpara_max)=$(round.(extrema(w_vpara_abscissa),sigdigits=4)) $(W_abscissas_units[1][w_vpara_ind+1])") # +1 because get_vpara_abscissa() returns index shifted by 1, due to not including measurement bin centers abscissa as function input
        w_vperp_abscissa, w_vperp_ind = get_vperp_abscissa(W_abscissas[1][2:end],W_abscissas_units[1][2:end]; returnExtra=true); nwvpe = length(w_vperp_abscissa)
        verbose && println("------> (vperp_min,vperp_max)=$(round.(extrema(w_vperp_abscissa),sigdigits=4)) $(W_abscissas_units[1][w_vperp_ind+1])") # +1 because get_vperp_abscissa() returns index shifted by 1, due to not including measurement bin centers abscissa as function input
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
                global min_E; min_E = E
            end
            if E>max_E
                global max_E; max_E = E
            end
            if p<min_p
                global min_p; min_p = p
            end
            if p>max_p
                global max_p; max_p = p
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
            f_SD_Ep = slowing_down_function(E_inj, p_inj, E_SD_array, p_SD_array, n_e, T_e, FI_species, regularization_thermal_ion_species, thermal_ion_densities_at_rho, thermal_ion_temperatures_at_rho; dampen=true, damp_type=:erfc, sigma=5.0)
            f_SD_VEL, vpara_SD_array, vperp_SD_array = Ep2VparaVperp(E_SD_array, p_SD_array, f_SD_Ep; my_gcp=getGCP(FI_species)(0.0,0.0,0.0,0.0), needJac=true, returnAbscissas=true)
            local nodes; nodes = (vpara_SD_array, vperp_SD_array)
            local itp; itp = Interpolations.interpolate(nodes,f_SD_VEL,Gridded(Linear()))
            local etp; etp = Interpolations.extrapolate(itp,Interpolations.Flat())
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
else
    verbose && println("No collision physics regularization included.")
end

append!(timestamps,time()) # The timestamp when the collision physics basis functions have been computed
dictionary_of_sections[prnt-1] = ("Computing collision physics basis functions (if requested)",diff(timestamps)[end])
###########################################################################################################
# SECTION: ENFORCE THE NOISE FLOOR
println("------------------------------------------------ SECTION $(prnt) ------------------------------------------------"); prnt += 1

verbose && println("Enforcing a noise floor of $(noise_floor_factor). All error values smaller than $(noise_floor_factor)*maximum(S) will be set to $(noise_floor_factor)*maximum(S)... ")
for (i,s) in enumerate(S)
    noise_floor_i = maximum(s)*noise_floor_factor
    verbose && print("---> Diagnostic $(i) has a noise floor of $(noise_floor_i)... ")
    below_noise_floor_inds = findall(x-> x<noise_floor_i, S_errs[i])
    if isempty(below_noise_floor_inds)
        verbose && println("no error values affected.")
    else
        println("")
        verbose && println("------> Found $(length(below_noise_floor_inds)) ($(round(100*length(below_noise_floor_inds)/length(S_errs[i]),digits=1))% of all) error values smaller than $(noise_floor_i). Re-setting... ")
        S_errs[i][below_noise_floor_inds] .= noise_floor_i
    end
end

append!(timestamps,time()) # The timestamp when the noise floor has been enforced
dictionary_of_sections[prnt-1] = ("Enforcing a noise floor on the measurement uncertainties",diff(timestamps)[end])
###########################################################################################################
# SECTION: EXCLUDE UNWANTED MEASUREMENT BINS
println("------------------------------------------------ SECTION $(prnt) ------------------------------------------------"); prnt += 1

if exclude_zero_measurements
    verbose && println("Performing exclusion of measurements equal to zero... ")
    for i in eachindex(S)
        verbose && println("---> For diagnostic $(i)... ")
        W_i = W[i]
        S_i = S[i]
        S_errs_i = S_errs[i]
        S_abscissas_i = S_abscissas[i]
        zero_inds = findall(x-> x==0.0, S_i)
        if !isempty(zero_inds)
            verbose && println("------> Found $(length(zero_inds)) ($(round(100*length(zero_inds)/length(S_i),digits=1))% of all) measurements equal to 0. Excluding... ")
            good_inds_i = filter(xx -> !(xx in zero_inds), collect(1:length(S_abscissas_i))) # The indices of all measurement bins to keep
            S[i] = S_i[good_inds_i] # Update measurements
            S_errs[i] = S_errs_i[good_inds_i] # Update errors/uncertainties
            S_abscissas[i] = S_abscissas_i[good_inds_i] # Update measurement abscissa
            W[i] = W_i[good_inds_i,:] # Update weight function matrix
        end
    end
end

if !(length(excluded_measurement_intervals)==length(excluded_measurement_units)==length(filepaths_S)) && !isempty(excluded_measurement_units) && !isempty(excluded_measurement_intervals)
    @warn "The lengths of input variables 'excluded_measurement_intervals' and 'excluded_measurement_units' were not equal to 'filepaths_S'. Cannot execute exclusion of unwanted measurement bins. All measurement bins will be included for all diagnostics."
else
    if !isempty(excluded_measurement_units) && !isempty(excluded_measurement_intervals)
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
                bad_inds_ii = findall(xx-> (xx>=interval_min && xx<=interval_max),S_abscissas_i)
                verbose && println("------> Identified $(length(bad_inds_ii)) measurements in excluded region (min,max)=($(interval_min),$(interval_max)) $(S_abscissas_units_i)") 
                append!(bad_inds_i, bad_inds_ii)
            end
            bad_inds_i = unique(bad_inds_i) # The indices of the excluded measurement intervals
            good_inds_i = filter(xx -> !(xx in bad_inds_i), collect(1:length(S_abscissas_i))) # The indices of all measurement bins to keep
            verbose && println("---> A total of $(length(bad_inds_i)) measurement bin centers will be excluded for diagnostic $(i)...")
            S[i] = S_i[good_inds_i] # Update measurements
            S_errs[i] = S_errs_i[good_inds_i] # Update errors/uncertainties
            S_abscissas[i] = S_abscissas_i[good_inds_i] # Update measurement abscissa
            W[i] = W_i[good_inds_i,:] # Update weight function matrix
        end
    else
        verbose && println("Input variables 'excluded_measurement_intervals' and 'excluded_measurement_units' not specified (empty). All measurements will be included when solving the inverse problem.")
    end
end

append!(timestamps,time()) # The timestamp when unwanted measurements have been excluded
dictionary_of_sections[prnt-1] = ("Excluding unwanted measurement bins",diff(timestamps)[end])
###########################################################################################################
# SECTION: INITIALIZE REGULARIZATION AND PRIOR SKELETON (EMPTY ARRAYS). IF REQUESTED, INCLUDE REGULARIZATION 
# MATRICES. THESE INCLUDE (COLLISION PHYSICS INCLUDED IN SECTION ABOVE):
#   - 0th ORDER TIKHONOV
#   - 1st ORDER TIKHONOV
#   - ICRF PHYSICS INFORMATION
# THIS SECTION WILL BE CHANGED IN THE FUTURE IF MORE FORMS OF REGULARIZATIONS AND PRIORS ARE TO BE INCLUDED
println("------------------------------------------------ SECTION $(prnt) ------------------------------------------------"); prnt += 1

rec_space_abscissas = W_abscissas[1][2:end] # The abscissas of the space of the solution to the inverse problem. All w abscissas are the same (because of sections 4 and 5), hence [1]. First w abscissa is the diagnostic measurement bins, hence [2:end].
rec_space_abscissas_units = W_abscissas_units[1][2:end] # The units of the rec_space_abscissas
rec_space_size = Tuple(length.(rec_space_abscissas)) # The size of the reconstruction space, e.g. (33,52,12).
rec_space_coords = CartesianIndices(rec_space_size) # All reconstruction space coordinates. E.g. [(1,1,1), (1,1,2), ..., (N1,N2,N3)] where Ni is the number of grid points in the i:th reconstruction space dimension
verbose && println("Creating list of reconstruction space coordinates (length=$(length(rec_space_coords))=$(reduce(*,map(x-> "x$(x)",rec_space_size))[2:end]))... ")

verbose && println("Creating inverse problem regularization and prior information skeleton (empty arrays)... ")
L = []
priors = []
priors_name = []

if "zerotikhonov" in lowercase.(String.(regularization))
    verbose && println(":ZEROTIKHONOV included in 'regularization' input variable. Adding 0th order Tikhonov regularization to the inverse solving algorithm... ")
    append!(L, [sparse(1.0(SparseArrays.I),length(rec_space_coords),length(rec_space_coords))]) # Sparse identity matrix
    append!(priors, [zeros(length(rec_space_coords))])
    append!(priors_name, ["0th order Tikhonov"])
end

first_tikhonov_reg = false
if "firsttikhonov" in lowercase.(String.(regularization))
    verbose && println(":FIRSTTIKHONOV included in 'regularization' input variable. Adding 1st order Tikhonov regularization to the inverse solving algorithm... ")
    first_tikhonov_reg = true
    L1 = forward_difference_matrix(rec_space_size) # The forward difference matrix of size (length(rec_space_size) * reduce(*,rec_space_size), reduce(*,rec_space_size))
    append!(L,[L1])
    append!(priors, [zeros(size(L1,1))])
    append!(priors_name, ["1st order Tikhonov"])
end

icrf_physics_reg = false
if "icrf" in lowercase.(String.(regularization))
    verbose && print(":ICRF included in 'regularization' input variable. Checking reconstruction space compatibility... ")
    # is_energy_pitch, is_vpara_vperp etc are functions in OWCF/extra/dependencies.jl
    Ep_bool, Ep_E_ind, Ep_p_ind = is_energy_pitch(rec_space_abscissas, rec_space_abscissas_units; returnExtra=true)
    VEL_bool, VEL_vpara_ind, VEL_vperp_ind = is_vpara_vperp(rec_space_abscissas, rec_space_abscissas_units; returnExtra=true)
    if length(rec_space_abscissas)==3 # (E,mu,Pphi) or (E,Lambda,Pphi_n)
        COM_noSigma_bool, COM_noSigma_E_ind, COM_noSigma_mu_ind, COM_noSigma_Pphi_ind = is_COM(rec_space_abscissas, rec_space_abscissas_units; returnExtra=true)
        normalized_COM_noSigma_bool, normalized_COM_noSigma_E_ind, normalized_COM_noSigma_Lambda_ind, normalized_COM_noSigma_Pphi_n_ind = is_normalized_COM(rec_space_abscissas, rec_space_abscissas_units; returnExtra=true)
    elseif length(rec_space_abscissas)==4 # (E,mu,Pphi,sigma) or (E,Lambda,Pphi_n,sigma)
        COM_bool, COM_E_ind, COM_mu_ind, COM_Pphi_ind, COM_sigma_ind = is_COM(rec_space_abscissas,rec_space_abscissas_units; returnExtra=true)
        normalized_COM_bool, normalized_COM_E_ind, normalized_COM_Lambda_ind, normalized_COM_Pphi_n_ind, normalized_COM_sigma_ind = is_normalized_COM(rec_space_abscissas,rec_space_abscissas_units; returnExtra=true)
    else # ADD MORE CHECKS HERE IN FUTURE VERSIONS, IF NEEDED <---------------------------------------------------------------------
    end
    COMPATIBILITY_ARRAY = [Ep_bool, VEL_bool, COM_noSigma_bool, normalized_COM_noSigma_bool, COM_bool, normalized_COM_bool] # An array of true's and false's
    COMPATIBLE_SPACE_INDICES = findall(x-> x,COMPATIBILITY_ARRAY) # Get the index pointing to which function in COMPATIBLE_OPTIONS returned true (if any of them)
    if isempty(COMPATIBLE_SPACE_INDICES)
        error(":ICRF was specified in 'regularization' input variable, but reconstruction space is not supported. Currently supported options include $(OWCF_ICRF_STREAMLINES_SUPPORTED_COORD_SYSTEMS). Please correct and re-try.")
    end
    if length(COMPATIBLE_SPACE_INDICES)>1 # If more than one of the functions in COMPATIBLE_OPTIONS returned true... 
        error("Compatibility check malfunction. Please post an issue with a screenshot of this error message at https://github.com/juliaFusion/OWCF/issues or try contacting henrikj@dtu.dk or anvalen@dtu.dk.")
    end
    RECONSTRUCTION_SPACE_ID = COMPATIBLE_SPACE_INDICES[1] # The identification number of the reconstruction space (see below)
    verbose && println("ok!")

    verbose && print("---> Loading magnetic equilibrium from $(regularization_equil_filepath)... ")
    verbose && println("Loading magnetic equilibrium... ")
    M, wall = nothing, nothing # Initialize global magnetic equilibrium variables
    try
        global M; global wall # Declare global scope
        M, wall = read_geqdsk(regularization_equil_filepath,clockwise_phi=false) # Assume the phi-direction is pointing counter-clockwise when tokamak is viewed from above. This is true for almost all coordinate systems used in the field of plasma physics
    catch # Otherwise, assume magnetic equilibrium is a saved .jld2 file
        global M; global wall; local myfile
        myfile = jldopen(regularization_equil_filepath,false,false,false,IOStream)
        M = myfile["S"]
        wall = myfile["wall"]
        close(myfile)
    end
    psi_axis, psi_bdry = psi_limits(M)
    verbose && println("ok!")

    icrf_physics_reg = true
    if !first_tikhonov_reg
        verbose && println("Computing finite difference matrix to be used in :ICRF regularization... ")
        L1 = forward_difference_matrix(rec_space_size) # The forward difference matrix of size (length(rec_space_size) * reduce(*,rec_space_size), reduce(*,rec_space_size))
    end
    
    if RECONSTRUCTION_SPACE_ID==1 # (E,p)
        verbose && print("---> Computing ICRF regularization matrix for (assumed) (E,p) reconstruction space (or permutations thereof)... ")
        L_ICRF = icrf_regularization_matrix(M, rec_space_abscissas, FI_species, regularization_wave_frequency, regularization_cyclotron_harmonic; 
                                            toroidal_mode_number=regularization_toroidal_mode_number, coord_system="(E,p)", coord_system_order=(Ep_E_ind,Ep_p_ind), 
                                            R_of_interest=units_conversion_factor(R_of_interest_units,"m")*R_of_interest, 
                                            z_of_interest=units_conversion_factor(z_of_interest_units,"m")*z_of_interest, L1=L1)
    elseif RECONSTRUCTION_SPACE_ID==2 # (vpara,vperp)
        verbose && print("---> Computing ICRF regularization matrix for (assumed) (vpara,vperp) reconstruction space (or permutations thereof)... ")
        L_ICRF = icrf_regularization_matrix(M, rec_space_abscissas, FI_species, regularization_wave_frequency, regularization_cyclotron_harmonic; 
                                            toroidal_mode_number=regularization_toroidal_mode_number, coord_system="(vpara,vperp)", coord_system_order=(VEL_vpara_ind,VEL_vperp_ind), 
                                            R_of_interest=units_conversion_factor(R_of_interest_units,"m")*R_of_interest, 
                                            z_of_interest=units_conversion_factor(z_of_interest_units,"m")*z_of_interest, L1=L1)
    elseif RECONSTRUCTION_SPACE_ID==3 # (E,mu,Pphi)
        verbose && print("---> Computing ICRF regularization matrix for (assumed) (E,mu,Pphi) reconstruction space (or permutations thereof)... ")
        L_ICRF = icrf_regularization_matrix(M, rec_space_abscissas, FI_species, regularization_wave_frequency, regularization_cyclotron_harmonic; 
                                            toroidal_mode_number=regularization_toroidal_mode_number, coord_system="(E,mu,Pphi)", 
                                            coord_system_order=(COM_noSigma_E_ind,COM_noSigma_mu_ind,COM_noSigma_Pphi_ind), L1=L1)
    elseif RECONSTRUCTION_SPACE_ID==4 # (E,Lambda,Pphi_n)
        verbose && print("---> Computing ICRF regularization matrix for (assumed) (E,Lambda,Pphi_n) reconstruction space (or permutations thereof)... ")
        L_ICRF = icrf_regularization_matrix(M, rec_space_abscissas, FI_species, regularization_wave_frequency, regularization_cyclotron_harmonic; 
                                            toroidal_mode_number=regularization_toroidal_mode_number, coord_system="(E,Lambda,Pphi_n)", 
                                            coord_system_order=(normalized_COM_noSigma_E_ind,normalized_COM_noSigma_Lambda_ind,normalized_COM_noSigma_Pphi_n_ind), L1=L1)
    elseif RECONSTRUCTION_SPACE_ID==5 # (E,mu,Pphi,sigma)
        verbose && print("---> Computing ICRF regularization matrix for (assumed) (E,mu,Pphi,sigma) reconstruction space (or permutations thereof)... ")
        L_ICRF = icrf_regularization_matrix(M, rec_space_abscissas, FI_species, regularization_wave_frequency, regularization_cyclotron_harmonic; 
                                            toroidal_mode_number=regularization_toroidal_mode_number, coord_system="(E,mu,Pphi,sigma)", 
                                            coord_system_order=(COM_E_ind,COM_mu_ind,COM_Pphi_ind,COM_sigma_ind), L1=L1)
    elseif RECONSTRUCTION_SPACE_ID==6 # (E,Lambda,Pphi_n,sigma)
        verbose && print("---> Computing ICRF regularization matrix for (assumed) (E,Lambda,Pphi_n,sigma) reconstruction space (or permutations thereof)... ")
        L_ICRF = icrf_regularization_matrix(M, rec_space_abscissas, FI_species, regularization_wave_frequency, regularization_cyclotron_harmonic; 
                                            toroidal_mode_number=regularization_toroidal_mode_number, coord_system="(E,Lambda,Pphi_n,sigma)", 
                                            coord_system_order=(normalized_COM_E_ind,normalized_COM_Lambda_ind,normalized_COM_Pphi_n_ind,normalized_COM_sigma_ind), L1=L1)
    else
        error("This error print should be impossible to reach. Please post an issue with a screenshot of this error message at https://github.com/juliaFusion/owCF/issues or try contacting henrikj@dtu.dk or anvalen@dtu.dk.")
    end
    verbose && println("ok!")

    append!(L,[L_ICRF])
    append!(priors, [zeros(size(L_ICRF,1))])
    append!(priors_name, ["ICRF"])
end

# The check below is only relevant in future versions where a general prior can be used!
!(sum(diff(length.(priors)))==0) && error("The number of reconstruction space points should be equal for all priors, i.e. [N,...,N]. Instead, $(length.(priors)) was found. Please correct and re-try.")

append!(timestamps,time()) # The timestamp when the regularization matrices have been initialized
dictionary_of_sections[prnt-1] = ("Initializing regularization matrices (if requested)",diff(timestamps)[end])
###########################################################################################################
# SECTION: EXCLUDE ALL RECONSTRUCTION SPACE POINTS FOR WHICH THE WEIGHT MATRICES OF ALL DIAGNOSTICS ARE ZERO
# PLEASE NOTE! THIS ALSO TAKES CARE OF INVALID ORBITS, SHOULD THE RECONSTRUCTION SPACE BE E.G. ORBIT SPACE 
# OR CONSTANTS-OF-MOTION SPACE.
println("------------------------------------------------ SECTION $(prnt) ------------------------------------------------"); prnt += 1

verbose && println("Finding reconstruction space points for which the weight matrix (column) is zero (no sensitivity), for all diagnostics... ")
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
if !isempty(rec_space_zero_dict)
    bad_inds = []
    for key in keys(rec_space_zero_dict)
        if rec_space_zero_dict[key]==length(W) # If the reconstruction space point has zero sensitivity, for all diagnostics
            verbose && println("---> The reconstruction space point (bad point) $([W_abscissas[1][i+1][key[i]] for i in 1:length(key)]) with units $([W_abscissas_units[1][i+1] for i in 1:length(key)]) has no sensitivity for all diagnostics. Excluding from reconstruction... ")
            append!(bad_inds, findfirst(x-> x==key, vec(rec_space_coords)))
        end
    end
    bad_inds = unique(bad_inds) # Should be redundant, but just in case
    good_inds = filter(x-> !(x in bad_inds), collect(1:length(rec_space_coords)))
    
    rec_space_coords = rec_space_coords[good_inds] # Update the list of reconstruction space coordinates, to be able to correctly inflate tomographic reconstruction vector to full N-dimensional array
    for (i,w) in enumerate(W)
        verbose && println("------> Excluding $(length(bad_inds)) bad points from the weight matrix of diagnostic $(i)... ")
        W[i] = w[:,good_inds]
    end
    if "collisions" in lowercase.(String.(regularization)) # Exclude bad points from slowing-down (collision physics) basis functions
        verbose && println("------> Excluding $(length(bad_inds)) bad points from the slowing-down (collision physics) basis functions... ")
        F_SD_safe = F_SD_safe[good_inds,:]
    end

    if !isempty(priors) # If not empty, exclude bad points from the regularization matrices and priors
        verbose && println("------> Excluding $(length(bad_inds)) bad points from the regularization matrices and priors... ")
        for (i,l) in enumerate(L)
            L[i] = l[good_inds,good_inds]
            priors[i] = (priors[i])[good_inds]
        end
    end
else
    verbose && println("---> No such reconstruction space points found.")
end

append!(timestamps,time()) # The timestamp when reconstruction space points with zero sensitivity for all diagnostics have been excluded from the problem
dictionary_of_sections[prnt-1] = ("Excluding reconstruction space points with zero sensitivity",diff(timestamps)[end])
###########################################################################################################
# SECTION: INCLUDE COLLISIONAL PHYSICS BY MULTIPLYING THE WEIGHT MATRICES WITH THE SLOWING-DOWN BASIS FUNCTIONS
# MODIFY THE OTHER REGULARIZATION MATRICES TO MATCH THE NUMBER OF SLOWING-DOWN BASIS FUNCTIONS 
println("------------------------------------------------ SECTION $(prnt) ------------------------------------------------"); prnt += 1

if "collisions" in lowercase.(String.(regularization)) # Don't need all the safety checks here, it was already taken care of in section 9
    verbose && println("Multiplying weight matrices with slowing-down (collision physics) matrix to include collision physics as prior information/regularization... ")
    for (i,w) in enumerate(W)
        verbose && print("---> For diagnostic $(i)... ")
        W[i] = w*F_SD_safe
        verbose && println("ok!")
    end

    verbose && println("Modifying other regularization matrices to work with collision physics regularization... ")
    for ir in eachindex(L)
        if priors_name[ir]=="1st order Tikhonov"
            L[ir] = forward_difference_matrix((size(F_SD_safe,2),))
            priors[ir] = zeros(size(F_SD_safe,2))
        end
        if priors_name[ir]=="0th order Tikhonov"
            L[ir] = diagm(ones(size(F_SD_safe,2)))
            priors[ir] = zeros(size(F_SD_safe,2))
        end
        if priors_name[ir]=="ICRF"
            @warn ":ICRF regularization is currently not supported together with :COLLISIONS regularization. Using 1st order Tikhonov regularization instead"
            L[ir] = forward_difference_matrix((size(F_SD_safe,2),))
            priors[ir] = zeros(size(F_SD_safe,2))
        end
    end
else
    verbose && println("No collision physics basis functions included.")
end

append!(timestamps,time()) # The timestamp when the weight matrices have been transformed using collision physics
dictionary_of_sections[prnt-1] = ("Multiplying weight matrices by collision physics basis function (if requested)",diff(timestamps)[end])
###########################################################################################################
# SECTION: RESCALE WEIGHT FUNCTIONS
println("------------------------------------------------ SECTION $(prnt) ------------------------------------------------"); prnt += 1

if rescale_W
    verbose && println("Rescaling weight matrices using computed rescale factors... ")
    for (i,r) in enumerate(rescale_W_factors)
        verbose && print("---> Rescaling weight matrix for diagnostic $(i) with a factor of $(round(r,sigdigits=4))... ")
        W[i] = r .*W[i]
        verbose && println("ok!")
    end
else
    verbose && println("No rescaling of weight matrices.")
end

append!(timestamps,time()) # The timestamp when the weight matrices have been rescaled
dictionary_of_sections[prnt-1] = ("Rescaling weight matrices (if requested)",diff(timestamps)[end])
###########################################################################################################
# SECTION: NORMALIZE ALL WEIGHT MATRICES AND MEASUREMENTS BY THE MEASUREMENT ERRORS
println("------------------------------------------------ SECTION $(prnt) ------------------------------------------------"); prnt += 1

verbose && println("Normalizing all weight matrices and measurements by the measurement uncertainties... ")
W_hat = Vector{Matrix{Float64}}(undef,length(filepaths_W))
S_hat = Vector{Vector{Float64}}(undef,length(filepaths_S))
for (i,w) in enumerate(W)
    errs_i = S_errs[i]
    numOZeroErrs = findall(x-> x<=0.0, errs_i)
    if !isempty(numOZeroErrs)
        if !allow_measurements_with_zero_uncertainty
            error("Found $(length(numOZeroErrs)) measurements with zero uncertainty for diagnostic $(i). This is not allowed (input variable 'allow_measurements_with_zero_uncertainty' was set to false). Set input variable 'allow_measurements_with_zero_uncertainty' to true, or correct and re-try")
        end
        verbose && println("---> Found $(length(numOZeroErrs)) measurements with 0 uncertainty for diagnostic $(i). Giving maximum importance... ")
        minNonZeroErr = Inf
        for err in errs_i
            if 0<err<minNonZeroErr
                minNonZeroErr = err
            end
        end
        if !isfinite(minNonZeroErr)
            verbose && println("------> No non-zero uncertainty found for diagnostic $(i). All measurements will be given equal importance... ")
            minNonZeroErr = 1.0 # Use 1.0 to leave measurements unaffected
        else
            verbose && println("------> The minimum non-zero uncertainty is $(minNonZeroErr). The following measurements had 0 uncertainty and will be given that uncertainty value instead: ")
            if verbose
                for i_err in findall(x-> x<=0.0, errs_i)
                    println("---------> $(errs_i[i_err]) $(S_errs_units[i])")
                end
            end
        end
        errs_i[findall(x-> x<=0.0, errs_i)] .= minNonZeroErr
    end
    W_hat[i] = w ./ errs_i
    S_hat[i] = S[i] ./ errs_i
end

append!(timestamps,time()) # The timestamp when the weight matrices and measurements have been normalized
dictionary_of_sections[prnt-1] = ("Normalize measurements by uncertainties",diff(timestamps)[end])
###########################################################################################################
# SECTION: CONCATENATE ALL MEASUREMENTS AND BUILD TOTAL WEIGHT MATRIX
println("------------------------------------------------ SECTION $(prnt) ------------------------------------------------"); prnt += 1

verbose && println("Concatenating all diagnostic measurements into one long S vector... ")
S_hat_long = vcat(S_hat...)

verbose && println("Concatenating all weight matrices into one big W matrix... ")
W_hat_long = vcat(W_hat...)

append!(timestamps,time()) # The timestamp when the total weight matrices and measurement vector have been built
dictionary_of_sections[prnt-1] = ("Concatenating all weight matrices and measurements",diff(timestamps)[end])
###########################################################################################################
# SECTION: NORMALIZE WEIGHT MATRIX, SIGNAL AND REGULARIZATION MATRICES TO ORDER UNITY
println("------------------------------------------------ SECTION $(prnt) ------------------------------------------------"); prnt += 1

verbose && println("Normalizing weight matrix W and measurements S to order unity... ")
Whm = maximum(W_hat_long)
Shm = maximum(S_hat_long)
W_hh = W_hat_long ./Whm
S_hh = S_hat_long ./Shm

verbose && println("Normalizing regularization matrices L to order unity... ")
for (il, l) in enumerate(L)
    L[il] = l ./maximum(l)
end

append!(timestamps,time()) # The timestamp when the weight matrices and signals have been normalized to order unity
dictionary_of_sections[prnt-1] = ("Normalizing weight matrix and measurements to order unity",diff(timestamps)[end])
###########################################################################################################
# SECTION: SOLVE THE INVERSE PROBLEM
println("------------------------------------------------ SECTION $(prnt) ------------------------------------------------"); prnt += 1

if !(length(nr_array)==length(regularization))
    @warn "length(nr_array) ($(length(nr_array))) was not equal to length(regularization) ($(length(regularization))). Elements in 'nr_array' can therefore not be related to the elements in 'regularization'. A default number of 20 grid points will therefore be used for all regularization types ($(regularization))"
    nr_array = Int64.(20 .*ones(length(regularization)))
end
if !(length(regularization_loglambda_mins)==length(regularization))
    @warn "length(regularization_loglambda_mins) ($(length(regularization_loglambda_mins))) was not equal to length(regularization) ($(length(regularization))). Elements in 'regularization_loglambda_mins' can therefore not be related to the elements in 'regularization'. A default lower bound of log10(λ)=-3 will therefore be used for all regularization types ($(regularization))"
    regularization_loglambda_mins = (-3) .*ones(length(regularization))
end
if !(length(regularization_loglambda_maxs)==length(regularization))
    @warn "length(regularization_loglambda_maxs) ($(length(regularization_loglambda_maxs))) was not equal to length(regularization) ($(length(regularization))). Elements in 'regularization_loglambda_maxs' can therefore not be related to the elements in 'regularization'. A default upper bound of log10(λ)=3 will therefore be used for all regularization types ($(regularization))"
    regularization_loglambda_maxs = 3 .*ones(length(regularization))
end

# nr_array[i] should specify the number of λ grid points for regularization[i]
# regularization_loglambda_mins[i] should specify the lower bound of log10(λ) for regularization[i]
# regularization_loglambda_maxs[i] should specify the upper bound of log10(λ) for regularization[i]
# The spacing between λ grid points is logarithmic
nop = length(priors)
s = nop>1 ? "s" : ""
verbose && println("Preparing $(nop) regularization strength (hyper-parameter) logarithmic range$(s)... ")
hyper_arrays = []
for ir in eachindex(priors) # For all regularizations/priors that require a hyper-parameter...
    nr = 0
    loglambda_min = 0
    loglambda_max = 0
    if priors_name[ir]=="1st order Tikhonov"
        iir = findfirst(x-> x=="firsttikhonov", lowercase.(String.(regularization)))
        nr = nr_array[iir]
        loglambda_min = regularization_loglambda_mins[iir]
        loglambda_max = regularization_loglambda_maxs[iir]
    end
    if priors_name[ir]=="0th order Tikhonov"
        iir = findfirst(x-> x=="zerotikhonov", lowercase.(String.(regularization)))
        nr = nr_array[iir]
        loglambda_min = regularization_loglambda_mins[iir]
        loglambda_max = regularization_loglambda_maxs[iir]
    end
    if priors_name[ir]=="ICRF"
        iir = findfirst(x-> x=="icrf", lowercase.(String.(regularization)))
        nr = nr_array[iir]
        loglambda_min = regularization_loglambda_mins[iir]
        loglambda_max = regularization_loglambda_maxs[iir]
    end
    append!(hyper_arrays,[10 .^(range(Float64(loglambda_min),stop=Float64(loglambda_max),length=Int64(nr)))])
    verbose && println("---> $(priors_name[ir]) with $(nr) λ grid points between 10^$(loglambda_min) and 10^$(loglambda_max)... ")
end
hyper_points = Iterators.product(hyper_arrays...)

verbose && println("Concatenating all priors into one long prior vector... ")
prior = vcat(priors...)

verbose && println("Pre-allocating solution arrays and L-curve quantities... ")
sols = zeros(size(W_hh,2), length(hyper_points))
sols_rawnorm = zeros(length(hyper_points))
sols_datafit = zeros(length(hyper_points))

nonneg = false
if "nonneg" in lowercase.(String.(constraints))
    verbose && println(":NONNEG found in input variable 'constraints'. Inverse problem will be solved with non-negativity constraint... ")
    nonneg = true
end

collision_physics_reg = false
if "collisions" in lowercase.(String.(regularization))
    collision_physics_reg = true # Just a bool to speed up computations
end
verbose && println("")

if verbose
    verbose_onlyForSolve = true # If verbose is true, so is verbose_onlyForSolve
end

nip = length(hyper_points); s = nip>1 ? "s" : ""
verbose_onlyForSolve && println("Solving $(nip) inverse problem$(s) (each with $(size(W_hh,2)) unknowns and $(size(W_hh,1)) measurements)... ")
for (i,hp) in enumerate(hyper_points)
    verbose_onlyForSolve && println("---> ($(i)/$(nip)) for regularization hyper-parameter point $(round.(vcat(hp...),sigdigits=4)) ($([pn for pn in priors_name]))... ")
    L_i = Matrix(vcat((vcat(hp...) .*L)...)) # Multiply hyper-parameter values with the corresponding regularization matrices. Then, concatenate into one big regularization matrix
    x = Convex.Variable(size(W_hh,2))

    # In future versions, this will be coded more elegantly
    if !nonneg && !collision_physics_reg # Case 1
        problem = Convex.minimize(Convex.sumsquares(vcat(W_hh,L_i) * x - vcat(S_hh,prior)), [])
    end
    if nonneg && !collision_physics_reg # Case 2
        problem = Convex.minimize(Convex.sumsquares(vcat(W_hh,L_i) * x - vcat(S_hh,prior)), [x >= 0])
    end
    if !nonneg && collision_physics_reg # Case 3 (identical to Case 1)
        problem = Convex.minimize(Convex.sumsquares(vcat(W_hh,L_i) * x - vcat(S_hh,prior)), [])
    end
    if nonneg && collision_physics_reg # Case 4
        problem = Convex.minimize(Convex.sumsquares(vcat(W_hh,L_i) * x - vcat(S_hh,prior)), [F_SD_safe*x >= 0])
    end

    solve!(problem, SCS.Optimizer; silent=true)

    sol_i = dropdims(x.value,dims=2) # One redundant dimension is returned by the solve!() method

    sols[:,i] = (Shm/Whm) .*sol_i
    sols_rawnorm[i] = norm(Matrix(vcat(L...))*sol_i) # ||Lx||
    sols_datafit[i] = norm(W_hh*sol_i - S_hh) # ||WF* - S||
end

# Collect hyper_points iterators, for element indexing purposes in plot and save sections
hyper_points = collect(hyper_points)

append!(timestamps,time()) # The timestamp when the inverse problem has been solved
dictionary_of_sections[prnt-1] = ("Solving the inverse problem",diff(timestamps)[end])
###########################################################################################################
# SECTION: DETERMINE INDEX OF BALANCED SOLUTION
println("------------------------------------------------ SECTION $(prnt) ------------------------------------------------"); prnt += 1

verbose && println("Determining index of balanced solution... ")
verbose && println("---> by computing the point of maximum curvature on the L-curve (default)... ")
sols_datafit_itp = Interpolations.interpolate(sols_datafit, BSpline(Cubic())) # Create a Cubic spline interpolation, so that you can...
sols_rawnorm_itp = Interpolations.interpolate(sols_rawnorm, BSpline(Cubic())) # -||-
sols_datafit_prim = map(x-> Vector(Interpolations.gradient(sols_datafit_itp,x))[1],eachindex(sols_datafit)) # ...easily compute the derivative
sols_rawnorm_prim = map(x-> Vector(Interpolations.gradient(sols_rawnorm_itp,x))[1],eachindex(sols_rawnorm)) # -||-
sols_datafit_prim_itp = Interpolations.interpolate(sols_datafit_prim, BSpline(Cubic())) # Create a Cubic spline interpolation of the derivative, so that you can...
sols_rawnorm_prim_itp = Interpolations.interpolate(sols_rawnorm_prim, BSpline(Cubic())) # -||-
sols_datafit_primprim = map(x-> Vector(Interpolations.gradient(sols_datafit_prim_itp,x))[1],eachindex(sols_datafit)) # ...easily compute the second derivative
sols_rawnorm_primprim = map(x-> Vector(Interpolations.gradient(sols_rawnorm_prim_itp,x))[1],eachindex(sols_rawnorm)) # -||-

sdf = sols_datafit; sdfp = sols_datafit_prim; sdfpp = sols_datafit_primprim # Create abbreviations for easy overview of for-loop below
srn = sols_rawnorm; srnp = sols_rawnorm_prim; srnpp = sols_rawnorm_primprim # -||-

gamma = zeros(size(sols,2)) # The curvature of the L-curve, pre-allocated
for i in eachindex(gamma)
    gamma[i] = abs(srnpp[i]*sdfp[i] - srnp[i]*sdfpp[i])/((srnp[i]*srnp[i] + sdfp[i]*sdfp[i])^(3/2)) # Curvature formula
end
ilm = argmax(gamma) # Index of L-curve maximum curvature, i.e. the default index of the balanced solution
ilm_default = true # Start by assuming that the computed ilm value is ok

# Check for non-robustness in gamma
diff_gamma = diff(gamma[1:ilm]) # The differences of all curvatures up to the maximum curvature index
diff_gamma_neg_inds = findall(x-> x<0, diff_gamma) # Find all negative differences
if !isempty(diff_gamma_neg_inds) || (ilm in collect(1:Int64(round(0.1*length(gamma))))) || (ilm in collect(Int64(round(0.9*length(gamma))):length(gamma))) # If there are any negative differences, or if ilm is in the first (last) 10% of gamma indices
    gamma_i_mean = sum(eachindex(gamma) .*(inv(sum(gamma)) .*gamma)) # Compute the mean index of gamma
    if !((ilm-1) < gamma_i_mean < (ilm+1)) # If the mean index is NOT within ilm +/- 1
        verbose && println("---> Signs of non-robustness found when computing L-curve curvature! Using minimum of abs(||S-WF|| - ||L*F||) instead... ")
        # The curvature curve is likely unstable and the maximum curvature computation is not robust enough to identify 'ilm'
        # Therefore, use the index of the point of minium P = abs(sols_datafit - sols_rawnorm) instead
        sdf_n = sdf ./maximum(sdf) # Normalize quantities to be able to identify usable intersection
        srn_n = srn ./maximum(srn) # Normalize quantities to be able to identify usable intersection
        ilm = argmin(abs.(sdf_n - srn_n))
        gamma = inv.(abs.(sdf_n - srn_n)) # Curvature proxy
        ilm_default = false # The ilm value was not computed in the default way (the point of L-curve maximum curvature), but via the minimum of abs(||S-WF|| - ||L*F||)
    end
end

append!(timestamps,time()) # The timestamp when the index of the balanced solution has been determined
dictionary_of_sections[prnt-1] = ("Determining index of balanced solution",diff(timestamps)[end])
###########################################################################################################
# SECTION: PLOT SOLUTIONS AND RELATED QUANTITIES, IF SPECIFIED
println("------------------------------------------------ SECTION $(prnt) ------------------------------------------------"); prnt += 1

date_and_time = split("$(Dates.now())","T")[1]*"at"*split("$(Dates.now())","T")[2][1:5] # Might be needed immidiately below, and definitely when saving output data
if plot_solutions || gif_solutions
    verbose && plot_solutions && println("Creating .png files of solutions... ")
    verbose && gif_solutions && println("Creating .gif file of solutions... ")

    # Pre-compute these quantities, needed for plots
    sdf0, sdfN = extrema(sdf); sdf_diff = abs(sdfN-sdf0) # sdf = sols_datafit ||S-WF|| (please see above)
    srn0, srnN = extrema(srn); srn_diff = abs(srnN-srn0) # srn = sols_rawnorm ||L*F|| (please see above)

    rec_space_dims = tuple(eachindex(rec_space_size)...) # The dimensions of the reconstruction space numbered, e.g. (1,2,3)
    rec_space_diff = map(x-> x[1],diff.(W_abscissas[1][2:end])) # The grid spacings of the reconstruction space, e.g. [5.0,0.01,0.04,0.05]. Assume equidistant grid

    verbose && println("---> Pre-computing plot quantities and layouts... ")
    if length(W)==1
        layout_WF = (1,1)
    elseif length(W)==2
        layout_WF = (1,2)
    elseif length(W)==3
        layout_WF = (1,3)
    elseif length(W)>3 && mod(length(W),3)==0
        layout_WF = (Int64(length(W)/3),3)
    elseif length(W)>3 && !(mod(length(W),3)==0)
        layout_WF = (ceil(length(W)/3),3)
    else
        error("This should be impossible to reach! (unless you have no diagnostics...?)")
    end
    if length(rec_space_size)==1 || length(rec_space_size)==2 # The 1- and 2-dimensional cases can be plotted in a single plot (no subplots)
        layout_F = (1,1)
    elseif length(rec_space_size)==3 # If the dimensionality of the reconstruction space is 3.. Currently, the two supported cases are (E,pm,Rm) and (E,mu,Phi). Neither requires a Jacobian when integrating
        layout_F = (1,3)
        plot_coords = [(1,2),(1,3),(2,3)] # Plot the 2-combinations of all indices, excluding duplicates and self-combinations
    elseif length(rec_space_size)==4 # If the dimensionality of the reconstruction space is 4... Currently, the only supported case is (E,p,R,z). More to be added in future versions
        layout_F = (2,3)
        if is_EpRz(W_abscissas[1], W_abscissas_units[1])
            dummy_bool, w_energy_ind, w_pitch_ind, w_R_ind, w_z_ind, R_array, z_array = is_EpRz(W_abscissas[1],W_abscissas_units[1]; returnExtra=true)
            w_energy_ind -= 1; w_pitch_ind -=1; w_R_ind -= 1; w_z_ind -=1 # All indices need to be adjusted by 1, since return of is_EpRz() is given relative to the elements of W_abscissas, whose first element always is the diagnostics measurement bins
            dE = rec_space_diff[w_energy_ind]; dp = rec_space_diff[w_pitch_ind]; dR = rec_space_diff[w_R_ind]; dz = rec_space_diff[w_z_ind]
            two_pi_R_jacobian = 2*pi .*reshape(R_array,tuple([(xx-> xx==w_R_ind ? length(R_array) : 1) for ir in 1:4]...)) # The R grid points, but reshaped as a 4D array with the R values placed in the correct dimension. All other dimensions have length 1
            plot_no_jac_coords = [(w_energy_ind, w_R_ind), (w_pitch_ind, w_R_ind), (w_R_ind, w_z_ind)] # These will not need the jacobian when plotting
            plot_jac_coords = [(w_energy_ind, w_pitch_ind),(w_energy_ind, w_z_ind),(w_pitch_ind, w_z_ind)] # These will need the jacobian when plotting (since we integrate over R)
        elseif false # <--- ADD MORE 4D COORDINATE CASES HERE IN FUTURE VERSIONS
        else
            @warn "Currently, only (E,p,R,z) coordinates are supported for plotting of 4D fast-ion distributions. An empty subfigure will be plotted instead."
            layout_F = (1,1)
        end
    elseif false # <--- HIGHER-DIMENSIONAL CASES WILL BE SUPPORTED IN FUTURE VERSIONS
    else
        @warn "Plotting of a 5D or 6D fast-ion distribution is currently not supported. An empty subfigure will be plotted instead."
        layout_F = (1,1)
    end

    verbose && println("---> Plotting solutions... ")
    lcurveproxystring = ilm_default ? "" : " (proxy)" # To show if the L-curve curvature plot is a proxy or not
    plot_inds = 1:length(hyper_points)
    anim = @animate for i=plot_inds
        sol_i = sols[:,i]
        hp = hyper_points[i]

        # CREATE THE L-CURVE SUBFIGURE, which will be at the top of the plot stack
        l_mag = sqrt(((sdf[i]-sdf0)^2)/(sdf_diff^2)+((srn[i]-srn0)^2)/(srn_diff^2))
        plt_L_curvature = Plots.plot(axis=([],false)) # Dummy plot, for white space
        plt_L_curvature = Plots.plot(gamma, label="", title="L-curve curvature$(lcurveproxystring)")
        plt_L_curvature = Plots.scatter!([ilm],[gamma[ilm]], label="Max. curvature")
        plt_L_curvature = Plots.scatter!([i], [gamma[i]], label="")

        plt_L_curve = Plots.plot([0,sdf[i]],[0,srn[i]],color=:black,label="")
        plt_L_curve = Plots.plot!(plt_L_curve, sdf, srn,xlims=(minimum(sdf)-0.05*sdf_diff, maximum(sdf)+0.05*sdf_diff),ylims=(minimum(srn)-0.05*srn_diff, maximum(srn)+0.05*srn_diff),label="")
        plt_L_curve = Plots.scatter!(plt_L_curve, [sdf[1]],[srn[1]],label="Start")
        plt_L_curve = Plots.scatter!(plt_L_curve, [sdf[end]],[srn[end]],label="End")
        plt_L_curve = Plots.scatter!(plt_L_curve, [sdf[ilm]],[srn[ilm]],label="Balanced F*")
        plt_L_curve = Plots.scatter!(plt_L_curve, [sdf[i]],[srn[i]],label="", title="$(i): log10(λ)=$(round.(log10.(vcat(hp...)),sigdigits=4))")
        plt_L_curve = Plots.plot!(xlabel="||WF-S||", ylabel="||L*F||")

        subfig_L_curve = Plots.plot(plt_L_curve, plt_L_curvature, layout=(1,2))

        # CREATE THE S vs WF SUBFIGURE, which will be in the middle of the plot stack
        WF_plots = Vector{Any}(undef,length(W))
        for (iw,w) in enumerate(W)
            WF = w*sol_i
            wf_plot = Plots.plot(S_abscissas[iw], WF, title="Diagnostic $(iw)",label="WF*", linewidth=2.5)
            wf_plot = Plots.scatter!(wf_plot, S_abscissas[iw], S[iw], label="S")
            wf_plot = Plots.yerror!(wf_plot, S_abscissas[iw], S[iw]; label="", yerror=S_errs[iw])
            wf_plot = Plots.plot!(wf_plot, xlabel=S_abscissas_units[iw], ylabel=S_units[iw])
            WF_plots[iw] = wf_plot
        end
        subfig_WF = Plots.plot(WF_plots..., layout=layout_WF)

        # CREATE THE F* SUBFIGURE, which will be at the bottom of the plot stack
        # If slowing-down basis function coefficients, compute fast-ion distribution
        if collision_physics_reg
            F_vec = F_SD_safe*sol_i
        else
            F_vec = sol_i
        end
        # Reshape into inflated fast-ion distribution
        F = zeros(rec_space_size) # Pre-allocate the inflated fast-ion distribution
        for (ic,c) in enumerate(rec_space_coords)
            F[c] = F_vec[ic]
        end
        if length(rec_space_size)==1 # If the dimensionality of the reconstruction space is 1
            F_plots = [Plots.plot(W_abscissas[1][2],F,xlabel=W_abscissas_units[1][2],ylabel="F* [$(units_inverse(W_abscissas_units[1][2]))]", title="F* [$(units_inverse(W_abscissas_units[1][2]))]",label="")]
        elseif length(rec_space_size)==2 # If the dimensionality of the reconstruction space is 2
            F_plots = [Plots.heatmap(W_abscissas[1][2],W_abscissas[1][3],transpose(F),xlabel=W_abscissas_units[1][2],ylabel=W_abscissas_units[1][3], title="F* [$(units_multiply(units_inverse(W_abscissas_units[1][2]), units_inverse(W_abscissas_units[1][3])))]",fillcolor=cgrad([:white, :darkblue, :green, :yellow, :orange, :red]))]
        elseif length(rec_space_size)==3 # -||- is 3. Currently, the two supported cases are (E,pm,Rm) and (E,mu,Phi). Neither requires a Jacobian when integrating
            F_plots = Vector{Any}(undef,length(plot_coords)) # Pre-allocate vector with plots
            for (ipc, plot_coord) in enumerate(plot_coords)
                ix = plot_coord[1]; x = W_abscissas[1][ix+1]; x_units = W_abscissas_units[1][ix+1] # +1 because W_abscissas and W_abscissas_units has a first dimension of diagnostic measurement bins
                iy = plot_coord[2]; y = W_abscissas[1][iy+1]; y_units = W_abscissas_units[1][iy+1] # +1 because W_abscissas and W_abscissas_units has a first dimension of diagnostic measurement bins
                iz = filter(xx-> !(xx==ix || xx==iy), rec_space_dims)
                F_sub = dropdims(sum(rec_space_diff[iz...] .*F,dims=iz),dims=iz)
                F_plot = Plots.heatmap(x,y,F_sub; xlabel=x_units, ylabel=y_units, title="F* [$(units_multiply(units_inverse(x_units), units_inverse(y_units)))]", fillcolor=cgrad([:white, :darkblue, :green, :yellow, :orange, :red]))
                F_plots[ipc] = F_plot
            end
        elseif length(rec_space_size)==4 # -||- is 4. Currently, the only supported case is (E,p,R,z). Suitable Jacobian is computed prior to plot for-loop
            if is_EpRz(W_abscissas[1], W_abscissas_units[1])
                F_plots = Vector{Any}(undef,length(plot_no_jac_coords)+length(plot_jac_coords)) # Pre-allocate vector with plots
                global ipc; ipc = 1 # To keep track of plots
                for plot_no_jac_coord in plot_no_jac_coords
                    global ipc
                    ix = plot_coord[1]; x = W_abscissas[1][ix+1]; x_units = W_abscissas_units[1][ix+1] # +1 because W_abscissas and W_abscissas_units has a first dimension of diagnostic measurement bins
                    iy = plot_coord[2]; y = W_abscissas[1][iy+1]; y_units = W_abscissas_units[1][iy+1] # +1 because W_abscissas and W_abscissas_units has a first dimension of diagnostic measurement bins
                    izs = filter(xx-> !(xx==ix || xx==iy), rec_space_dims)
                    F_sub = dropdims(sum(reduce(*,rec_space_diffs[[izs...]]) .*F,dims=izs),dims=izs)
                    F_plot = Plots.heatmap(x,y,F_sub; xlabel=x_units, ylabel=y_units, title="F* [$(units_multiply(units_inverse(x_units), units_inverse(y_units)))]", fillcolor=cgrad([:white, :darkblue, :green, :yellow, :orange, :red]))
                    F_plots[ipc] = F_plot; ipc += 1
                end
                for plot_jac_coord in plot_jac_coords
                    global ipc
                    ix = plot_jac_coord[1]; x = W_abscissas[1][ix+1]; x_units = W_abscissas_units[1][ix+1] # +1 because W_abscissas and W_abscissas_units has a first dimension of diagnostic measurement bins
                    iy = plot_jac_coord[2]; y = W_abscissas[1][iy+1]; y_units = W_abscissas_units[1][iy+1] # +1 because W_abscissas and W_abscissas_units has a first dimension of diagnostic measurement bins
                    izs = filter(xx-> !(xx==ix || xx==iy), rec_space_dims)
                    F_sub = dropdims(sum(reduce(*,rec_space_diffs[[izs...]]) .*two_pi_R_jacobian .*F,dims=izs),dims=izs)
                    F_plot = Plots.heatmap(x,y,F_sub; xlabel=x_units, ylabel=y_units, title="F* [$(units_multiply(units_inverse(x_units), units_inverse(y_units)))]", fillcolor=cgrad([:white, :darkblue, :green, :yellow, :orange, :red]))
                    F_plots[ipc] = F_plot; ipc += 1
                end
            elseif false # <--- ADD MORE 4D COORDINATE OPTIONS HERE IN FUTURE VERSIONS
            else
                plt_dummy = Plots.plot(axis=([],false)) # Dummy plot, for white space
                F_plots = [plt_dummy]
            end
        elseif false # <--- ADD MORE 5D- and 6D-DIMENSIONAL OPTIONS HERE IN FUTURE VERSIONS
        else
            plt_dummy = Plots.plot(axis=([], false)) # Dummy plot, for white space
            F_plots = [plt_dummy]
        end
        subfig_F = Plots.plot(F_plots..., layout=layout_F)
        
        lwf = layout_WF[1]
        lf = layout_F[1]
        ltot = Plots.@layout([a; b; c]) 
        fig_tot = Plots.plot(subfig_L_curve, subfig_WF, subfig_F, size=(800,400+400*lwf+400*lf), layout=ltot, left_margin=6Plots.mm, dpi=100)
        display(fig_tot)

        if plot_solutions
            filename_out_fig_default = "solInvProb_$(date_and_time)_resultsPlot_$(i)_of_$(length(hyper_points))"
            # Check if user (correctly) specified a specific output file name
            if !(typeof(filename_out)==String)
                @warn "The 'filename_out' input variable was not specified correctly (typeof(filename_out)==String evaluates to false). Default OWCF output .png file naming convention will be used"
                filename_out_fig = filename_out_fig_default
            else
                if filename_out==""
                    filename_out_fig = filename_out_fig_default
                else
                    filename_out_fig = filename_out*"_$(i)_of_$(length(hyper_points))"
                end
            end
            verbose && println("------> Saving .png file for plot $(i) of $(length(hyper_points))... ")
            png(fig_tot,folderpath_out*filename_out_fig)
        end
    end
    if gif_solutions
        filename_out_gif_default = "solInvProb_$(date_and_time)"
        # Check if user (correctly) specified a specific output file name
        if !(typeof(filename_out)==String)
            @warn "The 'filename_out' input variable was not specified correctly (typeof(filename_out)==String evaluates to false). Default OWCF output .gif file naming convention will be used"
            filename_out_gif = filename_out_gif_default
        else
            if filename_out==""
                filename_out_gif = filename_out_gif_default
            else
                filename_out_gif = filename_out
            end
        end
        verbose && println("---> Saving .gif file... ")
        gif(anim,folderpath_out*filename_out_gif*".gif",fps=2)
    end
else
    verbose && println("Plots (.png or .gif output files) of the solutions have NOT been requested.")
end

append!(timestamps,time()) # The timestamp when the results figures have been plotted
dictionary_of_sections[prnt-1] = ("Plotting solutions",diff(timestamps)[end])
###########################################################################################################
# SECTION: SAVE RECONSTRUCTION AND RELATED QUANTITIES
println("------------------------------------------------ SECTION $(prnt) ------------------------------------------------"); prnt += 1

verbose && println("Saving inverse problem solution and related quantities... ")
verbose && println("---> Inflating solutions to save-ready (N+R)-dimensional format... ")
hyper_space_size = ()
tot_space_size = deepcopy(rec_space_size) # Total space = reconstruction space + hyper-parameter space
for hyper_array in hyper_arrays
    global tot_space_size; global hyper_space_size
    hyper_space_size = tuple(hyper_space_size[:]..., length(hyper_array))
    tot_space_size = tuple(tot_space_size[:]..., length(hyper_array))
end
sols_inflated = zeros(tot_space_size) # Pre-allocate the save-ready (N+R)-dimensional format
ilm_hyper_coord = zeros(Int64, length(hyper_points[1])) # The hyper-point coordinate of the maximum L-curve curvature
ilm_hyper_point = zeros(length(hyper_points[1])) # The hyper-point of the maximum L-curve curvature
sol_ilm = zeros(rec_space_size) # The inflated solution at the maximum L-curve curvature
for i=1:length(hyper_points) # For each solution. Not the most efficient. But ok, since the total number of hyper-points/solutions will always be relatively small (<1000)
    verbose && println("------> Inflating solution $(i) of $(length(hyper_points))... ")
    sol_i = sols[:,i]
    # If slowing-down basis function coefficients, compute fast-ion distribution
    if collision_physics_reg
        sol_i = F_SD_safe*sol_i
    end
    hp = hyper_points[i]
    hp = vcat(hp...) # Tuple to Vector
    hp_indices = zeros(Int64, length(hp))
    for (ip,reg_strength) in enumerate(hp)
        hp_indices[ip] = findfirst(x-> x==reg_strength,hyper_arrays[ip])
    end
    if i==ilm # If it's the hyper-point with maximum L-curve curvature (happens only once)
        ilm_hyper_coord[:] = hp_indices # Save its coordinate (in hyper space)
        ilm_hyper_point[:] = hp # Save it
        for (ic,c) in enumerate(rec_space_coords)
            totc = vcat(Tuple(c)...,hp_indices...)
            sols_inflated[totc...] = sol_i[ic] # Add to the inflated ((N+R)-dimensional) array of solutions (.= must be used, otherwise Julia breakes. Even though no broadcasting will ever be used here)
            sol_ilm[c] = sol_i[ic] # Construct the inflated (N-dimensional) solution
        end
    else # Otherwise.. (note: This should be the most efficient approach, to have the i==ilm check be performed only length(hyper_points) times, instead of length(hyper_points)*length(rec_space_coords) times)
        for (ic,c) in enumerate(rec_space_coords)
            totc = vcat(Tuple(c)...,hp_indices...)
            sols_inflated[totc...] = sol_i[ic] # Add to the inflated ((N+R)-dimensional) array of solutions (.= must be used, otherwise Julia breakes. Even though no broadcasting will ever be used here)
        end
    end
end
if collision_physics_reg
    verbose && println("---> Inflating collision physics basis functions to N-dimensional format... ")
    F_SD_safe_inflated = zeros(tuple(rec_space_size[:]...,size(F_SD_safe,2)))
    for i=1:size(F_SD_safe,2)
        verbose && println("------> Inflating basis function $(i) of $(size(F_SD_safe,2))... ")
        for (ic,c) in enumerate(rec_space_coords)
            totc = vcat(Tuple(c)...,i)
            F_SD_safe_inflated[totc...] = F_SD_safe[ic,i]
        end
    end
end
verbose && println("---> Creating unique output data file name... ")
# Create default output file name
recSpaceGridSize = reduce(*,map(x-> "x$(x)",rec_space_size))[2:end]
hyperParamGridSize = reduce(*,map(x-> "x$(x)",hyper_space_size))[2:end]
filename_out_default = "solInvProb_$(date_and_time)_$(tokamak)_S_$(length(S))_$(length(S_hat_long))_W_$(length(W))_$(recSpaceGridSize)_$(hyperParamGridSize)"

# Check if user (correctly) specified a specific output file name
if !(typeof(filename_out)==String)
    @warn "The 'filename_out' input variable was not specified correctly (typeof(filename_out)==String evaluates to false). Default OWCF output file naming convention will be used"
    filename_out = filename_out_default
else
    if filename_out==""
        filename_out = filename_out_default
    end
end

# Adding _(1), _(2) etc if necessary
filepath_output_orig = folderpath_out*filename_out
filepath_output = deepcopy(filepath_output_orig)
C = 1
while isfile(filepath_output*".jld2") || isfile(filepath_output*".hdf5") # To take care of not overwriting files. Add _(1), _(2) etc
    global filepath_output; global C
    filepath_output = filepath_output_orig*"_($(Int64(C)))"
    C += 1
end

# Opening .jld2 or .h5 file
verbose && println("---> Opening save file... ")
if saveInHDF5format
    verbose && println("------> Output data will be saved to file $(filepath_output*".hdf5")... ")
    myfile = h5open(filepath_output*".hdf5","w")
else
    verbose && println("------> Output data will be saved to file $(filepath_output*".jld2")... ")
    myfile = jldopen(filepath_output*".jld2",true,true,false,IOStream)
end

# Saving data
verbose && println("---> Saving output data... ")
write(myfile,"sols",sols_inflated)
for (i,ab) in enumerate(W_abscissas[1][2:end])
    if saveInHDF5format
        ab = Float64.(ab) # Convert from Real to Float64, since .hdf5 files cannot handle Real types
    end
    write(myfile,"sols_abscissa_$(i)",ab)
    write(myfile,"sols_abscissa_units_$(i)",W_abscissas_units[1][i+1])
end
write(myfile, "sols_hyper_points", hyper_points)
for (i,hp_ab) in enumerate(hyper_arrays)
    write(myfile,"sols_hyper_abscissa_$(i)",hp_ab)
end
write(myfile,"FI_species",FI_species)
for (j,s) in enumerate(S)
    write(myfile,"S_$(j)",s)
    write(myfile,"S_units_$(j)",S_units[j])
    write(myfile,"S_err_$(j)",S_errs[j])
    write(myfile,"S_err_units_$(j)",S_errs_units[j])
    write(myfile,"S_abscissa_$(j)",S_abscissas[j])
    write(myfile,"S_abscissa_units_$(j)",S_abscissas_units[j])
end
write(myfile,"l_curve_x",sdf)
write(myfile,"l_curve_y",srn)
write(myfile,"l_curve_opt_index",ilm)
write(myfile,"l_curve_opt_hyper_point",ilm_hyper_point)
write(myfile,"l_curve_opt_hyper_coord",ilm_hyper_coord)
write(myfile,"l_curve_opt_sol",sol_ilm)
write(myfile,"filepath_start",filepath_start)
if rescale_W
    write(myfile,"rescale_W_factors",rescale_W_factors)
end
if collision_physics_reg || icrf_physics_reg
    write(myfile,"regularization_equil_filepath",regularization_equil_filepath)
end
if collision_physics_reg
    write(myfile,"coll_phys_basis",F_SD_safe_inflated)
    write(myfile,"coll_phys_thermal_species",regularization_thermal_ion_species)
    write(myfile,"coll_phys_Ti",thermal_ion_temperatures)
    write(myfile,"coll_phys_ni",thermal_ion_densities)
    write(myfile,"coll_phys_Te",thermal_electron_temperatures)
    write(myfile,"coll_phys_ne",thermal_electron_densities)
    write(myfile,"coll_phys_rho",rho_of_interests)
end
close(myfile)

append!(timestamps,time()) # The timestamp when the output data has been saved
dictionary_of_sections[prnt-1] = ("Saving output data",diff(timestamps)[end])
###########################################################################################################
# FINAL SECTION: PRINT COMPLETION STATEMENT AND STATS
println("------------------------------------------------ SECTION $(prnt) ------------------------------------------------")

timesum = sum(diff(timestamps))

println("------------------------------ solveInverseProblem.jl completed successfully! ------------------------------")
println("---> Total script run time: $(round(timesum,sigdigits=5)) seconds, where")
for key in sort([k for k in keys(dictionary_of_sections)]) # To get the keys sorted correctly
    println("------> SECTION $(key). $(dictionary_of_sections[key][1]) took $(round(dictionary_of_sections[key][2],sigdigits=3)) seconds ($(round(100*dictionary_of_sections[key][2]/timesum,sigdigits=3)) %)")
end
println("------------------------------------------------------------------------------------------------------------")