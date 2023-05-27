###################################### extractNullOrbits.jl #########################################

#### Description:
# This script examines a diagnostic's signal and notices where the measurements are zero (i.e. 
# below the signalNullThreshold). It will then proceed to examine the orbit weight function for that 
# specific diagnostic's energy bin, and note for which valid orbits the orbit sensitivity is 
# non-zero (above orbNullThreshold). Those orbits will be labelled as null orbits. By logic, 
# these orbits must be unpopulated since otherwise they would result in a non-zero measurement 
# signal (because the orbit sensitivity is non-zero). Those null orbits are then saved as '1.0' 
# in a vector equal in length to size(W,2) (and F_os_raw from ps2os.jl, if the same orbit grid was used).

#### Inputs (units given when defined in the script):
# Please see start_extrNullOrbs_template.jl

#### Outputs
# -

#### Saved files
# nullOrbs_[tokamak]_[TRANSP ID]_at[timepoint]s_[diagnostic]_[FI species]_[nE]x[npm]x[nRm].jld2
# This saved file will have the fields:
#   nullOrbs - A 2D array where '1.0' means a null orbit, and '0.0' means not. Same column ordering of orbits 
#              as that of F_os_raw from ps2os.jl - Array{Float64,2}
#   Ed_array - The measurement bins of the diagnostic spectrum - Array{Float64,1}
#   E_array - The fast-ion energy grid array used for orbit space - Array{Float64,1}
#   pm_array - The fast-ion pm grid array used for orbit space - Array{Float64,1}
#   Rm_array - The fast-ion Rm grid array used for orbit space - Array{Float64,1}
#   filepath_W - The filepath to the orbit weights file used - String
#   filepath_S - The filepath to the signal file used - String 
# if include2Dto4D, the saved file will also have the key
#   nullOrbs_indices - A vector containing cartesian 4D indices, marking the orbit-sapce location of the null orbits for each orbit weight function. size: length(Ed_array)xlength(E_array)xlength(pm_array)xlength(Rm_array) - Vector{CartesianIndex{4}}

# Script written by Henrik Järleblad. Last maintained 2022-09-20.
#####################################################################################################

## ------
# Loading packages
verbose && println("Loading Julia packages... ")
@everywhere begin
    using JLD2
    using FileIO
    using SparseArrays
    include(folderpath_OWCF*"misc/species_func.jl")
end

if include2Dto4D
    verbose && println("Loading extra Julia packages for include2Dto4D... ")
    @everywhere begin
        using EFIT
        using Equilibrium
        using GuidingCenterOrbits
        using OrbitTomography
    end
end

## ------
if include2Dto4D
    # Loading tokamak equilibrium
    verbose && println("Loading equilibrium... ")
    if ((split(filepath_equil,"."))[end] == "eqdsk") || ((split(filepath_equil,"."))[end] == "geqdsk")
        M, wall = read_geqdsk(filepath_equil,clockwise_phi=false) # Assume counter-clockwise phi-direction
        jdotb = M.sigma # The sign of the dot product between the plasma current and the magnetic field

        # Extract timepoint information from .eqdsk/.geqdsk file
        eqdsk_array = split(filepath_equil,".")
        XX = (split(eqdsk_array[end-2],"-"))[end] # Assume format ...-XX.YYYY.eqdsk where XX are the seconds and YYYY are the decimals
        YYYY = eqdsk_array[end-1] # Assume format ...-XX.YYYY.eqdsk where XX are the seconds and YYYY are the decimals
        timepoint = (timepoint == nothing ? XX*","*YYYY : timepoint) # Format XX,YYYY to avoid "." when including in filename of saved output
    else # Otherwise, assume magnetic equilibrium is a saved .jld2 file
        myfile = jldopen(filepath_equil,false,false,false,IOStream)
        M = myfile["S"]
        wall = myfile["wall"]
        close(myfile)
        jdotb = (M.sigma_B0)*(M.sigma_Ip)
        timepoint = (timepoint == nothing ? "00,0000" : timepoint)
    end
end

filepath_W_array = split(filepath_W,"_")
time_string = filepath_W_array[end-3] # Assume this to start with
if (time_string[1]=='a') && (time_string[2]=='t') && (time_string[end]=='s')
    # Correct format identified. Superseeds user input and/or timepoint found via equilibrium file
    timepoint = time_string[3:end-1] # First to elements are "at" and last element is "s"
else
    timepoint = (timepoint == nothing ? "00,0000" : timepoint)
end

## ------
# Loading weights file
verbose && println("Loading weight matrix data from filepath_W... ")
myfile = jldopen(filepath_W,false,false,false,IOStream)
if haskey(myfile,"W")
    W2D = myfile["W"]
else
    error(".jld2 key for weight matrix not recognized. Please change to 'W' and use .jld2 file as input again.")
end
E_array = myfile["E_array"]
pm_array = myfile["pm_array"]
Rm_array = myfile["Rm_array"]
Ed_array = myfile["Ed_array"]
close(myfile)

## ------
# Loading signal file
verbose && println("Loading signal data from filepath_S... ")
myfile = jldopen(filepath_S,false,false,false,IOStream)
if haskey(myfile,"spec")
    S = myfile["spec"]
elseif haskey(myfile,"S")
    S = myfile["S"]
elseif haskey(myfile,"s")
    S = myfile["s"]
else
    error(".jld2 key for diagnostic signal not recognized. Please change to 'S', 's' or 'spec' and use .jld2 file as input again.")
end
Ed_array_S = myfile["Ed_array"]
close(myfile)

## ------
verbose && println("Double-checking number of diagnostic energy bins... ")
if !(length(Ed_array_S)==length(Ed_array))
    error("Number of diagnostic energy bins in W does not match S. Please correct and re-try.")
end

## ------
# Extract info from filepath_W file path
verbose && println("Extracting various info from filepath_W file path... ")
fW_array = split(filepath_W,"_")
tokamak = fW_array[end-4]
TRANSP_id = fW_array[end-3]
diagnostic = fW_array[end-2]
reaction = fW_array[end-1]

## ------
# Printing script info and inputs
println("")
println("-------------------------------------------------extractNullOrbits.jl-------------------------------------------------")
println("Null orbits will be extracted for a $(length(E_array))x$(length(pm_array))x$(length(Rm_array)) orbit grid from a $(size(W2D,1))x$(size(W2D,2)) orbit weight matrix")
println("for the "*diagnostic*" diagnostic.")
println("")
println("The signal threshold for null signals will be below $(signalNullThreshold).")
println("The orbit threshold for null orbits will be above $(orbNullThreshold).")
println("")
println("Weights file specified: ")
println("--- : "*filepath_W)
println("")
println("Signal file specified: ")
println("--- : "*filepath_S)
println("")
println("Fast-ion species specified:")
println("--- :"*FI_species)
println("")
if include2Dto4D
    println("Orbit-space indices of null orbits will be computed.")
    println("")
    println("Magnetic equilibrium file specified:")
    println("--- : "*filepath_equil)
end
println("")
println("Extra keyword arguments specified: ")
println(extra_kw_args)
println("")
println("Written by Henrik Järleblad. Last maintained 2022-09-20.")
println("----------------------------------------------------------------------------------------------------------------------")
println("")

if !(@isdefined og) && include2Dto4D # Could be running via some other script that's already computed orbit-grid og
    verbose && println("Computing orbit grid... ")
    og_orbs, og = OrbitTomography.orbit_grid(M, E_array, pm_array, Rm_array; amu=getSpeciesAmu(FI_species), q=getSpeciesEcu(FI_species), wall=wall, extra_kw_args...)
    if !(length(og_orbs)==size(W2D,2))
        println("Computation of orbit grid did not produce same number of valid orbits as that of loaded weight matrix.")
        println("Please ensure that you are using the same orbit algorithm inputs as when you computed the weight matrix and try again.")
        error("extractNullOrbits.jl must terminate.")
    end
    og_orbs = nothing # Minimize memory
end

## ------
# Compute null orbits
# Assume (E,pm,Rm) dimensions order
verbose && println("Extracting null orbits... ")
null_orbits_2D = zeros(size(W2D))
for i=1:length(S)
    if S[i] <= signalNullThreshold*maximum(S)
        max_W = maximum(W2D[i,:])
        global null_orbits_2D[i,:] = map(x -> x > orbNullThreshold*max_W ? 1.0 : 0.0, W2D[i,:])
    end
end
null_orbits_2D = map(x -> x > 0.0 ? 1.0 : 0.0, null_orbits_2D) # Unnecessary? Let this line be.

if include2Dto4D
    verbose && println("Extracting cartesian orbit-space indices for null orbits... ")
    null_orbits_4D = zeros(length(S),length(E_array),length(pm_array),length(Rm_array))
    for i=1:size(null_orbits_4D,1)
        null_orbits_4D[i,:,:,:] = [oi == 0 ? zero(null_orbits_2D[1,1]) : null_orbits_2D[i,oi] for oi in og.orbit_index]
    end
    nz_inds = findall(!iszero, null_orbits_4D) # Store only the non-zero indices. For efficient storage
end

## ------
# Saving
verbose && println("Saving... ")
global filepath_output_orig = folderpath_o*"nullOrbits_"*tokamak*"_"*TRANSP_id*"_at"*timepoint*"s_"*diagnostic*"_"*reaction*"_$(length(E_array))x$(length(pm_array))x$(length(Rm_array))"
global filepath_output = deepcopy(filepath_output_orig)
global count = 1
while isfile(filepath_output*".jld2") # To take care of not overwriting files. Add _(1), _(2) etc
    global filepath_output = filepath_output_orig*"_($(Int64(count)))"
    global count += 1 # global scope, to surpress warnings
end
global filepath_output = filepath_output*".jld2"
myfile = jldopen(filepath_output,true,true,false,IOStream)
write(myfile,"nullOrbs",null_orbits_2D)
if include2Dto4D
    write(myfile,"nullOrbs_indices",nz_inds)
end
write(myfile,"Ed_array",Ed_array)
write(myfile,"E_array",E_array)
write(myfile,"pm_array",pm_array)
write(myfile,"Rm_array",Rm_array)
write(myfile,"filepath_W",filepath_W)
write(myfile,"filepath_S",filepath_S)
close(myfile)

println("~~~extractNullOrbits.jl done!~~~")