################################### orbWeights_2Dto4D.jl ##########################################

#### Description:
# This script will convert a 2D matrix (representing an orbit weight function) to it's inflated
# 4D form. The 2D matrix .jld2-file could be from calcOrbWeights.jl for example, and have rows
# corresponding to channels and columns corresponding to (valid) orbits. The orbits are ordered according to
# the scheme by L. Stagner's OrbitTomography.jl Julia package (specifically, as in orbits.jl/OrbitGrid.orbit_index).

#### Inputs (units given when defined in the script):
# Please see the start_2Dto4D_template.jl file.

#### Outputs
# -

#### Saved files
# orbWeights4D_[tokamak]_[TRANSP_id]_at[timepoint]s_[diagnostic_name]_[reaction/FI_species]_[nEd]x[nE]x[npm]x[nRm].jld2
#   This saved file will have the fields:
#   Wtot - The computed NES orbit weights. Dimensions are (channels, nE, npm, nRm) - Array{Float64,2}
#   E_array - The fast-ion energy grid array used for orbit space. Length nE - Array{Float64,1}
#   pm_array - The fast-ion pm grid array used for orbit space. Length npm - Array{Float64,1}
#   Rm_array - The fast-ion Rm grid array used for orbit space. Length nRm - Array{Float64,1}
#   Ed_array - The diagnostic energy bin centers. Length nEd - Array{Float64,1}
# If the compact form of the fusion reaction (a-b) is defined, the saved file will also include the key
#   reaction - The reactants of the fusion reaction for the orbit weight functions. 'a-b', where a is thermal, b is fast - String
# If the full form of the fusion reaction a(b,c)d is defined, the saved file will also include the key
#   reaction_full - The fusion reaction a(b,c)d. a is thermal, b is fast, c is emitted particle and d is the product nucleus - String
# If analical orbit weight functions were inflated, the saved file will also include the key
#   analyticalOWs - If true, then the orbit weight functions were computed using projected velocities - Bool

### Other
# You should delete the orbWeights4D file once you're done with it. These files are massive and take
# up an unnessecary amount of space in general.

# Script written by Henrik Järleblad. Last maintained 2022-10-09.
###################################################################################################

## ------
# Loading packages
println("Loading packages... ")
@everywhere begin
    cd(folderpath_OWCF)
    using EFIT # For calculating magnetic equilibrium quantities
    using Equilibrium # For loading flux function data, tokamak geometry data etc.
    using GuidingCenterOrbits # For calculating guiding-center orbits
    using OrbitTomography # For computing orbit grid structs
    using JLD2 # To write/open .jld2-files (Julia files, basically)
    using FileIO # To write/open files in general
    using SparseArrays # To enable utilization of sparse matrices/vectors
    include(folderpath_OWCF*"misc/species_func.jl") # To convert species label to particle mass
    include(folderpath_OWCF*"misc/availReacts.jl") # To check reaction availability and extract fast-ion and thermal species
    include(folderpath_OWCF*"misc/rewriteReacts.jl") # To rewrite a fusion reaction from the A(b,c)D format to the A-b=c-D format
end

## ------
# Loading tokamak equilibrium
verbose && println("Loading tokamak equilibrium... ")
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

filepath_W_array = split(filepath_W,"_")
time_string = filepath_W_array[end-3] # Assume this to start with
if (time_string[1]=='a') && (time_string[2]=='t') && (time_string[end]=='s')
    # Correct format identified. Superseeds user input and/or timepoint found via equilibrium file
    timepoint = time_string[3:end-1] # First to elements are "at" and last element is "s"
else
    timepoint = (timepoint == nothing ? "00,0000" : timepoint)
end

## ------
# Loading orbit  weights and data
verbose && println("Loading orbit weights... ")
myfile = jldopen(filepath_W,false,false,false,IOStream)
if haskey(myfile,"Wtot")
    W2D = myfile["Wtot"]
else
    W2D = myfile["W"]
end
E_array = myfile["E_array"]
pm_array = myfile["pm_array"]
Rm_array = myfile["Rm_array"]
if haskey(myfile,"En_array")
    Ed_array = myfile["En_array"]
elseif haskey(myfile,"Ed_array")
    Ed_array = myfile["Ed_array"]
else
    error("orbWeights_2Dto4D did not recognize known diagnostic data array type in provide 'filepath_W'. Please re-try another file.")
end
if haskey(myfile,"reaction")
    reaction = myfile["reaction"]
end
if haskey(myfile,"reaction_full")
    reaction_full = myfile["reaction_full"]
end
if haskey(myfile,"analyticalOWs")
    analyticalOWs = true
else
    analyticalOWs = false
end
close(myfile)

## ---------------------------------------------------------------------------------------------
# Determine fast-ion species from reaction
if (@isdefined reaction_full)
    thermal_species, FI_species = checkReaction(reaction_full; verbose=verbose, projVelocity=analyticalOWs)
    @everywhere FI_species = $FI_species # Transfer variable to all external processes
elseif (@isdefined reaction)
    FI_species = (split(reaction,"-"))[1] # Assume first species specified in reaction to be the fast-ion species. For example, in 'p-t' the 'p' will be assumed the fast ion species.
    @everywhere FI_species = $FI_species # Transfer variable to all external processes
else
    FI_species = FI_species # If no reaction has been specified, assume fast-ion species to be deuterium
end

## ---------------------------------------------------------------------------------------------
# Create modified Ed_array from Ed_array (loaded from filepath_W) and Ed_energies
verbose && println("Determining which diagnostic energy bins to extract orbit weight functions from... ")
if !(@isdefined Ed_energies)
    Ed_energies = nothing
end
if (Ed_energies===nothing)
    Edi_array = collect(1:length(Ed_array)) # Use all diagnostic energy bins
elseif typeof(Ed_energies)==Int
    if !(Ed_energies<=floor(length(Ed_array)/2))
        error("Ed_energies > half of number of diagnostic energy bins. Please correct and re-try!")
    end
    Edi = 1
    Edi_array = []
    while Edi<=length(Ed_array) # Not efficient. But it's ok, because array lengths are expected to be small
        global Edi # Declare that you are talking about the global Edi variable
        push!(Edi_array,Edi) # Use only diagnostic energy bins that are evenly spaced along the indices of Ed_array
        global Edi = Edi + floor(length(Ed_array)/Ed_energies) # 'global' declaration not required in Julia 1.5 and above
    end
else
    # Checking typeof Ed_energies
    if !(typeof(Ed_energies)<:AbstractVector)
        error("Ed_energies not specified correctly. Please correct and re-try.")
    end
    # Again, this code is not efficient, but it's ok thanks to assumed small array lengths
    Edi_array = []
    for Ed_energy in Ed_energies
        Ed_energy_array = Ed_energy .* ones(size(Ed_array))
        Edi = findmin(abs.(Ed_energy_array-Ed_array))[2] # Find closest value to Ed_energy in Ed_array

        if !(Edi in Edi_array) # If it hasn't already been added to Edi_array
            push!(Edi_array,Edi) # Use only diagnostic energy bins that are sufficiently close to Ed_energies
        end
    end
end
Edi_array = Int.(Edi_array) # Just make sure that all elements in Edi_array are Ints

## ---------------------------------------------------------------------------------------------
# Printing script info and inputs
println("")
println("---------------------------------orbWeights_2Dto4D.jl---------------------------------")
print("Magnetic equilibrium file specified: ")
print(filepath_equil)
print("    ")
println("Diagnostic name: "*diagnostic_name)
println("")
if (@isdefined reaction)
    println("Fusion reaction: "*reaction)
end
if (@isdefined reaction_full)
    println("Fusion reaction: "*reaction_full)
end
println("Fast-ion species specified: "*FI_species)
println("")
println("Orbit weight functions will be inflated for a $(length(E_array))x$(length(pm_array))x$(length(Rm_array)) orbit-space grid with")
println("Energy: [$(minimum(E_array)),$(maximum(E_array))] keV")
println("Pitch maximum: [$(minimum(pm_array)),$(maximum(pm_array))]")
println("Radius maximum: [$(minimum(Rm_array)),$(maximum(Rm_array))] m")
println("")
if distributed
    println("This will be done using $(nprocs()) CPU cores.")
    println("")
end
if !(Ed_energies==nothing)
    println("The diagnostic measurement bins of interest are: $(Ed_energies)")
    println("")
end
println("The resulting 4D matrix will have the dimensions $(length(Edi_array))x$(length(E_array))x$(length(pm_array))x$(length(Rm_array))")
println("")
print("Tokamak: "*tokamak)
print("     ")
println("TRANSP ID: "*TRANSP_id)
println("")
println("Extra keyword arguments specified: ")
println(extra_kw_args)
println("")
println("If you would like to change any settings, please edit the start_2Dto4D_template.jl file.")
println("Written by Henrik Järleblad. Last maintained 2022-10-09.")
println("--------------------------------------------------------------------------------------")
println("")

## ---------------------------------------------------------------------------------------------
# Computing orbit orbit grid
verbose && println("Computing orbit grid... ")
og_orbs, og = orbit_grid(M, E_array, pm_array, Rm_array; q=getSpeciesEcu(FI_species), amu=getSpeciesAmu(FI_species), wall=wall,extra_kw_args...)
og_orbs=nothing

## ---------------------------------------------------------------------------------------------
# Converting orbit weight function
verbose && println("Converting orbit weight functions... ")
dE = abs(E_array[2]-E_array[1]) # Assume equidistant orbit-space grid
dpm = abs(pm_array[2]-pm_array[1]) # Assume equidistant orbit-space grid
dRm = abs(Rm_array[2]-Rm_array[1]) # Assume equidistant orbit-space grid
sort!(Edi_array) # Sort the extracted diagnostic energy bin indices, to ensure W4D is consistent
W4D = zeros(length(Edi_array),length(E_array),length(pm_array),length(Rm_array))
for (i,Edi) in enumerate(Edi_array)
    verbose && println("Converting orbit weight function for diagnostic energy bin with energy $(round(Ed_array[Edi],digits=2)) keV ($(i)/$(length(Edi_array)))... ")
    W3D = map_orbits(og, Vector(W2D[Edi,:]),true) * (dE*dpm*dRm) # Assume equidistant orbit-space grid. For now.
    W4D[i,:,:,:] .= W3D
end

## ---------------------------------------------------------------------------------------------
# Saving result
input_folder = split(filepath_W,"/")
output_folder = ""
for i=1:length(input_folder)-1
    global output_folder = output_folder*input_folder[i] # Put the output in the same folder as the weights input
    global output_folder = output_folder*"/" # Put the output in the same folder as the weights input
end

verbose && println("Saving inflated 4D orbit weight functions... ")
if (@isdefined reaction_full)
    global filepath_output_orig = output_folder*"orbWeights4D_"*tokamak*"_"*TRANSP_id*"_at"*timepoint*"s_"*diagnostic_name*"_"*pretty2scpok(reaction_full; projVel = analyticalOWs)*"_$(length(Ed_array[Edi_array]))x$(length(E_array))x$(length(pm_array))x$(length(Rm_array))"
elseif (@isdefined reaction)
    global filepath_output_orig = output_folder*"orbWeights4D_"*tokamak*"_"*TRANSP_id*"_at"*timepoint*"s_"*diagnostic_name*"_"*reaction*"_$(length(Ed_array[Edi_array]))x$(length(E_array))x$(length(pm_array))x$(length(Rm_array))"
else
    global filepath_output_orig = output_folder*"orbWeights4D_"*tokamak*"_"*TRANSP_id*"_at"*timepoint*"s_"*diagnostic_name*"_"*FI_species*"_$(length(Ed_array[Edi_array]))x$(length(E_array))x$(length(pm_array))x$(length(Rm_array))"
end
global filepath_output = deepcopy(filepath_output_orig)
global count = 1
while isfile(filepath_output*".jld2") # To take care of not overwriting files. Add _(1), _(2) etc
    global filepath_output = filepath_output_orig*"_($(Int64(count)))"
    global count += 1 # global scope, to surpress warnings
end
global filepath_output = filepath_output*".jld2"
myfile = jldopen(filepath_output,true,true,false,IOStream)
write(myfile, "Wtot", W4D)
write(myfile, "E_array", E_array)
write(myfile, "pm_array", pm_array)
write(myfile, "Rm_array", Rm_array)
write(myfile, "Ed_array", Ed_array[Edi_array]) # Save only diagnostic energy bins that have been used
if (@isdefined reaction)
    write(myfile, "reaction", reaction)
end
if (@isdefined reaction_full)
    write(myfile, "reaction_full", reaction_full)
end
if analyticalOWs
    write(myfile, "analyticalOWs", analyticalOWs)
end
close(myfile)
println("Done!")
