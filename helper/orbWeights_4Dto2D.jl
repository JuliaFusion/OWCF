################################### orbWeights_4Dto2D.jl ##########################################

#### Description:
# This script will convert a 4D matrix representing an orbit weight function, to it's compressed
# 2D form. The 4D matrix .jld2-file could be from e.g. orbWeights_2Dto4D.jl, and have the dimensions 
# (nEd)x(nE)x(npm)x(nRm) where nEd is the number of diagnostic measurement bins, nE is the number
# of orbit-space energy grid points, npm is the number of orbit-space pitch maximum grid points and
# nRm is the number of orbit-space radius maximum grid points. The 2D matrix will have the dimensions 
# (nEd)x(valid orbits) where the valid orbits are ordered according to the scheme by L. Stagner's 
# OrbitTomography.jl Julia package (specifically, as in orbits.jl/OrbitGrid.orbit_index).

#### Inputs (units given when defined in the script):
# Please see the start_4Dto2D_template.jl file.

#### Outputs
# -

#### Saved files
# orbWeights_[tokamak]_[TRANSP_id]_at[timepoint]s_[diagnostic]_[reaction/FI species]_[nE]x[npm]x[nRm].jld2
# This saved file will have the fields:
#   Wtot - The computed orbit weights. Dimensions are (channels,valid orbits) - Array{Float64,2}
#   E_array - The fast-ion energy grid array used for orbit space - Array{Float64,1}
#   pm_array - The fast-ion pm grid array used for orbit space - Array{Float64,1}
#   Rm_array - The fast-ion Rm grid array used for orbit space - Array{Float64,1}
#   Ed_array - The diagnostic energy bin centers - Array{Float64,1}
# If the compact form of the fusion reaction (a-b) is defined, the saved file will also include the key
#   reaction - The reactants of the fusion reaction for the orbit weight functions. 'a-b', where a is thermal, b is fast - String
# If the full form of the fusion reaction a(b,c)d is defined, the saved file will also include the key
#   reaction_full - The fusion reaction a(b,c)d. a is thermal, b is fast, c is emitted particle and d is the product nucleus - String

### Other
# You should delete the orbWeights4D file once you're done with it. These files are massive and take
# up an unnessecary amount of space, in general. Or save it externally on a USB memory etc.

# Script written by Henrik Järleblad. Last maintained 2022-09-09.
###################################################################################################

verbose && println("Loading Julia packages... ")
@everywhere begin
    cd(folderpath_OWCF)
    using JLD2 # To write/open .jld2-files (Julia data files, basically)
    using EFIT # For calculating magnetic equilibrium quantities
    using Equilibrium # For loading flux function data, tokamak geometry data etc.
    using GuidingCenterOrbits # For calculating guiding-center orbits
    using OrbitTomography # This is what all this is about!
    include(folderpath_OWCF*"extra/dependencies.jl") # To load dependencies functions
    include(folderpath_OWCF*"misc/species_func.jl") # To convert species label to particle mass
    include(folderpath_OWCF*"misc/availReacts.jl") # To check reaction availability and extract fast-ion and thermal species
    include(folderpath_OWCF*"misc/rewriteReacts.jl") # To rewrite a fusion reaction from the A(b,c)D format to the A-b=c-D format
end

verbose && println("Loading orbit weight function and orbit-space grid arrays... ")
myfile = jldopen(filepath_W,false,false,false,IOStream)
Wtot = myfile["Wtot"]
E_array = myfile["E_array"]
pm_array = myfile["pm_array"]
Rm_array = myfile["Rm_array"]
if haskey(myfile,"Ed_array")
    Ed_array = myfile["Ed_array"]
elseif haskey(myfile,"En_array")
    Ed_array = myfile["En_array"]
else
    error("orbWeights_4Dto2D did not recognize known diagnostic data array type in provide 'filepath_W'. Please re-try another file.")
end
if haskey(myfile,"reaction")
    reaction = myfile["reaction"]
end
if haskey(myfile,"reaction_full")
    reaction_full = myfile["reaction_full"]
end
close(myfile)

## ---------------------------------------------------------------------------------------------
# Determine fast-ion species from reaction
if (@isdefined reaction_full)
    thermal_species, FI_species = checkReaction(reaction_full)
    @everywhere FI_species = $FI_species # Transfer variable to all external processes
elseif (@isdefined reaction)
    FI_species = (split(reaction,"-"))[1] # Assume first species specified in reaction to be the fast-ion species. For example, in 'p-t' the 'p' will be assumed the thermal species.
    @everywhere FI_species = $FI_species # Transfer variable to all external processes
else
    FI_species = FI_species # If no reaction has been specified, assume fast-ion species to be what is specified in input file
end

## ------
verbose && println("Loading magnetic equilibrium... ")
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
# Printing script info and inputs
println("")
println("--------------------orbWeights_4Dto2D.jl--------------------")
print("Magnetic equilibrium file specified: ")
print(filepath_equil)
print("    ")
println("Diagnostic name: "*diagnostic)
println("")
if (@isdefined reaction)
    println("Fusion reaction: "*reaction)
end
if (@isdefined reaction_full)
    println("Fusion reaction: "*reaction_full)
end
println("Fast-ion species specified: "*FI_species)
println("")
println("Orbit weight functions will be compressed for a $(length(Ed_array))x$(length(E_array))x$(length(pm_array))x$(length(Rm_array)) orbit weight functions with")
println("Energy: [$(minimum(E_array)),$(maximum(E_array))] keV")
println("Pitch maximum: [$(minimum(pm_array)),$(maximum(pm_array))]")
println("Radius maximum: [$(minimum(Rm_array)),$(maximum(Rm_array))] m")
println("")
if distributed
    println("This will be done using $(nprocs()) CPU cores.")
    println("")
end
println("The resulting 2D matrix will have the dimensions $(length(Ed_array))x(valid orbits)")
println("")
print("Tokamak: "*tokamak)
print("     ")
println("TRANSP ID: "*TRANSP_id)
println("")
println("Extra keyword arguments specified: ")
println(extra_kw_args)
println("")
println("If you would like to change any settings, please edit the start_4Dto2D_template.jl file.")
println("")
println("Written by Henrik Järleblad. Last maintained 2022-09-09.")
println("------------------------------------------------------------")
println("")

## ------
# Computing orbit orbit grid
verbose && println("Computing orbit grid... ")
og_orbs, og = orbit_grid(M, E_array, pm_array, Rm_array; q=getSpeciesEcu(FI_species), amu=getSpeciesAmu(FI_species), wall=wall, extra_kw_args...)

## ------
verbose && println("Number of valid orbits: $(length(og_orbs))")

## ------
Wtot2D = zeros(length(Ed_array),length(og_orbs))
dE = abs(E_array[2]-E_array[1]) # Assume equidistant orbit-space grid
dpm = abs(pm_array[2]-pm_array[1]) # Assume equidistant orbit-space grid
dRm = abs(Rm_array[2]-Rm_array[1]) # Assume equidistant orbit-space grid
for Edi=1:length(Ed_array)
    verbose && println("Unmapping $(Edi) of $(length(Ed_array))...")
    Wtot2D[Edi,:] .= unmap_orbits(og, Wtot[Edi,:,:,:]) ./ (dE*dpm*dRm) # Assume equidistant orbit-space grid
end

## ------
# Saving result
input_folder = split(filepath_W,"/")
output_folder = ""
for i=1:length(input_folder)-1
    global output_folder = output_folder*input_folder[i] # Put the output in the same folder as the weights input
    global output_folder = output_folder*"/" # Put the output in the same folder as the weights input
end

verbose && println("Saving compressed 2D orbit weight functions... ")
if (@isdefined reaction_full)
    global filepath_output_orig = output_folder*"orbWeights2Dfr4D_"*tokamak*"_"*TRANSP_id*"_at"*timepoint*"s_"*diagnostic*"_"*pretty2scpok(reaction_full)*"_$(length(E_array))x$(length(pm_array))x$(length(Rm_array))"
elseif (@isdefined reaction)
    global filepath_output_orig = output_folder*"orbWeights2Dfr4D_"*tokamak*"_"*TRANSP_id*"_at"*timepoint*"s_"*diagnostic*"_"*reaction*"_$(length(E_array))x$(length(pm_array))x$(length(Rm_array))"
else
    global filepath_output_orig = output_folder*"orbWeights2Dfr4D_"*tokamak*"_"*TRANSP_id*"_at"*timepoint*"s_"*diagnostic*"_"*FI_species*"_$(length(E_array))x$(length(pm_array))x$(length(Rm_array))"
end
global filepath_output = deepcopy(filepath_output_orig)
global count = 1
while isfile(filepath_output*".jld2") # To take care of not overwriting files. Add _(1), _(2) etc
    global filepath_output = filepath_output_orig*"_($(Int64(count)))"
    global count += 1 # global scope, to surpress warnings
end
global filepath_output = filepath_output*".jld2"
myfile = jldopen(filepath_output,true,true,false,IOStream)
write(myfile,"Wtot", Wtot2D)
write(myfile,"E_array",E_array)
write(myfile,"pm_array",pm_array)
write(myfile,"Rm_array",Rm_array)
write(myfile,"Ed_array",Ed_array)
if (@isdefined reaction)
    write(myfile, "reaction", reaction)
end
if (@isdefined reaction_full)
    write(myfile, "reaction_full", reaction_full)
end
close(myfile)

println("~~~~~~Done!!!~~~~~~")
