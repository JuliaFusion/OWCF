####################################################  os2com.jl ####################################################

#### Description:
# This script will transform a quantity from (E,pm,Rm) coordinates to 
# (E,Λ,Pϕ_n;σ) coordinates. If the quantity is given per phase-space volume, a Jacobian is needed when transforming. 
# This is taken care of automatically, via automatic differentiation using dual numbers. Λ is the normalized magnetic 
# moment given by Λ=μ*B0/E where μ is the magnetic moment, B0=B(mag. axis) and E is the energy. Pϕ_n is the normalized 
# toroidal canonical momentum, given by Pϕ_n=Pϕ/(q*|Ψ_w|) where q is the charge of the fast ion and Ψ_w is the poloidal 
# mangetic flux at the last closed flux surfaces (LCFS). If Ψ(LCFS)==0 for some reason (e.g. due to some convention),
# Ψ_w=Ψ(mag. axis) is assumed instead. σ is a binary coordinate. σ=-1 (Julia index 1) corresponds to counter-current 
# orbits. σ=+1 (Julia index 2) corresponds to co-current orbits.
#
# Please note that, if os2com.jl is used for orbit weight functions, it expects the input
# to be in the 4D inflated format. That is, an output file from orbWeights_2Dto4D.jl or similar.
# In short, 3D and 4D formats are accepted as input. 2D formats are NOT.

#### Inputs (Units given when defined in script)
# Given via input file start_os2com_template.jl, for example. 'template' should be replaced by whatever.

#### Outputs
# -

#### Saved files
# os2com_output_[tokamak]_[FI species]_at[timepoint]s_[length(E_array)]x[length(Lambda_array)]x[length(Pphi_n_array)].jld2
#   The saved file will have different keys depending on the mapped quantity. For example, if a topological
#   map has been mapped from (E,pm,Rm) to (E,Λ,Pϕ_n;σ), the key of the mapped quantity will be
#   'topoMap'. The pertaining Λ- and Pϕ_n-arrays will have the keys 'Lambda_array_topoMap' and 'Pphi_n_array_topoMap'. 
#   If the mapping from (E,pm,Rm) to (E,Λ,Pϕ_n;σ) was instead done for an inflated orbit weight matrix, the keys will be 
#   'Lambda_array_W' and 'Pphi_n_array_W' etc. If the mapping was instead done for a fast-ion distribution, the keys 
#   will be 'Lambda_array_F' and 'Pphi_n_array_F' etc. To see all data keys for a .jld2 file, do 'myfile = jldopen(...)' 
#   followed by 'keys(myfile)'. 
# In addition to the COM keys, the saved file will have the keys
#   E_array - The energy grid points of the (E,Λ,Pϕ_n;σ) grid (and the (E,pm,Rm) grid) - Array{Float64,1}
#   pm_array - The pitch maximum grid points of the (E,pm,Rm) grid - Array{Float64,1}
#   Rm_array - The radius maximum grid points of the (E,pm,Rm) grid - Array{Float64,1}
#   FI_species - The fast-ion species. E.g. 'D', 'T', '3he' etc - String
#   filepath_equil - The path to the magnetic equilibrium file used for the mapping - String
#   tokamak - The tokamak for the mapped quantity - String
#   extra_kw_args - A dictionary with the extra keyword arguments used for orbit integration - Dict
# In addition, if defined, the file will also have the key
#   TRANSP_id - The TRANSP run identification number of the quantities mapped - String

### Other
#

# Script written by Henrik Järleblad. Last maintained 2025-01-22.
######################################################################################################################

## ---------------------------------------------------------------------------------------------
verbose && println("Loading Julia packages... ")
@everywhere begin
    using EFIT # For calculating magnetic equilibrium quantities
    using Equilibrium # For loading flux function data, tokamak geometry data etc.
    using GuidingCenterOrbits # For calculating guiding-center orbits
    using OrbitTomography # To compute orbit grid structs
    using ProgressMeter # To display computational progress during parallel computations
    using JLD2 # To write/open .jld2 files (Julia data files, basically)
    using FileIO # To write/open files in general
    include(folderpath_OWCF*"misc/species_func.jl") # To convert species labels to particle mass
    include(folderpath_OWCF*"misc/availReacts.jl") # To check reaction availability and extract fast-ion and thermal species
    include(folderpath_OWCF*"extra/dependencies.jl") # To be able to use orbit_grid() without progress bar
end
## ---------------------------------------------------------------------------------------------

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
    timepoint = (timepoint == nothing ? XX*","*YYYY : timepoint) # Format XX,YYYY to avoid "." when including in filename of saved output
else # Otherwise, assume magnetic equilibrium is a saved .jld2 file
    myfile = jldopen(filepath_equil,false,false,false,IOStream)
    M = myfile["S"]
    wall = myfile["wall"]
    close(myfile)
    jdotb = (M.sigma_B0)*(M.sigma_Ip)
    timepoint = (timepoint == nothing ? "00,0000" : timepoint)
end
## ---------------------------------------------------------------------------------------------

## ---------------------------------------------------------------------------------------------
# Load orbit grid vectors and quantities to be mapped from (E,pm,Rm) to (E,Λ,Pϕ_n;σ) filepath_Q
analyticalOWs = false # Assume by default that, if orbit weight functions, the orbit weight functions are not analytical
Q_dict = Dict() # A dictionary containing all the quantities to be mapped
Q_message = "" # A message that will be printed informing the user of the quantitie(s) to be mapped
if isfile(filepath_Q)
    verbose && println("Loading from data file... ")
    myfile = jldopen(filepath_Q,false,false,false,IOStream)
    E_array = collect(myfile["E_array"])
    pm_array = collect(myfile["pm_array"])
    Rm_array = collect(myfile["Rm_array"])
    if haskey(myfile,"W")
        Ed_array = myfile["Ed_array"]
        Q_dict["W"] = myfile["W"]
        Q_message *= "An (E,pm,Rm) orbit weight matrix of size $(size(Q_dict["W"])) will be mapped to (E,Λ,Pϕ_n;σ) constants-of-motion space."
    elseif haskey(myfile,"W2D")
        error("It looks like you have provided an output file from calcOrbWeights.jl as input? Please provide an output file from orbWeights_2Dto4D.jl as input instead!")
    elseif haskey(myfile,"topoMap")
        Q_dict["topoMap"] = myfile["topoMap"]
        Q_message *= "An (E,pm,Rm) topological map of size $(size(Q_dict["topoMap"])) will be mapped to (E,Λ,Pϕ_n;σ) constants-of-motion space."
        if haskey(myfile,"polTransTimes")
            Q_dict["polTransTimes"] = myfile["polTransTimes"]
            Q_dict["torTransTimes"] = myfile["torTransTimes"]
            Q_message *= "\n Maps of poloidal and toroidal transit times will be included."
        end
    elseif haskey(myfile,"topoBounds")
        Q_dict["topoBounds"] = myfile["topoBounds"]
        Q_message *= "Topological boundaries in (E,pm,Rm) will be mapped to (E,Λ,Pϕ_n;σ) constants-of-motion space."
    elseif haskey(myfile,"nullOrbs_indices")
        Ed_array = myfile["Ed_array"]
        nullOrbs = zeros(length(Ed_array), length(E_array), length(pm_array), length(Rm_array))
        nullOrbs[myfile["nullOrbs_indices"]] = 1.0
        Q_dict["nullOrbs"] = nullOrbs
        Q_message *= "Null orbits in (E,pm,Rm) will be mapped to (E,Λ,Pϕ_n;σ) constants-of-motion space."
    elseif haskey(myfile,"F_os_3D")
        Q_dict["F"] = myfile["F_os_3D"]
        Q_message *= "An (E,pm,Rm) fast-ion distribution will be transformed to (E,Λ,Pϕ_n;σ) constants-of-motion space."
    else
        error("filepath_Q did not contain weights functions, topological map, topological boundaries or indices for null orbits (via include2Dto4D when executing extractNullOrbits.jl). Please correct and re-try.")
    end
    if haskey(myfile,"reaction")
        reaction = myfile["reaction"]
    end
    if haskey(myfile,"reaction_full")
        reaction_full = myfile["reaction_full"]
    end
    if haskey(myfile,"FI_species")
        FI_species = myfile["FI_species"]
    end
    if haskey(myfile,"extra_kw_args")
        extra_kw_args = myfile["extra_kw_args"]
    end
    if haskey(myfile,"analyticalOWs")
        analyticalOWs = true
    end
    close(myfile)
else
    error("filepath_Q not specified correctly. Please correct and re-try.")
end
## ---------------------------------------------------------------------------------------------

## ---------------------------------------------------------------------------------------------
# Determine fast-ion species from reaction
verbose && println("Identifying fast-ion species... ")
if (@isdefined reaction_full)
    thermal_species, FI_species = checkReaction(reaction_full; verbose=verbose, projVelocity=analyticalOWs)
    @everywhere FI_species = $FI_species # Transfer variable to all external processes
elseif (@isdefined reaction)
    FI_species = (split(reaction,"-"))[1] # Assume first species specified in reaction to be the fast-ion species. For example, in 'p-t' the 'p' will be assumed the thermal species.
    @everywhere FI_species = $FI_species # Transfer variable to all external processes
else
    # Backwards compatibility with output files from older versions of ps2os.jl and F_os_1Dto3D.jl
    if !(@isdefined FI_species)
        # Must be the case of output file from older version of ps2os.jl or F_os_1Dto3D.jl
        # Try to determine FI_species from filename
        if Sys.iswindows()
            filename_Q = split(filepath_Q,"\\")[end] # Special Windows case with backwards slash
        else
            filename_Q = split(filepath_Q,"/")[end]
        end
        FI_species = split(filename_Q,"_")[end-1] # For both F_os_1Dto3D.jl and ps2os.jl output, the FI_species should be found second to last in the filename
    end
    FI_species = FI_species # Already defined
end
## ---------------------------------------------------------------------------------------------

## ---------------------------------------------------------------------------------------------
# Printing script info and inputs
println("")
println("-------------------------------------------------os2com.jl-------------------------------------------------")
println("Fast-ion species: "*FI_species)
print("Tokamak: "*tokamak)
print("          ")
if @isdefined TRANSP_id
    println("TRANSP ID: "*TRANSP_id)
end
println("")
print("Equilibrium and geometry file specified: ")
println(filepath_equil)
println("")
println(Q_message)
println("")
println("(E,pm,Rm) grid has dimensions $(length(E_array))x$(length(pm_array))x$(length(Rm_array))")
println("(E,Λ,Pϕ_n;σ) grid will have same dimensions $(length(E_array))x$(nmu)x$(nPphi)x2.")
println("")
println("Extra orbit integration algorithm keywords specified: ")
println(extra_kw_args)
println("")
println("Please remove previously saved files with the same file name (if any) prior to script completion. Quickly!")
println("")
println("If you would like to change any settings, please edit the start_os2com_template.jl file or similar.")
println("")
println("Written by Henrik Järleblad. Last maintained 2025-01-22.")
println("------------------------------------------------------------------------------------------------------------------")
println("")
## ---------------------------------------------------------------------------------------------

## ---------------------------------------------------------------------------------------------
Q_dict_COM = Dict()
for key in keys(Q_dict)
    Q_dict_COM[key] = nothing # Initialize hashmap element
end

for key in keys(Q_dict)
    verbose && println("Mapping "*key*" to (E,Λ,Pϕ_n;σ)... ")
    if key=="W"
        good_coords = nothing # We don't know anything...
    elseif key=="F"
        good_coords = nothing # We don't know anything...
    elseif key=="topoMap"
        good_coords = findall(x-> x!=9.0,Q_dict[key])
    elseif key=="polTransTimes"
        good_coords = findall(x-> (x!=9.0) && (x!=7.0), Q_dict["topoMap"]) # We know the key 'topoMap' must exist, if 'polTransTimes' exists
    elseif key=="torTransTimes"
        good_coords = findall(x-> (x!=9.0) && (x!=7.0), Q_dict["topoMap"]) # We know the key 'topoMap' must exist, if 'torTransTimes' exists
    elseif key=="topoBounds"
        good_coords = findall(x-> x==1.0, Q_dict[key])
    elseif key=="nullOrbs"
        good_coords = nothing # We don't know anything. Every orbit weight function will likely have different null orbits
    else
        error("This should be impossible to reach, given the load check above.")
    end
    Q_dict_COM[key], E_array[:], Q_dict_COM["Lambda_array_"*key], Q_dict_COM["Pphi_n_array_"*key] = os2COM(M, Q_dict[key], Vector(E_array), Vector(pm_array), Vector(Rm_array), FI_species; nl=nmu, npp=nPphi, isTopoMap=(key=="topoMap" ? true : false), needJac=(key=="F" ? true : false), verbose=verbose, good_coords=good_coords, wall=wall, extra_kw_args=extra_kw_args)
end
## ---------------------------------------------------------------------------------------------

## ---------------------------------------------------------------------------------------------
verbose && println("Success! Saving to file... ")
global filepath_output_orig = folderpath_o*"os2com_output_"*tokamak*"_"*FI_species*"_at"*timepoint*"s_"*"$(length(E_array))x$(nmu)x$(nPphi)x2"
global filepath_output = deepcopy(filepath_output_orig)
global count = 1
while isfile(filepath_output*".jld2") # To take care of not overwriting files. Add _(1), _(2) etc
    global filepath_output = filepath_output_orig*"_($(Int64(count)))"
    global count += 1 # global scope, to surpress warnings
end
global filepath_output = filepath_output*".jld2"
myfile = jldopen(filepath_output,true,true,false,IOStream)
for key in keys(Q_dict_COM)
    write(myfile,key,Q_dict_COM[key])
end
if @isdefined Ed_array
    write(myfile,"Ed_array",Ed_array)
end
write(myfile,"E_array",E_array)
write(myfile,"pm_array",pm_array)
write(myfile,"Rm_array",Rm_array)
if @isdefined FI_species
    write(myfile,"FI_species",FI_species)
end
write(myfile,"filepath_equil",filepath_equil)
write(myfile,"tokamak",tokamak)
if @isdefined TRANSP_id
    write(myfile,"TRANSP_id",TRANSP_id)
end
write(myfile,"extra_kw_args",extra_kw_args)
close(myfile)

verbose && println("~~~os2com.jl finished!~~~")