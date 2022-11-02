#########################################  calcOrbGrid.jl #########################################

#### Description:
# This script will compute an orbit grid and the pertaining valid orbits. It will then save the 
# results to a .jld2 file, which can be easily read and provided as input to other scripts in the 
# OWCF. The orbit grid and valid orbits are computed using the (E,pm,Rm) coordinates, just as
# everything else in the OWCF.

#### Inputs (Units given when defined in script)
# Given via input file start_calcOG_template.jl, for example. 'template' should be replaced by whatever.

#### Outputs
# -

#### Saved files
# orbGrid_[tokamak]_[TRANSP_id]_at[timepoint]s_[FI species]_[nE]x[npm]x[nRm].jld2
# This saved file will have the fields:
# og - The OrbitGrid struct computed by the calcOrbGrid.jl script. Please see OrbitTomography.jl/orbits.jl for more info - OrbitGrid
# og_orbs - A vector containing Orbit structs; the valid orbits for the orbit grid og. Please see GuidingCenterOrbits.jl/orbit.jl for more info - Vector{Orbit}
# FI_species - The fast-ion species of the guiding-center orbits in og_orbs. "D", "T" etc - String
# tokamak - The tokamak of the magnetic equilibrium used to compute the orbits in og and og_orbs - String
# TRANSP_id - The TRANSP identification number of the TRANSP run corresponding to the magnetic equilibrium. If used. Otherwise, it will be "" - String
# filepath_equil - The path to the (.eqdsk/.geqdsk/.jld2) file used as magnetic equilibrium data for the orbit computations - String
# extra_kw_args - The extra input arguments used as keyword arguments for the (relativistic) equations-of-motion integration algorithm - Dictionary

### Other
# Please note, in Julia, a Vector{T} type is the same type as Array{T,1}. That is, they are
# equivalent. Vector{T} is an alias for Array{T,1}, where T can be any type.

# Script written by Henrik Järleblad. Last maintained 2022-09-09.
################################################################################################

## ---------------------------------------------------------------------------------------------
verbose && println("Loading Julia packages... ")
@everywhere begin
    cd(folderpath_OWCF) # Necessary to move all the workers to the correct folder
    using EFIT # For calculating magnetic equilibrium quantities
    using Equilibrium # For loading flux function data, tokamak geometry data etc.
    using GuidingCenterOrbits # For calculating guiding-center orbits
    using OrbitTomography # This is what all this is about!
    using ProgressMeter # To display computational progress during parallel computations
    using JLD2 # To write/open .jld2 files (Julia files, basically)
    using FileIO # To write/open files in general
    include(folderpath_OWCF*"misc/species_func.jl") # To convert species labels to particle mass
    include(folderpath_OWCF*"extra/dependencies.jl") # To be able to use orbit_grid() without progress bar
end

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
        timepoint = "00,0000" # Unknown timepoint for magnetic equilibrium
    end
end
## ---------------------------------------------------------------------------------------------

## ---------------------------------------------------------------------------------------------
# Defining orbit grid vectors
verbose && println("Defining orbit grid vectors... ")
if !(E_array == nothing)
    Emin = minimum(E_array)
    Emax = maximum(E_array)
    nE = length(E_array)
else
    E_array = range(Emin,stop=Emax,length=nE)
end
if !(pm_array == nothing)
    pm_min = minimum(pm_array)
    pm_max = maximum(pm_array)
    npm = length(pm_array)
else
    pm_array = range(pm_min, stop=pm_max, length=npm)
end
if !(Rm_array == nothing)
    Rm_min = minimum(Rm_array)
    Rm_max = maximum(Rm_array)
    nRm = length(Rm_array)
else
    if (Rm_min==nothing) || (Rm_max==nothing)
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
## ---------------------------------------------------------------------------------------------

## ---------------------------------------------------------------------------------------------
# Printing script info and inputs
println("")
println("-------------------------------------------------calcOrbGrid.jl-------------------------------------------------")
print("Fast-ion species: "*FI_species)
print("          ")
print("Tokamak: "*tokamak)
print("          ")
println("TRANSP ID: "*TRANSP_id)
println("")
println("A $(length(E_array))x$(length(pm_array))x$(length(Rm_array)) orbit grid and the pertaining valid orbits will be computed for the following (E,pm,Rm) arrays: ")
println("E-array: [$(minimum(E_array)),$(maximum(E_array))] keV. length(E_array): $(length(E_array))")
println("pm-array: [$(minimum(pm_array)),$(maximum(pm_array))]. length(pm_array): $(length(pm_array))")
println("Rm-array: [$(minimum(Rm_array)),$(maximum(Rm_array))] m. length(Rm_array): $(length(Rm_array))")
println("")
println("Equilibrium and geometry file specified: ")
println(filepath_equil)
println("")
println("Extra orbit integration algorithm keywords specified: ")
println(extra_kw_args)
println("")
println("Please remove previously saved files with the same file name (if any) prior to script completion. Quickly!")
println("")
println("If you would like to change any settings, please edit the start_calcOG_template.jl file or similar.")
println("Written by Henrik Järleblad. Last maintained 2022-09-09.")
println("------------------------------------------------------------------------------------------------------------------")
println("")

## ---------------------------------------------------------------------------------------------
# Calculating orbit grid
verbose && debug && println("")
verbose && println("Calculating the orbit grid... ")
og_orbs, og = orbit_grid(M, visualizeProgress, E_array, pm_array, Rm_array; amu=getSpeciesAmu(FI_species), q=getSpeciesEcu(FI_species), wall=wall, extra_kw_args...)

verbose && println("Success! Saving to file... ")
myfile = jldopen(folderpath_o*"orbGrid_"*tokamak*"_"*TRANSP_id*"_at"*timepoint*"s_"*FI_species*"_$(length(E_array))x$(length(pm_array))x$(length(Rm_array)).jld2",true,true,false,IOStream)
write(myfile,"og",og)
write(myfile,"og_orbs",og_orbs)
write(myfile,"FI_species",FI_species)
write(myfile,"tokamak",tokamak)
write(myfile,"TRANSP_id",TRANSP_id)
write(myfile,"filepath_equil",filepath_equil)
write(myfile,"extra_kw_args",extra_kw_args)
close(myfile)

verbose && println("~~~calcOrbGrid.jl finished!~~~")