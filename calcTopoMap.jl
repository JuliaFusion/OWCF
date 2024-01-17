############################################### calcTopoMap.jl #########################################

#### Description:
# This script computes the topological map for an orbit grid. It creates a 3D
# matrix where each element in the matrix corresponds to one grid point in orbit space
# (with the same indices etc.). Each matrix element will have one of the following integer
# values:
# 1 - This integer value is given to all stagnation orbits. Red.
# 2 - This integer value is given to all trapped orbits. Blue.
# 3 - This integer value is given to all co-passing orbits. Green.
# 4 - This integer value is given to all ctr-passing orbits. Purple.
# 5 - This integer value is given to all potato orbits. Orange.
# 8 - This integer value is given to all counter-stagnation orbits. Pink.
# 9 - This integer value is given to all invalid, incomplete and lost orbits. Gray.
#
# If distinguishLost==true, then all lost orbits will instead be given the
# integer value 7, which corresponds to the color brown.
#
# If distinguishIncomplete==true, then all incomplete orbits will instead be given the
# integer value 6, which corresponds to the color yellow.
#
# The color scheme used is the 'Set1_9' and can be found at: https://docs.juliaplots.org/latest/generated/colorschemes/
#
# If isfile(filepath_W), make sure to provide the correct path to the weights file.
# If the weights file is a .jld2-file, make sure to set weightsFileJLD2 to true.
# Please run, for example, calcOrbWeights.jl prior to running this script, to acquire
# a file with orbit weights in .jld2-format.

#### Inputs (units given when defined in the script)
# Given in start_calcTopoMap_template.jl

#### Outputs
# -

#### Saved files
# topoMap_[tokamak]_[TRANSP_id]_at[timepoint]s_[FI species]_[nE]x[npm]x[nRm].jld2
# If the input distinguishIncomplete was set to true, the output file will instead have the name
#   topoMap_[tokamak]_[TRANSP_id]_at[timepoint]s_[FI species]_[nE]x[npm]x[nRm]_wIncomp.jld2
# If the input distinguishLost was set to true, the output file will instead have the name
#   topoMap_[tokamak]_[TRANSP_id]_at[timepoint]s_[FI species]_[nE]x[npm]x[nRm]_wLost.jld2
# If both distinguishIncomplete and distinguishLost were set to true, the output file will instead have the name
#   topoMap_[tokamak]_[TRANSP_id]_at[timepoint]s_[FI species]_[nE]x[npm]x[nRm]_wIncomp_wLost.jld2
# The saved file will have the fields:
#   topoMap - The saved topological map. Dimensions are size(orbit grid) - Array{Float64,3}
#   E_array - The fast-ion energy grid array used for orbit space - Array{Float64,1}
#   pm_array - The fast-ion pm grid array used for orbit space - Array{Float64,1}
#   Rm_array - The fast-ion Rm grid array used for orbit space - Array{Float64,1}
#   FI_species - The fast-ion species of the orbit grid - String
#   distinguishIncomplete - The boolean specified as input by the user - Bool
#   distinguishLost - The boolean specified as input by the user - Bool
#   extra_kw_args - The extra input arguments used as keyword arguments for the (relativistic) equations-of-motion integration algorithm - Dictionary
# If the input saveTransitTimeMaps was set to true, the output file will also include the keys
#   polTransTimes - The poloidal transit time for every valid orbit. size(orbit grid) - Array{Float64,3}
#   torTransTimes - The toroidal transit time for every valid orbit. size(orbit grid) - Array{Float64,3}

#### Other
# Please note, because of numerical reasons, all integers (1,2,3,4,5,6,7,8,9) will actually be saved as 
# their Float64 counterparts (1.0, 2.0 etc).

# Script written by Henrik Järleblad. Last maintained 2023-01-09.
#######################################################################################################

## ------
# Loading packages
verbose && println("Loading packages (usually takes approx 30 secs/processor)... ")
@everywhere begin
    using EFIT
    using Equilibrium
    using GuidingCenterOrbits
    using JLD2
    using FileIO
    using ProgressMeter
    using Base.Iterators
    using HDF5
    using NetCDF
    include("misc/species_func.jl")
    include("misc/availReacts.jl") # To check reaction availability and extract fast-ion and thermal species
end

## ------
# Loading tokamak equilibrium
verbose && println("Loading tokamak equilibrium... ")
if ((split(filepath_equil,"."))[end] == "eqdsk") || ((split(filepath_equil,"."))[end] == "geqdsk")
    M, wall = read_geqdsk(filepath_equil,clockwise_phi=false) # Assume counter-clockwise phi-direction
    jdotb = M.sigma # The sign of the dot product between the plasma current and the magnetic field

    # Extract timepoint information from .eqdsk/.geqdsk file
    eqdsk_array = split(filepath_equil,".")
    if length(eqdsk_array)>2
        XX = (split(eqdsk_array[end-2],"-"))[end] # Assume format ...-XX.YYYY.eqdsk where XX are the seconds and YYYY are the decimals
        YYYY = eqdsk_array[end-1] # Assume format ...-XX.YYYY.eqdsk where XX are the seconds and YYYY are the decimals
        timepoint = XX*","*YYYY # Format XX,YYYY to avoid "." when including in filename of saved output
    else
        timepoint = "00,0000"
    end
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

## ------
# Defining orbit space basis vectors
analyticalOWs = false # Define this for its own sake
if isfile(filepath_W)
    verbose && println("Loading orbit-space grid points from weight file... ")
    if weightsFileJLD2
        myfile = jldopen(filepath_W,false,false,false,IOStream)
        E_array = myfile["E_array"]
        pm_array = myfile["pm_array"]
        Rm_array = myfile["Rm_array"]
        if haskey(myfile,"reaction")
            reaction = myfile["reaction"]
        end
        if haskey(myfile,"reaction_full")
            reaction_full = myfile["reaction_full"]
        end
        if haskey(myfile, "FI_species")
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
        myfile = h5open(filepath_W)
        DIAGNOSTIC = myfile[keyname_diagHDF5]
        E_array = read(DIAGNOSTIC["E_array"])
        pm_array = read(DIAGNOSTIC["pm_array"])
        Rm_array = read(DIAGNOSTIC["Rm_array"])
        close(myfile)
    end
elseif isfile(filepath_OG)
    verbose && println("Loading orbit-space grid points from .jld2 orbit grid file... ")
    myfile = jldopen(filepath_OG,false,false,false,IOStream)
    og = myfile["og"]
    extra_kw_args = myfile["extra_kw_args"]
    close(myfile)
    E_array = og.energy
    pm_array = og.pitch
    Rm_array = og.r
    og = nothing # Remove unnecessary memory usage
else
    if (E_array) == nothing
        E_array = collect(range(Emin,stop=Emax,length=nE))
    end
    if (pm_array) == nothing
        pm_array = collect(range(pm_min,stop=pm_max,length=npm))
    end
    if (Rm_array) == nothing
        if (Rm_min==nothing) || (Rm_max==nothing)
            if inclPrideRockOrbs
                # 4/5 of the distance from the HFS wall to the magnetic axis is usually enough to capture all the Pride Rock orbits
                Rm_array = range((4*magnetic_axis(M)[1]+minimum(wall.r))/5, stop=maximum(wall.r), length=nRm)
            else
                Rm_array = range(magnetic_axis(M)[1], stop=maximum(wall.r), length=nRm)
            end
        else
            Rm_array = range(Rm_min, stop=Rm_max, length=nRm)
        end
    end
end

## ---------------------------------------------------------------------------------------------
# Determine fast-ion plasma species from reaction
if (@isdefined reaction_full)
    thermal_species, FI_species = checkReaction(reaction_full; verbose=verbose, projVelocity=analyticalOWs)
    @everywhere FI_species = $FI_species # Transfer variable to all external processes
elseif (@isdefined reaction)
    FI_species = (split(reaction,"-"))[2] # Assume first species specified in reaction to be the fast-ion species. For example, in 'p-t' the 'p' will be assumed the thermal species.
    @everywhere FI_species = $FI_species # Transfer variable to all external processes
else
    FI_species = FI_species # Already defined
end

## ------
# Printing script info and inputs
println("")
println("----------------------------------------calcTopoMap.jl----------------------------------------")
print("Tokamak specified: "*tokamak*"        ")
print("TRANSP_id specified: "*TRANSP_id*"        ")
println("Timepoint: $(timepoint)")
println("Magnetic equilibrium file specified: ")
println("--- "*filepath_equil)
if (@isdefined reaction)
    println("")
    println("Fusion reaction: "*reaction)
    println("")
elseif (@isdefined reaction_full)
    println("")
    println("Fusion reaction: "*reaction_full)
    println("")
else
end
if isfile(filepath_W)
    println("Weights file specified: ")
    println("--- "*filepath_W)
elseif isfile(filepath_OG)
    println("Orbit-grid file specified: ")
    println("--- "*filepath_OG)
else
end
println("")
if distributed
    println("$(nprocs()) processors will be used for parallel computing.")
    println("")
else
    println("Single-threaded computing the topological map.")
    println("")
end
if distinguishLost
    println("Lost orbits will be given the integer 7 (brown).")
    println("")
end
if distinguishIncomplete
    println("Incomplete orbits will be given the integer 6 (yellow).")
    println("")
end
println("Topological map will be computed for a $(length(E_array))x$(length(pm_array))x$(length(Rm_array)) orbit-space grid with")
println("--- Energy: [$(minimum(E_array)),$(maximum(E_array))] keV")
println("--- Pitch maximum: [$(minimum(pm_array)),$(maximum(pm_array))]")
if (Rm_min==nothing) || (Rm_max==nothing)
    if inclPrideRockOrbs
        println("--- Radius maximum: [$((4*magnetic_axis(M)[1]+minimum(wall.r))/5),$(maximum(wall.r))] m")
    else
        println("--- Radius maximum: [$(magnetic_axis(M)[1]),$(maximum(wall.r))] m")
    end
else
    println("--- Radius maximum: [$(Rm_min),$(Rm_max)] m")
end
println("")
println("Fast-ion species specified: "*FI_species)
println("")
if saveTransitTimeMaps
    println("Poloidal and toroidal transit times will be saved as .jld2-file keys 'polTransTimes' and 'torTransTimes'.")
    println("")
end
if includeExtractTopoBounds
    println("Boundaries between topological regions will be extracted (and saved) upon completition.")
    println("")
end
println("Extra keyword arguments specified: ")
println(extra_kw_args)
println("")
println("If you would like to change any settings, please edit the start_calcTopoMap_template.jl file.")
println("Written by Henrik Järleblad. Last maintained 2023-01-09.")
println("-----------------------------------------------------------------------------------------------")
println("")

## ------
# Calculating the topological map
npoints = length(pm_array)*length(Rm_array) # Number of points for one energy slice
pm_rep = repeat(pm_array,inner=length(Rm_array)) # To get all points
Rm_rep = repeat(Rm_array,outer=length(pm_array)) # To get all points
pmRm_array = collect(take(zip(pm_rep,Rm_rep),length(pm_rep))) # All points in an energy slice, in one long vector
topoMap_tottot = zeros(length(E_array),length(pm_array),length(Rm_array)) # The total 3D topological map
polTimeMap_tottot = zeros(length(E_array),length(pm_array),length(Rm_array)) # The total 3D poloidal transit time map
torTimeMap_tottot = zeros(length(E_array),length(pm_array),length(Rm_array)) # The total 3D toroidal transit time map
for Ei=1:length(E_array) ###########################################
E = E_array[Ei]
verbose && println("Creating topological map for: $(round(E,digits=2)) keV... ")
if distributed
    if visualizeProgress # if you want the progress to be visualized...
        p = Progress(npoints) # Define the progress bar
        channel = RemoteChannel(()->Channel{Bool}(npoints),1) # Define the channel from which the progress bar draws data
        dataMap_tot = fetch(@sync begin # Start the distributed computational process, fetch result when done
            @async while take!(channel) # An asynchronous process, with no need for sync, since it simply displays the progress bar
                ProgressMeter.next!(p)
            end
            @async begin # No internal syncronization needed here either, only external sync needed
                dataMap = @distributed (+) for i=1:length(pmRm_array) # Compute one result, and reduce (add) it to a resulting matrix
                    pmRm = pmRm_array[i]
                    pm = pmRm[1]
                    Rm = pmRm[2]
                    topoMap_i = zeros(length(pm_array),length(Rm_array))
                    polTimeMap_i = zeros(length(pm_array),length(Rm_array))
                    torTimeMap_i = zeros(length(pm_array),length(Rm_array))
                    EPRc = EPRCoordinate(M,E,pm,Rm;amu=getSpeciesAmu(FI_species), q=getSpeciesEcu(FI_species))
                    if false # for de-bugging purposes. Should it be needed. If so, change to 'true'.
                        o = get_orbit(M,EPRc;wall=wall,store_path=false, verbose=verbose, extra_kw_args...)
                    else
                        o = get_orbit(M,EPRc;wall=wall,store_path=false, extra_kw_args...)
                    end
                    pmi = first(findall(x-> x==pm,pm_array))
                    Rmi = first(findall(x-> x==Rm,Rm_array))
                    polTimeMap_i[pmi,Rmi] = o.tau_p
                    torTimeMap_i[pmi,Rmi] = o.tau_t
                    if (o.class == :lost) && distinguishLost
                        topoMap_i[pmi,Rmi] = 7
                        polTimeMap_i[pmi,Rmi] = 0.0
                        torTimeMap_i[pmi,Rmi] = 0.0
                    elseif (o.class == :incomplete) && distinguishIncomplete
                        topoMap_i[pmi,Rmi] = 6
                        polTimeMap_i[pmi,Rmi] = 0.0
                        torTimeMap_i[pmi,Rmi] = 0.0
                    elseif o.class == :trapped
                        topoMap_i[pmi,Rmi] = 2
                    elseif o.class == :co_passing
                        topoMap_i[pmi,Rmi] = 3
                    elseif (o.class == :stagnation && o.coordinate.r>=magnetic_axis(M)[1]) # Regular stagnation orbit
                        topoMap_i[pmi,Rmi] = 1
                    elseif (o.class == :stagnation && o.coordinate.r<magnetic_axis(M)[1]) # Counterstagnation orbit
                        topoMap_i[pmi,Rmi] = 8
                    elseif o.class == :potato
                        topoMap_i[pmi,Rmi] = 5
                    elseif o.class == :ctr_passing
                        topoMap_i[pmi,Rmi] = 4
                    else
                        topoMap_i[pmi,Rmi] = 9
                        polTimeMap_i[pmi,Rmi] = 0.0
                        torTimeMap_i[pmi,Rmi] = 0.0
                    end

                    dataMap_i = zeros(3,length(pm_array),length(Rm_array))
                    dataMap_i[1,:,:] = topoMap_i
                    dataMap_i[2,:,:] = polTimeMap_i
                    dataMap_i[3,:,:] = torTimeMap_i
                    put!(channel,true) # Update the progress bar
                    dataMap_i # Declare dataMap_i as result to add to topoMap
                end
                put!(channel,false) # Update progress bar
                dataMap # Delcare dataMap as done/result, so it can be fetched
            end
        end)
    else
        dataMap_tot = @distributed (+) for i=1:length(pmRm_array)
            pmRm = pmRm_array[i]
            pm = pmRm[1]
            Rm = pmRm[2]
            topoMap_i = zeros(length(pm_array),length(Rm_array))
            polTimeMap_i = zeros(length(pm_array),length(Rm_array))
            torTimeMap_i = zeros(length(pm_array),length(Rm_array))
            EPRc = EPRCoordinate(M,E,pm,Rm;amu=getSpeciesAmu(FI_species), q=getSpeciesEcu(FI_species))
            if false # for de-bugging purposes. Should it be needed. If so, change to 'true'.
                o = get_orbit(M,EPRc;wall=wall,store_path=false, verbose=verbose, extra_kw_args...)
            else
                o = get_orbit(M,EPRc;wall=wall,store_path=false, extra_kw_args...)
            end
            pmi = first(findall(x-> x==pm,pm_array))
            Rmi = first(findall(x-> x==Rm,Rm_array))
            polTimeMap_i[pmi,Rmi] = o.tau_p
            torTimeMap_i[pmi,Rmi] = o.tau_t
            if (o.class == :lost) && distinguishLost
                topoMap_i[pmi,Rmi] = 7
                polTimeMap_i[pmi,Rmi] = 0.0
                torTimeMap_i[pmi,Rmi] = 0.0
            elseif (o.class == :incomplete) && distinguishIncomplete
                topoMap_i[pmi,Rmi] = 6
                polTimeMap_i[pmi,Rmi] = 0.0
                torTimeMap_i[pmi,Rmi] = 0.0
            elseif o.class == :trapped
                topoMap_i[pmi,Rmi] = 2
            elseif o.class == :co_passing
                topoMap_i[pmi,Rmi] = 3
            elseif (o.class == :stagnation && o.coordinate.r>=magnetic_axis(M)[1]) # Regular stagnation orbit
                topoMap_i[pmi,Rmi] = 1
            elseif (o.class == :stagnation && o.coordinate.r<magnetic_axis(M)[1]) # Counterstagnation orbit
                topoMap_i[pmi,Rmi] = 8
            elseif o.class == :potato
                topoMap_i[pmi,Rmi] = 5
            elseif o.class == :ctr_passing
                topoMap_i[pmi,Rmi] = 4
            else
                topoMap_i[pmi,Rmi] = 9
                polTimeMap_i[pmi,Rmi] = 0.0
                torTimeMap_i[pmi,Rmi] = 0.0
            end

            dataMap_i = zeros(3,length(pm_array),length(Rm_array))
            dataMap_i[1,:,:] = topoMap_i
            dataMap_i[2,:,:] = polTimeMap_i
            dataMap_i[3,:,:] = torTimeMap_i
            dataMap_i # Declare topoMap_i as result to add to dataMap_tot
        end
    end
else # ... good luck
    topoMap_i = zeros(length(pm_array),length(Rm_array))
    polTimeMap_i = zeros(length(pm_array),length(Rm_array))
    torTimeMap_i = zeros(length(pm_array),length(Rm_array))
    count = 1
    for pmRm in pmRm_array # Compute one result, and reduce (add) it to a resulting matrix
        pm = pmRm[1]
        Rm = pmRm[2]
        verbose && println("$(count)/$(length(pmRm_array)):    (E,pm,Rm) = ($(E),$(pm),$(Rm))")
        EPRc = EPRCoordinate(M,E,pm,Rm;amu=getSpeciesAmu(FI_species), q=getSpeciesEcu(FI_species))
        if false # for de-bugging purposes. Should it be needed. If so, change to 'true'.
            o = get_orbit(M,EPRc;wall=wall,store_path=false, verbose=verbose, extra_kw_args...)
        else
            o = get_orbit(M,EPRc;wall=wall,store_path=false, extra_kw_args...)
        end
        pmi = first(findall(x-> x==pm,pm_array))
        Rmi = first(findall(x-> x==Rm,Rm_array))
        polTimeMap_i[pmi,Rmi] = o.tau_p
        torTimeMap_i[pmi,Rmi] = o.tau_t
        if (o.class == :lost) && distinguishLost
            topoMap_i[pmi,Rmi] = 7
            polTimeMap_i[pmi,Rmi] = 0.0
            torTimeMap_i[pmi,Rmi] = 0.0
        elseif (o.class == :incomplete) && distinguishIncomplete
            topoMap_i[pmi,Rmi] = 6
            polTimeMap_i[pmi,Rmi] = 0.0
            torTimeMap_i[pmi,Rmi] = 0.0
        elseif o.class == :trapped
            topoMap_i[pmi,Rmi] = 2
        elseif o.class == :co_passing
            topoMap_i[pmi,Rmi] = 3
        elseif (o.class == :stagnation && o.coordinate.r>=magnetic_axis(M)[1]) # Regular stagnation orbit
            topoMap_i[pmi,Rmi] = 1
        elseif (o.class == :stagnation && o.coordinate.r<magnetic_axis(M)[1]) # Counterstagnation orbit
            topoMap_i[pmi,Rmi] = 8
        elseif o.class == :potato
            topoMap_i[pmi,Rmi] = 5
        elseif o.class == :ctr_passing
            topoMap_i[pmi,Rmi] = 4
        else
            topoMap_i[pmi,Rmi] = 9
            polTimeMap_i[pmi,Rmi] = 0.0
            torTimeMap_i[pmi,Rmi] = 0.0
        end
        count = count + 1
    end
    dataMap_tot = zeros(3,length(pm_array),length(Rm_array))
    dataMap_tot[1,:,:] = topoMap_i
    dataMap_tot[2,:,:] = polTimeMap_i
    dataMap_tot[3,:,:] = torTimeMap_i
end

global topoMap_tottot[Ei,:,:] = dataMap_tot[1,:,:] # One energy slice of topological map
global polTimeMap_tottot[Ei,:,:] = dataMap_tot[2,:,:] # One energy slice of poloidal transit times
global torTimeMap_tottot[Ei,:,:] = dataMap_tot[3,:,:] # One energy slice of toroidal transit times
end ########################################### <----------- PLEASE NOTE THE 'end'
## ------
# Save the results
verbose && println("Saving... ")
ident = ""
if distinguishIncomplete
    ident *= "_wIncomp"
end
if distinguishLost
    ident *= "_wLost"
end
# global declaration to include extractTopoBounds.jl scope
global filepath_tm_orig = folderpath_o*"topoMap_"*tokamak*"_"*TRANSP_id*"_at"*timepoint*"s_"*FI_species*"_$(length(E_array))x$(length(pm_array))x$(length(Rm_array))"*ident
global filepath_tm = deepcopy(filepath_tm_orig)
global count = 1
while isfile(filepath_tm*".jld2") # To take care of not overwriting files. Add _(1), _(2) etc
    global filepath_tm = filepath_tm_orig*"_($(Int64(count)))"
    global count += 1 # global scope, to surpress warnings
end
global filepath_tm = filepath_tm*".jld2"
myfile = jldopen(filepath_tm,true,true,false,IOStream)
write(myfile,"topoMap",topoMap_tottot)
write(myfile,"E_array",collect(E_array))
write(myfile,"pm_array",collect(pm_array))
write(myfile,"Rm_array",collect(Rm_array))
write(myfile,"FI_species",FI_species)
write(myfile,"distinguishIncomplete",distinguishIncomplete)
write(myfile,"distinguishLost",distinguishLost)
if saveTransitTimeMaps
    write(myfile,"polTransTimes",polTimeMap_tottot)
    write(myfile,"torTransTimes",torTimeMap_tottot)
end
write(myfile,"extra_kw_args",extra_kw_args)
close(myfile)

if includeExtractTopoBounds
    verbose && println("Initiating extraction of topoBounds for saved topoMap... ")
    include("helper/extractTopoBounds.jl")
end
println("~~~~~~ calcTopoMap.jl done! ~~~~~~")
