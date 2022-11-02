################################ extractTopoBounds.jl #########################################

#### Description:
# This script takes a file with a computed topological map and extracts the boundaries of the
# topological regions. The topological map is a three-dimensional matrix, and the result will be
# a three-dimensional matrix with the same dimensions. The boundaries are marked with the
# Float64 '1.0', while the other elements of the matrix are set to Float64 '0.0'.
# Please run calcTopoMap.jl prior to running this script, to acquire a topological map saved in
# .jld2-format.

# Prior to running this script, please make sure you have run calcTopoMap.jl, and noted
# the path to the saved file (filepath_tm).

#### Inputs (units given when defined in the script):
# folderpath_OWCF - The path to the OWCF folder on your computer - String
# filepath_tm - The path to the .jld2-file containing the topological map (and more) - String
# tokamak - The tokamak for which the topological map is created - String
# TRANSP_id - The TRANSP_id of the tokamak shot - Int64
# timepoint - The tokamak shot timepoint of the magnetic equilibrium. Format XX,YYYY where XX are seconds and YYYY are decimals - String
# folderpath_o - The path to the desired output folder - String
# FI_species - The fast-ion species of the orbits. "D", "T" etc - String
# verbose - If true, then the script will talk a lot! - Bool


#### Outputs
# -

#### Saved files
# topoBounds_[tokamak]_[TRANSP_id]_at[timepoint]s_[nE]x[npm]x[nRm].jld2
#   This saved file will have the fields:
#   topoBounds - The saved topological bounds. Dimensions are size(orbit grid) - Array{Float64,3}
#   E_array - The fast-ion energy grid array used for orbit space - Array{Float64,1}
#   pm_array - The fast-ion pm grid array used for orbit space - Array{Float64,1}
#   Rm_array - The fast-ion Rm grid array used for orbit space - Array{Float64,1}
# If extra_kw_args is defined, then the saved file will also have the key
#   extra_kw_args - The keyword arguments fed into the equations-of-motion integration algorithm - Dictionary

# Script written by Henrik Järleblad. Last maintained 2022-09-01.
###############################################################################################

## ------
# Inputs
partOfCalcTopoMap = true # Set this to false if you are running this script independently. Remember to reset it to 'true' afterwards
if !(partOfCalcTopoMap) # If this is run as a standalone script, please specify inputs below, and set partOfcalcTopoMap to false
    folderpath_OWCF = ""
    filepath_tm = "" # The name of the .jld2 file containing the topological map
    timepoint = nothing # Format "XX,YYYY" where XX are seconds and YYYY are decimals. If unknown, left as nothing. The algorithm will try to figure it out automatically
    tokamak = "" # The name of the tokamak. Leave unspecified if unknown/not applicable
    TRANSP_id = "00000X00" # The TRANSP shot id. Leave unspecified if unknown/not applicable
    folderpath_o = "" # The path to the desired output folder
    FI_species = "" # D, t, p etc
    verbose = true

    # These three lines are only executed if extractTopoBounds.jl is not executed as a part of calcTopoMap.jl
    cd(folderpath_OWCF)
    using Pkg
    Pkg.activate(".")
end

## ------
# Loading packages
verbose && println("Loading Julia packages... ")
using JLD2
using FileIO

verbose && println("Loading topological map file... ")
myfile = jldopen(filepath_tm,false,false,false,IOStream)
topoMap = myfile["topoMap"]
E_array = myfile["E_array"]
pm_array = myfile["pm_array"]
Rm_array = myfile["Rm_array"]
distinguishIncomplete = myfile["distinguishIncomplete"]
distinguishLost = myfile["distinguishLost"]
if haskey(myfile,"extra_kw_args")
    extra_kw_args = myfile["extra_kw_args"]
end
close(myfile)

## ------
# Try to figure out timepoint from input data
if (timepoint==nothing)
    filepath_tm_array = split(filepath_tm,"_")
    time_string = filepath_tm_array[end-3] # Assume this to start with
    if (time_string[1]=='a') && (time_string[2]=='t') && (time_string[end]=='s')
        # Correct format identified. Superseeds user input and/or timepoint found via equilibrium file
        timepoint = time_string[3:end-1] # First to elements are "at" and last element is "s"
    else
        timepoint = "00,0000"
    end
end

## ------
# Printing script info and inputs
println("")
println("----------------------------------------extractTopoBounds.jl----------------------------------------")
print("Tokamak specified: "*tokamak*"     ")
println("TRANSP_id specified: $(TRANSP_id)     Timepoint: "*timepoint*"       FI_species: "*FI_species)
println("Topological map (.jld2) file specified: ")
println(filepath_tm)
println("")
println("Topological boundaries will be extracted for a $(length(E_array))x$(length(pm_array))x$(length(Rm_array)) orbit-space grid with")
println("Energy: [$(minimum(E_array)),$(maximum(E_array))] keV")
println("Pitch maximum: [$(minimum(pm_array)),$(maximum(pm_array))]")
println("Radius maximum: [$(minimum(Rm_array)),$(maximum(Rm_array))] m")
println("")
println("If you would like to change any settings, please edit the extractTopoBounds.jl file.")
println("")
println("Written by Henrik Järleblad. Last maintained 2022-08-22.")
println("----------------------------------------------------------------------------------------------------")
println("")

## ------
# Compute topological boundaries
# Assume (E,pm,Rm) dimensions order
verbose && println("Extracting topological boundaries... ")
topoBounds = zeros(size(topoMap))
for Ei=1:size(topoMap,1)
for Rmi=1:size(topoMap,3)
    for type_i=[1.0,2.0,3.0,4.0,5.0,8.0] # 1.0=stagnation, 2.0=trapped, 3.0=co-passing, 4.0=ctr-passing, 5.0=potato, 8.0=counter-stagnation
        # first and last element index
        first_ind = findfirst(x-> x==type_i,topoMap[Ei,:,Rmi])
        last_ind = findlast(x-> x==type_i,topoMap[Ei,:,Rmi])
        if !(typeof(first_ind) == Nothing)
            topoBounds[Ei,first_ind,Rmi] = 1.0
            topoBounds[Ei,last_ind,Rmi] = 1.0
        end

        # check whether to mark vertical boundaries
        if !(Rmi==1 || Rmi==size(topoMap,3))
            all_ind_prev = findall(x-> x==type_i,topoMap[Ei,:,Rmi-1])
            all_ind_next = findall(x-> x==type_i,topoMap[Ei,:,Rmi+1])
            all_ind = findall(x-> x==type_i,topoMap[Ei,:,Rmi])
            if((length(all_ind_prev)==0 || length(all_ind_next)==0) && length(all_ind)!=0)
                topoBounds[Ei,all_ind,Rmi] .= 1.0
            end
        end
    end
end
for pmi=1:size(topoMap,2)
    for type_i=[1.0,2.0,3.0,4.0,5.0,8.0] # 1.0=stagnation, 2.0=trapped, 3.0=co-passing, 4.0=ctr-passing, 5.0=potato, 8.0=counter-stagnation
        # first and last element index
        first_ind = findfirst(x-> x==type_i,topoMap[Ei,pmi,:])
        last_ind = findlast(x-> x==type_i,topoMap[Ei,pmi,:])
        if !(typeof(first_ind) == Nothing)
            topoBounds[Ei,pmi,first_ind] = 1.0
            topoBounds[Ei,pmi,last_ind] = 1.0
        end

        # check whether to mark vertical boundaries
        if !(pmi==1 || pmi==size(topoMap,2))
            all_ind_prev = findall(x-> x==type_i,topoMap[Ei,pmi-1,:])
            all_ind_next = findall(x-> x==type_i,topoMap[Ei,pmi+1,:])
            all_ind = findall(x-> x==type_i,topoMap[Ei,pmi,:])
            if((length(all_ind_prev)==0 || length(all_ind_next)==0) && length(all_ind)!=0)
                topoBounds[Ei,pmi,all_ind] .= 1.0
            end
        end
    end
end
end

## ------
# Saving
ident = ""
if distinguishIncomplete
    ident *= "_wIncomp"
end
if distinguishLost
    ident *= "wLost"
end
verbose && println("Saving... ")
timepoint = (timepoint == nothing ? "00,0000" : timepoint)
myfile = jldopen(folderpath_o*"topoBounds_"*tokamak*"_"*TRANSP_id*"_"*FI_species*"_$(length(E_array))x$(length(pm_array))x$(length(Rm_array))"*ident*".jld2",true,true,false,IOStream)
write(myfile,"topoBounds",topoBounds)
write(myfile,"E_array",collect(E_array))
write(myfile,"pm_array",collect(pm_array))
write(myfile,"Rm_array",collect(Rm_array))
if @isdefined extra_kw_args
    write(myfile,"extra_kw_args",extra_kw_args)
end
close(myfile)
println("~~~~~~ extractTopoBounds.jl done! ~~~~~~")
