################################ calcEpRzTopoMap.jl #########################################

#### Description:
# This script calculates the topological map for an (E,p,R,z) grid. It creates a 4D
# matrix where each element in the matrix corresponds to one grid point in particle space
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
# If you set useDistrFile to true, make sure to provide the correct path to the distribution file.
# If distribution file is in .jld2 file format, set distrFileJLD2 to true.

#### Inputs (units given when defined in the script)
# Given in start_calcEpRzTopoMap_template.jl

#### Outputs
# -

#### Saved files
# EpRzTopoMap_[tokamak]_[TRANSP_id]_at[timepoint]s_[FI species]_[nE]x[np]x[nR]x[nz].jld2
# If distinguishLost==true, the name of the output file will instead be
#   EpRzTopoMap_[tokamak]_[TRANSP_id]_at[timepoint]s_[FI species]_[nE]x[np]x[nR]x[nz]_wLost.jld2
# If distinguishIncomplete==true, the name of the output file will instead be
#   EpRzTopoMap_[tokamak]_[TRANSP_id]_at[timepoint]s_[FI species]_[nE]x[np]x[nR]x[nz]_wIncomp.jld2
# If distinguishLost==distinguishIncomplete==true, the name of the output file will instead be
#   EpRzTopoMap_[tokamak]_[TRANSP_id]_at[timepoint]s_[FI species]_[nE]x[np]x[nR]x[nz]_wIncomp_wLost.jld2
# This saved file will have the fields:
#   topoMap - The saved topological map. Dimensions are size(particle space) - Array{Float64,4}
#   E_array - The fast-ion energy grid array used for particle space - Array{Float64,1}
#   p_array - The fast-ion p grid array used for particle space - Array{Float64,1}
#   R_array - The fast-ion R grid array used for particle space - Array{Float64,1}
#   z_array - The fast-ion z grid array used for particle space - Array{Float64,1}
#   distinguishIncomplete - The boolean specified as input by the user - Bool
#   distinguishLost - The boolean specified as input by the user - Bool
#   FI_species - The fast-ion species. "D", "3he", "T" etc - String
#   filepath_distr - If used, the filepath to the fast-ion distribution .jld2 file - String
# If saveTransitTimeMaps==true
#   polTransTimes - The poloidal transit time for all orbits of the particle-space grid - Array{Float64,4}
#   torTransTimes - The toroidal transit time for all orbits of the particle-space grid - Array{Float64,4}

# Script written by Henrik Järleblad. Last maintained 2022-10-05.
########################################################################################

## ------
# Loading packages
verbose && println("Loading packages (usually takes approx 30 secs/processor)... ")
@everywhere begin
    using EFIT
    using Equilibrium
    using GuidingCenterOrbits
    using JLD2
    using HDF5
    using FileIO
    using ProgressMeter
    using Base.Iterators
    include("misc/species_func.jl")
end

## ------
# Loading tokamak equilibrium
verbose && println("Loading equilibrium... ")
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

## ------
# Defining orbit space basis vectors
if useDistrFile
    verbose && println("Loading from distribution file... ")
    if distrFileJLD2
        myfile = jldopen(filepath_distr,false,false,false,IOStream)
        E_array = myfile["energy"]
        p_array = myfile["pitch"]
        R_array = myfile["R"]
        if haskey(myfile,"Z")
            z_array = myfile["Z"]
        else
            z_array = myfile["z"]
        end
        close(myfile)
    else
        myfile = h5open(filepath_distr)
        E_array = read(myfile["energy"])
        p_array = read(myfile["pitch"])
        R_array = read(myfile["R"])
        if haskey(myfile,"Z")
            z_array = read(myfile["Z"])
        else
            z_array = read(myfile["z"])
        end
        close(myfile)
    end
else
    if (E_array) == nothing
        E_array = collect(range(Emin,stop=Emax,length=nE))
    end
    if (p_array) == nothing
        p_array = collect(range(p_min,stop=p_max,length=np))
    end
    if (R_array) == nothing
        if (R_min==nothing) || (R_max==nothing)
            R_array = range(minimum(wall.r), stop=maximum(wall.r), length=nR)
        else
            R_array = range(R_min, stop=R_max, length=nR)
        end
    end
    if (z_array) == nothing
        if (z_min==nothing) || (z_max==nothing)
            z_array = range(minimum(wall.z), stop=maximum(wall.z), length=nz)
        else
            z_array = range(z_min, stop=z_max, length=nz)
        end
    end
end

## ------
# Convert R and z to meters, if necessary
if maximum(R_array)>100.0
    R_array = R_array ./100
end
if maximum(z_array)>100.0
    z_array = z_array ./100
end

## ------
# Printing script info and inputs
println("")
println("----------------------------------------calcEpRzTopoMap.jl----------------------------------------")
print("Tokamak specified: "*tokamak*"        ")
print("TRANSP_id specified: "*TRANSP_id*"        ")
println("Fast-ion species specified: "*FI_species)
println("Equilibrium file specified: ")
println("--- "*filepath_equil)
if filepath_distr != ""
    println("Distribution file specified: ")
    println("--- "*filepath_distr)
end
println("")
if distributed
    println("$(nprocs()) processors will be used for parallel computing.")
    println("")
else
    println("Single-threaded computing the EpRz topological map...")
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
println("Topological map will be computed for a $(length(E_array))x$(length(p_array))x$(length(R_array))x$(length(z_array)) particle-space grid with")
println("--- Energy: [$(minimum(E_array)),$(maximum(E_array))] keV")
println("--- Pitch: [$(minimum(p_array)),$(maximum(p_array))]")
println("--- R: [$(minimum(R_array)),$(maximum(R_array))] meters")
println("--- z: [$(minimum(z_array)),$(maximum(z_array))] meters")
println("")
if saveTransitTimeMaps
    println("Poloidal and toroidal transit times will be saved as .jld2-file keys 'polTransTimes' and 'torTransTimes'.")
    println("")
end
println("Extra keyword arguments specified: ")
println(extra_kw_args)
println("")
println("If you would like to change any settings, please edit the start_calcEpRzTopoMap_template.jl file.")
println("Written by Henrik Järleblad. Last maintained 2022-10-05.")
println("------------------------------------------------------------------------------------------------")
println("")

## ------
# Calculating the topological map
npoints = length(p_array)*length(R_array)*length(z_array) # Number of points for one energy level
pRz_array = Array{Tuple{Float64,Float64,Float64}}(undef,npoints) # To be shared among available processor cores
point = 1
for ip=1:length(p_array)
    for iR=1:length(R_array)
        for iz=1:length(z_array)
            pRz_array[point] = (p_array[ip],R_array[iR],z_array[iz]) # Create a tuple
            global point += 1
        end
    end
end
topoMap_tottot = zeros(length(E_array),length(p_array),length(R_array),length(z_array)) # The total 4D topological map
polTimeMap_tottot = zeros(length(E_array),length(p_array),length(R_array),length(z_array)) # The total 4D poloidal transit time map
torTimeMap_tottot = zeros(length(E_array),length(p_array),length(R_array),length(z_array)) # The total 4D toroidal transit time map
count = 1 # For single-threaded computational progress visualization
for iE=1:length(E_array) ###########################################
E = E_array[iE]
verbose && println("Creating topological map for: $(round(E,digits=2)) keV... ")
if distributed
    if visualizeProgress # if you want the progress to be visualized...
        prg = Progress(npoints) # Define the progress bar
        channel = RemoteChannel(()->Channel{Bool}(npoints),1) # Define the channel from which the progress bar draws data
        dataMap_tot = fetch(@sync begin # Start the distributed computational process, fetch result when done
            @async while take!(channel) # An asynchronous process, with no need for sync, since it simply displays the progress bar
                ProgressMeter.next!(prg)
            end
            @async begin # No internal syncronization needed here either, only external sync needed
                dataMap = @distributed (+) for i=1:length(pRz_array) # Compute one result, and reduce (add) it to a resulting matrix
                    p = pRz_array[i][1]
                    R = pRz_array[i][2]
                    z = pRz_array[i][3]

                    topoMap_i = zeros(length(p_array),length(R_array),length(z_array))
                    polTimeMap_i = zeros(length(p_array),length(R_array),length(z_array))
                    torTimeMap_i = zeros(length(p_array),length(R_array),length(z_array))

                    my_gcp = getGCP(FI_species)

                    if false # for de-bugging purposes. Should it be needed. If so, change to 'true'.
                        o = get_orbit(M,my_gcp(E,p,R,z);wall=wall,store_path=false, verbose=verbose, extra_kw_args...)
                    else
                        o = get_orbit(M,my_gcp(E,p,R,z);wall=wall,store_path=false, extra_kw_args...)
                    end
                    ip = first(findall(x-> x==p,p_array))
                    iR = first(findall(x-> x==R,R_array))
                    iz = first(findall(x-> x==z,z_array))
                    polTimeMap_i[ip,iR,iz] = o.tau_p
                    torTimeMap_i[ip,iR,iz] = o.tau_t

                    if (o.class == :lost) && distinguishLost
                        topoMap_i[ip,iR,iz] = 7
                        polTimeMap_i[ip,iR,iz] = 0.0
                        torTimeMap_i[ip,iR,iz] = 0.0
                    elseif (o.class == :incomplete) && distinguishIncomplete
                        topoMap_i[ip,iR,iz] = 6
                        polTimeMap_i[ip,iR,iz] = 0.0
                        torTimeMap_i[ip,iR,iz] = 0.0
                    elseif o.class == :trapped
                        topoMap_i[ip,iR,iz] = 2
                    elseif o.class == :co_passing
                        topoMap_i[ip,iR,iz] = 3
                    elseif (o.class == :stagnation && o.coordinate.r>=M.axis[1]) # Regular stagnation orbit
                        topoMap_i[ip,iR,iz] = 1
                    elseif (o.class == :stagnation && o.coordinate.r<M.axis[1]) # Counterstagnation orbit
                        topoMap_i[ip,iR,iz] = 8
                    elseif o.class == :potato
                        topoMap_i[ip,iR,iz] = 5
                    elseif o.class == :ctr_passing
                        topoMap_i[ip,iR,iz] = 4
                    else
                        topoMap_i[ip,iR,iz] = 9
                        polTimeMap_i[ip,iR,iz] = 0.0
                        torTimeMap_i[ip,iR,iz] = 0.0
                    end
                    dataMap_i = zeros(3,length(p_array),length(R_array),length(z_array))
                    dataMap_i[1,:,:,:] = topoMap_i
                    dataMap_i[2,:,:,:] = polTimeMap_i
                    dataMap_i[3,:,:,:] = torTimeMap_i
                    put!(channel,true) # Update the progress bar
                    dataMap_i # Declare dataMap_i as result to add to topoMap
                end
                put!(channel,false) # Update progress bar
                dataMap # Delcare dataMap as done/result, so it can be fetched
            end
        end)
    else
        dataMap_tot = @distributed (+) for i=1:length(pRz_array)
            p = pRz_array[i][1]
            R = pRz_array[i][2]
            z = pRz_array[i][3]

            topoMap_i = zeros(length(p_array),length(R_array),length(z_array))
            polTimeMap_i = zeros(length(p_array),length(R_array),length(z_array))
            torTimeMap_i = zeros(length(p_array),length(R_array),length(z_array))

            my_gcp = getGCP(FI_species)

            if false # for de-bugging purposes. Should it be needed. If so, change to 'true'.
                o = get_orbit(M,my_gcp(E,p,R,z);wall=wall,store_path=false, verbose=verbose, extra_kw_args...)
            else
                o = get_orbit(M,my_gcp(E,p,R,z);wall=wall,store_path=false, extra_kw_args...)
            end
            ip = first(findall(x-> x==p,p_array))
            iR = first(findall(x-> x==R,R_array))
            iz = first(findall(x-> x==z,z_array))
            polTimeMap_i[ip,iR,iz] = o.tau_p
            torTimeMap_i[ip,iR,iz] = o.tau_t

            if (o.class == :lost) && distinguishLost
                topoMap_i[ip,iR,iz] = 7
                polTimeMap_i[ip,iR,iz] = 0.0
                torTimeMap_i[ip,iR,iz] = 0.0
            elseif (o.class == :incomplete) && distinguishIncomplete
                topoMap_i[ip,iR,iz] = 6
                polTimeMap_i[ip,iR,iz] = 0.0
                torTimeMap_i[ip,iR,iz] = 0.0
            elseif o.class == :trapped
                topoMap_i[ip,iR,iz] = 2
            elseif o.class == :co_passing
                topoMap_i[ip,iR,iz] = 3
            elseif (o.class == :stagnation && o.coordinate.r>=M.axis[1]) # Regular stagnation orbit
                topoMap_i[ip,iR,iz] = 1
            elseif (o.class == :stagnation && o.coordinate.r<M.axis[1]) # Counterstagnation orbit
                topoMap_i[ip,iR,iz] = 8
            elseif o.class == :potato
                topoMap_i[ip,iR,iz] = 5
            elseif o.class == :ctr_passing
                topoMap_i[ip,iR,iz] = 4
            else
                topoMap_i[ip,iR,iz] = 9
                polTimeMap_i[ip,iR,iz] = 0.0
                torTimeMap_i[ip,iR,iz] = 0.0
            end
            
            dataMap_i = zeros(3,length(p_array),length(R_array),length(z_array))
            dataMap_i[1,:,:,:] = topoMap_i
            dataMap_i[2,:,:,:] = polTimeMap_i
            dataMap_i[3,:,:,:] = torTimeMap_i
            dataMap_i # Declare dataMap_i as result to add to dataMap_tot
        end
    end
else # ... good luck
    topoMap_i = zeros(length(p_array),length(R_array),length(z_array))
    polTimeMap_i = zeros(length(p_array),length(R_array),length(z_array))
    torTimeMap_i = zeros(length(p_array),length(R_array),length(z_array))
    for pRz in pRz_array # Compute one result, and reduce (add) it to a resulting matrix
        p = pRz[1]
        R = pRz[2]
        z = pRz[3]

        verbose && println("$(count)/$(length(topoMap_tottot)):    (E,p,R,z) = ($(E),$(p),$(R),$(z))")

        my_gcp = getGCP(FI_species)

        if false # for de-bugging purposes. Should it be needed. If so, change to 'true'.
            o = get_orbit(M,my_gcp(E,p,R,z);wall=wall,store_path=false, verbose=verbose, extra_kw_args...)
        else
            o = get_orbit(M,my_gcp(E,p,R,z);wall=wall,store_path=false, extra_kw_args...)
        end
        ip = first(findall(x-> x==p,p_array))
        iR = first(findall(x-> x==R,R_array))
        iz = first(findall(x-> x==z,z_array))
        polTimeMap_i[ip,iR,iz] = o.tau_p
        torTimeMap_i[ip,iR,iz] = o.tau_t

        if (o.class == :lost) && distinguishLost
            topoMap_i[ip,iR,iz] = 7
            polTimeMap_i[ip,iR,iz] = 0.0
            torTimeMap_i[ip,iR,iz] = 0.0
        elseif (o.class == :incomplete) && distinguishIncomplete
            topoMap_i[ip,iR,iz] = 6
            polTimeMap_i[ip,iR,iz] = 0.0
            torTimeMap_i[ip,iR,iz] = 0.0
        elseif o.class == :trapped
            topoMap_i[ip,iR,iz] = 2
        elseif o.class == :co_passing
            topoMap_i[ip,iR,iz] = 3
        elseif (o.class == :stagnation && o.coordinate.r>=M.axis[1]) # Regular stagnation orbit
            topoMap_i[ip,iR,iz] = 1
        elseif (o.class == :stagnation && o.coordinate.r<M.axis[1]) # Counterstagnation orbit
            topoMap_i[ip,iR,iz] = 8
        elseif o.class == :potato
            topoMap_i[ip,iR,iz] = 5
        elseif o.class == :ctr_passing
            topoMap_i[ip,iR,iz] = 4
        else
            topoMap_i[ip,iR,iz] = 9
            polTimeMap_i[ip,iR,iz] = 0.0
            torTimeMap_i[ip,iR,iz] = 0.0
        end

        global count = count + 1
    end
    dataMap_tot = zeros(3,length(p_array),length(R_array),length(z_array))
    dataMap_tot[1,:,:,:] = topoMap_i
    dataMap_tot[2,:,:,:] = polTimeMap_i
    dataMap_tot[3,:,:,:] = torTimeMap_i
end

global topoMap_tottot[iE,:,:,:] = dataMap_tot[1,:,:,:] # One energy slice of topological map
global polTimeMap_tottot[iE,:,:,:] = dataMap_tot[2,:,:,:] # One energy slice of poloidal transit times
global torTimeMap_tottot[iE,:,:,:] = dataMap_tot[3,:,:,:] # One energy slice of toroidal transit times
end ###########################################

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

filepath_tm = folderpath_o*"EpRzTopoMap_"*tokamak*"_"*TRANSP_id*"_at"*timepoint*"s_"*FI_species*"_$(length(E_array))x$(length(p_array))x$(length(R_array))x$(length(z_array))"*ident*".jld2"
myfile = jldopen(filepath_tm,true,true,false,IOStream)
write(myfile,"topoMap",topoMap_tottot)
write(myfile,"E_array",collect(E_array))
write(myfile,"p_array",collect(p_array))
write(myfile,"R_array",collect(R_array))
write(myfile,"z_array",collect(z_array))
write(myfile,"distinguishIncomplete",distinguishIncomplete)
write(myfile,"distinguishLost",distinguishLost)
write(myfile,"FI_species",FI_species)
if saveTransitTimeMaps
    write(myfile,"polTransTimes",polTimeMap_tottot)
    write(myfile,"torTransTimes",torTimeMap_tottot)
end
if useDistrFile
    write(myfile,"filepath_distr",filepath_distr)
end
close(myfile)

println("~~~~~~ calcEpRzTopoMap.jl done! ~~~~~~")
