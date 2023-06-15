################################ EpRzWebApp.jl #########################################

#### Description:
# This script provides an application to visualize (guiding-centre) orbits in a tokamak 
# in an interactive and intuitive manner. The topological regions of particle space (E,p,R,z)
# are visualized as well (in terms of orbits). Also, it can visualize a fast-ion distribution 
# and superimpose the topological boundaries to examine fast-ion distribution in terms of orbits.
# In addition, if computed, it can visualize maps of poloidal and toroidal transit times as well.
# Prior to running this script, please make sure you have run the calcEpRzTopoMap.jl script (or equivalent). 
# Also, please make sure that you have noted the path to the saved topological map.

#### Inputs (units given when defined in the script):
# folderpath_OWCF - The path to the OWCF folder on your computer. Needed for correct loading - String
# filepath_distr - (Optional) The path to a .h5/.jld2-file containing the (E,p,R,z) fast-ion distribution - String
# filepath_equil - The path to the .eqdsk-file (or .jld2-file) with the tokamak magnetic equilibrium and geometry - String
# filepath_tm - The path to the .jld2-file containing the topological map (and more) - String
# filepath_tb - (Optional) The path to a .jld2-file containing topological boundaries in (E,pm,Rm) - String
# FI_species - The fast-ion species, e.g. "D", "T", "alpha", "3he" etc - String
# verbose - If set to true, the app will talk a lot! - Bool
# port - An integer. It will correspond to the internet port of your computer through which the app is run - Int64

#### Outputs
# -

#### Saved files
# -

# Script written by Henrik Järleblad and Andrea Valentini. Last maintained 2022-08-25.
#########################################################################################

## --------------------------------------------------------------------------
# Please specify the OWCF folder and let the app change directory to the 
# OWCF folder when orbitsWebApp.jl is executed. This is to be able to load the
# correct versions of the Julia packages as specified in the Project.toml and 
# Manifest.toml files.
folderpath_OWCF = "" # Finish with '/'
cd(folderpath_OWCF)
using Pkg
Pkg.activate(".")

## --------------------------------------------------------------------------
# Must load JLD2 package first, to be able to check filepath_tm and filepath_W for 'extra_kw_args'
using JLD2

## ------
# Inputs
filepath_distr = ""
filepath_equil = "" 
filepath_tm = "" 
filepath_tb = "" # Energy (E) grid points in (E,pm,Rm) must match E grid points of topoMap in (E,p,R,z)
isfile(filepath_tb) && (nR_for_mu=500) # If an (E,pm,Rm) topoBounds file has been specified, specify how many R grid points should be used to resolve the magnetic moment across the plasma
FI_species = "" # Example deuterium: "D"
verbose = true
port = 5432

# EXTRA KEYWORD ARGUMENTS BELOW (these will go into the get_orbit() function from GuidingCenterOrbits.jl)
myfile = jldopen(filepath_tm,false,false,false,IOStream)
if haskey(myfile,"extra_kw_args")
    extra_kw_args = myfile["extra_kw_args"]
else
    extra_kw_args = Dict(:limit_phi => true, :max_tries => 0)
    # limits the number of toroidal turns for orbits
    # The orbit integration algorithm will try progressively smaller timesteps these number of times
end
close(myfile)

## ------
# Loading packages
verbose && println("Loading packages... ")
using Interact
using EFIT
using Equilibrium
using GuidingCenterOrbits
using Plots
using HDF5
using FileIO
using Mux
using WebIO
include(folderpath_OWCF*"extra/dependencies.jl")
include(folderpath_OWCF*"misc/species_func.jl")

verbose && println("Defining helper function (1)... ")
"""
    orbitClass2Int(o)
    orbitClass2Int(o; distinguishIncomplete=true, distinguishLost=true)

Please see apps/calcTopoMap.jl for info.
"""
function orbitClass2Int(M::AbstractEquilibrium, o::Orbit; distinguishIncomplete=false, distinguishLost=false)
    if (o.class == :lost) && distinguishLost
        return 7
    elseif (o.class == :incomplete) && distinguishIncomplete
        return 6
    elseif o.class == :trapped
        return 2
    elseif o.class == :co_passing
        return 3
    elseif (o.class == :stagnation && o.coordinate.r>=magnetic_axis(M)[1]) # Regular stagnation orbit
        return 1
    elseif (o.class == :stagnation && o.coordinate.r<magnetic_axis(M)[1]) # Counterstagnation orbit
        return 8
    elseif o.class == :potato
        return 5
    elseif o.class == :ctr_passing
        return 4
    else
        return 9
    end
end

verbose && println("Defining helper function (2)... ")
function extractTopoBounds(topoMap::Array{Float64,2}) # The (E,p) part of the full 4D topoMap

    topoBounds = zeros(size(topoMap,1),size(topoMap,2))

    for Ei=1:size(topoMap,1)
        for type_i=[1.0,2.0,3.0,4.0,5.0,8.0] # 1.0=stagnation, 2.0=trapped, 3.0=co-passing, 4.0=ctr-passing, 5.0=potato, 8.0=counter-stagnation
            # first and last element index
            first_ind = findfirst(x-> x==type_i,topoMap[Ei,:])
            last_ind = findlast(x-> x==type_i,topoMap[Ei,:])
            if !(typeof(first_ind) == Nothing)
                topoBounds[Ei,first_ind] = 1.0
                topoBounds[Ei,last_ind] = 1.0
            end
    
            # check whether to mark vertical boundaries
            if !(Ei==1 || Ei==size(topoMap,1))
                all_ind_prev = findall(x-> x==type_i,topoMap[Ei-1,:])
                all_ind_next = findall(x-> x==type_i,topoMap[Ei+1,:])
                all_ind = findall(x-> x==type_i,topoMap[Ei,:])
                if((length(all_ind_prev)==0 || length(all_ind_next)==0) && length(all_ind)!=0)
                    topoBounds[Ei,all_ind] .= 1.0
                end
            end
        end
    end

    for pi=1:size(topoMap,2)
        for type_i=[1.0,2.0,3.0,4.0,5.0,8.0] # 1.0=stagnation, 2.0=trapped, 3.0=co-passing, 4.0=ctr-passing, 5.0=potato, 8.0=counter-stagnation
            # first and last element index
            first_ind = findfirst(x-> x==type_i,topoMap[:,pi])
            last_ind = findlast(x-> x==type_i,topoMap[:,pi])
            if !(typeof(first_ind) == Nothing)
                topoBounds[first_ind,pi] = 1.0
                topoBounds[last_ind,pi] = 1.0
            end
    
            # check whether to mark vertical boundaries
            if !(pi==1 || pi==size(topoMap,2))
                all_ind_prev = findall(x-> x==type_i,topoMap[:,pi-1])
                all_ind_next = findall(x-> x==type_i,topoMap[:,pi+1])
                all_ind = findall(x-> x==type_i,topoMap[:,pi])
                if((length(all_ind_prev)==0 || length(all_ind_next)==0) && length(all_ind)!=0)
                    topoBounds[all_ind,pi] .= 1.0
                end
            end
        end
    end
    return topoBounds
end

## ------
verbose && println("Loading tokamak equilibrium... ")
if ((split(filepath_equil,"."))[end] == "eqdsk") || ((split(filepath_equil,"."))[end] == "geqdsk") 
    M, wall = read_geqdsk(filepath_equil,clockwise_phi=false) # Assume counter-clockwise phi-direction
    jdotb = M.sigma # The sign of the dot product between the plasma current and the magnetic field
else # Otherwise, assume magnetic equilibrium is a saved .jld2 file
    myfile = jldopen(filepath_equil,false,false,false,IOStream)
    M = myfile["S"]
    wall = myfile["wall"]
    close(myfile)
    jdotb = (M.sigma_B0)*(M.sigma_Ip)
end

#########################################################################################
verbose && println("Computing flux function on 100x100 (R,z)-grid (to plot flux surfaces)... ")
flux_r = range(extrema(wall.r)...,length=100)
flux_z = range(extrema(wall.z)...,length=100)
inds = CartesianIndices((length(flux_r),length(flux_z)))
psi_rz = [M(flux_r[ind[1]], flux_z[ind[2]]) for ind in inds]
psi_mag, psi_bdry = psi_limits(M)

## ------
# Read the .jld2-file for displaying the topology of orbit space
verbose && println("Loading topological map... ")
myfile = jldopen(filepath_tm,false,false,false,IOStream)
topoMap = myfile["topoMap"]
distinguishIncomplete = myfile["distinguishIncomplete"]
distinguishLost = myfile["distinguishLost"]
E_array = myfile["E_array"]
p_array = myfile["p_array"]
R_array = myfile["R_array"]
z_array = myfile["z_array"]
poltor = false
if haskey(myfile,"polTransTimes")
    verbose && println("-- Maps of poloidal and toroidal transit times found! Including... ")
    polTransTimes = myfile["polTransTimes"]
    torTransTimes = myfile["torTransTimes"]
    poltor = true
end
jac = false
if haskey(myfile,"jacobian")
    verbose && println("-- Jacobian from (x,y,z,vx,vy,vz) to (E,p,R,z) found! Including... ")
    jacobian = myfile["jacobian"]
    jac = true
end
close(myfile)

## ------
# If provided, read the .jld2-file for displaying the topological boundaries in (E,pm,Rm)
topobo = false
if isfile(filepath_tb)
    verbose && println("Loading (E,pm,Rm) topological boundaries... ")
    myfile = jldopen(filepath_tb,false,false,false,IOStream)
    topoBounds = myfile["topoBounds"]
    E_array_tb = myfile["E_array"]
    pm_array = myfile["pm_array"]
    Rm_array = myfile["Rm_array"]
    close(myfile)

    if !(E_array==E_array_tb)
        error("Energy (E) grid points of (E,pm,Rm) topological boundaries did not match E grid points of topological map in (E,p,R,z). Please correct and re-try!")
    end

    topobo = true

    ### WRITE CODE TO PRE-COMPUTE R_grid, B_array, ... etc FOR mu_func(). 
    ### USE INPUT ARGUMENT nR_for_mu TO DETERMINE LENGTH OF R_grid, B_array, ... etc
end

## ------
# Check the topological map for incomplete orbits
# If any are found, attempt to re-compute them
verbose && print("Checking topological map for incomplete orbits... ")
incompleteInds = findall(x-> x==6,topoMap)
if length(incompleteInds)>0
    verbose && println("Found $(length(incompleteInds))!")
    verbose && print("Attempting to re-compute the incomplete orbits... ")

    numOsavedOrbs = 0
    for incompleteInd in incompleteInds
        iE = incompleteInd[1]
        ip = incompleteInd[2]
        iR = incompleteInd[3]
        iz = incompleteInd[4]

        my_gcp = getGCP(FI_species)
        o = get_orbit(M,my_gcp(E_array[iE],p_array[ip],R_array[iR],z_array[iz]); wall=wall, extra_kw_args...)

        orbInt = orbitClass2Int(M, o; distinguishIncomplete=distinguishIncomplete, distinguishLost=distinguishLost)
        if !(orbInt==6)
            global numOsavedOrbs += 1
        end
        topoMap[iE,ip,iR,iz] = orbInt
    end
    verbose && println("Done! (successfully re-computed: $(numOsavedOrbs))")
else
    verbose && println("Found none!")
end

if isfile(filepath_distr)
    verbose && println("Loading fast-ion distribution... ")
    # Determine fast-ion distribution file type
    fileext_distr = (split(filepath_distr,"."))[end] # Assume last part after final '.' is the file extension
    if fileext_distr=="h5" || fileext_distr=="hdf5"
        myfile = h5open(filepath_distr,"r")
        F = read(myfile["f"])
        E_array_FI = read(myfile["energy"])
        p_array_FI = read(myfile["pitch"])
        R_array_FI = read(myfile["R"])
        if haskey(myfile,"Z")
            z_array_FI = read(myfile["Z"])
        else
            z_array_FI = read(myfile["z"])
        end
        close(myfile)

        if !(size(F,1)==length(E_array_FI) && size(F,2)==length(p_array_FI) && size(F,3)==length(R_array_FI) && size(F,4)==length(z_array_FI))
            println("Correctly rotating loaded .h5 fast-ion distribution (dimensions were in wrong order)... ")
            f = zeros(length(E_array_FI),length(p_array_FI),length(R_array_FI),length(z_array_FI))
            for i=1:size(f,1)
                for j=1:size(f,2)
                    for k=1:size(f,3)
                        for l=1:size(f,4)
                            f[i,j,k,l] = F[l,k,j,i]
                        end
                    end
                end
            end
            F = f
        end
    elseif fileext_distr=="jld2"
        myfile = jldopen(filepath_distr,false,false,false,IOStream)
        if haskey(myfile,"F_ps")
            F = myfile["F_ps"]
        else
            F = myfile["f"] # Otherwise, assume has 'f' key
        end
        E_array_FI = myfile["energy"]
        p_array_FI = myfile["pitch"]
        R_array_FI = myfile["R"]
        if haskey(myfile,"Z")
            z_array_FI = myfile["Z"]
        else
            z_array_FI = myfile["z"]
        end
        close(myfile)
    else
        error("Unknown fast-ion distribution file extensions. Please correct and instead use a .h5, .hdf5 or .jld2 file.")
    end

    if maximum(R_array_FI)>100.0 # Convert from centimeter (cm) to meter (m)
        verbose && println("Re-scaling fast-ion data from centimeters to meters... ")
        R_array_FI = R_array_FI ./100
        z_array_FI = z_array_FI ./100 # Assume both R and z need conversion from centimeter to meter
        F = F .* reshape(100 .*ones(length(R_array_FI)),(1,1,length(R_array_FI),1)) # To account for Jacobian (R) cm->m scaling
        F = (100*100) .* F # To account for dR and dz cm->m scaling
    end
    if !(E_array==E_array_FI && p_array==p_array_FI && R_array==R_array_FI && z_array==z_array_FI)
        # We need interpolation
        verbose && println("Grid points of topological map and loaded fast-ion data did not match. Interpolating... ")
        F = interpFps(F,E_array_FI,p_array_FI,R_array_FI,z_array_FI,E_array,p_array,R_array,z_array)
    end

    dE4D, dp4D, dR4D, dz4D = get4DDiffs(E_array, p_array,R_array,z_array) 
    F_Rz = dropdims(sum(dE4D .* dp4D .* F,dims=(1,2)),dims=(1,2)) # Will need for (R,z) heatmap plot
end

verbose && println("Building web application... ")
R_hfs = minimum(wall.r) # R-coord of high-field side wall
R_lfs = maximum(wall.r) # R-coord of low-field side wall
phi = collect(0:1:359).*(2*pi/180.0) # Toroidal angle
topview_R_hfs_x = (R_hfs).*cos.(phi)
topview_R_hfs_y = (R_hfs).*sin.(phi)
topview_R_lfs_x = (R_lfs).*cos.(phi)
topview_R_lfs_y = (R_lfs).*sin.(phi)

E_array = vec(collect(E_array)) # Ensure type Array{Float64,1}
p_array = vec(collect(p_array)) # Ensure type Array{Float64,1}
R_array = vec(collect(R_array)) # Ensure type Array{Float64,1}
z_array = vec(collect(z_array)) # Ensure type Array{Float64,1}

## ------
# The web application
verbose && println("--- You can access the EpRzWebApp via an internet web browser when you see 'Task (runnable)...' ")
verbose && println("--- When 'Task (runnable)...' has appeared, please visit the website localhost:$(port) ---")
verbose && println("--- Remember: It might take a minute or two to load the webpage. Please be patient. ---")
function app(req) # Here is where the app starts!
    @manipulate for tokamak_wall = Dict("on" => true, "off" => false), E=E_array, p=p_array, R=R_array, z=z_array, study = Dict("Valid orbits" => :orbs, "Fast-ion distribution" => :FI, "Poloidal times" => :tpol, "Toroidal times" => :ttor, "(E,pm,Rm) at (R,z)" => :OS), save_plots = Dict("on" => true, "off" => false)
        
        my_gcp = getGCP(FI_species)

        o = get_orbit(M,my_gcp(E,p,R,z); wall=wall, extra_kw_args...)
        
        topview_o_x = cos.(o.path.phi).*(o.path.r)
        topview_o_y = sin.(o.path.phi).*(o.path.r)

        orb_color = :black
        orb_linestyle = :solid
        if o.class==:invalid
            orb_color = :gray
            orb_linestyle = :dash
        elseif o.class == :lost
            orb_color = :brown
        elseif o.class == :incomplete # If this happens, you are in trouble. Because it will likely take forever to calculate. Please just re-start the app instead.
            orb_color = :yellow
        elseif o.class == :trapped
            orb_color = :blue
        elseif o.class == :co_passing
            orb_color = :green
        elseif (o.class == :stagnation && o.coordinate.r>=magnetic_axis(M)[1]) # Regular stagnation orbit
            orb_color = :red
        elseif o.class == :potato
            orb_color = :orange
        elseif o.class == :ctr_passing
            orb_color = :purple
        elseif (o.class == :stagnation && o.coordinate.r<magnetic_axis(M)[1]) # Counter-stagnation orbit
            orb_color = :pink
        else
            error("Something's gone wrong!!! Orbit class unknown!")
        end

        #cross-sectional (R,z) plot
        if study==:FI && !(filepath_distr=="")
            plt_crs = Plots.heatmap(R_array, z_array, F_Rz', fillcolor=cgrad([:white, :darkblue, :green, :yellow, :orange, :red]))
            plt_crs = Plots.scatter!([magnetic_axis(M)[1]],[magnetic_axis(M)[2]],label="Magnetic axis", mc=:gray, aspect_ratio=:equal, xlabel="R [m]", ylabel="z [m]",title="Poloidal cross-section")
        else
            plt_crs = Plots.scatter([magnetic_axis(M)[1]],[magnetic_axis(M)[2]],label="Magnetic axis", mc=:gray, aspect_ratio=:equal, xlabel="R [m]", ylabel="z [m]",title="Poloidal cross-section")
        end
        if tokamak_wall
            plt_crs = Plots.contour!(flux_r,flux_z,psi_rz',levels=collect(range(psi_mag,stop=psi_bdry,length=5)),color=:gray, linestyle=:dot,linewidth=1.5, label="",colorbar=false)
            wall_dR = maximum(wall.r)-minimum(wall.r)
            plt_crs = Plots.plot!(wall.r,wall.z, label="Tokamak wall", color=:black, linewidth=1.5,xaxis=[minimum(wall.r)-wall_dR/10,maximum(wall.r)+wall_dR])
        end
        plt_crs = Plots.plot!(o.path.r,o.path.z, label="$(o.class) orbit", color=orb_color, linestyle=orb_linestyle, linewidth=1.5)
        plt_crs = Plots.vline!([R],linestyle=:dot,color=:gray,linewidth=1.0,label="")
        plt_crs = Plots.hline!([z],linestyle=:dot,color=:gray,linewidth=1.0,label="")
        plt_crs = Plots.scatter!([R],[z], mc=orb_color, label="(R,z)", grid=false)
        if save_plots
            plt_crs = Plots.plot!(dpi=600)
            png(plt_crs, "plt_crs_$(round(E, digits=2))_$(round(p, digits=2))_$(round(R,digits=2))_$(round(z,digits=2))")
        end

        #topological (E,p) plot
        Rci = argmin(abs.(collect(R_array) .- R))
        zci = argmin(abs.(collect(z_array) .- z))
        Eci = argmin(abs.(collect(E_array) .- E))
        pci = argmin(abs.(collect(p_array) .- p))

        if study==:FI && !(filepath_distr=="")
            topoBounds = extractTopoBounds(topoMap[:,:,Rci,zci])
            ones_carinds = findall(x-> x==1.0,topoBounds)
            E_scatvals_tb = zeros(length(ones_carinds))
            p_scatvals_tb = zeros(length(ones_carinds))
            for (ind,carinds) in enumerate(ones_carinds)
                E_scatvals_tb[ind] = E_array[carinds[1]]
                p_scatvals_tb[ind] = p_array[carinds[2]]
            end
            plt_topo = Plots.heatmap(E_array,p_array, F[:,:,Rci,zci]', fillcolor=cgrad([:white, :darkblue, :green, :yellow, :orange, :red]), xlabel="E [keV]", ylabel="p [-]", title="f(E,p) at (R,z)", right_margin=3Plots.mm)
            ms = save_plots ? 2.6 : 1.8
            #plt_weights = Plots.scatter!(E_scatvals_tb,p_scatvals_tb,markersize=ms,leg=false,markercolor=:black)
        elseif study==:tpol && poltor
            topoBounds = extractTopoBounds(topoMap[:,:,Rci,zci])
            ones_carinds = findall(x-> x==1.0,topoBounds)
            E_scatvals_tb = zeros(length(ones_carinds))
            p_scatvals_tb = zeros(length(ones_carinds))
            for (ind,carinds) in enumerate(ones_carinds)
                E_scatvals_tb[ind] = E_array[carinds[1]]
                p_scatvals_tb[ind] = p_array[carinds[2]]
            end
            pTT_microsecs = polTransTimes[:,:,Rci,zci] ./(1.0e-6) # Convert from seconds to microseconds
            nz_coords = findall(x-> x>0.0,pTT_microsecs) # Find the 2D matrix coordinates of all non-zero elements
            my_coords = length(nz_coords) > 1 ? nz_coords : CartesianIndices(size(pTT_microsecs)) # Are there actually more than one non-zero element? If not, use all elements
            min_pol, max_pol = extrema(pTT_microsecs[my_coords]) # Find minimum and maximum values
            min_OOM, max_OOM = (floor(log10(min_pol)),ceil(log10(max_pol))) # The orders of magnitude of the minimum and maximum values
            if !((max_OOM-min_OOM)==0.0) && (length(nz_coords) > 1) # All values NOT within same order of magnitude AND more than one non-zero element. Use logarithmic colorbar
                plt_topo = Plots.heatmap(E_array,p_array, pTT_microsecs', title="tau_pol(E,p) at (R,z) [microseconds] \n tau_pol($(round(E,sigdigits=3)),$(round(p,sigdigits=2)))=$(round(o.tau_p /(1.0e-6),sigdigits=3)) microseconds", fillcolor=cgrad([:white, :darkblue, :green, :yellow, :orange, :red]), xlabel="E [keV]", ylabel="p [-]", colorbar=true, colorbar_scale=:log10, clims = (10^min_OOM, 10^max_OOM), top_margin=3Plots.mm) # Get nice powers-of-ten limits for the colorbar
            else # Else, use linear colorbar 
                plt_topo = Plots.heatmap(E_array,p_array, pTT_microsecs', title="tau_pol(E,p) at (R,z) [microseconds] \n tau_pol($(round(E,sigdigits=3)),$(round(p,sigdigits=2)))=$(round(o.tau_p /(1.0e-6),sigdigits=3)) microseconds", fillcolor=cgrad([:white, :darkblue, :green, :yellow, :orange, :red]), xlabel="E [keV]", ylabel="p [-]", colorbar=true, top_margin=3Plots.mm)
            end
            ms = save_plots ? 2.6 : 1.8
            #plt_weights = Plots.scatter!(E_scatvals_tb,p_scatvals_tb,markersize=ms,label="",markercolor=:black,leg=false)
        elseif study==:ttor && poltor
            topoBounds = extractTopoBounds(topoMap[:,:,Rci,zci])
            ones_carinds = findall(x-> x==1.0,topoBounds)
            E_scatvals_tb = zeros(length(ones_carinds))
            p_scatvals_tb = zeros(length(ones_carinds))
            for (ind,carinds) in enumerate(ones_carinds)
                E_scatvals_tb[ind] = E_array[carinds[1]]
                p_scatvals_tb[ind] = p_array[carinds[2]]
            end
            tTT_microsecs = torTransTimes[:,:,Rci,zci] ./(1.0e-6)
            nz_coords = findall(x-> x>0.0,tTT_microsecs) # Find the 2D matrix coordinates of all non-zero elements
            my_coords = length(nz_coords) > 1 ? nz_coords : CartesianIndices(size(tTT_microsecs)) # Are there actually more than one non-zero element? If not, use all elements
            min_pol, max_pol = extrema(tTT_microsecs[my_coords]) # Find minimum and maximum values
            min_OOM, max_OOM = (floor(log10(min_pol)),ceil(log10(max_pol))) # The orders of magnitude of the minimum and maximum values
            if !((max_OOM-min_OOM)==0.0) && (length(nz_coords) > 1) # All values NOT within same order of magnitude AND more than one non-zero element. Use logarithmic colorbar
                plt_topo = Plots.heatmap(E_array,p_array, tTT_microsecs', title="tau_tor(E,p) at (R,z) [microseconds] \n tau_tor($(round(E,sigdigits=3)),$(round(p,sigdigits=2)))=$(round(o.tau_t /(1.0e-6),sigdigits=3)) microseconds", fillcolor=cgrad([:white, :darkblue, :green, :yellow, :orange, :red]), xlabel="E [keV]", ylabel="p [-]", colorbar=true, colorbar_scale=:log10, clims = (10^min_OOM, 10^max_OOM),top_margin=3Plots.mm) # Get nice powers-of-ten limits for the colorbar
            else # Else, use linear colorbar 
                plt_topo = Plots.heatmap(E_array,p_array, tTT_microsecs', title="tau_tor(E,p) at (R,z) [microseconds] \n tau_tor($(round(E,sigdigits=3)),$(round(p,sigdigits=2)))=$(round(o.tau_t /(1.0e-6),sigdigits=3)) microseconds", fillcolor=cgrad([:white, :darkblue, :green, :yellow, :orange, :red]), xlabel="E [keV]", ylabel="p [-]", colorbar=true,top_margin=3Plots.mm)
            end
            ms = save_plots ? 2.6 : 1.8
            #plt_weights = Plots.scatter!(E_scatvals_tb,p_scatvals_tb,markersize=ms,label="",markercolor=:black,leg=false)
        elseif study==:OS && topobo
            ones_carinds = findall(x-> x==1.0,topoBounds[Eci,:,:])
            pm_scatvals_tb = zeros(length(ones_carinds))
            Rm_scatvals_tb = zeros(length(ones_carinds))
            for (ind,carinds) in enumerate(ones_carinds)
                pm_scatvals_tb[ind] = pm_array[carinds[1]]
                Rm_scatvals_tb[ind] = Rm_array[carinds[2]]
            end

            array_o_hcTopo_tuples = map(i-> (HamiltonianCoordinate(M,my_gcp(E,p_array[i],R,z)), topoMap[Eci,i,Rci,zci]), 1:length(p_array))
            array_o_EPRcTopo_tuples = map(x-> (EPRCoordinate4(; R_grid=R_grid, B_array=B_on_R_grid, ...),x[2]), array_o_hcTopo_tuples)

        else
            topoMap_raw = topoMap[:,:,Rci,zci] # might not include all integers between 1 and 9 (needed for correct Set1_9 heatmap coloring)
            topoMap_ext = ones(size(topoMap_raw,1),size(topoMap_raw,2)+1) # We ensure 1.0 in data by using ones() function. Add one extra column
            topoMap_ext[1:size(topoMap_raw,1),1:size(topoMap_raw,2)] .= topoMap_raw # put in the true topoMap to be plotted
            topoMap_ext[end,end] = 9.0 # Get the 9.0. This is such an ugly solution (to acquire correct orbit type coloring)...
            p_array_ext = vcat(p_array,2*p_array[end]-p_array[end-1]) # Extend p_array by one dummy element
            plt_topo = Plots.heatmap(E_array,p_array_ext,topoMap_ext',color=:Set1_9,legend=false,xlabel="E [keV]", ylabel="p [-]", title="(E,p) orbit topology at (R,z)", ylims=(minimum(p_array),maximum(p_array)))
        end
        plt_topo = Plots.vline!([E],linestyle=:dot,color=:gray,linewidth=1.0)
        plt_topo = Plots.hline!([p],linestyle=:dot,color=:gray,linewidth=1.0)
        plt_topo = Plots.scatter!([E],[p],markershape=:circle,mc=orb_color,legend=false,markersize=6, grid=false)
        if save_plots
            plt_topo = Plots.plot!(dpi=600)
            png(plt_topo, "plt_topo_$(round(E, digits=2))_$(round(p, digits=2))_$(round(R,digits=2))_$(round(z,digits=2))")
        end
        
        # Orbit fractions plot
        # Please see signalWebApp.jl for more info
        o_types=[1.0,2.0,3.0,4.0,5.0,8.0]
        valid_orb_inds = findall(x-> x in o_types,topoMap[:,:,Rci,zci])
        if study==:FI && !(filepath_distr=="")
            F_Ep = F[:,:,Rci,zci]
            sum_o_FI_distr = sum(F_Ep[valid_orb_inds])
            stagnation_fraction = sum(F_Ep[findall(x-> x==1.0,topoMap[:,:,Rci,zci])]) / sum_o_FI_distr
            trapped_fraction = sum(F_Ep[findall(x-> x==2.0,topoMap[:,:,Rci,zci])]) / sum_o_FI_distr
            copassing_fraction = sum(F_Ep[findall(x-> x==3.0,topoMap[:,:,Rci,zci])]) / sum_o_FI_distr
            counterpassing_fraction = sum(F_Ep[findall(x-> x==4.0,topoMap[:,:,Rci,zci])]) / sum_o_FI_distr
            potato_fraction = sum(F_Ep[findall(x-> x==5.0,topoMap[:,:,Rci,zci])]) / sum_o_FI_distr
            counterstagnation_fraction = sum(F_Ep[findall(x-> x==8.0,topoMap[:,:,Rci,zci])]) / sum_o_FI_distr

            plt_orb_fracs = Plots.plot([1],[stagnation_fraction],color=:red,linetype=:bar,label="Stagnation")
            plt_orb_fracs = Plots.plot!([2],[trapped_fraction],color=:blue,linetype=:bar,label="Trapped")
            plt_orb_fracs = Plots.plot!([3],[copassing_fraction],color=:green,linetype=:bar,label="Co-passing")
            plt_orb_fracs = Plots.plot!([4],[counterpassing_fraction],color=:purple,linetype=:bar,label="Counter-passing")
            plt_orb_fracs = Plots.plot!([5],[potato_fraction],color=:orange,linetype=:bar,label="Potato")
            plt_orb_fracs = Plots.plot!([6],[counterstagnation_fraction],color=:pink,linetype=:bar,label="Counter-stagnation", xticks=false, title="Orbit types distribution of f(E,p) at (R,z)",ylabel="Fraction",ylims=(0.0, 1.0),xlabel="Orbit types")
        else
            sum_o_topoMap = sum(topoMap[valid_orb_inds,Rci,zci])
            stagnation_fraction = sum(topoMap[findall(x-> x==1.0,topoMap[:,:,Rci,zci]),Rci,zci]) / sum_o_topoMap
            trapped_fraction = sum(topoMap[findall(x-> x==2.0,topoMap[:,:,Rci,zci]),Rci,zci]) / sum_o_topoMap
            copassing_fraction = sum(topoMap[findall(x-> x==3.0,topoMap[:,:,Rci,zci]),Rci,zci]) / sum_o_topoMap
            counterpassing_fraction = sum(topoMap[findall(x-> x==4.0,topoMap[:,:,Rci,zci]),Rci,zci]) / sum_o_topoMap
            potato_fraction = sum(topoMap[findall(x-> x==5.0,topoMap[:,:,Rci,zci]),Rci,zci]) / sum_o_topoMap
            counterstagnation_fraction = sum(topoMap[findall(x-> x==8.0,topoMap[:,:,Rci,zci]),Rci,zci]) / sum_o_topoMap

            plt_orb_fracs = Plots.plot([1],[stagnation_fraction],color=:red,linetype=:bar,label="Stagnation")
            plt_orb_fracs = Plots.plot!([2],[trapped_fraction],color=:blue,linetype=:bar,label="Trapped")
            plt_orb_fracs = Plots.plot!([3],[copassing_fraction],color=:green,linetype=:bar,label="Co-passing")
            plt_orb_fracs = Plots.plot!([4],[counterpassing_fraction],color=:purple,linetype=:bar,label="Counter-passing")
            plt_orb_fracs = Plots.plot!([5],[potato_fraction],color=:orange,linetype=:bar,label="Potato")
            plt_orb_fracs = Plots.plot!([6],[counterstagnation_fraction],color=:pink,linetype=:bar,label="Counter-stagnation", xticks=false, title="Valid orbit types at (R,z)",ylabel="Fraction",ylims=(0.0, 1.0),xlabel="Orbit types")
        end
        if save_plots
            plt_orb_fracs = Plots.plot!(dpi=600)
            png(plt_orb_fracs, "plt_orb_fracs_$(round(E, digits=2))_$(round(p, digits=2))_$(round(R,digits=2))_$(round(z,digits=2))")
        end

        #pitch visualization plot
        x_array = collect(range(0.0,stop=o.tau_p,length=length(o.path.pitch)))
        my_xlabel = "Poloidal time [microseconds]"
        plt_pitc = Plots.plot(x_array ./(1.0e-6), o.path.pitch,color=orb_color,title="Pitch along orbit path",label="",linewidth=1.5, xlabel=my_xlabel, right_margin=4Plots.mm)
        if save_plots
            plt_pitc = Plots.plot!(dpi=600)
            png(plt_pitc, "plt_pitc_$(round(E, digits=2))_$(round(p, digits=2))_$(round(R,digits=2))_$(round(z,digits=2))")
        end

        vbox(vskip(1em),
            hbox(Plots.plot(plt_crs),Plots.plot(plt_orb_fracs)),
            hbox(Plots.plot(plt_topo), Plots.plot(plt_pitc))
        )
    end
end
webio_serve(page("/",app), port)
