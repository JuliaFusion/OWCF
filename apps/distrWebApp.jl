################################### distrWebApp.jl ##########################################

#### Description:
# This script provides an application to visualize orbit-space distribution functions in an 
# interactive and intuitive manner. It can also be used to visualize other 3D orbit-space 
# quantities, such as transformed flux surfaces to orbit space etc (for orbit weight functions,
# please use weightsWebApp.jl or weightWebApp.jl. For topological maps and maps of the poloidal
# and toroidal transit times, please use orbitsWebApp.jl).
#
# It visualizes distribution functions as a function of orbit space (2D slices of 3D space). 
#
# A second orbit-space distribution (or other quantity) can be specified. In that case, the app
# will show to slices side-by-side so that the two 3D quantities can be compared as the user
# scans through fast-ion energy slices.
#
# Prior to running this script, please make sure you have run the following scripts:
#   - ps2os.jl (or equivalent)
#   - F_os_1Dto3D.jl (or equivalent)
#   - extractTopoBounds.jl
# And that you have noted the paths to the outputs. You will also need the following: 
#   - .eqdsk/.jld2-file containing tokamak magnetic equilibrium and geometry
#
# Run the script by starting the Julia command prompt, and type 'include("app/distrWebApp.jl")'
# when standing in the OWCF folder (please see the howToOWCF.pdf document for more info).
# NOTE! Make sure to set the inputs in the ## ------ # Inputs section first.
# After running the above command, the script will create the web application (might take several
# minutes). When you see "Task runnable...", open a web browser of your choice and type 
# 'localhost:'+ the port number you have set with the port variable.
# The web application will then load (might take a minute or two). Et voilà, have fun!
#
# ALSO NOTE! It is recommended to ensure that the distribution function and the topological boundaries
# have exactly the same dimensions and ranges, by using the same orbit-space grid to calculate the 
# topological map in getTopoMap.jl in the first place.

#### Inputs (units given when defined in the script):
# port - The I/O port on which to host the web application (0001-9999) - Int64
# verbose - If set to true, the app will talk a lot! - Bool
# filepath_tb - The path to the .jld2-file containing the topological boundaries - String
# filepath_equil - The path to the .eqdsk-file with the tokamak magnetic equilibrium and geometry - String
# filepath_distr - The path to the 3D orbit-space fast-ion distribution to be visualized - String
# filepath_distr_2 - The path to a second 3D orbit-space distribution (/quantity). To be compared to the first - String
# filepath_no - The path to the 3D orbit-space fast-ion null-measurement region boundaries. Enabled via showNo - String
#
# THere are also more inputs further down in the script. These inputs can be ignored, unless
# you are specifically visualizing flux surfaces in orbit space. The inputs are
# psi_value - The value between 0.0 and 1.0 corresponding to the flux surface you are visualizing in orbit space.

#### Outputs
# -

#### Saved files
# - 

#### Other
# Warning! Please note! For orbit-space grids containing more than approximately 150 000 valid orbits (e.g. 20x100x100),
# you should NOT use distrWebApp.jl (or any other interactive app). As of OWCF version 1.0, the 
# web interface simply becomes too slow. Please do instead plot the energy slices manually instead (in scratch.jl for example).
# You do this by, for example, coding
#
# folderpath_OWCF = "/the/path/to/the/OWCF/folder/"
# cd(folderpath_OWCF)
# using Pkg
# Pkg.activate(".")
# using JLD2
# using Plots
#
# myfile = jldopen(filepath_distr, false, false, false, IOStream)
# F_os_3D = myfile["F_os_3D"]
# E_array = myfile["E_array"]
# pm_array = myfile["pm_array"]
# Rm_array = myfile["Rm_array"]
# close(myfile)
#
# myfile = jldopen(filepath_tb,false,false,false,IOStream)
# topoBounds = myfile["topoBounds"]
# close(myfile)
#
# E = 150.0 # Example of 150 keV
# iE = argmin(abs.(E_array - E)) # Find the closest value to E in E_array
#
# ones_carinds = findall(x-> x==1.0,topoBounds[iE,:,:])
# pm_scatvals_tb = zeros(length(ones_carinds))
# Rm_scatvals_tb = zeros(length(ones_carinds))
# for (ind,carinds) in enumerate(ones_carinds)
#     pm_scatvals_tb[ind] = pm_array[carinds[1]]
#     Rm_scatvals_tb[ind] = Rm_array[carinds[2]]
# end
#
# Plots.heatmap(Rm_array, pm_array, (F_os_3D[iE,:,:])./maximum(F_os_3D[iE,:,:]), colorbar=true, title="Fast-ion distribution slice  ($(round(maximum(F_os_3D[iE,:,:]), sigdigits=4)) = 1.0)",fillcolor=cgrad([:white, :darkblue, :green, :yellow, :orange, :red]))
# Plots.scatter!(Rm_scatvals_tb,pm_scatvals_tb,markersize=ms,leg=false,markercolor=:black, xlabel="Rm [m]", ylabel="pm")

# Script written by Henrik Järleblad. Last maintained 2025-08-01.
###############################################################################################

## --------------------------------------------------------------------------
# Please specify the OWCF folder and let the app change directory to the 
# OWCF folder when distrWebApp.jl is executed. This is to be able to load the
# correct versions of the Julia packages as specified in the Project.toml and 
# Manifest.toml files.
folderpath_OWCF = "" # Finish with '/'
cd(folderpath_OWCF)
using Pkg
Pkg.activate(".")

## ------
# Loading JLD2 before reading the input variables, to enable input variable selection
using JLD2

## ------
# Inputs
port = 3333
verbose = true
FI_species = "" # Example deuterium: "D"
filepath_tb = ""
filepath_equil = ""
filepath_distr = ""
filepath_distr_2 = ""
(showNo = false) && (filepath_no = "")

# If you want to visualize output from e.g. a reconstruction, it might be in 4D format.
# The first dimension then denotes different reconstructions for different regularization parameter values
# Use 'id' to select which reconstruction to visualize, if filepath_distr is a reconstruction output
# Use 'id_2' to select which reconstruction to visualize, if filepath_distr_2 is a reconstruction output
id = 1
id_2 = 1

# EXTRA KEYWORD ARGUMENTS BELOW (these will go into the orbit_grid() and get_orbit() functions from GuidingCenterOrbits.jl)
myfile = jldopen(filepath_tb,false,false,false,IOStream)
if haskey(myfile,"extra_kw_args")
    extra_kw_args = myfile["extra_kw_args"]
else
    extra_kw_args = Dict(:limit_phi => true, :max_tries => 0)
    # limits the number of toroidal turns for orbits
    # The orbit integration algorithm will try progressively smaller timesteps these number of times
end
close(myfile)

###############################################################################################
###############################################################################################
###############################################################################################
###############################################################################################
###############################################################################################
###############################################################################################
######################################### START OF APP ######################################## 
######################################### DO NOT CHANGE #######################################
######################################### ANYTHING BELOW ######################################
###############################################################################################
###############################################################################################
###############################################################################################
###############################################################################################
###############################################################################################
###############################################################################################

## ------
# Loading packages
verbose && println("Loading packages... ")
using Interact
using EFIT
using Equilibrium
using GuidingCenterOrbits
using OrbitTomography
using SparseArrays
using Plots
using FileIO
using Mux
using WebIO
include(folderpath_OWCF*"misc/species_func.jl")
include(folderpath_OWCF*"extra/dependencies.jl")

## ------
# Read the .jld2-file for displaying the topological boundaries of orbit space.
# Pre-calculated with getTopoMap.jl followed by extractTopoBounds.jl.
println("Loading topological boundaries... ")
myfile = jldopen(filepath_tb,false,false,false,IOStream)
topoBounds = myfile["topoBounds"]
close(myfile)

## ------
# Read the .jld2-file for displaying the null-region boundaries of orbit space.
# Pre-calculated with extractnullOrbs.jl.
if showNo
    if !isfile(filepath_no)
        error("Null-measurement boundaries set for plot (showNo=true). But filepath_no is invalid. Please correct and re-try.")
    end
    verbose && println("Loading null-region boundaries... ")
    myfile = jldopen(filepath_no,false,false,false,IOStream)
    if haskey(myfile,"nullOrbs_indices")
        verbose && println("Found null-orbit (4D) indices (pressumably computed with extractNullOrbs.jl). Loading... ")
        nullOrbs_indices = myfile["nullOrbs_indices"]
    elseif haskey(myfile,"nullOrbs_1D") && haskey(myfile,"og")
        verbose && println("Found null-orbit vector and orbit grid (pressumably computed with solveInverseProblem.jl). Loading... ")
        nullOrbs_1D = myfile["nullOrbs_1D"]
        og_nb = myfile["og"]
    else
        @warn "Input null-orbits file filepath_no did not have known null-orbit keys. Please re-compute extractNullOrbits.jl with input variable 'include2Dto4D' set to true, and use its output file as new filepath_no. Or re-compute solveInverse.jl with autoDetectAndUseNullOrbs set to true and use its output file as new filepath_no."
    end
    close(myfile)
end

## ------
# Load the 3D orbit-space fast-ion distribution, the slices of which will be visualized
println("Loading 3D orbit-space fast-ion distribution... ")
myfile = jldopen(filepath_distr,false,false,false,IOStream)
mykeys = keys(myfile)
psi3D = false # Please ignore. This will be switched on automatically.
psi_value = 0.6
if "F_os_3D" in mykeys
    F_os_3D = myfile["F_os_3D"]
elseif "psi_os_3D" in mykeys
    F_os_3D = myfile["psi_os_3D"]
    Rz_array = myfile["Rz_array"]
    psi3D = true
elseif "F_reconstructions_3D" in mykeys
    F_os_3D = myfile["F_reconstructions_3D"]
else
    error("Error! Unknown orbit-space 3D object in 'filepath_distr'. Please manually check distrWebApp.jl for accepted keys.")
end
if haskey(myfile,"og")
    og = myfile["og"]
    E_array = vec(og.energy)
    pm_array = vec(og.pitch)
    Rm_array = vec(og.r)
else
    E_array = myfile["E_array"]
    pm_array = myfile["pm_array"]
    Rm_array = myfile["Rm_array"]
end
if haskey(myfile,"FI_species")
    FI_species = myfile["FI_species"]
end

close(myfile)

## ------
# Trim distribution data if necessary
if length(size(F_os_3D))==4
    F_os_3D = F_os_3D[id,:,:,:]
end

## ------
# Load the other 3D orbit-space fast-ion distribution, the slices of which will be visualized
if !(filepath_distr_2===nothing)
    println("Loading the second 3D orbit-space fast-ion distribution... ")
    myfile = jldopen(filepath_distr_2,false,false,false,IOStream)
    mykeys = keys(myfile)
    psi3D_2 = false # Please ignore. This will be switched on automatically.
    psi_value_2 = 0.6
    if "F_os_3D" in mykeys
        F_os_3D_2 = myfile["F_os_3D"]
    elseif "psi_os_3D" in mykeys
        F_os_3D_2 = myfile["psi_os_3D"]
        Rz_array_2 = myfile["Rz_array"]
        psi3D_2 = true
    elseif "F_reconstructions_3D" in mykeys
        F_os_3D_2 = myfile["F_reconstructions_3D"]
    else
        error("Error! Unknown orbit-space 3D object in 'filepath_distr_2'. Please manually check distrWebApp.jl for accepted keys.")
    end
    if haskey(myfile,"og")
        og_2 = myfile["og"]
        E_array_2 = vec(og_2.energy)
        pm_array_2 = vec(og_2.pitch)
        Rm_array_2 = vec(og_2.r)
    else
        E_array_2 = myfile["E_array"]
        pm_array_2 = myfile["pm_array"]
        Rm_array_2 = myfile["Rm_array"]
    end
    if haskey(myfile,"FI_species")
        FI_species = myfile["FI_species"]
    end
    close(myfile)
end

## ----------
# Double-check that both loaded fast-ion distributions have the same E-, pm- and Rm-arrays
if !(filepath_distr_2===nothing)
    # Trim second distribution data if necessary
    if length(size(F_os_3D_2))==4
        F_os_3D_2 = F_os_3D_2[id_2,:,:,:]
    end
    if !((E_array==E_array_2) && (pm_array==pm_array_2) && (Rm_array==Rm_array_2))
        error("Orbit-space grid arrays do not match for the two specified fast-ion distributions. Please correct and re-try.")
    end
end
# Check that the size of the null-orbit orbit grid matches the size of F_os_3D (and, by extension, the size of F_os_3D_2, if defined)
if showNo && (@isdefined og_nb)
    if !(size(og_nb.orbit_index)==size(F_os_3D))
        error("Sizes of orbit grids used for null orbits and fast-ion distribution(s) do not match. Please correct and re-try.")
    end
end

## ----------
# Loading tokamak equilibrium
verbose && println("Loading magnetic equilibrium... ")
M, wall, jdotb = nothing, nothing, nothing # Initialize global magnetic equilibrium variables
try
    global M; global wall; global jdotb; global timepoint#; global timepoint_source # Declare global scope
    M, wall = read_geqdsk(filepath_equil,clockwise_phi=false) # Assume counter-clockwise phi-direction
    jdotb = M.sigma # The sign of the dot product between the plasma current and the magnetic field
catch # Otherwise, assume magnetic equilibrium is a saved .jld2 file
    global M; global wall; global jdotb; global timepoint#; global timepoint_source
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

## ----------
# Check that the dimensions of topoBounds and F_os_3D match, and match the E-, pm- and Rm-arrays.
println("Checking that dimensions are consistent... ")
if !(size(F_os_3D)==size(topoBounds)) || !(size(F_os_3D,1) == length(E_array)) || !(size(F_os_3D,2) == length(pm_array)) || !(size(F_os_3D,3) == length(Rm_array))
    println("Size(F_os_3D/psi_os_3D): $(size(F_os_3D))")
    println("Size(topoBounds): $(size(topoBounds))")
    println("Length(E_array): $(length(E_array))")
    println("Length(pm_array): $(length(pm_array))")
    println("Length(Rm_array): $(length(Rm_array))")
    error("Dimensions of topoBounds and/or F_os_3D don't match given E-, pm- and/or Rm-arrays. Please correct and re-try.")
end

if psi3D
    println("Unzipping Rz_array... ")
    psiR_array = zeros(size(Rz_array))
    psiz_array = zeros(size(Rz_array))
    for i=1:length(Rz_array)
        psiR_array[i] = (Rz_array[i])[1]
        psiz_array[i] = (Rz_array[i])[2]
    end
end

if !(filepath_distr_2===nothing)
    if psi3D_2
        println("Unzipping second Rz_array... ")
        psiR_array_2 = zeros(size(Rz_array_2))
        psiz_array_2 = zeros(size(Rz_array_2))
        for i=1:length(Rz_array_2)
            psiR_array_2[i] = (Rz_array_2[i])[1]
            psiz_array_2[i] = (Rz_array_2[i])[2]
        end
    end
end

verbose && println("Size of pm_array: $(size(pm_array))")
verbose && println("Size of Rm_array: $(size(Rm_array))")
verbose && println("Size of E_array: $(size(E_array))")
verbose && println("Size of F_os_3D(/psi_os_3D): $(size(F_os_3D))")
verbose && println("Size of topoBounds: $(size(topoBounds))")

## ------
# Make use of the loaded null-orbits data and prepare it for app visualization format
if showNo && (@isdefined nullOrbs_indices)
    verbose && println("Creating 3D null-orbit matrix from loaded 4D null-orbit data... ")
    nullOrbs = zeros(size(F_os_3D))
    for i=1:length(nullOrbs_indices)
        nullOrbs[nullOrbs_indices[i][2],nullOrbs_indices[i][3],nullOrbs_indices[i][4]] = 1.0
    end
end
if showNo && (@isdefined og_nb)
    verbose && println("Creating 3D null-orbit matrix from loaded 1D null-orbit data and pertaining orbit grid... ")
    nullOrbs = map_orbits(og_nb, nullOrbs_1D, true)
end

E_array = vec(collect(E_array)) # Ensure type Array{Float64,1}
pm_array = vec(collect(pm_array)) # Ensure type Array{Float64,1}
Rm_array = vec(collect(Rm_array)) # Ensure type Array{Float64,1}


## ------
# Pre-compute f(E), if needed
if !(filepath_distr_2===nothing)
    verbose && println("Pre-computing f(E) for visualization and comparison of the two fast-ion distributions... ")
    dE3D, dpm3D, dRm3D = OrbitTomography.get3DDiffs(E_array,pm_array,Rm_array)
    fE_orig = dropdims(sum(F_os_3D .*dpm3D .*dRm3D,dims=(2,3)),dims=(2,3))
    fE_2_orig = dropdims(sum(F_os_3D_2 .*dpm3D .*dRm3D,dims=(2,3)),dims=(2,3))

    fE = fE_orig ./ (maximum(fE_orig)== 0.0 ? 1.0 : maximum(fE_orig))
    fE_2 = fE_2_orig ./ (maximum(fE_2_orig)== 0.0 ? 1.0 : maximum(fE_2_orig))

    fE_min, fE_max = extrema(fE)
    fE_2_min, fE_2_max = extrema(fE_2)
    if (round(log10(fE_max))-round(log10(fE_min)))!=0.0 || (round(log10(fE_2_max))-round(log10(fE_2_min)))!=0.0
        y_scale = :log10
        gi_fE = findall(x-> x>0.0, fE)
        gi_fE_2 = findall(x-> x>0.0, fE_2)
    else
        y_scale = :identity
        gi_fE = 1:length(fE)
        gi_fE_2 = 1:length(fE_2)
    end
end
## ------
# Minimizing memory usage
og_nb = nothing
og = nothing
og_2 = nothing
# Add more variables that won't be used again here...

## ------
# The web application
R_hfs = minimum(wall.r) # R-coord of high-field side wall
R_lfs = maximum(wall.r) # R-coord of low-field side wall
phi = collect(0:1:359).*(2*pi/180.0) # Toroidal angle
topview_R_hfs_x = (R_hfs).*cos.(phi)
topview_R_hfs_y = (R_hfs).*sin.(phi)
topview_R_lfs_x = (R_lfs).*cos.(phi)
topview_R_lfs_y = (R_lfs).*sin.(phi)
verbose && println("--- You can access the distrWebApp via an internet web browser when you see 'Task (runnable)...' ")
verbose && println("--- When 'Task (runnable)...' has appeared, please visit the website localhost:$(port) ---")
verbose && println("--- Remember: It might take a minute or two to load the webpage. Please be patient. ---")
function app(req)
    @manipulate for tokamak_wall = Dict("on" => true, "off" => false), colorbar_scale = Dict("0.0-1.0" => "zero2one", "0.0-0.1" => "zero2aTenth"), compare_distributions = Dict("On" => true, "Off" => false), E=E_array, pm=pm_array, Rm=Rm_array, save_plots = Dict("on" => true, "off" => false)

        EPRc = EPRCoordinate(M, E, pm, Rm; amu=getSpeciesAmu(FI_species), q=getSpeciesEcu(FI_species))
        o = get_orbit(M,EPRc; wall=wall, extra_kw_args...)
        topview_o_x = cos.(o.path.phi).*(o.path.r)
        topview_o_y = sin.(o.path.phi).*(o.path.r)

        orb_color = :black
        orb_linestyle = :solid
        if o.class==:invalid
            orb_color = :gray
            orb_linestyle = :dash
        elseif o.class == :lost
            orb_color = :brown
        elseif o.class == :incomplete # If this happens, you are unfortunately in trouble. Because it will likely take forever to calculate. Please just re-start the app instead.
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
        elseif (o.class == :stagnation && o.coordinate.r<magnetic_axis(M)[1]) # Pinch (HFS stagnation) orbit
            orb_color = :pink
        else
            error("Something's gone wrong!!! Orbit class unknown!")
        end

        # Cross-sectional plot
        plt_crs = Plots.scatter([magnetic_axis(M)[1]],[magnetic_axis(M)[2]],label="Magnetic axis", mc=:grey, aspect_ratio=:equal, xlabel="R [m]", ylabel=" z[m]")
        plt_crs = Plots.plot!(o.path.r,o.path.z, label="$(o.class) orbit", color=orb_color, linestyle=orb_linestyle, linewidth=1.5,title="E: $(round(E, digits=4)) keV  pm: $(round(o.coordinate.pitch, digits=2))  Rm: $(round(o.coordinate.r,digits=2))")
        if tokamak_wall
            plt_crs = Plots.contour!(flux_r,flux_z,psi_rz',levels=collect(range(psi_mag,stop=psi_bdry,length=5)),color=:gray, linestyle=:dot,linewidth=1.5, label="",colorbar=false)
            plt_crs = Plots.plot!(wall.r,wall.z, label="JET wall", color=:black, linewidth=1.5,xaxis=[1.5,5.5])
        end
        plt_crs = Plots.scatter!([o.coordinate.r],[o.coordinate.z], mc=orb_color, label="(Rm,zm)")
        if psi3D
            plt_crs = Plots.scatter!(psiR_array,psiz_array, label="Flux surface (Ψ=$(psi_value))", mc=:gray,ms=1.5)
        end
        if !(filepath_distr_2===nothing)
            if psi3D_2 && compare_distributions
                plt_crs = Plots.scatter!(psiR_array_2,psiz_array_2, label="Flux surface (Ψ=$(psi_value))", mc=:gray13,ms=1.5)
            end
        end
        if save_plots
            png(plt_crs, "plt_crs_$(round(E, digits=2))_$(round(o.coordinate.pitch, digits=2))_$(round(o.coordinate.r,digits=2))")
        end
            
        # Find the E in E_array. For fast-ion slice plot.
        Ei = (findall(x-> x==E,E_array))[1] # Should only be 1 exactly equal element

        # Extract correct topological boundaries and convert to vectors (for scatter plot)
        ones_carinds = findall(x-> x==1.0,topoBounds[Ei,:,:])
        pm_scatvals_tb = zeros(length(ones_carinds))
        Rm_scatvals_tb = zeros(length(ones_carinds))
        for (ind,carinds) in enumerate(ones_carinds)
            pm_scatvals_tb[ind] = pm_array[carinds[1]]
            Rm_scatvals_tb[ind] = Rm_array[carinds[2]]
        end

        # Extract correct null boundaries and convert to vectors (for scatter plot)
        if showNo && (@isdefined nullOrbs)
            ones_carinds = findall(x-> x==1.0,nullOrbs[Ei,:,:])
            pm_scatvals_nb = zeros(length(ones_carinds))
            Rm_scatvals_nb = zeros(length(ones_carinds))
            for (ind,carinds) in enumerate(ones_carinds)
                pm_scatvals_nb[ind] = pm_array[carinds[1]]
                Rm_scatvals_nb[ind] = Rm_array[carinds[2]]
            end
        end

        if colorbar_scale == "zero2one"
            clims = (0.0,1.0)
        else
            clims = (0.0,0.1)
        end

        #topview plot
        plt_top = Plots.plot(topview_R_lfs_x,topview_R_lfs_y, label="JET wall", color=:black, linewidth=1.5)
        plt_top = Plots.plot!(topview_R_hfs_x,topview_R_hfs_y, label="", color=:black,linewidth=1.5, aspect_ratio=:equal, title="Top view")
        plt_top = Plots.plot!(topview_o_x,topview_o_y,label="$(o.class) orbit", color=orb_color, linestyle=orb_linestyle, linewidth=1.5)
        if save_plots
            png(plt_top, "plt_top_$(round(E, digits=2))_$(round(o.coordinate.pitch, digits=2))_$(round(o.coordinate.r,digits=2))")
        end

        # distribution plot
        plt_distr = Plots.heatmap(Rm_array, pm_array, (F_os_3D[Ei,:,:])./maximum(F_os_3D[Ei,:,:]), colorbar=true, title="Fast-ion distribution slice  ($(round(maximum(F_os_3D[Ei,:,:]), sigdigits=4)) = 1.0)", clims=clims,fillcolor=cgrad([:white, :darkblue, :green, :yellow, :orange, :red]))
        plt_distr = Plots.scatter!(Rm_scatvals_tb,pm_scatvals_tb,markersize=1.8,leg=false,markercolor=:white, xlabel="Rm [m]", ylabel="pm")
        if showNo && (@isdefined nullOrbs)
            plt_distr = Plots.scatter!(Rm_scatvals_nb,pm_scatvals_nb,markersize=5,leg=false,markercolor=:teal, markershape=:circle)
        end
        plt_distr = Plots.scatter!([Rm],[pm],markershape=:circle,mc=orb_color,markersize=5.0) # orbit coordinate marker
        if save_plots
            png(plt_distr, "plt_distr_$(round(E, digits=2))_$(round(o.coordinate.pitch, digits=2))_$(round(o.coordinate.r,digits=2))")
        end

        #pitch visualization plot, or second distribution plot. Depending on input
        if !(filepath_distr_2===nothing) && compare_distributions
            # distribution plot
            plt_pitc = Plots.heatmap(Rm_array, pm_array, (F_os_3D_2[Ei,:,:])./maximum(F_os_3D_2[Ei,:,:]), colorbar=true, title="Fast-ion reference slice  ($(round(maximum(F_os_3D_2[Ei,:,:]), sigdigits=4)) = 1.0)", clims=clims,fillcolor=cgrad([:white, :darkblue, :green, :yellow, :orange, :red]))
            plt_pitc = Plots.scatter!(Rm_scatvals_tb,pm_scatvals_tb,markersize=1.8,leg=false,markercolor=:white, xlabel="Rm [m]", ylabel="pm")
            if showNo && (@isdefined nullOrbs)
                plt_pitc = Plots.scatter!(Rm_scatvals_nb,pm_scatvals_nb,markersize=5,leg=false,markercolor=:teal, markershape=:circle)
            end
            plt_pitc = Plots.scatter!([Rm],[pm],markershape=:circle,mc=orb_color,markersize=5.0) # orbit coordinate marker
            if save_plots
                png(plt_pitc, "plt_distr2_$(round(E, digits=2))_$(round(o.coordinate.pitch, digits=2))_$(round(o.coordinate.r,digits=2))")
            end
        else
            # pitch visualization plot
            xmin = (pm >= 0.0) ? 0.0 : pm
            xp = pm
            yp = sqrt(1-pm^2)
            plt_pitc = Plots.quiver([0.0],[0.0],quiver=([1.0],[0.0]), label="B field", markerstrokewidth=4.0, color=:purple)
            plt_pitc = Plots.plot!([0.0,1.0],[0.0,0.0], label="B field", linewidth=3, color=:purple, annotations=(0.98,0.05,Plots.text("B",:right)))
            plt_pitc = Plots.quiver!([0.0],[0.0],quiver=([xp],[yp]), label="v", markerstrokewidth=4.0, color=:black)
            plt_pitc = Plots.plot!([0.0,xp],[0.0,yp], label="v", linewidth=3, color=:black, annotations=(xp+0.02,yp+0.03,Plots.text("v")))
            plt_pitc = Plots.quiver!([0.0],[0.0],quiver=([0.0],[yp]), label="vperp", markerstrokewidth=4.0, color=:red)
            plt_pitc = Plots.plot!([0.0,0.0],[0.0,yp], label="vperp", linewidth=3, color=:red, annotations=(0.02,yp+0.03,Plots.text("vperp")))
            plt_pitc = Plots.quiver!([0.0],[0.0],quiver=([xp],[0.0]), label="vpara", markerstrokewidth=4.0, color=:blue)
            plt_pitc = Plots.plot!([0.0,xp],[0.0,0.0], label="vpara", xaxis=[xmin,1.0], yaxis=[0.0,1.0], linewidth=3, legend=true, color=:blue, fontsize=14.0, fontweight="bold", annotations=(xp,0.05,Plots.text("vpara")))
            plt_pitc = Plots.plot!([0.0,xp],[yp,yp],linestyle=:dash,label="",color=:black)
            plt_pitc = Plots.plot!([xp,xp],[0.0,yp],linestyle=:dash,label="",color=:black, title="pm: $(round(pm, digits=3))")
            if save_plots
                png(plt_pitc, "plt_pitc_$(round(E, digits=2))_$(round(o.coordinate.pitch, digits=2))_$(round(o.coordinate.r,digits=2))")
            end
        end

        if !(filepath_distr_2===nothing) && compare_distributions
            plt_fE = Plots.plot(E_array[gi_fE],fE[gi_fE],linewidth=3.0,color=:blue,label="f(E)", yaxis=y_scale)
            plt_fE = Plots.plot!(E_array[gi_fE_2],fE_2[gi_fE_2],linewidth=3.0,color=:black,label="f(E) reference",xlabel="Energy [keV]",ylabel="f(E) [keV^-1]")
            plt_fE = Plots.scatter!([E_array[Ei]],[fE[Ei]],ms=5.0,mc=:blue,label="$(round(fE_orig[Ei],digits=2))")
            plt_fE = Plots.scatter!([E_array[Ei]],[fE_2[Ei]],ms=5.0,mc=:black,label="$(round(fE_2_orig[Ei],digits=2))", title="Normalized f(E) distributions")
            plt_fE = Plots.plot!(size=(1200,400),left_margin=5Plots.mm,bottom_margin=5Plots.mm)
        end
        
        if !(filepath_distr_2===nothing) && compare_distributions
            vbox(vskip(1em),
                plt_fE,
                hbox(plt_distr,plt_pitc),
                hbox(plt_top,plt_crs)
            )
        else
            vbox(vskip(1em),
                hbox(plt_distr,plt_pitc),
                hbox(plt_top,plt_crs)
            )
        end
    end
end

webio_serve(page("/",app), port)