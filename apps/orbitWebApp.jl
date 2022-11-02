################################ orbitWebApp.jl #########################################

#### Description:
# This script provides an application to visualize a single (guiding-centre) orbit and its
# path around the tokamak in an interactive and intuitive manner. The (E,pm,Rm) coordinate
# is specified and the orbit can be visualize via sliders. 

# Prior to running this script, please make sure you have run the 
# calcTopoMap.jl script (or equivalent). Also, please make sure that you have noted the 
# path to the saved topological map.

#### Inputs (units given when defined in the script):
# folderpath_OWCF - The path to the OWCF folder on your computer. Needed for correct loading - String
# port - The I/O internet port through which this app will be accessed - Int64
# filepath_equil - The path to the .eqdsk-file (or .jld2-file) with the tokamak magnetic equilibrium and geometry - String
# FI_species - The species of the particle being simulated (deuterium, tritium etc) - String
# verbose - If set to true, the app will talk a lot! - Bool

#### Outputs
# -

#### Saved files
# -

# Script written by Henrik Järleblad. Last maintained 2022-08-25.
#########################################################################################

## --------------------------------------------------------------------------
# Please specify the OWCF folder and let the app change directory to the 
# OWCF folder when orbitWebApp.jl is executed. This is to be able to load the
# correct versions of the Julia packages as specified in the Project.toml and 
# Manifest.toml files.
folderpath_OWCF = "" # Finish with '/'
cd(folderpath_OWCF)
using Pkg
Pkg.activate(".")

## ------
# Inputs
port = 4444 # Should not be changed lightly
filepath_equil = ""
FI_species = "" # Example deuterium: "D"
verbose = true

extra_kw_args = Dict(:toa => true, :limit_phi => true, :maxiter => 0)
# toa is 'try only adaptive'
# limits the number of toroidal turns for orbits
# The orbit integration algorithm will try progressively smaller timesteps these number of times

## ------
# Loading packages
verbose && println("Loading packages... ")
using Interact
using EFIT
using Equilibrium
using GuidingCenterOrbits
using Plots
using JLD2
using FileIO
using Mux
using WebIO
include(folderpath_OWCF*"misc/species_func.jl")

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
# The web application
verbose && println("Building web application... ")
R_hfs = minimum(wall.r) # R-coord of high-field side wall
R_lfs = maximum(wall.r) # R-coord of low-field side wall
phi = collect(0:1:359).*(2*pi/180.0) # Toroidal angle
topview_R_hfs_x = (R_hfs).*cos.(phi)
topview_R_hfs_y = (R_hfs).*sin.(phi)
topview_R_lfs_x = (R_lfs).*cos.(phi)
topview_R_lfs_y = (R_lfs).*sin.(phi)
verbose && println("--- You can access the orbitWebApp via an internet web browser when you see 'Task (runnable)...' ")
verbose && println("--- When 'Task (runnable)...' has appeared, please visit the website localhost:$(port) ---")
verbose && println("--- Remember: It might take a minute or two to load the webpage. Please be patient. ---")
function app(req)
    @manipulate for tokamak_wall = Dict("on" => true, "off" => false), E=88.0, pm=0.5, Rm=(maximum(wall.r)+magnetic_axis(M)[1])/2, i=1:1:500, save_plots = Dict("on" => true, "off" => false)

        EPRc = EPRCoordinate(M, E, pm, Rm, amu=getSpeciesAmu(FI_species), q=getSpeciesEcu(FI_species))
        o = get_orbit(M,EPRc; wall=wall, interp_dt=1.0e-10, max_length=500, extra_kw_args...) # interp_dt is set to ridiculously small value, to ensure orbit path length of 500
        
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

        #topview plot
        plt_top = Plots.plot(topview_o_x[1:i],topview_o_y[1:i],label="$(o.class) orbit", color=orb_color, linestyle=orb_linestyle, linewidth=1.5)
        plt_top = Plots.plot!(topview_R_lfs_x,topview_R_lfs_y, label="JET wall", color=:black, linewidth=1.5)
        plt_top = Plots.plot!(topview_R_hfs_x,topview_R_hfs_y, label="", color=:black,linewidth=1.5, aspect_ratio=:equal, title="Top view")
        plt_top = Plots.scatter!([topview_o_x[i]],[topview_o_y[i]],label="",mc=orb_color, xlabel="x [m]", ylabel="y [m]")
        if save_plots
            png(plt_top, "plt_top_$(round(E, digits=2))_$(round(o.coordinate.pitch, digits=2))_$(round(o.coordinate.r,digits=2))")
        end

        #cross-sectional plot
        plt_crs = Plots.scatter([magnetic_axis(M)[1]],[magnetic_axis(M)[2]],label="Magnetic axis", mc=:gray, aspect_ratio=:equal, xlabel="R [m]", ylabel=" z[m]",title="E: $(round(E,digits=2)) keV  pm: $(round(o.coordinate.pitch, digits=2))  Rm: $(round(o.coordinate.r,digits=2))")
        plt_crs = Plots.plot!(o.path.r[1:i],o.path.z[1:i], label="$(o.class) orbit", color=orb_color, linestyle=orb_linestyle, linewidth=1.5)
        if tokamak_wall
            plt_crs = Plots.contour!(flux_r,flux_z,psi_rz',levels=collect(range(psi_mag,stop=psi_bdry,length=5)),color=:gray, linestyle=:dot,linewidth=1.5, label="",colorbar=false)
            plt_crs = Plots.plot!(wall.r,wall.z, label="JET wall", color=:black, linewidth=1.5)
        end
        plt_crs = Plots.scatter!([o.path.r[i]],[o.path.z[i]],mc=orb_color,label="")
        if save_plots
            png(plt_crs, "plt_crs_$(round(E, digits=2))_$(round(o.coordinate.pitch, digits=2))_$(round(o.coordinate.r,digits=2))")
        end

        #pitch visualization plot
        t_array = collect(range(0.0,stop=((o.tau_p)*(i/length(o.path.pitch))), length=i)) ./(1.0e-6) # Microseconds
        plt_pitc = Plots.plot(t_array, o.path.pitch[1:i],color=orb_color,title="Pitch (p) along orbit path \n p=$(round(o.path.pitch[i],sigdigits=3))",label="", xlabel="Poloidal time [microseconds]")
        plt_pitc = Plots.scatter!([t_array[end]],[o.path.pitch[i]],color=orb_color,label="", ylabel="pitch [-]")
        if save_plots
            png(plt_pitc, "plt_pitc_$(round(E, digits=2))_$(round(o.coordinate.pitch, digits=2))_$(round(o.coordinate.r,digits=2))")
        end

        vbox(vskip(1em),
            hbox(Plots.plot(plt_top),Plots.plot(plt_crs)),
            hbox(Plots.plot(plt_pitc))
        )
    end
end
webio_serve(page("/",app), port)
