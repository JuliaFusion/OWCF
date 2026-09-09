### A Pluto.jl notebook ###
# v1.0.3

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ f24dc2c0-a139-11f1-309d-21cee6dd4826
md"""
# orbitWebApp

## Description:
This notebook provides an application to visualize a single (guiding-centre) orbit and its path around the tokamak in an interactive and intuitive manner. The (E,pm,Rm) coordinate is specified and the orbit can be visualized via sliders. 

At the bottom of the web application, there is a cell that allows the user to save a .gif file of an animation that shows how the particle moves along the orbit.

## Inputs:
- folderpath_OWCF - The path to the OWCF folder on your computer. Needed for correct loading - String
- filepath_equil - The path to the .eqdsk-file (or .jld2-file) with the tokamak magnetic equilibrium and tokamak wall geometry. A .jld2 file with a magnetic equilibrium and tokamak wall geometry can be computed using the OWCF/extra/createCustomMagneticEquilibrium.jl script - String
- FI\_species - The species of the particle being simulated (deuterium, tritium etc). Specified as "D", "T" etc. Please see the keys of the OWCF\_ChemElemSymbol\_to\_Z Dict object in the file OWCF/misc/species_func.jl for a full list of all available particle species in the OWCF - String
- verbose - If set to true, the app will talk a lot! - Bool

## Outputs:
# -

## Saved files:
# -

### Original script written by Henrik Järleblad, henrikj@dtu.dk
### Last maintained 2026-09-04
"""

# ╔═╡ ca2c21a7-429e-4be4-ae40-9e1d0bf35092
begin
	# SPECIFY THE INPUTS IN THIS CELL

	# Please specify the OWCF folder and let the notebook change directory to the 
	# OWCF folder when the cell below is executed. This is to be able to load the
	# correct versions of the Julia packages as specified in the Project.toml and 
	# Manifest.toml files.
	folderpath_OWCF = "REPLACE-THIS-TEXT-WITH-THE-PATH-TO-THE-OWCF-FOLDER-ON-YOUR-COMPUTER" # Finish with '/'
	
	filepath_equil = folderpath_OWCF*"equilibrium/JET/g99971/g99971_474-48.9.eqdsk" # Or .jld2. The g99971_474-48.9.eqdsk file is an example.
	FI_species = "D" # Example deuterium: "D"
	verbose = true
	extra_kw_args = Dict(:toa => true, :limit_phi => true) # Extra keyword arguments for the orbit-integration algorithm. toa is 'try only adaptive' and limit_phi limits the number of toroidal turns for orbits
end

# ╔═╡ 1b10b234-a468-4772-a4a4-7771b1f2b445
begin
	# Load the OWCF environment and all necessary Julia packages
	cd(folderpath_OWCF)
	using Pkg
	Pkg.activate(".")
	
	verbose && println("Loading packages... ")
	using EFIT
	using Equilibrium
	using GuidingCenterOrbits
	using Plots
	using PlutoUI
	using JLD2
	using FileIO
	using Suppressor
	include(folderpath_OWCF*"misc/species_func.jl")
	include(folderpath_OWCF*"extra/gui.jl") # For orbit movie
	println()
end

# ╔═╡ 88a35068-f22c-451a-8a8d-0b35a8d5ed6f
begin
	verbose && println("Loading magnetic equilibrium... ")
    M, wall, jdotb = nothing, nothing, nothing # Initialize global magnetic equilibrium variables
    try
        global M; global wall; global jdotb # Declare global scope
        M, wall = read_geqdsk(filepath_equil,clockwise_phi=false) # Assume counter-clockwise phi-direction
        jdotb = M.sigma # The sign of the dot product between the plasma current and the magnetic field
    catch # Otherwise, assume magnetic equilibrium is a saved .jld2 file
        global M; global wall; global jdotb; local myfile # Declare global scope and local scope for variables
        myfile = jldopen(filepath_equil,false,false,false,IOStream)
        M = myfile["S"]
        wall = myfile["wall"]
        close(myfile)
        jdotb = (M.sigma_B0)*(M.sigma_Ip)
    end
    
    verbose && println("Computing flux function on 100x100 (R,z)-grid (to plot flux surfaces)... ")
    flux_r = range(extrema(wall.r)...,length=100)
    flux_z = range(extrema(wall.z)...,length=100)
    inds = CartesianIndices((length(flux_r),length(flux_z)))
    psi_rz = [M(flux_r[ind[1]], flux_z[ind[2]]) for ind in inds]
    psi_mag, psi_bdry = psi_limits(M)
    
    verbose && println("Defining necessary quantities... ")
    R_hfs = minimum(wall.r) # R-coord of high-field side wall
    R_lfs = maximum(wall.r) # R-coord of low-field side wall
    phi = collect(0:1:359).*(2*pi/180.0) # Toroidal angle
    topview_R_hfs_x = (R_hfs).*cos.(phi)
    topview_R_hfs_y = (R_hfs).*sin.(phi)
    topview_R_lfs_x = (R_lfs).*cos.(phi)
    topview_R_lfs_y = (R_lfs).*sin.(phi)
end

# ╔═╡ 8c351874-c009-4462-97d3-2df97de824ea
function plot_inputs()

	return PlutoUI.combine() do Child

		inputs=[
			md""" Show tokamak wall: $(
				Child("wall", Switch(default=true))
			)""",
			md""" Energy (keV): $(
				Child("E", Slider(1:1:1000.0, default=100.0))
			)""",
			md""" Pitch maximum, pm (-): $(
				Child("pm", Slider(-1.0:0.01:1.0, default=0.7))
			)""",
			md""" Radius maximum, Rm (m): $(
				Child("Rm", Slider(((4*M.axis[1]+minimum(wall.r))/5):0.01:maximum(wall.r), default=(maximum(wall.r)+magnetic_axis(M)[1])/2))
			)""",
			md""" Orbit path index, i (-): $(
				Child("i", Slider(1:1:500, default=480))
			)"""
		]

		md"""
		#### Plot controls
		$(inputs)
		"""
	end
end

# ╔═╡ fa5b0903-83d8-40da-8f97-97f711ea60af
@bind my_plot_inputs plot_inputs()

# ╔═╡ 26eee712-0e9e-4415-9534-65db1e734d34
let
	@suppress begin
		tokamak_wall = my_plot_inputs.wall
		E = my_plot_inputs.E
		pm = my_plot_inputs.pm
		Rm = my_plot_inputs.Rm
		i = my_plot_inputs.i
		
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
		elseif o.class == :incomplete 
			# If this happens, you are in trouble. Because it will likely take forever to calculate. Please just re-start the app instead.
		    orb_color = :yellow
		elseif o.class == :trapped
		    orb_color = :blue
		elseif o.class == :co_passing
		    orb_color = :green
		elseif (o.class == :stagnation && o.coordinate.r>=magnetic_axis(M)[1]) 
			# Regular stagnation orbit
		    orb_color = :red
		elseif o.class == :potato
		    orb_color = :orange
		elseif o.class == :ctr_passing
		    orb_color = :purple
		elseif (o.class == :stagnation && o.coordinate.r<magnetic_axis(M)[1]) 
			# Counter-stagnation orbit
		    orb_color = :pink
		else
		    error("Something's gone wrong!!! Orbit class unknown!")
		end

		# topview plot
		plt_top = Plots.plot(topview_o_x[1:i],topview_o_y[1:i],label="$(o.class) orbit", color=orb_color, linestyle=orb_linestyle, linewidth=1.5)
		plt_top = Plots.plot!(topview_R_lfs_x,topview_R_lfs_y, label="JET wall", color=:black, linewidth=1.5)
		plt_top = Plots.plot!(topview_R_hfs_x,topview_R_hfs_y, label="", color=:black,linewidth=1.5, aspect_ratio=:equal, title="Top view")
		plt_top = Plots.scatter!([topview_o_x[i]], [topview_o_y[i]], label="", mc=orb_color, xlabel="x [m]", ylabel="y [m]")

		# cross-sectional plot
		plt_crs = Plots.scatter([magnetic_axis(M)[1]],[magnetic_axis(M)[2]],label="Magnetic axis", mc=:gray, aspect_ratio=:equal, xlabel="R [m]", ylabel=" z[m]", title="E: $(round(E,digits=2)) keV  pm: $(round(o.coordinate.pitch, digits=2))  Rm: $(round(o.coordinate.r,digits=2))")
		plt_crs = Plots.plot!(o.path.r[1:i],o.path.z[1:i], label="$(o.class) orbit", color=orb_color, linestyle=orb_linestyle, linewidth=1.5)
		if tokamak_wall
		    plt_crs = Plots.contour!(flux_r,flux_z,psi_rz',levels=collect(range(psi_mag,stop=psi_bdry,length=5)),color=:gray, linestyle=:dot,linewidth=1.5, label="",colorbar=false)
		    plt_crs = Plots.plot!(wall.r,wall.z, label="JET wall", color=:black, linewidth=1.5)
		end
		plt_crs = Plots.scatter!([o.path.r[i]],[o.path.z[i]],mc=orb_color,label="")

		# pitch visualization plot
		t_array = collect(range(0.0,stop=((o.tau_p)*(i/length(o.path.pitch))), length=i)) ./(1.0e-6) # Microseconds
		plt_pitc = Plots.plot(t_array, o.path.pitch[1:i],color=orb_color, title="Pitch (p) along orbit path \n p=$(round(o.path.pitch[i],sigdigits=3))", label="", xlabel="Poloidal time [microseconds]")
		plt_pitc = Plots.scatter!([t_array[end]],[o.path.pitch[i]],color=orb_color,label="", ylabel="pitch [-]")

		plt_dum = Plots.plot(axis=([],false), aspect_ratio=:equal)

		Plots.plot(plt_crs, plt_pitc, plt_top, plt_dum, layout=(2,2), size=(1000, 800))
	end
end

# ╔═╡ 704e7216-b629-4975-97f1-2d4141256cbc
let
	# AN EXTRA CELL THAT LET'S YOU SAVE AN ANIMATION OF HOW THE PARTICLE ACTUALLY MOVES ALONG THE ORBIT. TO SAVE THE ANIMATION
	# - Set the (E,pm,Rm) values below
	# - Set the 'save_animation' variable to 'true'
	# - Remember to set the 'save_animation' variable to 'false' once you are done with the app, or when you don't want to save the animation
	E = 100.0 # keV
	pm = 0.3
	Rm = 3.47 # meters
	save_animation = true # Change this to 'true' to save animation when cell is run

	@suppress begin
		plot_orbit_movie(M, E, pm, Rm; FI_species=FI_species, wall=wall, save_anim=save_animation)
	end
end

# ╔═╡ Cell order:
# ╠═f24dc2c0-a139-11f1-309d-21cee6dd4826
# ╠═ca2c21a7-429e-4be4-ae40-9e1d0bf35092
# ╠═1b10b234-a468-4772-a4a4-7771b1f2b445
# ╠═88a35068-f22c-451a-8a8d-0b35a8d5ed6f
# ╠═8c351874-c009-4462-97d3-2df97de824ea
# ╠═fa5b0903-83d8-40da-8f97-97f711ea60af
# ╠═26eee712-0e9e-4415-9534-65db1e734d34
# ╠═704e7216-b629-4975-97f1-2d4141256cbc
