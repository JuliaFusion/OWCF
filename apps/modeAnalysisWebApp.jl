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

# ╔═╡ e1984f36-a6a7-11f1-16fc-a9fcbbb60876
md"""
# modeAnalysisWebApp

## Description:
This notebook provides a web application to visualize resonances between fast ions and magnetohydrodynamic (MHD) Alfvén eigenmodes in an interactive and intuitive manner. The resonance criterion is based on the equation

ω = n ω\_ϕ + m ω\_θ + ω\_rot

where ω is the mode (angular) frequency, n is the toroidal mode number, ω\_ϕ is the toroidal transit (angular) frequency of the fast ion, m is the poloidal mode number, ω\_θ is the poloidal transit (angular) frequency of the fast ion and ω\_rot is the plasma rotation (angular) frequency. ω\_ϕ and ω\_θ are computed by integrating the guiding-centre equations of motion, and then evaluating the toroidal (τ\_ϕ) and poloidal (τ\_θ) transit times. ω, n and m can be interactively changed by the user using the sliders in this web application.

τ\_ϕ and τ\_θ are related to ω\_ϕ and ω\_θ as ω\_ϕ = 2\*pi/τ\_ϕ and ω\_θ = 2\*pi/τ\_θ, respectively.

The resonance between mode and fast ion is then plotted in (E,pm,Rm) (or, via the press of a button in the app, in (E,Λ,Pϕ\_n;σ) CURRENTLY OUT-OF-ORDER). The resonance is visualized as a heatmap where the values (A) of the figure colorbar are computed according to the equation

A = - log10(|ω - 2πn/τ\_ϕ - 2πm/τ\_θ - ω\_rot|)

where larger values of A thus corresponds to greater resonance, and smaller values of A corresponds to less resonance.

### PLEASE NOTE! 
To be able to use the app, a topological map containing toroidal (τ\_ϕ) and poloidal (τ_θ) transit time data has to be provided. This type of data is to be provided as a .jld2 file with the keys 'torTransTimes' and 'polTransTimes'. Such a file can be easily obtained by using the OWCF/calcTopoMap.jl tool.

## Inputs:
- enable\_COM - If true, (E,pm,Rm) -> (E,Λ,Pϕ\_n;σ) can be performed via a toggle button. Set to false to minimize computation time. CURRENTLY OUT-OF-ORDER - Bool
- folderpath\_OWCF - The path to the OWCF folder on your computer. Needed for correct loading - String
- filepath\_equil - The path to the file with the tokamak magnetic equilibrium and geometry. Can be either an .eqdsk file or an output of the OWCF/extra/createCustomMagneticEquilibrium.jl tool - String
- filepath\_tm - The path to the .jld2 file containing the topological map (and toroidal and poloidal transit times). The .jld2 file is an output of the OWCF/calcTopoMap.jl script - String
- FI\_species - The species of the particle being simulated, e.g. "D" or "T" (deuterium, tritium etc). Please see OWCF/misc/species_func.jl for a list of all available particle species in the OWCF - String
- ω\_init - The mode frequency to investigate (in kHz). This can be interactively changed later in the app - Float64
- ω\_rot\_init - The plasma rotation frequency (in kHz). This can be interactively changed later in the app - Float64
- n - The toroidal mode number to investigate. This can be interactively changed later in web application - Int64
- m - The toroidal mode number to investigate. This can be interactively changed later in web application - Int64
- verbose - If set to true, the app will talk a lot! - Bool

## Outputs:
# -

## Saved files:
# -

### Notebook written by Henrik Järleblad, henrikj@dtu.dk
### Last maintained 2026-09-02
"""

# ╔═╡ 4e590375-f4fd-47ca-b89c-fe46dbdc22c2
begin
	# SPECIFY THE INPUTS IN THIS CELL

	# Please specify the OWCF folder and let the notebook change directory to the 
	# OWCF folder when the cell below is executed. This is to be able to load the
	# correct versions of the Julia packages as specified in the Project.toml and 
	# Manifest.toml files.
	folderpath_OWCF = "REPLACE-THIS-TEXT-WITH-THE-PATH-TO-THE-OWCF-FOLDER-ON-YOUR-COMPUTER" # Finish with '/'
	
	enable_COM = false # CURRENTLY OUT-OF-ORDER!!! DO NOT SET TO 'true'!!! Set to false for large grids (>10x100x100 in (E,pm,Rm)). The computation time simply becomes too large. Or...
	if enable_COM
	    filepath_tm_COM = "" # ...please specify an output of the os2com.jl script (which contains the key "topoMap")(and possibly "polTransTimes" and "torTransTimes"). Leave unspecified if orbitsWebApp.jl should compute the (E,pm,Rm) -> (E,Λ,Pϕ_n;σ) map itself.
	end
	filepath_equil = folderpath_OWCF*"equilibrium/JET/g99971/g99971_474-48.9.eqdsk"  # Example JET shot 96100 at 13s (53 minus 40): g96100/g96100_0-53.0012.eqdsk" #
	filepath_tm = folderpath_OWCF*"apps/example_data/topoMap_JET_99971L72_at48,9s_D_6x101x102_wLost.jld2"
	FI_species = "D" # Example deuterium: "D"
	ω_init = 250.0 # Mode frequency. kHz
	ω_rot_init = 2.0 # Plasma rotation frequency. kHz
	n_init = 1 # Toroidal mode number
	m_init = 1 # Poloidal mode number
	verbose = true
	
	# You can specify additional entries in the dictionary below, to include more keyword arguments in 
	# the guiding-centre orbit integration algorithm
	extra_kw_args = Dict(:limit_phi => true, :maxiter => 0, :toa => true)
	# limits the number of toroidal turns for orbits
	# The orbit integration algorithm will try progressively smaller timesteps these number of times
end

# ╔═╡ c30b75f7-34e1-4b60-9c52-2efe4d57f803
begin
	# Load the OWCF environment and all necessary Julia packages
	cd(folderpath_OWCF)
	using Pkg
	Pkg.activate(".")
	
	verbose && println("Loading packages... ")
	using EFIT
	using Equilibrium
	using JLD2
	using GuidingCenterOrbits
	using Plots
	using PlutoUI
	using FileIO
	include(folderpath_OWCF*"extra/dependencies.jl")
end

# ╔═╡ 753914e7-52f4-4df3-a97a-ca2c179e5b8d
md"""
### Warning! Please note! 
For topological maps containing more than approximately 150 000 valid orbits (e.g. 20x100x100),
you should NOT use this web app (or any other OWCF web apps). As of the current OWCF version, the interactive web interface simply becomes too slow. Please do instead plot the energy slices manually. This can be done as shown in the first Julia cell below this cell (if you correctly specify the folderpath\_OWCF, filepath\_tm and E variables in the cell below, the cell will be able to run and plot):
"""

# ╔═╡ e4c79031-c111-4cd2-8d8f-f763b6c2bed1
let
	n = 1 # Toroidal mode number
	m = 1 # Poloidal mode number
	ω = 250 # kHz Wave frequency in kHz
	ω_rot = 2 # kHz Plasma rotation frequency in kHz
	E = 5 # keV Energy in keV

	
	myfile = jldopen(filepath_tm, false, false, false, IOStream)
	topoMap = myfile["topoMap"]
	E_array = myfile["E_array"]
	pm_array = myfile["pm_array"]
	Rm_array = myfile["Rm_array"]
	torTransTimes = myfile["torTransTimes"]
	polTransTimes = myfile["polTransTimes"]
	close(myfile)

	bad_coords = findall(x-> (x==9.0 || x==6.0 || x==7.0), topoMap) # Invalid, incomplete or lost orbits
    
    # Taking care of bad coordinates for the poloidal and toroidal transit times. NaNs will lead to those pixels not being plotted
    polTransTimes[bad_coords] .= NaN
    torTransTimes[bad_coords] .= NaN
    
    ω_θ = map(x-> !isnan(x) ? 0.001*2*pi/x : x, polTransTimes) # kHz
    ω_ϕ = map(x-> !isnan(x) ? 0.001*2*pi/x : x, torTransTimes) # kHz
	
	iE = argmin(abs.(E_array .- E)) # Find the closest value to E in E_array
	E = E_array[iE]
	npm = length(pm_array)
	nRm = length(Rm_array)
	
	resonance = -log10.(abs.(
		ω .*ones(npm,nRm) 
		- n .*ω_ϕ[iE,:,:] 
		- m .*ω_θ[iE,:,:]
		- ω_rot .*ones(npm,nRm) 
							))
	
	Plots.heatmap(Rm_array,pm_array,resonance,xlabel="Rm [m]", ylabel="pm", title="E: $(round(E,digits=3)) keV", fillcolor=cgrad([:white, :darkblue, :green, :yellow, :orange, :red]))
end

# ╔═╡ 84693850-008b-433d-b2d8-4fc493167189
begin
	# Process input data and pre-compute necessary quantities for web application

	# Loading tokamak equilibrium
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
    
    # Read the .jld2-file for displaying the topology of orbit space
    verbose && println("Loading topological map... ")
    myfile = jldopen(filepath_tm,false,false,false,IOStream)
    topoMap = myfile["topoMap"]
    E_array = myfile["E_array"]
    pm_array = myfile["pm_array"]
    Rm_array = myfile["Rm_array"]
    if haskey(myfile,"polTransTimes")
        polTransTimes = myfile["polTransTimes"]
        torTransTimes = myfile["torTransTimes"]
    else
        error("Toroidal and poloidal transit time data required for modeAnalysisWebApp.jl. Please re-run calcTopoMap.jl with saveTransitTimeMaps set to true, and use output file as filepath_tm input to modeAnalysisWebApp.jl.")
    end
    close(myfile)
    
    # Mapping topological map to (E,Λ,Pϕ_n;σ)
    if enable_COM && !(isfile(filepath_tm_COM))
        verbose && println(">>>>>> Mapping topological map from (E,pm,Rm) to (E,Λ,Pϕ_n;σ) <<<<<<... ")
        topoMap_COM, E_array, Λ_array, Pϕ_n_array = os2COM(M, topoMap, Vector(E_array), pm_array, Rm_array, FI_species; nl=2*length(pm_array), npp=2*length(Rm_array), isTopoMap=true, verbose=verbose)
    elseif enable_COM && isfile(filepath_tm_COM)
        verbose && println("Loading topological map in (E,Λ,Pϕ_n;σ) coordinates from filepath_tm_COM... ")
        myfile = jldopen(filepath_tm_COM,false,false,false,IOStream)
        topoMap_COM = myfile["topoMap"]
        E_array_COM = myfile["E_array"]
        Λ_array = myfile["Lambda_array_topoMap"]
        Pϕ_n_array = myfile["Pphi_n_array_topoMap"]
        close(myfile)
        if !(E_array==E_array_COM)
            error("Energy grid points in (E,pm,Rm) do not match energy grid points in (E,Λ,Pϕ_n;σ). Please correct and re-try.")
        end
    else
        verbose && println("Switching (E,pm,Rm) -> (E,Λ,Pϕ_n;σ) will not be possible.")
    end
    
    # Mapping maps of the poloidal and toroidal transit times to (E,Λ,Pϕ_n;σ), if available
    if enable_COM && !(isfile(filepath_tm_COM))
        verbose && println(">>>>>> Mapping poloidal transit times from (E,pm,Rm) to (E,Λ,Pϕ_n;σ) <<<<<<... ")
        valid_orbit_indices = findall(x-> (x!=9.0) && (x!=7.0), topoMap) # 9 and 7 are the integers representing invalid and lost orbits in the calcTopoMap.jl script, respectively. We don't want them.
        polTransTimes_COM, E_array_pol, Λ_array_pol, Pϕ_n_array_pol = os2COM(M, valid_orbit_indices, polTransTimes, E_array, pm_array, Rm_array, FI_species; nl=2*length(pm_array), npp=2*length(Rm_array), verbose=verbose)
        verbose && println(">>>>>> Mapping toroidal transit times from (E,pm,Rm) to (E,Λ,Pϕ_n;σ) <<<<<<... ")
        torTransTimes_COM, E_array_tor, Λ_array_tor, Pϕ_n_array_tor = os2COM(M, valid_orbit_indices, torTransTimes, E_array, pm_array, Rm_array, FI_species; nl=2*length(pm_array), npp=2*length(Rm_array), verbose=verbose)
    elseif enable_COM && isfile(filepath_tm_COM)
        verbose && println("Loading maps of poloidal and toroidal transit times in (E,Λ,Pϕ_n;σ) coordinates from filepath_tm_COM... ")
        myfile = jldopen(filepath_tm_COM,false,false,false,IOStream)
        polTransTimes_COM = myfile["polTransTimes"]
        E_array_pol = myfile["E_array"] # Very silly
        Λ_array_pol = myfile["Lambda_array_polTransTimes"]
        Pϕ_n_array_pol = myfile["Pphi_n_array_polTransTimes"]
        torTransTimes_COM = myfile["torTransTimes"]
        E_array_tor = myfile["E_array"] # Even more silly
        Λ_array_tor = myfile["Lambda_array_torTransTimes"]
        Pϕ_n_array_tor = myfile["Pphi_n_array_torTransTimes"]
        close(myfile)
    else
    end
    
    # Safety-checking
    if enable_COM
        verbose && println("Checking that the (E,pm,Rm) -> (E,Λ,Pϕ_n;σ) mapping yielded roughly the same E-, Λ- and Pϕ_n-grid points for topoMap, polTransTimes and torTransTimes... ")
        if !isapprox(E_array, E_array_pol)
            error("Energy grid points of topological map were not (roughly) equal to those of the poloidal transit time map. This should be impossible. orbitsWebApp.jl needs to be corrected manually.")
        end
        if !isapprox(E_array, E_array_tor)
            error("Energy grid points of topological map were not (roughly) equal to those of the toroidal transit time map. This should be impossible. orbitsWebApp.jl needs to be corrected manually.")
        end
    end
    if enable_COM
        verbose && println("Checking that (E,Λ,Pϕ_n;σ) grid points coincide for τ_ϕ and τ_θ... ")
        if !isapprox(Λ_array_pol,Λ_array_tor)
            error("Λ grid points did not approximately coincide for τ_ϕ and τ_θ. Something has gone terribly wrong. Please contact henrikj@dtu.dk or anvalen@dtu.dk.")
        end
        if !isapprox(Pϕ_n_array_pol,Pϕ_n_array_tor)
            error("Pϕ_n grid points did not approximately coincide for τ_ϕ and τ_θ. Something has gone terribly wrong. Please contact henrikj@dtu.dk or anvalen@dtu.dk.")
        end
    end
    
    # Double-checking so nothing has NaNs in it
    if !(sum(isnan.(topoMap))==0)
        topoMap = map(x-> isnan(x) ? 6.0 : x, topoMap)
        @warn "Topological map in (E,pm,Rm) had NaNs in it! orbitsWebApp.jl automatically set them to 6.0 (incomplete)."
    end
    if enable_COM
        if !(sum(isnan.(topoMap_COM))==0)
            topoMap_COM = map(x-> isnan(x) ? 6.0 : x, topoMap_COM)
            @warn "Topological map in (E,Λ,Pϕ_n;σ) had NaNs in it! orbitsWebApp.jl automatically set them to 6.0 (incomplete)."
        end
    end
    
    # Inverting poloidal and toroidal transit time data, to enable computation of resonance
    # Then multiply with 2*pi to get ω_ϕ and ω_θ
    verbose && println("Inverting poloidal and toroidal transit time data, to compute ω_θ and ω_ϕ, respectively... ")
    
    bad_coords = findall(x-> (x==9.0 || x==6.0 || x==7.0), topoMap) # Invalid, incomplete or lost orbits
    enable_COM && (bad_coords_COM = findall(x-> (x==9.0 || x==6.0 || x==7.0), topoMap_COM))
    
    # Taking care of bad coordinates for the poloidal and toroidal transit times. NaNs will lead to those pixels not being plotted
    polTransTimes[bad_coords] .= NaN
    torTransTimes[bad_coords] .= NaN
    enable_COM && (polTransTimes_COM[bad_coords_COM] .= NaN)
    enable_COM && (torTransTimes_COM[bad_coords_COM] .= NaN)
    
    ω_θ = map(x-> !isnan(x) ? 0.001*2*pi/x : x, polTransTimes) # kHz
    ω_ϕ = map(x-> !isnan(x) ? 0.001*2*pi/x : x, torTransTimes) # kHz
    enable_COM && (ω_θ_COM = map(x-> !isnan(x) ? 0.001*2*pi/x : x, polTransTimes_COM)) # kHz
    enable_COM && (ω_ϕ_COM = map(x-> !isnan(x) ? 0.001*2*pi/x : x, torTransTimes_COM)) # kHz
    
    verbose && println("Computing flux function on 100x100 (R,z)-grid... ")
    flux_r = range(extrema(wall.r)...,length=100)
    flux_z = range(extrema(wall.z)...,length=100)
    inds = CartesianIndices((length(flux_r),length(flux_z)))
    psi_rz = [M(flux_r[ind[1]], flux_z[ind[2]]) for ind in inds]
    psi_mag, psi_bdry = psi_limits(M)
    
    # Compute necessary quantities for COM tool
    if enable_COM
        verbose && println("Computing necessary quantities for constants-of-motion option... ")
        B0 = norm(Equilibrium.Bfield(M,magnetic_axis(M)...)) # Tesla
        q = getSpeciesCharge(FI_species) # Coulomb
        if psi_bdry==0
            @warn "The magnetic flux at the last closed flux surface (LCFS) is found to be 0 for the magnetic equilibrium in $(filepath_equil). Pϕ_n=Pϕ/(q*|Ψ_w|) where Ψ_w=Ψ(mag. axis) is assumed instead of Ψ_w=Ψ(LCFS)."
            Ψ_w_norm = abs(psi_axis)
        else
            Ψ_w_norm = abs(psi_bdry)
        end
    end
    
    # Geometrical quantities
    R_hfs = minimum(wall.r) # R-coord of high-field side wall
    R_lfs = maximum(wall.r) # R-coord of low-field side wall
    phi = collect(0:1:359).*(2*pi/180.0) # Toroidal angle
    topview_R_hfs_x = (R_hfs).*cos.(phi)
    topview_R_hfs_y = (R_hfs).*sin.(phi)
    topview_R_lfs_x = (R_lfs).*cos.(phi)
    topview_R_lfs_y = (R_lfs).*sin.(phi)

    # Ensure compatible data types for grid points
    E_array = vec(collect(E_array)) # Ensure type Array{Float64,1} (vector)
    pm_array = vec(collect(pm_array)) # Ensure type Array{Float64,1} (vector)
    Rm_array = vec(collect(Rm_array)) # Ensure type Array{Float64,1} (vector)
end

# ╔═╡ f73b5e13-4b8d-4770-afec-51df471bb7c2
function plot_inputs()

	return PlutoUI.combine() do Child

		inputs=[
			md""" Coordinate space: $(
				Child("coordspace", Select([:OS => "(E,pm,Rm)", :COM => "(E,Λ,Pϕ;σ)"]))
			)""",
			md""" Show tokamak wall: $(
				Child("wall", Switch(default=true))
			)""",
			md""" Show coordinate: $(
				Child("showcoord", Switch(default=true))
			)""",
			md""" Mode frequency (kHz): $(
				Child("omega", NumberField(10:10:1000, default=ω_init))
			)""",
			md""" Plasma rotation (kHz): $(
				Child("prot", NumberField(0:1:100, default=round(ω_rot_init)))
			)""",
			md""" Toroidal mode number (-): $(
				Child("n", Scrubbable(0:1:40, default=n_init))
			)""",
			md""" Poloidal mode number (-): $(
				Child("m", Scrubbable(0:1:40, default=m_init))
			)""",
			md""" Energy (keV): $(
				Child("E", Slider(E_array, default=E_array[1]))
			)""",
			md""" Pitch maximum, pm (-): $(
				Child("pm", Slider(pm_array, default=pm_array[Int64(round(3*length(pm_array)/4))]))
			)""",
			md""" Radius maximum, Rm (m): $(
				Child("Rm", Slider(Rm_array, default=Rm_array[Int64(round(length(Rm_array)/2))]))
			)"""
		]

		md"""
		#### Plot controls
		$(inputs)
		"""
	end
end

# ╔═╡ e0515346-ae9c-4bb6-8475-91c9699ea8c9
@bind my_plot_inputs plot_inputs()

# ╔═╡ 066df4bc-ff34-482b-82ba-2a3bec2490f7
let
	E = my_plot_inputs.E 
	pm = my_plot_inputs.pm
	Rm = my_plot_inputs.Rm 
	ω = my_plot_inputs.omega 
	ω_rot = my_plot_inputs.prot
	n = my_plot_inputs.n 
	m = my_plot_inputs.m 
	phase_space = my_plot_inputs.coordspace 
	tokamak_wall = my_plot_inputs.wall 
	show_coordinate = my_plot_inputs.showcoord

	EPRc = EPRCoordinate(M, E, pm, Rm; amu=getSpeciesAmu(FI_species), q=getSpeciesEcu(FI_species))
    o = get_orbit(M,EPRc; wall=wall, extra_kw_args...)
    if (phase_space==:COM) && enable_COM
        myHc = HamiltonianCoordinate(M, EPRc)
        μ = myHc.mu; Λ = μ*B0/E_joule
        Pϕ = myHc.p_phi; Pϕ_n = Pϕ/(q*Ψ_w_norm)
    end
    
    topview_o_x = cos.(o.path.phi).*(o.path.r)
    topview_o_y = sin.(o.path.phi).*(o.path.r)

    orb_color = :black
    orb_linestyle = :solid
    if o.class==:invalid
        orb_color = :gray
        orb_linestyle = :dash
    elseif o.class == :lost
        orb_color = :brown
    elseif o.class == :incomplete # If this happens, you are scr*w*d. Because it will likely take forever to calculate.
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
    plt_top = Plots.plot(topview_o_x,topview_o_y,label="$(o.class) orbit", color=orb_color, linestyle=orb_linestyle, linewidth=1.5)
    plt_top = Plots.plot!(topview_R_lfs_x,topview_R_lfs_y, label="Tokamak wall", color=:black, linewidth=1.5, xlabel="x [m]", ylabel="y [m]")
    plt_top = Plots.plot!(topview_R_hfs_x,topview_R_hfs_y, label="", color=:black,linewidth=1.5, aspect_ratio=:equal, title="Top view")

    #cross-sectional plot
    plt_crs = Plots.plot()
    if tokamak_wall
        wall_dR = maximum(wall.r)-minimum(wall.r)
        plt_crs = Plots.contour!(flux_r, flux_z, psi_rz', levels=collect(range(psi_mag, stop=psi_bdry,length=5)), color=:gray, linestyle=:dot,linewidth=1.5, label="",colorbar=false)
        plt_crs = Plots.plot!(wall.r, wall.z, label="Tokamak wall", color=:black, linewidth=1.5,xaxis=[minimum(wall.r)-wall_dR/10,maximum(wall.r)+wall_dR])
    end
    plt_crs = Plots.plot!(o.path.r,o.path.z, label="$(o.class) orbit", color=orb_color, linestyle=orb_linestyle, linewidth=1.5)
    if (phase_space==:OS) || !enable_COM
        plt_crs = Plots.plot!(title="E: $(round(E,digits=2)) keV  pm: $(round(o.coordinate.pitch, digits=2))  Rm: $(round(o.coordinate.r,digits=2))")
    else # phase_space==:COM (OS = orbit space, COM = constants-of-motion)
        plt_crs = Plots.plot!(title="E: $(round(E,digits=2)) keV  Λ: $(round(Λ, sigdigits=2))  Pϕ_n: $(round(Pϕ_n,sigdigits=2))")
    end
    plt_crs = Plots.scatter!([magnetic_axis(M)[1]],[magnetic_axis(M)[2]],label="Magnetic axis", mc=:gray, aspect_ratio=:equal, xlabel="R [m]", ylabel=" z [m]")
    plt_crs = Plots.scatter!([o.coordinate.r],[o.coordinate.z], mc=orb_color, label="(Rm,zm)")

    #topological plot(s)
    Eci = argmin(abs.(collect(E_array) .- E))
    if (phase_space==:OS) || !enable_COM
        pm_array_ext = vcat(pm_array[1]-(diff(pm_array))[1],pm_array) # Extend pm_array one row below
        topoMap_ext = vcat(vcat([1,2,3,4,5,6,7,8,9],ones(length(Rm_array)-9))', topoMap[Int64(Eci),:,:]) # Extend topoMap one row below, to ensure correct colormapping for orbit types (please see calcTopoMap.jl for more info). The y-limits (ylims) will make sure the extra row is not visible in the plot
        plt_topo = Plots.heatmap(Rm_array,pm_array_ext,topoMap_ext,color=:Set1_9,legend=false,xlabel="Rm [m]", ylabel="pm", title="E: $(round(E,digits=3)) keV", ylims=extrema(pm_array), xlims=extrema(Rm_array))
        if tokamak_wall && (maximum(Rm_array) >= maximum(wall.r))
            plt_topo = Plots.plot!(maximum(wall.r).*ones(length(pm_array)), pm_array, color=:black, linewidth=2)
        end
    else
        Λ_array_ext = vcat(Λ_array[1]-(diff(Λ_array))[1],Λ_array) # Extend Λ_array one row below
        if pm<0.0
            topoMap_COM_ext = vcat(vcat([1,2,3,4,5,6,7,8,9],ones(length(Pϕ_n_array)-9))', topoMap_COM[Int64(Eci),:,:,1]) # Extend topoMap_COM one row below, to ensure correct colormapping for orbit types (please see calcTopoMap.jl for more info). The y-limits (ylims) will make sure the extra row is not visible in the plot
            plt_topo = Plots.heatmap(Pϕ_n_array,Λ_array_ext,topoMap_COM_ext,color=:Set1_9,legend=false,xlabel="Pϕ_n", ylabel="Λ", title="E: $(round(E_array[Int64(Eci)],digits=3)) keV", ylims=extrema(Λ_array), xlims=extrema(Pϕ_n_array))
        else
            topoMap_COM_ext = vcat(vcat([1,2,3,4,5,6,7,8,9],ones(length(Pϕ_n_array)-9))', topoMap_COM[Int64(Eci),:,:,2]) # Extend topoMap_COM one row below, to ensure correct colormapping for orbit types (please see calcTopoMap.jl for more info). The y-limits (ylims) will make sure the extra row is not visible in the plot
            plt_topo = Plots.heatmap(Pϕ_n_array,Λ_array_ext,topoMap_COM_ext,color=:Set1_9,legend=false,xlabel="Pϕ_n", ylabel="Λ", title="E: $(round(E_array[Int64(Eci)],digits=3)) keV", ylims=extrema(Λ_array), xlims=extrema(Pϕ_n_array))
        end
    end
    if show_coordinate
        if (phase_space==:OS) || !enable_COM
            plt_topo = Plots.scatter!([Rm],[pm],markershape=:circle,mc=orb_color,legend=false,markersize=6)
        else
            plt_topo = Plots.scatter!([Pϕ_n],[Λ],markershape=:circle,mc=orb_color,legend=false,markersize=6)
        end
    end

    # Resonance plot
    if (phase_space==:OS) || !enable_COM
        npm = length(pm_array)
        nRm = length(Rm_array)
        resonance = -log10.(abs.(
			ω .*ones(npm,nRm) 
			- n .*ω_ϕ[Int64(Eci),:,:] 
			- m .*ω_θ[Int64(Eci),:,:]
			- ω_rot .*ones(npm,nRm)
		))
        plt_res = Plots.heatmap(Rm_array,pm_array,resonance,xlabel="Rm [m]", ylabel="pm", title="E: $(round(E,digits=3)) keV", fillcolor=cgrad([:white, :darkblue, :green, :yellow, :orange, :red]), ylims=extrema(pm_array), xlims=extrema(Rm_array))
    else
        nl = length(Λ_array_tor) # tor or pol should not matter
        npp = length(Pϕ_n_array_pol) # tor or pol should not matter
        if pm<0.0
            resonance = -log10.(abs.(
				ω .*ones(nl,npp) 
				- n .*ω_ϕ_COM[Int64(Eci),:,:,1] 
				- m .*ω_θ_COM[Int64(Eci),:,:,1] 
				- ω_rot .*ones(nl,npp)
			))
        else
            resonance = -log10.(abs.(
				ω .*ones(nl,npp) 
				- n .*ω_ϕ_COM[Int64(Eci),:,:,2] 
				- m .*ω_θ_COM[Int64(Eci),:,:,2] 
				- ω_rot .*ones(nl,npp)
			))
        end

        plt_res = Plots.heatmap(Pϕ_n_array_tor[Int64(Eci),:],Λ_array_tor[Int64(Eci),:],resonance,legend=false,xlabel="Pϕ_n", ylabel="Λ", title="E: $(round(E,digits=3)) keV", fillcolor=cgrad([:white, :darkblue, :green, :yellow, :orange, :red]), ylims=extrema(Λ_array_tor), xlims=extrema(Pϕ_n_array_tor))
    end
    if show_coordinate
        if (phase_space==:OS) || !enable_COM
            plt_res = Plots.scatter!([Rm],[pm],markershape=:circle,mc=orb_color,legend=false,markersize=6)
        else
            plt_res = Plots.scatter!([Pϕ_n],[Λ],markershape=:circle,mc=orb_color,legend=false,markersize=6)
        end
    end

	Plots.plot(plt_res, plt_topo, plt_top, plt_crs, layout=(2,2), size=(800, 800))
end

# ╔═╡ 8572b2a4-34c8-4fe5-b106-759d78d20d33


# ╔═╡ f0823b81-2340-4495-8f1b-544f5e0ec9e3


# ╔═╡ 26c118fb-774d-4c2f-a2ee-2c2f202e3e56


# ╔═╡ 115b02e7-3660-4d20-9c18-11fb6e4b7718


# ╔═╡ 714de67a-7f48-4c7d-90c9-5759aeb45258


# ╔═╡ 09d757ea-9472-48a9-b22d-43876492a889


# ╔═╡ 3ed81523-c733-492e-9bca-19771db89d7a


# ╔═╡ 35b0b0bc-98a5-4adb-a29b-2d7fd3bf55b3


# ╔═╡ 52bba0f0-e584-4e5d-aca4-8f75d22cfd97


# ╔═╡ Cell order:
# ╠═e1984f36-a6a7-11f1-16fc-a9fcbbb60876
# ╠═4e590375-f4fd-47ca-b89c-fe46dbdc22c2
# ╠═c30b75f7-34e1-4b60-9c52-2efe4d57f803
# ╠═753914e7-52f4-4df3-a97a-ca2c179e5b8d
# ╠═e4c79031-c111-4cd2-8d8f-f763b6c2bed1
# ╠═84693850-008b-433d-b2d8-4fc493167189
# ╠═f73b5e13-4b8d-4770-afec-51df471bb7c2
# ╠═e0515346-ae9c-4bb6-8475-91c9699ea8c9
# ╠═066df4bc-ff34-482b-82ba-2a3bec2490f7
# ╠═8572b2a4-34c8-4fe5-b106-759d78d20d33
# ╠═f0823b81-2340-4495-8f1b-544f5e0ec9e3
# ╠═26c118fb-774d-4c2f-a2ee-2c2f202e3e56
# ╠═115b02e7-3660-4d20-9c18-11fb6e4b7718
# ╠═714de67a-7f48-4c7d-90c9-5759aeb45258
# ╠═09d757ea-9472-48a9-b22d-43876492a889
# ╠═3ed81523-c733-492e-9bca-19771db89d7a
# ╠═35b0b0bc-98a5-4adb-a29b-2d7fd3bf55b3
# ╠═52bba0f0-e584-4e5d-aca4-8f75d22cfd97
