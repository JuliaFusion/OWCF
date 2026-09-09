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

# ╔═╡ a761bf9e-a159-11f1-2ec7-711778cfac0d
md"""
# orbitsWebApp

## Description:
This notebook provides an application to visualize (guiding-centre) orbits in a tokamak in an interactive and intuitive manner. Also, the topological regions of orbit space (E,pm,Rm) are visualized. In addition, if computed, maps of the poloidal and toroidal transit times are visualized as well. Prior to running this notebook, please make sure you have run the OWCF/calcTopoMap.jl tool (or equivalent). Also, please make sure that you have noted the path to the saved topological map output file. More specifically, If the filepath\_tm .jld2-file (see below) has keys named 'polTransTimes' and 'torTransTimes', the user will also be able to visualize maps of the poloidal and toroidal transit times. In the app below, this is done automatically.

## Inputs:
- enable\_COM - If true, (E,pm,Rm) -> (E,Λ,Pϕ_n;σ) can be performed via a toggle button. Set to false to minimize computation time. THIS IS TEMPORARILY OUT-OF-ORDER - Bool
- folderpath\_OWCF - The path to the OWCF folder on your computer. Needed for correct loading - String
- filepath\_equil - The path to the file with the tokamak magnetic equilibrium and geometry - String
- filepath\_tm - The path to the .jld2-file containing the topological map (and more). Can be computed using the OWCF/calcTopoMap.jl script - String
- FI\_species - The species of the particle being simulated (deuterium, tritium etc). Specified as "D", "T" etc. Please see the keys of the OWCF\_ChemElemSymbol\_to\_Z Dict object in the file OWCF/misc/species_func.jl for a full list of all available particle species in the OWCF - String
- verbose - If set to true, the app will talk a lot! - Bool

## Outputs:
# -

## Saved files:
# -

### Notebook written by Henrik Järleblad, henrikj@dtu.dk
### Last maintained 2026-09-09
"""

# ╔═╡ cd6a92a9-7113-4281-8596-a520642c5cbe
begin
    # SPECIFY THE INPUTS IN THIS CELL

    # Please specify the OWCF folder and let the notebook change directory to the 
    # OWCF folder when the cell below is executed. This is to be able to load the
    # correct versions of the Julia packages as specified in the Project.toml and 
    # Manifest.toml files.
    folderpath_OWCF = "REPLACE-THIS-TEXT-WITH-THE-PATH-TO-THE-OWCF-FOLDER-ON-YOUR-COMPUTER" # Finish with '/'
    
    enable_COM = false # TEMPORARILY OUT-OF-ORDER!!! DO NOT SET TO 'true'!!! Set to false for large grids (>10x100x100 in (E,pm,Rm)). The computation time simply becomes too large. Or...
    if enable_COM
        filepath_tm_COM = "" # ...please specify an output of the os2com.jl script (which contains the key "topoMap")(and possibly "polTransTimes" and "torTransTimes"). This can also have been produced automatically from orbitsWebApp.jl. Leave unspecified if orbitsWebApp.jl should compute the (E,pm,Rm) -> (E,Λ,Pϕ_n;σ) map itself.
    end
    filepath_equil = folderpath_OWCF*"equilibrium/JET/g99971/g99971_474-48.9.eqdsk"  # Example JET shot 96100 at 13s (53 minus 40): g96100/g96100_0-53.0012.eqdsk" #
    filepath_tm = folderpath_OWCF*"apps/example_data/topoMap_JET_99971L72_at48,9s_D_6x101x102_wLost.jld2"
    FI_species = "D" # Example deuterium: "D"
    verbose = true
    
    # You can specify additional entries in the dictionary below, to include more keyword arguments in the guiding-centre orbit integration algorithm:
    extra_kw_args = Dict(:limit_phi => true, :maxiter => 0, :toa => true)
    # 'limit_phi' limits the number of toroidal turns for orbits. 
    # 'maxiter' will let the orbit integration algorithm try progressively smaller timesteps these number of times
    # 'toa' stands for "try only adaptive" orbit trajectory integration algorithm
end

# ╔═╡ 29d500df-9ab5-449c-b1a2-a081e3cd770f
begin
	# Load OWCF environment and Julia packages
	cd(folderpath_OWCF)
	using Pkg
	Pkg.activate(".")
	
	verbose && println("Loading packages... ")
	using EFIT
	using Equilibrium
	using GuidingCenterOrbits
	using JLD2
	using LinearAlgebra
	using Plots
	using PlutoUI
	using FileIO
	using Dates
	using Statistics
	include(folderpath_OWCF*"extra/dependencies.jl")
	println("")
end

# ╔═╡ 23772204-992e-4bd3-ad1f-f1d86c33f6aa
md"""
### Warning! Please note! 
For topological maps containing more than approximately 150 000 valid orbits (e.g. 20x100x100),
you should NOT use this web app (or any other OWCF web apps). As of the current OWCF version, the interactive web interface simply becomes too slow. Please do instead plot the energy slices manually. This can be done as shown in the first Julia cell below this cell (if you correctly specify the folderpath\_OWCF, filepath\_tm and E variables in the cell below, the cell will be able to run and plot):
"""

# ╔═╡ dfacf685-d70b-4433-9d7c-e24fac4cf794
let
	E = 150.0 # keV Energy of the (pm,Rm) slice of interest
	
	myfile = jldopen(filepath_tm, false, false, false, IOStream)
	topoMap = myfile["topoMap"]
	E_array = myfile["E_array"]
	pm_array = myfile["pm_array"]
	Rm_array = myfile["Rm_array"]
	close(myfile)
	
	iE = argmin(abs.(E_array .- E)) # Find the closest value to E in E_array
	Plots.heatmap(Rm_array,pm_array,topoMap[iE,:,:],color=:Set1_9,legend=false,xlabel="Rm [m]", ylabel="pm", title="E: $(round(E,digits=3)) keV")
end

# ╔═╡ 206460f9-93e3-48cc-bd65-fbf0ce0a9967
begin
	## Read and load data
    
    # Loading tokamak equilibrium
    verbose && println("Loading magnetic equilibrium... ")
    M, wall, jdotb = nothing, nothing, nothing # Initialize global magnetic equilibrium variables
    try
        global M; global wall; global jdotb # Declare global scope
        M, wall = read_geqdsk(filepath_equil,clockwise_phi=false) # Assume counter-clockwise phi-direction
        jdotb = M.sigma # The sign of the dot product between the plasma current and the magnetic field
    catch # Otherwise, assume magnetic equilibrium is a saved .jld2 file
        global M; global wall; global jdotb; local myfile # Declare global scope 
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
    poltor = false
    if haskey(myfile,"polTransTimes")
        println("-- Maps of poloidal and toroidal transit times found! Including... ")
        polTransTimes = myfile["polTransTimes"]
        torTransTimes = myfile["torTransTimes"]
        poltor = true
    end
    close(myfile)

    # Map from (E.pm,Rm) to (E,Λ,Pϕ_n;σ), if required
	# Mapping topological map to (E,Λ,Pϕ_n;σ)
    if enable_COM && !(isfile(filepath_tm_COM))
        verbose && println(">>>>>> Mapping topological map from (E,pm,Rm) to (E,Λ,Pϕ_n;σ) <<<<<<... ")
        topoMap_COM, E_dummy, Λ_array, Pϕ_n_array = os2COM(M, topoMap, Vector(E_array), Vector(pm_array), Vector(Rm_array), FI_species; nl=2*length(pm_array), npp=2*length(Rm_array), isTopoMap=true, verbose=verbose)
    elseif enable_COM && isfile(filepath_tm_COM)
        verbose && println("Loading topological map in (E,Λ,Pϕ_n;σ) coordinates from filepath_tm_COM... ")
        myfile_tm_COM = jldopen(filepath_tm_COM,false,false,false,IOStream)
        topoMap_COM = myfile_tm_COM["topoMap"]
        E_array_COM = myfile_tm_COM["E_array"]
        Λ_array = myfile_tm_COM["Lambda_array_topoMap"]
        Pϕ_n_array = myfile_tm_COM["Pphi_n_array_topoMap"]
        close(myfile_tm_COM)
        if !(E_array==E_array_COM)
            error("Energy grid points in (E,pm,Rm) do not match energy grid points in (E,Λ,Pϕ_n;σ). Please correct and re-try.")
        end
    else
        verbose && println("Switching (E,pm,Rm) -> (E,Λ,Pϕ_n;σ) will not be possible.")
    end
    
    # Mapping maps of the poloidal and toroidal transit times to (E,Λ,Pϕ_n;σ), if available
    if poltor
        if enable_COM && !(isfile(filepath_tm_COM))
            verbose && println(">>>>>> Mapping poloidal transit times from (E,pm,Rm) to (E,Λ,Pϕ_n;σ) <<<<<<... ")
            valid_orbit_indices = findall(x-> (x!=9.0) && (x!=7.0), topoMap) # 9 and 7 are the integers representing invalid and lost orbits in the calcTopoMap.jl script, respectively. We don't want them.
            polTransTimes_COM, E_array_pol, Λ_array_pol, Pϕ_n_array_pol = os2COM(M, valid_orbit_indices, polTransTimes, E_array, pm_array, Rm_array, FI_species; nl=2*length(pm_array), npp=2*length(Rm_array), verbose=verbose)
            verbose && println(">>>>>> Mapping toroidal transit times from (E,pm,Rm) to (E,Λ,Pϕ_n;σ) <<<<<<... ")
            torTransTimes_COM, E_array_tor, Λ_array_tor, Pϕ_n_array_tor = os2COM(M, valid_orbit_indices, torTransTimes, E_array, pm_array, Rm_array, FI_species; nl=2*length(pm_array), npp=2*length(Rm_array), verbose=verbose)
        elseif enable_COM && isfile(filepath_tm_COM)
            verbose && println("Loading maps of poloidal and toroidal transit times in (E,Λ,Pϕ_n;σ) coordinates from filepath_tm_COM... ")
            myfile_tm_COM = jldopen(filepath_tm_COM,false,false,false,IOStream)
            polTransTimes_COM = myfile_tm_COM["polTransTimes"]
            E_array_pol = myfile_tm_COM["E_array"] # Very silly
            Λ_array_pol = myfile_tm_COM["Lambda_array_polTransTimes"]
            Pϕ_n_array_pol = myfile_tm_COM["Pphi_n_array_polTransTimes"]
            torTransTimes_COM = myfile_tm_COM["torTransTimes"]
            E_array_tor = myfile_tm_COM["E_array"] # Even more silly
            Λ_array_tor = myfile_tm_COM["Lambda_array_torTransTimes"]
            Pϕ_n_array_tor = myfile_tm_COM["Pphi_n_array_torTransTimes"]
            close(myfile)
        else
        end
    end
    
    if enable_COM && !(isfile(filepath_tm_COM))
        verbose && println("Saving topological map in (E,Λ,Pϕ_n;σ) format... ")
        nmu = length(Λ_array)
        nPphi = length(Pϕ_n_array)
        date_and_time = split("$(Dates.now())","T")[1]*"at"*split("$(Dates.now())","T")[2][1:5]
        filepath_output_orig = folderpath_OWCF*"orbitsWebApp_COM_data_$(length(E_array))x$(nmu)x$(nPphi)x2_"*date_and_time
        global filepath_output = deepcopy(filepath_output_orig)
        global C = 1
        while isfile(filepath_output*".jld2") # To take care of not overwriting files. Add _(1), _(2) etc
            global filepath_output
            global C
            filepath_output = filepath_output_orig*"_($(Int64(C)))"
            C += 1 # global scope, to surpress warnings
        end
        filepath_output = filepath_output*".jld2"
        myfile_COMsave = jldopen(filepath_output,true,true,false,IOStream)
        write(myfile_COMsave,"topoMap",topoMap_COM)
        write(myfile_COMsave,"E_array",E_array)
        write(myfile_COMsave,"Lambda_array_topoMap", Λ_array)
        write(myfile_COMsave,"Pphi_n_array_topoMap", Pϕ_n_array)
        if poltor
            verbose && println("Saving poloidal and toroidal transit times in (E,Λ,Pϕ_n;σ) format... ")
            write(myfile_COMsave,"polTransTimes",polTransTimes_COM)
            write(myfile_COMsave,"Lambda_array_polTransTimes",Λ_array_pol)
            write(myfile_COMsave,"Pphi_n_array_polTransTimes",Pϕ_n_array_pol)
            write(myfile_COMsave,"torTransTimes",torTransTimes_COM)
            write(myfile_COMsave,"Lambda_array_torTransTimes",Λ_array_tor)
            write(myfile_COMsave,"Pphi_n_array_torTransTimes",Pϕ_n_array_tor)
        end
        close(myfile_COMsave)
        verbose && println("----> NEXT TIME orbitsWebApp.jl IS RUN WITH THE SAME INPUTS, set the 'filepath_tm_COM' input variable to "*filepath_output*", to avoid having to re-do the (E,pm,Rm) -> (E,Λ,Pϕ_n;σ) mapping!")
        verbose && println("----> PLEASE DELETE "*filepath_output*" MANUALLY IF NOT NEEDED.")
    end

    # Safety checks
	## Safety-checking energy grid points
    if poltor && enable_COM
        verbose && println("Checking that the (E,pm,Rm) -> (E,Λ,Pϕ_n;σ) mapping yielded roughly the same Λ- and Pϕ_n-grid points for topoMap, polTransTimes and torTransTimes... ")
        if !isapprox(E_array, E_array_pol)
            error("Energy grid points of topological map were not (roughly) equal to those of the poloidal transit time map. This should be impossible. orbitsWebApp.jl needs to be corrected manually.")
        end
        if !isapprox(E_array, E_array_tor)
            error("Energy grid points of topological map were not (roughly) equal to those of the toroidal transit time map. This should be impossible. orbitsWebApp.jl needs to be corrected manually.")
        end
    end
    
    ## Double-checking so nothing has NaNs in it
    if !(sum(isnan.(topoMap))==0)
        topoMap = map(x-> isnan(x) ? 0.0 : x, topoMap)
        @warn "Topological map in (E,pm,Rm) had NaNs in it! orbitsWebApp.jl automatically set them to 0.0."
    end
    if enable_COM
        if !(sum(isnan.(topoMap_COM))==0)
            topoMap_COM = map(x-> isnan(x) ? 0.0 : x, topoMap_COM)
            @warn "Topological map in (E,Λ,Pϕ_n;σ) had NaNs in it! orbitsWebApp.jl automatically set them to 0.0."
        end
    end
    if poltor
        if !(sum(isnan.(polTransTimes))==0)
            polTransTimes = map(x-> isnan(x) ? 0.0 : x, polTransTimes)
            @warn "Map of poloidal transit times in (E,pm,Rm) had NaNs in it! orbitsWebApp.jl automatically set them to 0.0."
        end
        if !(sum(isnan.(torTransTimes))==0)
            torTransTimes = map(x-> isnan(x) ? 0.0 : x, torTransTimes)
            @warn "Map of toroidal transit times in (E,pm,Rm) had NaNs in it! orbitsWebApp.jl automatically set them to 0.0."
        end
        if enable_COM
            if !(sum(isnan.(polTransTimes_COM))==0)
                polTransTimes_COM = map(x-> isnan(x) ? 0.0 : x, polTransTimes_COM)
                @warn "Map of poloidal transit times in (E,Λ,Pϕ_n;σ) had NaNs in it! orbitsWebApp.jl automatically set them to 0.0."
            end
            if !(sum(isnan.(torTransTimes_COM))==0)
                torTransTimes_COM = map(x-> isnan(x) ? 0.0 : x, torTransTimes_COM)
                @warn "Map of toroidal transit times in (E,Λ,Pϕ_n;σ) had NaNs in it! orbitsWebApp.jl automatically set them to 0.0."
            end
        end
    end

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
    
    # Necessary quantities for plotting the web application
    R_hfs = minimum(wall.r) # R-coord of high-field side wall
    R_lfs = maximum(wall.r) # R-coord of low-field side wall
    phi = collect(0:1:359).*(2*pi/180.0) # Toroidal angle
    topview_R_hfs_x = (R_hfs).*cos.(phi)
    topview_R_hfs_y = (R_hfs).*sin.(phi)
    topview_R_lfs_x = (R_lfs).*cos.(phi)
    topview_R_lfs_y = (R_lfs).*sin.(phi)
    
    E_array = vec(collect(E_array)) # Ensure type Array{Float64,1} (vector)
    pm_array = vec(collect(pm_array)) # Ensure type Array{Float64,1} (vector)
    Rm_array = vec(collect(Rm_array)) # Ensure type Array{Float64,1} (vector)
end

# ╔═╡ e9cd967d-c461-4324-99e0-6110bf7d43fc
function plot_inputs()

	return PlutoUI.combine() do Child

		inputs=[
			md""" Show tokamak wall: $(
				Child("wall", Switch(default=true))
			)""",
			md""" Coordinate space: $(
				Child("coordspace", Select([:OS => "(E,pm,Rm)", :COM => "(E,Λ,Pϕ;σ)"]))
			)""",
			md""" Show coordinate: $(
				Child("showcoord", Switch(default=true))
			)""",
			md""" Energy (keV): $(
				Child("E", Slider(E_array, default=E_array[Int64(round(length(E_array)/2))]))
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

# ╔═╡ a37e90da-5be0-4b14-8613-1e146dca5ec0
@bind my_plot_inputs plot_inputs()

# ╔═╡ 684923b0-dcf4-4ff3-9d50-8291ddb0b593
let
    phase_space = my_plot_inputs.coordspace
    tokamak_wall = my_plot_inputs.wall
    show_coordinate = my_plot_inputs.showcoord
    E = my_plot_inputs.E
    pm = my_plot_inputs.pm
    Rm = my_plot_inputs.Rm
    
	EPRc = EPRCoordinate(M, E, pm, Rm; amu=getSpeciesAmu(FI_species), q=getSpeciesEcu(FI_species))
    o = get_orbit(M,EPRc; wall=wall, extra_kw_args...)
    if (phase_space==:COM) && enable_COM
        E_joule = (1000*GuidingCenterOrbits.e0)*E
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
    elseif o.class == :incomplete # If this happens, you are scr*w*ed. Because it will likely take forever to calculate.
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
        plt_crs = Plots.contour!(flux_r,flux_z,psi_rz',levels=collect(range(psi_mag,stop=psi_bdry,length=5)),color=:gray, linestyle=:dot,linewidth=1.5, label="",colorbar=false)
        plt_crs = Plots.plot!(wall.r,wall.z, label="Tokamak wall", color=:black, linewidth=1.5,xaxis=[minimum(wall.r)-wall_dR/10,maximum(wall.r)+wall_dR])
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

    # Poloidal and toroidal transit time maps
    if poltor
        if (phase_space==:OS) || !enable_COM
            pTT_microsecs = polTransTimes[Int64(Eci),:,:] ./(1.0e-6) # Convert from seconds to microseconds
            nz_coords = findall(x-> x>0.0,pTT_microsecs) # Find the 2D matrix coordinates of all non-zero elements
            my_coords = length(nz_coords) > 1 ? nz_coords : CartesianIndices(size(pTT_microsecs)) # Are there actually more than one non-zero element? If not, use all elements
            min_pol, max_pol = extrema(pTT_microsecs[my_coords]) # Find minimum and maximum values
            min_OOM, max_OOM = (floor(log10(min_pol)),ceil(log10(max_pol))) # The orders of magnitude of the minimum and maximum values
            if !((max_OOM-min_OOM)==0.0) && (length(nz_coords) > 1) && ((max_pol/min_pol) > 10) # If all values are NOT within same order of magnitude AND more than one non-zero element, use logarithmic colorbar
                plt_pol = Plots.heatmap(Rm_array, pm_array, pTT_microsecs, xlabel="Rm [m]", ylabel="pm", title="tau_pol(Rm,pm) [microseconds] \n tau_pol($(round(Rm,digits=2)),$(round(pm,digits=2)))=$(round(o.tau_p /(1.0e-6),sigdigits=3)) microseconds", fillcolor=cgrad([:white, :darkblue, :green, :yellow, :orange, :red]), colorbar=true, colorbar_scale=:log10, clims = (10^min_OOM, 10^max_OOM), top_margin=3Plots.mm) # Get nice powers-of-ten limits for the colorbar
            else
                plt_pol = Plots.heatmap(Rm_array, pm_array, pTT_microsecs, xlabel="Rm [m]", ylabel="pm", title="tau_pol(Rm,pm) [microseconds] \n tau_pol($(round(Rm,digits=2)),$(round(pm,digits=2)))=$(round(o.tau_p /(1.0e-6),sigdigits=3)) microseconds", fillcolor=cgrad([:white, :darkblue, :green, :yellow, :orange, :red]), colorbar=true, top_margin=3Plots.mm)
            end
            tTT_microsecs = torTransTimes[Int64(Eci),:,:] ./(1.0e-6) # Convert from seconds to microseconds
            nz_coords = findall(x-> x>0.0,tTT_microsecs) # Find the 2D matrix coordinates of all non-zero elements
            my_coords = length(nz_coords) > 1 ? nz_coords : CartesianIndices(size(tTT_microsecs)) # Are there actually more than one non-zero element? If not, use all elements
            min_tor, max_tor = extrema(tTT_microsecs[my_coords]) # Find minimum and maximum values
            min_OOM, max_OOM = (floor(log10(min_tor)),ceil(log10(max_tor))) # The orders of magnitude of the minimum and maximum values
            if !((max_OOM-min_OOM)==0.0) && (length(nz_coords) > 1) && ((max_tor/min_tor) > 10) # If all values are NOT within same order of magnitude AND more than one non-zero element, use logarithmic colorbar
                plt_tor = Plots.heatmap(Rm_array, pm_array, tTT_microsecs, xlabel="Rm [m]", ylabel="pm", title="tau_tor(Rm,pm) [microseconds] \n tau_tor($(round(Rm,digits=2)),$(round(pm,digits=2)))=$(round(o.tau_t /(1.0e-6),sigdigits=3)) microseconds", fillcolor=cgrad([:white, :darkblue, :green, :yellow, :orange, :red]), colorbar=true, colorbar_scale=:log10, clims = (10^min_OOM, 10^max_OOM), top_margin=3Plots.mm) # Get nice powers-of-ten limits for the colorbar
            else # Else, use linear colorbar
                plt_tor = Plots.heatmap(Rm_array, pm_array, tTT_microsecs, xlabel="Rm [m]", ylabel="pm", title="tau_tor(Rm,pm) [microseconds] \n tau_tor($(round(Rm,digits=2)),$(round(pm,digits=2)))=$(round(o.tau_t /(1.0e-6),sigdigits=3)) microseconds", fillcolor=cgrad([:white, :darkblue, :green, :yellow, :orange, :red]), colorbar=true, top_margin=3Plots.mm)
            end
            if show_coordinate
                plt_pol = Plots.scatter!(plt_pol, [Rm],[pm],markershape=:circle,mc=orb_color,legend=false,markersize=6)
                plt_tor = Plots.scatter!(plt_tor, [Rm],[pm],markershape=:circle,mc=orb_color,legend=false,markersize=6)
            end
        else
            if pm<0.0
                pTT_microsecs = polTransTimes_COM[Int64(Eci),:,:,1] ./(1.0e-6)
                tTT_microsecs = torTransTimes_COM[Int64(Eci),:,:,1] ./(1.0e-6)
            else
                pTT_microsecs = polTransTimes_COM[Int64(Eci),:,:,2] ./(1.0e-6)
                tTT_microsecs = torTransTimes_COM[Int64(Eci),:,:,2] ./(1.0e-6)
            end
            nz_coords = findall(x-> x>0.0,pTT_microsecs) # Find the 2D matrix coordinates of all non-zero elements
            my_coords = length(nz_coords) > 1 ? nz_coords : CartesianIndices(size(pTT_microsecs)) # Are there actually more than one non-zero element? If not, use all elements
            min_pol, max_pol = extrema(pTT_microsecs[my_coords]) # Find minimum and maximum values
            min_OOM, max_OOM = (floor(log10(min_pol)),ceil(log10(max_pol))) # The orders of magnitude of the minimum and maximum values
            if !((max_OOM-min_OOM)==0.0) && (length(nz_coords) > 1) && ((max_pol/min_pol) > 10) # If all values are NOT within same order of magnitude AND more than one non-zero element, use logarithmic colorbar
                plt_pol = Plots.heatmap(Pϕ_n_array_pol,Λ_array_pol, pTT_microsecs,xlabel="Pϕ_n", ylabel="Λ", title="tau_pol(Pϕ_n,Λ) [microseconds] \n tau_pol($(round(Pϕ_n,sigdigits=2)),$(round(Λ,sigdigits=2)))=$(round(o.tau_p /(1.0e-6),sigdigits=3)) microseconds", fillcolor=cgrad([:white, :darkblue, :green, :yellow, :orange, :red]), colorbar=true, colorbar_scale=:log10, clims = (10^min_OOM, 10^max_OOM), top_margin=3Plots.mm, ylims=extrema(Λ_array_pol), xlims=extrema(Pϕ_n_array_pol)) # Get nice powers-of-ten limits for the colorbar
            else
                plt_pol = Plots.heatmap(Pϕ_n_array_pol,Λ_array_pol, pTT_microsecs,xlabel="Pϕ_n", ylabel="Λ", title="tau_pol(Pϕ_n,Λ) [microseconds] \n tau_pol($(round(Pϕ_n,sigdigits=2)),$(round(Λ,sigdigits=2)))=$(round(o.tau_p /(1.0e-6),sigdigits=3)) microseconds", fillcolor=cgrad([:white, :darkblue, :green, :yellow, :orange, :red]), colorbar=true, top_margin=3Plots.mm, ylims=extrema(Λ_array_pol), xlims=extrema(Pϕ_n_array_pol))
            end

            nz_coords = findall(x-> x>0.0,tTT_microsecs) # Find the 2D matrix coordinates of all non-zero elements
            my_coords = length(nz_coords) > 1 ? nz_coords : CartesianIndices(size(tTT_microsecs)) # Are there actually more than one non-zero element? If not, use all elements
            min_tor, max_tor = extrema(tTT_microsecs[my_coords]) # Find minimum and maximum values
            min_OOM, max_OOM = (floor(log10(min_tor)),ceil(log10(max_tor))) # The orders of magnitude of the minimum and maximum values
            if !((max_OOM-min_OOM)==0.0) && (length(nz_coords) > 1) && ((max_tor/min_tor) > 10) # If all values are NOT within same order of magnitude AND more than one non-zero element, use logarithmic colorbar
                plt_tor = Plots.heatmap(Pϕ_n_array_tor,Λ_array_tor, tTT_microsecs,xlabel="Pϕ_n", ylabel="Λ", title="tau_tor(Pϕ_n,Λ) [microseconds] \n tau_tor($(round(Pϕ_n,sigdigits=2)),$(round(Λ,sigdigits=2)))=$(round(o.tau_t /(1.0e-6),sigdigits=3)) microseconds", fillcolor=cgrad([:white, :darkblue, :green, :yellow, :orange, :red]), colorbar=true, colorbar_scale=:log10, clims = (10^min_OOM, 10^max_OOM), top_margin=3Plots.mm, ylims=extrema(Λ_array_tor), xlims=extrema(Pϕ_n_array_tor)) # Get nice powers-of-ten limits for the colorbar
            else
                plt_tor = Plots.heatmap(Pϕ_n_array_tor,Λ_array_tor, tTT_microsecs,xlabel="Pϕ_n", ylabel="Λ", title="tau_tor(Pϕ_n,Λ) [microseconds] \n tau_tor($(round(Pϕ_n,sigdigits=2)),$(round(Λ,sigdigits=2)))=$(round(o.tau_t /(1.0e-6),sigdigits=3)) microseconds", fillcolor=cgrad([:white, :darkblue, :green, :yellow, :orange, :red]), colorbar=true, top_margin=3Plots.mm, ylims=extrema(Λ_array_tor), xlims=extrema(Pϕ_n_array_tor))
            end
            if show_coordinate
                plt_pol = Plots.scatter!(plt_pol, [Pϕ_n],[Λ],markershape=:circle,mc=orb_color,label="",markersize=6)
                plt_tor = Plots.scatter!(plt_tor, [Pϕ_n],[Λ],markershape=:circle,mc=orb_color,label="",markersize=6)
            end
        end
    end

    #pitch maximum visualization plot
    xmin = (pm >= 0.0) ? 0.0 : pm
    xp = pm
    yp = sqrt(1-pm^2)
    plt_pm = Plots.quiver([0.0],[0.0],quiver=([1.0],[0.0]), label="B field", markerstrokewidth=4.0, color=:purple)
    plt_pm = Plots.plot!([0.0,1.0],[0.0,0.0], label="B field", linewidth=3, color=:purple, annotations=(0.98,0.05,Plots.text("B",:right)))
    plt_pm = Plots.quiver!([0.0],[0.0],quiver=([xp],[yp]), label="v", markerstrokewidth=4.0, color=:black)
    plt_pm = Plots.plot!([0.0,xp],[0.0,yp], label="v", linewidth=3, color=:black, annotations=(xp+0.02,yp+0.03,Plots.text("v")))
    plt_pm = Plots.quiver!([0.0],[0.0],quiver=([0.0],[yp]), label="vperp", markerstrokewidth=4.0, color=:red)
    plt_pm = Plots.plot!([0.0,0.0],[0.0,yp], label="vperp", linewidth=3, color=:red, annotations=(0.02,yp+0.03,Plots.text("vperp")))
    plt_pm = Plots.quiver!([0.0],[0.0],quiver=([xp],[0.0]), label="vpara", markerstrokewidth=4.0, color=:blue)
    plt_pm = Plots.plot!([0.0,xp],[0.0,0.0], label="vpara", xaxis=[xmin,1.0], yaxis=[0.0,1.0], linewidth=3, legend=true, color=:blue, fontsize=14.0, fontweight="bold", annotations=(xp,0.05,Plots.text("vpara")))
    plt_pm = Plots.plot!([0.0,xp],[yp,yp],linestyle=:dash,label="",color=:black)
    plt_pm = Plots.plot!([xp,xp],[0.0,yp],linestyle=:dash,label="",color=:black, title="pm: $(round(pm, digits=3))")

    if !poltor
        myplt = Plots.plot(plt_topo, plt_pm, plt_top, plt_crs, layout=(2,2), size=(800,800))
    else
        myplt = Plots.plot(plt_topo, plt_pol, plt_tor, plt_pm, plt_top, plt_crs, layout=(3,2), size=(800, 1200), left_margin=5Plots.mm)
    end
    Plots.plot(myplt)
end

# ╔═╡ d3c783ea-2867-4e1e-968a-a6fe982cfa74


# ╔═╡ 09af6e22-0663-4e5d-8f0c-b06055270e0e


# ╔═╡ 4419b398-8ebf-4c15-a45a-6ff0d287d585


# ╔═╡ 34acf9b7-8e20-4e2e-84ca-e9bcc051609a


# ╔═╡ Cell order:
# ╠═a761bf9e-a159-11f1-2ec7-711778cfac0d
# ╠═cd6a92a9-7113-4281-8596-a520642c5cbe
# ╠═29d500df-9ab5-449c-b1a2-a081e3cd770f
# ╠═23772204-992e-4bd3-ad1f-f1d86c33f6aa
# ╠═dfacf685-d70b-4433-9d7c-e24fac4cf794
# ╠═206460f9-93e3-48cc-bd65-fbf0ce0a9967
# ╠═e9cd967d-c461-4324-99e0-6110bf7d43fc
# ╠═a37e90da-5be0-4b14-8613-1e146dca5ec0
# ╠═684923b0-dcf4-4ff3-9d50-8291ddb0b593
# ╠═d3c783ea-2867-4e1e-968a-a6fe982cfa74
# ╠═09af6e22-0663-4e5d-8f0c-b06055270e0e
# ╠═4419b398-8ebf-4c15-a45a-6ff0d287d585
# ╠═34acf9b7-8e20-4e2e-84ca-e9bcc051609a
