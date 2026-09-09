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

# ╔═╡ 23423083-edd9-4c93-8eaf-39038f8447a9
begin
	# SPECIFY THE INPUTS IN THIS CELL
	# Please specify the OWCF folder and let the notebook change directory to the 
	# OWCF folder when the cell below is executed. This is to be able to load the
	# correct versions of the Julia packages as specified in the Project.toml and 
	# Manifest.toml files.
	folderpath_OWCF = "REPLACE-THIS-TEXT-WITH-THE-PATH-TO-THE-OWCF-FOLDER-ON-YOUR-COMPUTER" # Finish with '/'
	
	filepath_distr = "REPLACE-THIS-TEXT-WITH-THE-PATH-TO-YOUR-OWN-FAST-ION-DISTRIBUTION-FILE-ON-YOUR-COMPUTER"
	filepath_equil = folderpath_OWCF*"equilibrium/JET/g99971/g99971_474-48.9.eqdsk" 
	filepath_tm = folderpath_OWCF*"apps/example_data/EpRzTopoMap_JET_99971L72_at48,9s_D_51x41x11x13_wLost.jld2" 
	filepath_tb = folderpath_OWCF*"apps/example_data/topoBounds_JET_99971L72_at48,9s_D_51x52x53wLost.jld2" # Energy (E) grid points in (E,pm,Rm) must match E grid points of topoMap in (E,p,R,z)
	isfile(filepath_tb) && (nR_for_mu=500) # If an (E,pm,Rm) topoBounds file has been specified, specify how many R grid points should be used to resolve the magnetic moment across the plasma
	FI_species = "D" # Example deuterium: "D"
	verbose = true

	# h5_is_rowmajor - If the 'filepath_EpRz' file is an .h5/.hdf5 file and it was saved using a row-major programming language (e.g. C/C++ or NumPy in Python, please see https://en.wikipedia.org/wiki/Row-_and_column-major_order for more info), this input variable should be set to true - Bool
	h5_is_rowmajor = false
	
	# You can specify additional entries in the dictionary below, to include more keyword arguments in 
	# the guiding-centre orbit integration algorithm
	extra_kw_args = Dict(:limit_phi => true, :maxiter => 0, :toa => true)
	# limits the number of toroidal turns for orbits
	# The orbit integration algorithm will try progressively smaller timesteps these number of times
	# toa stands for "try only adaptive" orbit trajectory integration
end

# ╔═╡ 4556c934-451a-4392-8c15-3cbc32c464e2
begin
	cd(folderpath_OWCF)
	using Pkg
	Pkg.activate(".")
	
	# Loading packages
	verbose && println("Loading packages... ")
	using EFIT
	using Equilibrium
	using GuidingCenterOrbits
	using JLD2
	using Plots
	using HDF5
	using FileIO
	using PlutoUI
	include(folderpath_OWCF*"extra/dependencies.jl")
end

# ╔═╡ a82b82d8-a6d7-11f1-0495-6174c0501847
md"""
# EpRzWebApp

## Description:
This notebook provides a web application that can be used to visualize (E,p,R,z) phase space and (guiding-centre) orbits in a tokamak in an interactive and intuitive manner. The orbit topological regions of (E,p,R,z) space are visualized and explored interactively. The web app can also be used to visualize a fast-ion distribution 
and superimpose topological boundaries to examine the distributions in terms of orbits in (E,p,R,z) space. In addition, if computed, the app can also be used to visualize maps of poloidal and toroidal transit times.

Prior to running this notebook, please make sure you have run the calcEpRzTopoMap.jl script (or equivalent). Also, please make sure that you have noted the path to the saved topological map output data.

## Inputs:
- folderpath_OWCF - The path to the OWCF folder on your computer. Needed for correct loading - String
- filepath_distr - (Optional) The path to a .h5/.jld2-file containing the (E,p,R,z) fast-ion distribution. Please see the functions h5to4D() and JLD2to4D() in the OWCF/extra/dependencies.jl library of functions, for a description of necessary file data keys - String
- filepath_equil - The path to the .eqdsk-file (or .jld2-file) with the tokamak magnetic equilibrium and geometry. A .jld2 file with a magnetic equilibrium can be computed with the OWCF/extra/createCustomMagneticEquilibrium.jl tool - String
- filepath_tm - The path to the .jld2-file containing the topological map (and more). Should be an output file of the OWCF/calcEpRzTopoMap.jl script - String
- filepath_tb - (Optional) The path to a .jld2-file containing topological boundaries in (E,pm,Rm). Should be an output file of the OWCF/helper/extractTopoBounds.jl tool or the OWCF/calcTopoMap.jl script - String
- FI_species - The fast-ion species, e.g. "D", "T", "alpha", "3he" etc. Please see the OWCF/misc/species_func.jl for an overview of all available particle species in the OWCF - String
- verbose - If set to true, the app will talk a lot! - Bool

## Outputs:
# -

## Saved files:
# -

### Notebook written by Henrik Järleblad (henrikj@dtu.dk) and Andrea Valentini (anvalen@dtu.dk)
### Last maintained 2026-09-02.
"""

# ╔═╡ 97aff330-3f29-47e1-a5fc-ac284a0ebbab
begin
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
    close(myfile)
    
    # If provided, read the .jld2-file for displaying the topological boundaries in (E,pm,Rm)
    topobo = false
    if isfile(filepath_tb)
        verbose && println("Loading (E,pm,Rm) topological boundaries... ")
        myfile = jldopen(filepath_tb,false,false,false,IOStream)
        topoBounds_OS = myfile["topoBounds"]
        E_array_tb = myfile["E_array"]
        pm_array = myfile["pm_array"]
        Rm_array = myfile["Rm_array"]
        close(myfile)
    
        if !(E_array==E_array_tb)
            error("Energy (E) grid points of (E,pm,Rm) topological boundaries did not match E grid points of topological map in (E,p,R,z). Please correct and re-try!")
        end
    
        topobo = true
        println("topobo is true")
    end
    
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
    
            o = get_orbit(M,getGCP(FI_species;E=E_array[iE],p=p_array[ip],R=R_array[iR],z=z_array[iz]); wall=wall, extra_kw_args...)
    
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
            F, E_array_FI, p_array_FI, R_array_FI, z_array_FI = h5to4D(filepath_distr; rowmajor=h5_is_rowmajor)
        elseif fileext_distr=="jld2"
            F, E_array_FI, p_array_FI, R_array_FI, z_array_FI = JLD2to4D(filepath_distr)
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
            verbose && println("Grid points of topological map and loaded fast-ion data did not match.")
            verbose && println("---> nE_topo: $(length(E_array))   nE_FI: $(length(E_array_FI))   extrema(E_topo): $(extrema(E_array))   extrema(E_FI): $(extrema(E_array_FI))")
            verbose && println("---> np_topo: $(length(p_array))   np_FI: $(length(p_array_FI))   extrema(p_topo): $(extrema(p_array))   extrema(p_FI): $(extrema(p_array_FI))")
            verbose && println("---> nR_topo: $(length(R_array))   nR_FI: $(length(R_array_FI))   extrema(R_topo): $(extrema(R_array))   extrema(R_FI): $(extrema(R_array_FI))")
            verbose && println("---> nz_topo: $(length(z_array))   nz_FI: $(length(z_array_FI))   extrema(z_topo): $(extrema(z_array))   extrema(z_FI): $(extrema(z_array_FI))")
                               
            verbose && println("Interpolating... ")
            F = interpFps(F,E_array_FI,p_array_FI,R_array_FI,z_array_FI,E_array,p_array,R_array,z_array)
        end
    
        dE4D, dp4D, dR4D, dz4D = get4DDiffs(E_array, p_array,R_array,z_array) 
        F_Rz = dropdims(sum(dE4D .* dp4D .* F,dims=(1,2)),dims=(1,2)) # Will need for (R,z) heatmap plot
    end
    
    verbose && println("Preparing final data for web application... ")
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
    dE = abs.(diff(E_array))[1] # Assume equidistant
    dp = abs.(diff(p_array))[1] # Assume equidistant
    dR = abs.(diff(R_array))[1] # Assume equidistant
    dz = abs.(diff(z_array))[1] # Assume equidistant
    
    keV = 1000*(GuidingCenterOrbits.e0)
    
    c_array = [:red, :blue, :green, :purple, :orange, :yellow, :brown, :pink, :gray]
end

# ╔═╡ f1684069-7a13-4fab-9a38-b2f1f525dddc
function plot_inputs()

	return PlutoUI.combine() do Child

		inputs=[
			md""" Show tokamak wall: $(
				Child("wall", Switch(default=true))
			)""",
			md""" Study: $(
				Child("study", Select([:orbs => "Orbit topology", :FI => "Fast-ion distribution", :tpol => "Poloidal transit times", :ttor => "Toroidal transit times", :OS => "(E,pm,Rm) at (R,z)"]))
			)""",
			md""" Coordinate space: $(
				Child("coordspace", Select([:EP => "(E,p)", :VEL => "(v_para,v_perp)"]))
			)""",
			md""" Energy, E (keV): $(
				Child("E", Slider(E_array, default=E_array[Int64(round(length(E_array)/2))]))
			)""",
			md""" Pitch, p (-): $(
				Child("p", Slider(p_array, default=p_array[Int64(round(3*length(p_array)/4))]))
			)""",
			md""" Major radius, R (m): $(
				Child("R", Slider(R_array, default=R_array[Int64(round(length(R_array)/2))]))
			)""",
			md""" Vertical pos., z (m): $(
				Child("z", Slider(z_array, default=z_array[Int64(round(length(z_array)/2))]))
			)"""
		]

		md"""
		#### Plot controls
		$(inputs)
		"""
	end
end

# ╔═╡ d8211649-132a-45e3-b993-546bbd565a0d
@bind my_plot_inputs plot_inputs()

# ╔═╡ 72fff0c5-83d1-4a96-b04c-8590275b8bca
let
	tokamak_wall = my_plot_inputs.wall
	study = my_plot_inputs.study
	space = my_plot_inputs.coordspace
	E = my_plot_inputs.E
	p = my_plot_inputs.p
	R = my_plot_inputs.R
	z = my_plot_inputs.z

	my_gcp = getGCP(FI_species,E=E,p=p,R=R,z=z)

    o = get_orbit(M,my_gcp; wall=wall, extra_kw_args...)
    
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
    if study==:FI && isfile(filepath_distr)
        plt_crs = Plots.heatmap(R_array, z_array, F_Rz', fillcolor=cgrad([:white, :darkblue, :green, :yellow, :orange, :red]))
        plt_crs = Plots.scatter!([magnetic_axis(M)[1]],[magnetic_axis(M)[2]],label="Magnetic axis", mc=:gray, aspect_ratio=:equal, xlabel="R [m]", ylabel="z [m]",title="Poloidal cross-section")
    else
        plt_crs = Plots.scatter([magnetic_axis(M)[1]],[magnetic_axis(M)[2]],label="Magnetic axis", mc=:gray, aspect_ratio=:equal, xlabel="R [m]", ylabel="z [m]",title="Poloidal cross-section")
    end
    if tokamak_wall
        plt_crs = Plots.contour!(flux_r,flux_z,psi_rz',levels=collect(range(psi_mag,stop=psi_bdry,length=5)),color=:gray, linestyle=:dot,linewidth=2.5, label="",colorbar=false)
        wall_dR = maximum(wall.r)-minimum(wall.r)
        plt_crs = Plots.plot!(wall.r,wall.z, label="Tokamak wall", color=:black, linewidth=1.5,xaxis=[minimum(wall.r)-wall_dR/10,maximum(wall.r)+wall_dR])
    end
    plt_crs = Plots.plot!(o.path.r,o.path.z, label="$(o.class) orbit", color=orb_color, linestyle=orb_linestyle, linewidth=2.5)
    plt_crs = Plots.vline!([R],linestyle=:dot,color=:gray,linewidth=1.0,label="")
    plt_crs = Plots.hline!([z],linestyle=:dot,color=:gray,linewidth=1.0,label="")
    plt_crs = Plots.scatter!([R],[z], mc=orb_color, label="(R,z)", grid=false)

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
        if space==:EP # (E,p) plot
            plt_topo = Plots.heatmap(E_array,p_array, F[:,:,Rci,zci]', fillcolor=cgrad([:white, :darkblue, :green, :yellow, :orange, :red]), xlabel="E [keV]", ylabel="p [-]", title="f(E,p) at (R,z)", right_margin=3Plots.mm)
        else # :VEL => (v_para, v_perp) plot
            F_VEL, vpara_array, vperp_array = Ep2VparaVperp(E_array, p_array, F[:,:,Rci,zci]; my_gcp=my_gcp, needJac=true, returnAbscissas=true)
            plt_topo = Plots.heatmap(vpara_array, vperp_array, F_VEL', fillcolor=cgrad([:white, :darkblue, :green, :yellow, :orange, :red]), xlabel="v_para [m/s]", ylabel="v_perp [m/s]", title="f(v_para,v_perp) at (R,z)", right_margin=3Plots.mm )
        end
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
        if space==:EP
            z_matrix, x_array, y_array,  = pTT_microsecs, E_array, p_array
            title = "tau_pol(E,p) at (R,z) [microseconds] \n tau_pol($(round(E,sigdigits=3)),$(round(p,sigdigits=2)))=$(round(o.tau_p /(1.0e-6),sigdigits=3)) microseconds"
            xlabel = "E [keV]"
            ylabel = "p [-]"
        else
            z_matrix, x_array, y_array  = Ep2VparaVperp(E_array, p_array, pTT_microsecs; my_gcp=my_gcp, returnAbscissas=true)
            title = "tau_pol(vpara,vperp) at (R,z) [microseconds] \n tau_pol($(round(p*sqrt(2*keV*E/my_gcp.m),sigdigits=3)),$(round(sqrt(1-p^2)*sqrt(2*keV*E/my_gcp.m),sigdigits=3)))=$(round(o.tau_p /(1.0e-6),sigdigits=3)) microseconds"
            xlabel = "vpara [m/s]"
            ylabel = "vperp [m/s]"
        end
        if !((max_OOM-min_OOM)==0.0) && (length(nz_coords) > 1) # All values NOT within same order of magnitude AND more than one non-zero element. Use logarithmic colorbar
            plt_topo = Plots.heatmap(x_array,y_array, z_matrix', title=title, fillcolor=cgrad([:white, :darkblue, :green, :yellow, :orange, :red]), xlabel=xlabel, ylabel=ylabel, colorbar=true, colorbar_scale=:log10, clims = (10^min_OOM, 10^max_OOM), top_margin=3Plots.mm) # Get nice powers-of-ten limits for the colorbar
        else # Else, use linear colorbar 
            plt_topo = Plots.heatmap(x_array,y_array, z_matrix', title=title, fillcolor=cgrad([:white, :darkblue, :green, :yellow, :orange, :red]), xlabel=xlabel, ylabel=ylabel, colorbar=true, top_margin=3Plots.mm)
        end
        ms = 1.8
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
        if space==:EP
            z_matrix, x_array, y_array = tTT_microsecs, E_array, p_array
            title = "tau_tor(E,p) at (R,z) [microseconds] \n tau_tor($(round(E,sigdigits=3)),$(round(p,sigdigits=2)))=$(round(o.tau_t /(1.0e-6),sigdigits=3)) microseconds"
            xlabel = "E [keV]"
            ylabel = "p [-]"
        else
            z_matrix, x_array, y_array = Ep2VparaVperp(E_array, p_array, tTT_microsecs; my_gcp=my_gcp, returnAbscissas=true)
            title = "tau_tor(vpara,vperp) at (R,z) [microseconds] \n tau_tor($(round(p*sqrt(2*keV*E/my_gcp.m),sigdigits=3)),$(round(sqrt(1-p^2)*sqrt(2*keV*E/my_gcp.m),sigdigits=3)))=$(round(o.tau_t /(1.0e-6),sigdigits=3)) microseconds"
            xlabel = "vpara [m/s]"
            ylabel = "vperp [m/s]"
        end
        if !((max_OOM-min_OOM)==0.0) && (length(nz_coords) > 1) # All values NOT within same order of magnitude AND more than one non-zero element. Use logarithmic colorbar
            plt_topo = Plots.heatmap(x_array,y_array, z_matrix', title=title, fillcolor=cgrad([:white, :darkblue, :green, :yellow, :orange, :red]), xlabel=xlabel, ylabel=ylabel, colorbar=true, colorbar_scale=:log10, clims = (10^min_OOM, 10^max_OOM),top_margin=3Plots.mm) # Get nice powers-of-ten limits for the colorbar
        else # Else, use linear colorbar 
            plt_topo = Plots.heatmap(x_array,y_array, z_matrix', title=title, fillcolor=cgrad([:white, :darkblue, :green, :yellow, :orange, :red]), xlabel=xlabel, ylabel=ylabel, colorbar=true,top_margin=3Plots.mm)
        end
        ms = 1.8
        #plt_weights = Plots.scatter!(E_scatvals_tb,p_scatvals_tb,markersize=ms,label="",markercolor=:black,leg=false)
    elseif study==:OS && topobo
        ms = 1.8
        ones_carinds = findall(x-> x==1.0,topoBounds_OS[Eci,:,:])
        pm_scatvals_tb = zeros(length(ones_carinds))
        Rm_scatvals_tb = zeros(length(ones_carinds))
        for (ind,carinds) in enumerate(ones_carinds)
            pm_scatvals_tb[ind] = pm_array[carinds[1]]
            Rm_scatvals_tb[ind] = Rm_array[carinds[2]]
        end

        array_o_orb_tuples = map(i-> (get_orbit(M,getGCP(FI_species,E=E,p=p_array[i],R=R,z=z); wall=wall, toa=true), topoMap[Eci,i,Rci,zci]), 1:length(p_array)) # Genius line <3
        
        plt_topo = Plots.scatter(Rm_scatvals_tb,pm_scatvals_tb,markersize=ms,leg=false,markercolor=:black, xlabel="Rm [m]", ylabel="pm")
        for orb_tuple in array_o_orb_tuples
            oo, ic = orb_tuple
            if !(o.class == :lost)
                plt_topo = Plots.scatter!([oo.coordinate.r],[oo.coordinate.pitch],mc=c_array[Int64(ic)])
            else
                ### WRITE A SPECIAL ALGORITHM TO DETERMINE Rm FOR LOST ORBITS. If possible/makes sense? Some lost orbits might not go through midplane?
            end
        end
    else
        topoMap_raw = topoMap[:,:,Rci,zci] # might not include all integers between 1 and 9 (needed for correct Set1_9 heatmap coloring)
        if space==:EP
            topoMap_ext = ones(size(topoMap_raw,1),size(topoMap_raw,2)+1) # We ensure 1.0 in data by using ones() function. Add one extra column
            topoMap_ext[1:size(topoMap_raw,1),1:size(topoMap_raw,2)] .= topoMap_raw # put in the true topoMap to be plotted
            topoMap_ext[end,end] = 9.0 # Get the 9.0. This is such an ugly solution (to acquire correct orbit type coloring)...
            p_array_ext = vcat(p_array,2*p_array[end]-p_array[end-1]) # Extend p_array by one dummy element
            plt_topo = Plots.heatmap(E_array,p_array_ext,topoMap_ext',color=:Set1_9,legend=false,xlabel="E [keV]", ylabel="p [-]", title="(E,p) orbit topology at (R,z)", ylims=(minimum(p_array),maximum(p_array)))
        else
            topoMap_raw_VEL, vpara_array, vperp_array = Ep2VparaVperp(E_array, p_array, topoMap_raw; my_gcp=my_gcp, isTopoMap=true, returnAbscissas=true)
            topoMap_ext_VEL = ones(size(topoMap_raw_VEL,1),size(topoMap_raw_VEL,2)+1) # We ensure 1.0 in data by using ones() function. Add one extra column
            topoMap_ext_VEL[1:size(topoMap_raw_VEL,1),1:size(topoMap_raw_VEL,2)] .= topoMap_raw_VEL # put in the true topoMap to be plotted
            topoMap_ext_VEL[end,end] = 9.0 # Get the 9.0. This is such an ugly solution (to acquire correct orbit type coloring)...
            vperp_array_ext = vcat(vperp_array,2*vperp_array[end]-vperp_array[end-1]) # Extend vperp_array by one dummy element
            plt_topo = Plots.heatmap(vpara_array,vperp_array_ext,topoMap_ext_VEL',color=:Set1_9,legend=false,xlabel="vpara [m/s]", ylabel="vperp [m/s]", title="(vpara,vperp) orbit topology at (R,z)", ylims=(minimum(vperp_array),maximum(vperp_array)), xlims=(minimum(vpara_array),maximum(vpara_array)))
        end
    end
    if !(study==:OS && topobo)
        if space==:EP
            plt_topo = Plots.vline!([E],linestyle=:dot,color=:gray,linewidth=1.0)
            plt_topo = Plots.hline!([p],linestyle=:dot,color=:gray,linewidth=1.0)
            plt_topo = Plots.scatter!([E],[p],markershape=:circle,mc=orb_color,legend=false,markersize=6, grid=false)
        else
            plt_topo = Plots.vline!([p*sqrt(2*keV*E/my_gcp.m)],linestyle=:dot,color=:gray,linewidth=1.0)
            plt_topo = Plots.hline!([sqrt(1-p^2)*sqrt(2*keV*E/my_gcp.m)],linestyle=:dot,color=:gray,linewidth=1.0)
            plt_topo = Plots.scatter!([p*sqrt(2*keV*E/my_gcp.m)],[sqrt(1-p^2)*sqrt(2*keV*E/my_gcp.m)],markershape=:circle,mc=orb_color,legend=false,markersize=6, grid=false)
        end
    else
        if !(o.class == :lost)
            plt_topo = Plots.vline!([o.coordinate.r],linestyle=:dot,color=:gray,linewidth=1.0)
            plt_topo = Plots.hline!([o.coordinate.pitch],linestyle=:dot,color=:gray,linewidth=1.0)
            plt_topo = Plots.scatter!([o.coordinate.r],[o.coordinate.pitch],mc=:black, markersize=6)
        else
            ### WRITE A SPECIAL ALGORITHM TO DETERMINE Rm FOR LOST ORBITS. If possible/makes sense? Some lost orbits might not go through midplane?
        end
    end
    
    # Orbit fractions plot
    # Please see signalWebApp.jl for more info
    o_types=[1.0,2.0,3.0,4.0,5.0,8.0]
    valid_orb_inds = findall(x-> x in o_types,topoMap[:,:,Rci,zci])
    if study==:FI && !(filepath_distr=="")
        F_Ep = F[:,:,Rci,zci]
        sum_o_FI_distr = sum(F_Ep[valid_orb_inds])
        #n_FI = sum((2*pi*R*dE*dp*dR*dz) .*F_Ep[valid_orb_inds])
        #println(n_FI)
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

    #pitch visualization plot
    x_array = collect(range(0.0,stop=o.tau_p,length=length(o.path.pitch)))
    my_xlabel = "Poloidal time [microseconds]"
    plt_pitc = Plots.plot(x_array ./(1.0e-6), o.path.pitch,color=orb_color,title="Pitch along orbit path",label="",linewidth=1.5, xlabel=my_xlabel, right_margin=4Plots.mm)

	Plots.plot(plt_crs, plt_orb_fracs, plt_topo, plt_pitc, layout=(2,2), size=(800, 800))
end

# ╔═╡ 0bcaa411-0c18-4db6-afc6-02de1586ee1d


# ╔═╡ a61e1bb3-16b6-491a-8700-1a8526a2c431


# ╔═╡ 6a1260d5-ccab-49b4-b52b-d1f4f5f03b32


# ╔═╡ Cell order:
# ╠═a82b82d8-a6d7-11f1-0495-6174c0501847
# ╠═23423083-edd9-4c93-8eaf-39038f8447a9
# ╠═4556c934-451a-4392-8c15-3cbc32c464e2
# ╠═97aff330-3f29-47e1-a5fc-ac284a0ebbab
# ╠═f1684069-7a13-4fab-9a38-b2f1f525dddc
# ╠═d8211649-132a-45e3-b993-546bbd565a0d
# ╠═72fff0c5-83d1-4a96-b04c-8590275b8bca
# ╠═0bcaa411-0c18-4db6-afc6-02de1586ee1d
# ╠═a61e1bb3-16b6-491a-8700-1a8526a2c431
# ╠═6a1260d5-ccab-49b4-b52b-d1f4f5f03b32
