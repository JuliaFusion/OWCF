################################ modeAnalysisWebApp.jl ########################################################

#### Description:
# This script provides a web application to visualize resonances between fast ions and magnetohydrodynamic (MHD)
# Alfvén eigenmodes in an interactive and intuitive manner. The resonance criterion is based on the equation
#
# ω = n ω_ϕ + m ω_θ + ω_rot
# 
# where ω is the mode (angular) frequency, n is the toroidal mode number, ω_ϕ is the toroidal transit (angular) frequency 
# of the fast ion, m is the poloidal mode number, ω_θ is the poloidal transit (angular) frequency of the fast ion and ω_rot
# is the plasma rotation (angular) frequency.
# ω_ϕ and ω_θ are computed by integrating the guiding-centre equations of motion, and then evaluating the toroidal (τ_ϕ)
# and poloidal (τ_θ) transit times. ω, n and m can be interactively changed by the user using the sliders in the web application.
#
# τ_ϕ and τ_θ are related to ω_ϕ and ω_θ as ω_ϕ = 2*pi/τ_ϕ and ω_θ = 2*pi/τ_θ, respectively.
#
# The resonance between mode and fast ion is then plotted in (E,pm,Rm) (or, via the press of a button in the app, in (E,Λ,Pϕ_n;σ))
# and is visualized as a heatmap where the values (A) of the figure colorbar are computed according to the equation
#
# A = - log10(|ω - n ω_ϕ - m ω_θ - ω_rot|)
#
# where larger values of A thus corresponds to greater resonance, and smaller values of A corresponds to less resonance.
#
# PLEASE NOTE! To be able to use the app, a topological map containing toroidal (tau_t) and poloidal (tau_p) transit time data 
# has to be provided. This type of data is to be provided as a .jld2 file with the keys 'torTransTimes' and 'polTransTimes'.
# Such a file can be easily obtained by using the calcTopoMap.jl OWCF script.

#### Inputs (units given when defined in the script):
# enable_COM - If true, (E,pm,Rm) -> (E,Λ,Pϕ_n;σ) can be performed via a toggle button. Set to false to minimize computation time - Bool
# folderpath_OWCF - The path to the OWCF folder on your computer. Needed for correct loading - String
# filepath_equil - The path to the file with the tokamak magnetic equilibrium and geometry - String
# filepath_tm - The path to the .jld2-file containing the topological map (and toroidal and poloidal transit times) - String
# FI_species - The species of the particle being simulated (deuterium, tritium etc) - String
# ω_init - The mode frequency to investigate (in kHz). This can be interactively changed later in the app - Float64
# ω_rot_init - The plasma rotation frequency (in kHz). This can be interactively changed later in the app - Float64
# n - The toroidal mode number to investigate. This can be interactively changed later in web application - Int64
# m - The toroidal mode number to investigate. This can be interactively changed later in web application - Int64
# verbose - If set to true, the app will talk a lot! - Bool
# port - An integer. It will correspond to the internet port of your computer through which the app is run - Int64

#### Outputs
# -

#### Saved files
# -

#### Other
# Warning! Please note! For topological maps containing more than approximately 150 000 valid orbits (e.g. 20x100x100),
# you should NOT use modeAnalysisWebApp.jl (or any other interactive app). As of OWCF version 1.2.1, the 
# web interface simply becomes too slow. Please do instead plot resonances manually (in scratch.jl for example).
# You do this by, e.g. for an (n,m)=(1,2) mode with ω=150 kHz, coding the following
#
# folderpath_OWCF = "/the/path/to/the/OWCF/folder/"
# cd(folderpath_OWCF)
# using Pkg
# Pkg.activate(".")
# using JLD2
# using Plots
#
# myfile = jldopen(filepath_tm, false, false, false, IOStream)
# topoMap = myfile["topoMap"]
# E_array = myfile["E_array"]
# pm_array = myfile["pm_array"]
# Rm_array = myfile["Rm_array"]
# torTransTimes = myfile["torTransTimes"]
# polTransTimes = myfile["polTransTimes"]
# close(myfile)
#
# n = 1
# m = 2
# ω = 150 # kHz
# E = 150 keV # Example energy 'slice' of 150 keV
# iE = argmin(abs.(E_array - E)) # Find the closest value to E in E_array
# npm = length(pm_array)
# nRm = length(Rm_array)
# 
# resonance = -log10.((2*pi/(ω*1000)) .*ones(size(npm,nRm)) - n .*torTransTimes[iE,:,:] - m .*polTransTimes[iE,:,:])
#
# Plots.heatmap(Rm_array,pm_array,resonance,legend=false,xlabel="Rm [m]", ylabel="pm", title="E: $(round(E,digits=3)) keV", fillcolor=cgrad([:white, :darkblue, :green, :yellow, :orange, :red]))

# Script written by Henrik Järleblad. Last maintained 2025-01-23.
########################################################################################################

#########################################################################################
# Please specify the OWCF folder and let the app change directory to the 
# OWCF folder when orbitsWebApp.jl is executed. This is to be able to load the
# correct versions of the Julia packages as specified in the Project.toml and 
# Manifest.toml files.
folderpath_OWCF = "" # Finish with '/'
cd(folderpath_OWCF)
using Pkg
Pkg.activate(".")

#########################################################################################
# Must load JLD2 package first, to be able to check filepath_tm and filepath_W for 'extra_kw_args'
using JLD2

#########################################################################################
# Inputs
enable_COM = false # Set to false for large grids (>10x100x100 in (E,pm,Rm)). The computation time simply becomes too large. Or...
if enable_COM
    filepath_tm_COM = "" # ...please specify an output of the os2com.jl script (which contains the key "topoMap")(and possibly "polTransTimes" and "torTransTimes"). Leave unspecified if orbitsWebApp.jl should compute the (E,pm,Rm) -> (E,Λ,Pϕ_n;σ) map itself.
end
filepath_equil = ""  # Example JET shot 96100 at 13s (53 minus 40): g96100/g96100_0-53.0012.eqdsk" #
filepath_tm = ""
FI_species = "" # Example deuterium: "D"
ω_init = 0.0 # Mode frequency. kHz
ω_rot_init = 0.0 # Plasma rotation frequency. kHz
n_init = 1 # Toroidal mode number
m_init = 2 # Poloidal mode number
verbose = true
port = 2424 # The port number hosting the web app. orbitsWebApp will be accessible via your web browser at the web address 'localhost:port' where you replace 'port' with the integer number you have specified
# EXTRA KEYWORD ARGUMENTS BELOW (these will go into the orbit_grid() and get_orbit() functions from GuidingCenterOrbits.jl)
myfile = jldopen(filepath_tm,false,false,false,IOStream)
if haskey(myfile,"extra_kw_args")
    extra_kw_args = myfile["extra_kw_args"]
else
    verbose && println("No extra keyword arguments for orbit integration found in topological map file. Assuming :toa, :limit_phi and :max_tries=0")
    extra_kw_args = Dict(:limit_phi => true, :max_tries => 0)
    # limit_phi: limits the number of toroidal turns for orbits
    # max_tries: The orbit integration algorithm will try progressively smaller timesteps these number of times
end
close(myfile)
if haskey(extra_kw_args, :max_tries) # To allow compatibility with older OWCF files
    delete!(extra_kw_args, :max_tries)
end

#########################################################################################
# Loading packages
verbose && println("Loading packages... ")
using Interact
using EFIT
using Equilibrium
using GuidingCenterOrbits
using Plots
using FileIO
using Mux
using WebIO
include(folderpath_OWCF*"misc/species_func.jl")
include(folderpath_OWCF*"extra/dependencies.jl")

#########################################################################################
# Loading tokamak equilibrium
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

#########################################################################################
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

#########################################################################################
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

#########################################################################################
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

#########################################################################################
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
#if !(sum(isnan.(polTransTimes))==0)
#    polTransTimes = map(x-> isnan(x) ? 0.0 : x, polTransTimes)
#    @warn "Map of poloidal transit times in (E,pm,Rm) had NaNs in it! orbitsWebApp.jl automatically set them to 0.0."
#end
#if !(sum(isnan.(torTransTimes))==0)
#    torTransTimes = map(x-> isnan(x) ? 0.0 : x, torTransTimes)
#    @warn "Map of toroidal transit times in (E,pm,Rm) had NaNs in it! orbitsWebApp.jl automatically set them to 0.0."
#end
#if enable_COM
#    if !(sum(isnan.(polTransTimes_COM))==0)
#        polTransTimes_COM = map(x-> isnan(x) ? 0.0 : x, polTransTimes_COM)
#        @warn "Map of poloidal transit times in (E,Λ,Pϕ_n;σ) had NaNs in it! orbitsWebApp.jl automatically set them to 0.0."
#    end
#    if !(sum(isnan.(torTransTimes_COM))==0)
#        torTransTimes_COM = map(x-> isnan(x) ? 0.0 : x, torTransTimes_COM)
#        @warn "Map of toroidal transit times in (E,Λ,Pϕ_n;σ) had NaNs in it! orbitsWebApp.jl automatically set them to 0.0."
#    end
#end

#########################################################################################
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

ω_θ = map(x-> !isnan(x) ? 2*pi/x : x, polTransTimes)
ω_ϕ = map(x-> !isnan(x) ? 2*pi/x : x, torTransTimes)
enable_COM && (ω_θ_COM = map(x-> !isnan(x) ? 2*pi/x : x, polTransTimes_COM))
enable_COM && (ω_ϕ_COM = map(x-> !isnan(x) ? 2*pi/x : x, torTransTimes_COM))

#########################################################################################
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
#########################################################################################
## --------------------------------------------------------------------------
# The web application
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

verbose && println("--- You can access the modeAnalysisWebApp via an internet web browser when you see 'Task (runnable)...' ")
verbose && println("--- When 'Task (runnable)...' has appeared, please visit the website localhost:$(port) ---")
verbose && println("--- Remember: It might take a minute or two to load the webpage. Please be patient. ---")
function app(req) # Start the Interact.jl app
    @manipulate for E=E_array, pm=pm_array, Rm=Rm_array, ω=ω_init, ω_rot=ω_rot_init, n=n_init, m=m_init, phase_space = Dict("(E,Λ,Pϕ_n;σ)" => :COM, "(E,pm,Rm)" => :OS), tokamak_wall = Dict("on" => true, "off" => false), show_coordinate = Dict("on" => true, "off" => false), save_plots = Dict("on" => true, "off" => false)

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
        if save_plots
            if (phase_space==:OS) || !enable_COM
                png(plt_top, "plt_top_$(round(E, digits=2))_$(round(o.coordinate.pitch, digits=2))_$(round(o.coordinate.r,digits=2))")
            else
                png(plt_top, "plt_top_$(round(E, digits=2))_$(round(Λ, sigdigits=2))_$(round(Pϕ_n,sigdigits=2))")
            end
        end

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
        if save_plots
            if (phase_space==:OS) || !enable_COM
                png(plt_crs, "plt_crs_$(round(E, digits=2))_$(round(o.coordinate.pitch, digits=2))_$(round(o.coordinate.r,digits=2))")
            else
                png(plt_crs, "plt_crs_$(round(E, digits=2))_$(round(Λ, sigdigits=2))_$(round(Pϕ_n,sigdigits=2))")
            end
        end

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
        if save_plots
            plt_topo = Plots.plot!(dpi=600)
            if (phase_space==:OS) || !enable_COM
                png(plt_topo, "plt_topo_$(round(E, digits=2))_$(round(o.coordinate.pitch, digits=2))_$(round(o.coordinate.r,digits=2))")
            else
                png(plt_topo, "plt_topo_$(round(E, digits=2))_$(round(Λ, sigdigits=2))_$(round(Pϕ_n,sigdigits=2))")
            end
        end

        # Resonance plot
        if (phase_space==:OS) || !enable_COM
            npm = length(pm_array)
            nRm = length(Rm_array)
            resonance = -log10.(abs.((ω*1000) .*ones(npm,nRm) - n .*ω_ϕ[Eci,:,:] - m .*ω_θ[Eci,:,:] - (ω_rot*1000) .*ones(npm,nRm)))
            plt_res = Plots.heatmap(Rm_array,pm_array,resonance,xlabel="Rm [m]", ylabel="pm", title="E: $(round(E,digits=3)) keV", fillcolor=cgrad([:white, :darkblue, :green, :yellow, :orange, :red]), ylims=extrema(pm_array), xlims=extrema(Rm_array))
        else
            nl = length(Λ_array_tor) # tor or pol should not matter
            npp = length(Pϕ_n_array_pol) # tor or pol should not matter
            if pm<0.0
                resonance = -log10.(abs.((ω*1000) .*ones(nl,npp) - n .*ω_ϕ_COM[Int64(Eci),:,:,1] - m .*ω_θ_COM[Int64(Eci),:,:,1] - (ω_rot*1000) .*ones(nl,npp)))
            else
                resonance = -log10.(abs.((ω*1000) .*ones(nl,npp) - n .*ω_ϕ_COM[Int64(Eci),:,:,2] - m .*ω_θ_COM[Int64(Eci),:,:,2] - (ω_rot*1000) .*ones(nl,npp)))
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
        if save_plots
            plt_res = Plots.plot!(dpi=600)
            if (phase_space==:OS) || !enable_COM
                png(plt_res, "plt_res_$(round(E, digits=2))_$(round(o.coordinate.pitch, digits=2))_$(round(o.coordinate.r,digits=2))")
            else
                png(plt_res, "plt_res_$(round(E, digits=2))_$(round(Λ, sigdigits=2))_$(round(Pϕ_n,sigdigits=2))")
            end
        end

        # Put all the plots together
        vbox(vskip(1em),
        hbox(Plots.plot(plt_res),Plots.plot(plt_topo)),
        hbox(Plots.plot(plt_top), Plots.plot(plt_crs))
        )
    end
end
webio_serve(page("/",app), port)