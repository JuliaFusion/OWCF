### A Pluto.jl notebook ###
# v1.0.3

using Markdown
using InteractiveUtils

# ╔═╡ b910b53a-a83c-11f1-1c4d-290ec76846d3
using Pkg; Pkg.activate("."); using JLD2; using Interpolations; using Plots; using Equilibrium; include(".../OWCF/extra/dependencies.jl")

# ╔═╡ 29c89906-87d4-463f-924b-7df6bd2ee15c
begin
	folderpath_OWCF = "REPLACE-THIS-TEXT-WITH-THE-PATH-TO-THE-OWCF-FOLDER-ON-YOUR-COMPUTER"
myfile = jldopen(folderpath_OWCF*"apps/example_data/topoMap_JET_99971L72_at48,9s_D_6x101x102_wLost.jld2", false, false, false, IOStream)
	q = Dict()
	for key in keys(myfile)
		q[key] = myfile[key]
	end
	close(myfile)
end

# ╔═╡ 98ff9a87-b7d0-4c65-9eb9-49b2bdc04bc8
q

# ╔═╡ b045305a-13ce-4e7f-8354-0e102b04444e
M, wall = read_geqdsk(folderpath_OWCF*"equilibrium/JET/g99971/g99971_474-48.9.eqdsk", clockwise_phi=false)

# ╔═╡ 5b19fb38-70aa-4095-ac2f-ce4f3a42ec11
begin
	npm = length(q["pm_array"])
	nRm = length(q["Rm_array"])
	topoMap_COM, E_array_COM, Λ_array, Pϕ_n_array = os2COM(M, q["topoMap"], q["E_array"], q["pm_array"], q["Rm_array"], "D"; nl=2*npm, npp=2*nRm, isTopoMap=true, verbose=true, wall=wall)
end

# ╔═╡ 5cb14a84-a43a-4d91-8f66-3d90270c8fb9
let
	E = 100.0 # keV
	iE = argmin(abs.(q["E_array"] .- E)) # Find the closest value to E in E_array
	E = q["E_array"][iE]
	myplt = Plots.heatmap(q["Rm_array"],q["pm_array"],q["topoMap"][iE,:,:],color=:Set1_9,legend=false,xlabel="Rm [m]", ylabel="pm", title="E: $(round(E,digits=3)) keV")
	myplt
end

# ╔═╡ fb55786a-4ab4-4c7c-ba16-ac92054bea62
let
	E = 100.0 # keV
	iE = argmin(abs.(q["E_array"] .- E)) # Find the closest value to E in E_array
	E = q["E_array"][iE]
	myplt = Plots.heatmap(Pϕ_n_array, Λ_array, topoMap_COM[iE,:,:,2],color=:Set1_9,legend=false,xlabel="Pϕ_n [-]", ylabel="Λ [-]", title="E: $(round(E,digits=3)) keV")
	myplt
end

# ╔═╡ 2842c577-bd0d-4bf8-b214-8093d8f7216a
let
	E_array = q["E_array"]
	Λ_array_ext = vcat(Λ_array[1]-(diff(Λ_array))[1],Λ_array) # Extend Λ_array one row below
	
	E = 100.0 # keV
	iE = argmin(abs.(q["E_array"] .- E)) # Find the closest value to E in E_array
	topoMap_COM_ext = vcat(vcat([1,2,3,4,5,6,7,8,9],ones(length(Pϕ_n_array)-9))', topoMap_COM[Int64(iE),:,:,1]) # Extend topoMap_COM one row below, to ensure correct colormapping for orbit types (please see calcTopoMap.jl for more info). The y-limits (ylims) will make sure the extra row is not visible in the plot
    plt_topo = Plots.heatmap(Pϕ_n_array, Λ_array_ext, topoMap_COM_ext, color=:Set1_9, legend=false, xlabel="Pϕ_n", ylabel="Λ", title="E: $(round(E_array[Int64(iE)],digits=3)) keV", ylims=extrema(Λ_array), xlims=extrema(Pϕ_n_array))
end

# ╔═╡ b81d816e-6cd2-4cd1-8362-43c7ae78c5ad
let
    polTransTimes = q["polTransTimes"]
    Rm_array = q["Rm_array"]
    pm_array = q["pm_array"]
    E_array = q["E_array"]
    pm = 0.7
    E = 100.0 # keV
    Eci = argmin(abs.(q["E_array"] .- E))
    
	pTT_microsecs = polTransTimes[Int64(Eci),:,:] ./(1.0e-6) # Convert from seconds to microseconds
    nz_coords = findall(x-> x>0.0,pTT_microsecs) # Find the 2D matrix coordinates of all non-zero elements
    my_coords = length(nz_coords) > 1 ? nz_coords : CartesianIndices(size(pTT_microsecs)) # Are there actually more than one non-zero element? If not, use all elements
    min_pol, max_pol = extrema(pTT_microsecs[my_coords]) # Find minimum and maximum values
    min_OOM, max_OOM = (floor(log10(min_pol)),ceil(log10(max_pol))) # The orders of magnitude of the minimum and maximum values
    println(min_OOM)
    println(max_OOM)
    if !((max_OOM-min_OOM)==0.0) && (length(nz_coords) > 1) && ((max_pol/min_pol) > 10) # If all values are NOT within same order of magnitude AND more than one non-zero element, use logarithmic colorbar
        plt_pol = Plots.heatmap(Rm_array, pm_array, pTT_microsecs, xlabel="Rm [m]", ylabel="pm", title="tau_pol(Rm,pm) [microseconds]", fillcolor=cgrad([:white, :darkblue, :green, :yellow, :orange, :red]), colorbar=true, colorbar_scale=:log10, clims = (10^min_OOM, 10^max_OOM), top_margin=3Plots.mm) # Get nice powers-of-ten limits for the colorbar
    else
        plt_pol = Plots.heatmap(Rm_array, pm_array, pTT_microsecs, xlabel="Rm [m]", ylabel="pm", title="tau_pol(Rm,pm) [microseconds]", fillcolor=cgrad([:white, :darkblue, :green, :yellow, :orange, :red]), colorbar=true, top_margin=3Plots.mm)
    end
    plt_pol
end

# ╔═╡ c0227554-0103-4b21-94e1-4285c068eab6
begin
	good_coords = findall(x-> x!=9.0 && x!=7.0, q["topoMap"])
	println(length(good_coords))
end

# ╔═╡ 900d2f8a-4a80-40d9-a903-b3e439559d29
begin
	polTransTimes_COM, _, _, _ = os2COM(M, q["polTransTimes"], q["E_array"], q["pm_array"], q["Rm_array"], "D"; nl=2*npm, npp=2*nRm, isTopoMap=false, good_coords=good_coords, verbose=true, wall=wall)
end

# ╔═╡ aee17650-1da7-4487-b124-182730164831
let
    pm = 0.7
    E = 100.0 # keV
    Eci = argmin(abs.(q["E_array"] .- E))
    
	if pm<0.0
        pTT_microsecs = polTransTimes_COM[Int64(Eci),:,:,1] ./(1.0e-6)
        #tTT_microsecs = torTransTimes_COM[Int64(Eci),:,:,1] ./(1.0e-6)
    else
        pTT_microsecs = polTransTimes_COM[Int64(Eci),:,:,2] ./(1.0e-6)
        #tTT_microsecs = torTransTimes_COM[Int64(Eci),:,:,2] ./(1.0e-6)
    end
    nz_coords = findall(x-> x>0.0,pTT_microsecs) # Find the 2D matrix coordinates of all non-zero elements
    my_coords = length(nz_coords) > 1 ? nz_coords : CartesianIndices(size(pTT_microsecs)) # Are there actually more than one non-zero element? If not, use all elements
    min_pol, max_pol = extrema(pTT_microsecs[my_coords]) # Find minimum and maximum values
    min_OOM, max_OOM = (floor(log10(min_pol)),ceil(log10(max_pol))) # The orders of magnitude of the minimum and maximum values
    min_OOM = 0.0
    max_OOM = 3.0
    if !((max_OOM-min_OOM)==0.0) && (length(nz_coords) > 1) && ((max_pol/min_pol) > 10) # If all values are NOT within same order of magnitude AND more than one non-zero element, use logarithmic colorbar
        plt_pol = Plots.heatmap(Pϕ_n_array,Λ_array, pTT_microsecs,xlabel="Pϕ_n", ylabel="Λ", title="tau_pol(Pϕ_n,Λ) [microseconds]", fillcolor=cgrad([:white, :darkblue, :green, :yellow, :orange, :red]), colorbar=true, colorbar_scale=:log10, clims = (10^min_OOM, 10^max_OOM), top_margin=3Plots.mm, ylims=extrema(Λ_array), xlims=extrema(Pϕ_n_array)) # Get nice powers-of-ten limits for the colorbar
    else
        plt_pol = Plots.heatmap(Pϕ_n_array,Λ_array, pTT_microsecs,xlabel="Pϕ_n", ylabel="Λ", title="tau_pol(Pϕ_n,Λ) [microseconds]", fillcolor=cgrad([:white, :darkblue, :green, :yellow, :orange, :red]), colorbar=true, top_margin=3Plots.mm, ylims=extrema(Λ_array), xlims=extrema(Pϕ_n_array))
    end
    plt_pol
end

# ╔═╡ e8cf668b-0173-4e77-9d19-27635cc44a2a
let
    pm = -0.7
    E = 100.0 # keV
    Eci = argmin(abs.(q["E_array"] .- E))
    
	if pm<0.0
        pTT_microsecs = polTransTimes_COM[Int64(Eci),:,:,1] ./(1.0e-6)
        #tTT_microsecs = torTransTimes_COM[Int64(Eci),:,:,1] ./(1.0e-6)
    else
        pTT_microsecs = polTransTimes_COM[Int64(Eci),:,:,2] ./(1.0e-6)
        #tTT_microsecs = torTransTimes_COM[Int64(Eci),:,:,2] ./(1.0e-6)
    end
    nz_coords = findall(x-> x>0.0,pTT_microsecs) # Find the 2D matrix coordinates of all non-zero elements
    my_coords = length(nz_coords) > 1 ? nz_coords : CartesianIndices(size(pTT_microsecs)) # Are there actually more than one non-zero element? If not, use all elements
    min_pol, max_pol = extrema(pTT_microsecs[my_coords]) # Find minimum and maximum values
    min_OOM, max_OOM = (floor(log10(min_pol)),ceil(log10(max_pol))) # The orders of magnitude of the minimum and maximum values
    min_OOM = 0.0
    max_OOM = 3.0
    if !((max_OOM-min_OOM)==0.0) && (length(nz_coords) > 1) && ((max_pol/min_pol) > 10) # If all values are NOT within same order of magnitude AND more than one non-zero element, use logarithmic colorbar
        plt_pol = Plots.heatmap(Pϕ_n_array,Λ_array, pTT_microsecs,xlabel="Pϕ_n", ylabel="Λ", title="tau_pol(Pϕ_n,Λ) [microseconds]", fillcolor=cgrad([:white, :darkblue, :green, :yellow, :orange, :red]), colorbar=true, colorbar_scale=:log10, clims = (10^min_OOM, 10^max_OOM), top_margin=3Plots.mm, ylims=extrema(Λ_array), xlims=extrema(Pϕ_n_array)) # Get nice powers-of-ten limits for the colorbar
    else
        plt_pol = Plots.heatmap(Pϕ_n_array,Λ_array, pTT_microsecs,xlabel="Pϕ_n", ylabel="Λ", title="tau_pol(Pϕ_n,Λ) [microseconds]", fillcolor=cgrad([:white, :darkblue, :green, :yellow, :orange, :red]), colorbar=true, top_margin=3Plots.mm, ylims=extrema(Λ_array), xlims=extrema(Pϕ_n_array))
    end
    plt_pol
end

# ╔═╡ 65861658-a9b9-4efd-9827-d41a0ff510ee
let
    torTransTimes = q["torTransTimes"]
    Rm_array = q["Rm_array"]
    pm_array = q["pm_array"]
    E_array = q["E_array"]
    pm = 0.7
    E = 100.0 # keV
    Eci = argmin(abs.(q["E_array"] .- E))
    
	tTT_microsecs = torTransTimes[Int64(Eci),:,:] ./(1.0e-6) # Convert from seconds to microseconds
    nz_coords = findall(x-> x>0.0,tTT_microsecs) # Find the 2D matrix coordinates of all non-zero elements
    my_coords = length(nz_coords) > 1 ? nz_coords : CartesianIndices(size(tTT_microsecs)) # Are there actually more than one non-zero element? If not, use all elements
    min_pol, max_pol = extrema(tTT_microsecs[my_coords]) # Find minimum and maximum values
    min_OOM, max_OOM = (floor(log10(min_pol)),ceil(log10(max_pol))) # The orders of magnitude of the minimum and maximum values
    println(min_OOM)
    println(max_OOM)
    if !((max_OOM-min_OOM)==0.0) && (length(nz_coords) > 1) && ((max_pol/min_pol) > 10) # If all values are NOT within same order of magnitude AND more than one non-zero element, use logarithmic colorbar
        plt_tor = Plots.heatmap(Rm_array, pm_array, tTT_microsecs, xlabel="Rm [m]", ylabel="pm", title="tau_tor(Rm,pm) [microseconds]", fillcolor=cgrad([:white, :darkblue, :green, :yellow, :orange, :red]), colorbar=true, colorbar_scale=:log10, clims = (10^min_OOM, 10^max_OOM), top_margin=3Plots.mm) # Get nice powers-of-ten limits for the colorbar
    else
        plt_tor = Plots.heatmap(Rm_array, pm_array, tTT_microsecs, xlabel="Rm [m]", ylabel="pm", title="tau_tor(Rm,pm) [microseconds]", fillcolor=cgrad([:white, :darkblue, :green, :yellow, :orange, :red]), colorbar=true, top_margin=3Plots.mm)
    end
    plt_tor
end

# ╔═╡ 726367b7-c7a5-4e04-be19-9d79e03e42d3
begin
	torTransTimes_COM, _, _, _ = os2COM(M, q["torTransTimes"], q["E_array"], q["pm_array"], q["Rm_array"], "D"; nl=2*npm, npp=2*nRm, isTopoMap=false, good_coords=good_coords, verbose=true, wall=wall)
end

# ╔═╡ 38aca531-a5be-49c0-b54a-35f7683044a4
let
    pm = 0.7
    E = 100.0 # keV
    Eci = argmin(abs.(q["E_array"] .- E))
    
	if pm<0.0
        #pTT_microsecs = polTransTimes_COM[Int64(Eci),:,:,1] ./(1.0e-6)
        tTT_microsecs = torTransTimes_COM[Int64(Eci),:,:,1] ./(1.0e-6)
    else
        #pTT_microsecs = polTransTimes_COM[Int64(Eci),:,:,2] ./(1.0e-6)
        tTT_microsecs = torTransTimes_COM[Int64(Eci),:,:,2] ./(1.0e-6)
    end
    nz_coords = findall(x-> x>0.0,tTT_microsecs) # Find the 2D matrix coordinates of all non-zero elements
    my_coords = length(nz_coords) > 1 ? nz_coords : CartesianIndices(size(tTT_microsecs)) # Are there actually more than one non-zero element? If not, use all elements
    min_pol, max_pol = extrema(tTT_microsecs[my_coords]) # Find minimum and maximum values
    min_OOM, max_OOM = (floor(log10(min_pol)),ceil(log10(max_pol))) # The orders of magnitude of the minimum and maximum values
    min_OOM = 0.0
    max_OOM = 5.0
    if !((max_OOM-min_OOM)==0.0) && (length(nz_coords) > 1) && ((max_pol/min_pol) > 10) # If all values are NOT within same order of magnitude AND more than one non-zero element, use logarithmic colorbar
        plt_tor = Plots.heatmap(Pϕ_n_array,Λ_array, tTT_microsecs,xlabel="Pϕ_n", ylabel="Λ", title="tau_tor(Pϕ_n,Λ) [microseconds]", fillcolor=cgrad([:white, :darkblue, :green, :yellow, :orange, :red]), colorbar=true, colorbar_scale=:log10, clims = (10^min_OOM, 10^max_OOM), top_margin=3Plots.mm, ylims=extrema(Λ_array), xlims=extrema(Pϕ_n_array)) # Get nice powers-of-ten limits for the colorbar
    else
        plt_tor = Plots.heatmap(Pϕ_n_array,Λ_array, tTT_microsecs,xlabel="Pϕ_n", ylabel="Λ", title="tau_tor(Pϕ_n,Λ) [microseconds]", fillcolor=cgrad([:white, :darkblue, :green, :yellow, :orange, :red]), colorbar=true, top_margin=3Plots.mm, ylims=extrema(Λ_array), xlims=extrema(Pϕ_n_array))
    end
    plt_tor
end

# ╔═╡ 80efb786-e28b-4b6e-9285-ae23cc1b71ce
let
    pm = -0.7
    E = 100.0 # keV
    Eci = argmin(abs.(q["E_array"] .- E))
    
	if pm<0.0
        #pTT_microsecs = polTransTimes_COM[Int64(Eci),:,:,1] ./(1.0e-6)
        tTT_microsecs = torTransTimes_COM[Int64(Eci),:,:,1] ./(1.0e-6)
    else
        #pTT_microsecs = polTransTimes_COM[Int64(Eci),:,:,2] ./(1.0e-6)
        tTT_microsecs = torTransTimes_COM[Int64(Eci),:,:,2] ./(1.0e-6)
    end
    nz_coords = findall(x-> x>0.0,tTT_microsecs) # Find the 2D matrix coordinates of all non-zero elements
    my_coords = length(nz_coords) > 1 ? nz_coords : CartesianIndices(size(tTT_microsecs)) # Are there actually more than one non-zero element? If not, use all elements
    min_pol, max_pol = extrema(tTT_microsecs[my_coords]) # Find minimum and maximum values
    min_OOM, max_OOM = (floor(log10(min_pol)),ceil(log10(max_pol))) # The orders of magnitude of the minimum and maximum values
    min_OOM = 0.0
    max_OOM = 5.0
    if !((max_OOM-min_OOM)==0.0) && (length(nz_coords) > 1) && ((max_pol/min_pol) > 10) # If all values are NOT within same order of magnitude AND more than one non-zero element, use logarithmic colorbar
        plt_tor = Plots.heatmap(Pϕ_n_array,Λ_array, tTT_microsecs,xlabel="Pϕ_n", ylabel="Λ", title="tau_tor(Pϕ_n,Λ) [microseconds]", fillcolor=cgrad([:white, :darkblue, :green, :yellow, :orange, :red]), colorbar=true, colorbar_scale=:log10, clims = (10^min_OOM, 10^max_OOM), top_margin=3Plots.mm, ylims=extrema(Λ_array), xlims=extrema(Pϕ_n_array)) # Get nice powers-of-ten limits for the colorbar
    else
        plt_tor = Plots.heatmap(Pϕ_n_array,Λ_array, tTT_microsecs,xlabel="Pϕ_n", ylabel="Λ", title="tau_tor(Pϕ_n,Λ) [microseconds]", fillcolor=cgrad([:white, :darkblue, :green, :yellow, :orange, :red]), colorbar=true, top_margin=3Plots.mm, ylims=extrema(Λ_array), xlims=extrema(Pϕ_n_array))
    end
    plt_tor
end

# ╔═╡ f98ae570-5e51-4939-9e09-5a931d62fc82
begin
	myfile_f = jldopen(folderpath_OWCF*"apps/example_data/F_os_3D_JET_99971L72_at48,9s_D_6x101x102.jld2", false, false, false, IOStream)
	q_f = Dict()
	for key in keys(myfile_f)
		q_f[key] = myfile_f[key]
	end
	close(myfile_f)
end

# ╔═╡ 5a36bf1a-7bea-4d67-9c16-fed7f4a0cf76
q_f

# ╔═╡ 4e520b7a-9972-4614-9a32-7435a1bf7973
begin
	F_COM, _, _, _ = os2COM(M, q_f["F_os_3D"], q["E_array"], q["pm_array"], q["Rm_array"], "D"; nl=2*npm, npp=2*nRm, isTopoMap=false, verbose=true, good_coords=good_coords, needJac=true, wall=wall)
end

# ╔═╡ 381c45c0-2bd3-4365-8669-ceaeecff1166
let
    F_os_3D = q_f["F_os_3D"]
    Rm_array = q["Rm_array"]
    pm_array = q["pm_array"]
    E_array = q["E_array"]
    pm = 0.7
    E = 100.0 # keV
    Eci = argmin(abs.(q["E_array"] .- E))
    
	F_os_2D = F_os_3D[Int64(Eci),:,:] 
    nz_coords = findall(x-> x>0.0,F_os_2D) # Find the 2D matrix coordinates of all non-zero elements
    my_coords = length(nz_coords) > 1 ? nz_coords : CartesianIndices(size(F_os_2D)) # Are there actually more than one non-zero element? If not, use all elements
    min_pol, max_pol = extrema(F_os_2D[my_coords]) # Find minimum and maximum values
    min_OOM, max_OOM = (floor(log10(min_pol)),ceil(log10(max_pol))) # The orders of magnitude of the minimum and maximum values
    println(maximum(F_os_2D))
    if !((max_OOM-min_OOM)==0.0) && (length(nz_coords) > 1) && ((max_pol/min_pol) > 10) && false # If all values are NOT within same order of magnitude AND more than one non-zero element, use logarithmic colorbar
        plt_f = Plots.heatmap(Rm_array, pm_array, F_os_2D, xlabel="Rm [m]", ylabel="pm", title="f(Rm,pm) at $(q["E_array"][Eci]) keV", fillcolor=cgrad([:white, :darkblue, :green, :yellow, :orange, :red]), colorbar=true, colorbar_scale=:log10, clims = (10^min_OOM, 10^max_OOM), top_margin=3Plots.mm) # Get nice powers-of-ten limits for the colorbar
    else
        plt_f = Plots.heatmap(Rm_array, pm_array, F_os_2D, xlabel="Rm [m]", ylabel="pm", title="f(Rm,pm) at $(q["E_array"][Eci]) keV", fillcolor=cgrad([:white, :darkblue, :green, :yellow, :orange, :red]), colorbar=true, top_margin=3Plots.mm, right_margin=6Plots.mm, clims=(0.0, maximum(F_os_2D)))
    end
    plt_f
end

# ╔═╡ 4b778295-8701-40bb-9a4a-c104326f50d0
let
    pm = 0.7
    E = 100.0 # keV
    Eci = argmin(abs.(q["E_array"] .- E))
    
	if pm<0.0
        F_COM_2D = F_COM[Int64(Eci),:,:,1]
    else
        F_COM_2D = F_COM[Int64(Eci),:,:,2]
    end
    nz_coords = findall(x-> x>0.0,F_COM_2D) # Find the 2D matrix coordinates of all non-zero elements
    my_coords = length(nz_coords) > 1 ? nz_coords : CartesianIndices(size(F_COM_2D)) # Are there actually more than one non-zero element? If not, use all elements
    min_pol, max_pol = extrema(F_COM_2D[my_coords]) # Find minimum and maximum values
    min_OOM, max_OOM = (floor(log10(min_pol)),ceil(log10(max_pol))) # The orders of magnitude of the minimum and maximum values
    min_OOM = 16.0
    max_OOM = 20.0
    cmax = 2e19
    if !((max_OOM-min_OOM)==0.0) && (length(nz_coords) > 1) && ((max_pol/min_pol) > 10) && false # If all values are NOT within same order of magnitude AND more than one non-zero element, use logarithmic colorbar
        plt_f_COM = Plots.heatmap(Pϕ_n_array,Λ_array, F_COM_2D,xlabel="Pϕ_n", ylabel="Λ", title="f(Pϕ_n,Λ) at $(q["E_array"][Eci]) keV", fillcolor=cgrad([:white, :darkblue, :green, :yellow, :orange, :red]), colorbar=true, colorbar_scale=:log10, clims = (10^min_OOM, 10^max_OOM), top_margin=3Plots.mm, ylims=extrema(Λ_array), xlims=extrema(Pϕ_n_array)) # Get nice powers-of-ten limits for the colorbar
    else
        plt_f_COM = Plots.heatmap(Pϕ_n_array,Λ_array, F_COM_2D,xlabel="Pϕ_n", ylabel="Λ", title="f(Pϕ_n,Λ) at $(q["E_array"][Eci]) keV", fillcolor=cgrad([:white, :darkblue, :green, :yellow, :orange, :red]), colorbar=true, top_margin=3Plots.mm, ylims=extrema(Λ_array), xlims=extrema(Pϕ_n_array), clims=(0.0, cmax))
    end
    Plots.plot(plt_f_COM, right_margin=6Plots.mm)
end

# ╔═╡ 6e10eed9-d1ea-439c-8fcb-bde1de85038e
let
    pm = -0.7
    E = 100.0 # keV
    Eci = argmin(abs.(q["E_array"] .- E))
    
	if pm<0.0
        F_COM_2D = F_COM[Int64(Eci),:,:,1]
    else
        F_COM_2D = F_COM[Int64(Eci),:,:,2]
    end
    nz_coords = findall(x-> x>0.0,F_COM_2D) # Find the 2D matrix coordinates of all non-zero elements
    my_coords = length(nz_coords) > 1 ? nz_coords : CartesianIndices(size(F_COM_2D)) # Are there actually more than one non-zero element? If not, use all elements
    min_pol, max_pol = extrema(F_COM_2D[my_coords]) # Find minimum and maximum values
    min_OOM, max_OOM = (floor(log10(min_pol)),ceil(log10(max_pol))) # The orders of magnitude of the minimum and maximum values
    min_OOM = 16.0
    max_OOM = 20.0
    cmax = 2e19
    if !((max_OOM-min_OOM)==0.0) && (length(nz_coords) > 1) && ((max_pol/min_pol) > 10) && false # If all values are NOT within same order of magnitude AND more than one non-zero element, use logarithmic colorbar
        plt_f_COM = Plots.heatmap(Pϕ_n_array,Λ_array, F_COM_2D,xlabel="Pϕ_n", ylabel="Λ", title="f(Pϕ_n,Λ) at $(q["E_array"][Eci]) keV", fillcolor=cgrad([:white, :darkblue, :green, :yellow, :orange, :red]), colorbar=true, colorbar_scale=:log10, clims = (10^min_OOM, 10^max_OOM), top_margin=3Plots.mm, ylims=extrema(Λ_array), xlims=extrema(Pϕ_n_array)) # Get nice powers-of-ten limits for the colorbar
    else
        plt_f_COM = Plots.heatmap(Pϕ_n_array,Λ_array, F_COM_2D,xlabel="Pϕ_n", ylabel="Λ", title="f(Pϕ_n,Λ) at $(q["E_array"][Eci]) keV", fillcolor=cgrad([:white, :darkblue, :green, :yellow, :orange, :red]), colorbar=true, top_margin=3Plots.mm, ylims=extrema(Λ_array), xlims=extrema(Pϕ_n_array), clims=(0.0, cmax))
    end
    Plots.plot(plt_f_COM, right_margin=6Plots.mm)
end

# ╔═╡ 3ea233ad-4532-4fab-9466-3083f1464895
let
	dpm = diff(q["pm_array"])[1]
	dRm = diff(q["Rm_array"])[1]
	fE_OS = dropdims((dpm*dRm) .*sum(q_f["F_os_3D"],dims=(2,3)),dims=(2,3))
	nz_coords_OS = findall(x-> x>0.0, fE_OS)

	dΛ = diff(Λ_array)[1]
	dPϕ_n = diff(Pϕ_n_array)[1]
	fE_COM = dropdims((dΛ*dPϕ_n) .*sum(F_COM,dims=(2,3,4)),dims=(2,3,4))
	nz_coords_COM = findall(x-> x>0.0, fE_COM)

	
	myplt = Plots.plot(q["E_array"][nz_coords_OS], fE_OS[nz_coords_OS], yaxis=:log10, yticks=10 .^collect(14:18), label="f(E) OS")

	myplt = Plots.plot!(myplt, q["E_array"][nz_coords_COM], fE_COM[nz_coords_COM], label="f(E) COM")

	dE = diff(q["E_array"])[1]
	N_OS = dE * sum(fE_OS)
	N_COM = dE * sum(fE_COM)

	println("N_OS: $(N_OS)")
	println("N_COM: $(N_COM)")
	ndiff = round((N_COM - N_OS)/(0.01*N_OS),sigdigits=3)
	println("Diff: $(ndiff) %")
	
	myplt
end

# ╔═╡ Cell order:
# ╠═b910b53a-a83c-11f1-1c4d-290ec76846d3
# ╠═29c89906-87d4-463f-924b-7df6bd2ee15c
# ╠═98ff9a87-b7d0-4c65-9eb9-49b2bdc04bc8
# ╠═b045305a-13ce-4e7f-8354-0e102b04444e
# ╠═5b19fb38-70aa-4095-ac2f-ce4f3a42ec11
# ╠═5cb14a84-a43a-4d91-8f66-3d90270c8fb9
# ╠═fb55786a-4ab4-4c7c-ba16-ac92054bea62
# ╠═2842c577-bd0d-4bf8-b214-8093d8f7216a
# ╠═b81d816e-6cd2-4cd1-8362-43c7ae78c5ad
# ╠═c0227554-0103-4b21-94e1-4285c068eab6
# ╠═900d2f8a-4a80-40d9-a903-b3e439559d29
# ╠═aee17650-1da7-4487-b124-182730164831
# ╠═e8cf668b-0173-4e77-9d19-27635cc44a2a
# ╠═65861658-a9b9-4efd-9827-d41a0ff510ee
# ╠═726367b7-c7a5-4e04-be19-9d79e03e42d3
# ╠═38aca531-a5be-49c0-b54a-35f7683044a4
# ╠═80efb786-e28b-4b6e-9285-ae23cc1b71ce
# ╠═f98ae570-5e51-4939-9e09-5a931d62fc82
# ╠═5a36bf1a-7bea-4d67-9c16-fed7f4a0cf76
# ╠═4e520b7a-9972-4614-9a32-7435a1bf7973
# ╠═381c45c0-2bd3-4365-8669-ceaeecff1166
# ╠═4b778295-8701-40bb-9a4a-c104326f50d0
# ╠═6e10eed9-d1ea-439c-8fcb-bde1de85038e
# ╠═3ea233ad-4532-4fab-9466-3083f1464895
