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

# ╔═╡ 2f9bc6d2-aa8f-11f1-1d73-81e4088d833b
md"""
# signalWebApp

## Description:
This notebook provides an easy way of visualizing a synthetic diagnostic signal (WF) and its contributions in terms of orbit types. All you need to specify is the output file from the OWCF ps2WF.jl tool. HOWEVER, PLEASE NOTE, you need to have run ps2WF.jl with the input variable 'calcWFOs' set to true.

The interactive web application below will also allow you to visualize the fast-ion distribution and the orbit-space sensitivity.

## Inputs:
- folderpath_OWCF - The path to the OWCF folder on your computer - String
- filepath\_ps2WF_output - The output file from the ps2WF.jl script - String
- energy\_ticks - The x-axis ticks when visualizing f(E), w(E) and w(E)f(E) - Array{Union{Float64,Int64},1}
- pm\_ticks - The x-axis ticks when visualizing f(pm), w(pm) and w(pm)f(pm) - Array{Union{Float64,Int64},1}
- Rm\_ticks - The x-axis ticks when visualizing f(Rm), w(Rm) and w(Rm)f(Rm) - Array{Union{Float64,Int64},1}
- verbose - If set to true, the app will talk a lot! - Bool

## Outputs:
# -

## Saved files:
# -

## Other:
Please note, you CANNOT multiply w(E) with f(E) and expect to get the correct w(E)f(E) (and equivalent for pm and Rm). Because all the different weights have been binned in energy, but they actually vary for pm and Rm as well. You can sum all the f(E) to get the total number of fast ions however (and equivalent for pm and Rm). And you can sum all the w(E)f(E) and get the correct WF signal for that particular diagnostic energy (and equivalent for pm and Rm).

### Notebook written by Henrik Järleblad, henrikj@dtu.dk
### Last maintained 2026-09-07
"""

# ╔═╡ 8ffe494b-f5a0-472e-a912-16d9f5009b8b
begin
	# SPECIFY THE INPUTS IN THIS CELL

	# Please specify the OWCF folder and let the notebook change directory to the 
	# OWCF folder when the cell below is executed. This is to be able to load the
	# correct versions of the Julia packages as specified in the Project.toml and 
	# Manifest.toml files.
	folderpath_OWCF = "REPLACE-THIS-TEXT-WITH-THE-PATH-TO-THE-OWCF-FOLDER-ON-YOUR-COMPUTER" # Finish with '/'
	
	filepath_ps2WF_output = folderpath_OWCF*"apps/example_data/ps2WF_results_JET_99971L72_at8,9s_AB_T-D--n-4He_126x101x104.jld2" # The path to a .jld2 file that is the output file os the ps2WF.jl script
	energy_ticks = [0, 25, 50, 75, 100, 150, 200, 250] # Example energy ticks between 0 and 175 keV: [0,10,20,30,40,50,60,80,100,125,150,175]
	pm_ticks = [-1.0, -0.5, 0.0, 0.5, 1.0] # Example pitch maximum ticks between -1.0 and 1.0: [-1.0,-0.5,0.0,0.5,1.0]
	Rm_ticks = [2.75, 3.00, 3.25, 3.50, 3.75] # Example radius maximum ticks between 2.75 and 3.75 meters (good for JET): [2.75, 3.00, 3.25, 3.50, 3.75]
	verbose = true
end

# ╔═╡ 3c464d4c-7d47-4a2c-a0df-b8ab1f19d069
begin
	# Change directory to the OWCF folder, load the Pkg package manager and activate the OWCF environment, to get access to all Julia packages
	cd(folderpath_OWCF)
	using Pkg
	Pkg.activate(".")
	
	# Loading packages
	verbose && println("Loading packages... ")
	using Plots
	using PlutoUI
	using JLD2
	using FileIO
end

# ╔═╡ 21db47b1-fec7-439e-a465-17d0c4ba1670
begin
	# Load data och pre-process for app usage
	
	# Read the .jld2-file for the required data
	verbose && println("Loading data from file... ")
	myfile = jldopen(filepath_ps2WF_output,false,false,false,IOStream)
	S_WF = myfile["S_WF"]
	Ed_array = myfile["Ed_array"]
	E_array = myfile["E_array"]
	pm_array = myfile["pm_array"]
	Rm_array = myfile["Rm_array"]
	
	WO_E = myfile["WO_E"]
	WO_pm = myfile["WO_pm"]
	WO_Rm = myfile["WO_Rm"]
	
	FO_E = myfile["FO_E"]
	FO_pm = myfile["FO_pm"]
	FO_Rm = myfile["FO_Rm"]
	if length(size(FO_E))==3
	    FO_E = FO_E[1,:,:]
	end
	if length(size(FO_pm))==3
	    FO_pm = FO_pm[1,:,:]
	end
	if length(size(FO_Rm))==3
	    FO_Rm = FO_Rm[1,:,:]
	end
	
	WFO_E = myfile["WFO_E"]
	WFO_pm = myfile["WFO_pm"]
	WFO_Rm = myfile["WFO_Rm"]
	close(myfile)
	
	# Determining diagnostic from file
	verbose && println("Determining diagnostic... ")
	diagnostic = lowercase((split(split(filepath_ps2WF_output,"/")[end],"_"))[6]) # Know that the sixth element will be the diagnostic
	if (diagnostic=="tofor") || (diagnostic=="km14")
	    sig_color = :green3
	elseif (diagnostic=="ab") || (diagnostic=="ne213") || (diagnostic=="km15")
	    sig_color = :red1
	elseif (diagnostic=="mpru")
	    sig_color = :purple
	else
	    println("Diagnostic: "*diagnostic)
	    @warn "Unknown diagnostic for loaded data. Signal color will be gray."
	    sig_color = :gray
	end
end

# ╔═╡ 005854d6-097d-4d43-907d-e07626e7dd4c
function plot_inputs()

	return PlutoUI.combine() do Child

		inputs=[
			md""" Density: $(
				Child("density", Select(["orbits", "total"]))
			)""",
			md""" Mode: $(
				Child("mode", Select(["absolute", "fraction"]))
			)""",
			md""" Y-scale: $(
				Child("y_scale", Select([:identity => "linear", :log10 => "logarithmic"]))
			)""",
			md""" Legend position: $(
				Child("legend_position", Select([:topleft => "Top-left", :topright => "Top-right", :bottomleft => "Bottom-left", :bottomright => "Bottom-right"]))
			)""",
			md""" Coordinate: $(
				Child("coord", Select(["E", "pm", "Rm"]))
			)""",
			md""" Diagnostic measurement bin (Ed): $(
				Child("Ed", Slider(Ed_array, default=Ed_array[Int64(round(length(Ed_array)/2))]))
			)"""
		]

		md"""
		#### Plot controls
		$(inputs)
		"""
	end
end

# ╔═╡ 07070349-5108-40fc-92e6-193d6c0b8f30
@bind my_plot_inputs plot_inputs()

# ╔═╡ 8d776aee-5c9d-4103-a518-a6f47d1d612d
let
    density = my_plot_inputs.density
    mode = my_plot_inputs.mode
    y_scale = my_plot_inputs.y_scale
    legend_position = my_plot_inputs.legend_position
    coord = my_plot_inputs.coord
    Ed = my_plot_inputs.Ed
    
	iEd = findfirst(x-> x==Ed, Ed_array)

    plt_S_WF = Plots.plot(Ed_array, S_WF, xlabel="Diagnostic measurement bin", ylabel="Signal counts [particles/(s*bin)]",linewidth=2.5, color=sig_color, title="WF")
    plt_S_WF = Plots.scatter!([Ed_array[iEd]],[S_WF[iEd]],mc=sig_color, legend=false, ms=5.0)

    # Compute the orbit fractions of the total signal at Ed
    stagnation_fraction = sum(WFO_E[iEd,:,1]) / sum(WFO_E[iEd,:,:])
    trapped_fraction = sum(WFO_E[iEd,:,2]) / sum(WFO_E[iEd,:,:])
    copassing_fraction = sum(WFO_E[iEd,:,3]) / sum(WFO_E[iEd,:,:])
    counterpassing_fraction = sum(WFO_E[iEd,:,4]) / sum(WFO_E[iEd,:,:])
    potato_fraction = sum(WFO_E[iEd,:,5]) / sum(WFO_E[iEd,:,:])
    counterstagnation_fraction = sum(WFO_E[iEd,:,6]) / sum(WFO_E[iEd,:,:])

    # linetype=:bar gives bar plot
    plt_orb_fracs = Plots.plot([1],[stagnation_fraction],color=:red,linetype=:bar,label="Stagnation")
    plt_orb_fracs = Plots.plot!([2],[trapped_fraction],color=:blue,linetype=:bar,label="Trapped")
    plt_orb_fracs = Plots.plot!([3],[copassing_fraction],color=:green,linetype=:bar,label="Co-passing")
    plt_orb_fracs = Plots.plot!([4],[counterpassing_fraction],color=:purple,linetype=:bar,label="Counter-passing")
    plt_orb_fracs = Plots.plot!([5],[potato_fraction],color=:orange,linetype=:bar,label="Potato")
    plt_orb_fracs = Plots.plot!([6],[counterstagnation_fraction],color=:pink,linetype=:bar,label="Counter-stagnation", xticks=false, title="Signal origin",ylabel="Fraction",ylims=(0.0, 1.0),xlabel="Orbit types")

    if coord=="E"
        if density=="orbits"
            if mode=="absolute"
                denom_WF = ones(size(WFO_E[iEd,:,1]))
                ylabel_WF = "Signal [Signal/keV]"
                title_WF = "w(E)f(E)"

                denom_F = ones(size(FO_E[:,1]))
                ylabel_F = "Distribution [Fast ions/keV]"
                title_F = "f(E)"

                denom_W = ones(size(WO_E[iEd,:,1]))
                ylabel_W = "Sensitivity [Signal/fast ion]"
                title_W = "w(E)"
            else
                denom_WF = sum(WFO_E[iEd,:,:],dims=2)
                ylabel_WF = "Signal [a.u.]"
                title_WF = "w(E)f(E) fractions"

                denom_F = sum(FO_E[:,:],dims=2)
                ylabel_F = "Distribution [a.u.]"
                title_F = "f(E) fractions"

                denom_W = sum(WO_E[iEd,:,:],dims=2)
                ylabel_W = "Sensitivity [a.u.]"
                title_W = "w(E) fractions"
            end
            if y_scale==:identity
                gi_WF_1 = 1:length(WFO_E[iEd,:,1])
                gi_WF_2 = 1:length(WFO_E[iEd,:,2])
                gi_WF_3 = 1:length(WFO_E[iEd,:,3])
                gi_WF_4 = 1:length(WFO_E[iEd,:,4])
                gi_WF_5 = 1:length(WFO_E[iEd,:,5])
                gi_WF_6 = 1:length(WFO_E[iEd,:,6])

                gi_F_1 = 1:length(FO_E[:,1])
                gi_F_2 = 1:length(FO_E[:,2])
                gi_F_3 = 1:length(FO_E[:,3])
                gi_F_4 = 1:length(FO_E[:,4])
                gi_F_5 = 1:length(FO_E[:,5])
                gi_F_6 = 1:length(FO_E[:,6])

                gi_W_1 = 1:length(WO_E[iEd,:,1])
                gi_W_2 = 1:length(WO_E[iEd,:,2])
                gi_W_3 = 1:length(WO_E[iEd,:,3])
                gi_W_4 = 1:length(WO_E[iEd,:,4])
                gi_W_5 = 1:length(WO_E[iEd,:,5])
                gi_W_6 = 1:length(WO_E[iEd,:,6])
            else
                gi_WF_1 = findall(x -> x > 0.0, WFO_E[iEd,:,1])
                gi_WF_2 = findall(x -> x > 0.0, WFO_E[iEd,:,2])
                gi_WF_3 = findall(x -> x > 0.0, WFO_E[iEd,:,3])
                gi_WF_4 = findall(x -> x > 0.0, WFO_E[iEd,:,4])
                gi_WF_5 = findall(x -> x > 0.0, WFO_E[iEd,:,5])
                gi_WF_6 = findall(x -> x > 0.0, WFO_E[iEd,:,6])

                gi_F_1 = findall(x -> x > 0.0, FO_E[:,1])
                gi_F_2 = findall(x -> x > 0.0, FO_E[:,2])
                gi_F_3 = findall(x -> x > 0.0, FO_E[:,3])
                gi_F_4 = findall(x -> x > 0.0, FO_E[:,4])
                gi_F_5 = findall(x -> x > 0.0, FO_E[:,5])
                gi_F_6 = findall(x -> x > 0.0, FO_E[:,6])

                gi_W_1 = findall(x -> x > 0.0, WO_E[iEd,:,1])
                gi_W_2 = findall(x -> x > 0.0, WO_E[iEd,:,2])
                gi_W_3 = findall(x -> x > 0.0, WO_E[iEd,:,3])
                gi_W_4 = findall(x -> x > 0.0, WO_E[iEd,:,4])
                gi_W_5 = findall(x -> x > 0.0, WO_E[iEd,:,5])
                gi_W_6 = findall(x -> x > 0.0, WO_E[iEd,:,6])
            end
            plt_WF = Plots.plot(E_array[gi_WF_1], (WFO_E[iEd,:,1] ./denom_WF)[gi_WF_1], label="Stagnation", color=:red, linewidth=2.5)
            plt_WF = Plots.plot!(E_array[gi_WF_2], (WFO_E[iEd,:,2] ./denom_WF)[gi_WF_2], label="Trapped", color=:blue, linewidth=2.5)
            plt_WF = Plots.plot!(E_array[gi_WF_3], (WFO_E[iEd,:,3] ./denom_WF)[gi_WF_3], label="Co-passing", color=:green, linewidth=2.5)
            plt_WF = Plots.plot!(E_array[gi_WF_4], (WFO_E[iEd,:,4] ./denom_WF)[gi_WF_4], label="Counter-passing", color=:purple, linewidth=2.5)
            plt_WF = Plots.plot!(E_array[gi_WF_5], (WFO_E[iEd,:,5] ./denom_WF)[gi_WF_5], label="Potato", color=:orange, linewidth=2.5)
            plt_WF = Plots.plot!(E_array[gi_WF_6], (WFO_E[iEd,:,6] ./denom_WF)[gi_WF_6], label="Counter-stagnation", color=:pink, linewidth=2.5, xlabel="Fast-ion energy [keV]", ylabel=ylabel_WF, title=title_WF)
            plt_WF = Plots.plot!(size=(1200,400),left_margin=5Plots.mm,bottom_margin=5Plots.mm)
            plt_WF = Plots.vline!([maximum(E_array)],color=:gray,linewidth=1.0,label="",linestyle=:dash)
            plt_WF = Plots.plot!(xlims=(minimum(E_array)-0.02*(maximum(E_array) - minimum(E_array)),maximum(E_array) + 0.3333*(maximum(E_array) - minimum(E_array))))
            plt_WF = Plots.plot!(xticks=energy_ticks, yscale=y_scale, legend=legend_position)

            plt_F = Plots.plot(E_array[gi_F_1], (FO_E[:,1] ./denom_F)[gi_F_1], label="Stagnation", color=:red, linewidth=2.5)
            plt_F = Plots.plot!(E_array[gi_F_2], (FO_E[:,2] ./denom_F)[gi_F_2], label="Trapped", color=:blue, linewidth=2.5)
            plt_F = Plots.plot!(E_array[gi_F_3], (FO_E[:,3] ./denom_F)[gi_F_3], label="Co-passing", color=:green, linewidth=2.5)
            plt_F = Plots.plot!(E_array[gi_F_4], (FO_E[:,4] ./denom_F)[gi_F_4], label="Counter-passing", color=:purple, linewidth=2.5)
            plt_F = Plots.plot!(E_array[gi_F_5], (FO_E[:,5] ./denom_F)[gi_F_5], label="Potato", color=:orange, linewidth=2.5)
            plt_F = Plots.plot!(E_array[gi_F_6], (FO_E[:,6] ./denom_F)[gi_F_6], label="Counter-stagnation", color=:pink, linewidth=2.5, xlabel="Fast-ion energy [keV]", ylabel=ylabel_F, title=title_F)
            plt_F = Plots.plot!(size=(1200,400),left_margin=5Plots.mm,bottom_margin=5Plots.mm)
            plt_F = Plots.vline!([maximum(E_array)],color=:gray,linewidth=1.0,label="",linestyle=:dash)
            plt_F = Plots.plot!(xlims=(minimum(E_array)-0.02*(maximum(E_array) - minimum(E_array)),maximum(E_array) + 0.3333*(maximum(E_array) - minimum(E_array))))
            plt_F = Plots.plot!(xticks=energy_ticks, yscale=y_scale, legend=legend_position)

            plt_W = Plots.plot(E_array[gi_W_1], (WO_E[iEd,:,1] ./denom_W)[gi_W_1], label="Stagnation", color=:red, linewidth=2.5)
            plt_W = Plots.plot!(E_array[gi_W_2], (WO_E[iEd,:,2] ./denom_W)[gi_W_2], label="Trapped", color=:blue, linewidth=2.5)
            plt_W = Plots.plot!(E_array[gi_W_3], (WO_E[iEd,:,3] ./denom_W)[gi_W_3], label="Co-passing", color=:green, linewidth=2.5)
            plt_W = Plots.plot!(E_array[gi_W_4], (WO_E[iEd,:,4] ./denom_W)[gi_W_4], label="Counter-passing", color=:purple, linewidth=2.5)
            plt_W = Plots.plot!(E_array[gi_W_5], (WO_E[iEd,:,5] ./denom_W)[gi_W_5], label="Potato", color=:orange, linewidth=2.5)
            plt_W = Plots.plot!(E_array[gi_W_6], (WO_E[iEd,:,6] ./denom_W)[gi_W_6], label="Counter-stagnation", color=:pink, linewidth=2.5, xlabel="Fast-ion energy [keV]", ylabel=ylabel_W, title=title_W)
            plt_W = Plots.plot!(size=(1200,400),left_margin=5Plots.mm,bottom_margin=5Plots.mm)
            plt_W = Plots.vline!([maximum(E_array)],color=:gray,linewidth=1.0,label="",linestyle=:dash)
            plt_W = Plots.plot!(xlims=(minimum(E_array)-0.02*(maximum(E_array) - minimum(E_array)),maximum(E_array) + 0.3333*(maximum(E_array) - minimum(E_array))))
            plt_W = Plots.plot!(xticks=energy_ticks, yscale=y_scale, legend=legend_position)
        else
            WF = dropdims(sum(WFO_E[iEd,:,:],dims=2), dims=2) # Sum over all orbit types
            F = dropdims(sum(FO_E[:,:],dims=2), dims=2) # Sum over all orbit types
            W = dropdims(sum(WO_E[iEd,:,:],dims=2), dims=2) # Sum over all orbit types
            if mode=="absolute"
                denom_WF = ones(size(WF))
                ylabel_WF = "Signal [Signal/keV]"
                title_WF = "w(E)f(E)"

                denom_F = ones(size(F))
                ylabel_F = "Distribution [Fast ions/keV]"
                title_F = "f(E)"

                denom_W = ones(size(W))
                ylabel_W = "Sensitivity [Signal/fast ion]"
                title_W = "w(E)"
            else
                denom_WF = maximum(WF)
                ylabel_WF = "Signal [a.u.]"
                title_WF = "w(E)f(E) normalized"

                denom_F = maximum(F)
                ylabel_F = "Distribution [a.u.]"
                title_F = "f(E) normalized"

                denom_W = maximum(W)
                ylabel_W = "Sensitivity [a.u.]"
                title_W = "w(E) normalized"
            end
            if y_scale==:identity
                gi_WF = 1:length(WF)
                gi_F = 1:length(F)
                gi_W = 1:length(W)
            else
                gi_WF = findall(x -> x>0.0,WF)
                gi_F = findall(x -> x>0.0,F)
                gi_W = findall(x -> x>0.0,W)
            end
            plt_WF = Plots.plot(E_array[gi_WF], (WF ./denom_WF)[gi_WF], label="Total", color=:black, linewidth=2.5, xlabel="Fast-ion energy [keV]", ylabel=ylabel_WF, title=title_WF)
            plt_WF = Plots.plot!(size=(1200,400),left_margin=5Plots.mm,bottom_margin=5Plots.mm)
            plt_WF = Plots.vline!([maximum(E_array)],color=:gray,linewidth=1.0,label="",linestyle=:dash)
            plt_WF = Plots.plot!(xlims=(minimum(E_array)-0.02*(maximum(E_array) - minimum(E_array)),maximum(E_array) + 0.3333*(maximum(E_array) - minimum(E_array))))
            plt_WF = Plots.plot!(xticks=energy_ticks, yscale=y_scale, legend=legend_position)

            plt_F = Plots.plot(E_array[gi_F], (F ./denom_F)[gi_F], label="Total", color=:black, linewidth=2.5, xlabel="Fast-ion energy [keV]", ylabel=ylabel_F, title=title_F)
            plt_F = Plots.plot!(size=(1200,400),left_margin=5Plots.mm,bottom_margin=5Plots.mm)
            plt_F = Plots.vline!([maximum(E_array)],color=:gray,linewidth=1.0,label="",linestyle=:dash)
            plt_F = Plots.plot!(xlims=(minimum(E_array)-0.02*(maximum(E_array) - minimum(E_array)),maximum(E_array) + 0.3333*(maximum(E_array) - minimum(E_array))))
            plt_F = Plots.plot!(xticks=energy_ticks, yscale=y_scale, legend=legend_position)

            plt_W = Plots.plot(E_array[gi_W], (W ./denom_W)[gi_W], label="Total", color=:black, linewidth=2.5, xlabel="Fast-ion energy [keV]", ylabel=ylabel_W, title=title_W)
            plt_W = Plots.plot!(size=(1200,400),left_margin=5Plots.mm,bottom_margin=5Plots.mm)
            plt_W = Plots.vline!([maximum(E_array)],color=:gray,linewidth=1.0,label="",linestyle=:dash)
            plt_W = Plots.plot!(xlims=(minimum(E_array)-0.02*(maximum(E_array) - minimum(E_array)),maximum(E_array) + 0.3333*(maximum(E_array) - minimum(E_array))))
            plt_W = Plots.plot!(xticks=energy_ticks, yscale=y_scale, legend=legend_position)
        end
    elseif coord=="pm"
        if density=="orbits"
            if mode=="absolute"
                denom_WF = ones(size(WFO_pm[iEd,:,1]))
                ylabel_WF = "Signal [Signal]"
                title_WF = "w(pm)f(pm)"

                denom_F = ones(size(FO_pm[:,1]))
                ylabel_F = "Distribution [Fast ions]"
                title_F = "f(pm)"

                denom_W = ones(size(WO_pm[iEd,:,1]))
                ylabel_W = "Sensitivity [Signal/fast ion]"
                title_W = "w(pm)"
            else
                denom_WF = sum(WFO_pm[iEd,:,:],dims=2)
                ylabel_WF = "Signal [a.u.]"
                title_WF = "w(pm)f(pm) fractions"

                denom_F = sum(FO_pm[:,:],dims=2)
                ylabel_F = "Distribution [a.u.]"
                title_F = "f(pm) fractions"

                denom_W = sum(WO_pm[iEd,:,:],dims=2)
                ylabel_W = "Sensitivity [a.u.]"
                title_W = "w(pm) fractions"
            end
            if y_scale==:identity
                gi_WF_1 = 1:length(WFO_pm[iEd,:,1])
                gi_WF_2 = 1:length(WFO_pm[iEd,:,2])
                gi_WF_3 = 1:length(WFO_pm[iEd,:,3])
                gi_WF_4 = 1:length(WFO_pm[iEd,:,4])
                gi_WF_5 = 1:length(WFO_pm[iEd,:,5])
                gi_WF_6 = 1:length(WFO_pm[iEd,:,6])

                gi_F_1 = 1:length(FO_pm[:,1])
                gi_F_2 = 1:length(FO_pm[:,2])
                gi_F_3 = 1:length(FO_pm[:,3])
                gi_F_4 = 1:length(FO_pm[:,4])
                gi_F_5 = 1:length(FO_pm[:,5])
                gi_F_6 = 1:length(FO_pm[:,6])

                gi_W_1 = 1:length(WO_pm[iEd,:,1])
                gi_W_2 = 1:length(WO_pm[iEd,:,2])
                gi_W_3 = 1:length(WO_pm[iEd,:,3])
                gi_W_4 = 1:length(WO_pm[iEd,:,4])
                gi_W_5 = 1:length(WO_pm[iEd,:,5])
                gi_W_6 = 1:length(WO_pm[iEd,:,6])
            else
                gi_WF_1 = findall(x -> x > 0.0, WFO_pm[iEd,:,1])
                gi_WF_2 = findall(x -> x > 0.0, WFO_pm[iEd,:,2])
                gi_WF_3 = findall(x -> x > 0.0, WFO_pm[iEd,:,3])
                gi_WF_4 = findall(x -> x > 0.0, WFO_pm[iEd,:,4])
                gi_WF_5 = findall(x -> x > 0.0, WFO_pm[iEd,:,5])
                gi_WF_6 = findall(x -> x > 0.0, WFO_pm[iEd,:,6])

                gi_F_1 = findall(x -> x > 0.0, FO_pm[:,1])
                gi_F_2 = findall(x -> x > 0.0, FO_pm[:,2])
                gi_F_3 = findall(x -> x > 0.0, FO_pm[:,3])
                gi_F_4 = findall(x -> x > 0.0, FO_pm[:,4])
                gi_F_5 = findall(x -> x > 0.0, FO_pm[:,5])
                gi_F_6 = findall(x -> x > 0.0, FO_pm[:,6])

                gi_W_1 = findall(x -> x > 0.0, WO_pm[iEd,:,1])
                gi_W_2 = findall(x -> x > 0.0, WO_pm[iEd,:,2])
                gi_W_3 = findall(x -> x > 0.0, WO_pm[iEd,:,3])
                gi_W_4 = findall(x -> x > 0.0, WO_pm[iEd,:,4])
                gi_W_5 = findall(x -> x > 0.0, WO_pm[iEd,:,5])
                gi_W_6 = findall(x -> x > 0.0, WO_pm[iEd,:,6])
            end
            plt_WF = Plots.plot(pm_array[gi_WF_1], (WFO_pm[iEd,:,1] ./denom_WF)[gi_WF_1], label="Stagnation", color=:red, linewidth=2.5)
            plt_WF = Plots.plot!(pm_array[gi_WF_2], (WFO_pm[iEd,:,2] ./denom_WF)[gi_WF_2], label="Trapped", color=:blue, linewidth=2.5)
            plt_WF = Plots.plot!(pm_array[gi_WF_3], (WFO_pm[iEd,:,3] ./denom_WF)[gi_WF_3], label="Co-passing", color=:green, linewidth=2.5)
            plt_WF = Plots.plot!(pm_array[gi_WF_4], (WFO_pm[iEd,:,4] ./denom_WF)[gi_WF_4], label="Counter-passing", color=:purple, linewidth=2.5)
            plt_WF = Plots.plot!(pm_array[gi_WF_5], (WFO_pm[iEd,:,5] ./denom_WF)[gi_WF_5], label="Potato", color=:orange, linewidth=2.5)
            plt_WF = Plots.plot!(pm_array[gi_WF_6], (WFO_pm[iEd,:,6] ./denom_WF)[gi_WF_6], label="Counter-stagnation", color=:pink, linewidth=2.5, xlabel="Pitch maximum [-]", ylabel=ylabel_WF, title=title_WF)
            plt_WF = Plots.plot!(size=(1200,400),left_margin=5Plots.mm,bottom_margin=5Plots.mm)
            plt_WF = Plots.vline!([maximum(pm_array)],color=:gray,linewidth=1.0,label="",linestyle=:dash)
            plt_WF = Plots.plot!(xlims=(minimum(pm_array)-0.02*(maximum(pm_array) - minimum(pm_array)),maximum(pm_array) + 0.3333*(maximum(pm_array) - minimum(pm_array))))
            plt_WF = Plots.plot!(xticks=pm_ticks, yscale=y_scale, legend=legend_position)

            plt_F = Plots.plot(pm_array[gi_F_1], (FO_pm[:,1] ./denom_F)[gi_F_1], label="Stagnation", color=:red, linewidth=2.5)
            plt_F = Plots.plot!(pm_array[gi_F_2], (FO_pm[:,2] ./denom_F)[gi_F_2], label="Trapped", color=:blue, linewidth=2.5)
            plt_F = Plots.plot!(pm_array[gi_F_3], (FO_pm[:,3] ./denom_F)[gi_F_3], label="Co-passing", color=:green, linewidth=2.5)
            plt_F = Plots.plot!(pm_array[gi_F_4], (FO_pm[:,4] ./denom_F)[gi_F_4], label="Counter-passing", color=:purple, linewidth=2.5)
            plt_F = Plots.plot!(pm_array[gi_F_5], (FO_pm[:,5] ./denom_F)[gi_F_5], label="Potato", color=:orange, linewidth=2.5)
            plt_F = Plots.plot!(pm_array[gi_F_6], (FO_pm[:,6] ./denom_F)[gi_F_6], label="Counter-stagnation", color=:pink, linewidth=2.5, xlabel="Pitch maximum [-]", ylabel=ylabel_F, title=title_F)
            plt_F = Plots.plot!(size=(1200,400),left_margin=5Plots.mm,bottom_margin=5Plots.mm)
            plt_F = Plots.vline!([maximum(pm_array)],color=:gray,linewidth=1.0,label="",linestyle=:dash)
            plt_F = Plots.plot!(xlims=(minimum(pm_array)-0.02*(maximum(pm_array) - minimum(pm_array)),maximum(pm_array) + 0.3333*(maximum(pm_array) - minimum(pm_array))))
            plt_F = Plots.plot!(xticks=pm_ticks, yscale=y_scale, legend=legend_position)

            plt_W = Plots.plot(pm_array[gi_W_1], (WO_pm[iEd,:,1] ./denom_W)[gi_W_1], label="Stagnation", color=:red, linewidth=2.5)
            plt_W = Plots.plot!(pm_array[gi_W_2], (WO_pm[iEd,:,2] ./denom_W)[gi_W_2], label="Trapped", color=:blue, linewidth=2.5)
            plt_W = Plots.plot!(pm_array[gi_W_3], (WO_pm[iEd,:,3] ./denom_W)[gi_W_3], label="Co-passing", color=:green, linewidth=2.5)
            plt_W = Plots.plot!(pm_array[gi_W_4], (WO_pm[iEd,:,4] ./denom_W)[gi_W_4], label="Counter-passing", color=:purple, linewidth=2.5)
            plt_W = Plots.plot!(pm_array[gi_W_5], (WO_pm[iEd,:,5] ./denom_W)[gi_W_5], label="Potato", color=:orange, linewidth=2.5)
            plt_W = Plots.plot!(pm_array[gi_W_6], (WO_pm[iEd,:,6] ./denom_W)[gi_W_6], label="Counter-stagnation", color=:pink, linewidth=2.5, xlabel="Pitch maximum [-]", ylabel=ylabel_W, title=title_W)
            plt_W = Plots.plot!(size=(1200,400),left_margin=5Plots.mm,bottom_margin=5Plots.mm)
            plt_W = Plots.vline!([maximum(pm_array)],color=:gray,linewidth=1.0,label="",linestyle=:dash)
            plt_W = Plots.plot!(xlims=(minimum(pm_array)-0.02*(maximum(pm_array) - minimum(pm_array)),maximum(pm_array) + 0.3333*(maximum(pm_array) - minimum(pm_array))))
            plt_W = Plots.plot!(xticks=pm_ticks, yscale=y_scale, legend=legend_position)
        else
            WF = dropdims(sum(WFO_pm[iEd,:,:],dims=2), dims=2) # Sum over all orbit types
            F = dropdims(sum(FO_pm[:,:],dims=2), dims=2) # Sum over all orbit types
            W = dropdims(sum(WO_pm[iEd,:,:],dims=2), dims=2) # Sum over all orbit types
            if mode=="absolute"
                denom_WF = ones(size(WF))
                ylabel_WF = "Signal [Signal]"
                title_WF = "w(pm)f(pm)"

                denom_F = ones(size(F))
                ylabel_F = "Distribution [Fast ions]"
                title_F = "f(pm)"

                denom_W = ones(size(W))
                ylabel_W = "Sensitivity [Signal/fast ion]"
                title_W = "w(pm)"
            else
                denom_WF = maximum(WF)
                ylabel_WF = "Signal [a.u.]"
                title_WF = "w(pm)f(pm) normalized"

                denom_F = maximum(F)
                ylabel_F = "Distribution [a.u.]"
                title_F = "f(pm) normalized"

                denom_W = maximum(W)
                ylabel_W = "Sensitivity [a.u.]"
                title_W = "w(pm) normalized"
            end
            if y_scale==:identity
                gi_WF = 1:length(WF)
                gi_F = 1:length(F)
                gi_W = 1:length(W)
            else
                gi_WF = findall(x -> x>0.0,WF)
                gi_F = findall(x -> x>0.0,F)
                gi_W = findall(x -> x>0.0,W)
            end
            plt_WF = Plots.plot(pm_array[gi_WF], (WF ./denom_WF)[gi_WF], label="Total", color=:black, linewidth=2.5, xlabel="Pitch maximum [-]", ylabel=ylabel_WF, title=title_WF)
            plt_WF = Plots.plot!(size=(1200,400),left_margin=5Plots.mm,bottom_margin=5Plots.mm)
            plt_WF = Plots.vline!([maximum(pm_array)],color=:gray,linewidth=1.0,label="",linestyle=:dash)
            plt_WF = Plots.plot!(xlims=(minimum(pm_array)-0.02*(maximum(pm_array) - minimum(pm_array)),maximum(pm_array) + 0.3333*(maximum(pm_array) - minimum(pm_array))))
            plt_WF = Plots.plot!(xticks=pm_ticks, yscale=y_scale, legend=legend_position)

            plt_F = Plots.plot(pm_array[gi_F], (F ./denom_F)[gi_F], label="Total", color=:black, linewidth=2.5, xlabel="Pitch maximum [-]", ylabel=ylabel_F, title=title_F)
            plt_F = Plots.plot!(size=(1200,400),left_margin=5Plots.mm,bottom_margin=5Plots.mm)
            plt_F = Plots.vline!([maximum(pm_array)],color=:gray,linewidth=1.0,label="",linestyle=:dash)
            plt_F = Plots.plot!(xlims=(minimum(pm_array)-0.02*(maximum(pm_array) - minimum(pm_array)),maximum(pm_array) + 0.3333*(maximum(pm_array) - minimum(pm_array))))
            plt_F = Plots.plot!(xticks=pm_ticks, yscale=y_scale, legend=legend_position)

            plt_W = Plots.plot(pm_array[gi_W], (W ./denom_W)[gi_W], label="Total", color=:black, linewidth=2.5, xlabel="Pitch maximum [-]", ylabel=ylabel_W, title=title_W)
            plt_W = Plots.plot!(size=(1200,400),left_margin=5Plots.mm,bottom_margin=5Plots.mm)
            plt_W = Plots.vline!([maximum(pm_array)],color=:gray,linewidth=1.0,label="",linestyle=:dash)
            plt_W = Plots.plot!(xlims=(minimum(pm_array)-0.02*(maximum(pm_array) - minimum(pm_array)),maximum(pm_array) + 0.3333*(maximum(pm_array) - minimum(pm_array))))
            plt_W = Plots.plot!(xticks=pm_ticks, yscale=y_scale, legend=legend_position)
        end
    else # Else, it's Rm!
        if density=="orbits"
            if mode=="absolute"
                denom_WF = ones(size(WFO_Rm[iEd,:,1]))
                ylabel_WF = "Signal [Signal/m]"
                title_WF = "w(Rm)f(Rm)"

                denom_F = ones(size(FO_Rm[:,1]))
                ylabel_F = "Distribution [Fast ions/meter]"
                title_F = "f(Rm)"

                denom_W = ones(size(WO_Rm[iEd,:,1]))
                ylabel_W = "Sensitivity [Signal/fast ion]"
                title_W = "w(Rm)"
            else
                denom_WF = sum(WFO_Rm[iEd,:,:],dims=2)
                ylabel_WF = "Signal [a.u.]"
                title_WF = "w(Rm)f(Rm) fractions"

                denom_F = sum(FO_Rm[:,:],dims=2)
                ylabel_F = "Distribution [a.u.]"
                title_F = "f(Rm) fractions"

                denom_W = sum(WO_Rm[iEd,:,:],dims=2)
                ylabel_W = "Sensitivity [a.u.]"
                title_W = "w(Rm) fractions"
            end
            if y_scale==:identity
                gi_WF_1 = 1:length(WFO_Rm[iEd,:,1])
                gi_WF_2 = 1:length(WFO_Rm[iEd,:,2])
                gi_WF_3 = 1:length(WFO_Rm[iEd,:,3])
                gi_WF_4 = 1:length(WFO_Rm[iEd,:,4])
                gi_WF_5 = 1:length(WFO_Rm[iEd,:,5])
                gi_WF_6 = 1:length(WFO_Rm[iEd,:,6])

                gi_F_1 = 1:length(FO_Rm[:,1])
                gi_F_2 = 1:length(FO_Rm[:,2])
                gi_F_3 = 1:length(FO_Rm[:,3])
                gi_F_4 = 1:length(FO_Rm[:,4])
                gi_F_5 = 1:length(FO_Rm[:,5])
                gi_F_6 = 1:length(FO_Rm[:,6])

                gi_W_1 = 1:length(WO_Rm[iEd,:,1])
                gi_W_2 = 1:length(WO_Rm[iEd,:,2])
                gi_W_3 = 1:length(WO_Rm[iEd,:,3])
                gi_W_4 = 1:length(WO_Rm[iEd,:,4])
                gi_W_5 = 1:length(WO_Rm[iEd,:,5])
                gi_W_6 = 1:length(WO_Rm[iEd,:,6])
            else
                gi_WF_1 = findall(x -> x > 0.0, WFO_Rm[iEd,:,1])
                gi_WF_2 = findall(x -> x > 0.0, WFO_Rm[iEd,:,2])
                gi_WF_3 = findall(x -> x > 0.0, WFO_Rm[iEd,:,3])
                gi_WF_4 = findall(x -> x > 0.0, WFO_Rm[iEd,:,4])
                gi_WF_5 = findall(x -> x > 0.0, WFO_Rm[iEd,:,5])
                gi_WF_6 = findall(x -> x > 0.0, WFO_Rm[iEd,:,6])

                gi_F_1 = findall(x -> x > 0.0, FO_Rm[:,1])
                gi_F_2 = findall(x -> x > 0.0, FO_Rm[:,2])
                gi_F_3 = findall(x -> x > 0.0, FO_Rm[:,3])
                gi_F_4 = findall(x -> x > 0.0, FO_Rm[:,4])
                gi_F_5 = findall(x -> x > 0.0, FO_Rm[:,5])
                gi_F_6 = findall(x -> x > 0.0, FO_Rm[:,6])

                gi_W_1 = findall(x -> x > 0.0, WO_Rm[iEd,:,1])
                gi_W_2 = findall(x -> x > 0.0, WO_Rm[iEd,:,2])
                gi_W_3 = findall(x -> x > 0.0, WO_Rm[iEd,:,3])
                gi_W_4 = findall(x -> x > 0.0, WO_Rm[iEd,:,4])
                gi_W_5 = findall(x -> x > 0.0, WO_Rm[iEd,:,5])
                gi_W_6 = findall(x -> x > 0.0, WO_Rm[iEd,:,6])
            end
            plt_WF = Plots.plot(Rm_array[gi_WF_1], (WFO_Rm[iEd,:,1] ./denom_WF)[gi_WF_1], label="Stagnation", color=:red, linewidth=2.5)
            plt_WF = Plots.plot!(Rm_array[gi_WF_2], (WFO_Rm[iEd,:,2] ./denom_WF)[gi_WF_2], label="Trapped", color=:blue, linewidth=2.5)
            plt_WF = Plots.plot!(Rm_array[gi_WF_3], (WFO_Rm[iEd,:,3] ./denom_WF)[gi_WF_3], label="Co-passing", color=:green, linewidth=2.5)
            plt_WF = Plots.plot!(Rm_array[gi_WF_4], (WFO_Rm[iEd,:,4] ./denom_WF)[gi_WF_4], label="Counter-passing", color=:purple, linewidth=2.5)
            plt_WF = Plots.plot!(Rm_array[gi_WF_5], (WFO_Rm[iEd,:,5] ./denom_WF)[gi_WF_5], label="Potato", color=:orange, linewidth=2.5)
            plt_WF = Plots.plot!(Rm_array[gi_WF_6], (WFO_Rm[iEd,:,6] ./denom_WF)[gi_WF_6], label="Counter-stagnation", color=:pink, linewidth=2.5, xlabel="Radius maximum [m]", ylabel=ylabel_WF, title=title_WF)
            plt_WF = Plots.plot!(size=(1200,400),left_margin=5Plots.mm,bottom_margin=5Plots.mm)
            plt_WF = Plots.vline!([maximum(Rm_array)],color=:gray,linewidth=1.0,label="",linestyle=:dash)
            plt_WF = Plots.plot!(xlims=(minimum(Rm_array)-0.02*(maximum(Rm_array) - minimum(Rm_array)),maximum(Rm_array) + 0.3333*(maximum(Rm_array) - minimum(Rm_array))))
            plt_WF = Plots.plot!(xticks=Rm_ticks, yscale=y_scale, legend=legend_position)

            plt_F = Plots.plot(Rm_array[gi_F_1], (FO_Rm[:,1] ./denom_F)[gi_F_1], label="Stagnation", color=:red, linewidth=2.5)
            plt_F = Plots.plot!(Rm_array[gi_F_2], (FO_Rm[:,2] ./denom_F)[gi_F_2], label="Trapped", color=:blue, linewidth=2.5)
            plt_F = Plots.plot!(Rm_array[gi_F_3], (FO_Rm[:,3] ./denom_F)[gi_F_3], label="Co-passing", color=:green, linewidth=2.5)
            plt_F = Plots.plot!(Rm_array[gi_F_4], (FO_Rm[:,4] ./denom_F)[gi_F_4], label="Counter-passing", color=:purple, linewidth=2.5)
            plt_F = Plots.plot!(Rm_array[gi_F_5], (FO_Rm[:,5] ./denom_F)[gi_F_5], label="Potato", color=:orange, linewidth=2.5)
            plt_F = Plots.plot!(Rm_array[gi_F_6], (FO_Rm[:,6] ./denom_F)[gi_F_6], label="Counter-stagnation", color=:pink, linewidth=2.5, xlabel="Radius maximum [m]", ylabel=ylabel_F, title=title_F)
            plt_F = Plots.plot!(size=(1200,400),left_margin=5Plots.mm,bottom_margin=5Plots.mm)
            plt_F = Plots.vline!([maximum(Rm_array)],color=:gray,linewidth=1.0,label="",linestyle=:dash)
            plt_F = Plots.plot!(xlims=(minimum(Rm_array)-0.02*(maximum(Rm_array) - minimum(Rm_array)),maximum(Rm_array) + 0.3333*(maximum(Rm_array) - minimum(Rm_array))))
            plt_F = Plots.plot!(xticks=Rm_ticks, yscale=y_scale, legend=legend_position)

            plt_W = Plots.plot(Rm_array[gi_W_1], (WO_Rm[iEd,:,1] ./denom_W)[gi_W_1], label="Stagnation", color=:red, linewidth=2.5)
            plt_W = Plots.plot!(Rm_array[gi_W_2], (WO_Rm[iEd,:,2] ./denom_W)[gi_W_2], label="Trapped", color=:blue, linewidth=2.5)
            plt_W = Plots.plot!(Rm_array[gi_W_3], (WO_Rm[iEd,:,3] ./denom_W)[gi_W_3], label="Co-passing", color=:green, linewidth=2.5)
            plt_W = Plots.plot!(Rm_array[gi_W_4], (WO_Rm[iEd,:,4] ./denom_W)[gi_W_4], label="Counter-passing", color=:purple, linewidth=2.5)
            plt_W = Plots.plot!(Rm_array[gi_W_5], (WO_Rm[iEd,:,5] ./denom_W)[gi_W_5], label="Potato", color=:orange, linewidth=2.5)
            plt_W = Plots.plot!(Rm_array[gi_W_6], (WO_Rm[iEd,:,6] ./denom_W)[gi_W_6], label="Counter-stagnation", color=:pink, linewidth=2.5, xlabel="Radius maximum [m]", ylabel=ylabel_W, title=title_W)
            plt_W = Plots.plot!(size=(1200,400),left_margin=5Plots.mm,bottom_margin=5Plots.mm)
            plt_W = Plots.vline!([maximum(Rm_array)],color=:gray,linewidth=1.0,label="",linestyle=:dash)
            plt_W = Plots.plot!(xlims=(minimum(Rm_array)-0.02*(maximum(Rm_array) - minimum(Rm_array)),maximum(Rm_array) + 0.3333*(maximum(Rm_array) - minimum(Rm_array))))
            plt_W = Plots.plot!(xticks=Rm_ticks, yscale=y_scale, legend=legend_position)
        else
            WF = dropdims(sum(WFO_Rm[iEd,:,:],dims=2), dims=2) # Sum over all orbit types
            F = dropdims(sum(FO_Rm[:,:],dims=2), dims=2) # Sum over all orbit types
            W = dropdims(sum(WO_Rm[iEd,:,:],dims=2), dims=2) # Sum over all orbit types
            if mode=="absolute"
                denom_WF = ones(size(WF))
                ylabel_WF = "Signal [Signal/m]"
                title_WF = "w(Rm)f(Rm)"

                denom_F = ones(size(F))
                ylabel_F = "Distribution [Fast ion/meter]"
                title_F = "f(Rm)"

                denom_W = ones(size(W))
                ylabel_W = "Sensitivity [Signal/fast ion]"
                title_W = "w(Rm)"
            else
                denom_WF = maximum(WF)
                ylabel_WF = "Signal [a.u.]"
                title_WF = "w(Rm)f(Rm) normalized"

                denom_F = maximum(F)
                ylabel_F = "Distribution [a.u.]"
                title_F = "f(Rm) normalized"

                denom_W = maximum(W)
                ylabel_W = "Sensitivity [a.u.]"
                title_W = "w(Rm) normalized"
            end
            if y_scale==:identity
                gi_WF = 1:length(WF)
                gi_F = 1:length(F)
                gi_W = 1:length(W)
            else
                gi_WF = findall(x -> x>0.0,WF)
                gi_F = findall(x -> x>0.0,F)
                gi_W = findall(x -> x>0.0,W)
            end
            plt_WF = Plots.plot(Rm_array[gi_WF], (WF ./denom_WF)[gi_WF], label="Total", color=:black, linewidth=2.5, xlabel="Radius maximum [m]", ylabel=ylabel_WF, title=title_WF)
            plt_WF = Plots.plot!(size=(1200,400),left_margin=5Plots.mm,bottom_margin=5Plots.mm)
            plt_WF = Plots.vline!([maximum(Rm_array)],color=:gray,linewidth=1.0,label="",linestyle=:dash)
            plt_WF = Plots.plot!(xlims=(minimum(Rm_array)-0.02*(maximum(Rm_array) - minimum(Rm_array)),maximum(Rm_array) + 0.3333*(maximum(Rm_array) - minimum(Rm_array))))
            plt_WF = Plots.plot!(xticks=Rm_ticks, yscale=y_scale, legend=legend_position)
            
            plt_F = Plots.plot(Rm_array[gi_F], (F ./denom_F)[gi_F], label="Total", color=:black, linewidth=2.5, xlabel="Radius maximum [m]", ylabel=ylabel_F, title=title_F)
            plt_F = Plots.plot!(size=(1200,400),left_margin=5Plots.mm,bottom_margin=5Plots.mm)
            plt_F = Plots.vline!([maximum(Rm_array)],color=:gray,linewidth=1.0,label="",linestyle=:dash)
            plt_F = Plots.plot!(xlims=(minimum(Rm_array)-0.02*(maximum(Rm_array) - minimum(Rm_array)),maximum(Rm_array) + 0.3333*(maximum(Rm_array) - minimum(Rm_array))))
            plt_F = Plots.plot!(xticks=Rm_ticks, yscale=y_scale, legend=legend_position)

            plt_W = Plots.plot(Rm_array[gi_W], (W ./denom_W)[gi_W], label="Total", color=:black, linewidth=2.5, xlabel="Radius maximum [m]", ylabel=ylabel_W, title=title_W)
            plt_W = Plots.plot!(size=(1200,400),left_margin=5Plots.mm,bottom_margin=5Plots.mm)
            plt_W = Plots.vline!([maximum(Rm_array)],color=:gray,linewidth=1.0,label="",linestyle=:dash)
            plt_W = Plots.plot!(xlims=(minimum(Rm_array)-0.02*(maximum(Rm_array) - minimum(Rm_array)),maximum(Rm_array) + 0.3333*(maximum(Rm_array) - minimum(Rm_array))))
            plt_W = Plots.plot!(xticks=Rm_ticks, yscale=y_scale, legend=legend_position)
        end
    end
    plt_dummy = Plots.plot(axis=([], false)) # Dummy plot
    plt_upper = Plots.plot(plt_S_WF, plt_orb_fracs, plt_dummy, layout=(1,3),size=(1200,400), left_margin=4Plots.mm)
    plt_lower = Plots.plot(plt_WF, plt_F, plt_W, layout=(3,1), size=(1200,1200), left_margin=4Plots.mm)
    Plots.plot(plt_upper, plt_lower, layout = grid(2, 1, heights=[0.25 , 0.75]), left_margin=4Plots.mm)
end

# ╔═╡ Cell order:
# ╠═2f9bc6d2-aa8f-11f1-1d73-81e4088d833b
# ╠═8ffe494b-f5a0-472e-a912-16d9f5009b8b
# ╠═3c464d4c-7d47-4a2c-a0df-b8ab1f19d069
# ╠═21db47b1-fec7-439e-a465-17d0c4ba1670
# ╠═005854d6-097d-4d43-907d-e07626e7dd4c
# ╠═07070349-5108-40fc-92e6-193d6c0b8f30
# ╠═8d776aee-5c9d-4103-a518-a6f47d1d612d
