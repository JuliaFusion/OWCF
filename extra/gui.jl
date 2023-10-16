######################################### gui.jl ##################################################
# This script defines a collection of functions that can be used to easily plot quantities computed
# with the OWCF. While the web applications are interactive, these functions can be used if a simple
# plot is desired. As of this version of the OWCF, there are more plot functions planned for gui.jl.
#
# Please note, not all plotting functionalities from the web applications are included in gui.jl.
# In fact, most are unique to the web applications. You can always use the 'save_plots' toggle 
# button in the web applications to save the current figures that you are observing. Be sure to
# un-toggle the 'save_plots' button when you have saved the figures that you want. Otherwise, the 
# web application will keep saving all figures all the time.
#
# Prior to loading this collection of functions (for example via 'include(extra/gui.jl)' when standing
# in the OWCF-folder), you should have performed the usual activation of the OWCF Julia environment.
# This is done as follows
#
# folderpath_OWCF = "/path/to/the/OWCF/folder/"
# using Distributed
# @everywhere begin
#     using Pkg
#     cd(folderpath_OWCF)
#     Pkg.activate(".")
# end
# 
# This is performed in e.g. every OWCF start template file.
#
# Written by H. Järleblad. Last updated 2023-09-26.
###################################################################################################

println("Loading the Julia packages for the grahpical user interfaces (GUI) of the OWCF... ")
using Plots, FileIO, HDF5, JLD2, GuidingCenterOrbits, LaTeXStrings
include("dependencies.jl")

"""
    plotSigComp(Ed_S, S, Ed_WF, WF)
    plotSigComp(-||- ; lc_S = :black, lc_WF = :green3, lw = 2.5, xlabel = "Diagnostic energy [keV]", ylabel = "Counts [s^-1]", normalize = false)

Plot the signal S and the WF signal in the same plot. The corresponding measurement bins are Ed_S and Ed_WF, respectively. By default, use green color for the WF and
don't normalize the signals.
"""
function plotSigComp(Ed_S::AbstractVector, S::AbstractVector, Ed_WF::AbstractVector, WF::AbstractVector; lc_S = :black, lc_WF = :green3, lw::Float64 = 2.5, xlabel::String = "Diagnostic energy [keV]", ylabel::String = "Counts [s^-1]", normalize::Bool = false)
    if normalize 
        S = S ./ maximum(S)
        WF = WF ./ maximum(WF)
    end

    Plots.plot(Ed_S, S, color=lc_S, lw = lw, label="S")
    Plots.plot!(Ed_WF, WF, color=lc_WF, lw = lw, xlabel=xlabel, ylabel=ylabel, title="S vs WF", label="WF")
end

"""
    plot_fE_comp(F_os_3D, E_array, pm_array, Rm_array, filepath_distr)
    plot_fE_comp(-||-; verbose=false, logplot=true, os_equidistant=true)

Given a fast-ion orbit-space distribution and its grid points, defined by F_os_3D, E_array, pm_array and Rm_array, respectively, compute the
integrated energy dependance f(E) and compare it with the energy dependence of the fast-ion (E,p,R,z) distribution given by the data in 
filepath_distr.

By default, assume
- no verbose printing
- log-y plot wanted
- that the grid in orbit-space is equidistant
"""
function plot_fE_comp(F_os_3D::Array{Float64,3}, E_array::AbstractArray, pm_array::AbstractArray, Rm_array::AbstractArray, filepath_distr::String; verbose::Bool=false, logplot::Bool=true, os_equidistant::Bool=true)

    # Determine filepath_distr file extension
    fileext_distr = (split(filepath_distr,"."))[end] # Assume last part after final '.' is the file extension
    if (lowercase(fileext_distr) == "h5") || (lowercase(fileext_distr) == "hdf5")
        verbose && println("Loading .h5 file and extracting f, E, p, R and z... ")
        F_EpRz, energy, pitch, R, z = h5to4D(filepath_distr; verbose = verbose) # Load the (E,p,R,z) fast-ion distribution. For energy comparison.
    elseif lowercase(fileext_distr) == "jld2"
        F_EpRz, energy, pitch, R, z = JLD2to4D(filepath_distr; verbose = verbose) # Load the (E,p,R,z) fast-ion distribution. For energy comparison.
    else
        error("Unknown file extension for fast-ion distribution file. Please input a .h5 or .jld2 file.")
    end

    verbose && println("Calculating the 4D differences... ")
    dE4D, dp4D, dR4D, dZ4D = get4DDiffs(energy, pitch, R, z)

    verbose && println("Calculting fE_ps... ")
    fr = F_EpRz .* reshape(R,(1,1,length(R),1)) # Multiply by R in the R-dimension, to account for toroidal symmetry
    fE_ps = dropdims(sum(fr.*dp4D.*dR4D.*dZ4D, dims=(2,3,4))*(2*pi), dims=(2,3,4)) # Integrate out all dimensions but the energy.

    verbose && println("Calculting fE_os... ")
    if os_equidistant
        dpm = abs(pm_array[2]-pm_array[1])
        dRm = abs(Rm_array[2]-Rm_array[1])

        fE_os = dropdims(sum(F_os_3D, dims=(2,3))*dpm*dRm,dims=(2,3))
    else
        throw(ArgumentError("ERROR. os_equidistant=false not implemented yet."))
        # TO BE ADDED TO THE OWCF. MAYBE.
    end

    verbose && println("Plotting... ")
    good_inds_ps = findall(x-> x>0, fE_ps) # Because of log-plot, find all non-zero elements. Only plot those.
    good_inds_os = findall(x-> x>0, fE_os) # Because of log-plot, find all non-zero elements. Only plot those.
    if logplot
        Plots.plot(energy[good_inds_ps], fE_ps[good_inds_ps], label="f(E) [E, p, R, z]", xlabel="E [keV]", yaxis=:log10, ylabel= "f(E)", lw = 2.5, linecolor="green", legend=:bottomleft)
        Plots.plot!(E_array[good_inds_os], fE_os[good_inds_os], label="f(E) [E, pm, Rm]", xlabel="E [keV]", yaxis=:log10, ylabel="f(E)", lw = 2.5, linecolor="blue", legend=:bottomleft)
    else
        Plots.plot(energy[good_inds_ps], fE_ps[good_inds_ps], label="f(E) [E, p, R, z]", xlabel="E [keV]", ylabel= "f(E)", lw = 2.5, linecolor="green")
        Plots.plot!(E_array[good_inds_os], fE_os[good_inds_os], label="f(E) [E, pm, Rm]", xlabel="E [keV]", ylabel="f(E)", lw = 2.5, linecolor="blue")
    end
end

"""
    plot_fE_comp(F_os, og, filepath_distr)
    plot_fE_comp(-||-, verbose=false, os_equidistant=true, logplot=true)

Take the fast-ion orbit-space distribution (compressed into a 1D vector format), inflate it to its 3D form and plot f(E) versus
the (E,p,R,z) fast-ion distribution data in filepath_distr.

By default, assume
- no verbose printing
- that the orbit-space grid is equidistant
"""
function plot_fE_comp(F_os::AbstractArray, og::OrbitGrid, filepath_distr::String; verbose::Bool=false, os_equidistant::Bool=true, kwargs... )

    verbose && println("Mapping orbits to grid... ")
    F_os_3D = map_orbits_OWCF(og, F_os, os_equidistant) # Map the 1D orbit-space distribution vector to the full 3D grid
    E_array = og.energy
    pm_array = og.pitch
    Rm_array = og.r

    plot_fE_comp(F_os_3D, E_array, pm_array, Rm_array, filepath_distr; verbose=verbose, os_equidistant=os_equidistant, kwargs... )
end

"""
    plot_fpm(F_os_3D, E_array, pm_array, Rm_array)
    plot_fpm(-||-; verbose=false, os_equidistant=true, logplot=false, lgd=:topright)

Given a fast-ion orbit-space distribution and its grid points, defined by F_os_3D, E_array, pm_array and Rm_array, respectively, compute the
integrated pm dependance f(pm) and plot it.

By default, assume
- no verbose printing
- that the grid in orbit-space is equidistant
- log-y plot not wanted
- that the figure legend should be in the top-right corner of the plot
"""
function plot_fpm(F_os_3D::Array{Float64,3}, E_array::AbstractArray, pm_array::AbstractArray, Rm_array::AbstractArray; verbose::Bool=false, os_equidistant::Bool=true, logplot::Bool=false,lgd=:topright)
    if os_equidistant
        dE = abs(E_array[2]-E_array[1])
        dRm = abs(Rm_array[2]-Rm_array[1])

        fpm_os = dropdims(sum(F_os_3D, dims=(1,3))*dE*dRm,dims=(1,3))
    else
        throw(ArgumentError("ERROR. os_equidistant=false not implemented yet."))
        # TO BE WRITTEN
    end

    verbose && println("Plotting... ")
    if logplot
        good_inds_os = findall(x-> x>0, fpm_os) # Because of log-plot, find all non-zero elements. Only plot those.
        Plots.plot(pm_array[good_inds_os], fpm_os[good_inds_os], label="f(pm) [E, pm, Rm]", xlabel="pm [-]", yaxis=:log, ylabel="f(pm)", lw = 2.5, linecolor="blue", legend=lgd)
    else
        Plots.plot(pm_array, fpm_os, label="f(pm) [E, pm, Rm]", xlabel="pm [-]", ylabel="f(pm)", lw = 2.5, linecolor="blue", legend=lgd)
    end
end

"""
    plot_fpm(F_os, og)
    plot_fpm(-||-; verbose=false, os_equidistant=true, logplot=false)

Take the fast-ion orbit-space distribution (compressed into a 1D vector format), inflate it to its 3D form and plot f(pm).

By default, assume
- no verbose printing
- that the orbit-space grid is equidistant
"""
function plot_fpm(F_os::AbstractArray, og::OrbitGrid; verbose::Bool=false, os_equidistant::Bool=true, kwargs... )
    verbose && println("Mapping orbits to grid... ")
    F_os_3D = map_orbits_OWCF(og, F_os, os_equidistant) # Map the 1D orbit-space distribution vector to the full 3D grid
    plot_fpm(F_os_3D, og.energy, og.pitch, og.r; verbose=verbose, os_equidistant=os_equidistant, kwargs... )
end

"""
    plot_fRm(F_os_3D, E_array, pm_array, Rm_array)
    plot_fRm(-||-; verbose=false, os_equidistant=true, logplot=false, lgd=:topright)

Given a fast-ion orbit-space distribution and its grid points, defined by F_os_3D, E_array, pm_array and Rm_array, respectively, compute the
integrated Rm dependance f(Rm) and plot it.

By default, assume
- no verbose printing
- that the grid in orbit-space is equidistant
- log-y plot not wanted
- that the figure legend should be in the top-right corner of the plot
- that the f(Rm) curve should not be normalized
"""
function plot_fRm(F_os_3D::Array{Float64,3}, E_array::AbstractArray, pm_array::AbstractArray, Rm_array::AbstractArray; verbose::Bool=false, os_equidistant::Bool=true, logplot::Bool=false,lgd=:topright, normalize::Bool=false)
    if os_equidistant
        dE = abs(E_array[2]-E_array[1])
        dpm = abs(pm_array[2]-pm_array[1])

        fRm_os = dropdims(sum(F_os_3D, dims=(1,2))*dE*dpm,dims=(1,2))
    else
        throw(ArgumentError("ERROR. os_equidistant=false not implemented yet."))
        # TO BE WRITTEN
    end

    verbose && println("Plotting... ")
    if logplot
        good_inds_os = findall(x-> x>0, fRm_os) # Because of log-plot, find all non-zero elements. Only plot those.
        fRm_os = fRm_os[good_inds_os]
        if normalize
            fRm_os = fRm_os ./ maximum(fRm_os)
        end
        Plots.plot(Rm_array[good_inds_os], fRm_os, label="f(Rm) [E, pm, Rm]", xlabel="Rm [m]", yaxis=:log, ylabel="f(Rm) [m^-1]", lw = 2.5, linecolor="blue", legend=lgd)
    else
        if normalize
            fRm_os = fRm_os ./ maximum(fRm_os)
        end
        Plots.plot(Rm_array, fRm_os, label="f(Rm) [E, pm, Rm]", xlabel="Rm [m]", ylabel="f(Rm) [m^-1]", lw = 2.5, linecolor="blue", legend=lgd)
    end
end

"""
    plot_fRm(F_os, og)
    plot_fRm(-||-; verbose=false, os_equidistant=true, logplot=false)

Take the fast-ion orbit-space distribution (compressed into a 1D vector format), inflate it to its 3D form and plot f(Rm).

By default, assume
- no verbose printing
- that the orbit-space grid is equidistant
"""
function plot_fRm(F_os::AbstractArray, og::OrbitGrid; verbose::Bool=false, os_equidistant::Bool=true, kwargs...)
    verbose && println("Mapping orbits to grid... ")
    F_os_3D = map_orbits_OWCF(og, F_os, os_equidistant) # Map the 1D orbit-space distribution vector to the full 3D grid

    plot_fRm(F_os_3D, og.energy, og.pitch, og.r; verbose=verbose, os_equidistant=os_equidistant, kwargs... )
end

"""
    plot_F_Ep(F_EpRz, energy, pitch, R, z)
    plot_F_Ep(-||-, verbose=false, xaxis=extrema(energy), kwargs... )

Take the F_EpRz (E,p,R,z) fast-ion distribution, integrate over (R,z) and plot the E,p-dependance, i.e. f(E,p).
"""
function plot_F_Ep(F_EpRz::Array{Float64,4}, energy::Array{Float64,1}, pitch::Array{Float64,1}, R::Array{Float64,1}, z::Array{Float64,1}; verbose::Bool=false, xaxis=extrema(energy), kwargs...)

    verbose && println("Calculating the 4D differences... ")
    dE4D, dp4D, dR4D, dZ4D = get4DDiffs(energy, pitch, R, z)

    F_R = F_EpRz .* reshape(R,(1,1,length(R),1)) # Multiply by R in the R-dimension, to account for toroidal symmetry
    F_Ep = dropdims(sum(F_R.*dR4D.*dZ4D, dims=(3,4))*(2*pi), dims=(3,4)) # Integrate out R,z

    Plots.heatmap(energy,pitch, fEp', xlabel="Energy [keV]", ylabel="Pitch [-]", xaxis=xaxis, fillcolor=cgrad([:white, :darkblue, :green, :yellow, :orange, :red]), kwargs... )
end

"""
    plot_fEp(filepath_distr)
    plot_fEp(filepath_distr,verbose=true, kwargs... )

Load the F_EpRz (E,p,R,z) fast-ion distribution from file (.h5/.hdf5/.jld2), integrate over (R,z) and plot the E,p-dependance, i.e. f(E,p).
"""
function plot_F_Ep(filepath_distr::String; verbose::Bool=false, kwargs...)
    # Determine filepath_distr file extension
    fileext_distr = (split(filepath_distr,"."))[end] # Assume last part after final '.' is the file extension
    if (lowercase(fileext_distr) == "h5") || (lowercase(fileext_distr) == "hdf5")
        verbose && println("Loading .h5 file and extracting f, E, p, R and z... ")
        F_EpRz, energy, pitch, R, z = h5to4D(filepath_distr; verbose = verbose) # Load the (E,p,R,z) fast-ion distribution. For energy comparison.
    elseif lowercase(fileext_distr) == "jld2"
        verbose && println("Loading .jld2 file and extracting f, E, p, R and z... ")
        F_EpRz, energy, pitch, R, z = JLD2to4D(filepath_distr; verbose = verbose) # Load the (E,p,R,z) fast-ion distribution. For energy comparison.
    else
        error("Unknown file extension for fast-ion distribution file. Please input a .h5 or .jld2 file.")
    end

    plot_F_Ep(F_EpRz, energy, pitch, R, z; verbose=verbose, kwargs...)
end

"""
    plot_F_Rz(F_EpRz, energy, pitch, R, z)
    plot_F_Rz(-||-, wall=wall, verbose=true, kwargs...)

Take the (E,p,R,z) fast-ion distribution F_EpRz, integrate over (E,p) and plot the R,z-dependance, i.e. f(R,z). Include the tokamak wall, if available.
"""
function plot_F_Rz(F_EpRz::Array{Float64,4}, energy::Array{Float64,1}, pitch::Array{Float64,1}, R::Array{Float64,1}, z::Array{Float64,1}; wall::Union{Nothing,Boundary{T}}=nothing, verbose::Bool=false, kwargs...) where {T}

    verbose && println("Calculating the 4D differences... ")
    dE4D, dp4D, dR4D, dZ4D = get4DDiffs(energy, pitch, R, z)

    F_Rz = dropdims(sum(F_EpRz.*dE4D.*dp4D, dims=(1,2)), dims=(1,2)) # Integrate out energy, pitch

    # Convert from cm to m
    if minimum(R)>100.0
        R = R./100.0
        z = z./100.0
    end

    verbose && println("Plotting... ")
    Plots.heatmap(R,z, F_Rz', xlabel="R [m]", ylabel="z[m]", aspect_ratio=:equal,fillcolor=cgrad([:white, :darkblue, :green, :yellow, :orange, :red]))
    if !(wall==nothing)
        Plots.plot!(wall.r,wall.z, color=:black, legend=false)
    end
end

"""
    plot_F_Rz(filepath_distr)
    plot_F_Rz(-||-, verbose=false, filepath_equil=nothing, clockwise_phi_false, kwargs... )

Load the (E,p,R,z) fast-ion distribution F_EpRz from file (.h5/.hdf5/.jld2), integrate over (E,p) and plot the R,z-dependance, i.e. f(R,z). 
Also load and include the tokamak wall, if available.
"""
function plot_F_Rz(filepath_distr::String; verbose::Bool=false, filepath_equil::Bool=nothing, clockwise_phi::Bool=false, kwargs...)
    # Determine filepath_distr file extension
    fileext_distr = (split(filepath_distr,"."))[end] # Assume last part after final '.' is the file extension
    if (lowercase(fileext_distr) == "h5") || (lowercase(fileext_distr) == "hdf5")
        verbose && println("Loading .h5 file and extracting f, E, p, R and z... ")
        F_EpRz, energy, pitch, R, z = h5to4D(filepath_distr; verbose = verbose) # Load the (E,p,R,z) fast-ion distribution. For energy comparison.
    elseif lowercase(fileext_distr) == "jld2"
        verbose && println("Loading .jld2 file and extracting f, E, p, R and z... ")
        F_EpRz, energy, pitch, R, z = JLD2to4D(filepath_distr; verbose = verbose) # Load the (E,p,R,z) fast-ion distribution. For energy comparison.
    else
        error("Unknown file extension for fast-ion distribution file. Please input a .h5 or .jld2 file.")
    end

    if !(filepath_equil==nothing)
        verbose && println("Loading tokamak equilibrium, extracting wall info... ")
        if ((split(filepath_equil,"."))[end] == "eqdsk") || ((split(filepath_equil,"."))[end] == "geqdsk")
            M, wall = read_geqdsk(filepath_equil,clockwise_phi=clockwise_phi) # Assume counter-clockwise phi-direction
        else # Otherwise, assume magnetic equilibrium is a saved .jld2 file
            myfile = jldopen(filepath_equil,false,false,false,IOStream)
            wall = myfile["wall"]
            close(myfile)
        end
        plot_F_Rz(F_EpRz, energy, pitch, R, z; wall=wall, verbose=verbose, kwargs... )
    else
        plot_F_Rz(F_EpRz, energy, pitch, R, z; verbose=verbose, kwargs... )
    end
end

################################################################################################
# Define the orbit pile plot function that will be able to plot split of diagnostic signals
# into their orbit-type constituents (as in Figure 10 and 11 in H. Järleblad et al, NF, 2022)
@userplot OrbPile

@recipe function f(of::OrbPile)
    weights, returns = of.args
    weights = cumsum(weights,dims=2)
    seriestype := :shape
    orbcolors = [:red :blue :green :purple :orange :pink] # Assume maximum 6 orbit types
    for c=1:size(weights,2)
        sx = vcat(weights[:,c], c==1 ? zeros(length(returns)) : reverse(weights[:,c-1]))
        sy = vcat(returns, reverse(returns))
        @series begin
            fillcolor := orbcolors[c]
            Shape(sy, sx)
        end
    end
end
################################################################################################

"""
    plotDistrOrbSplit(ps2WFoutputFilepath)
    plotDistrOrbSplit(ps2WFoutputFilepath; verbose=false, save_plot=false, normalize=false, log_plotting=false)

Plot the orbit constituents of the fast-ion distribution using an output file from ps2WF.jl when 'calcWFOs' was set to true.
The orbit constituents will be plotted as a bar plot. To visualize the energy, pm or Rm dependency of the orbit constituents, 
please use the apps/signalWebApp.jl.

verbose - If true, the function will talk a lot!
save_plot - If true, a .png file will be saved of the plot.
normalize - If true, the orbit-type constitutents will be normalized so that their sum is 1.
log_plotting - If true, manual override to use log-scale for y-axis. If false, manual override to user linear-scale for y-axis.
"""
function plotDistrOrbSplit(ps2WFoutputFilepath::String; verbose::Bool = false, save_plot::Bool = false, normalize::Bool=false, log_plotting::Union{Nothing,Bool}=nothing)
    myfile = jldopen(ps2WFoutputFilepath,false,false,false,IOStream)
    if !(haskey(myfile,"WFO_E"))
        error("Specified ps2WF output file did not have the data necessary for splitting signal into orbit types. Please try re-running ps2WF.jl with 'calcWFOs' set to true and use that output file as input to this function.")
    end

    verbose && println("Loading data from ps2WF output file... ")
    E_array = myfile["E_array"]
    FO_E = myfile["FO_E"]
    nfast = myfile["nfast_orig"]
    close(myfile)

    verbose && println("(a) Total number of fast ions (originally): $(nfast)")

    dE = abs(E_array[2]-E_array[1]) # Assume equidistant energy grid points
    if length(size(FO_E))==3
        FO_E = FO_E[1,:,:]
    end

    stagnation_part = sum(dE .*FO_E[:,1])
    trapped_part = sum(dE .*FO_E[:,2])
    copassing_part = sum(dE .*FO_E[:,3])
    counterpassing_part = sum(dE .*FO_E[:,4])
    potato_part = sum(dE .*FO_E[:,5])
    counterstagnation_part = sum(dE .*FO_E[:,6])
    constituents = [stagnation_part, trapped_part, copassing_part, counterpassing_part, potato_part, counterstagnation_part]
    verbose && println("(b) Sum of constituents: $(sum(constituents))")
    verbose && println("Any discrepancy between (a) and (b) is due to numerical inaccuracies.")

    title_F = "Orbit-type constituents of \n the fast-ion distribution"
    if !normalize
        ylabel_F = "Fast ions"
    else
        ylabel_F = "Fractions"
        constituents = constituents ./ norm(constituents)
    end

    min_y, max_y = extrema(constituents)
    if min_y==0.0
        p = sortperm(constituents) # Return an array p with indices that puts the elements of constituents in ascending order
        min_y = (constituents[p])[findfirst(x-> x>0.0, constituents[p])] # Only ok coding because constituents is so short. And efficiency is not an issue here
    end
    min_OOM, max_OOM = (floor(log10(min_y)),ceil(log10(max_y))) # The orders of magnitude of the minimum and maximum values
    if !((max_OOM-min_OOM)==0.0) && ((max_y/min_y) > 10)
        verbose && println("Log-plotting activated (between 10^$(min_OOM) and 10^$(max_OOM))... ")
        yaxis_scale_F = :log10
        ylims_F = 10 .^(min_OOM, max_OOM)
        yticks_F = 10 .^collect(min_OOM:1:max_OOM)
    else
        verbose && println("Linear plotting... ")
        yaxis_scale_F = :identity
        ylims_F = (min_y, max_y)
        yticks_F = round.(collect(range(min_y,stop=max_y,length=5)),sigdigits=3)
    end

    if typeof(log_plotting)==Bool
        if !log_plotting && yaxis_scale_F==:log10
            verbose && println("Log-plotting deactivated due to user input.")
            yaxis_scale_F = :identity
            ylims_F = (min_y, max_y)
            yticks_F = round.(collect(range(min_y,stop=max_y,length=5)),sigdigits=3)
        end
        if log_plotting && yaxis_scale_F==:identity
            verbose && println("Log-plotting activated due to user input.")
            yaxis_scale_F = :log10
            ylims_F = 10 .^(min_OOM, max_OOM)
            yticks_F = 10 .^collect(min_OOM:1:max_OOM)
        end
    end


    verbose && println("Plotting... ")
    dx = 0.2
    my_plt = Plots.plot(Shape([1.0-2*dx,1.0,1.0,1.0-2*dx],[1.0e-9,1.0e-9,constituents[1],constituents[1]]),color=:red,label="Stagnation")
    Plots.plot!(my_plt, Shape([2.0-2*dx,2.0,2.0,2.0-2*dx],[1.0e-9,1.0e-9,constituents[2],constituents[2]]),color=:blue,label="Trapped")
    Plots.plot!(my_plt, Shape([3.0-2*dx,3.0,3.0,3.0-2*dx],[1.0e-9,1.0e-9,constituents[3],constituents[3]]),color=:green,label="Co-passing")
    Plots.plot!(my_plt, Shape([4.0-2*dx,4.0,4.0,4.0-2*dx],[1.0e-9,1.0e-9,constituents[4],constituents[4]]),color=:purple,label="Counter-passing")
    Plots.plot!(my_plt, Shape([5.0-2*dx,5.0,5.0,5.0-2*dx],[1.0e-9,1.0e-9,constituents[5],constituents[5]]),color=:orange,label="Potato")
    Plots.plot!(my_plt, Shape([6.0-2*dx,6.0,6.0,6.0-2*dx],[1.0e-9,1.0e-9,constituents[6],constituents[6]]),color=:pink,label="Counter-stagnation")
    Plots.plot!(my_plt, yaxis=yaxis_scale_F, ylims=ylims_F, yticks=yticks_F, xticks=false, title=title_F,ylabel=ylabel_F,xlabel="Orbit types")

    Plots.display(my_plt)

    if save_plot
        verbose && println("Saving plot... ")
        my_plt = Plots.plot!(my_plt, dpi=600)
        sadd = (normalize ? "norm" : "abs")
        saddd = (yaxis_scale_F==:log10 ? "_log10" : "")
        png(my_plt, "plotDistrOrbSplits_"*sadd*saddd)
    end
end

"""
    plotMeasurementOrbSplits(Ed_WF, measurement_bin_of_interest, WFO_E, E_array)
    plotMeasurementOrbSplits(-||-; verbose=false, save_plot=false, normalize=false, ylabel = (normalize ? "Normalized signal [a.u.]" : "Counts [(keV*s)^-1]"), title = (normalize ? "Orbit-split of diagnostic measurement \n (normalized)" : "Orbit-split of diagnostic measurement"), debug=false)

Plot bar plot showing the orbit constituents of the diagnostic measurement at the measurement bin of interest.
Ed_WF are the measurement bins. measurement_bin_of_interest if the measurement bin of interest. WFO_E is the WFO_E quantity from the ps2WF.jl output file, run with calcWFOs set to true.
E_array are the fast-ion energy grid points. 

If verbose, then the function will talk a lot!
If save_plot, then the plot will be saved in .png format.
If normalize, then the bars will be plotted as fractions of the total measurement.
ylabel sets the label of the y-axis.
title sets the title of the bar plot.
debug turns on debug mode.
"""
function plotMeasurementOrbSplits(Ed_WF::Union{Vector{Float64},Vector{Int64}}, measurement_bin_of_interest::Union{Float64,Int64}, WFO_E::Array{Float64,3}, E_array::AbstractVector; verbose::Bool=false, save_plot::Bool=false, normalize::Bool=false, ylabel = (normalize ? "Normalized signal [a.u.]" : "Counts [(keV*s)^-1]"), title = (normalize ? "Orbit-split of diagnostic measurement \n (normalized)" : "Orbit-split of diagnostic measurement"), debug=false)
    verbose && println("Finding measurement bin closest to requested bin... ")
    iEd = closest_index(Ed_WF, measurement_bin_of_interest)
    dE = abs(diff(E_array)[1]) # Assume equidistant
    verbose && println("Computing orbit contributions to measurement... ")
    stagnation_contribution = sum(dE .*WFO_E[iEd,:,1])
    trapped_contribution = sum(dE .*WFO_E[iEd,:,2])
    copassing_contribution = sum(dE .*WFO_E[iEd,:,3])
    counterpassing_contribution = sum(dE .*WFO_E[iEd,:,4])
    potato_contribution = sum(dE .*WFO_E[iEd,:,5])
    counterstagnation_contribution = sum(dE .*WFO_E[iEd,:,6])

    verbose && println("Determining plot axis data... ")
    contributions = [stagnation_contribution, trapped_contribution, copassing_contribution, counterpassing_contribution, potato_contribution, counterstagnation_contribution]
    if normalize
        contributions = contributions ./norm(contributions)
    end
    min_y, max_y = extrema(contributions)
    if min_y==0.0
        p = sortperm(contributions) # Return an array p with indices that puts the elements of contributions in ascending order
        min_y = (contributions[p])[findfirst(x-> x>0.0, contributions[p])] # Only ok coding because contributions is so short. And efficiency is not an issue here
    end
    min_OOM, max_OOM = (floor(log10(min_y)),ceil(log10(max_y))) # The orders of magnitude of the minimum and maximum values
    if !((max_OOM-min_OOM)==0.0) && ((max_y/min_y) > 10)
        verbose && println("Log-plotting activated... ")
        yaxis_scale = :log10
        ylims = (min_OOM, max_OOM)
        yticks = collect(min_OOM:1:max_OOM)
    else
        verbose && println("Linear plotting... ")
        yaxis_scale = :identity
        ylims = (min_y, max_y)
        yticks = round.(collect(range(min_y,stop=max_y,length=5)),sigdigits=3)
    end

    verbose && println("Plotting... ")
    dx = 0.2
    my_plt = Plots.plot(Shape([1.0-2*dx,1.0,1.0,1.0-2*dx],[1.0e-4,1.0e-4,contributions[1],contributions[1]]),color=:red,label="Stagnation")
    Plots.plot!(my_plt, Shape([2.0-2*dx,2.0,2.0,2.0-2*dx],[1.0e-4,1.0e-4,contributions[2],contributions[2]]),color=:blue,label="Trapped")
    Plots.plot!(my_plt, Shape([3.0-2*dx,3.0,3.0,3.0-2*dx],[1.0e-4,1.0e-4,contributions[3],contributions[3]]),color=:green,label="Co-passing")
    Plots.plot!(my_plt, Shape([4.0-2*dx,4.0,4.0,4.0-2*dx],[1.0e-4,1.0e-4,contributions[4],contributions[4]]),color=:purple,label="Counter-passing")
    Plots.plot!(my_plt, Shape([5.0-2*dx,5.0,5.0,5.0-2*dx],[1.0e-4,1.0e-4,contributions[5],contributions[5]]),color=:orange,label="Potato")
    Plots.plot!(my_plt, Shape([6.0-2*dx,6.0,6.0,6.0-2*dx],[1.0e-4,1.0e-4,contributions[6],contributions[6]]),color=:pink,label="Counter-stagnation")
    Plots.plot!(my_plt, yaxis=yaxis_scale, ylims=(10 .^ylims),yticks=(10 .^yticks), xticks=false, title=title,ylabel=ylabel,xlabel="Orbit types", top_margin=3Plots.mm, xlims=(0.4,7.0))

    Plots.display(my_plt)

    if save_plot
        verbose && println("Saving plot... ")
        my_plt = Plots.plot!(my_plt, dpi=600)
        sadd = (normalize ? "norm" : "abs")
        png(my_plt, "plotMeasurementOrbSplits_"*sadd)
    end
end

"""
    plotSignalOrbSplits(S_WF, Ed_WF, WFO_E, WFO_pm, WFO_Rm, E_array, pm_array, Rm_array)
    plotSignalOrbSplits(-||-; normalize=false, verbose=false, save_plot=false, i0=1, iend=0, xlabel="Deposited energy [keV]", ylabel = (normalize ? "Normalized signal [a.u.]" : "Counts [(keV*s)^-1]"), legend_margin = 39, title = (normalize ? "Orbit-split of diagnostic signal" : "Orbit-split of diagnostic signal (normalized)"), WF_marker_color = :green3)

Take the ps2WF.jl results data and plot an orbit split of the synthetic diagnostic signal. 

By default, assume
- The user does NOT want a normalized orbit split (Figure 10 vs 11 in H. Järleblad et al, NF, 2022)
- The figure should not be saved
- The split should start from the first diagnostic measurement bin (i0=1)
- The split should finish at the last diagnostic measurement bin (iend=0)(increasing 'iend' means cut more and more of upper diagnostic measurement interval)
- The measurement spectrum x-axis corresponds to deposited energy in keV
- The measurement spectrum y-axis corresponds to particle counts per keV per second (or arbitrary units [a.u.] in the case of normalized split)
- To create 39 pixels of margin for the figure legend
- The title to be *please see the title keyword alternatives*
- The marker color of the WF scatter plot should be bright green (:green3)
"""
function plotSignalOrbSplits(S_WF::AbstractVector, Ed_WF::AbstractVector, WFO_E::Array{Float64,3}, WFO_pm::Array{Float64,3}, WFO_Rm::Array{Float64,3}, E_array::AbstractVector, pm_array::AbstractVector, Rm_array::AbstractVector; normalize::Bool=false, verbose::Bool=false, save_plot::Bool=false, i0=1, iend=0, xlabel="Deposited energy [keV]", ylabel = (normalize ? "Normalized signal [a.u.]" : "Counts [(keV*s)^-1]"), legend_margin = 39, title = (normalize ? "Orbit-split of diagnostic signal (normalized)" : "Orbit-split of diagnostic signal"), WF_marker_color = :green3)
    dE = abs(diff(E_array)[1]) # Assume equidistant
    dpm = abs(diff(pm_array)[1]) # Assume equidistant
    dRm = abs(diff(Rm_array)[1]) # Assume equidistant

    EdvsOrbs_WFOE = dropdims(dE .*sum(WFO_E, dims=2),dims=2)
    EdvsOrbs_WFOpm = dropdims(dpm .*sum(WFO_pm, dims=2),dims=2) # Not really needed, should be identical to EdvsOrbs_WFOE if safety plot looks fine
    EdvsOrbs_WFORm = dropdims(dRm .*sum(WFO_Rm, dims=2),dims=2) # Not really needed, should be identical to EdvsOrbs_WFOE if safety plot looks fine

    verbose && println("Computing plot weights for signal split plot... ")
    weights = zeros(size(EdvsOrbs_WFOE))
    for i=1:size(EdvsOrbs_WFOE,1)
        rowsum = sum(EdvsOrbs_WFOE[i,:])
        if rowsum != 0.0
            rowsum = (normalize ? rowsum : 1.0)
            weights[i,:] = EdvsOrbs_WFOE[i,:] ./rowsum
        end
    end

    verbose && print("Plotting... ")
    verbose && normalize && print("(normalized)... ")
    verbose && println("")

    dx = (maximum(Ed_WF)-minimum(Ed_WF))*0.005
    my_plt = Plots.plot(OrbPile((weights[i0:end-iend,:], Ed_WF[i0:end-iend])), labels=["Stagnation" "Trapped" "Co-passing" "Counter-passing" "Potato" "Counter-stagnation"], legend=:topright, xlabel=xlabel, ylabel=ylabel, xlims=(minimum(Ed_WF)-dx, maximum(Ed_WF)+legend_margin*dx))
    S_WF_norm = (normalize ? maximum(S_WF[i0:end-iend]) : 1.0)
    my_plt = Plots.scatter!(my_plt, Ed_WF[i0:end-iend], S_WF[i0:end-iend] ./S_WF_norm , markerstrokealpha=1.0, markerstrokecolor=:white, markercolor=WF_marker_color, label="WF", markerstrokewidth=1.5, title=title)
    Plots.display(my_plt)

    if save_plot
        my_plt = Plots.plot!(my_plt, dpi=600)
        sadd = (normalize ? "norm" : "abs")
        png(my_plt, "plotSignalOrbSplits_"*sadd)
    end
end

"""
    plotSafetyPlot(S_WF, Ed_WF, WFO_E, WFO_pm, WFO_Rm)
    plotSafetyPlot(-||-; save_plot=false, verbose=false)

Take the provided input data, which is part of the ps2WF.jl output data, and plot a safety plot to confirm the
validity of the data.

By default, assume
- That the plot should NOT be saved
- That the function should NOT be talkative
"""
function plotSafetyPlot(S_WF::AbstractVector, Ed_WF::AbstractVector, WFO_E::Array{Float64,3}, WFO_pm::Array{Float64,3}, WFO_Rm::Array{Float64,3}, E_array::AbstractVector, pm_array::AbstractVector, Rm_array::AbstractVector; save_plot::Bool = false, verbose::Bool=false)
    dE = abs(diff(E_array)[1]) # Assume equidistant
    dpm = abs(diff(pm_array)[1]) # Assume equidistant
    dRm = abs(diff(Rm_array)[1]) # Assume equidistant
    
    S_WF_WFOE = dropdims(dE .*sum(WFO_E,dims=(2,3)),dims=(2,3)) # Sum over all energies and orbit types
    S_WF_WFOpm = dropdims(dpm .*sum(WFO_pm,dims=(2,3)),dims=(2,3)) # Sum over all pm values and orbit types
    S_WF_WFORm = dropdims(dRm .*sum(WFO_Rm,dims=(2,3)),dims=(2,3)) # Sum over all Rm values and orbit types

    verbose && println("Plotting safety plot... ")
    my_plt = Plots.plot(Ed_WF,S_WF,label="S_WF",linewidth=8)
    my_plt = Plots.plot!(my_plt, Ed_WF,S_WF_WFOE,label="S_WF_WFOE",linewidth=6)
    my_plt = Plots.plot!(my_plt, Ed_WF,S_WF_WFOpm,label="S_WF_WFOpm",linewidth=4)
    my_plt = Plots.plot!(my_plt, Ed_WF,S_WF_WFORm,label="S_WF_WFORm",linewidth=2)
    Plots.display(my_plt)

    if save_plot
        png(my_plt, "safety_plot.png")
    end
end

"""
    plotSignalOrbSplits(ps2WFoutputFilepath)
    plotSignalOrbSplits(-||-; verbose=false, scenario=:all, kwargs... )

Load the ps2WF.jl results data stored in the ps2WFoutputFilepath .jld2-file. Depending on the scenario keyword input,
plot either the absolute or normalized orbit splits of the WF signal, or the safety plot, or all of them. Scenarios:
- :all
- :abs
- :norm
- :safety

By default, assume:
- scenario = :all
- Do NOT print verbose print statements

If 'specific_measurement_bin_of_interest' is specified (a Float64 or Int64), then a bar plot showing the orbit-type constituents for that specific measurement bin of
interest will be plotted. 
"""
function plotSignalOrbSplits(ps2WFoutputFilepath::String; verbose::Bool = false, scenario = :all, save_plot::Bool = false, specific_measurement_bin_of_interest::Union{Nothing,Float64,Int64} = nothing, kwargs...)
    myfile = jldopen(ps2WFoutputFilepath,false,false,false,IOStream)
    if !(haskey(myfile,"WFO_E"))
        error("Specified ps2WF output file did not have the data necessary for splitting signal into orbit types. Please try re-running ps2WF.jl with 'calcWFOs' set to true and use that output file as input instead.")
    end

    verbose && println("Loading data from ps2WF output file... ")
    S_WF = myfile["S_WF"]
    E_array = myfile["E_array"]
    pm_array = myfile["pm_array"]
    Rm_array = myfile["Rm_array"]
    Ed_WF = myfile["Ed_array"]
    WFO_E = myfile["WFO_E"]
    WFO_pm = myfile["WFO_pm"]
    WFO_Rm = myfile["WFO_Rm"]
    nfast = myfile["nfast_orig"]
    close(myfile)

    if typeof(specific_measurement_bin_of_interest)==Nothing
        if scenario == :all
            plotSafetyPlot(S_WF, Ed_WF, WFO_E, WFO_pm, WFO_Rm, E_array, pm_array, Rm_array; verbose = verbose, save_plot=save_plot)
            plotSignalOrbSplits(S_WF, Ed_WF, WFO_E, WFO_pm, WFO_Rm, E_array, pm_array, Rm_array; verbose=verbose, save_plot=save_plot, kwargs... ) # Absolute values. Corresponding to Figure 10 in H. Järleblad et al, NF, 2022
            plotSignalOrbSplits(S_WF, Ed_WF, WFO_E, WFO_pm, WFO_Rm, E_array, pm_array, Rm_array; verbose=verbose, save_plot=save_plot, normalize=true, kwargs... ) # Normalized values. Corresponding to Figure 11 in H. Järleblad et al, NF, 2022
        elseif scenario == :abs
            plotSignalOrbSplits(S_WF, Ed_WF, WFO_E, WFO_pm, WFO_Rm, E_array, pm_array, Rm_array; verbose=verbose, save_plot=save_plot, kwargs... ) # Absolute values. Corresponding to Figure 10 in H. Järleblad et al, NF, 2022
        elseif scenario == :norm
            plotSignalOrbSplits(S_WF, Ed_WF, WFO_E, WFO_pm, WFO_Rm, E_array, pm_array, Rm_array; verbose=verbose, save_plot=save_plot, normalize=true, kwargs... ) # Normalized values. Corresponding to Figure 11 in H. Järleblad et al, NF, 2022
        elseif scenario == :safety
            plotSafetyPlot(S_WF, Ed_WF, WFO_E, WFO_pm, WFO_Rm, E_array, pm_array, Rm_array; verbose = verbose, save_plot=save_plot)
        else
            error("Unknown scenario. Accepted values are :all, :abs, :norm and :safety.")
        end
    else
        if scenario == :all
            plotSafetyPlot(S_WF, Ed_WF, WFO_E, WFO_pm, WFO_Rm, E_array, pm_array, Rm_array; verbose = verbose, save_plot=save_plot)
            plotMeasurementOrbSplits(Ed_WF, specific_measurement_bin_of_interest, WFO_E, E_array; verbose=verbose, save_plot=save_plot, kwargs...)
            plotMeasurementOrbSplits(Ed_WF, specific_measurement_bin_of_interest, WFO_E, E_array; verbose=verbose, save_plot=save_plot, normalize=true, kwargs...)
        elseif scenario == :abs
            plotMeasurementOrbSplits(Ed_WF, specific_measurement_bin_of_interest, WFO_E, E_array; verbose=verbose, save_plot=save_plot, kwargs...)
        elseif scenario == :norm
            plotMeasurementOrbSplits(Ed_WF, specific_measurement_bin_of_interest, WFO_E, E_array; verbose=verbose, save_plot=save_plot, normalize=true, kwargs...)
        elseif scenario == :safety
            plotSafetyPlot(S_WF, Ed_WF, WFO_E, WFO_pm, WFO_Rm, E_array, pm_array, Rm_array; verbose = verbose, save_plot=save_plot)
        else
            error("Unknown scenario. Accepted values are :all, :abs, :norm and :safety.")
        end
    end
end


"""
    plot_orbit_movie(M::AbstractEquilibrium, o::Orbit; wall=nothing, top_view=true, verbose=false, anim_numOtau_p = 3, save_anim=false)

Plot an animation of the motion of the guiding-center orbit 'o'. If 'wall' is specified, include tokamak wall in animation. If 'top_view' is true, show the animation 
from a top view of the tokamak, otherwise show the animation from a poloidal cross-sectional view. Let the animation last for 'anim_numOtau_p' number of poloidal transit times.
If 'save_anim' is set to true, save the animation after plotting.
"""
function plot_orbit_movie(M::AbstractEquilibrium, o::Orbit; wall::Union{Nothing,Boundary{T}}=nothing, top_view::Bool=true, verbose::Bool=false, anim_numOtau_p = 3, save_anim::Bool=false) where {T}
    if !(wall==nothing)
        verbose && println("Creating tokamak wall data for topview plot... ")
        R_hfs = minimum(wall.r) # R-coord of high-field side wall
        R_lfs = maximum(wall.r) # R-coord of low-field side wall
        phi = collect(0:1:359).*(2*pi/180.0) # Toroidal angle
        topview_R_hfs_x = (R_hfs).*cos.(phi)
        topview_R_hfs_y = (R_hfs).*sin.(phi)
        topview_R_lfs_x = (R_lfs).*cos.(phi)
        topview_R_lfs_y = (R_lfs).*sin.(phi)
    end

    verbose && println("Defining orbit plot color... ")
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

    if !(wall==nothing)
        verbose && println("Computing magnetic flux surfaces for plotting... ")
        # Define magnetic flux (R,z) grid range
        flux_r = range(minimum(wall.r),stop=maximum(wall.r),length=33)
        flux_z = range(minimum(wall.z),stop=maximum(wall.z),length=33)
        inds = CartesianIndices((length(flux_r),length(flux_z)))
        psi_rz = [M(flux_r[ind[1]], flux_z[ind[2]]) for ind in inds]
        psi_mag, psi_bdry = psi_limits(M)
    end

    if (typeof(o.coordinate)==EPRCoordinate{Float64}) || (typeof(o.coordinate)==EPRCoordinate{Int64})# Float64/Int64 needs to be included. Float/Int sensitive
        verbose && println("(E,pm,Rm) found in orbit. Using...")
        myEPRc = o.coordinate
    else
        verbose && print("(E,mu,Pphi;sigma) found in orbit. Computing (E,pm,Rm)... ")
        # Otherwise, must be HamiltonianCoordinate
        sigma = 1
        if (orb_color == :purple) || (orb_color == :pink)
            sigma = -1
        end
        myEPRc = EPRCoordinate(M, o.coordinate; sigma=sigma, wall=wall, verbose=verbose)
    end

    verbose && println("Computing orbit trajectory data with high temporal resolution... ")
    mygcp = GCParticle(myEPRc) # Define the guiding-center (GC) particle object
    # Integrate. Get the orbit path (path) and the orbit status (stat) objects
    path, stat = integrate(M, mygcp; one_transit=false, r_callback=false, wall=wall, interp_dt=1.0e-10, max_length=5000, tmax=anim_numOtau_p*o.tau_p)
    gcvalid = gcde_check(M, mygcp, path) # Check if usage of the guiding-center approximation was valid, given the length scale of the variation of the magnetic field

    rmax = stat.rm
    if stat.class != :incomplete && stat.class != :lost # If the orbit was not incomplete, nor lost...
        if rmax > myEPRc.r && (false || !isapprox(rmax,myEPRc.r,rtol=1e-4)) # If it ended up outside of its initial R coordinate, or it's not exact...
            stat.class = :invalid # It's invalid!
        end
    else
        stat.tau_p=zero(stat.tau_p) # Else, it's an incomplete or lost orbits, and has no poloidal...
        stat.tau_t=zero(stat.tau_t) # ...and toroidal transit times
    end
    o_long = Orbit(myEPRc,o.class,stat.tau_p,stat.tau_t,path,gcvalid)

    mv_o_x = cos.(o_long.path.phi).*(o_long.path.r) # Compute x-coordinates for the animation orbit trajectory
    mv_o_y = sin.(o_long.path.phi).*(o_long.path.r) # Compute y-coordinates for the animation orbit trajectory

    if top_view
        verbose && println("Creating top-view animation of orbit trajectory... ")
        plt_anim = plot((mv_o_x)[1:2],(mv_o_y)[1:2],color=orb_color,alpha=0.3)
        plt_anim = plot!((mv_o_x)[1:2],(mv_o_y)[1:2],color=orb_color,alpha=0.8)
        plt_anim = scatter!([(mv_o_x)[1]],[(mv_o_y)[1]],color=orb_color)
        if !(wall==nothing)
            plt_anim = plot!(topview_R_hfs_x, topview_R_hfs_y, color=:black)
            plt_anim = plot!(topview_R_lfs_x, topview_R_lfs_y, color=:black)
        end
        # To achieve fast and efficient animation plotting, trick 1 is to not plot the axis ticks...
        plt_anim = plot!(plt_anim,aspect_ratio=:equal, legend=false, xlabel="x [m]", ylabel="y [m]",xticks=true,yticks=true,left_margin=4Plots.mm)
        # ...and trick 2 is to define the needed indexes outside of the @animate loop
        A,B = -500:1:length(mv_o_x)-500, -50:1:length(mv_o_x)-50 # We want negative indexes in for the first 50, 500, etc ones...
        A,B = ifelse.(A .> 0, A, 1), ifelse.(B .> 0, B, 1) # ...so that we can easily put them to 1 in this next line.

        plt_movie = @animate for j=1:30:length(mv_o_x)
            # With the above plt_anim initialization, it is faster now to simply update the coordinates of the already existing plots
            plt_anim[1][1][:x], plt_anim[1][1][:y] = (mv_o_x)[A[j]:j], (mv_o_y)[A[j]:j]
            plt_anim[1][2][:x], plt_anim[1][2][:y] = (mv_o_x)[B[j]:j], (mv_o_y)[B[j]:j]
            plt_anim[1][3][:x], plt_anim[1][3][:y] = [(mv_o_x)[j]], [(mv_o_y)[j]]
        end every 1 # And create a "movie object" 'plt_movie' from every frame ('every 1')
        my_gif = gif(plt_movie,"plot_orbit_movie_top_view.gif",fps=15)
        Plots.display(my_gif)
        if !(save_anim)
            rm("plot_orbit_movie_top_view.gif",force=true)
            verbose && println("'save_anim' set to false. Saved animation automatically removed.")
        end
    else
        verbose && println("Creating poloidal view animation of orbit trajectory... ")
        mv_o_r = o_long.path.r # Compute r-coordinates for the animation orbit trajectory
        mv_o_z = o_long.path.z # Compute z-coordinates for the animation orbit trajectory
        plt_anim = plot(mv_o_r, mv_o_z,color=orb_color,linewidth=2.5)
        plt_anim = scatter!([(mv_o_r)[1]],[(mv_o_z)[1]],color=orb_color)
        if !(wall==nothing)
            plt_crs = Plots.contour!(flux_r,flux_z,psi_rz',levels=collect(range(psi_mag,stop=psi_bdry,length=5)),color=:gray, linestyle=:dot,linewidth=1.5, label="",colorbar=false)
            plt_crs = Plots.plot!(wall.r,wall.z, label="Tokamak wall", color=:black, linewidth=1.5)
        end
        plt_anim = scatter!([magnetic_axis(M)[1]],[magnetic_axis(M)[2]],mc=:gray,label="")
        # To achieve fast and efficient animation plotting, trick 1 is to not plot the axis ticks...
        plt_anim = plot!(plt_anim,aspect_ratio=:equal, legend=false, xlabel="R [m]", ylabel="z [m]",xticks=true,yticks=true,left_margin=4Plots.mm)
        # ...and trick 2 is to define the needed indexes outside of the @animate loop
        A,B = -500:1:length(mv_o_r)-500, -50:1:length(mv_o_z)-50 # We want negative indexes in for the first 50, 500, etc ones...
        A,B = ifelse.(A .> 0, A, 1), ifelse.(B .> 0, B, 1) # ...so that we can easily put them to 1 in this next line.

        plt_movie = @animate for j=1:30:length(mv_o_r)
            plt_anim[1][2][:x], plt_anim[1][2][:y] = [(mv_o_r)[j]], [(mv_o_z)[j]]
        end every 1 # And create a "movie object" 'plt_movie' from every frame ('every 1')
        my_gif = gif(plt_movie,"plot_orbit_movie_crs_view.gif",fps=15)
        Plots.display(my_gif)
        if !(save_anim)
            rm("plot_orbit_movie_crs_view.gif",force=true)
            verbose && println("'save_anim' set to false. Saved animation automatically removed.")
        end
    end
end

"""
    plot_orbit_movie(M, myEPRc; wall=nothing, extra_kw_args=Dict(:toa => true, :limit_phi => true, :maxiter => 0), verbose=false, kwargs...)

Please see documentation for other versions of plot_orbit_movie().
"""
function plot_orbit_movie(M::AbstractEquilibrium, myEPRc::EPRCoordinate; wall::Union{Boundary{T},Nothing}=nothing, extra_kw_args=Dict(:toa => true, :limit_phi => true, :maxiter => 0), verbose=false, kwargs...) where {T}
    verbose && println("Computing orbit from (E,pm,Rm) coordinate... ")
    o = get_orbit(M,myEPRc; wall=wall, interp_dt=1.0e-10, max_length=500, extra_kw_args...) # interp_dt is set to ridiculously small value, to ensure orbit path length of 500
    plot_orbit_movie(M, o; wall=wall, verbose=verbose, kwargs... )
end

"""
    plot_orbit_movie(M, E, pm, Rm; wall=nothing, FI_species="D", extra_kw_args=Dict(:toa => true, :limit_phi => true, :maxiter => 0), verbose=false, anim_numOtau_p = 3, save_anim=false)

Given the magnetic equilibrium and (E,pm,Rm) triplet, compute the resulting orbit and plot an animation of its motion. 'wall' specifies the tokamak first wall, 'FI_species' specifies the fast-ion species (e.g. "D", "T" etc),
'extra_kw_args' specifies extra keyword arguments for the equations-of-motion integration, 'verbose' specifies extra function print output activity, 'anim_numOtau_p' specifies length of animation and 'save_anim' specifies
whether to save the animation to .gif file.
"""
function plot_orbit_movie(M::AbstractEquilibrium, E::Number, pm::Number, Rm::Number; FI_species::String="D", verbose::Bool=false, kwargs...)
    verbose && println("Computing EPRCoordinate from (E,pm,Rm) input... ")
    myEPRc = EPRCoordinate(M, E, pm, Rm; amu=getSpeciesAmu(FI_species), q=getSpeciesEcu(FI_species))
    plot_orbit_movie(M, myEPRc; verbose=verbose, kwargs...)
end