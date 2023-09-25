######################################################################################
# This insignificant appendix preserves a useful PyPlot function utility.
# Because the Plots.jl and PyPlot.jl packages do NOT talk very well with each other,
# they should NEVER be loaded simultaneously.
#
# Therefore, this plot function can NOT be included in gui.jl
#
# It has been preserved in this gui appendix, for future care.
#
# The plot_orbit_distribution() function was originally coded by L. Stagner in 2019.
######################################################################################

using PyPlot, OrbitTomography

"""
    plot_orbit_distribution(E_array, pm_array, Rm_array, F_os_3D)
    plot_orbit_distribution(-||-; save_fig_name="", plot_font_size="medium", vert=false, cmax=[], fig=nothing, ax=nothing, normalize=false)

Take the 3D fast-ion orbit-space distribution F_os_3D and plot f(E,pm), f(pm,Rm) and f(E,Rm). plot_font_size defines font-size.
vert defines plot orientation. cmax defines color-scale maximums. fig defines an already existing figure, for adding to already existing plot.
ax defines the axes of that plot, if already existing.
"""
function plot_orbit_distribution(E_array::AbstractVector, pm_array::AbstractVector, Rm_array::AbstractVector, F_os_3D::Array{T,3}; save_fig_name::String="", plot_font_size::String = "medium", vert::Bool=false, cmax=[], fig=nothing, ax=nothing, normalize::Bool=false) where {T<:Number}

    if length(size(F_os_3D)) != 3
        error("Quantity to be plotted is not 3D. Please input 3D quantity.")
    end
    dr = abs(Rm_array[2]-Rm_array[1]) # Assume equidistant os-grid
    dE = abs(E_array[2]-E_array[1]) # Assume equidistant os-grid
    dp = abs(pm_array[2]-pm_array[1]) # Assume equidistant os-grid

    ep_pdf = dropdims(sum(F_os_3D,dims=3),dims=3)*dr # Integrate out the Rm dimension
    er_pdf = dropdims(sum(F_os_3D,dims=2),dims=2)*dp # Integrate out the pm dimension
    rp_pdf = dropdims(sum(F_os_3D,dims=1),dims=1)*dE # Integrate out the E dimension

    if normalize
        ep_pdf = ep_pdf./maximum(ep_pdf)
        er_pdf = er_pdf./maximum(er_pdf)
        rp_pdf = rp_pdf./maximum(rp_pdf)
    end
    if length(cmax) == 0 # If no maximum for the color scale is defined, use the maximum of the functions
        cmax = zeros(3)
        cmax[1] = maximum(ep_pdf)
        cmax[2] = maximum(rp_pdf)
        cmax[3] = maximum(er_pdf)
    end

    if fig == nothing && ax == nothing # If no figure and axes submitted, create new ones
        if vert # Use vertical plot orientation
            fig, ax = plt.subplots(nrows=3)
            fig.set_size_inches(4,6)
        else # Don't
            fig, ax = plt.subplots(ncols=3)
            fig.set_size_inches(13.5,3)
        end
    end

    c = ax[1].pcolormesh(E_array, pm_array, ep_pdf',vmax=cmax[1],cmap=:inferno,clim=[0,cmax[1]], shading="auto") # Plot it! Use the inferno colormap

    c.set_edgecolor("face") # Set the edgecolor category to face
    #ax[1].set_xticks([70,795,1520,2245,2970]) # Manually set xticks. Change these if needed
    #ax[1].set_yticks([-1.0,-0.5,0.0,0.5,1.0]) # Manually set yticks. Change these if needed
    ax[1].set_ylabel(L"\rm{Pitch\;at\;R_{m}}",labelpad=-4,fontsize=plot_font_size) # Labels
    ax[1].set_xlabel(L"\rm{Energy\;[keV]}",fontsize=plot_font_size) # Labels
    ax[1].set_ylim(extrema(pm_array)) # Set the y-plotlimits to the extrema (minimum,maximum)
    ax[1].set_xlim(extrema(E_array)) # -||-
    cb1 = fig.colorbar(c,ax=ax[1],pad=0.01)
    #cb1.clim(0,cmax[1])
    cb1.formatter.set_powerlimits((0,0))
    cb1.ax.yaxis.set_offset_position("left")
    cb1.update_ticks()

    c2 = ax[2].pcolormesh(Rm_array, pm_array, rp_pdf,vmax=cmax[2],cmap=:inferno,clim=[0,cmax[2]], shading="auto")

    c2.set_edgecolor("face")
    #ax[2].set_yticks([-1.0,-0.5,0.0,0.5,1.0])
    #ax[2].set_xticks([3.125,3.25,3.375,3.5,3.625,3.75])
    ax[2].set_ylabel(L"\rm{Pitch\;at\;R_{m}}",labelpad=-4,fontsize=plot_font_size)
    ax[2].set_xlabel(L"\rm{R_{m}\;[m]}",fontsize=plot_font_size)
    ax[2].set_xlim(extrema(Rm_array))
    ax[2].set_ylim(extrema(pm_array))
    cb2 = fig.colorbar(c2,ax=ax[2],pad=0.01)
    #cb2.clim(0,cmax[2])
    cb2.formatter.set_powerlimits((0,0))
    cb2.ax.yaxis.set_offset_position("left")
    cb2.update_ticks()

    c3 = ax[3].pcolormesh(E_array, Rm_array, er_pdf',vmax=cmax[3],cmap=:inferno,clim=[0,cmax[3]], shading="auto")

    c3.set_edgecolor("face")
    #ax[3].set_yticks([3.125,3.25,3.375,3.5,3.625,3.75])
    #ax[3].set_xticks([70,795,1520,2245,2970])
    ax[3].set_ylabel(L"\rm{R_{m}\;[m]}",fontsize=plot_font_size,labelpad=1)
    ax[3].set_xlabel(L"\rm{Energy\;[keV]}",fontsize=plot_font_size)
    ax[3].set_ylim(extrema(Rm_array))
    ax[3].set_xlim(extrema(E_array))
    cb3 = fig.colorbar(c3,ax=ax[3],pad=0.01)
    #cb3.clim((0,cmax[3]))
    cb3.formatter.set_powerlimits((0,0))
    cb3.ax.yaxis.set_offset_position("left")
    cb3.update_ticks.()
    if save_fig_name != "" # If a save_fig_name is specified, save the plot
        fig.savefig(save_fig_name,dpi=300,bbox_inches="tight") # Use tight border-box, for nice plot. Dpi defines figure quality
    end
    return cmax
end

"""
    plot_orbit_distribution(grid, F)
    plot_orbit_distribution(-||-, save_fig_name="myfile.png",plot_font_size="small", vert=true, cmax=[1e16,2e16,3e16], fig=myfig, ax=myax)

Plots the orbit distribution (or orbit weight function) projected onto 2D planes in the 3 orthogonal directions.
If its3D=true, then the mapping to 3D will not be performed, because the provided input 'f' is already 3D. 
"""
function plot_orbit_distribution(grid::OrbitGrid, F::Union{Array{Float64,3}, Vector{Float64}}; kwargs... )

    if length(size(F))==1
        F_os_3D = map_orbits(grid,F)
    else
        F_os_3D = F
    end

    return plot_orbit_distribution(grid.energy, grid.pitch, grid.r, F_os_3D; kwargs... )

end