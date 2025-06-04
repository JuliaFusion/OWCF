################################ createCustomMagneticEquilibrium.jl ######################################################
# This script will create a Solov'ev magnetic equilibrium given the inputs. If no wall is specified (nothing), the script 
# will assume the wall to be given by the Equilibrium.boundary() function (please see the solovev.jl file of the 
# Equilibrium.jl package, likely located at C:/Users/[your name]/.julia/packages/Equilibrium/) or type 'using Pkg; 
# Pkg.add("https://github.com/JuliaFusion/Equilibrium.jl"); using Equilibrium; ?Equilibrium.boundary()). Basically, the wall 
# will then be right outside the plasma given by the Solov'ev equilibrium. If a .jld2 file is specified, the script will 
# just try to load the 'R' and 'z' data and use that as wall. It is assumed that 'R' and 'z' are vectors of Float64, where each 
# major radius coordinate in 'R' has a corresponding vertical coordinate in 'z'. The units are assumed to 
# be meters.

### The inputs are as follows:
# Please see the OWCF/templates/start_createCustomMagneticEquilibrium_template.jl file for information.

### The output file will be named
# solovev_equilibrium_[DATE].jld2
# And it will have the fields
#   S - The Solov'ev Equilbrium.jl object
#   wall - The wall data

### Other
# To find a suitable total plasma β-value, you can try creating a Solov'ev equilibrium with a trial α value.
# You can then check the resulting total β-value by doing
# 
# using JLD2
# using Equilibrium
# myfile = jldopen("solovev_equilibrium_[DATE].jld2",false,false,false,IOStream)
# S = myfile["S"]
# close(myfile)
# βt = S.beta_t # Toroidal plasma β
# βp = S.beta_p # Poloidal plasma β
# β = βt*βp / (βp+βt)
#
# You can then increase/decrease your α value, until you find the desired total plasma β-value.

# Script written by Henrik Järleblad. Last maintained 2025-04-14.
##########################################################################################################################

## --------------------------------------------------------------------------
# Loading Julia packages
verbose && println("Loading Julia packages... ")
using Equilibrium
using JLD2
using Dates
plot_equilibrium && (using Plots)

## --------------------------------------------------------------------------
# Creating Solov'ev equilibrium from inputs
verbose && println("Creating Solov'ev equilibrium... ")
x_point = (diverted ? (R0*(1-1.1*δ*ϵ),-R0*1.1*κ*ϵ) : nothing)
S = Equilibrium.solovev(B0, R0, ϵ, δ, κ, α, qstar; B0_dir = B0_dir, Ip_dir = Ip_dir, diverted = diverted, x_point = x_point, symmetric = isnothing(x_point))

## --------------------------------------------------------------------------
# Checking wall file for data. If not available, create default wall for Solov'ev equilibrium
if !(isfile(filepath_wall))
    verbose && println("No wall data found in 'filepath_wall' (or could not be loaded). Creating default tokamak wall... ")
    wall = Equilibrium.boundary(S)
else
    verbose && println("Loadable file detected at "*filepath_wall*"... ")
    verbose && print("Attempting to load wall data... ")
    myfile = jldopen(filepath_wall,false,false,false,IOStream)
    R_array = myfile["R"]
    z_array = myfile["z"]
    close(myfile)
    verbose && println("Success!")
    verbose && println("Creating wall from data... ")
    wall = Equilibrium.Boundary(R_array, z_array)
end

## --------------------------------------------------------------------------
# Determine output file name
verbose && println("Determining output file name... ")
date_and_time = split("$(Dates.now())","T")[1]*"at"*split("$(Dates.now())","T")[2][1:5]

if !(filename_out=="") # If the 'filename_out' input variable has been specified... 
    filepath_output_orig = folderpath_out*filename_out # Use the 'filename_out' input variable to name the output data file
else # Otherwise, if the 'filename_out' input variable was left unspecified (default), use the default file name format
    filepath_output_orig = folderpath_out*"solovev_equilibrium_"*date_and_time
end
filepath_output = deepcopy(filepath_output_orig)

count = 1
while isfile(filepath_output*".jld2") # To take care of not overwriting files. Add _(1), _(2) etc
    global filepath_output; global count # Declare global scope of specified variables
    filepath_output = filepath_output_orig*"_($(Int64(count)))"
    count += 1 # global scope, to surpress warnings
end

## --------------------------------------------------------------------------
# Plot the magnetic equilibrium (and wall), if requested
if plot_equilibrium
    flux_r = range(extrema(wall.r)...,length=100)
    flux_z = range(extrema(wall.z)...,length=100)
    inds = CartesianIndices((length(flux_r),length(flux_z)))
    psi_rz = [S(flux_r[ind[1]], flux_z[ind[2]]) for ind in inds]
    psi_mag, psi_bdry = psi_limits(S)

    wall_dR = maximum(wall.r)-minimum(wall.r)
    plot_font = "Computer Modern"
    Plots.default(fontfamily=plot_font)
    plt_crs = Plots.contour(flux_r,flux_z,psi_rz',levels=collect(range(psi_mag,stop=psi_bdry,length=5)),color=:gray, α=0.75, linewidth=1.5, label="",colorbar=false)
    plt_crs = Plots.plot!(wall.r,wall.z,label="Tokamak first wall",linewidth=2.5,color=:black)
    plt_crs = Plots.scatter!([magnetic_axis(S)[1]],[magnetic_axis(S)[2]],label="Mag. axis",markershape=:xcross,markercolor=:black,markerstrokewidth=4)
    plt_crs = Plots.plot!(aspect_ratio=:equal,xlabel="R [m]",ylabel="z [m]", xlims=(minimum(wall.r)-0.1*wall_dR,maximum(wall.r)+wall_dR))
    plt_crs = Plots.plot!(xtickfontsize=14,ytickfontsize=14,xguidefontsize=16,yguidefontsize=16)
    plt_crs = Plots.plot!(legend=:bottomright,legendfontsize=13)
    plt_crs = Plots.plot!(title="Mag. equil. (pol. proj) $(date_and_time)",titlefontsize=14)
    plt_crs = Plots.plot!(dpi=200)
    display(plt_crs)

    if save_equil_plot
        verbose && println("Saving magnetic equilibrium plot in .png file format... ")
        png(plt_crs, filepath_output)
    end
end

## --------------------------------------------------------------------------
# Save the data
filepath_output = filepath_output*".jld2"
myfile = jldopen(filepath_output,true,true,false,IOStream)
write(myfile,"S",S)
write(myfile,"wall",wall)
close(myfile)
verbose && println("Saved output file at "*filepath_output)
println("~~~createCustomMagneticEquilibrium.jl completed successfully!~~~")