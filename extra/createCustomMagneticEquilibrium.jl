################################ createCustomMagneticEquilibrium.jl ######################################################
# This script will create a Solov'ev magnetic equilibrium given the inputs. If no wall is specified, the script 
# will assume the wall to be given by the Equilibrium.boundary() function (please see the solovev.jl file of the 
# Equilibrium.jl package, likely located at C:/Users/[your name]/.julia/packages/Equilibrium/) or type 'using Pkg; 
# Pkg.add("https://github.com/JuliaFusion/Equilibrium.jl"); using Equilibrium; ?Equilibrium.boundary()). Basically, the wall 
# will then be right outside the plasma given by the Solov'ev equilibrium. Please see the 
# OWCF/templates/start_createCustomMagneticEquilibrium_template.jl template for more information.

### The inputs are as follows:
# Please see the OWCF/templates/start_createCustomMagneticEquilibrium_template.jl file for information.

### The output file will be named
# solovev_equilibrium_[DATE AND TIME].jld2
# And it will have the keys
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

# Script written by Henrik Järleblad. Last maintained 2026-09-24.
##########################################################################################################################

## --------------------------------------------------------------------------
# Loading Julia packages
verbose && println("Loading Julia packages... ")
using Equilibrium
using JLD2
using Dates
using ProgressMeter
plot_equilibrium && (using Plots)
debug = debug
date_and_time = split("$(Dates.now())","T")[1]*"at"*split("$(Dates.now())","T")[2][1:5]

## --------------------------------------------------------------------------
# Creating Solov'ev equilibrium from inputs
verbose && println("Creating Solov'ev equilibrium (qstar=$(round(qstar,sigdigits=4)))... ")
if β isa Number
    verbose && println("---> By using β=$(β)")
    x_point = (diverted ? (R0*(1-1.1*δ*ϵ),-R0*1.1*κ*ϵ) : nothing)
    alpha_array = collect(range(-10.0, stop=10.0, length=1001)) # Array of possible α values
    beta_array = zeros(length(alpha_array))
    @showprogress desc="Finding optimal α given β=$(β)" for (i,alpha) in enumerate(alpha_array)
        local S
        S = Equilibrium.solovev(B0, R0, ϵ, δ, κ, alpha, qstar; B0_dir = B0_dir, Ip_dir = Ip_dir, diverted = diverted, x_point = x_point, symmetric = isnothing(x_point))
        βt = S.beta_t # Toroidal plasma β
        βp = S.beta_p # Poloidal plasma β
        beta_array[i] = βt*βp / (βp+βt) # Total plasma β
    end
    noNaN_inds = findall(x-> !isnan(x), beta_array)
    debug && println("DEBUG - Num o NaNs: $(length(beta_array) - length(noNaN_inds))")
    alpha_array = alpha_array[noNaN_inds] # Remove NaNs
    beta_array = beta_array[noNaN_inds] # Remove NaNs
    diff_beta_array = abs.(beta_array .- β)
    p = sortperm(diff_beta_array)
    best_betas = beta_array[p][1:5] # The elements in beta_array closest to β, in sorted order. Keep only the 5 best points
    best_alphas = alpha_array[p][1:5] # The corresponding elements in alpha_array. Keep only the 5 best points
    best_beta = best_betas[1] # The value in beta_array closest to β
    best_alpha = best_alphas[1] # The corresponding alpha value
    if debug
        println("DEBUG ------> Plotting β(α)... ")
        myplt = Plots.plot(alpha_array, beta_array, xlabel="α", ylabel="β", label="", linewidth=2.5)
        myplt = Plots.vline!(myplt, [best_alpha], linestyle=:dash, color=:gray, label="")
        myplt = Plots.hline!(myplt, [β], linestyle=:dash, color=:gray, label="β=$(β)")
        myplt = Plots.scatter!(myplt, vcat(best_alphas), vcat(best_betas), markershape=:xcross, label="β(α) approx. $(β)")
        myplt = Plots.scatter!(myplt, [best_alpha], [best_beta], label="β($(round(best_alpha,sigdigits=3))) = $(round(best_beta, sigdigits=3))", ylims=(0,1))
        println("DEBUG - Saving β(α) plot in .png file format... ")
        png(myplt, folderpath_out*"createCustomMagneticEquilibrium_beta_of_alpha_"*date_and_time)
    end
    S = Equilibrium.solovev(B0, R0, ϵ, δ, κ, best_alpha, qstar; B0_dir = B0_dir, Ip_dir = Ip_dir, diverted = diverted, x_point = x_point, symmetric = isnothing(x_point))
else
    verbose && println("---> By using α=$(α)")
    x_point = (diverted ? (R0*(1-1.1*δ*ϵ),-R0*1.1*κ*ϵ) : nothing)
    S = Equilibrium.solovev(B0, R0, ϵ, δ, κ, α, qstar; B0_dir = B0_dir, Ip_dir = Ip_dir, diverted = diverted, x_point = x_point, symmetric = isnothing(x_point))
end

## --------------------------------------------------------------------------
# Checking wall file for data. If not available, create default wall for Solov'ev equilibrium
if !(isfile(filepath_wall))
    if lowercase(filepath_wall)=="jet"
        verbose && println("Using JET wall geometry as wall data... ")
        _, wall = read_geqdsk(folderpath_OWCF*"equilibrium/JET/g99971/g99971_474-48.9.eqdsk", clockwise_phi=false)
    elseif lowercase(filepath_wall)=="iter"
        verbose && println("Using ITER wall geometry as wall data... ")
        _, wall = read_geqdsk(folderpath_OWCF*"equilibrium/ITER/test/80MW_equilibrium.geqdsk", clockwise_phi=false)
    else
        verbose && println("No wall data found in 'filepath_wall' (or could not be loaded). Creating default tokamak wall... ")
        wall = Equilibrium.boundary(S)
    end
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

if !(filename_out=="") # If the 'filename_out' input variable has been specified... 
    filepath_output_orig = folderpath_out*filename_out # Use the 'filename_out' input variable to name the output data file
else # Otherwise, if the 'filename_out' input variable was left unspecified (default), use the default file name format
    filepath_output_orig = folderpath_out*"createCustomMagneticEquilibrium_solovev_"*date_and_time
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

    verbose && println("Saving magnetic equilibrium plot in .png file format... ")
    png(plt_crs, filepath_output)
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