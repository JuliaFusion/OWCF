################################ comWebApp.jl #########################################

#### Description:
# This script provides an application to visualize orbits in a tokamak in an interactive
# and intuitive manner. The orbits are visualized as a function of the classical constants
# of motion (COM). They are energy (E), magnetic moment (mu) and toroidal canonical momentum (Pphi).
# The degeneracy for certain combination of EmuPphi is not adressed. Thus, for some EmuPphi-coordinates
# there will be two orbits plotted.
#
# The magnetic equilibrium is created from scratch with the help of the Solov'ev equilibrium formula
# from the book by Goedbloed. The parameters that can be specified in real-time are inverse aspect ratio,
# plasma elongation, triangularity, minor radius (in meters) and Bϕ on axis.
#
# Please note. Of all the apps in the OWCF, this might be the least developed one. It works, but 
# there are some (hard-coded) lines into the application the cause limitations. For example,
# on line 80 and 84 are hard-coded initial limit values for the mu and Pphi sliders.

#### Inputs (units given when defined in the script):
# m - Mass of particle. In kg. - Float64
# q - Charge of particle. In Coloumb. - Float64
# Emin - The minimum fast-ion energy value to be accessible. keV - Float64
# Emax - The maximum fast-ion energy value to be accessible. keV - Float64
# verbose - If set to true, the app will talk a lot! - Bool
# port - The I/O through which the app will be accessed via a web browser - Int64

#### Outputs
# -

#### Saved files
# -

# Script written by Henrik Järleblad. Last maintained 2022-08-25.
#########################################################################################

## --------------------------------------------------------------------------
# Please specify the OWCF folder and let the app change directory to the 
# OWCF folder when comWebApp.jl is executed. This is to be able to load the
# correct versions of the Julia packages as specified in the Project.toml and 
# Manifest.toml files.
folderpath_OWCF = "" # Finish with '/'
cd(folderpath_OWCF)
using Pkg
Pkg.activate(".")

## ------
# Inputs
m = 0.0 # Mass in kg. Change manually here. Example value alpha particle: 4.001506179127 .* 1.660539e-27
q = 0.0 # Charge in Coulombs. Change manually here. Example value alpha particle: 2.0 * (1.602176e-19)
Emin = 0.0 # keV
Emax = 0000.0 # KeV
Bmin = 2.0 # Minimum magnetic field value in the tokamak. Tesla. Example JET approx.: 2.0
verbose = true
port = 18888

## ------
# Loading packages
verbose && println("Loading packages... ")
using Interact
using Plots
using LinearAlgebra
using Mux
using WebIO
using Contour

## ------
verbose && println("Defining μ function... ")
μ_func(E, B, Pϕ, Ψ, RBϕ) = begin
    res = E/B - (B/(2*m)) * ((Pϕ-q*Ψ)/RBϕ)^2
    (res > 0.0) ? res : 0.0
end
verbose && println("Done!")

## ------
## Approximate ranges. Need better method of determining them.
println("Calculating range of E, μ and Pϕ values... ")
E_array = round.(collect(range(Emin,stop=Emax, length=500)),sigdigits=3) # keV
vmax = sqrt(maximum(E_array)*1000*q/m) # 1000*q is to convert from keV to Joule. Then sqrt(E/m) to obtain speed.
Bmin = Bmin
μ_max = (m*(vmax^2))/(2*minimum(Bmin))
μ_array = round.(collect(range(0.0,stop=μ_max,length=50)),sigdigits=3)

Pϕ_max = q * 2.6 + ((1/0.4)+1.0) * m * vmax # 2.6 is assumed for Ψ_LCFS. 
Pϕ_min = -((1/0.4)+1.0) * m * vmax # To be optimized
Pϕ_array = round.(collect(range(Pϕ_min,stop=Pϕ_max,length=50)),sigdigits=3)

println("Pϕ_min: $(Pϕ_min)")
println("Pϕ_max: $(Pϕ_max)")

## ------
# The web application
verbose && println("Building web application... ")

verbose && println("--- You can access the comWebApp via an internet web browser when you see 'Task (runnable)...' ")
verbose && println("--- When 'Task (runnable)...' has appeared, please visit the website localhost:$(port) ---")
verbose && println("--- Remember: It might take a minute or two to load the webpage. Please be patient. ---")
function app(req)
    @manipulate for inverse_aspect_ratio=0.4, plasma_elongation=1.4, triangularity=2.0, minor_radius=1.0, vacuum_Bϕ_on_axis=2.6, E=E_array, μ=μ_array, Pϕ=Pϕ_array

        # Re-written from Mirko Salewski's solovev.m MATLAB script
        maxx=1.2   #plasma in x in [-1;1], plot a little wider
        cellsize=0.02
        x=-maxx:cellsize:maxx
        y=-plasma_elongation*maxx:cellsize:plasma_elongation*maxx
        X_mesh = x' .* ones(length(y)) # Meshgrid
        Y_mesh = ones(length(x))' .* y # Meshgrid

        R0 = minor_radius/inverse_aspect_ratio

        Rvec=R0*ones(length(x)) .+ minor_radius .* collect(x)
        zvec=minor_radius .* collect(y)

        R_mesh = Rvec' .*ones(length(zvec)) # Meshgrid
        z_mesh = ones(length(Rvec))' .* zvec # Meshgrid

        # CALCULATE PSI. The Solov'ev formula.
        Ψ_norm = (X_mesh - 0.5*inverse_aspect_ratio .*(ones(size(X_mesh)) - X_mesh.^2)).^2 +(1 - 0.25*inverse_aspect_ratio^2) .*(ones(size(X_mesh)) + inverse_aspect_ratio*triangularity .* X_mesh .*(2 .*ones(size(X_mesh)) + inverse_aspect_ratio .*X_mesh)) .* (Y_mesh ./plasma_elongation).^2
        Ψ_scale = (minor_radius^2)*vacuum_Bϕ_on_axis /(2*pi)
        Ψ_matrix = Ψ_scale .* Ψ_norm

        # CALCULATE MAGNETIC FIELD
        dψdR = zeros(size(Ψ_matrix))
        dψdR[1:end-1,:] .= diff(Ψ_matrix,dims=1) # dims=1 => R
        dψdR[end,:] .= dψdR[end-1,:] # Set last row to same as second last
        dψdR = (1/(abs(Rvec[2]-Rvec[1]))) .* dψdR # Scale
        dψdz = zeros(size(Ψ_matrix))
        dψdz[:,1:end-1] .= diff(Ψ_matrix,dims=2) # dims=2 => z
        dψdz[:,end] .= dψdz[:,end-1] # Set last column to same as second last
        dψdz = (1/(abs(zvec[2]-zvec[1]))) .* dψdz # Scale

        BR = (-1) .* dψdz ./ R_mesh
        Bz = dψdR ./ R_mesh
        RBϕ_matrix = (vacuum_Bϕ_on_axis * R0) .* ones(size(R_mesh))
        Bϕ = RBϕ_matrix ./ R_mesh
        B_matrix = sqrt.(BR.^2 + Bϕ.^2 + Bz.^2)

        Pϕ_matrix = Pϕ .* ones(size(B_matrix))
        E_matrix = (E*1000*(GuidingCenterOrbits.e0)) .* ones(size(B_matrix)) # 1000*q is to convert from keV to Joule
        print("Extrema of B_matrix: "); verbose && println(extrema(B_matrix))
        print("Extrema of Ψ_matrix: "); verbose && println(extrema(Ψ_matrix))
        μ_matrix = map(μ_func, E_matrix, B_matrix, Pϕ_matrix, Ψ_matrix, RBϕ_matrix)
        print("Extrema of μ_matrix: "); verbose && println(extrema(μ_matrix))
        verbose && println("")

        #cross-sectional plot
        cl = Contour.contour(Rvec,zvec,μ_matrix',μ)
        l = Contour.lines(cl)
        good_Rms = []
        good_zms = []
        for li=1:length(l)
            ll = l[li]
            Rll, zll = Contour.coordinates(ll)
            closed_curve_measure = sqrt((Rll[end]-Rll[1])^2  +(zll[end]-zll[1])^2)
            closed_curve_ref = sqrt((Rll[end]-Rll[end-1])^2  +(zll[end]-zll[end-1])^2)
            if isapprox(closed_curve_measure, closed_curve_ref, atol=1e-1)
                Rmi = argmax(Rll)
                Rm = Rll[Rmi]
                zm = zll[Rmi]
                push!(good_Rms,Rm)
                push!(good_zms,zm)
            end
        end
        plt_crs = Plots.contour(Rvec,zvec,μ_matrix, levels=[μ], label="", aspect_ratio=:equal, xlabel="R [m]", ylabel="z [m]", linewidth=2.5, colorbar=false)
        plt_crs = Plots.contour!(Rvec,zvec,Ψ_matrix, levels=50, label="", linestyle=:dot, color=:grey, colorbar=false)
        plt_crs = Plots.scatter!(good_Rms,good_zms,mc=:red, label="(Rm,zm) points", colorbar=false)

        vbox(vskip(1em),
            hbox(Plots.plot(plt_crs))
        )
    end
end
webio_serve(page("/",app), port)
