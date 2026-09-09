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

# ╔═╡ d8114058-2553-46cf-a28e-a4920c3501b9
begin
	# SPECIFY THE INPUTS IN THIS CELL

	# Please specify the OWCF folder and let the cell below change directory to the 
	# OWCF folder when the cell is run. This is to be able to load the
	# correct versions of the Julia packages as specified in the Project.toml and 
	# Manifest.toml files.
	folderpath_OWCF = "REPLACE-THIS-TEXT-WITH-THE-PATH-TO-THE-OWCF-FOLDER-ON-YOUR-COMPUTER" # Finish with '/'
	
	m = 2.0*1.660539e-27 # Mass in kg. Change manually here. Example value alpha particle: 4.001506179127 .* 1.660539e-27
	q = 2.0 * (1.602176e-19) # Charge in Coulombs. Change manually here. Example value alpha particle: 2.0 * (1.602176e-19)
	Emin = 10.0 # keV
	Emax = 2000.0 # KeV
	Bmin = 2.0 # Minimum magnetic field value in the tokamak. Tesla. Example JET approx.: 2.0
	verbose = true
end

# ╔═╡ f271690d-fa0a-426f-8c6c-6ea59ea15966
begin
	cd(folderpath_OWCF)
	using Pkg
	Pkg.activate(".")
	verbose && println("Loading packages... ")
	using GuidingCenterOrbits
	using Interact
	using Plots
	using LinearAlgebra
	using Contour
	using PlutoUI
end

# ╔═╡ 3bcc3ca4-137b-4ba6-995c-81088feb065e
md"""
# comWebApp

## Description:
This notebook provides an application to visualize orbits in a tokamak in an interactive and intuitive manner. The orbits are visualized as a function of the classical constants-of-motion (COM). They are energy (E), magnetic moment (mu) and toroidal canonical momentum (Pphi). The degeneracy for certain combination of EmuPphi is not adressed. Thus, for some EmuPphi-coordinates there will be two orbits plotted. The magnetic equilibrium is created from scratch with the help of the Solov'ev equilibrium formula from the book by Goedbloed. The parameters that can be specified in real-time are the inverse aspect ratio, plasma elongation, triangularity, minor radius (in meters) and Bϕ on axis.

The app should be shown as an interactivate plot at the bottom of this notebook. If it is not, please post an issue on the Github www.github.com/JuliaFusion/OWCF.
- The COM-triplet (E,mu,Pphi) can be interactively varied via the slider above the plot. The plot should then update automatically in real-time.
- The magnetic equilibrium parameters can be interactivaley varied via the number fields above the COM-sliders. When they are changed, the plot should update automatically in real-time.

### Please note
Of all the apps in the OWCF, this might be the least developed one. It works, but there are some (hard-coded) lines into the application the cause limitations. For example, there are are hard-coded limit values for the mu and Pphi sliders.

## Inputs:
- m - Mass of particle. In kg. - Float64
- q - Charge of particle. In Coloumb. - Float64
- Emin - The minimum fast-ion energy value to be accessible. keV - Float64
- Emax - The maximum fast-ion energy value to be accessible. keV - Float64
- verbose - If set to true, the app will talk a lot! - Bool

## Outputs
#### -

## Saved files
#### -

## Notebook written by Henrik Järleblad, henrikj@dtu.dk
## Last maintained 2026-08-26
"""

# ╔═╡ ed858360-360c-42c0-a945-a7e2b114b61a
begin
	verbose && println("Defining μ function... ")
    μ_func(E, B, Pϕ, Ψ, RBϕ) = begin
        res = E/B - (B/(2*m)) * ((Pϕ-q*Ψ)/RBϕ)^2
        (res > 0.0) ? res : 0.0
    end
    verbose && println("Done!")
end

# ╔═╡ 7b7ebbca-4452-45c3-b65f-91e0fbd316f3
begin
	## Approximate ranges. Need better method of determining them.
	println("Calculating range of E, μ and Pϕ values... ")
	E_array = round.(collect(range(Emin,stop=Emax, length=500)),sigdigits=3) # keV
	vmax = sqrt(maximum(E_array)*1000*q/m) # 1000*q is to convert from keV to Joule. Then sqrt(E/m) to obtain speed.
	μ_max = (m*(vmax^2))/(2*minimum(Bmin))
	μ_array = round.(collect(range(0.0,stop=μ_max,length=50)),sigdigits=3)
	
	Pϕ_max = q * 2.6 + ((1/0.4)+1.0) * m * vmax # 2.6 is assumed for Ψ_LCFS. 
	Pϕ_min = -((1/0.4)+1.0) * m * vmax # To be optimized
	Pϕ_array = round.(collect(range(Pϕ_min,stop=Pϕ_max,length=50)),sigdigits=3)
	
	println("Pϕ_min: $(Pϕ_min)")
	println("Pϕ_max: $(Pϕ_max)")
end

# ╔═╡ 6f6eb7cb-a799-460e-b609-050e21062ca5
function COM_input()
	
	return PlutoUI.combine() do Child
		
		inputs = [
			md""" Energy (keV): $(
				Child("E", Slider(E_array, default=900))
			)""",
			md""" μ (A*m^2): $(
				Child("mu", Slider(μ_array, default=4.25e-14))
			)""",
			md""" Pϕ (kg*m^2/s): $(
				Child("pphi", Slider(Pϕ_array, default=5.09e-20))
			)"""
		]
		
		md"""
		#### Constants-of-motion
		$(inputs)
		"""
	end
end

# ╔═╡ 8a4f0c9b-9af1-46c9-9f96-5ca45b6e305f
function plasma_input()

	return PlutoUI.combine() do Child

		inputs=[
			md""" Vacuum Bϕ on-axis (T): $(
				Child("B0", NumberField(0.1:0.1:10.0, default=2.6))
			)""",
			md""" Minor radius (m): $(
				Child("a", NumberField(0.1:0.1:5.0, default=1.0))
			)""",
			md""" Triangularity (-): $(
				Child("triangularity", NumberField(-5.0:0.1:5.0, default=2.0))
			)""",
			md""" Plasma elongation (-): $(
				Child("kappa", NumberField(0.1:0.1:5.0, default=1.4))
			)""",
			md""" Inverse aspect ratio (-): $(
				Child("a_over_R0", NumberField(0.05:0.05:1.0, default=0.4))
			)"""
		]

		md"""
		#### Plasma
		$(inputs)
		"""
	end
end

# ╔═╡ 30636474-2a78-4cf9-8f20-4f919b8d35fd
@bind plasma plasma_input()

# ╔═╡ bf30c698-30b2-4040-b73e-6877d7aa72ea
@bind com COM_input()

# ╔═╡ db3988a1-9f33-415d-aa01-f270ca1cb905
let
    E = com.E
    μ = com.mu
    Pϕ = com.pphi
    plasma_elongation = plasma.kappa
    minor_radius = plasma.a
    inverse_aspect_ratio = plasma.a_over_R0
    triangularity = plasma.triangularity
    vacuum_Bϕ_on_axis = plasma.B0
    
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
    μ_matrix = map(μ_func, E_matrix, B_matrix, Pϕ_matrix, Ψ_matrix, RBϕ_matrix)

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

    Plots.plot(plt_crs)
end

# ╔═╡ Cell order:
# ╠═3bcc3ca4-137b-4ba6-995c-81088feb065e
# ╠═d8114058-2553-46cf-a28e-a4920c3501b9
# ╠═f271690d-fa0a-426f-8c6c-6ea59ea15966
# ╠═ed858360-360c-42c0-a945-a7e2b114b61a
# ╠═7b7ebbca-4452-45c3-b65f-91e0fbd316f3
# ╠═6f6eb7cb-a799-460e-b609-050e21062ca5
# ╠═8a4f0c9b-9af1-46c9-9f96-5ca45b6e305f
# ╠═30636474-2a78-4cf9-8f20-4f919b8d35fd
# ╠═bf30c698-30b2-4040-b73e-6877d7aa72ea
# ╠═db3988a1-9f33-415d-aa01-f270ca1cb905
