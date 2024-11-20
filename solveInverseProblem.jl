###########################################################
# THIS SCRIPT IS CURRENTLY UNDER CONSTRUCTION!!!
# THIS SCRIPT IS CURRENTLY UNDER CONSTRUCTION!!!
# THIS SCRIPT IS CURRENTLY UNDER CONSTRUCTION!!!
#
#
# When completed, it will be able to solve the S=WF inverse
# problem. That is, given a set of measurements S (either from 
# calcSpec.jl or experimental data), use 0th or 1st order 
# Tikhonov together with the weight functions (computed by
# calcOrbWeights.jl or calc2DWeights.jl or even calc4DWeights.jl)
# to reconstruct the fast-ion distribution F. This is done for 
# many regularization parameter values, and returned as a 
# multi-dimensional array.
#
# THIS SCRIPT IS CURRENTLY UNDER CONSTRUCTION!!!
# THIS SCRIPT IS CURRENTLY UNDER CONSTRUCTION!!!
# THIS SCRIPT IS CURRENTLY UNDER CONSTRUCTION!!!
###########################################################

# Load packages
folderpath_OWCF = "/home/henrikj/Codes/OWCF/"
cd(folderpath_OWCF)
using Pkg
Pkg.activate(".")
using JLD2
using HDF5
using FileIO
using Plots
###########################################################
# Load the experimental data, including MPRu, and try with that
# KM14
filepath_KM14_data = "/home/henrikj/Documents/dissemination/papers/99971_2Drec_paper/data/KM14/99971_data_KM14_t=48.5_49.5s.dat"
datafile = open(filepath_KM14_data,"r")
lines = readlines(datafile)
Ed_array_S_KM14 = zeros(length(lines)-1)
spec_KM14 = zeros(length(lines)-1)
spec_err_KM14 = zeros(length(lines)-1)
for i in eachindex(lines)
    if i>1
        line_array = split(lines[i],"\t")
        Ed_array_S_KM14[i-1] = parse(Float64,line_array[2])
        spec_KM14[i-1] = parse(Float64,line_array[3])
        spec_err_KM14[i-1] = parse(Float64,line_array[4])
    end
end
close(datafile)
gi_KM14 = findall(x-> x>=7.75 && x<=9.5,Ed_array_S_KM14)
Ed_array_S_KM14_orig = Ed_array_S_KM14[gi_KM14]
S_exp_KM14_orig = spec_KM14[gi_KM14]
err_exp_KM14_orig = spec_err_KM14[gi_KM14] 
###########################################################
# KM15
filepath_KM15_data = "/home/henrikj/Documents/dissemination/papers/99971_2Drec_paper/data/KM15/99971_data_KM15_t=48.5_49.5s.dat"
datafile = open(filepath_KM15_data,"r")
lines = readlines(datafile)
Ed_array_S_KM15 = zeros(length(lines)-1)
spec_KM15 = zeros(length(lines)-1)
spec_err_KM15 = zeros(length(lines)-1)
for i in eachindex(lines)
    if i>1
        line_array = split(lines[i]," ")
        Ed_array_S_KM15[i-1] = parse(Float64,line_array[2])
        spec_KM15[i-1] = parse(Float64,line_array[3])
        spec_err_KM15[i-1] = parse(Float64,line_array[4])
    end
end
close(datafile)
gi_KM15 = findall(x-> x>=7.75 && x<=9.5,Ed_array_S_KM15)
Ed_array_S_KM15_orig = Ed_array_S_KM15[gi_KM15]
S_exp_KM15_orig = spec_KM15[gi_KM15]
err_exp_KM15_orig = spec_err_KM15[gi_KM15] 
###########################################################
# MPRu
using HDF5
myfile = h5open("/home/henrikj/Documents/dissemination/papers/99971_2Drec_paper/data/MPRu/99971_MPRu_exp_data.h5","r")
Ed_array_S_MPRu = read(myfile["Xpos_array"])
spec_err_MPRu = read(myfile["erro_array"])
spec_MPRu = read(myfile["spec_array"])
close(myfile)
gi_MPRu = findall(x-> x>=20.0 && x<=60.0,Ed_array_S_MPRu)
Ed_array_S_MPRu_orig = Ed_array_S_MPRu[gi_MPRu]
S_exp_MPRu_orig = spec_MPRu[gi_MPRu]
err_exp_MPRu_orig = spec_err_MPRu[gi_MPRu] 
###########################################################
pltKM14nKM15 = Plots.scatter(Ed_array_S_KM14_orig,S_exp_KM14_orig; label="KM14",markershape=:diamond,markersize=4,markercolor=:green)
pltKM14nKM15 = Plots.yerror!(Ed_array_S_KM14_orig,S_exp_KM14_orig; yerror=err_exp_KM14_orig)
pltKM14nKM15 = Plots.plot!(xlabel="Deposited energy [MeV]",ylabel="Counts",yaxis=:identity)
pltKM14nKM15 = Plots.plot!(title="Experimental data",xtickfontsize=13,xguidefontsize=14)
pltKM14nKM15 = Plots.plot!(ytickfontsize=13,yguidefontsize=14)
pltKM14nKM15 = Plots.plot!(titlefontsize=15)
pltKM14nKM15 = Plots.scatter!(Ed_array_S_KM15_orig,S_exp_KM15_orig; label="KM15",markershape=:diamond,markersize=4,markercolor=:red)
pltKM14nKM15 = Plots.yerror!(Ed_array_S_KM15_orig,S_exp_KM15_orig; yerror=err_exp_KM15_orig)

pltMPRu = Plots.scatter(Ed_array_S_MPRu_orig,S_exp_MPRu_orig; label="MPRu",markershape=:diamond,markersize=4,mc=:purple)
pltMPRu = Plots.yerror!(Ed_array_S_MPRu_orig,S_exp_MPRu_orig; yerror=err_exp_MPRu_orig)
pltMPRu = Plots.plot!(xlabel="Proton impact position [cm]",ylabel="Counts",yaxis=:identity)
pltMPRu = Plots.plot!(title="Experimental data",xtickfontsize=13,xguidefontsize=14)
pltMPRu = Plots.plot!(ytickfontsize=13,yguidefontsize=14)
pltMPRu = Plots.plot!(titlefontsize=15)
###########################################################
# Load KM14 weight functions
using JLD2
using HDF5
using FileIO
filepath_W = "/home/henrikj/Codes/OWCF_results/cycle31/17/weights/velWeights_JET_99971L72_at8,9s_KM14_T-d=n-4He.hdf5"
myfile = h5open(filepath_W,"r")
W_KM14 = read(myfile["W"])
Ed_array_W_KM14 = read(myfile["Ed_array"])
Ed_array_W_KM14 = Ed_array_W_KM14 ./(1000.0) # keV to MeV
E_array = read(myfile["E_array"])
p_array = read(myfile["p_array"])
R_of_interest = read(myfile["R"])
z_of_interest = read(myfile["z"])
En_array_W_KM14 = read(myfile["Ed_array_raw"])
close(myfile)
###########################################################
# Define coarse grid vectors. Compute interpolations of 
# weight functions onto coarse grid.
E_coarse = collect(range(minimum(E_array),stop=250.0,length=100))
p_coarse = collect(range(extrema(p_array)...,length=35))
using Interpolations
nodes = (Ed_array_W_KM14, E_array, p_array)

# Create interpolation object
itp = interpolate(nodes, W_KM14, Gridded(Linear()))

W_coarse_KM14 = zeros(length(Ed_array_S_KM14_orig),length(E_coarse),length(p_coarse)) # Pre-allocate 4D array
numObadInds = 0
for Edi=1:length(Ed_array_S_KM14_orig)
    for Ei=1:length(E_coarse)
        for pi=1:length(p_coarse)
            try
                W_coarse_KM14[Edi,Ei,pi] = itp(Ed_array_S_KM14_orig[Edi],E_coarse[Ei],p_coarse[pi]) # Interpolate at 3D query point
            catch
                numObadInds += 1
                debug && println("(Edi: $(Edi), Ei: $(Ei), pi: $(pi)) <--- Interpolation failed for this index") # Print if failed (should not happen)
            end
        end
    end
end
W_2D_KM14_orig = reshape(W_coarse_KM14,(length(Ed_array_S_KM14_orig),length(E_coarse)*length(p_coarse)))
###########################################################
for iEd in axes(W_coarse_KM14,1)
    display(Plots.heatmap(E_coarse,p_coarse,W_coarse_KM14[iEd,:,:]'))
end
###########################################################
# Load the KM15 weight functions
filepath_W_KM15 = "/home/henrikj/Codes/OWCF_results/cycle31/17/weights/velWeights_JET_99971L72_at8,9s_KM15_T-d=n-4He.hdf5"
myfile = h5open(filepath_W_KM15,"r")
W_KM15 = read(myfile["W"])
W_KM15_raw = read(myfile["W_raw"])
Ed_array_W_KM15 = read(myfile["Ed_array"])
Ed_array_W_KM15 = Ed_array_W_KM15 ./(1000.0) # keV to MeV
Ed_array_W_KM15_raw = read(myfile["Ed_array_raw"])
close(myfile)
###########################################################
# Do an equivalent interpolation for the KM15, as for the KM14
nodes = (Ed_array_W_KM15, E_array, p_array)
itp = interpolate(nodes, W_KM15, Gridded(Linear()))
W_coarse_KM15 = zeros(length(Ed_array_S_KM15_orig),length(E_coarse),length(p_coarse)) # Pre-allocate 4D array
numObadInds = 0
for Edi=1:length(Ed_array_S_KM15_orig)
    for Ei=1:length(E_coarse)
        for pi=1:length(p_coarse)
            try
                W_coarse_KM15[Edi,Ei,pi] = itp(Ed_array_S_KM15_orig[Edi],E_coarse[Ei],p_coarse[pi]) # Interpolate at 3D query point
            catch
                numObadInds += 1
                debug && println("(Edi: $(Edi), Ei: $(Ei), pi: $(pi)) <--- Interpolation failed for this index") # Print if failed (should not happen)
            end
        end
    end
end
W_2D_KM15_orig = reshape(W_coarse_KM15,(length(Ed_array_S_KM15_orig),length(E_coarse)*length(p_coarse)))
###########################################################
for iEd in axes(W_coarse_KM15,1)
    display(Plots.heatmap(E_coarse,p_coarse,W_coarse_KM15[iEd,:,:]'))
end
###########################################################
# Compute the MPRu weight functions by multiplying the KM15_raw weight functions with the MPRu instrumental response
# Then, interpolate onto a nice grid
using PyCall
py"""
import numpy as np

En_array_transf = np.loadtxt("/home/henrikj/Codes/OWCF/vc_data/MPRu/En_keV.txt")
transf_mat = np.loadtxt("/home/henrikj/Codes/OWCF/vc_data/MPRu/matrix.txt")
X_array_transf = np.loadtxt("/home/henrikj/Codes/OWCF/vc_data/MPRu/Xaxis_cm.txt")
"""
En_array_transf = py"En_array_transf"
transf_mat = py"transf_mat"
X_array_transf = py"X_array_transf"
lo = findfirst(x-> x>minimum(Ed_array_W_KM15_raw),En_array_transf)
hi = findlast(x-> x<maximum(Ed_array_W_KM15_raw),En_array_transf)
W_MPRu_orig = zeros(length(X_array_transf),length(E_array),length(p_array))
transf_mat_transp = transf_mat'
for (iE,E) in enumerate(E_array), (ip,p) in enumerate(p_array)
    nodes = (Ed_array_W_KM15_raw,)
    itp = interpolate(nodes, W_KM15_raw[:,iE,ip], Gridded(Linear()))
    W_KM15_raw_i = itp.(En_array_transf[lo:hi])
    W_MPRu_orig[:,iE,ip] = (transf_mat_transp[:,lo:hi])*W_KM15_raw_i
end

nodes = (X_array_transf, E_array, p_array)
itp = interpolate(nodes, W_MPRu_orig, Gridded(Linear()))
W_coarse_MPRu = zeros(length(Ed_array_S_MPRu_orig),length(E_coarse),length(p_coarse)) # Pre-allocate 4D array
numObadInds = 0
for Edi=1:length(Ed_array_S_MPRu_orig)
    for Ei=1:length(E_coarse)
        for pi=1:length(p_coarse)
            try
                W_coarse_MPRu[Edi,Ei,pi] = itp(Ed_array_S_MPRu_orig[Edi],E_coarse[Ei],p_coarse[pi]) # Interpolate at 3D query point
            catch
                numObadInds += 1
                debug && println("(Xpos: $(Edi), Ei: $(Ei), pi: $(pi)) <--- Interpolation failed for this index") # Print if failed (should not happen)
            end
        end
    end
end
W_2D_MPRu_orig = reshape(W_coarse_MPRu,(length(Ed_array_S_MPRu_orig),length(E_coarse)*length(p_coarse)))
###########################################################
for iEd in axes(W_coarse_MPRu,1)
    display(Plots.heatmap(E_coarse,p_coarse,W_coarse_MPRu[iEd,:,:]'))
end
###########################################################
# Compute test Maxwellian fast-ion distribution
# Compute test Maxwellian fast-ion distribution
# Compute test Maxwellian fast-ion distribution
E0 = 200.0 # keV 
δE = 50
p0 = 0.5 # -
δp = 0.15
f_test = zeros(length(E_coarse),length(p_coarse))
for (iE,E) in enumerate(E_coarse)
    for (ip,p) in enumerate(p_coarse)
        ΔE = E - E0
        Δp = p - p0
        f_test[iE,ip] = exp(-(ΔE/δE)^2) * exp(-(Δp/δp)^2)
    end
end
maximum(f_test)
f_test = (1.0e19/sum(diff(E_coarse)[1]*diff(p_coarse)[1] .*f_test)) .*f_test
f_test = map(x-> x<0.001*maximum(f_test) ? 0.0 : x, f_test)
maximum(f_test)
###########################################################
# OR, load f_test as the TRANSP distribution, and interpolate onto grid of interest
# OR, load f_test as the TRANSP distribution, and interpolate onto grid of interest
# OR, load f_test as the TRANSP distribution, and interpolate onto grid of interest
include("/home/henrikj/Codes/OWCF/misc/convert_units.jl")
myfile = jldopen("/home/henrikj/Data/JET/TRANSP/99971/L72/99971L72_fi_2.jld2",false,false,false,IOStream)
F_EpRz_TRANSP = myfile["F_ps"]
E_array_TRANSP = myfile["energy"]
p_array_TRANSP = myfile["pitch"]
R_array_TRANSP = myfile["R"]
z_array_TRANSP = myfile["z"]
close(myfile)
iR_TRANSP = argmin(abs.(R_array_TRANSP .- 100*R_of_interest))
iz_TRANSP = argmin(abs.(z_array_TRANSP .- 100*z_of_interest))
f_test_orig = units_conversion_factor("cm^-3","m^-3") .*F_EpRz_TRANSP[:,:,iR_TRANSP,iz_TRANSP] # TRANSP distribution is in cm^-3. Convert to m^-3
nodes_TRANSP = (E_array_TRANSP,p_array_TRANSP)
itp = Interpolations.interpolate(nodes_TRANSP,f_test_orig,Gridded(Linear()))
etp = Interpolations.extrapolate(itp,Interpolations.Flat())
f_test = zeros(length(E_coarse),length(p_coarse))
for (iE,E) in enumerate(E_coarse)
    for (ip,p) in enumerate(p_coarse)
        f_test[iE,ip] = etp(E,p)
    end
end
###########################################################
# Plot test distribution (Maxwellian)
Plots.heatmap(E_coarse,p_coarse,f_test',title="Ground truth")
###########################################################
# Compute and plot synthetic signal
f_1D = reshape(f_test,(length(E_coarse)*length(p_coarse),1))
S_synth_KM14_clean = dropdims(W_2D_KM14_orig * f_1D,dims=2)
Plots.plot(Ed_array_S_KM14_orig,S_synth_KM14_clean)
###########################################################
# Add noise to the signal. Superimpose it on signal plot
k = 0.1
b = 0.05
include("extra/dependencies.jl")
S_synth_KM14_noisy, err_synth_KM14 = add_noise(S_synth_KM14_clean, b; k=k)
Plots.plot!(Ed_array_S_KM14_orig,S_synth_KM14_noisy)
###########################################################
S_synth_KM15_clean = dropdims(W_2D_KM15_orig * f_1D,dims=2)
Plots.plot!(Ed_array_S_KM15_orig, S_synth_KM15_clean)
S_synth_KM15_noisy, err_synth_KM15 = add_noise(S_synth_KM15_clean, b; k=k)
Plots.plot!(Ed_array_S_KM15_orig,S_synth_KM15_noisy)
###########################################################
S_synth_MPRu_clean = dropdims(W_2D_MPRu_orig * f_1D,dims=2)
Plots.plot(Ed_array_S_MPRu_orig,S_synth_MPRu_clean)
S_synth_MPRu_noisy, err_synth_MPRu = add_noise(S_synth_MPRu_clean, b; k=k)
Plots.plot!(Ed_array_S_MPRu_orig,S_synth_MPRu_noisy)
###########################################################
# Remnant code from earlier version
function nonunique2(x::AbstractArray{T}) where T
    xs = sort(x)
    duplicatevector = T[]
    for i=2:length(xs)
        if (isequal(xs[i],xs[i-1]) && (length(duplicatevector)==0 || !isequal(duplicatevector[end],xs[i])))
            push!(duplicatevector,xs[i])
        end
    end
    duplicatevector
end
###########################################################
using NetCDF
"""
read_ncdf(filepath; vars=nothing)

# Function description here
"""
function read_ncdf(filepath::String; wanted_keys=nothing)

    d = Dict()
    d["err"] = 1
    if isfile(filepath)
        d["err"] = 0
        NetCDF.open(filepath) do nc
            cdf_variables = nc.vars
            if !(wanted_keys==nothing)
                for wanted_key in wanted_keys
                    if wanted_key in keys(cdf_variables)
                        values = NetCDF.readvar(nc,wanted_key)
                        if () == size(values) # If the size of values is 0 (i.e. values is a scalar)
                            d[wanted_key] = first(values) # Take the 'first' element, parse a float and store it in d with key 'wanted_key'
                        else
                            if typeof(values) <: Vector{NetCDF.ASCIIChar}
                                values = reduce(*,map(x-> "$(x)",values)) # Concatenate all ASCII characters into one string
                            end
                            d[wanted_key] = values # Parse all elements in values as floats and store them as an array in d with key 'wanted_key'
                        end
                    end
                end
            else
                for (key,_) in cdf_variables
                    values = NetCDF.readvar(nc,key)
                    if () == size(values) # If the size of values is 0 (i.e values is a scalar)
                        d[key] = first(values) # Take the 'first' element, parse a float and store it in d with key 'wanted_key'
                    else
                        if typeof(values) <: Vector{NetCDF.ASCIIChar}
                            values = reduce(*,map(x-> "$(x)",values)) # Concatenate all ASCII characters into one string
                        end
                        d[key] = values # Parse all elements in values as floats and store them as an array in d with key 'wanted_key'
                    end
                end
            end
        end
    else
        error("FILE DOES NOT EXIST: "*filepath)
    end
    return d
end
###########################################################
# Create SD basis functions
include("extra/dependencies.jl")
include("misc/temp_n_dens.jl")
my_dict = read_ncdf("/home/henrikj/Data/JET/TRANSP/99971/L72/99971L72_fi_2.cdf",wanted_keys=["TIME"])
t_mid = my_dict["TIME"]
R_of_interest = 2.8989
z_of_interest = 0.27561
M, wall = read_geqdsk("/home/henrikj/Data/JET/eqdsk/99971/g99971_474-48.9.eqdsk",clockwise_phi=false)
psi_axis, psi_bdry = psi_limits(M)
rho_of_interest = sqrt((M(R_of_interest,z_of_interest)-psi_axis)/(psi_bdry-psi_axis))
ne_itp = getDensProfileFromTRANSP(t_mid,"/home/henrikj/Data/JET/TRANSP/99971/L72/99971L72.cdf","e")
nD_itp = getDensProfileFromTRANSP(t_mid,"/home/henrikj/Data/JET/TRANSP/99971/L72/99971L72.cdf","D")
nT_itp = getDensProfileFromTRANSP(t_mid,"/home/henrikj/Data/JET/TRANSP/99971/L72/99971L72.cdf","T")
n_e = ne_itp(rho_of_interest)
n_D = nD_itp(rho_of_interest)
n_T = nT_itp(rho_of_interest)
Te_itp = getTempProfileFromTRANSP(t_mid,"/home/henrikj/Data/JET/TRANSP/99971/L72/99971L72.cdf","e")
TD_itp = getTempProfileFromTRANSP(t_mid,"/home/henrikj/Data/JET/TRANSP/99971/L72/99971L72.cdf","D")
TT_itp = getTempProfileFromTRANSP(t_mid,"/home/henrikj/Data/JET/TRANSP/99971/L72/99971L72.cdf","T")
T_e = Te_itp(rho_of_interest)
T_D = TD_itp(rho_of_interest)
T_T = TT_itp(rho_of_interest)
###########################################################
# Test upgraded Debye length function
λ_D = debye_length(n_e, T_e, ["T","D"], [n_T,n_D], [T_T, T_D])
###########################################################
include("extra/dependencies.jl")
E0_array = E_coarse[1:1:end-1] .+ diff(E_coarse)[1]/2
p0_array = p_coarse[2:1:end-1]
F_SD = zeros(length(E_coarse)*length(p_coarse),length(E0_array)*length(p0_array))
F_SD_c_inds = CartesianIndices((length(E0_array),length(p0_array)))
i_SD = 1
E_oi = 250.0
iE_oi = argmin(abs.(E0_array .- E_oi))
p_oi = 0.0
ip_oi = argmin(abs.(p0_array .- p_oi))
for (ic,c) in enumerate(F_SD_c_inds)
    iE = c[1]; ip = c[2]; E_0 = E0_array[iE]; p_0 = p0_array[ip]
    #println("$(round(100*i_SD/(length(E_coarse)*length(p_coarse)),digits=3)) %")
    f_SD = slowing_down_function(E_0, p_0, E_coarse, p_coarse, n_e, T_e, "D", ["T","D"], [n_T, n_D], [T_T, T_D]; dampen=true, damp_type=:erfc, sigma=2.0)
    #println(extrema(f_SD))
    dE = diff(E_coarse)[1]
    dp = diff(p_coarse)[1]
    f_SD = f_SD ./sum((dE*dp) .*f_SD) # Normalize the basis function so they integrate to 1.0
    F_SD[:,ic] .= reshape(f_SD,(length(E_coarse)*length(p_coarse),1))
    i_SD += 1
    if iE==iE_oi && ip==ip_oi
        my_plt = Plots.heatmap(E_coarse,p_coarse,f_SD',fillcolor=cgrad([:white, :yellow, :orange, :red, :black]))
        display(my_plt)
    end
end
###########################################################
anim = @animate for i in 1:25:size(F_SD,2)
    f_SD = reshape(F_SD[:,i],(length(E_coarse),length(p_coarse)))
    my_plt = Plots.heatmap(E_coarse,p_coarse,f_SD',fillcolor=cgrad([:white, :yellow, :orange, :red, :black]))
    display(my_plt)
end
#gif(anim,"test.gif",fps=2)
#rm("test.gif"; force=true)
###########################################################
function testyMcTestface(x::T) where T<:Real
    return 2*x
end
###########################################################
# Filter out all basis functions that are NOT identically zero and do NOT contain NaNs 
# hcat them back into a matrix
F_SD_new = reduce(hcat,filter(col-> sum(col)!=0.0 && sum(isnan.(col))==0, eachcol(F_SD)))
###########################################################
for i in 1:50:size(F_SD_new,2)
    f_SD_new = reshape(F_SD_new[:,i],(length(E_coarse),length(p_coarse)))
    my_plt = Plots.heatmap(E_coarse,p_coarse,f_SD_new',fillcolor=cgrad([:white, :yellow, :orange, :red, :black]))
    display(my_plt)
end
###########################################################
SD_prior = true
if !SD_prior
    # Without SD prior
    W_SD_KM14_orig = W_2D_KM14_orig
    W_SD_KM15_orig = W_2D_KM15_orig
    W_SD_MPRu_orig = W_2D_MPRu_orig
else
    # With SD prior
    W_SD_KM14_orig = W_2D_KM14_orig*F_SD_new
    W_SD_KM15_orig = W_2D_KM15_orig*F_SD_new
    W_SD_MPRu_orig = W_2D_MPRu_orig*F_SD_new
end
###########################################################
# Rescale weight functions, if necessary
# In final OWCF version, code an intelligent check, to 
# determine if re-scaling is necessary
# It could be a flag like: 
# (car_W = true) && (typical_N_FI = 1.0e19)
# 'car_W' stands for check and rescale weights
# A sample of Maxwellians will then be created. They each have 'typical_N_FI'
# as the value of their zeroth moment. Their signals will be computed and 
# compared with the experimental data. If none of the signals are within an 
# order of magnitude from the maximum of any of the experimental signals, 
# the weight functions should be re-scaled as below
rescale_func = maximum # mean
W_2D_KM14_a = (rescale_func(S_exp_KM14_orig)/rescale_func(S_synth_KM14_clean)) .*W_SD_KM14_orig
W_2D_KM15_a = (rescale_func(S_exp_KM15_orig)/rescale_func(S_synth_KM15_clean)) .*W_SD_KM15_orig
W_2D_MPRu_a = (rescale_func(S_exp_MPRu_orig)/rescale_func(S_synth_MPRu_clean)) .*W_SD_MPRu_orig
###########################################################
# Try using the non-rescaled weight functions on the experimental data
#W_2D_KM14_a = W_SD_KM14_orig#W_2D_KM14_orig
#W_2D_KM15_a = W_SD_KM15_orig#W_2D_KM15_orig
#W_2D_MPRu_a = W_SD_MPRu_orig#W_2D_MPRu_orig
###########################################################
# Prepare experimental data problem first
noise_floor_factor = 1.0e-4
noise_floor_KM14 = maximum(S_exp_KM14_orig)*noise_floor_factor
noise_floor_KM15 = maximum(S_exp_KM15_orig)*noise_floor_factor
noise_floor_MPRu = maximum(S_exp_MPRu_orig)*noise_floor_factor
err_exp_KM14_a = replace(x-> (x<noise_floor_KM14) ? noise_floor_KM14 : x, err_exp_KM14_orig)
err_exp_KM15_a = replace(x-> (x<noise_floor_KM15) ? noise_floor_KM15 : x, err_exp_KM15_orig)
err_exp_MPRu_a = replace(x-> (x<noise_floor_MPRu) ? noise_floor_MPRu : x, err_exp_MPRu_orig)
###########################################################
using Base.Sort

function smallestn(a, n)
  sort(a; alg=Sort.PartialQuickSort(n))[1:n]
end
###########################################################
# Exclude unwanted measurement bins from reconstruction
good_inds_KM14 = findall(x-> !(x>=8.16 && x<=8.56),Ed_array_S_KM14_orig)
good_inds_KM15 = findall(x-> !(x>=8.36 && x<=8.76),Ed_array_S_KM15_orig)
good_inds_MPRu = findall(x-> !(x>= 23.7 && x<= 27.7),Ed_array_S_MPRu_orig)
zero_inds_MPRu = findall(x-> iszero(x), S_exp_MPRu_orig)
good_inds_MPRu = filter(x-> !(x in zero_inds_MPRu), good_inds_MPRu)

W_2D_KM14_b = W_2D_KM14_a[good_inds_KM14,:]
W_2D_KM15_b = W_2D_KM15_a[good_inds_KM15,:]
W_2D_MPRu_b = W_2D_MPRu_a[good_inds_MPRu,:]
err_exp_KM14_b = err_exp_KM14_a[good_inds_KM14]
err_exp_KM15_b = err_exp_KM15_a[good_inds_KM15]
err_exp_MPRu_b = err_exp_MPRu_a[good_inds_MPRu]
S_exp_KM14_a = S_exp_KM14_orig[good_inds_KM14]
S_exp_KM15_a = S_exp_KM15_orig[good_inds_KM15]
S_exp_MPRu_a = S_exp_MPRu_orig[good_inds_MPRu]
# Put together one large weight matrix and measurement vector
W_hat = vcat(W_2D_KM14_b ./err_exp_KM14_b, W_2D_KM15_b ./err_exp_KM15_b, W_2D_MPRu_b ./err_exp_MPRu_b) #W_hat = vcat(W_2D_KM14 ./err_exp_KM14, W_2D_KM15 ./err_exp_KM15, W_2D_MPRu ./err_exp_MPRu)
s_hat = vcat(S_exp_KM14_a ./err_exp_KM14_b, S_exp_KM15_a ./err_exp_KM15_b, S_exp_MPRu_a ./err_exp_MPRu_b)
###########################################################
# Add 0th order Tikhonov
L0 = zeros(size(W_hat,2),size(W_hat,2)) #L0 = zeros(length(f_1D),length(f_1D))
for i=1:size(W_hat,2)#length(f_1D)
    for j=1:size(W_hat,2)#length(f_1D)
        if i==j
            L0[i,j] = 1.0
        end
    end
end
###########################################################
# Add 1st order Tikhonov
#L1 = zeros(length(W_hat),length(W_hat))
#for i=1:length(W_hat)
#    for j=1:length(W_hat)
#        if (i-j)==1
#            L1[i,j] = 1.0
#        end
#        if (i-j)==-1
#            L1[i,j] = -1.0
#        end
#    end
#end
#L1[1,:] .= 0.0
#L1[:,end] .= 0.0
#spy(L1)
###########################################################
#Normalize everything
Whm = maximum(W_hat)
W_hh = W_hat ./Whm
shm = maximum(s_hat)
s_hh = s_hat ./shm
L_orig = L0
f0 = zeros(size(W_hat,2))#zeros(size(L_orig,1))
###########################################################
# SOLVE THE EXPERIMENTAL DATA PROBLEM
using SCS, Convex
using LinearAlgebra
lambda_values = 10 .^(range(-3,stop=3.0,length=20))
x_sols = zeros(length(lambda_values))
p_sols = zeros(length(lambda_values))
F_sols = zeros(length(f0),length(lambda_values))
for (il,lambda) in enumerate(lambda_values)
    println("lambda: $(lambda)")
    L = lambda*L_orig
    x = Convex.Variable(length(f0))

    problem = Convex.minimize(Convex.sumsquares(vcat(W_hh,L) * x - vcat(s_hh,f0)), [x >= 0]) # FOR SD PRIOR, THIS SHOULD BE F_SD_new*x >= 0 !!!

    solve!(problem, SCS.Optimizer)

    x_sols[il] = norm(dropdims(x.value,dims=2))
    p_sols[il] = norm(W_hh*dropdims(x.value,dims=2) - s_hh)
    F_sols[:,il] = (shm/Whm) .*dropdims(x.value,dims=2)
end
###########################################################
# Plot L-curve
Plots.plot(p_sols,x_sols)
Plots.scatter!([p_sols[1]],[x_sols[1]],label="Start")
Plots.scatter!([p_sols[end]],[x_sols[end]],label="End")
###########################################################
using Interpolations
p_itp = Interpolations.interpolate(p_sols, BSpline(Cubic()))
x_itp = Interpolations.interpolate(x_sols, BSpline(Cubic()))
p_prim = map(x-> Vector(Interpolations.gradient(p_itp,x))[1],eachindex(p_sols))
x_prim = map(x-> Vector(Interpolations.gradient(x_itp,x))[1],eachindex(x_sols))
p_prim_itp = Interpolations.interpolate(p_prim, BSpline(Cubic()))
x_prim_itp = Interpolations.interpolate(x_prim, BSpline(Cubic()))
p_primprim = map(x-> Vector(Interpolations.gradient(p_prim_itp,x))[1],eachindex(p_sols))
x_primprim = map(x-> Vector(Interpolations.gradient(x_prim_itp,x))[1],eachindex(x_sols))

Plots.plot(p_itp)
Plots.plot!(p_prim)
Plots.plot!(p_primprim)
Plots.plot(x_itp)
Plots.plot!(x_prim)
Plots.plot!(x_primprim)
gamma = zeros(length(p_sols))
for i in eachindex(p_sols)
    gamma[i] = (x_primprim[i]*p_prim[i] - x_prim[i]*p_primprim[i])/((x_prim[i]*x_prim[i] + p_prim[i]*p_prim[i])^(3/2))
end
ilm = argmax(gamma)
Plots.plot(gamma)
###########################################################
plot_inds = 1:length(lambda_values) #15:33 # [24] #
anim = @animate for i=plot_inds
    Sr_KM14 = W_2D_KM14_b*F_sols[:,i]
    myplt_KM14 = Plots.plot(Ed_array_S_KM14_orig[good_inds_KM14], Sr_KM14,title="$(i): log10(λ)=$(round(log10(lambda_values[i]),sigdigits=4))",label="WF*")
    myplt_KM14 = Plots.scatter!(myplt_KM14, Ed_array_S_KM14_orig,S_exp_KM14_orig ,label="S")
    Sr_KM15 = W_2D_KM15_b*F_sols[:,i]
    myplt_KM15 = Plots.plot(Ed_array_S_KM15_orig[good_inds_KM15], Sr_KM15,title="$(i): log10(λ)=$(round(log10(lambda_values[i]),sigdigits=4))",label="WF*")
    myplt_KM15 = Plots.scatter!(myplt_KM15, Ed_array_S_KM15_orig,S_exp_KM15_orig ,label="S")
    Sr_MPRu = W_2D_MPRu_b*F_sols[:,i]
    myplt_MPRu = Plots.plot(Ed_array_S_MPRu_orig[good_inds_MPRu],Sr_MPRu ,title="$(i): log10(λ)=$(round(log10(lambda_values[i]),sigdigits=4))",label="WF*")
    myplt_MPRu = Plots.scatter!(myplt_MPRu, Ed_array_S_MPRu_orig,S_exp_MPRu_orig ,label="S")
    myplt_tot_1 = Plots.plot(myplt_KM14, myplt_KM15, myplt_MPRu, layout=(1,3))

    if SD_prior
        Fr_2D = reshape(F_SD_new*F_sols[:,i],length(E_coarse),length(p_coarse))#reshape(F_sols[:,i],size(f_test))
    else
        Fr_2D = reshape(F_sols[:,i],length(E_coarse),length(p_coarse))
    end

    normf_max_Fr_2D = floor(log10(maximum(Fr_2D)))
    normf_max_F_sols = floor(log10(maximum(F_sols[:,i])))
    gi = findall(x-> x<250.0,E_coarse)
    myplt = Plots.heatmap(E_coarse[gi],p_coarse,Fr_2D[gi,:]' ./(10^normf_max_Fr_2D),title="$(i): log10(λ)=$(round(log10(lambda_values[i]),sigdigits=4)) [x 10^$(Int(normf_max_Fr_2D))]",fillcolor=cgrad([:white, :yellow, :orange, :red, :black]))
    #myplt0 = Plots.plot(F_sols[:,i] ./(10^normf_max_F_sols))
    myplt0 = Plots.heatmap(E0_array,p0_array,reshape(F_sols[:,i] ./(10^normf_max_F_sols),(length(E0_array),length(p0_array)))',title="$(i): SD coefficients value [x 10^$(Int(normf_max_F_sols))]",fillcolor=cgrad([:white, :darkblue, :green, :yellow, :orange, :red]),xticks=[50,100,150,200],yticks=[-1.0,-0.5,0.0,0.5,1.0],ylims=(-1.0,1.0))
    p0, pN = extrema(p_sols); diffP = abs(pN-p0)
    x0, xN = extrema(x_sols); diffX = abs(xN-x0)
    l_mag = sqrt((p_sols[i]-p0)*(p_sols[i]-p0)/(diffP*diffP)+(x_sols[i]-x0)*(x_sols[i]-x0)/(diffX*diffX))
    myplt1 = Plots.plot([0,p_sols[i]],[0,x_sols[i]],color=:black,title="||L_n||=$(round(l_mag,sigdigits=4))",label="")
    myplt1 = Plots.plot!(p_sols,x_sols,xlims=(minimum(p_sols)-0.05*diffP, maximum(p_sols)+0.05*diffP),ylims=(minimum(x_sols)-0.05*diffX, maximum(x_sols)+0.05*diffX),label="")
    myplt1 = Plots.scatter!(myplt1, [p_sols[1]],[x_sols[1]],label="Start")
    myplt1 = Plots.scatter!(myplt1, [p_sols[end]],[x_sols[end]],label="End")
    myplt1 = Plots.scatter!(myplt1, [p_sols[ilm]],[x_sols[ilm]],label="gamma_max")
    myplt1 = Plots.scatter!(myplt1, [p_sols[i]],[x_sols[i]],label="$(round(log10(lambda_values[i]),sigdigits=4))")
    myplt_tot_2 = Plots.plot(left_margin=8Plots.mm, myplt, myplt0, myplt1,layout=(1,3))
    myplt_tot = Plots.plot(myplt_tot_1, myplt_tot_2, layout=(2,1),size=(1500,900))
    display(myplt_tot)  
    #display(myplt)
end
gif(anim,"test_SD_rescale_with_new_TRANSP_WF.gif",fps=1)
###########################################################
Ni = plot_inds[26]
Fr_2D = reshape(F_SD_new*F_sols[:,Ni],length(E_coarse),length(p_coarse))#reshape(F_sols[:,i],size(f_test))
gi = findall(x-> x<250.0,E_coarse)
E_array_for_saving = E_coarse[gi]
p_array_for_saving = p_coarse
lambda_for_saving = lambda_values[Ni]
F_2D_for_saving = Fr_2D[gi,:]
myfile = jldopen("99971_rec_SD_withRot_c3.jld2",true,true,false,IOStream)
write(myfile,"F",F_2D_for_saving)
write(myfile,"E_array",E_array_for_saving)
write(myfile,"p_array",p_array_for_saving)
write(myfile,"lambda",lambda_for_saving)
close(myfile)
###########################################################
# Redo the reconstruction. But for the synthetic case
# For comparison and sanity check
replace!(x-> (x==0.0) ? 1.0 : x, err_synth_KM14)
replace!(x-> (x==0.0) ? 1.0 : x, err_synth_KM15)
replace!(x-> (x==0.0) ? 1.0 : x, err_synth_MPRu)
# W_2D_KM14, W_2D_KM15 and W_2D_MPRu are f*cked here 
# Because of normalize to fit experimental data
W_hat = vcat(W_SD_KM14_orig ./err_exp_KM14, W_SD_KM15_orig ./err_exp_KM15, W_SD_MPRu_orig ./err_exp_MPRu) #vcat(W_2D_KM14_orig ./err_synth_KM14, W_2D_KM15_orig ./err_synth_KM15, W_2D_MPRu_orig ./err_synth_MPRu)
s_hat = vcat(S_synth_KM14_noisy ./err_synth_KM14, S_synth_KM15_noisy ./err_synth_KM15, S_synth_MPRu_noisy ./err_synth_MPRu)
Whm = maximum(W_hat)
W_hh = W_hat ./Whm
shm = maximum(s_hat)
s_hh = s_hat ./shm
f0 = zeros(size(L_orig,1))
x_sols_synth = zeros(length(lambda_values))
p_sols_synth = zeros(length(lambda_values))
F_sols_synth = zeros(length(f0),length(lambda_values))
for (il,lambda) in enumerate(lambda_values)
    println("lambda: $(lambda)")
    L = lambda*L0
    x = Convex.Variable(length(f0))

    problem = Convex.minimize(Convex.sumsquares(vcat(W_hh,L) * x - vcat(s_hh,f0)), [x >= 0])

    solve!(problem, SCS.Optimizer)

    x_sols_synth[il] = norm(dropdims(x.value,dims=2))
    p_sols_synth[il] = norm(W_hh*dropdims(x.value,dims=2) - s_hh)
    F_sols_synth[:,il] = (shm/Whm) .*dropdims(x.value,dims=2)
end
###########################################################
# Check the synthetic reconstruction
for i=1:length(lambda_values)
    Sr = W_hh*F_sols_synth[:,i]#W_hh*F_sols_synth[:,i]
    myplt = Plots.plot(Sr ./maximum(Sr),title="$(i): log10(λ)=$(log10(lambda_values[i]))",label="ŴF*")
    myplt = Plots.plot!(s_hh ./maximum(s_hh),label="Ŝ")
    display(myplt)
end
for i=1:length(lambda_values)
    Sr_KM14 = W_2D_KM14*F_SD_new*F_sols_synth[:,i]
    myplt = Plots.plot(Ed_array_S_KM14,Sr_KM14 ./maximum(Sr_KM14),title="$(i): log10(λ)=$(log10(lambda_values[i]))",label="WF*")
    myplt = Plots.plot!(Ed_array_S_KM14,S_synth_KM14_noisy ./maximum(S_synth_KM14_noisy),label="S")
    display(myplt)
end
for i=1:length(lambda_values)
    Sr_KM15 = W_2D_KM15*F_SD_new*F_sols_synth[:,i]
    myplt = Plots.plot(Ed_array_S_KM15,Sr_KM15 ./maximum(Sr_KM15),title="$(i): log10(λ)=$(log10(lambda_values[i]))",label="WF*")
    myplt = Plots.plot!(Ed_array_S_KM15,S_synth_KM15_noisy ./maximum(S_synth_KM15_noisy),label="S")
    display(myplt)
end
for i=1:length(lambda_values)
    Sr_MPRu = W_2D_MPRu*F_SD_new*F_sols_synth[:,i]
    myplt = Plots.plot(Ed_array_S_MPRu,Sr_MPRu ./maximum(Sr_MPRu),title="$(i): log10(λ)=$(log10(lambda_values[i]))",label="WF*")
    myplt = Plots.plot!(Ed_array_S_MPRu,S_synth_MPRu_noisy ./maximum(S_synth_MPRu_noisy),label="S")
    display(myplt)
end
for i=1:length(lambda_values)
    Fr_2D = reshape(F_SD_new*F_sols_synth[:,i],size(f_test))
    myplt = Plots.heatmap(E_coarse,p_coarse,Fr_2D',title="$(i): log10(λ)=$(log10(lambda_values[i]))")
    p0, pN = extrema(p_sols_synth); diffP = abs(pN-p0)
    x0, xN = extrema(x_sols_synth); diffX = abs(xN-x0)
    l_mag = sqrt((p_sols_synth[i]-p0)*(p_sols_synth[i]-p0)/(diffP*diffP)+(x_sols_synth[i]-x0)*(x_sols_synth[i]-x0)/(diffX*diffX))
    myplt1 = Plots.plot([0,p_sols_synth[i]],[0,x_sols_synth[i]],color=:black,title="||L_n||=$(round(l_mag,sigdigits=4))",label="")
    myplt1 = Plots.plot!(p_sols_synth,x_sols_synth,xlims=(minimum(p_sols_synth)-0.05*diffP, maximum(p_sols_synth)+0.05*diffP),ylims=(minimum(x_sols_synth)-0.05*diffX, maximum(x_sols_synth)+0.05*diffX),label="")
    myplt1 = Plots.scatter!(myplt1, [p_sols_synth[1]],[x_sols_synth[1]],label="Start")
    myplt1 = Plots.scatter!(myplt1, [p_sols_synth[end]],[x_sols_synth[end]],label="End")
    #myplt1 = Plots.scatter!(myplt1, [p_sols_synth[ilm]],[x_sols_synth[ilm]],label="gamma_max")
    myplt1 = Plots.scatter!(myplt1, [p_sols_synth[i]],[x_sols_synth[i]],label="$(round(log10(lambda_values[i]),sigdigits=4))")
    myplt_tot = Plots.plot(myplt,myplt1,layout=(1,2),size=(900,400))
    display(myplt_tot)  
end