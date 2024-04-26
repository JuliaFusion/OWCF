###########################################################
# This script is currently under construction.
# When completed, it will be able to solve the S=WF inverse
# problem. That is, given a set of measurements S (either from 
# calcSpec.jl or experimental data), use 0th or 1st order 
# Tikhonov together with the weight functions (computed by
# calcOrbWeights.jl or calc2DWeights.jl or even calc4DWeights.jl)
# to reconstruct the fast-ion distribution F. This is done for 
# many regularization parameter values, and returned as a 
# multi-dimensional array.
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
Ed_array_S_KM14 = Ed_array_S_KM14[gi_KM14]
S_exp_KM14 = spec_KM14[gi_KM14]
err_exp_KM14 = spec_err_KM14[gi_KM14] 
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
Ed_array_S_KM15 = Ed_array_S_KM15[gi_KM15]
S_exp_KM15 = spec_KM15[gi_KM15]
err_exp_KM15 = spec_err_KM15[gi_KM15] 
###########################################################
# MPRu
using HDF5
myfile = h5open("/home/henrikj/Documents/dissemination/papers/99971_2Drec_paper/data/MPRu/99971_MPRu_exp_data.h5","r")
Ed_array_S_MPRu = read(myfile["Xpos_array"])
spec_err_MPRu = read(myfile["erro_array"])
spec_MPRu = read(myfile["spec_array"])
close(myfile)
gi_MPRu = findall(x-> x>=20.0 && x<=60.0,Ed_array_S_MPRu)
Ed_array_S_MPRu = Ed_array_S_MPRu[gi_MPRu]
S_exp_MPRu = spec_MPRu[gi_MPRu]
err_exp_MPRu = spec_err_MPRu[gi_MPRu] 
###########################################################
pltKM14nKM15 = Plots.scatter(Ed_array_S_KM14,S_exp_KM14; label="KM14",markershape=:diamond,markersize=4,markercolor=:green)
pltKM14nKM15 = Plots.yerror!(Ed_array_S_KM14,S_exp_KM14; yerror=err_exp_KM14)
pltKM14nKM15 = Plots.plot!(xlabel="Deposited energy [MeV]",ylabel="Counts",yaxis=:identity)
pltKM14nKM15 = Plots.plot!(title="Experimental data",xtickfontsize=13,xguidefontsize=14)
pltKM14nKM15 = Plots.plot!(ytickfontsize=13,yguidefontsize=14)
pltKM14nKM15 = Plots.plot!(titlefontsize=15)
pltKM14nKM15 = Plots.scatter!(Ed_array_S_KM15,S_exp_KM15; label="KM15",markershape=:diamond,markersize=4,markercolor=:red)
pltKM14nKM15 = Plots.yerror!(Ed_array_S_KM15,S_exp_KM15; yerror=err_exp_KM15)

pltMPRu = Plots.scatter(Ed_array_S_MPRu,S_exp_MPRu; label="MPRu",markershape=:diamond,markersize=4,mc=:purple)
pltMPRu = Plots.yerror!(Ed_array_S_MPRu,S_exp_MPRu; yerror=err_exp_MPRu)
pltMPRu = Plots.plot!(xlabel="Proton impact position [cm]",ylabel="Counts",yaxis=:identity)
pltMPRu = Plots.plot!(title="Experimental data",xtickfontsize=13,xguidefontsize=14)
pltMPRu = Plots.plot!(ytickfontsize=13,yguidefontsize=14)
pltMPRu = Plots.plot!(titlefontsize=15)
###########################################################
# Load KM14 weight functions
using JLD2
using HDF5
using FileIO
filepath_W = "/home/henrikj/Codes/OWCF_results/cycle31/13/weights/velWeights_JET_99971L72_at8,9s_KM14_T-d=n-4He.hdf5"
myfile = h5open(filepath_W,"r")
W_KM14 = read(myfile["W"])
Ed_array_W_KM14 = read(myfile["Ed_array"])
Ed_array_W_KM14 = Ed_array_W_KM14 ./(1000.0) # keV to MeV
E_array = read(myfile["E_array"])
p_array = read(myfile["p_array"])
En_array_W_KM14 = read(myfile["Ed_array_raw"])
close(myfile)
###########################################################
# Define coarse grid vectors. Compute interpolations of 
# weight functions onto coarse grid.
E_coarse = collect(range(minimum(E_array),stop=500.0,length=100))
p_coarse = collect(range(extrema(p_array)...,length=35))
using Interpolations
nodes = (Ed_array_W_KM14, E_array, p_array)

# Create interpolation object
itp = interpolate(nodes, W_KM14, Gridded(Linear()))

W_coarse_KM14 = zeros(length(Ed_array_S_KM14),length(E_coarse),length(p_coarse)) # Pre-allocate 4D array
numObadInds = 0
for Edi=1:length(Ed_array_S_KM14)
    for Ei=1:length(E_coarse)
        for pi=1:length(p_coarse)
            try
                W_coarse_KM14[Edi,Ei,pi] = itp(Ed_array_S_KM14[Edi],E_coarse[Ei],p_coarse[pi]) # Interpolate at 3D query point
            catch
                numObadInds += 1
                debug && println("(Edi: $(Edi), Ei: $(Ei), pi: $(pi)) <--- Interpolation failed for this index") # Print if failed (should not happen)
            end
        end
    end
end
W_2D_KM14_orig = reshape(W_coarse_KM14,(length(Ed_array_S_KM14),length(E_coarse)*length(p_coarse)))
###########################################################
# Load the KM15 weight functions
filepath_W_KM15 = "/home/henrikj/Codes/OWCF_results/cycle31/13/weights/velWeights_JET_99971L72_at8,9s_KM15_T-d=n-4He.hdf5"
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
W_coarse_KM15 = zeros(length(Ed_array_S_KM15),length(E_coarse),length(p_coarse)) # Pre-allocate 4D array
numObadInds = 0
for Edi=1:length(Ed_array_S_KM15)
    for Ei=1:length(E_coarse)
        for pi=1:length(p_coarse)
            try
                W_coarse_KM15[Edi,Ei,pi] = itp(Ed_array_S_KM15[Edi],E_coarse[Ei],p_coarse[pi]) # Interpolate at 3D query point
            catch
                numObadInds += 1
                debug && println("(Edi: $(Edi), Ei: $(Ei), pi: $(pi)) <--- Interpolation failed for this index") # Print if failed (should not happen)
            end
        end
    end
end
W_2D_KM15_orig = reshape(W_coarse_KM15,(length(Ed_array_S_KM15),length(E_coarse)*length(p_coarse)))
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
W_coarse_MPRu = zeros(length(Ed_array_S_MPRu),length(E_coarse),length(p_coarse)) # Pre-allocate 4D array
numObadInds = 0
for Edi=1:length(Ed_array_S_MPRu)
    for Ei=1:length(E_coarse)
        for pi=1:length(p_coarse)
            try
                W_coarse_MPRu[Edi,Ei,pi] = itp(Ed_array_S_MPRu[Edi],E_coarse[Ei],p_coarse[pi]) # Interpolate at 3D query point
            catch
                numObadInds += 1
                debug && println("(Xpos: $(Edi), Ei: $(Ei), pi: $(pi)) <--- Interpolation failed for this index") # Print if failed (should not happen)
            end
        end
    end
end
W_2D_MPRu_orig = reshape(W_coarse_MPRu,(length(Ed_array_S_MPRu),length(E_coarse)*length(p_coarse)))
###########################################################
for iEd in axes(W_coarse_MPRu,1)
    display(Plots.heatmap(E_coarse,p_coarse,W_coarse_MPRu[iEd,:,:]'))
end
###########################################################
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
# Plot test distribution (Maxwellian)
Plots.heatmap(E_coarse,p_coarse,f_test',title="Ground truth")
###########################################################
# Compute and plot synthetic signal
f_1D = reshape(f_test,(length(E_coarse)*length(p_coarse),1))
S_synth_KM14_clean = dropdims(W_2D_KM14_orig * f_1D,dims=2)
Plots.plot(Ed_array_S_KM14,S_synth_KM14_clean)
###########################################################
# Add noise to the signal. Superimpose it on signal plot
k = 0.1
b = 0.05
include("extra/dependencies.jl")
S_synth_KM14_noisy, err_synth_KM14 = add_noise(S_synth_KM14_clean, b; k=k)
Plots.plot!(Ed_array_S_KM14,S_synth_KM14_noisy)
###########################################################
S_synth_KM15_clean = dropdims(W_2D_KM15_orig * f_1D,dims=2)
Plots.plot!(Ed_array_S_KM15, S_synth_KM15_clean)
S_synth_KM15_noisy, err_synth_KM15 = add_noise(S_synth_KM15_clean, b; k=k)
Plots.plot!(Ed_array_S_KM15,S_synth_KM15_noisy)
###########################################################
S_synth_MPRu_clean = dropdims(W_2D_MPRu_orig * f_1D,dims=2)
Plots.plot(Ed_array_S_MPRu,S_synth_MPRu_clean)
S_synth_MPRu_noisy, err_synth_MPRu = add_noise(S_synth_MPRu_clean, b; k=k)
Plots.plot!(Ed_array_S_MPRu,S_synth_MPRu_noisy)
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
# Create SD basis functions

###########################################################
# Rescale weight functions, if necessary
# In final OWCF version, code an intelligent check, to 
# determine if re-scaling is necessary
# It could be a flag like: 
# (caar_W = true) && (typical_N_FI = 1.0e19)
# 'caar_W' stands for check and auto-rescale weights
# A sample of Maxwellians will then be created. They each have 'typical_N_FI'
# as the value of their zeroth moment. Their signals will be computed and 
# compared with the experimental data. If none of the signals are within an 
# order of magnitude from the maximum of any of the experimental signals, 
# the weight functions should be re-scaled as below
W_2D_KM14 = (maximum(S_exp_KM14)/maximum(S_synth_KM14_clean)) .*W_2D_KM14_orig
W_2D_KM15 = (maximum(S_exp_KM15)/maximum(S_synth_KM15_clean)) .*W_2D_KM15_orig
W_2D_MPRu = (maximum(S_exp_MPRu)/maximum(S_synth_MPRu_clean)) .*W_2D_MPRu_orig
###########################################################
# Try using the non-rescaled weight functions on the experimental data
W_2D_KM14 = W_2D_KM14_orig
W_2D_KM15 = W_2D_KM15_orig
W_2D_MPRu = W_2D_MPRu_orig
###########################################################
# Prepare experimental data problem first
noise_floor_factor = 1.0e-4
noise_floor_KM14 = maximum(S_exp_KM14)*noise_floor_factor
noise_floor_KM15 = maximum(S_exp_KM15)*noise_floor_factor
noise_floor_MPRu = maximum(S_exp_MPRu)*noise_floor_factor
replace!(x-> (x<noise_floor_KM14) ? noise_floor_KM14 : x, err_exp_KM14)
replace!(x-> (x<noise_floor_KM15) ? noise_floor_KM15 : x, err_exp_KM15)
replace!(x-> (x<noise_floor_MPRu) ? noise_floor_MPRu : x, err_exp_MPRu)
W_hat = vcat(W_2D_KM14 ./err_exp_KM14, W_2D_KM15 ./err_exp_KM15, W_2D_MPRu ./err_exp_MPRu)
s_hat = vcat(S_exp_KM14 ./err_exp_KM14, S_exp_KM15 ./err_exp_KM15, S_exp_MPRu ./err_exp_MPRu)
###########################################################
# Add 0th order Tikhonov
L0 = zeros(length(f_1D),length(f_1D))
for i=1:length(f_1D)
    for j=1:length(f_1D)
        if i==j
            L0[i,j] = 1.0
        end
    end
end
###########################################################
# Add 1st order Tikhonov
L1 = zeros(length(f_1D),length(f_1D))
for i=1:length(f_1D)
    for j=1:length(f_1D)
        if (i-j)==1
            L1[i,j] = 1.0
        end
        if (i-j)==-1
            L1[i,j] = -1.0
        end
    end
end
L1[1,:] .= 0.0
L1[:,end] .= 0.0
#spy(L1)
###########################################################
#Normalize everything
Whm = maximum(W_hat)
W_hh = W_hat ./Whm
shm = maximum(s_hat)
s_hh = s_hat ./shm
L_orig = L0
f0 = zeros(size(L_orig,1))
###########################################################
# SOLVE THE EXPERIMENTAL DATA PROBLEM
using SCS, Convex
using LinearAlgebra
lambda_values = 10 .^(range(-2,stop=2.0,length=5))
x_sols = zeros(length(lambda_values))
p_sols = zeros(length(lambda_values))
F_sols = zeros(length(f0),length(lambda_values))
for (il,lambda) in enumerate(lambda_values)
    println("lambda: $(lambda)")
    L = lambda*L_orig
    x = Convex.Variable(length(f_1D))

    problem = Convex.minimize(Convex.sumsquares(vcat(W_hh,L) * x - vcat(s_hh,f0)), [x >= 0])

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
for i=1:length(lambda_values)
    Sr = W_hh*F_sols[:,i]
    myplt = Plots.plot(Sr ./maximum(Sr),title="$(i): log10(λ)=$(log10(lambda_values[i]))",label="ŴF*")
    myplt = Plots.plot!(s_hh,label="hat(Ŝ)")
    display(myplt)
end
for i=1:length(lambda_values)
    Sr_KM14 = W_2D_KM14*F_sols[:,i]
    myplt = Plots.plot(Ed_array_S_KM14,Sr_KM14 ./maximum(Sr_KM14),title="$(i): log10(λ)=$(log10(lambda_values[i]))",label="WF*")
    myplt = Plots.plot!(Ed_array_S_KM14,S_exp_KM14 ./maximum(S_exp_KM14),label="S")
    display(myplt)
end
for i=1:length(lambda_values)
    Sr_KM15 = W_2D_KM15*F_sols[:,i]
    myplt = Plots.plot(Ed_array_S_KM15,Sr_KM15 ./maximum(Sr_KM15),title="$(i): log10(λ)=$(log10(lambda_values[i]))",label="WF*")
    myplt = Plots.plot!(Ed_array_S_KM15,S_exp_KM15 ./maximum(S_exp_KM15),label="S")
    display(myplt)
end
for i=1:length(lambda_values)
    Sr_MPRu = W_2D_MPRu*F_sols[:,i]
    myplt = Plots.plot(Ed_array_S_MPRu,Sr_MPRu ./maximum(Sr_MPRu),title="$(i): log10(λ)=$(log10(lambda_values[i]))",label="WF*")
    myplt = Plots.plot!(Ed_array_S_MPRu,S_exp_MPRu ./maximum(S_exp_MPRu),label="S")
    display(myplt)
end
for i=1:length(lambda_values)
    Fr_2D = reshape(F_sols[:,i],size(f_test))
    myplt = Plots.heatmap(E_coarse,p_coarse,Fr_2D',title="$(i): log10(λ)=$(log10(lambda_values[i]))")
    p0, pN = extrema(p_sols); diffP = abs(pN-p0)
    x0, xN = extrema(x_sols); diffX = abs(xN-x0)
    l_mag = sqrt((p_sols[i]-p0)*(p_sols[i]-p0)/(diffP*diffP)+(x_sols[i]-x0)*(x_sols[i]-x0)/(diffX*diffX))
    myplt1 = Plots.plot([0,p_sols[i]],[0,x_sols[i]],color=:black,title="||L_n||=$(round(l_mag,sigdigits=4))",label="")
    myplt1 = Plots.plot!(p_sols,x_sols,xlims=(minimum(p_sols)-0.05*diffP, maximum(p_sols)+0.05*diffP),ylims=(minimum(x_sols)-0.05*diffX, maximum(x_sols)+0.05*diffX),label="")
    myplt1 = Plots.scatter!(myplt1, [p_sols[1]],[x_sols[1]],label="Start")
    myplt1 = Plots.scatter!(myplt1, [p_sols[end]],[x_sols[end]],label="End")
    myplt1 = Plots.scatter!(myplt1, [p_sols[ilm]],[x_sols[ilm]],label="gamma_max")
    myplt1 = Plots.scatter!(myplt1, [p_sols[i]],[x_sols[i]],label="$(round(log10(lambda_values[i]),sigdigits=4))")
    myplt_tot = Plots.plot(myplt,myplt1,layout=(1,2),size=(900,400))
    display(myplt_tot)  
end
###########################################################
# Redo the reconstruction. But for the synthetic case
# For comparison and sanity check
replace!(x-> (x==0.0) ? 1.0 : x, err_synth_KM14)
replace!(x-> (x==0.0) ? 1.0 : x, err_synth_KM15)
replace!(x-> (x==0.0) ? 1.0 : x, err_synth_MPRu)
# W_2D_KM14, W_2D_KM15 and W_2D_MPRu are f*cked here 
# Because of normalize to fit experimental data
W_hat = vcat(W_2D_KM14_orig ./err_synth_KM14, W_2D_KM15_orig ./err_synth_KM15, W_2D_MPRu_orig ./err_synth_MPRu)
s_hat = vcat(S_synth_KM14_noisy ./err_synth_KM14, S_synth_KM15_noisy ./err_synth_KM15, S_synth_MPRu_noisy ./err_synth_MPRu)
Whm = maximum(W_hat)
W_hh = W_hat ./Whm
shm = maximum(s_hat)
s_hh = s_hat ./shm
f0 = zeros(size(L,1))
x_sols_synth = zeros(length(lambda_values))
p_sols_synth = zeros(length(lambda_values))
F_sols_synth = zeros(length(f0),length(lambda_values))
for (il,lambda) in enumerate(lambda_values)
    println("lambda: $(lambda)")
    L = lambda*L0
    x = Convex.Variable(length(f_1D))

    problem = Convex.minimize(Convex.sumsquares(vcat(W_hh,L) * x - vcat(s_hh,f0)), [x >= 0])

    solve!(problem, SCS.Optimizer)

    x_sols_synth[il] = norm(dropdims(x.value,dims=2))
    p_sols_synth[il] = norm(W_hh*dropdims(x.value,dims=2) - s_hh)
    F_sols_synth[:,il] = (shm/Whm) .*dropdims(x.value,dims=2)
end
###########################################################
# Check the synthetic reconstruction
for i=1:length(lambda_values)
    Sr = W_hh*F_sols_synth[:,i]
    myplt = Plots.plot(Sr ./maximum(Sr),title="$(i): log10(λ)=$(log10(lambda_values[i]))",label="ŴF*")
    myplt = Plots.plot!(s_hh ./maximum(s_hh),label="Ŝ")
    display(myplt)
end
for i=1:length(lambda_values)
    Sr_KM14 = W_2D_KM14*F_sols_synth[:,i]
    myplt = Plots.plot(Ed_array_S_KM14,Sr_KM14 ./maximum(Sr_KM14),title="$(i): log10(λ)=$(log10(lambda_values[i]))",label="WF*")
    myplt = Plots.plot!(Ed_array_S_KM14,S_synth_KM14_noisy ./maximum(S_synth_KM14_noisy),label="S")
    display(myplt)
end
for i=1:length(lambda_values)
    Sr_KM15 = W_2D_KM15*F_sols_synth[:,i]
    myplt = Plots.plot(Ed_array_S_KM15,Sr_KM15 ./maximum(Sr_KM15),title="$(i): log10(λ)=$(log10(lambda_values[i]))",label="WF*")
    myplt = Plots.plot!(Ed_array_S_KM15,S_synth_KM15_noisy ./maximum(S_synth_KM15_noisy),label="S")
    display(myplt)
end
for i=1:length(lambda_values)
    Sr_MPRu = W_2D_MPRu*F_sols_synth[:,i]
    myplt = Plots.plot(Ed_array_S_MPRu,Sr_MPRu ./maximum(Sr_MPRu),title="$(i): log10(λ)=$(log10(lambda_values[i]))",label="WF*")
    myplt = Plots.plot!(Ed_array_S_MPRu,S_synth_MPRu_noisy ./maximum(S_synth_MPRu_noisy),label="S")
    display(myplt)
end
for i=1:length(lambda_values)
    Fr_2D = reshape(F_sols_synth[:,i],size(f_test))
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
#########################################################
# Define function to compute slowing-down basis function