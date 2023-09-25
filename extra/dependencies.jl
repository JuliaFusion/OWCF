######################################### dependencies.jl ##################################################
# This script defines a large collection of functions that is used by various scripts and apps of the OWCF.
#
# Prior to loading this collection of functions (for example via 'include(extra/dependencies.jl)' when standing
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
# This is performed in e.g. every OWCF start template file. It is also possible to skip the Distributed.jl
# package. This is then done as 
#
# folderpath_OWCF = "/path/to/the/OWCF/folder/"
# using Pkg
# cd(folderpath_OWCF)#
# Pkg.activate(".")
#
# Written by H. Järleblad. Last updated 2022-09-22.
###################################################################################################

println("Loading the Julia packages for the OWCF dependencies... ")
using Distributed
using EFIT
using Equilibrium
using GuidingCenterOrbits
using ForwardDiff
using LinearAlgebra
using OrbitTomography
using Optim
using JLD2
using HDF5
using FileIO
using Interpolations
using NearestNeighbors
using Statistics
using SharedArrays
using ProgressMeter
using Distributions
using VoronoiDelaunay
import OrbitTomography: orbit_grid

"""
    add_noise(S,0.05)
    add_noise(-||-, k=0.15)

A function that adds noise to a signal in a consistent way. The function adds both background noise as well as
noise to the signal itself. The level of the background noise and the signal noise is determined by the b
and k variables respectively. The levels are input as a fraction of the mean of the original signal S strength.
"""
function add_noise(s::AbstractVector, b::Union{Float64, AbstractVector}; k::Float64=0.1)
    sn = max.(s,0.0) .+ k.*(mean(sqrt.(abs.(s)))).*rand.(Normal.(0.0, max.(sqrt.(max.(s,0)), sqrt.(b))))
    err = k.*mean(sqrt.(abs.(s))).*max.(sqrt.(max.(s,0)), sqrt.(b))
    return sn, err
end

"""
    estimate_noise(S_tot)

A function that estimates the noise of a signal S. With the 'bad_inds' keyword argument,
estimate the noise of the S_tot-elements with the specific array indices in 'bad_inds'. 
Assume that the true signal 'S_tot' is a smooth function without the noise.
"""
function estimate_noise(S_tot::AbstractVector; bad_inds::Union{Nothing,AbstractVector}=nothing)
    noise_tot = zeros(length(S_tot))
    if bad_inds==nothing
        bad_inds = 1:length(S_tot)
    end
    for bad_ind in bad_inds # A simple noise estimation algorithm
        if bad_ind==1 # If the bad index is the very first element, we need to extrapolate twice
            S_bad_one_ahead_mean = (S_tot[bad_ind]+S_tot[bad_ind+1]+S_tot[bad_ind+2])/3
            S_extrap = 2*S_tot[bad_ind] - S_bad_one_ahead_mean
            S_worse_mean = (S_extrap+S_tot[bad_ind]+S_tot[bad_ind+1])/3
            S_extrap_extrap = 2*S_extrap - S_worse_mean 
            S_bad_mean = (S_extrap_extrap+S_extrap+S_tot[bad_ind]+S_tot[bad_ind+1]+S_tot[bad_ind+2])/5 # The average of all five points close to the point of interest
        elseif bad_ind==2 # If the bad index is the second element, we need to extrapolate once
            S_worse_mean = (S_tot[bad_ind-1]+S_tot[bad_ind]+S_tot[bad_ind+1])/3
            S_extrap = 2*S_tot[bad_ind-1] - S_worse_mean
            S_bad_mean = (S_extrap+S_tot[bad_ind-1]+S_tot[bad_ind]+S_tot[bad_ind+1]+S_tot[bad_ind+2])/5 # The average of all five points close to the point of interest
        elseif bad_ind==(length(S_tot)-1) # If the bad index is the next to last element, we need to extrapolate once
            S_worse_mean = (S_tot[bad_ind-1]+S_tot[bad_ind]+S_tot[bad_ind+1])/3
            S_extrap = 2*S_tot[bad_ind+1]-S_worse_mean
            S_bad_mean = (S_tot[bad_ind-2]+S_tot[bad_ind-1]+S_tot[bad_ind]+S_tot[bad_ind+1]+S_extrap)/5 # The average of all five points close to the point of interest
        elseif bad_ind==length(S_tot) # If the bad index is the very last element, we need to extrapolate twice
            S_bad_one_behind_mean = (S_tot[bad_ind-2]+S_tot[bad_ind-1]+S_tot[bad_ind])/3
            S_extrap = 2*S_tot[bad_ind]-S_bad_one_behind_mean
            S_worse_mean = (S_tot[bad_ind-1]+S_tot[bad_ind]+S_extrap)/3
            S_extrap_extrap = 2*S_extrap - S_worse_mean
            S_bad_mean = (S_tot[bad_ind-2]+S_tot[bad_ind-1]+S_tot[bad_ind]+S_extrap+S_extrap_extrap)/5 # The average of all five points close to the point of interest
        else
            S_bad_mean = (S_tot[bad_ind-2]+S_tot[bad_ind-1]+S_tot[bad_ind]+S_tot[bad_ind+1]+S_tot[bad_ind+2])/5 # The average of all five points close to the point of interest
        end
    
        noise_tot[bad_ind] = abs(S_tot[bad_ind] - S_bad_mean) # Noise (or error)
    end

    return noise_tot
end

"""
    apply_instrumental_response(S, Ed_array, instrumental_response_input, instrumental_response_output, instrumental_response_matrix)
    apply_instrumental_response(-||-; lo=nothing, hi=nothing, instrumental_response=true)

Apply instrumental response to the signal S via the matrix instrumental_response matrix. The input and output axes of the instrumental response
matrix are specified as vector instrumental_response_input and instrumental_response_output. The index of the element closest to minimum(Ed_array)
in instrumental_response_input can be specified via the lo keyword argument. The index of the element closest to maximum(Ed_array) in instrumental_response_input
can be specified via the hi keyword argument.
"""
function apply_instrumental_response(S::AbstractVector, Ed_array::AbstractVector, instrumental_response_input::AbstractVector, instrumental_response_output::AbstractVector, instrumental_response_matrix::AbstractMatrix; lo::Union{Int64,Nothing}=nothing, hi::Union{Int64,Nothing}=nothing)
    if isnothing(lo) # If no lower bound for the instrumental response input has been provided...
        lo = findfirst(x-> x>minimum(Ed_array),instrumental_response_input) # Find it
    end
    if isnothing(hi) # Similar as lo, but hi
        hi = findlast(x-> x<maximum(Ed_array),instrumental_response_input)
    end
    
    if isnothing(lo) || isnothing(hi) # If extrema(Ed_array) are (partially) outside of extrema(instrumental_response_input)...
        @warn "Instrumental response matrix input completely outside weight function measurement bin range. No diagnostic response will be applied."
        return S
    else
        if lo==1
            @warn "Lower bound of instrumental response matrix input might not be low enough to cover weight function measurement bin range. Diagnostic response representation might be inaccurate."
        end
        if hi==length(instrumental_response_input)
            @warn "Upper bound of instrumental response matrix input might not be high enough to cover weight function measurement bin range. Diagnostic response representation might be inaccurate."
        end
        instrumental_response_matrix = (instrumental_response_matrix[lo:hi,:])' # Take the transpose, to go from (input, output) shape to (output, input) shape
        itp = LinearInterpolation(Ed_array,S) # Create interpolation object for S
        S_itp = itp.(instrumental_response_input[lo:hi]) # Find interpolation values for S and instrumental_response_input points
        S_out = instrumental_response_matrix * S_itp # The instrumental response is applied
        return S_out
    end
end

"""
    apply_instrumental_response(Q, Ed_array, instrumental_response_input, instrumental_response_output, instrumental_response_matrix)
    apply_instrumental_response(-||-; lo=nothing, hi=nothing, instrumental_response=true)

Apply the instrumental response given by instrumental_response_matrix to each column of Q.
"""
function apply_instrumental_response(Q::Array{Float64,2}, Ed_array::AbstractVector, instrumental_response_input::AbstractVector, instrumental_response_output::AbstractVector, instrumental_response_matrix::AbstractMatrix; lo::Union{Int64,Nothing}=nothing, hi::Union{Int64,Nothing}=nothing)
    if isnothing(lo)
        lo = findfirst(x-> x>minimum(Ed_array),instrumental_response_input)
    end
    if isnothing(hi)
        hi = findlast(x-> x<maximum(Ed_array),instrumental_response_input)
    end
    
    if isnothing(lo) || isnothing(hi)
        @warn "Instrumental response matrix input completely outside weight function measurement bin range. No diagnostic response will be applied."
        return Q
    else
        Q_out = zeros(length(instrumental_response_output), size(Q,2))
        if lo==1
            @warn "Lower bound of instrumental response matrix input might not be low enough to cover weight function measurement bin range. Diagnostic response representation might be inaccurate."
        end
        if hi==length(instrumental_response_input)
            @warn "Upper bound of instrumental response matrix input might not be high enough to cover weight function measurement bin range. Diagnostic response representation might be inaccurate."
        end
        for i=1:size(Q,2) # Apply instrumental response to each column
            Q_out[:,i] = apply_instrumental_response(Q[:,i], Ed_array, instrumental_response_input, instrumental_response_output, instrumental_response_matrix; lo=lo, hi=hi)
        end
        return Q_out
    end
end
"""
    apply_instrumental_response(Q, Ed_array, instrumental_response_input, instrumental_response_output, instrumental_response_matrix)
    apply_instrumental_response(-||-; lo=nothing, hi=nothing, instrumental_response=true)

Apply the instrumental response given by instrumental_response_matrix to each matrix of Q.
"""
function apply_instrumental_response(Q::Array{Float64,3}, Ed_array::AbstractVector, instrumental_response_input::AbstractVector, instrumental_response_output::AbstractVector, instrumental_response_matrix::AbstractMatrix; lo::Union{Int64,Nothing}=nothing, hi::Union{Int64,Nothing}=nothing)
    if isnothing(lo)
        lo = findfirst(x-> x>minimum(Ed_array),instrumental_response_input)
    end
    if isnothing(hi)
        hi = findlast(x-> x<maximum(Ed_array),instrumental_response_input)
    end
    
    if isnothing(lo) || isnothing(hi)
        @warn "Instrumental response matrix input completely outside weight function measurement bin range. No diagnostic response will be applied."
        return Q
    else
        Q_out = zeros(length(instrumental_response_output), size(Q,2), size(Q,3))
        if lo==1
            @warn "Lower bound of instrumental response matrix input might not be low enough to cover weight function measurement bin range. Diagnostic response representation might be inaccurate."
        end
        if hi==length(instrumental_response_input)
            @warn "Upper bound of instrumental response matrix input might not be high enough to cover weight function measurement bin range. Diagnostic response representation might be inaccurate."
        end
        for i=1:size(Q,3) # Apply instrumental response to each matrix
            Q_out[:,:,i] = apply_instrumental_response(Q[:,:,i], Ed_array, instrumental_response_input, instrumental_response_output, instrumental_response_matrix; lo=lo, hi=hi)
        end
        return Q_out
    end
end
"""
    apply_instrumental_response(Q, Ed_array, instrumental_response_input, instrumental_response_output, instrumental_response_matrix)
    apply_instrumental_response(-||-; lo=nothing, hi=nothing, instrumental_response=true)

Apply the instrumental response given by instrumental_response_matrix to each 3D array of Q.
"""
function apply_instrumental_response(Q::Array{Float64,4}, Ed_array::AbstractVector, instrumental_response_input::AbstractVector, instrumental_response_output::AbstractVector, instrumental_response_matrix::AbstractMatrix; lo::Union{Int64,Nothing}=nothing, hi::Union{Int64,Nothing}=nothing)
    if isnothing(lo)
        lo = findfirst(x-> x>minimum(Ed_array),instrumental_response_input)
    end
    if isnothing(hi)
        hi = findlast(x-> x<maximum(Ed_array),instrumental_response_input)
    end
    
    if isnothing(lo) || isnothing(hi)
        @warn "Instrumental response matrix input completely outside weight function measurement bin range. No diagnostic response will be applied."
        return Q
    else
        Q_out = zeros(length(instrumental_response_output), size(Q,2), size(Q,3), size(Q,4))
        if lo==1
            @warn "Lower bound of instrumental response matrix input might not be low enough to cover weight function measurement bin range. Diagnostic response representation might be inaccurate."
        end
        if hi==length(instrumental_response_input)
            @warn "Upper bound of instrumental response matrix input might not be high enough to cover weight function measurement bin range. Diagnostic response representation might be inaccurate."
        end
        for i=1:size(Q,4) # Apply instrumental response to each 3D array
            Q_out[:,:,:,i] = apply_instrumental_response(Q[:,:,:,i], Ed_array, instrumental_response_input, instrumental_response_output, instrumental_response_matrix; lo=lo, hi=hi)
        end
        return Q_out
    end
end

"""
    apply_instrumental_response(Q, Ed_array, instrumental_response_input, instrumental_response_output, instrumental_response_matrix)
    apply_instrumental_response(-||-; lo=nothing, hi=nothing)

Apply the instrumental response given by instrumental_response_matrix to each 4D array of Q.
"""
function apply_instrumental_response(Q::Array{Float64,5}, Ed_array::AbstractVector, instrumental_response_input::AbstractVector, instrumental_response_output::AbstractVector, instrumental_response_matrix::AbstractMatrix; lo::Union{Int64,Nothing}=nothing, hi::Union{Int64,Nothing}=nothing)
    if isnothing(lo)
        lo = findfirst(x-> x>minimum(Ed_array),instrumental_response_input)
    end
    if isnothing(hi)
        hi = findlast(x-> x<maximum(Ed_array),instrumental_response_input)
    end
    
    if isnothing(lo) || isnothing(hi)
        @warn "Instrumental response matrix input completely outside weight function measurement bin range. No diagnostic response will be applied."
        return Q
    else
        Q_out = zeros(length(instrumental_response_output), size(Q,2), size(Q,3), size(Q,4), size(Q,5))
        if lo==1
            @warn "Lower bound of instrumental response matrix input might not be low enough to cover weight function measurement bin range. Diagnostic response representation might be inaccurate."
        end
        if hi==length(instrumental_response_input)
            @warn "Upper bound of instrumental response matrix input might not be high enough to cover weight function measurement bin range. Diagnostic response representation might be inaccurate."
        end
        for i=1:size(Q,5) # Apply instrumental response to each 4D array
            Q_out[:,:,:,:,i] = apply_instrumental_response(Q[:,:,:,:,i], Ed_array, instrumental_response_input, instrumental_response_output, instrumental_response_matrix; lo=lo, hi=hi)
        end
        return Q_out
    end
end

"""
    orbit_grid(M,false,eo,po,ro)
    orbit_grid(-||-, q=2, amu=OrbitTomography.p_amu, kwargs... )

This function is just like OrbitTomography.orbit_grid, but allows execution without progress bar.
Good for HPC batch job submission, to not overrun log files.
"""
function orbit_grid(M::AbstractEquilibrium, visualizeProgress::Bool, eo::AbstractVector, po::AbstractVector, ro::AbstractVector;
                    q::Int64 = 1, amu = OrbitTomography.H2_amu, kwargs...)

    nenergy = length(eo)
    npitch = length(po)
    nr = length(ro)

    orbit_index = zeros(Int,nenergy,npitch,nr)
    class = fill(:incomplete,(nenergy,npitch,nr))
    tau_t = zeros(Float64,nenergy,npitch,nr)
    tau_p = zeros(Float64,nenergy,npitch,nr)

    norbs = nenergy*npitch*nr
    subs = CartesianIndices((nenergy,npitch,nr))

    if visualizeProgress
        p = Progress(norbs)
        channel = RemoteChannel(()->Channel{Bool}(norbs), 1)
        orbs = fetch(@sync begin
            @async while take!(channel)
                ProgressMeter.next!(p)
            end
            @async begin
                orbs = @distributed (vcat) for i=1:norbs
                    ie,ip,ir = Tuple(subs[i])
                    c = GuidingCenterOrbits.EPRCoordinate(M,eo[ie],po[ip],ro[ir],q=q,amu=amu)
                    try
                        o = GuidingCenterOrbits.get_orbit(M, c; kwargs...)
                    catch
                        o = GuidingCenterOrbits.Orbit(EPRCoordinate(;q=q,amu=amu),:incomplete)
                    end

                    if o.class in (:incomplete,:invalid,:lost)
                        o = Orbit(o.coordinate,:incomplete)
                    end
                    put!(channel, true)
                    o
                end
                put!(channel, false)
                orbs
            end
        end)
    else
        orbs = @distributed (vcat) for i=1:norbs
            ie,ip,ir = Tuple(subs[i])
            c = GuidingCenterOrbits.EPRCoordinate(M,eo[ie],po[ip],ro[ir],q=q,amu=amu)
            try
                o = GuidingCenterOrbits.get_orbit(M, c; kwargs...)
            catch
                o = GuidingCenterOrbits.Orbit(EPRCoordinate(;q=q,amu=amu),:incomplete)
            end

            if o.class in (:incomplete,:invalid,:lost)
                o = Orbit(o.coordinate,:incomplete)
            end
            o
        end
    end


    for i=1:norbs
        class[subs[i]] = orbs[i].class
        tau_p[subs[i]] = orbs[i].tau_p
        tau_t[subs[i]] = orbs[i].tau_t
    end

    grid_index = filter(i -> orbs[i].class != :incomplete, 1:norbs)
    orbs = filter(x -> x.class != :incomplete, orbs)
    norbs = length(orbs)
    orbit_index[grid_index] = 1:norbs

    return orbs, OrbitTomography.OrbitGrid(eo,po,ro,fill(1,norbs),orbit_index,class,tau_p,tau_t)
end

"""
    mu_func()
    mu_func()

Bla bla bla
"""
function mu_func(energy::T, B::T, Pϕ::T, Ψ::T, RBϕ::T;m::Float64=GuidingCenterOrbits.H2_amu*GuidingCenterOrbits.mass_u, q::Float64=1*GuidingCenterOrbits.e0) where {T}
    res = energy/B - (B/(2*m)) * ((Pϕ-q*Ψ)/RBϕ)^2
    return (res > zero(energy)) ? res : zero(energy)
end

"""
    get_orbel_volume(myOrbitGrid, os_equidistant)

Get the orbit element volume of this specific orbit-grid.
If os_equidistant=false, then return a 3D array with all the orbit volumes.
Otherwise, just a scalar.
"""
function get_orbel_volume(og::OrbitGrid, os_equidistant::Bool)

    if os_equidistant
        dRm = abs(og.r[2]-og.r[1]) # The difference between the first and second element should be representative for the whole grid, due to equidistancy
        dE = abs(og.energy[2]-og.energy[1])
        dpm = abs(og.pitch[2]-og.pitch[1])
        dO = dE*dpm*dRm
        return dO
    else
        dO = zeros(length(og.energy), length(og.pitch), length(og.r)) # If not equidistant, create a 3D array with zeros
        for Ei=1:length(og.energy) # For all energies...
            if Ei==length(og.energy) # If at the edge of the orbit grid...
                # Assume edge orbit-element volume to be same as next-to-edge
                dEi = abs(og.energy[end]-og.energy[end-1])
            else
                dEi = abs(og.energy[Ei+1]-og.energy[Ei])
            end
            for pmi=1:length(og.pitch)
                if pmi==length(og.pitch)
                    dpmi = abs(og.pitch[end]-og.pitch[end-1])
                else
                    dpmi = abs(og.pitch[pmi+1]-og.pitch[pmi])
                end
                for Rmi=1:length(og.r)
                    if Rmi==length(og.r)
                        dRmi = abs(og.r[end]-og.r[end-1])
                    else
                        dRmi = abs(og.r[Rmi+1]-og.r[Rmi])
                    end

                    dO[Ei, pmi, Rmi] = dEi*dpmi*dRmi # Replace the zero-element in the 3D array with the resulting orbit-element volume
                end
            end
        end

        return dO
    end
end

"""
    closest_index(myArray, val)

Finds (the index of) the closest value to val in myArray. Just a function found on the internet as open source. It works.
"""
function closest_index(x::AbstractArray, val::Number)
    ibest = first(eachindex(x))
    dxbest = abs(x[ibest]-val)
    for I in eachindex(x)
        dx = abs(x[I]-val)
        if dx < dxbest
            dxbest = dx
            ibest = I
        end
    end
    ibest
end

"""
    interpFps(f,E,p,R,z,Eq,pq,Rq,zq)

The 4D particle-space fast-ion distribution f, with the corresponding
grid vectors in energy E, pitch p, radius R and vertical position z, are supplied as input. Evaluate the
fast-ion distribution at the query points (Eq, pq, Rq and zq) using linear interpolation in 4D.
Return the resulting 4D fast-ion distribution.
"""
function interpFps(f::AbstractArray, E::AbstractArray, p::AbstractArray, R::AbstractArray, z::AbstractArray, Eq::AbstractArray, pq::AbstractArray, Rq::AbstractArray, zq::AbstractArray; debug::Bool=false)

    # Assume equidistant elements
    nodes = (E, p, R, z)

    # Create interpolation object
    itp = interpolate(nodes, f, Gridded(Linear()))

    fq = zeros(length(Eq),length(pq),length(Rq),length(zq)) # Pre-allocate 4D array
    numObadInds = 0
    for Ei=1:length(Eq)
        for pi=1:length(pq)
            for Ri=1:length(Rq)
                for zi=1:length(zq)
                    try
                        fq[Ei,pi,Ri,zi] = itp(Eq[Ei],pq[pi],Rq[Ri],zq[zi]) # Interpolate at 4D query point
                    catch
                        numObadInds += 1
                        debug && println("(Ei: $(Ei), pi: $(pi), Ri: $(Ri), zi: $(zi)) <--- Interpolation failed for this index") # Print if failed (should not happen)
                    end
                end
            end
        end
    end

    if numObadInds > 0
        perc_frac = (numObadInds / length(fq))*100
        @warn "The interpolation algorithm failed to interpolate $(numObadInds) (E,p,R,z) points ($(round(perc_frac,sigdigits=3)) %)"
    end

    return fq
end

"""
    do4Dto2D(og, W)
    do4Dto2D(-||-, returnComplement=true)

This function converts orbit weights of the form (channel,E,pm,Rm) to the form
(channel,orbits). 'orbits' consists of the valid orbits for the (E,pm,Rm)-grid.
Usually, in orbit space, for a given orbit grid only about 67% of the (E,pm,Rm)-points
corresponds to valid orbits. So extracting the actual valid points as a 1D vector makes sense.

If returnComplement, then a 1D array with all the weight values
that should be zero is returned (C), with the corresponding coordinates (Ccoords) and indices (Cindices) as arrays of triplets.
The sum of the C array should be zero. If it is not, then invalid orbits have been given a non-zero weight. Good for debugging.
"""
function do4Dto2D(og::OrbitGrid, W4D::AbstractArray; returnComplement::Bool=false)

    W = zeros(size(W4D,1),length(og.counts)) # Pre-allocate the 2D (channel,orbits) matrix
    if returnComplement
        C = zeros(length(W4D)-length(og.counts)) #Pre-allocate
        Ccoords = Array{Tuple{Float64,Float64,Float64},1}(undef,length(W4D)-length(og.counts)) #Pre-allocate
        Cindices = Array{Tuple{Int64,Int64,Int64},1}(undef,length(W4D)-length(og.counts)) #Pre-allocate
        ii = 1
    end
    for c=1:size(W4D,1)
        for Ei=1:size(W4D,2)
            for pmi=1:size(W4D,3)
                for Rmi=1:size(W4D,4)
                    o_index = og.orbit_index[Ei,pmi,Rmi] # By the method of L. Stagner, every valid orbit corresponds to a non-zero integer: their index. The invalid and lost orbits are simply represented by a zero as their index.
                    if o_index==0 && returnComplement
                        C[ii] = W4D[c,Ei,pmi,Rmi]
                        Ccoords[ii] = (og.energy[Ei],og.pitch[pmi],og.r[Rmi])
                        Cindices[ii] = (Ei,pmi,Rmi)
                        ii = ii + 1
                    elseif o_index==0
                    else
                        W[c,o_index] = W4D[c,Ei,pmi,Rmi] # Put the orbit weight at the correct orbit position (column) in the 2D weight matrix
                    end
                end
            end
        end
    end
    if returnComplement
        return W,C,Ccoords,Cindices
    else
        return W
    end
end

"""
    flip_along_pm(W)
    flip_along_pm(W, with_channels=true)

Flip the elements of an orbit-space function along the pm axis. If with_channels=true,
then a 4D orbit-space function on the form (channel,E,pm,Rm) or (channel,Rm,pm,E) is assumed.
Otherwise a form of (E,pm,Rm) is assumed.
"""
function flip_along_pm(W::AbstractArray; with_channels::Bool=false)

    Wres = zeros(size(W))
    if with_channels
        for pmi=1:size(W,3)
            Wres[:,:,(end+1)-pmi,:] .= W[:,:,pmi,:]
        end
    else
        for pmi=1:size(W,2)
            Wres[:,(end+1)-pmi,:] .= W[:,pmi,:]
        end
    end

    return Wres
end

"""
    flip_along_p(F_EpRz)

This function flips the elements of a 4D fast-ion distribution along the pitch axis.
Assume the fast-ion distribution to be of the form (E,p,R,z)
"""
function flip_along_p(F_EpRz::AbstractArray)

    F_EpRz_res = zeros(size(F_EpRz))
    for pi=1:size(F_EpRz,2)
        F_EpRz_res[:,end+1-pi,:,:] .= F_EpRz[:,pi,:,:]
    end
    return F_EpRz_res
end

"""
    h5To4D(filepath_distr)
    h5To4D(filepath_distr, backwards=false, verbose=true)

Load and read an h5 file containing the data necessary to construct a 4D fast-ion distribution, with dimensions (E,p,R,z).
Automatically assume that the data will be loaded backwards, because that seems to be the case when exporting 4D
data from Python. Returns f, E, p, R, z
"""
function h5to4D(filepath_distr::AbstractString; backwards::Bool = true, verbose::Bool = false)

    if verbose
        println("Loading the 4D distribution from .h5 file... ")
    end
    distr_file = h5open(filepath_distr) # With first NES orbit weights
    fbackwards = read(distr_file,"f")
    energy = read(distr_file, "energy")
    pitch = read(distr_file, "pitch")
    R = read(distr_file, "R")
    if haskey(distr_file,"z")
        z = read(distr_file,"z")
    else
        z = read(distr_file, "Z")
    end

    if (size(fbackwards,1)==length(energy)) && (size(fbackwards,2)==length(pitch)) && (size(fbackwards,3)==length(R)) && (size(fbackwards,4)==length(z))
        backwards = false
    end

    if backwards
        verbose && println("Fast-ion distribution data is permutated. Attemping to fix... ")
        f = zeros(size(fbackwards,4),size(fbackwards,3),size(fbackwards,2),size(fbackwards,1))
        for i=1:size(f,1)
            for j=1:size(f,2)
                for k=1:size(f,3)
                    for l=1:size(f,4)
                        f[i,j,k,l] = fbackwards[l,k,j,i]
                    end
                end
            end
        end
    else
        f = fbackwards
    end

    return f, energy, pitch, R, z
end

"""
    jld2toh5(filepath_jld2)
    jld2toh5(filepath_jld2; verbose=false)

Convert a .jld2 file to a .hdf5 file.
"""
function jld2tohdf5(filepath_jld2::String; verbose::Bool=false)
    filepath_hdf5 = reduce(*, split(filepath_jld2,".")[1:end-1])*".hdf5"
    myfile_jld2 = jldopen(filepath_jld2,false,false,false,IOStream)
    myfile_hdf5 = h5open(filepath_hdf5,"w")

    for key in keys(myfile_jld2)
        verbose && println("Writing data: "*key)
        data = myfile_jld2[key]
        write(myfile_hdf5,key,data)
    end
    close(myfile_jld2)
    close(myfile_hdf5)
    verbose && println("All .jld2 file data written to "*filepath_hdf5)
    return filepath_hdf5
end

"""
    JLD2to4D(filepath_distr)
    JLD2to4D(filepath_distr; verbose=false) # default

Load the fast-ion distribution from a .jld2 file. Assume certain keys to access data. Print information if 'verbose'.
"""
function JLD2to4D(filepath_distr::String; verbose::Bool=false)
    myfile = jldopen(filepath_distr,false,false,false,IOStream)
    if haskey(myfile,"F_ps")
        verbose && println("Found F_ps key in .jld2 file.")
        F_EpRz = myfile["F_ps"]
    elseif haskey(myfile,"f")
        verbose && println("Found f key in .jld2 file.")
        F_EpRz = myfile["f"]
    elseif haskey(myfile,"F_EpRz")
        verbose && println("Found F_EpRz key in .jld2 file.")
        F_EpRz = myfile["F_EpRz"]
    else
        error("Fast-ion distribution .jld2 file did not have any expected keys for the distribution (F_ps, f or F_EpRz).")
    end
    E_array = myfile["energy"]
    p_array = myfile["pitch"]
    R_array = myfile["R"]
    if haskey(myfile,"z")
        verbose && println("Found z key in .jld2 file.")
        z_array = myfile["z"]
    else
        verbose && println("Found Z key in .jld2 file.")
        z_array = myfile["Z"]
    end
    close(myfile)

    return F_EpRz, E_array, p_array, R_array, z_array
end

"""
    h5toSampleReady(filepath_distr, interp)
    h5toSampleReady(-||-,nE_ps=201, np_ps=56, nR_ps=61, nz_ps=63, backwards=false, verbose=true)

Load the fast-ion distribution from a .h5 file and forward to function toSampleReady().

--- Input:
filepath_distr - The fast-ion distribution in .h5 file format - String
interp - If true, then fast-ion distribution will be interpolated onto grid with same end-points but dimension (nE_ps, np_ps, nR_ps, nz_ps) - Bool
nE_ps - Energy grid dimension to be interpolated onto, if interp - Int64
np_ps - Pitch grid dimension to be interpolated onto, if interp - Int64
nR_ps - R grid dimension to be interpolated onto, if interp - Int64
nz_ps - z grid dimension to be interpolated onto, if interp - Int64
backwards - If saved with Python, the .h5 file will have (z,R,p,E) orientation instead of (E,p,R,z). Set to false, if that is not that case - Bool
verbose - If set to true, you will get a talkative function indeed - Bool
"""
function h5toSampleReady(filepath_distr::AbstractString; interp::Bool=false, nE_ps::Int64=100, np_ps::Int64=51, nR_ps::Int64=60, nz_ps::Int64=61, backwards::Bool = true, verbose::Bool = false)
    verbose && println("Loading fast-ion distribution from .h5 file... ")
    F_EpRz, E_array, p_array, R_array, z_array = h5to4D(filepath_distr, backwards=backwards, verbose=verbose)

    return toSampleReady(F_EpRz, E_array, p_array, R_array, z_array; interp=interp, nE_ps=nE_ps, np_ps=np_ps, nR_ps=nR_ps, nz_ps=nz_ps, verbose=verbose)
end

"""
    jld2toSampleReady(filepath_distr)
    jld2toSampleReady(filepath_distr; verbose=false) # default

Load the fast-ion distribution from a .jld2 file and forward to function toSampleReady().

--- Input:
filepath_distr - The fast-ion distribution in .jld2 file format - String
verbose - If set to true, you will get a talkative function indeed - Bool
"""
function jld2toSampleReady(filepath_distr::AbstractString; verbose::Bool = false, kwargs...)
    verbose && println("Loading fast-ion distribution from .jld2 file... ")
    F_EpRz, E_array, p_array, R_array, z_array = JLD2to4D(filepath_distr; verbose=verbose)
    
    return toSampleReady(F_EpRz, E_array, p_array, R_array, z_array; verbose=verbose, kwargs...)
end

"""
    toSampleReady(F_EpRz, E_array, p_array, R_array, z_array)
    toSampleReady(-||-, nE_ps=202, np_ps=57, nR_ps=62, nz_ps=64, verbose=true)

Interpolate the fast-ion distribution and grid arrays if specified, then compute the necessary
quantities needed to be able to sample from the fast-ion distribution with the spaghettification
method (for example in calcSpec.jl). Then return those quantities.

--- Input:
F_EpRz - The 4D fast-ion distribution - Array{Float64,4}
E_array - The 1D array containing the energy grid points - Array{Float64,1}
p_array - The 1D array containing the pitch grid points - Array{Float64,1}
R_array - The 1D array containing the R grid points - Array{Float64,1}
z_array - The 1D array containing the z grid points - Array{Float64,1}
interp - If true, then fast-ion distribution will be interpolated onto grid with same end-points but dimension (nE_ps, np_ps, nR_ps, nz_ps) - Bool
nE_ps - Energy grid dimension to be interpolated onto, if interp - Int64
np_ps - Pitch grid dimension to be interpolated onto, if interp - Int64
nR_ps - R grid dimension to be interpolated onto, if interp - Int64
nz_ps - z grid dimension to be interpolated onto, if interp - Int64
verbose - If set to true, you will get a talkative function indeed - Bool
"""
function toSampleReady(F_EpRz::Array{Float64,4}, E_array::Array{Float64,1}, p_array::Array{Float64,1}, R_array::Array{Float64,1}, z_array::Array{Float64,1}; interp=false, nE_ps::Int64=100, np_ps::Int64=51, nR_ps::Int64=60, nz_ps::Int64=61, verbose::Bool = false)
    if interp
        energyq = range(E_array[1], E_array[end], length=nE_ps) # Set the energy query points
        pitchq = range(p_array[1], p_array[end], length=np_ps) # Etc
        Rq = range(R_array[1], R_array[end], length=nR_ps) # Etc
        zq = range(z_array[1], z_array[end], length=nz_ps) # Etc
        verbose && println("Interpolating Julia fast-ion distribution... ")
        F_EpRz = interpFps(F_EpRz, E_array, p_array, R_array, z_array, energyq, pitchq, Rq, zq) # Interpolate
        E_array = energyq # The energy query points are now the energy points
        p_array = pitchq # Etc
        R_array = Rq
        z_array = zq
    end

    fr = F_EpRz.*reshape(R_array,(1,1,length(R_array),1))
    verbose && println("Computing 4D vols... ")
    dvols = get4DVols(E_array, p_array, R_array, z_array)
    nfast = sum(fr.*dvols)*(2*pi)

    # Checking if units of R,z is meters or centimeters. Correct to meters if units is centimeters
    # PLEASE NOTE! IT IS VERY IMPORTANT TO DO THIS STEP AFTER YOU HAVE COMPUTED fr, dvols and nfast
    # OTHERWISE, THE TOTAL NUMBER OF FAST IONS WILL BE WRONG, BECAUSE THE ORIGINAL F_EpRz DATA WAS
    # IN UNITS OF PER CENTIMETER
    if maximum(R_array)>100.0 # Assume no tokamak has a major radius larger than 100 meters...
        verbose && println("Converting R- and z-arrays from centimeters to meters... ")
        R_array = R_array./100.0
        z_array = z_array./100.0
    end
    dE_array = vcat(abs.(diff(E_array)),abs(E_array[end]-E_array[end-1]))
    dp_array = vcat(abs.(diff(p_array)),abs(p_array[end]-p_array[end-1]))
    dR_array = vcat(abs.(diff(R_array)),abs(R_array[end]-R_array[end-1]))
    dz_array = vcat(abs.(diff(z_array)),abs(z_array[end]-z_array[end-1]))

    dims = size(fr) # Tuple
    verbose && println("Computing cumulative sum vector... ")
    frdvols_cumsum_vector = cumsum(vec(fr .*dvols)) # Vector. Reaaally long vector.
    subs = CartesianIndices(dims) # 4D matrix

    return frdvols_cumsum_vector, subs, E_array, p_array, R_array, z_array, dE_array, dp_array, dR_array, dz_array, nfast

end


"""
    get_orbel_volume(myOrbitGrid, equidistant)

Get the orbit element volume of this specific orbit-grid.
If equidistant=false, then return a 3D array with all the orbit volumes.
Otherwise, just a scalar.
"""
function get_orbel_volume(og::OrbitGrid, equidistant::Bool)

    if equidistant
        dRm = abs(og.r[2]-og.r[1])
        dE = abs(og.energy[2]-og.energy[1])
        dpm = abs(og.pitch[2]-og.pitch[1])
        dO = dE*dpm*dRm
        return dO
    else
        return get3DVols(og.energy, og.pitch, og.r)
    end
end


"""
    get3DDiffs(E, pm, Rm)

Calculate and return 3D-array of the diffs of orbit-space coordinates. Assume edge diff to be same as next-to-edge diff.
"""
function get3DDiffs(E::AbstractVector, pm::AbstractVector, Rm::AbstractVector)
    dE = vcat(abs.(diff(E)),abs(E[end]-E[end-1]))
    dE = reshape(dE,length(dE),1,1)
    dE3D = repeat(dE,1,length(pm),length(Rm))
    dpm = vcat(abs.(diff(pm)),abs(pm[end]-pm[end-1]))
    dpm = reshape(dpm,1,length(dpm),1)
    dpm3D = repeat(dpm,length(E),1,length(Rm))
    dRm = vcat(abs.(diff(Rm)),abs(Rm[end]-Rm[end-1]))
    dRm = reshape(dRm,1,1,length(dRm))
    dRm3D = repeat(dRm,length(E),length(pm),1)

    return dE3D, dpm3D, dRm3D
end


"""
    get3DVols(E, pm, Rm)

This function will calculate the volume of all the hyper-voxels pertaining to the 3D grid. It assumes the hyper-voxels pertaining to the
upper-end (edge) grid-points have the same volumes as the hyper-voxels just inside of them. It return a 3D array, containing all the hyper-voxel
volumes. The 3D array will have size()=(length(E), length(pm), length(Rm)). The function assumes a rectangular 3D grid.
"""
function get3DVols(E::AbstractVector, pm::AbstractVector, Rm::AbstractVector)

    # Safety-check to ensure vectors
    if !(1==length(size(E))==length(size(pm))==length(size(Rm)))
        throw(ArgumentError("E, pm, Rm inputs are not all vectors. Please correct and re-try."))
    end

    dE3D, dpm3D, dRm3D = get3DDiffs(E, pm, Rm)
    dvols = dE3D .*dpm3D .*dRm3D
    return dvols
end


"""
    get4DDiffs(E, p, R, z)

Calculate and return 4D-arrays of the diffs of particle-space coordinates. Assume edge diff to be same as next-to-edge diff.
"""
function get4DDiffs(E::AbstractVector, p::AbstractVector, R::AbstractVector, z::AbstractVector)
    dR = vcat(abs.(diff(R)),abs(R[end]-R[end-1]))
    dR = reshape(dR,1,1,length(dR),1)
    dR4D = repeat(dR,length(E),length(p),1,length(z))
    dz = vcat(abs.(diff(z)),abs(z[end]-z[end-1]))
    dz = reshape(dz,1,1,1,length(dz))
    dz4D = repeat(dz,length(E),length(p),length(R),1)
    dE = vcat(abs.(diff(E)),abs(E[end]-E[end-1]))
    dE = reshape(dE,length(dE),1,1,1)
    dE4D = repeat(dE,1,length(p),length(R),length(z))
    dp = vcat(abs.(diff(p)),abs(p[end]-p[end-1]))
    dp = reshape(dp,1,length(dp),1,1)
    dp4D = repeat(dp,length(E),1,length(R),length(z))

    return dE4D, dp4D, dR4D, dz4D
end

"""
    get4DVols(E, p, R, z)

This function will calculate the volume of all the hyper-voxels pertaining to the 4D grid. It assumes the hyper-voxels pertaining to the
upper-end (edge) grid-points have the same volumes as the hyper-voxels just inside of them. It return a 4D array, containing all the hyper-voxel
volumes. The 4D array will have size()=(length(E), length(p), length(R), length(z)). The function assumes a rectangular 4D grid.
"""
function get4DVols(E::AbstractVector, p::AbstractVector, R::AbstractVector, z::AbstractVector)

    # Safety-check to ensure vectors
    if !(1==length(size(E))==length(size(p))==length(size(R))==length(size(z)))
        throw(ArgumentError("E, p, R, z inputs are not all vectors. Please correct and re-try."))
    end

    dE4D, dp4D, dR4D, dz4D = get4DDiffs(E, p, R, z)
    dvols = dE4D .*dp4D .*dR4D .*dz4D
    return dvols
end

"""
    sample_f_OWCF(fr,dvols,energy,pitch,R,z,n=100_000)

This function is the OWCF version of the original function sample_f which already exists in the OrbitTomography.jl package.
This version has the ability to take non-equidistant 4D grid-points into consideration. w is energy, x is pitch, y is R and z is z.
"""
function sample_f_OWCF(fr::Array{T,N}, dvols, w::AbstractVector, x::AbstractVector, y::AbstractVector, z::AbstractVector; n::Int64=100_000) where {T,N}

    inds = OrbitTomography.sample_array(fr.*dvols,n)

    dw = vcat(abs.(diff(w)),abs(w[end]-w[end-1]))
    dx = vcat(abs.(diff(x)),abs(x[end]-x[end-1]))
    dy = vcat(abs.(diff(y)),abs(y[end]-y[end-1]))
    dz = vcat(abs.(diff(z)),abs(z[end]-z[end-1]))

    r = rand(N,n) .- 0.5
    xx = zeros(n)
    yy = zeros(n)
    zz = zeros(n)
    ww = zeros(n)

    @inbounds for i=1:n
        ww[i] = max(w[inds[1,i]] + r[1,i]*dw[inds[1,i]], 0.0)
        xx[i] = x[inds[2,i]] + r[2,i]*dx[inds[2,i]]
        yy[i] = y[inds[3,i]] + r[3,i]*dy[inds[3,i]]
        zz[i] = z[inds[4,i]] + r[4,i]*dz[inds[4,i]]
    end

    return ww, xx, yy, zz
end

"""
    class2int(o)
    class2int(o; plot=false, distinguishLost=true, distinguishIncomplete=true)

Take an orbit, examine its class (stagnation, trapped, co-passing etc) and return 
the appropriate integer. The integers are as follows
1 - stagnation
2 - trapped
3 - co-passing
4 - counter-passing
5 - potato
9 - invalid

Depending on the keyword arguments, give orbit appropriate integer
- If plot, give counter-stagnation orbits integer 8. Else, give 6
- If distinguishLost, give lost orbits integer 7. Else, give 9
- If distinguishIncomplete and plot, give incomplete orbits integer 6. Otherwise, give 9

Keyword argument 'plot' is meant to ensure correct color map when using output for OWCF apps.
If 'plot' is set to false, then you can simply feed all integers into a 6 element array and 
have the valid orbits (stagnation, trapped, co-passing, counter-passing, potato and counter-stagnation) in bins 1 to 6.
"""
function class2int(M::AbstractEquilibrium, o::GuidingCenterOrbits.Orbit;plot::Bool=true, distinguishLost::Bool=false, distinguishIncomplete::Bool=false)
    if (o.class == :lost) && distinguishLost
        return 7
    elseif (o.class == :incomplete) && distinguishIncomplete && plot
        return 6
    elseif o.class == :trapped
        return 2
    elseif o.class == :co_passing
        return 3
    elseif (o.class == :stagnation && o.coordinate.r>=M.axis[1]) # Regular stagnation orbit
        return 1
    elseif (o.class == :stagnation && o.coordinate.r<M.axis[1] && plot) # Counterstagnation orbit
        return 8
    elseif (o.class == :stagnation && o.coordinate.r<M.axis[1] && !plot) # Counterstagnation orbit, but treat as normal stagnation orbit
        return 6
    elseif o.class == :potato
        return 5
    elseif o.class == :ctr_passing
        return 4
    else
        return 9
    end
end

"""
    map_orbits_OWCF(og, f, equidistant)

This function is the OWCF version of the original function map_orbits which already exists in the OrbitTomography.jl package.
It has the ability to take non-equidistant 3D grid-points into account.
"""
function map_orbits_OWCF(grid::OrbitGrid, f::Vector, equidistant::Bool; weights::Bool=false)
    if equidistant
        return OrbitTomography.map_orbits(grid,f) # Use the normal map_orbits function
    else
        if length(grid.counts) != length(f)
            throw(ArgumentError("Incompatible sizes"))
        end

        if !(weights)
            dorb = get_orbel_volume(grid, equidistant) # Get the volumes of the voxels of the orbit grid
            return [i == 0 ? zero(f[1]) : f[i]/(grid.counts[i]*dorb[i]) for i in grid.orbit_index] # This line is complicated...
        else
            return [i == 0 ? zero(f[1]) : f[i]/(grid.counts[i]*1.0) for i in grid.orbit_index] # This line is complicated...
        end
    end
end

"""
    ps2os_streamlined(filepath_distr, filepath_equil, og, numOsamples)
    ps2os_streamlined(-||-, verbose=true, distr_dim=[99,43,54,56], btipsign=-1, GCP=GCProton)

This function is a streamlined version of sampling particle-space (PS) and converting those samples to orbit-space (OS).
It allows transformation from (E,p,R,z) to (E,pm,Rm) directly from filepath_distr
"""
function ps2os_streamlined(filepath_distr::AbstractString, filepath_equil::AbstractString, og::OrbitGrid; numOsamples::Int64, verbose::Bool=false, kwargs...)

    if verbose
        println("")
        println("Loading f, energy, pitch, R and z from file... ")
    end
    f, energy, pitch, R, z = h5to4D(filepath_distr, verbose = verbose) # Load fast-ion distribution

    return ps2os_streamlined(f,energy,pitch,R,z,filepath_equil, og; numOsamples=numOsamples, verbose=verbose, kwargs...)
end

"""
    ps2os_point()
    ps2os_point()

Bla bla bla
"""
function ps2os_point(M::AbstractEquilibrium, E::T, p::T, R::T, z::T, ) where {T}
    
end

"""
Continuation of the ps2os_streamlined function above.
"""
function ps2os_streamlined(F_EpRz::Array{Float64,4}, energy::AbstractVector, pitch::AbstractVector, R::AbstractVector, z::AbstractVector, filepath_equil::AbstractString, og::OrbitGrid; numOsamples::Int64, verbose::Bool=false, distr_dim = [], flip_F_EpRz_pitch::Bool=false, GCP = GCDeuteron, distributed::Bool=true, nbatch::Int64 = 1_000_000, clockwise_phi::Bool=false, kwargs...)

    verbose && println("Loading the magnetic equilibrium... ")
    if ((split(filepath_equil,"."))[end] == "eqdsk") || ((split(filepath_equil,"."))[end] == "geqdsk")
        M, wall = read_geqdsk(filepath_equil,clockwise_phi=clockwise_phi) # Assume counter-clockwise phi-direction
        jdotb = M.sigma # The sign of the dot product between the plasma current and the magnetic field
    else # Otherwise, assume magnetic equilibrium is a saved .jld2 file
        myfile = jldopen(filepath_equil,false,false,false,IOStream)
        M = myfile["S"]
        wall = myfile["wall"]
        close(myfile)
        jdotb = (M.sigma_B0)*(M.sigma_Ip)
    end

    if !isempty(distr_dim)
        if verbose
            println("Interpolating in 4D... ")
        end
        if length(distr_dim)!=4
            throw(ArgumentError("Length of distr_dim MUST be 4. Please correct and try again."))
        end

        energyq = range(energy[1], energy[end], length=distr_dim[1]) # Set the energy query points
        pitchq = range(pitch[1], pitch[end], length=distr_dim[2]) # Etc
        Rq = range(R[1], R[end], length=distr_dim[3]) # Etc
        zq = range(z[1], z[end], length=distr_dim[4]) # Etc
        F_EpRz = interpFps(F_EpRz, energy, pitch, R, z, energyq, pitchq, Rq, zq) # Interpolate

        energy = energyq # The energy query points are now the energy points
        pitch = pitchq # Etc
        R = Rq
        z = zq
    end

    # This operation will be necessary when working with old .h5-file distributions (16th April 2021 and earlier).
    # This is because a (M.sigma*) multiplication operation was removed in the get_orbit() calls in the helper functions (etc).
    if flip_F_EpRz_pitch
        F_EpRz = flip_along_p(F_EpRz)
    end

    return ps2os(M, wall, F_EpRz, energy, pitch, R, z, og; numOsamples=numOsamples, verbose=verbose, kwargs...)
end

"""
    ps2os(M, wall, F_EpRz, energy, pitch, R, z, og)
    ps2os(-||-; numOsamples, verbose=false, GCP = GCDeuteron, distributed=true, nbatch=1_000_000, saveProgress=true, performance=true, kwargs...)

Take the (E,p,R,z) fast-ion distribution and use Monte-Carlo sampling to transform it to (E,pm,Rm) orbit space.
- The number of Monte-Carlo samples is defined by numOsamples
- If verbose, the function will talk a lot!
- GCP defines the guiding-center particle species
- distributed controls whether to use multi-core processing when Monte-Carlo sampling
- The Monte-Carlo sampling will be saved every 'nbatch' number of samples. Good to use, in case something goes wrong and terminates early
- if saveProgress, a progress bar will be displayed when sampling
- If performance, then performance sampling will be used. Highly recommended! Normal sampling will be deprecated in later versions of the OWCF
"""
function ps2os(M::AbstractEquilibrium, wall::Boundary, F_EpRz::Array{Float64,4}, energy::AbstractVector, pitch::AbstractVector, R::AbstractVector, z::AbstractVector, og::OrbitGrid; numOsamples::Int64, verbose::Bool=false, GCP = GCDeuteron, distributed::Bool=true, nbatch::Int64 = 1_000_000, saveProgress::Bool=true, performance::Bool=true, kwargs...)

    if verbose
        println("Acquiring fr, dvols and nfast... ")
    end
    dR = vcat(abs.(diff(R)),abs(R[end]-R[end-1]))
    dR = reshape(dR,1,1,length(dR),1)
    dR4D = repeat(dR,length(energy),length(pitch),1,length(z))
    dz = vcat(abs.(diff(z)),abs(z[end]-z[end-1]))
    dz = reshape(dz,1,1,1,length(dz))
    dz4D = repeat(dz,length(energy),length(pitch),length(R),1)
    dE = vcat(abs.(diff(energy)),abs(energy[end]-energy[end-1]))
    dE = reshape(dE,length(dE),1,1,1)
    dE4D = repeat(dE,1,length(pitch),length(R),length(z))
    dp = vcat(abs.(diff(pitch)),abs(pitch[end]-pitch[end-1]))
    dp = reshape(dp,1,length(dp),1,1)
    dp4D = repeat(dp,length(energy),1,length(R),length(z))
    fr = F_EpRz.*reshape(R,(1,1,length(R),1))
    dvols = get4DVols(energy, pitch, R, z)
    nfast = sum(fr.*dE4D.*dp4D.*dR4D.*dz4D)*(2*pi)

    # Keeping memory usage to minimum
    F_EpRz = nothing
    dE = nothing
    dp = nothing
    dR = nothing
    dz = nothing
    dE4D = nothing
    dp4D = nothing
    dR4D = nothing
    dz4D = nothing

    # Checking if units of R,z is meters or centimeters. Correct to meters if units is centimeters
    if maximum(R)>100.0 # Assume no tokamak has a major radius larger than 100 meters...
        R = R./100.0
        z = z./100.0
    end

    #################################################################################
    # Handle re-start of ps2os-transformation process, if abruptly cancelled
    if isfile("ps2os_progress.jld2")
        myfile = jldopen("ps2os_progress.jld2",false,false,false,IOStream)
        numOsamples_sofar = deepcopy(myfile["numOsamples"])
        result_sofar = deepcopy(myfile["F_os"])
        numOsamples = numOsamples - numOsamples_sofar

        if length(og.counts) != length(result_sofar)
            error("Loaded orbit-space fast-ion distribution from progress-file is non-compatible with specified orbit grid. Please correct orbit grid, provide valid progress-file or remove progress-file from work directory.")
        end
        close(myfile)
    else
        numOsamples_sofar = 0
        result_sofar = zeros(size(og.counts))
    end
    #################################################################################

    #################################################################################
    verbose && print("Performance = $(performance) ==>  ")
    if performance
        verbose && println("Computing samples efficiently... ")
        dims = size(fr) # Tuple
        frdvols_cumsum_vector = cumsum(vec(fr .*dvols)) # Vector. Reaaally long vector.
        subs = CartesianIndices(dims) # 4D matrix
        fr = nothing # Memory efficiency
        dvols = nothing # Memory efficiency
        return ps2os_performance(M, wall, frdvols_cumsum_vector, subs, nfast, energy, pitch, R, z, og; numOsamples=numOsamples, numOsamples_sofar=numOsamples_sofar, result_sofar=result_sofar, distributed=distributed, GCP=GCP, saveProgress=saveProgress, verbose=verbose, kwargs...)
    end
    verbose && println("Computing samples the old way... ")
    #################################################################################

    # PLEASE NOTE! If performance = true, then the algorithm will never make it down here.
    #################################################################################
    if verbose
        println("Starting the sampling process... ")
    end
    #################################################################################
    if distributed # If paralell computing...
        subdivide = false # Needed to handle scenario when nbatch is larger than numOsamples from start
        while numOsamples > nbatch # Sample in chunks of nbatch, for safety
            subdivide = true
            numOsamples = numOsamples - nbatch
            if verbose
                println("Samples left: $(numOsamples)")
            end

            result_p = sample_helper(M, nbatch, fr, dvols, energy, pitch, R, z, og; wall=wall, GCP=GCP, kwargs...)
            result_sofar .+= result_p
            numOsamples_sofar += nbatch
            if saveProgress
                rm("ps2os_progress.jld2", force=true) #clear the previous file
                myfile = jldopen("ps2os_progress.jld2",true,true,false,IOStream)
                write(myfile,"F_os",result_sofar)
                write(myfile,"numOsamples",numOsamples_sofar)
                close(myfile)
            end
        end
        if verbose
            println("(Rest) Samples left: $(numOsamples)")
        end
        result_rest = sample_helper(M, numOsamples, fr, dvols, energy, pitch, R, z, og; wall=wall, GCP=GCP, kwargs...)
        numOsamples_rest = numOsamples

        if subdivide
            result = result_sofar + result_rest
            numOsamples = numOsamples_sofar + numOsamples_rest
        else
            result = result_rest
            numOsamples = numOsamples_rest
        end
    else # ...if not parallel computing... I wish you good luck.
        dE_os_end = abs((og.energy)[end]-(og.energy)[end-1])
        dE_os_1 = abs((og.energy)[2]-(og.energy)[1])
        
        for i=numOsamples_sofar+1:numOsamples
            if verbose
                println("Sample number: $(i)")
            end
            E_sample, p_sample, R_sample, z_sample = sample_f_OWCF(fr, dvols, energy, pitch, R, z, n=1) # Returns 4 1-element arrays

            # CHECK IF IT'S A GOOD SAMPLE
            good_sample = checkIfGoodSample(E_sample[1], p_sample[1], R_sample[1], z_sample[1], energy, pitch, R, z)

            if good_sample
                o = get_orbit(M,GCP(E_sample[1],p_sample[1],R_sample[1],z_sample[1]); store_path=false, wall=wall, kwargs...)
                if (o.coordinate.energy <= (og.energy[end]+dE_os_end/2) && o.coordinate.energy >= (og.energy[1]-dE_os_1/2)) # Make sure it's within the energy bounds (+one half grid cell)
                    F_os_i = bin_orbits(og,vec([o.coordinate]),weights=vec([1.0]))
                else
                    F_os_i = zeros(length(og.counts))
                end
            else
                F_os_i = zeros(length(og.counts))
            end

            result_sofar .+= F_os_i

            if (i%nbatch)==0 && saveProgress # Every nbatch sample, save
                rm("ps2os_progress.jld2", force=true) #clear the previous file
                myfile = jldopen("ps2os_progress.jld2",true,true,false,IOStream)
                write(myfile,"F_os",result_sofar)
                write(myfile,"numOsamples",i)
                close(myfile)
            end
        end
        result = result_sofar
    end
    if verbose
        println("Number of good samples/All samples: $(sum(result)/numOsamples)")
    end

    rm("ps2os_progress.jld2", force=true) # Finally, remove the file that is no longer needed

    return result, nfast
end

"""
    ps2os_performance(M, wall, fr, dvols, nfast, energy, pitch, R, z, og)
    ps2os_performance(-||-; numOsamples, numOsamples_sofar=0, result_sofar=zeros(size(og.counts))), distributed=true, GCP=GCDeuteron, saveProgess=true, nbatch=1_000_000, verbose=false, kwargs...)

The performance version of part of ps2os(). This function will likely completely replace ps2os() in the near future. It computes necessary quantities
once instead of for every sample (as ps2os() does). Such as the element-wise product of fr and dvols, and its cumulative sum.
M - The magnetic equilibrium struct
wall - The tokamak first wall
frdvols_cumsum_vector - A vector. The cumulative sum of the (E,p,R,z) fast-ion distribution
subs - Cartesian indices with the corresponding (E,p,R,z) points for the elements in frdvols_cumsum_vector
nfast - The total number of fast ions
energy - The energy grid points
pitch - The pitch grid points
R - The major radius grid points
z - The vertical grid points
og - The orbit grid struct
;
numOsamples - The total number of Monte-Carlo samples (includes numOsamples_sofar)
numOsamples_sofar - The number of Monte-Carlo samples sampled so far
result_sofar - The (E,pm,Rm) fast-ion distribution in 1D compressed vector format. So far, having completed numOsamples_sofar number of samples
distributed - If true, multi-core processing will be used
GCP - The guiding-center struct representing the fast-ion species
saveProgress - If true, Monte-Carlo sampling process will be saved every nbatch samples
nbatch - Compute the Monte-Carlo samples in batches, to optimize computation efficiency and enable subsequent progress saving
verbose - If true, the function will talk a lot
"""
function ps2os_performance(M::AbstractEquilibrium,
                            wall::Boundary,
                            frdvols_cumsum_vector::AbstractVector,
                            subs::CartesianIndices{4,NTuple{4,Base.OneTo{Int64}}},
                            nfast::Union{Int64,Float64},
                            energy::AbstractVector,
                            pitch::AbstractVector,
                            R::AbstractVector,
                            z::AbstractVector,
                            og::OrbitGrid;
                            numOsamples::Int64,
                            numOsamples_sofar::Int64=0,
                            result_sofar=zeros(size(og.counts)),
                            distributed::Bool=true,
                            GCP=GCDeuteron,
                            saveProgress::Bool=true,
                            nbatch::Int64 = 1_000_000,
                            verbose::Bool=false,
                            kwargs...)
    verbose && println("Pre-computing difference vectors... ")
    dE_vector = vcat(abs.(diff(energy)),abs(energy[end]-energy[end-1]))
    dp_vector = vcat(abs.(diff(pitch)),abs(pitch[end]-pitch[end-1]))
    dR_vector = vcat(abs.(diff(R)),abs(R[end]-R[end-1]))
    dz_vector = vcat(abs.(diff(z)),abs(z[end]-z[end-1]))

    verbose && println("Going into the distributed part... ")
    if distributed
        subdivide = false
        while numOsamples > nbatch
            subdivide = true
            verbose && println("Samples left: $(numOsamples)")
            numOsamples = numOsamples - nbatch
            result_p = performance_helper(M, nbatch, frdvols_cumsum_vector, subs, dE_vector, dp_vector, dR_vector, dz_vector, energy, pitch, R, z, og; wall=wall, GCP=GCP, kwargs...)
            result_sofar .+= result_p
            numOsamples_sofar += nbatch
            if saveProgress
                rm("ps2os_progress.jld2", force=true) #clear the previous file
                myfile = jldopen("ps2os_progress.jld2",true,true,false,IOStream)
                write(myfile,"F_os",result_sofar)
                write(myfile,"numOsamples",numOsamples_sofar)
                close(myfile)
            end
        end
        verbose && println("(Rest) Samples left: $(numOsamples)")
        result_rest = performance_helper(M, numOsamples, frdvols_cumsum_vector, subs, dE_vector, dp_vector, dR_vector, dz_vector, energy, pitch, R, z, og; wall=wall, GCP=GCP, kwargs...)
        numOsamples_rest = numOsamples

        if subdivide
            result = result_sofar + result_rest
            numOsamples = numOsamples_sofar + numOsamples_rest
        else
            result = result_rest
            numOsamples = numOsamples_rest
        end
    else # ...if not parallel computing... I wish you good luck.
        dE_os_end = abs((og.energy)[end]-(og.energy)[end-1])
        dE_os_1 = abs((og.energy)[2]-(og.energy)[1])

        for i=numOsamples_sofar+1:numOsamples
            if verbose
                println("Sample number: $(i)")
            end
            # Sample
            p = rand()*frdvols_cumsum_vector[end]
            j = searchsortedfirst(frdvols_cumsum_vector,p,Base.Order.Forward)
            inds = collect(Tuple(subs[j])) # First sample
            r = rand(4) .- 0.5 # 4 stands for the number of dimensions. 0.5 to sample within a hypercube
            E_sample = max(energy[inds[1]] + r[1]*dE_vector[inds[1]], 0.0)
            p_sample = pitch[inds[2]] + r[2]*dp_vector[inds[2]]
            R_sample = R[inds[3]] + r[3]*dR_vector[inds[3]]
            z_sample = z[inds[4]] + r[4]*dz_vector[inds[4]]

            # CHECK IF IT'S A GOOD SAMPLE
            good_sample = checkIfGoodSample(E_sample, p_sample, R_sample, z_sample, energy, pitch, R, z)

            if good_sample
                o = get_orbit(M,GCP(E_sample,p_sample,R_sample,z_sample); store_path=false, wall=wall, kwargs...)
                if (o.coordinate.energy <= (maximum(og.energy)+dE_os_end/2) && o.coordinate.energy >= (minimum(og.energy)-dE_os_1/2)) # Make sure it's within the energy bounds (+one half grid cell)
                    F_os_i = bin_orbits(og,vec([o.coordinate]),weights=vec([1.0]))
                else
                    F_os_i = zeros(length(og.counts))
                end
            else
                F_os_i = zeros(length(og.counts))
            end

            result_sofar .+= F_os_i

            if (i%nbatch)==0 && saveProgress # Every nbatch sample, save
                rm("ps2os_progress.jld2", force=true) #clear the previous file
                myfile = jldopen("ps2os_progress.jld2",true,true,false,IOStream)
                write(myfile,"F_os",result_sofar)
                write(myfile,"numOsamples",i)
                close(myfile)
            end
        end
        result = result_sofar
    end

    if verbose
        println("Number of good samples/All samples: $(sum(result)/numOsamples)")
    end
    rm("ps2os_progress.jld2", force=true) # As in ps2os(), remove the progress file that is no longer needed

    return result, nfast
end

"""
    sample_helper(M, numOsamples, fr, dvols, energy, pitch, R, z, og)

Help the function ps2os() with acquiring orbit samples when parallel computations are desired.
This is to enable the sampling process to be saved regularly when calculating a large number of samples.
If the sampling process is not saved, then progress will be lost when the super-user of the HPC terminates
the sampling process early, due to misinterpretation of Julia's way of distributed computing.
"""
function sample_helper(M::AbstractEquilibrium, numOsamples::Int64, fr::AbstractArray, dvols::AbstractArray, energy::AbstractVector, pitch::AbstractVector, R::AbstractVector, z::AbstractVector, og::OrbitGrid; wall::Union{Nothing,Boundary}, GCP=GCDeutron, visualizeProgress::Bool=false, kwargs...)

    dE_os_end = abs((og.energy)[end]-(og.energy)[end-1])
    dE_os_1 = abs((og.energy)[2]-(og.energy)[1])

    if numOsamples>0.0 # If there are actually a non-zero number of samples left to sample...
        if visualizeProgress
            p = Progress(numOsamples) # Define the progress bar
            channel = RemoteChannel(()->Channel{Bool}(numOsamples),1) # Define the channel from which the progress bar draws data
            result = fetch(@sync begin # Start the distributed computational process, fetch result when done
                @async while take!(channel) # An asynchronous process, with no need for sync, since it simply displays the progress bar
                    ProgressMeter.next!(p)
                end

                @async begin # No internal syncronization needed here either, only external sync needed
                    F_os = @distributed (+) for i=1:numOsamples # Compute one result, and reduce (add) it to a resulting vector F_os

                        # Sample
                        E_sample, p_sample, R_sample, z_sample = sample_f_OWCF(fr, dvols, energy, pitch, R, z, n=1) # Returns 4 1-element arrays
                        # CHECK IF IT'S A GOOD SAMPLE
                        good_sample = checkIfGoodSample(E_sample[1], p_sample[1], R_sample[1], z_sample[1], energy, pitch, R, z)
                        if good_sample
                            o = get_orbit(M,GCP(E_sample[1],p_sample[1],R_sample[1],z_sample[1]); store_path=false, wall=wall, kwargs...) # Calculate the orbit
                            if (o.coordinate.energy <= (maximum(og.energy)+dE_os_end/2) && o.coordinate.energy >= (minimum(og.energy)-dE_os_1/2)) # Make sure it's within the energy bounds (+one half grid cell)
                                F_os_i = bin_orbits(og,vec([o.coordinate]),weights=vec([1.0])) # Bin to the orbit grid
                            else
                                F_os_i = zeros(length(og.counts)) # Otherwise, zero
                            end
                        else
                            F_os_i = zeros(length(og.counts)) # Otherwise, zero
                        end
                        put!(channel,true) # Update the progress bar
                        F_os_i # Declare F_os_i as result to add to F_os
                    end
                    put!(channel,false) # Update progress bar
                    F_os # Delcare F_os as done/result, so it can be fetched
                end
            end)
        else
            result = @distributed (+) for i=1:numOsamples
                E_sample, p_sample, R_sample, z_sample = sample_f_OWCF(fr, dvols, energy, pitch, R, z, n=1) # Returns 4 1-element arrays
                good_sample = checkIfGoodSample(E_sample[1], p_sample[1], R_sample[1], z_sample[1], energy, pitch, R, z)
                if good_sample
                    o = get_orbit(M,GCP(E_sample[1],p_sample[1],R_sample[1],z_sample[1]); store_path=false, wall=wall, kwargs...) # Calculate the orbit
                    if (o.coordinate.energy <= (maximum(og.energy)+dE_os_end/2) && o.coordinate.energy >= (minimum(og.energy)-dE_os_1/2)) # Make sure it's within the energy bounds (+one half grid cell)
                        F_os_i = bin_orbits(og,vec([o.coordinate]),weights=vec([1.0])) # Bin to the orbit grid
                    else
                        F_os_i = zeros(length(og.counts)) # Otherwise, zero
                    end
                else
                    F_os_i = zeros(length(og.counts)) # Otherwise, zero
                end
                F_os_i # Declare F_os_i as result to add to result
            end
        end

        return result
    else # ...otherwise just return a zero vector with the length of the number of valid orbits
        return zeros(length(og.counts))
    end
end

"""
    performance_helper(M, nbatch, frdvols_cumsum_vector, subs, dE_vector, dp_vector, dR_vector, dz_vector, energy, pitch, R, z, og; wall=wall, GCP=GCP, kwargs...)

Help the function ps2os_performance() with acquiring orbit samples when parallel computations are desired.
This is to enable the sampling process to be saved regularly when calculating a large number of samples.
If the sampling process is not saved, then progress will be lost when the super-user of the HPC terminates
the sampling process early, due to misinterpretation of Julia's way of distributed computing.
"""
function performance_helper(M::AbstractEquilibrium, numOsamples::Int64, frdvols_cumsum_vector::AbstractVector, subs::CartesianIndices{4,NTuple{4,Base.OneTo{Int64}}}, dE_vector::AbstractVector, dp_vector::AbstractVector, dR_vector::AbstractVector, dz_vector::AbstractVector, energy::AbstractVector, pitch::AbstractVector, R::AbstractVector, z::AbstractVector, og::OrbitGrid; wall::Union{Nothing,Boundary{Float64}}, GCP=GCDeuteron, visualizeProgress::Bool=false, kwargs...)

    dE_os_end = abs((og.energy)[end]-(og.energy)[end-1])
    dE_os_1 = abs((og.energy)[2]-(og.energy)[1])

    if numOsamples>0 # If there are actually a non-zero number of samples left to sample...
        if visualizeProgress
            p = Progress(numOsamples) # Define the progress bar
            channel = RemoteChannel(()->Channel{Bool}(numOsamples),1) # Define the channel from which the progress bar draws data
            result = fetch(@sync begin # Start the distributed computational process, fetch result when done
                @async while take!(channel) # An asynchronous process, with no need for sync, since it simply displays the progress bar
                    ProgressMeter.next!(p)
                end

                @async begin # No internal syncronization needed here either, only external sync needed
                    F_os = @distributed (+) for i=1:numOsamples # Compute one result, and reduce (add) it to a resulting vector F_os

                        # Sample
                        p = rand()*frdvols_cumsum_vector[end]
                        j = searchsortedfirst(frdvols_cumsum_vector,p,Base.Order.Forward)
                        inds = collect(Tuple(subs[j])) # First sample
                        r = rand(4) .- 0.5 # 4 stands for the number of dimensions. 0.5 to sample within a hypercube
                        E_sample = max(energy[inds[1]] + r[1]*dE_vector[inds[1]], 0.0)
                        p_sample = pitch[inds[2]] + r[2]*dp_vector[inds[2]]
                        R_sample = R[inds[3]] + r[3]*dR_vector[inds[3]]
                        z_sample = z[inds[4]] + r[4]*dz_vector[inds[4]]

                        # CHECK IF IT'S A GOOD SAMPLE
                        good_sample = checkIfGoodSample(E_sample, p_sample, R_sample, z_sample, energy, pitch, R, z)
                        if good_sample
                            o = get_orbit(M,GCP(E_sample,p_sample,R_sample,z_sample); store_path=false, wall=wall, kwargs...) # Calculate the orbit
                            if (o.coordinate.energy <= (maximum(og.energy)+dE_os_end/2) && o.coordinate.energy >= (minimum(og.energy)-dE_os_1/2)) # Make sure it's within the energy bounds (+one half grid cell)
                                F_os_i = bin_orbits(og,vec([o.coordinate]),weights=vec([1.0])) # Bin to the orbit grid
                            else
                                F_os_i = zeros(length(og.counts)) # Otherwise, zero
                            end
                        else
                            F_os_i = zeros(length(og.counts)) # Otherwise, zero
                        end
                        put!(channel,true) # Update the progress bar
                        F_os_i # Declare F_os_i as result to add to F_os
                    end
                    put!(channel,false) # Update progress bar
                    F_os # Delcare F_os as done/result, so it can be fetched
                end
            end)
        else
            result = @distributed (+) for i=1:numOsamples
                # Sample
                p = rand()*frdvols_cumsum_vector[end]
                j = searchsortedfirst(frdvols_cumsum_vector,p,Base.Order.Forward)
                inds = collect(Tuple(subs[j])) # First sample
                r = rand(4) .- 0.5 # 4 stands for the number of dimensions. 0.5 to sample within a hypercube
                E_sample = max(energy[inds[1]] + r[1]*dE_vector[inds[1]], 0.0)
                p_sample = pitch[inds[2]] + r[2]*dp_vector[inds[2]]
                R_sample = R[inds[3]] + r[3]*dR_vector[inds[3]]
                z_sample = z[inds[4]] + r[4]*dz_vector[inds[4]]

                # CHECK IF IT'S A GOOD SAMPLE
                good_sample = checkIfGoodSample(E_sample, p_sample, R_sample, z_sample, energy, pitch, R, z)
                if good_sample
                    o = get_orbit(M,GCP(E_sample,p_sample,R_sample,z_sample); store_path=false, wall=wall, kwargs...) # Calculate the orbit
                    if (o.coordinate.energy <= (maximum(og.energy)+dE_os_end/2) && o.coordinate.energy >= (minimum(og.energy)-dE_os_1/2)) # Make sure it's within the energy bounds (+/- one half grid cell)
                        F_os_i = bin_orbits(og,vec([o.coordinate]),weights=vec([1.0])) # Bin to the orbit grid
                    else
                        F_os_i = zeros(length(og.counts)) # Otherwise, zero
                    end
                else
                    F_os_i = zeros(length(og.counts)) # Otherwise, zero
                end
                F_os_i # Declare F_os_i as result to add to result
            end
        end

        return result
    else # ...otherwise just return a zero vector with the length of the number of valid orbits
        return zeros(length(og.counts))
    end
end

"""

checkIfGoodSample(E_sample, p_sample, R_sample, z_sample, energy, pitch, R, z)

The function checks if a sample is within bounds. Returns true if that is the case. Otherwise false.
"""
function checkIfGoodSample(E_sample::Float64, p_sample::Float64, R_sample::Float64, z_sample::Float64, energy::AbstractVector, pitch::AbstractVector, R::AbstractVector, z::AbstractVector)

    dp_1 = pitch[2]-pitch[1]
    dp_end = pitch[end]-pitch[end-1]
    dR_1 = R[2]-R[1]
    dR_end = R[end]-R[end-1]
    dz_1 = z[2]-z[1]
    dz_end = z[end]-z[end-1]
    if E_sample<0.0 || p_sample>(pitch[end]+dp_end/2) || p_sample<(pitch[1]-dp_1/2) || R_sample>(R[end]+dR_end/2) || R_sample<(R[1]-dR_1/2) || z_sample>(z[end]+dz_end/2) || z_sample<(z[1]-dz_1/2)
        return false
    else
        return true
    end
end

"""

unmap_orbits(og, F_os_3D)

This function unmaps the orbits from a 3D array into a standard compressed 1D array (with the elements representing the valid orbits).
It is essentially the inverse of the OrbitTomography.map_orbits() function. The function takes a function in orbit-space of the form (E,pm,Rm) and transforms it into a 1D array,
where each element corresponds to a valid orbit in orbit-space, with the indexing from og.orbit_index.
NOTE! The edges of the fast-ion distribution has to correspond to the edges of the orbit-space.
"""
function unmap_orbits(og::OrbitGrid, F_os_3D::AbstractArray; verbose::Bool=false)

    if !(length(og.energy)==size(F_os_3D,1) && length(og.pitch)==size(F_os_3D,2) && length(og.r)==size(F_os_3D,3))
        throw(ArgumentError("Dimensions of fast-ion distribution does not correspond to orbit-grid dimensions!"))
    end

    dorbs = get_orbel_volume(og,false) # Should always be false for safety and applicability
    if verbose
        println("Unmapping the orbits... ")
    end
    F_os = zeros(length(og.counts))
    for i=1:length(F_os_3D)
        if og.orbit_index[i]==0
        else
            F_os[og.orbit_index[i]] = F_os_3D[i]*dorbs[i]*og.counts[og.orbit_index[i]]
        end
    end

    return F_os
end

"""
    os2ps(M, Σ_ff_inv, F_os, og_orbs, corr_lengths, Js, energy, pitch, R, z)
    os2ps(-||-, distributed=true, atol=1e-4, domain_check=(xx,yy) -> false, covariance=:global, checkpoint=false, warmstart=true, file="test.jld2", verbose=true)

This function is complicated. It transforms a fast-ion distribution in orbit-space (E,pm,Rm) to particle-space (E,p,R,z). It requires knowledge
of the magnetic field and the tokamak geometry (M), the covariance between the orbits of the orbit grid (Σ_ff_inv), a vector that is essentially
the fast-ion distribution in orbit space, corresponding to the valid orbits only (F_os), all the orbits of the orbit grid (og_orbs)(each element
of that vector is of the orbit class, defined in the Julia package GuidingCenterOrbits.jl), the hyper-parameters correlation lengths (corr_lenghts)
(essentially setting a length in (E,p,R,z) by which the fast-ion distribution can be expected to correlate), the Jacobian for the orbit grid
(Js)(calculated with the get_jacobian() function, and the energy, pitch, R and z array determining the 4D grid to which you wish to transform the fast-ion distribution).

The optional keyword arguments are plenty. The most important ones are
distributed : wether you want to perform the calculations using multi-core computation or not (HIGHLY RECOMMENDED. Unless you prefer to wait halfway to forever for the result)
atol : the tolerance for the calculations of the covariance between the orbits of the orbit grid and the orbits of a specific (R,z) point
domain_check: I don't know. I don't even know why I include this one. Don't bother setting this to anything by the default value
covariance: Set to :local by default. Keep it that way, I have never tried setting it to :global and to tell you the truth, I haven't really understood it yet. It has to do with how the covariance is defined w.r.t. orbit space
checkpoint: If you want to periodically save the computational process, in case the computation is interrupted unexpectedly
warmstart: If true, then the function will load a half-finished computational process from the file 'file'
file: the file, in .jld2-format, in which the half-finished computational process is saved. Don't worry, it doesn't take up much space.
verbose: Do you want continuous print-out updates on how the computation is going? (No.. because you will be sleeping while this is running)

Don't bother with the other keyword arguments. I don't.
"""
function os2ps(M::AbstractEquilibrium, Σ_ff_inv::AbstractArray, F_os::Vector, og_orbs, corr_lengths::AbstractVector, Js::AbstractVector, energy::Array{Float64,1}, pitch::Array{Float64,1}, r::Array{Float64,1}, z::Array{Float64,1};
    distributed=false, atol=1e-3, domain_check= (xx,yy) -> true,
    covariance=:local, norms=OrbitTomography.S3(1.0,1.0,1.0),
    checkpoint=true, warmstart=false,file="eprz_progress.jld2",
    GCP=GCDeuteron, wall, verbose=false, kwargs...)

    nenergy = length(energy)
    npitch = length(pitch)
    nr = length(r)
    nz = length(z)
    inds = CartesianIndices((nr,nz)) # Define a set of cartesian coordinates for the number of R and z points. For simplicity
    m = og_orbs[1].coordinate.m # The mass of the orbital particles
    q = og_orbs[1].coordinate.q # The charge of the orbital particles

    f_eprz = zeros(nenergy,npitch,nr,nz) # Pre-allocate the resulting fast-ion distribution in particle-space f(E,p,R,z)

    if warmstart && isfile(file) && (filesize(file) != 0) # Load a half-finished computational process, if one exists
        progress_file = jldopen(file,false,false,false,IOStream)
        f_eprz = progress_file["f_eprz"]
        last_ind = progress_file["last_ind"]
        close(progress_file)
    else
        last_ind = inds[1] # Else, start from the beginning
    end

    if checkpoint
        touch(file) # Just mark the file for saving
    end

    for (iI,I) in enumerate(inds) # For all the (R,z) points...

        (I != inds[1] && I < last_ind) && continue # If we are still not done, continue
        i = I[1]
        j = I[2]
        rr = r[i]
        zz = z[j]
        !(domain_check(rr,zz)) && continue # I don't know. Ask L. Stagner

        #### Get all the orbits for this specific (R,z) point ####
        verbose && println("Getting orbits for loop $(iI) of $(length(inds))... ")
        lorbs = reshape([get_orbit(M, GCP(energy[k],pitch[l],rr,zz); wall=wall, kwargs...) for k=1:nenergy,l=1:npitch],nenergy*npitch)

        #### Calculate the Jacobian for those orbits ####
        verbose && println("Calculating Jacobian for loop $(iI) of $(length(inds))... ")
        if distributed # If multi-core computing
            lJs = pmap(o->get_jacobian(M,o), lorbs, on_error = ex->zeros(2))
            #batch_size=round(Int, nenergy*npitch/(5*nprocs()))) # Don't bother
        else # ...good luck!
            lJs = [get_jacobian(M, o) for o in lorbs]
        end

        verbose && println("Calculating orbit-space covariance for loop $(iI) of $(length(inds))... ")
        if covariance == :local # If it's local
            Si = get_covariance_matrix(M, lorbs, og_orbs, corr_lengths, Js_1=lJs, Js_2=Js,
                    distributed=distributed, atol=atol) # Get the covariance between the (R,z) orbit and the orbits of the orbit grid
        else
            Si = get_global_covariance_matrix(lorbs, og_orbs, corr_lengths, norms=norms)
        end

        f_ep = reshape(max.(Si*(Σ_ff_inv*F_os),0.0),nenergy,npitch) # Perform the operation to obtain f(E,p) for this (R,z) point
        detJ = reshape([length(j) != 0 ? j[1] : 0.0 for j in lJs],nenergy,npitch) # Needed to be able to...
        f_ep .= f_ep./detJ # ... normalize by the determinant of the Jacobian
        f_ep[detJ .== 0.0] .= 0.0 # But we don't want NaNs
        w = reshape([o.class in (:lost, :incomplete, :unknown) for o in lorbs],nenergy, npitch) # For all the lost, incomplete or unknown orbits...
        f_ep[w] .= 0.0 # Set their corresponding values to zero
        f_eprz[:,:,i,j] .= f_ep # Add your f(E,p) to the 4D f(E,p,R,z)
        if checkpoint # If you want to save
            progress_file = jldopen(file, true,true,true,IOStream)
            write(progress_file,"f_eprz",f_eprz) # Save
            write(progress_file,"last_ind",last_ind) # Save
            close(progress_file)

            #pvis_file = jldopen("$(round(count/length(inds),digits=3)).jld2", true,true,true,IOStream)
            #write(pvis_file,"count",count)
            #close(pvis_filei)
        end
        last_ind = I
    end

    return EPRZDensity(f_eprz,energy,pitch,r,z) # Return as EPRZDensity object
end

#"""
#Just a simple object to hold the fast-ion distribution in orbit-space, together with the corresponding energy, pitch-maximum
#and radius-maximum arrays, defining the orbit grid.
#"""
#struct EpmRmDensity{T::AbstractArray,S::AbstractVector}
#    d::T
#    energy::S
#    pm::S
#    Rm::S
#end

"""
    wahba_estimator(OS; kwargs...)

Compute the Wahba expression G using the Tikhonov regularisation parameter in OS.alpha. The Wahba expression G can be written as

                   |Wx-S|^2
G = ---------------------------------------
    Trace(I-W*Inv(Wt*W+alpha*alpha*I)*Wt)^2

Where I is the identity matrix, Wt is the transpose of W and alpha is OS.alpha. Inv(A) is the inverse of the matrix A. This function
is heavily inspired by the marginal_loglike() function written by L. Stagner, which can be found in the OrbitTomography.jl 
package in the tomography.jl file.

WARNING! If the signal S has a lot of diagnostic energy bins, then this algorithm will require a lot of RAM memory.
"""
function wahba_estimator(OS; norm=1.0e18/size(OS.W,2),max_iter=30*size(OS.W,1))
    S = OS.d
    noise = OS.err
    Σ_S = Diagonal(noise.^2)
    Σ_S_inv = inv(Σ_S)

    W_scaled = norm*OS.W

    F_prior_scaled = OS.mu/norm

    # Scale covariance matrices
    Σ_ff_scaled = OS.alpha*OS.S
    Σ_ff_inv_scaled = inv(OS.alpha)*OS.S_inv
    Γ_scaled = sqrt(inv(OS.alpha))*Matrix(OS.G)

    F_scaled = (try 
        vec(nonneg_lsq(vcat(W_scaled./noise, Γ_scaled), vcat(S./noise, Γ_scaled*F_prior_scaled); alg=:fnnls, use_parallel=false, max_iter=max_iter))
    catch er
        if isa(er,InterruptException)
            rethrow(er)
        end
        @warn "Non-negative Least Squares failed. Using MAP estimate"
        println(err)

        Σ_inv = Σ_ff_inv_scaled .+ W_scaled'*Σ_S_inv*W_scaled
        Σ_inv \ (W_scaled'*(Σ_S_inv*S) + Σ_ff_inv_scaled*F_prior_scaled)
    end)

    F = vec((norm/size(W_scaled,2)) .*F_scaled) # Re-scale solution to correct size

    W = OS.W
    alpha = OS.alpha

    nom = (LinearAlgebra.norm(W*F-S))^2
    denom = tr(Diagonal(ones(size(W,1))) - W*inv(transpose(W)*W + (alpha*alpha) .*Diagonal(ones(size(W,2))))*transpose(W))^2

    return nom/denom
end

"""
    wahba_optimize_tikhonov_factor!(OS; kwargs...)
    wahba_optimize_tikhonov_factor!(OS; log_bounds=(-3,3), kwargs...)

Wahba lahba dub dub! This optimisation function attempts to find the optimal Tikhonov regularisation factor by minimizing
the expression G, as found at https://en.wikipedia.org/wiki/Tikhonov_regularization under 'Determination of the Tikhonov factor'.
Grace Wahba proved that the optimal Tikhonov regularisation parameter, in the sense of leave-one-out cross-validation, minimizes the expression G.
See for example Wahba, G. (1990). "Spline Models for Observational Data". CBMS-NSF Regional Conference Series in Applied Mathematics. Society for 
Industrial and Applied Mathematics. See also for example www.doi.org/10.1080/00401706.1979.10489751 (Golub, G.; Heath, M.; Wahba, G. (1979). Technometrics. 21 (2)).

The function structure is heavily inspired by the optimize_alpha!() function written by L. Stagner, which can be found in the OrbitTomography.jl
package in the tomography.jl file.

The optimization will print info every 2nd iterations by default (show_trace=true and show_every=2, respectively).

A time limit for the optimization is set to 15 minutes (15*60 seconds) by default.
"""
function wahba_optimize_tikhonov_factor!(OS; log_bounds=(-6,6), show_trace=true, show_every=2, iterations=1000, kwargs...)
    f = x -> begin
        OS.alpha = 10.0^x
        return -wahba_estimator(OS; kwargs...)
    end

    op = optimize(f, log_bounds[1], log_bounds[2], Brent(),rel_tol=1e-3, show_trace=show_trace, show_every=show_every, iterations=iterations)

    OS.alpha = 10.0^Optim.minimizer(op)

    return Optim.minimum(op)
end

"""
    loglike_optimize_tikhonov_factor!(OS)

Given the OrbitSystem struct OS, compute the optimal tikhonov regularization factor by maximizing the log
likelihood.
"""
function loglike_optimize_tikhonov_factor!(OS; show_trace=true, show_every=2, iterations=1000, log_bounds = (-6,6), kwargs...)
    f = x -> begin
        OS.alpha = 10.0^x
        return -marginal_loglike(OS; kwargs...)
    end

    op = optimize(f, log_bounds[1], log_bounds[2], Brent(),rel_tol=1e-3, show_trace=show_trace, show_every=show_every, iterations=iterations)

    OS.alpha = 10.0^Optim.minimizer(op)

    return Optim.minimum(op)
end

"""
    getGCP(species_identifier)

The same function as is provided by OWCF/misc/species_func.jl. However, since it's needed in the os2COM() function, which is also
a part of extra/dependencies.jl, it needs to be defined here as well. Sometimes, it's more efficient to only load misc/species_func.jl, 
since one does not always want all the packages needed for dependencies.jl. Therefore, the functions below are also loadable as a separate
script, i.e. misc/species_func.jl.
"""
function getGCP(species_identifier::AbstractString)
    if lowercase(species_identifier)=="d"
        return GuidingCenterOrbits.GCDeuteron
    elseif lowercase(species_identifier)=="t"
        return GuidingCenterOrbits.GCTriton
    elseif lowercase(species_identifier)=="3he"
        return GuidingCenterOrbits.GCHelium3
    elseif lowercase(species_identifier)=="p"
        return GuidingCenterOrbits.GCProton
    elseif lowercase(species_identifier)=="e"
        return GuidingCenterOrbits.GCElectron
    elseif (lowercase(species_identifier)=="alpha" || lowercase(species_identifier)=="4he")
        return GuidingCenterOrbits.GCAlpha
    else
        error("getSpeciesMass got unknown species as input. Please re-try.")
    end
end

function getSpeciesMass(species_identifier::AbstractString)
    return (getGCP(species_identifier)(0.0,0.0,0.0,0.0)).m # 0.0, 0.0, 0.0, 0.0 is just to activate the function
end

function getSpeciesAmu(species_identifier::AbstractString)
    return getSpeciesMass(species_identifier) / GuidingCenterOrbits.mass_u
end

function getSpeciesEcu(species_identifier::AbstractString)
    return (getGCP(species_identifier)(0.0,0.0,0.0,0.0)).q # 0.0, 0.0, 0.0, 0.0 is just to activate the function
end

function getSpeciesCharge(species_identifier::AbstractString)
    return getSpeciesEcu(species_identifier) * GuidingCenterOrbits.e0
end

"""
    getOSTopoMap(M, E, pmRm_inds, pm_array, Rm_array)
    getOSTopoMap(-||-; FI_species="D", wall=nothing, distinguishLost=false, distinguishIncomplete=false, extra_kw_args=Dict(:toa => true))

Compute the orbit-sapce topological map for the energy slice specified by E (keV). Use pmRm_inds to iterate through the pm and Rm coordinates. 
"""
function getOSTopoMap(M::AbstractEquilibrium, E::Float64, pmRm_inds_array::Vector{Tuple{CartesianIndex{1}, CartesianIndex{1}}}, pm_array::AbstractVector, Rm_array::AbstractVector; FI_species::String="D", wall::Union{Nothing,Boundary{Float64}}=nothing, distinguishLost::Bool=false, distinguishIncomplete::Bool=false, extra_kw_args=Dict(:toa => true))
    topoMap = @showprogress 1 "Computing topoMap for E=$(E) keV... " @distributed (+) for i=1:length(pmRm_inds_array)
        pmRm_inds = pmRm_inds_array[i]
        ipm = pmRm_inds[1]
        iRm = pmRm_inds[2]
        pm = pm_array[ipm]
        Rm = Rm_array[iRm]
        topoMap_i = zeros(length(pm_array),length(Rm_array))
        EPRc = EPRCoordinate(M,E,pm,Rm;amu=getSpeciesAmu(FI_species), q=getSpeciesEcu(FI_species))
        o = get_orbit(M,EPRc;wall=wall,store_path=false, extra_kw_args...)
        topoMap_i[ipm,iRm] = class2int(M, o; plot=true, distinguishLost=distinguishLost, distinguishIncomplete=distinguishIncomplete) # Plot=true ensures correct integer for c-stag orbits
        topoMap_i
    end

    return topoMap
end

"""
    getOSTopoMap(M, E_array, pm_array, Rm_array)
    getOSTopoMap(; kwargs... )

Compute the orbit-space topological map, given a grid spanned by the inputs E_array, pm_array and Rm_array.
"""
function getOSTopoMap(M::AbstractEquilibrium, E_array::AbstractVector, pm_array::AbstractVector, Rm_array::AbstractVector; kwargs... )
    pm_inds_rep = repeat(CartesianIndices(pm_array),inner=length(Rm_array)) # To get all points
    Rm_inds_rep = repeat(CartesianIndices(Rm_array),outer=length(pm_array)) # To get all points
    topoMap = zeros(length(E_array),length(pm_array),length(Rm_array)) # The total 3D topological map
    pmRm_inds_array = collect(zip(pm_inds_rep,Rm_inds_rep)) # All points in an energy slice, in one long vector
    for iE=1:length(E_array)
        topoMap[iE,:,:] = getOSTopoMap(M, E_array[iE], pmRm_inds_array, pm_array, Rm_array; kwargs... )
    end
    return topoMap
end

"""
    interpDelaunayTess(tess, data, x_array, y_array, vertices_hash)
    interpDelaunayTess(-||-; outside_value=0.0, nearest=false, verbose=false, vverbose=false)

Takes the points in 'x_array' and 'y_array', and find the Barycentrically interpolated value at (x,y) given the Delaunay tesselation 'tess' and the
'data'. Use the dictionary 'vertices_hash' to enable the Barycentric interpolation; relate the vertices of a triangle with the (x,y) point inside of it,
to the data value at the vertices. The 'vertices_hash' is an output of the getDelaunayTessVerticesHash() function, defined below.
"""
function interpDelaunayTess(tess::DelaunayTessellation2D{Point2D}, data::AbstractVector, x_array::AbstractVector, y_array::AbstractVector, vertices_hash::Dict; outside_value::Number=0.0, nearest::Bool=false, verbose::Bool=false, vverbose::Bool=false)
    data_interp = Array{Float64}(undef,length(x_array),length(y_array))
    verbose && println("Interpolating data onto all (x,y) query points... ")
    count = 1
    for (ix,xx) in enumerate(x_array), (iy,yy) in enumerate(y_array)
        vverbose && println("$(count)/$(length(x_array)*length(y_array)): (x,y)=($(xx),$(yy))")
        t = locate(tess, Point(xx, yy)) # Try to locate the query point inside of the tessellation

        if isexternal(t) == true # If the query point was outside of the tessellation...
            # HERE, THERE COULD BE SOME FANCY NEAREST NEIGHBOUR SCHEME THAT PRESERVES E.G. STAGNATION ORBITS BETTER
            data_interp[ix,iy] = outside_value
        else
            tID = "$(t._neighbour_a)_$(t._neighbour_b)_$(t._neighbour_c)" # Contruct the triangle ID from the order of the neighbours
            ia, ib, ic = vertices_hash[tID] # Look up the original indices of the vertices
            p = [xx, yy] # Query point in tri coordinates
            pa = [getx(geta(t)), gety(geta(t))] # Vertex a of triangle t
            pb = [getx(getb(t)), gety(getb(t))] # Vertex b of triangle t
            pc = [getx(getc(t)), gety(getc(t))] # Vertex c of triangle t

            μ = [(1/sum(abs2,p-pa)),(1/sum(abs2,p-pb)),(1/sum(abs2,p-pc))] 
            μ = (1/sum(μ)) .* μ # Barycentric weights. sum(μ)=1 must hold

            # Barycentric interpolation. All components of μ sum to one. 'Barycentric' simply means 'center-of-mass-like'. https://en.wikipedia.org/wiki/Barycentric_coordinate_system 
            vertices_data = [data[ia], data[ib], data[ic]] 
            interp_value = μ[1]*vertices_data[1]+μ[2]*vertices_data[2]+μ[3]*vertices_data[3]
            if nearest
                interp_value = vertices_data[closest_index(vertices_data,interp_value)] # closest_index() function found near top of dependencies.jl
            end
            data_interp[ix,iy] = interp_value
        end
        count += 1
    end
    return data_interp
end

"""
    getDelaunayTessVerticesHash(tess, iterator)
    getDelaunayTessVerticesHash(tess, iterator; verbose=true)

Take the Delaunay tesselation 'tess' and relate all points in 'iterator' to the triangle ID
of every triangle in the tesselation. The points in 'iterator' are the points at the vertices of the
Delaunay triangles making up the tesselation. The triangle ID should be unique for every triangle in the
tesselation. Thus, if we know the triangle, we can find the original indices in 'iterator' of the tree
vertices of the triangle, by using the output hashmap 'vertices'. 
"""
function getDelaunayTessVerticesHash(tess::DelaunayTessellation2D{Point2D}, iterator; verbose::Bool=false)
    vertices = Dict{String,Tuple{Int64,Int64,Int64}}() # A hashmap. The hashes are triangle IDs and the values are tuples with vertices indices in the same order as the a,b and c fields of the DelaunayTriangle object returned by locate(tess, point). The DelaunayTesselation struct provided no way to keep track of what indices of the (x,y) inputs that make up the vertices of the triangles in the tessellation... So I had to do it myself!
    count = 1
    for t in tess # For triangles in the tessellation
        verbose && print("Triangle $(count).     ")
        ind_a = findall(x-> x==(getx(geta(t)),gety(geta(t))), collect(iterator)) # Find the (x,y) data index of the first vertex to triangle t
        ind_b = findall(x-> x==(getx(getb(t)),gety(getb(t))), collect(iterator)) # Find the (x,y) data index of the second vertex to triangle t
        ind_c = findall(x-> x==(getx(getc(t)),gety(getc(t))), collect(iterator)) # Find the (x,y) data index of the third vertex to triangle t

        verbose && print("ind_a: $(ind_a)    ind_b: $(ind_b)     ind_c: $(ind_c).     ")

        triangle_ID = "$(t._neighbour_a)_$(t._neighbour_b)_$(t._neighbour_c)" # No other triangle in the tessellation will have this ID. Thus making it easily look-up-able later
        vertices[triangle_ID] = (first(ind_a), first(ind_b), first(ind_c)) # Essentially, we are mapping tesselation triangle ID to input data index ID. This is to be able to interpolate new values at a query point.
        verbose && println("Triangle ID: "*triangle_ID)
        count += 1
    end
    verbose && println()
    return vertices
end

"""
    pmRm_2_μPϕ(M, good_coords_pmRm, data, E, pm_array, Rm_array, FI_species)
    pmRm_2_μPϕ(-||- ; nμ=length(pm_array), nPϕ=length(Rm_array), isTopoMap=false, verbose=false)

Take the 2D array 'data' and map its elements from (pm,Rm) to (μ,Pϕ) for the energy E (given in keV) and ion species FI_species. 
The output will be a 3D array where the last dimension corresponds to the binary σ coordinate (+1 or -1) that keeps track of co- and 
counter-going orbits. Co-going orbits (co-passing, trapped, potato and stagnation) will have σ=+1. Counter-going orbits (counter-passing
and counter-stagnation) will have σ=-1. The good_coords_pmRm vector contains coordinate indices for the mappable (pm,Rm) coordinates 
(to avoid mapping invalid orbits).
"""
function pmRm_2_μPϕ(M::AbstractEquilibrium, good_coords_pmRm::Vector{CartesianIndex{2}}, data::Array{Float64,2}, E::Union{Float64,Int64}, pm_array::AbstractVector, Rm_array::AbstractVector, FI_species::AbstractString; nμ::Int64=length(pm_array), nPϕ::Int64=length(Rm_array), isTopoMap::Bool=false, needJac::Bool=false, transform = x -> x, verbose::Bool=false, vverbose::Bool=false, debug::Bool=false)

    verbose && println("Transforming (pm,Rm) coordinates into (μ,Pϕ) coordinates... ")
    μ_values = Array{Float64}(undef, length(good_coords_pmRm))
    Pϕ_values = Array{Float64}(undef, length(good_coords_pmRm))
    data_values = Array{Float64}(undef, length(good_coords_pmRm))
    pm_signs = Array{Float64}(undef, length(good_coords_pmRm))
    if needJac # If a Jacobian is needed for the (pm,Rm) -> (μ,Pϕ) mapping (i.e. for distributions)...
        # We need to convert our (E,pm,Rm) coordinates to Dual numbers (a+ϵb), where the dual part (b) is a unique orthonormal vector for each E-, pm- and Rm-direction
        verbose && println("Converting pmRm-, and pre-allocating μPϕ-, arrays as Dual-number arrays... ")
        E = ForwardDiff.Dual(E,(1.0,0.0,0.0))
        pm_array = map(x-> ForwardDiff.Dual(x,(0.0,1.0,0.0)), pm_array)
        Rm_array = map(x-> ForwardDiff.Dual(x,(0.0,0.0,1.0)), Rm_array)
        μ_values = Array{ForwardDiff.Dual{Nothing, Float64, 3}}(undef, length(good_coords_pmRm))
        Pϕ_values = Array{ForwardDiff.Dual{Nothing, Float64, 3}}(undef, length(good_coords_pmRm))
        pm_signs = Array{ForwardDiff.Dual{Nothing, Float64, 3}}(undef, length(good_coords_pmRm))
    end
    for (i_good_coord, good_coord) in enumerate(good_coords_pmRm)
        vverbose && println("Transforming (pm,Rm) coordinate $(i_good_coord) of $(length(good_coords_pmRm))... ")
        ipm, iRm = Tuple(good_coord)
        myEPRc = EPRCoordinate(M, E, pm_array[ipm], Rm_array[iRm]; amu=getSpeciesAmu(FI_species), q=getSpeciesEcu(FI_species))
        myHc = HamiltonianCoordinate(M, myEPRc)

        μ_values[i_good_coord] = myHc.mu
        Pϕ_values[i_good_coord] = myHc.p_phi
        data_values[i_good_coord] = data[ipm,iRm]
        pm_signs[i_good_coord] = sign(pm_array[ipm])
    end

    if needJac # Could, but should not, do this in the loop above. This way, we achieve optimized speed for pmRm_2_μPϕ() when needJac is false
        verbose && print("Scaling data with 1/|J|... ")
        for i=1:length(data_values)
            x = transform([E,μ_values[i],Pϕ_values[i]]) # Probably redundant
            detJac = max(abs(det(hcat((ForwardDiff.partials(xx) for xx in x)...))),0.0) # Extract the determinant of the Jacobi matrix (dEdμdPϕ/dEdpmdRm)
            data_values[i] = (1/detJac) * data_values[i] # Multiply by the inverse of the Jacobian to go from (E,pm,Rm) to (E,μ,Pϕ;σ). Please note! σ is automatically accounted for by the dual numbers
        end
        verbose && println("Success!")
        verbose && println("Reducing Dual-number arrays to real arrays... ")
        E = ForwardDiff.value(E) # Extract the real part of E (remove the dual part)
        pm_array = map(x-> ForwardDiff.value(x),pm_array) # Extract the real part of pm (remove the dual part)
        Rm_array = map(x-> ForwardDiff.value(x),Rm_array) # Extract the real part of Rm (remove the dual part)
        μ_values = map(x-> ForwardDiff.value(x),μ_values) # Extract the real part of μ (remove the dual part)
        Pϕ_values = map(x-> ForwardDiff.value(x),Pϕ_values) # Extract the real part of Pϕ (remove the dual part)
        pm_signs = map(x-> ForwardDiff.value(x),pm_signs) # Extract the real part of sign(pm) (remove the dual part)
    end

    # We will let σ=+1 correspond to all types of co-going orbits.
    # We will therefore need to treat them separately.
    # Orbits with pm>=0 are co-going. Let's start with them. 
    # 'cogo' is short for 'co-going'
    cogo_inds = findall(x-> x>=0.0, pm_signs)
    if isTopoMap
        lost_inds = findall(x-> x==7.0, data_values)
        lost_cogo_inds = filter(x-> x in lost_inds, cogo_inds)
        cogo_inds = filter(x-> !(x in lost_inds), cogo_inds) # Remove all lost orbits. They need to be treated separately
    end
    data_values_cogo = data_values[cogo_inds]
    μ_values_cogo = μ_values[cogo_inds]
    Pϕ_values_cogo = Pϕ_values[cogo_inds]

    # Since the μ- and Pϕ-values will be irregularly scattered in the (μ,Pϕ)-plane
    # (because we simply just computed them from (pm,Rm) points), we need to create 
    # a Delaunay tesselation to be able to interpolate onto our desired (nμ,nPϕ) 
    # output grid (which will be regularly spaced)
    tess = DelaunayTessellation(length(μ_values_cogo))

    # The Delaunay tesselation unfortunately only works if the nodes are
    # within (1+eps, 2-eps). We therefore need to map all our (μ,Pϕ)-points
    # into the box with vertices (1,1), (1,2), (2,1), (2,2).
    # min_coord = 1 + eps is automatically loaded from the DelaunayVoronoi.jl package
    # max_coord = 2 - eps is automatically loaded from the DelaunayVoronoi.jl package
    min_Pϕ, max_Pϕ = minimum(Pϕ_values), maximum(Pϕ_values) # Let's not use the 'cogo' arrays, but the total arrays. To enable the min/max values to be used for the counter-going orbits later as well
    min_μ, max_μ = minimum(μ_values), maximum(μ_values)
    Pϕ_tess(Pϕ) = min_coord + (max_coord-min_coord)*(Pϕ-min_Pϕ)/(max_Pϕ-min_Pϕ) # DelaunayTesselation expects all points to be within [min_coord,max_coord] (given by VoronoiDelaunay.min_coord etc). We therefore need a conversion function.
    Pϕ_tess_inv(Pϕ_tess) = min_Pϕ + (max_Pϕ-min_Pϕ)*(Pϕ_tess-min_coord)/(max_coord-min_coord) # We need a function to convert back to normal Pϕ values
    μ_tess(μ) = min_coord + (max_coord-min_coord)*(μ-min_μ)/(max_μ-min_μ) # DelaunayTesselation expects all points to be within [min_coord,max_coord] (given by VoronoiDelaunay.min_coord etc). We therefore need a conversion function.
    μ_tess_inv(μ_tess) = min_μ + (max_μ-min_μ)*(μ_tess-min_coord)/(max_coord-min_coord) # We need a function to convert back to normal μ values

    μPϕ_iterator_tess = zip(μ_tess.(μ_values_cogo),Pϕ_tess.(Pϕ_values_cogo)) # Create an iterator with all the co-going (μ,Pϕ) values as tuples
    a = Point2D[Point(mm_tess, pp_tess) for (mm_tess,pp_tess) in μPϕ_iterator_tess] # Put them all in a Point() object within a Point2D array
    push!(tess, a) # Feed all points to the Delaunay tessellation

    verbose && println("Mapping all co-going triangles to original vertex indices... ")
    vertices_hash = getDelaunayTessVerticesHash(tess, μPϕ_iterator_tess)

    μ_array = collect(range(min_μ,stop=max_μ,length=nμ)) # The μ grid points onto which to interpolate
    Pϕ_array = collect(range(min_Pϕ,stop=max_Pϕ,length=nPϕ)) # The Pϕ grid points onto which to interpolate
    if isTopoMap
        outside_value = 9.0 # If we are mapping a topological map and the query point is outside of the tesselation, it's an invalid orbit
    else
        outside_value = 0.0 # If we are mapping anything else (orbit weight function slice, transit time map etc.), all values should be zero outside of the tesselation
    end
    data_COM_cogo = interpDelaunayTess(tess, data_values_cogo, μ_tess.(μ_array), Pϕ_tess.(Pϕ_array), vertices_hash; outside_value=outside_value, nearest=isTopoMap, verbose=verbose, vverbose=vverbose)

    if isTopoMap && !(length(lost_inds)==0) # If we are transforming a topological map, and there are lost orbits included
        # We need to do the same for the lost (co-going) orbits as we just did for the co-going orbits
        vverbose && println("Mapping all co-going lost orbits to (μ, Pϕ)... ")
        data_values_lost_cogo = data_values[lost_cogo_inds]
        μ_values_lost_cogo = μ_values[lost_cogo_inds]
        Pϕ_values_lost_cogo = Pϕ_values[lost_cogo_inds]
        tess = DelaunayTessellation(length(μ_values_lost_cogo))
        μPϕ_iterator_tess = zip(μ_tess.(μ_values_lost_cogo),Pϕ_tess.(Pϕ_values_lost_cogo)) # Create an iterator with all the co-going (μ,Pϕ) values as tuples
        a = Point2D[Point(mm_tess, pp_tess) for (mm_tess,pp_tess) in μPϕ_iterator_tess] # Put them all in a Point() object within a Point2D array
        push!(tess, a) # Feed all points to the Delaunay tessellation
        vertices_hash = getDelaunayTessVerticesHash(tess, μPϕ_iterator_tess)
        data_COM_lost_cogo = interpDelaunayTess(tess, data_values_lost_cogo, μ_tess.(μ_array), Pϕ_tess.(Pϕ_array), vertices_hash; outside_value=outside_value, nearest=isTopoMap, verbose=verbose, vverbose=vverbose)
        data_COM_cogo = [((data_COM_cogo[sub] != 9.0) ? data_COM_cogo[sub] : data_COM_lost_cogo[sub]) for sub in CartesianIndices(data_COM_cogo)]
    end

    # Now we need to do the same for the counter-going orbits
    # We will let σ=-1 correspond to all types of counter-going orbits.
    # We will therefore need to treat them separately (just as for co-going) 
    # Orbits with pm<0 are counter-going. Let's continue with them now.
    # 'ctgo' is short for 'counter-going'
    ctgo_inds = findall(x-> x<0.0, pm_signs)
    if isTopoMap
        lost_ctgo_inds = filter(x-> x in lost_inds, ctgo_inds)
        ctgo_inds = filter(x-> !(x in lost_inds), ctgo_inds) # Remove all lost orbits. They need to be treated separately
    end
    data_values_ctgo = data_values[ctgo_inds]
    μ_values_ctgo = μ_values[ctgo_inds]
    Pϕ_values_ctgo = Pϕ_values[ctgo_inds]
    tess = DelaunayTessellation(length(μ_values_ctgo))
    μPϕ_iterator_tess = zip(μ_tess.(μ_values_ctgo),Pϕ_tess.(Pϕ_values_ctgo)) # Create an iterator with all the counter-going (μ,Pϕ) values as tuples
    a = Point2D[Point(mm_tess, pp_tess) for (mm_tess,pp_tess) in μPϕ_iterator_tess] # Put them all in a Point() object within a Point2D array
    push!(tess, a) # Feed all points to the Delaunay tessellation
    verbose && println("Mapping all counter-going triangles to original vertex indices... ")
    vertices_hash = getDelaunayTessVerticesHash(tess, μPϕ_iterator_tess)
    data_COM_ctgo = interpDelaunayTess(tess, data_values_ctgo, μ_tess.(μ_array), Pϕ_tess.(Pϕ_array), vertices_hash; outside_value=outside_value, nearest=isTopoMap, verbose=verbose, vverbose=vverbose)
    if isTopoMap && !(length(lost_inds)==0) # If we are transforming a topological map, and there are lost orbits included
        # We need to do the same for the lost (counter-going) orbits as we just did for the counter-going orbits
        vverbose && println("Mapping all counter-going lost orbits to (μ, Pϕ)... ")
        data_values_lost_ctgo = data_values[lost_ctgo_inds]
        μ_values_lost_ctgo = μ_values[lost_ctgo_inds]
        Pϕ_values_lost_ctgo = Pϕ_values[lost_ctgo_inds]
        tess = DelaunayTessellation(length(μ_values_lost_ctgo))
        μPϕ_iterator_tess = zip(μ_tess.(μ_values_lost_ctgo),Pϕ_tess.(Pϕ_values_lost_ctgo)) # Create an iterator with all the counter-going (μ,Pϕ) values as tuples
        a = Point2D[Point(mm_tess, pp_tess) for (mm_tess,pp_tess) in μPϕ_iterator_tess] # Put them all in a Point() object within a Point2D array
        push!(tess, a) # Feed all points to the Delaunay tessellation
        if debug
            # Write your own debug code here. For example:
            x, y = getplotxy(delaunayedges(tess))
            test_plt_1 = Plots.plot(Pϕ_tess_inv.(y),μ_tess_inv.(x); label="Delaunay tessellation",ylims=extrema(μ_tess_inv.(x)), xlims=extrema(Pϕ_tess_inv.(y)))
            display(test_plt_1)
        end
        vertices_hash = getDelaunayTessVerticesHash(tess, μPϕ_iterator_tess)
        data_COM_lost_ctgo = interpDelaunayTess(tess, data_values_lost_ctgo, μ_tess.(μ_array), Pϕ_tess.(Pϕ_array), vertices_hash; outside_value=outside_value, nearest=isTopoMap, verbose=verbose, vverbose=vverbose)
        data_COM_ctgo = [((data_COM_ctgo[sub] != 9.0) ? data_COM_ctgo[sub] : data_COM_lost_ctgo[sub]) for sub in CartesianIndices(data_COM_ctgo)]
    end
    data_COM = Array{Float64}(undef,nμ,nPϕ,2)
    #println(unique(data_COM_ctgo))
    data_COM[:,:,1] = data_COM_ctgo # Counter-going (first element, σ=-1)
    data_COM[:,:,2] = data_COM_cogo # Co-going (Second element, σ=+1)
    if debug
        return μ_values_cogo, μ_values_ctgo, Pϕ_values_cogo, Pϕ_values_ctgo, μ_array, Pϕ_array
    end

    return map(x-> isnan(x) ? 0.0 : x, data_COM), E, μ_array, Pϕ_array # Make sure to remove all NaNs. Sometimes, NaNs arise. It's inevitable.
end


"""
    os2COM(M, good_coords, data, E_array, pm_array, Rm_array, FI_species)
    os2COM(M, good_coords, data, E_array, pm_array, Rm_array, FI_species; nμ=length(pm_array), nPϕ=length(Rm_array), isTopoMap=false, verbose=false)

    This function maps an orbit-space (E,pm,Rm) quantity into COM-space (E, μ, Pϕ; σ). It assumes the quantity is given as 3D data and that 
    the (E,pm,Rm) coordinates suitable for mapping are given as good_coords, a Vector{CartesianIndex{3}}. The output data will be 4D. The last 
    dimension will correspond to the binary σ coordinate (+1 or -1) that keeps track of co- and counter-passing orbits with the same (E, μ, Pϕ) coordinate.
    
    PLEASE NOTE! Every 3D quantity corresponding to each index in the energy dimension of the 4D output (i.e. data_COM[1,:,:,:], data_COM[2,:,:,:] etc) has its own pair of μ- and 
    Pϕ-arrays (i.e. μ_matrix[1,:]/Pϕ_matrix[1,:], μ_matrix[2,:]/Pϕ_matrix[2,:] etc). This is because the minimum and maximum of μ and Pϕ scales with the energy.
"""
function os2COM(M::AbstractEquilibrium, good_coords::Vector{CartesianIndex{3}}, data::Array{Float64, 3}, E_array::AbstractVector, pm_array::AbstractVector, Rm_array::AbstractVector, FI_species::AbstractString; nμ::Int64=length(pm_array), nPϕ::Int64=length(Rm_array), verbose::Bool=false, kwargs...)
    
    data_COM = zeros(length(E_array),nμ, nPϕ, 2)
    μ_matrix = zeros(length(E_array),nμ) # Create matrix, because we need to save all possible μ-values for all possible energies (the min/max values of μ scale with the energy)
    Pϕ_matrix = zeros(length(E_array),nPϕ) # Create matrix, because we need to save all possible Pϕ-values for all possible energies (the min/max values of Pϕ scale with the energy)
    for (iE, E) in enumerate(E_array)
        verbose && println("Mapping energy slice $(iE) of $(length(E_array))... ")
        good_coords_iE = findall(x-> x[1]==iE, good_coords) # Returns a 1D vector with Integer elements
        good_coords_pmRm = map(x-> CartesianIndex(x[2],x[3]), good_coords[good_coords_iE]) # Returns a vector with CartesianIndex{2} elements
        data_COM[iE,:,:,:], E, μ_array, Pϕ_array = pmRm_2_μPϕ(M, good_coords_pmRm, data[iE,:,:], E, pm_array, Rm_array, FI_species; nμ=nμ, nPϕ=nPϕ, verbose=verbose, kwargs...)
        μ_matrix[iE,:] = μ_array # Save the array of possible μ-values FOR THE CURRENT ENERGY E in a matrix at row iE
        Pϕ_matrix[iE,:] = Pϕ_array # Save the array of possible Pϕ-values FOR THE CURRENT ENERGY E in a matrix at row iE
    end
    
    return data_COM, E_array, μ_matrix, Pϕ_matrix
end

"""
    os2COM(M, data, E_array, pm_array, Rm_array, FI_species)
    os2COM(M, data, E_array, pm_array, Rm_array, FI_species; nμ=length(pm_array), nPϕ=length(Rm_array), isTopoMap=false, verbose=false)

This function maps an orbit-space (E,pm,Rm) quantity into COM-space (E, μ, Pϕ; σ). If a 4D quantity is given as input data, the function will assume it's a collection of 3D (E,pm,Rm) quantities with
the last three dimensions corresponding to E, pm and Rm. Return transformed data in COM, along with E-array, μ-matrix and Pϕ-matrix. The output data will be either 4D or 5D. The last dimension will 
correspond to the binary σ coordinate (-1 or +1) that keeps track of counter- and co-going orbits with the same (E, μ, Pϕ) coordinate. The first Julia index (1)
corresponds to σ=-1 and the second Julia index (2) corresponds to σ=+1.

If the input data is not a topological map, then we need to identify the (E,pm,Rm) coordinates suitable for mapping by first computing the corresponding orbit grid.

PLEASE NOTE! Every 3D quantity corresponding to each index in the energy dimension of the 5D output (i.e. data_COM[iE,1,:,:,:], data_COM[iE,2,:,:,:] etc) has its own pair of μ- and 
Pϕ-arrays (i.e. μ_matrix[1,:]/Pϕ_matrix[1,:], μ_matrix[2,:]/Pϕ_matrix[2,:] etc). This is because the minimum and maximum of μ and Pϕ scales with the energy, so each energy slice needs
its own pair of μ and Pϕ grid points (given the input grid points in pm and Rm). You could, of course, interpolate all slices onto the same (μ,Pϕ) grid. BUT, since e.g. maximum(μ) scales
by a factor of 100 when you go from a 10 keV energy slice to a 1 MeV energy slice, it's not really practical. Then again, you could of course interpolate onto an energy normalized (μ,Pϕ) grid.
But I thought that maybe you, dear user, would like to do that yourself post-computation. I am giving you the freedom of choice here. That's why the (μ,Pϕ) output is given in matrix form, NOT
in vector form.
"""
function os2COM(M::AbstractEquilibrium, data::Union{Array{Float64, 3},Array{Float64, 4}}, E_array::AbstractVector, pm_array::AbstractVector, Rm_array::AbstractVector, FI_species::AbstractString; nμ::Int64=length(pm_array), nPϕ::Int64=length(Rm_array), isTopoMap::Bool=false, needJac::Bool=false, verbose::Bool=false, vverbose::Bool=false, good_coords=nothing, wall::Union{Nothing,Boundary{Float64}}=nothing, extra_kw_args=Dict(:toa => true, :limit_phi => true, :maxiter => 0))
    if !isTopoMap && !(typeof(good_coords)==Vector{CartesianIndex{3}})
        verbose && println("Input data is not a topological map, and keyword argument 'good_coords' has not been (correctly?) provided... ")
        verbose && println("--------> Computing orbit grid for (E,pm,Rm)-values to be able to deduce valid orbits for (E,pm,Rm)-grid... ")
        orbs, og = orbit_grid(M, verbose, E_array, pm_array, Rm_array; q=getSpeciesEcu(FI_species), amu=getSpeciesAmu(FI_species), wall=wall, extra_kw_args...)
    end

    if isTopoMap
        verbose && println("Identifying mappable (E,pm,Rm) coordinates from topological map... ")
        good_coords = findall(x-> x!=9.0,data) # If isTopoMap, then we know it's a 3D quantity. So good_coords will be a vector consisting of CartesianIndices{3}
    elseif (typeof(good_coords)==Vector{CartesianIndex{3}})
        verbose && println("Mappable (E,pm,Rm) coordinates provided via 'good_coords' arguments. os2COM will attempt to use them... ")
        good_coords = good_coords
    else
        verbose && print("Identifying mappable (E,pm,Rm) coordinates from orbit grid... ")
        good_coords = findall(x-> x>0.0, og.orbit_index) # Valid orbits (orbit_index > 0). So good_coords will be a vector consisting of CartesianIndices{3}
        verbose && println("($(length(good_coords)))")
    end

    input_4D = true
    if !(length(size(data))==4)
        verbose && println("Input data not 4D. Reshaping... ")
        input_4D = false
        data = reshape(data,(1,size(data)...)) # Reshape into general 4D shape. Sometimes, data might be 4D (many 3D OWs) and it's therefore best to always have the code handle 4D data.
    end

    data_COM = zeros(size(data,1),length(E_array),nμ, nPϕ, 2)
    μ_matrix = nothing # Just pre-define something
    Pϕ_matrix = nothing
    verbose && println("Mapping (E,pm,Rm) -> (E,μ,Pϕ;σ) for 3D quantity 1 of $(size(data,1))... ")
    data_COM[1, :, :, :, :], E_array, μ_matrix, Pϕ_matrix = os2COM(M, good_coords, data[1,:,:,:], E_array, pm_array, Rm_array, FI_species; nμ=nμ, nPϕ=nPϕ, isTopoMap=isTopoMap, needJac=needJac, verbose=verbose, vverbose=vverbose)
    if size(data,1)>1 # To avoid distributed errors
        data_COM_mp = @showprogress 1 "Mapping (E,pm,Rm) -> (E,μ,Pϕ;σ)... " @distributed (+) for iEd in collect(2:size(data,1))
            #verbose && println("Mapping (E,pm,Rm) -> (E,μ,Pϕ;σ) for 3D quantity $(iEd) of $(size(data,1))... ")
            data_COM_part = zeros(size(data,1),length(E_array),nμ, nPϕ, 2)
            # E_array_part, μ_matrix_part and Pϕ_matrix_part do not matter. But they need to be returned by os2COM.
            data_COM_part[iEd, :, :, :, :], E_array_part, μ_matrix_part, Pϕ_matrix_part = os2COM(M, good_coords, data[iEd,:,:,:], E_array, pm_array, Rm_array, FI_species; nμ=nμ, nPϕ=nPϕ, isTopoMap=isTopoMap, needJac=needJac, verbose=false, vverbose=false)
            data_COM_part # Declare ready for reduction (data_COM += data_COM_part but distributed over several CPU cores)
        end

        data_COM = data_COM + data_COM_mp # Maybe inefficient/non-optimal coding, but I don't care. It works!
    end

    if !input_4D
        return dropdims(data_COM,dims=1), E_array, μ_matrix, Pϕ_matrix 
    end
    return data_COM, E_array, μ_matrix, Pϕ_matrix
end


"""
    muPphi_2_pmRmNorm_funcGen(M, E, FI_species)
    muPphi_2_pmRmNorm_funcGen(-||-; amu=getSpeciesAmu(FI_species), q=getSpeciesEcu(FI_species), sigma=1, wall=nothing, nR=500, verbose=false, extra_kw_args=Dict(:max_tries => 0))

Return a function that returns sqrt(pm^2 + Rm^2) given a (mu,Pphi) coordinate as input. This function(al) is deprecated, but might be brought back in a future version of the OWCF.
"""
function muPphi_2_pmRmNorm_funcGen(M::AbstractEquilibrium, E::Float64, FI_species::String; amu::Number=getSpeciesAmu(FI_species), q::Number=getSpeciesEcu(FI_species), sigma::Int64=1, wall::Union{Nothing,Boundary{Float64}}=nothing, nR::Int64=500, verbose::Bool=false, extra_kw_args=Dict(:max_tries => 0))
    f = function func(x)
        myEPRc = EPRCoordinate4(M, HamiltonianCoordinate(E, x[1], x[2]; amu=amu, q=q); sigma=sigma, wall=wall, nR=nR, verbose=verbose, extra_kw_args=extra_kw_args)
        return (myEPRc.r > 0.0) ? sqrt((myEPRc.pitch)^2 + (myEPRc.r)^2) : Inf
    end
    return f
end

"""
    muPphiSigma_2_pmRm(...)
    muPphiSigma_2_pmRm(...)

Compute a (pm,Rm) coordinate, given a (mu,Pphi;sigma) coordinate as input. This function is deprecated, but might be brought back in a future version of the OWCF.
"""
function muPphiSigma_2_pmRm(M::AbstractEquilibrium, data::Array{Float64,3}, E::Union{Float64, Int64}, mu_array::AbstractVector, 
    Pphi_array::AbstractVector, FI_species::AbstractString, npm::Int64, nRm::Int64, wall::Union{GuidingCenterOrbits.Boundary{Float64},Nothing}, 
    pm_array::Union{AbstractVector,Nothing}, Rm_array::Union{AbstractVector,Nothing}, dataIsTopoMap::Bool, topoMap::Array{Float64,3}, needJac::Bool, 
    nR::Int64, n_stagopt::Int64, verbose::Bool, vverbose::Bool, debug::Bool, extra_kw_args::T) where {T<:Dict}
    ######################
end

function com2OS(M::AbstractEquilibrium, data::Array{Float64,4}, E_array::AbstractVector, mu_matrix::AbstractMatrix, Pphi_matrix::AbstractMatrix, FI_species::String; npm::Int64=length(mu_matrix[1,:]), nRm::Int64=length(Pphi_matrix[1,:]), wall::Union{Nothing,Boundary{Float64}}=nothing, pm_array::Union{Nothing,AbstractArray}=nothing, Rm_array::Union{Nothing,AbstractArray}=nothing, dataIsTopoMap::Bool=false, needJac::Bool=false, transform = x -> x, verbose::Bool=false, debug::Bool=false, kwargs...)
    if !(pm_array===nothing)
        npm = length(pm_array)
    else
        pm_array = range(-1.0,stop=1.0,length=npm)
    end
    if !(Rm_array===nothing)
        nRm = length(Rm_array)
    else
        if !(wall===nothing)
            Rm_array = vec(range((4*magnetic_axis(M)[1]+minimum(wall.r))/5, stop=maximum(wall.r), length=nRm))
        else
            R_lims, z_lims = limits(M)
            Rm_array = vec(range((4*magnetic_axis(M)[1]+R_lims[1])/5, stop=R_lims[2], length=nRm))
        end
    end

    # Before computing an (E,pm,Rm) topological map, assume that lost and incomplete orbits are not of interest
    distinguishLost = false
    distinguishIncomplete = false
    if dataIsTopoMap # If data is a topological map
        if 7 in data # If 7 is present, then we know that lost orbits have been distinguished, and are of interest
            distinguishLost = true
        end
        if 6 in data # If 6 is present, then we know that incomplete orbits have been distinguished, and are of interest
            distinguishIncomplete = true
        end
    end
    ## DON'T WANT TO DO THIS FOR EACH WEIGHT FUNCTION
    ## DON'T WANT TO DO THIS FOR EACH WEIGHT FUNCTION
    ## DON'T WANT TO DO THIS FOR EACH WEIGHT FUNCTION
    topoMap = getOSTopoMap(M, E_array, pm_array, Rm_array; FI_species=FI_species, wall=wall, distinguishLost=distinguishLost, distinguishIncomplete=distinguishIncomplete, kwargs... )
    if dataIsTopoMap # If data is topoMap, then we are done!
        return topoMap, E_array, pm_array, Rm_array
    end

    # If data is NOT a topological map, we need to extract valid (E,pm,Rm) coordinates to map to (E,mu,Pphi;sigma)
    ctrgo_inds = findall(x-> x<0, pm_array)
    cogo_inds = findall(x-> x>=0, pm_array)
    data_OS = zeros(length(E_array), length(pm_array), length(Rm_array))
    if debug
        if needJac
            Hc_OS = Array{HamiltonianCoordinate{ForwardDiff.Dual{Nothing, Float64, 3}}}(undef,size(data_OS))
        else
            Hc_OS = Array{HamiltonianCoordinate}(undef,size(data_OS))
        end
    end
    numObadInds = 0
    @showprogress 1 "Mapping energy slices... " for (iE,E) in enumerate(E_array)
        valid_OSinds = findall(x-> (x!=9) && (x!=7) && (x!=6), topoMap[iE,:,:])
        valid_ctrgo_OSinds = filter(x-> x[1] in ctrgo_inds, valid_OSinds)
        valid_cogo_OSinds = filter(x-> x[1] in cogo_inds, valid_OSinds)

        # Counter-going
        nodes = (mu_matrix[iE,:], Pphi_matrix[iE,:])
        itp = interpolate(nodes, data[iE,:,:,1], Gridded(Linear()))
        if needJac
            valid_ctrgo_OScoords = map(x-> (ForwardDiff.Dual(E,(1.0,0.0,0.0)), ForwardDiff.Dual(pm_array[x[1]],(0.0,1.0,0.0)), ForwardDiff.Dual(Rm_array[x[2]],(0.0,0.0,1.0))), valid_ctrgo_OSinds)
            #mu_values = Array{ForwardDiff.Dual{Nothing, Float64, 3}}(undef, length(valid_ctrgo_OScoords))
            #Pphi_values = Array{ForwardDiff.Dual{Nothing, Float64, 3}}(undef, length(valid_ctrgo_OScoords))
        else
            valid_ctrgo_OScoords = map(x-> (E, pm_array[x[1]], Rm_array[x[2]]), valid_ctrgo_OSinds)
            #mu_values = Array{Float64}(undef,length(valid_ctrgo_OScoords))
            #Pphi_values = Array{Float64}(undef,length(valid_ctrgo_OScoords))
        end

        for (i,OS_coord) in enumerate(valid_ctrgo_OScoords)
            myEPRc = EPRCoordinate(M, OS_coord[1], OS_coord[2], OS_coord[3]; amu=getSpeciesAmu(FI_species), q=getSpeciesEcu(FI_species))
            myHc = HamiltonianCoordinate(M, myEPRc)

            #mu_values[i] = myHc.mu
            #Pphi_values[i] = myHc.p_phi

            try
                x = transform([E,myHc.mu,myHc.p_phi]) # Probably redundant
                jac = (isempty(ForwardDiff.partials(myHc.mu))) ? 1.0 : max(abs(det(hcat((ForwardDiff.partials(xx) for xx in x)...))),0.0)
                data_OS[iE, (valid_ctrgo_OSinds[i])[1], (valid_ctrgo_OSinds[i])[2]] = jac * itp(ForwardDiff.value(myHc.mu), ForwardDiff.value(myHc.p_phi))
            catch
                numObadInds += 1
                debug && println("BAD INTERPOLATION COORDINATE ---> E: $(ForwardDiff.value(OS_coord[1]))  pm: $(ForwardDiff.value(OS_coord[2]))  Rm: $(ForwardDiff.value(OS_coord[3]))  mu: $(myHc.mu)  Pphi: $(myHc.p_phi)")
            end
            debug && (Hc_OS[iE, (valid_ctrgo_OSinds[i])[1], (valid_ctrgo_OSinds[i])[2]] = myHc )
        end

        # Co-going
        itp = interpolate(nodes, data[iE,:,:,2], Gridded(Linear()))
        if needJac
            valid_cogo_OScoords = map(x-> (ForwardDiff.Dual(E,(1.0,0.0,0.0)), ForwardDiff.Dual(pm_array[x[1]],(0.0,1.0,0.0)), ForwardDiff.Dual(Rm_array[x[2]],(0.0,0.0,1.0))), valid_cogo_OSinds)
            #mu_values = Array{ForwardDiff.Dual{Nothing, Float64, 3}}(undef, length(valid_ctrgo_OScoords))
            #Pphi_values = Array{ForwardDiff.Dual{Nothing, Float64, 3}}(undef, length(valid_ctrgo_OScoords))
        else
            valid_cogo_OScoords = map(x-> (E, pm_array[x[1]], Rm_array[x[2]]), valid_cogo_OSinds)
            #mu_values = Array{Float64}(undef,length(valid_ctrgo_OScoords))
            #Pphi_values = Array{Float64}(undef,length(valid_ctrgo_OScoords))
        end

        for (i,OS_coord) in enumerate(valid_cogo_OScoords)
            myEPRc = EPRCoordinate(M, OS_coord[1], OS_coord[2], OS_coord[3]; amu=getSpeciesAmu(FI_species), q=getSpeciesEcu(FI_species))
            myHc = HamiltonianCoordinate(M, myEPRc)

            #mu_values[i] = myHc.mu
            #Pphi_values[i] = myHc.p_phi

            try
                x = transform([E,myHc.mu,myHc.p_phi]) # Probably redundant
                jac = (isempty(ForwardDiff.partials(myHc.mu))) ? 1.0 : max(abs(det(hcat((ForwardDiff.partials(xx) for xx in x)...))),0.0)
                data_OS[iE, (valid_cogo_OSinds[i])[1], (valid_cogo_OSinds[i])[2]] = jac * itp(ForwardDiff.value(myHc.mu), ForwardDiff.value(myHc.p_phi))
            catch
                numObadInds += 1
                debug && println("BAD INTERPOLATION COORDINATE ---> E: $(ForwardDiff.value(OS_coord[1]))  pm: $(ForwardDiff.value(OS_coord[2]))  Rm: $(ForwardDiff.value(OS_coord[3]))  mu: $(myHc.mu)  Pphi: $(myHc.p_phi)")
            end
            debug && (Hc_OS[iE, (valid_cogo_OSinds[i])[1], (valid_cogo_OSinds[i])[2]] = myHc )
        end
    end

    if numObadInds > 0
        @warn "The interpolation algorithm failed to interpolate $(numObadInds) (E,mu,Pphi;sigma) points."
    end

    if debug
        return data_OS, E_array, pm_array, Rm_array, Hc_OS, valid_cogo_OScoords
    end

    return data_OS, E_array, pm_array, Rm_array
end

"""
    com2OS(M, data, E_array, mu_matrix, Pphi_matrix, FI_species)
    com2OS(M, data, E_array, mu_matrix, Pphi_matrix, FI_species; npm=length(mu_matrix[1,:]), nRm=length(Pphi_matrix[1,:]), wall=nothing, dataIsTopoMap=false, topoMap=nothing, verbose=false, kwargs...)

This function maps a constants-of-motion (E,mu,Pphi;sigma) quantity into orbit space (E,pm,Rm). The input data is given as a 5D quantity, and the function will assume it's a collection of 4D (E,mu,Pphi;sigma) quantities with
the last four dimensions corresponding to E, mu, Pphi and sigma, in that order. sigma is a binary coordinate (+1 or -1) that keeps track of co- (+1) and counter-passing (-1) orbits that has the same (E,mu,Pphi) triplet. The first Julia index (1) 
corresponds to sigma=-1 and the second Julia index (2) corresponds to sigma=+1. Return transformed data in OS, along with E-array, pm-array and Rm-array. The output will be 4D with the last three dimensions corresponding to
E, pm and Rm, in that order.

If you have data for your tokamak first wall, it needs to be provided via the 'wall' keyword input argument. From an EFIT magnetic equilibrium, you get the wall via 
M, wall = read_geqdsk(filepath_eqdsk,clockwise_phi=true/false)
From a Solov'ev equilibrium S, you get the wall via boundary(S). Please see Equilibrium.jl/boundary.jl for more info.

If the input data is a topological map, you need to pass the extra keyword argument 'dataIsTopoMap=true'.

If the input data is a fast-ion distribution, the mapping will need a Jacobian. For this, you pass the extra keyword argument needJac=true.

Other keyword arguments can also be passed. Please see below.

The mu_matrix and Pphi_matrix inputs are 2D arrays where the first dimension has length equal to length(E_array).
To explain via example, mu_matrix[3,:] and Pphi_matrix[3,:] correspond to the mu- and Pphi-grid points for the 
energy E_array[3]. This is to keep the number of (mu,Pphi)-grid points the same for all energies, while also minimizing
the number of invalid/impossible (E,mu,Pphi;sigma) coordinates.

kwargs... => dataIsTopoMap, needJac, vverbose, debug, extra_kw_args

dataIsTopoMap - Should be set to true if 'data' is a topological map.
needJac - Should be set to true if mapping a quantity that requires a Jacobian e.g. a fast-ion distribution.
vverbose - Should be set to true if you want the algorithm to be VERY talkative!
debug - If set to true, debug mode will be active.
extra_kw_args - A dictionary with extra keyword arguments for the orbit integration algorithm found in GuidingCenterOrbits.jl/orbit.jl, used to finalize (E,mu,Phi;sigma) -> (E,pm,Rm) transform.
"""
function com2OS(M::AbstractEquilibrium, data::Array{Float64, 5}, E_array::AbstractVector, mu_matrix::Array{Float64, 2}, Pphi_matrix::Array{Float64, 2}, FI_species::AbstractString; npm::Int64=length(mu_matrix[1,:]), nRm::Int64=length(Pphi_matrix[1,:]), wall::Union{Nothing,Boundary{Float64}}=nothing, verbose::Bool=false, kwargs...)

    if pm_array===nothing
        pm_array = range(-1.0,stop=1.0,length=npm)
    end
    if Rm_array===nothing
        if !(wall===nothing)
            Rm_array = vec(range((4*magnetic_axis(M)[1]+minimum(wall.r))/5, stop=maximum(wall.r), length=nRm))
        else
            R_lims, z_lims = limits(M)
            Rm_array = vec(range((4*magnetic_axis(M)[1]+R_lims[1])/5, stop=R_lims[2], length=nRm))
        end
    end

    data_OS = zeros(size(data,1),length(E_array),npm, nRm)
    verbose && println("Mapping (E,μ,Pϕ;σ) -> (E,pm,Rm) for 4D quantity 1 of $(size(data,1))... ")
    data_OS[1, :, :, :], E_array, pm_array, Rm_array = com2OS(M, data[1,:,:,:,:], E_array, mu_matrix, Pphi_matrix, FI_species; pm_array=pm_array, Rm_array=Rm_array, wall=wall, verbose=verbose, kwargs...)
    if size(data,1)>1 # To avoid distributed errors
        data_OS[2,:,:,:] = @showprogress 1 "Mapping (E,μ,Pϕ;σ) -> (E,pm,Rm)... "@distributed (+) for iEd in collect(2:size(data,1))
            data_OS_part = zeros(size(data,1),length(E_array),npm, nRm)
            # E_array_part, pm_array_part and Rm_array_part do not matter. But they need to be returned by com2OS.
            data_OS_part[iEd, :, :, :], E_array_part, pm_array_part, Rm_array_part = com2OS(M, data[iEd,:,:,:,:], E_array, mu_matrix, Pphi_matrix, FI_species; pm_array=pm_array, Rm_array=Rm_array, wall=wall, verbose=verbose, kwargs...)
            data_OS_part # Declare ready for reduction (data_OS += data_OS_part but distributed over several CPU cores)
        end
    end

    return data_OS, E_array, pm_array, Rm_array
end

"""
    Ep2VparaVperp(E_array, p_array, Q)
    Ep2VparaVperp(E_array, p_array, Q; my_gcp=GCDeuteron, needJac=false, isTopoMap=false, verbose=false)

Transform a 2D (E,p) quantity Q into (v_para,v_perp) space. If needJac, assume a jacobian is needed (e.g distribution).
If isTopoMap, assume Q is a topological map, implying special care is needed for successful transform.
PLEASE NOTE! E_array MUST be given in keV.
"""
function Ep2VparaVperp(E_array::Vector, p_array::Vector, Q::Matrix{Float64}; my_gcp::AbstractParticle{T}=GCDeuteron(0.0,0.0,0.0,0.0), needJac::Bool=false, isTopoMap::Bool=false, verbose::Bool=false) where {T}
    keV = 1000*(GuidingCenterOrbits.e0)
    min_keV_E_array = minimum(keV .*E_array)
    max_keV_E_array = maximum(keV .*E_array)
    min_p_array = minimum(p_array)
    max_p_array = maximum(p_array)

    vf = (GuidingCenterOrbits.c0)*sqrt(1-(1/((max_keV_E_array/((my_gcp.m)*(GuidingCenterOrbits.c0)^2))+1)^2))
    verbose && println("vf: $(vf) m/s")
    vpara_array = collect(range(-vf, stop=vf, length=Int64(round(sqrt(length(E_array)*length(p_array))))))
    vperp_array = collect(range(eps(), stop=vf, length=Int64(round(sqrt(length(E_array)*length(p_array))))))
    Q_VEL = zeros(length(vpara_array),length(vperp_array))
    nodes = (keV .*E_array, p_array) # Assume E_array in keV
    itp = interpolate(nodes, Q, Gridded(Linear()))
    itp = Interpolations.extrapolate(itp,isTopoMap ? 9.0 : 0.0)
    E_rep = zeros(length(E_array)*length(p_array))
    p_rep = zeros(length(E_array)*length(p_array))
    coords = CartesianIndices(Q)
    ind = 1
    for coord in coords
        E_rep[ind] = keV .*E_array[coord[1]]
        p_rep[ind] = p_array[coord[2]]
        ind += 1
    end

    brutetree = BruteTree(hcat(E_rep,p_rep)')
    for (ivpara, vpara) in enumerate(vpara_array), (ivperp, vperp) in enumerate(vperp_array)
        v = sqrt(vpara^2 + vperp^2)
        Eq = ((my_gcp.m)*(GuidingCenterOrbits.c0)^2)*(sqrt(1/(1-(v/(GuidingCenterOrbits.c0))^(2)))-1) # Relativistic energy
        pq = vpara/v
        
        #if (Eq < min_keV_E_array) || (Eq > max_keV_E_array) || (pq < min_p_array) || (pq > max_p_array) # If outside of bounds
            verbose && println("Point out-of-bounds found!")
            #idx, dist = nn(brutetree, [Eq, pq])
            #interp_value = Q[coords[idx]]
        #else
        interp_value = itp(Eq,pq)
        if isTopoMap # For topological maps, we need special treatment
            diffE = (keV .*E_array) .- Eq
            diffp = p_array .- pq
            sorted_diffE_inds = sortperm(abs.(diffE))
            sorted_diffp_inds = sortperm(abs.(diffp))
            ibest_E = sorted_diffE_inds[1]
            isecondbest_E = sorted_diffE_inds[2]
            ibest_p = sorted_diffp_inds[1]
            isecondbest_p = sorted_diffp_inds[2]
            closest_values = [Q[ibest_E,ibest_p], Q[ibest_E,isecondbest_p], Q[isecondbest_E,ibest_p], Q[isecondbest_E,isecondbest_p]]
            interp_value = closest_values[closest_index(closest_values,interp_value)]
        end
        #end
        jac = (needJac ? ((my_gcp.m) * vperp/v)/keV : 1.0)
        Q_VEL[ivpara,ivperp] = jac * interp_value # Including Jacobian
    end

    return vpara_array, vperp_array, Q_VEL
end