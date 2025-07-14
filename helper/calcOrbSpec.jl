################################################################ calcOrbSpec.jl ########################################################################

#### Description:
# This helper script is comprised of functions that take guiding-centre drift orbit Structs as input. The orbit Struct is defined in the orbit.jl script
# of the GuidingCenterOrbits.jl Julia package. The functions extract the (E,p,R,z) points of the orbit Struct, weight them and then sends them into 
# a forward model. The functions assume that the (E,p,R,z) points define charged guiding-centre particle trajectories inside a toroidally symmetric 
# magnetic confinement fusion device with
#
#   - E - The energy value(s) for the orbit, in keV
#   - p - The pitch values -||- (v_||/v where v_|| is the velocity component parallel to the magnetic field and v is the speed)
#   - R - The major radius -||-, in meters
#   - z - The vertical -||-, in meters
#
# The forward model then acts as a 'black box' and computes the expected signal for that orbit. It is then returned to the functions in 
# calcOrbSpec.jl which, in turn, return them to whatever script requested it.
#
# As of this version of the OWCF, the forward model used as input to calcOrbSpec() and calcOrbSpecs() needs to be on the form
#
#   forwardmodel(E, p, R, z, w)
#
# where the inputs are
#
#   - E - The energy points - AbstractVector with length(E)=N
#   - p - The pitch (p=v_||/v) points - AbstractVector with length(P)=N
#   - R - The major radius points - AbstractVector with length(R)=N
#   - z - The vertical points - AbstractVector with length(z)=N
#   - w - The weights of the points - AbstractVector with length(w)=N
#
# The output of the forwardmodel() function needs to be a vector of length M.

# Script written by Henrik Järleblad and Andrea Valentini. Last maintained 2025-07-11.
########################################################################################################################################################

# The Julia packages used in this file
using Distributed
using GuidingCenterOrbits
using OrbitTomography
using ProgressMeter
using SparseArrays

"""
    calcOrbSpec(o, n_fast, forwardmodel)
    calcOrbSpec(-||-; o_interp_length=500, debug=false)

Calculate the expected diagnostic spectrum for one (drift) orbit. The input arguments are as follows:

- o - The orbit Struct for which to compute the expected diagnostic spectrum. It is an Orbit struct from the Julia package GuidingCenterOrbits.jl.
- n_fast - The number of fast ions on the orbit. Set to 1, by default.
- forwardmodel - A function of the form forwardmodel = f(E,p,R,z,w) that outputs a Vector, as explained in the beginning of this file.

The keyword arguments are as follows:

- o_interp_length - The number of time points onto which the orbit trajectory will be interpolated. Equidistant points in time are used. - Int64
- debug - A boolean debug input variable. If set to true, the function will run in debug-mode.
"""
function calcOrbSpec(o::T where {T<:GuidingCenterOrbits.Orbit}, n_fast::T1 where {T1<:Real}, forwardmodel;
                     o_interp_length::Int64=200, debug::Bool=false)

    # Ensure that every orbit has o_interp_length number of EpRz points to serve as input to the spectrum calculator
    # We can assume that the EpRz points are equally distant in time, because this is ensured by GuidingCenterOrbits.get_orbit()
    o_spline = OrbitTomography.OrbitSpline(o, range(0.0, o.tau_p, length=length(o))) # Make a spline to enable interpolation onto o_interp_length points
    tpl = range(0.0, o.tau_p, length=o_interp_length) # The time 'grid' we want to interpolate onto
    o_energy = zeros(length(tpl)) # All the energy values of the points along the orbit. Vector.
    o_pitch = zeros(length(tpl)) # All the pitch values of the points along the orbit. Vector.
    o_R = zeros(length(tpl)) # All the R values of the points along the orbit. Vector.
    o_z = zeros(length(tpl)) # All the z values of the points along the orbit. Vector.
    for i in eachindex(tpl)
        Si = o_spline(tpl[i]) # Interpolate onto the tpl time points
        o_energy[i] = Si[1] # Store energy at point i
        o_pitch[i] = Si[2] # Store pitch at point i
        o_R[i] = Si[3] # Store R position at point i
        o_z[i] = Si[4] # Store z position at point i
    end
    o_w = eltype(o_R)[(tpl[min(i+1,o_interp_length)] - tpl[i]) for i=1:o_interp_length] # Incremental time (dt) weights for all EpRz points
    o_w = o_w .* (n_fast/o.tau_p) # Normalize the time weights by the poloidal transit time (and also multiply by some number of fast ions, n_fast. n_fast is usually set to 1.0)

    if debug
        # WRITE CODE HERE FOR DEBUGGING
    end
    
    return forwardmodel(o_energy, o_pitch, o_R, o_z, o_w)
end

"""
    calcOrbSpecs(orbs, orbs_n_fast, forwardmodel)
    calcOrbSpecs(-||-; distributed=false, visualizeProgress=false, verbose=false, kwargs... )

Like calcOrbSpec(), but with many orbits. Return all the spectra as a matrix. This function is essentially the script calcOrbWeights.jl, 
but re-written as a callable function.

Calculate the expected diagnostic spectra of the orbits in og_orbs. The inputs are as follows:

- orbs - The orbits for which to compute the expected synthetic spectra. It is a vector with orbit Structs from GuidingCenterOrbits.jl/orbit.jl.
- orbs_n_fast - A vector of length(orbs) where each element is the number of fast ions on the corresponding orbit in 'orbs'. The value of every element is set to 1.0, by default.
- forwardmodel - A function of the form forwardmodel = f(E,p,R,z,w) that outputs a Vector, as explained in the beginning of this file.

The keyword arguments are:

- distributed - If true, the orbit spectra will be computed using multiple CPU processes. If false, single-threaded computing will be used
- visualizeProgress - If true, a progress bar will be displayed when computing orbit spectra
- verbose - If true, the function execution will be very talkative

"""
function calcOrbSpecs(orbs::Vector{T} where {T<:GuidingCenterOrbits.Orbit}, orbs_n_fast::Vector{T1} where {T1<:Real}, forwardmodel;
                      distributed::Bool=false, visualizeProgress::Bool=false, verbose::Bool=false, kwargs...)

    n_orbs = length(orbs)
    if distributed # If parallel computating is desired...
        if visualizeProgress # if you want the progress to be visualized...
            p = Progress(n_orbs) # Create a progress bar that is 'n_orbs' long
            channel = RemoteChannel(()->Channel{Bool}(n_orbs), 1) # Utilize a channel
            Wtot = fetch(@sync begin
                @async while take!(channel)
                    ProgressMeter.next!(p)
                end
                @async begin
                    W = @distributed (+) for i=1:n_orbs
                        spec = calcOrbSpec(orbs[i], orbs_n_fast[i], forwardmodel; kwargs...) # Calculate the diagnostic energy spectrum for it
                        rows = append!(collect(1:length(spec)),length(spec)) # To be able to tell the sparse framework about the real size of the weight matrix
                        cols = append!(i .*ones(Int64, length(spec)), n_orbs) # To be able to tell the sparse framework about the real size of the weight matrix

                        # NOTE! Important to have this append!(spec,0.0) line AFTER the rows and cols assignment lines
                        append!(spec,0.0) # To be able to tell the sparse framework about the real size of the weight matrix
                        Wi = dropzeros(sparse(rows,cols,spec)) # Create a sparse weight matrix, with the current (i) column non-zero. Remove all redundant zeros.

                        put!(channel, true) # Update the progress bar
                        Wi # Declare this orbit weight to be added to the parallel computing reduction (@distributed (+))
                    end
                    put!(channel, false) # Update the progress bar
                    W # Declare the total weight function as ready for fetching
                end
            end)
        else
            Wtot = @distributed (+) for i=1:n_orbs
            spec = calcOrbSpec(orbs[i], orbs_n_fast[i], forwardmodel; kwargs...) # Calculate the diagnostic energy spectrum for it
                rows = append!(collect(1:length(spec)),length(spec)) # Please see similar line earlier in the script
                cols = append!(i .*ones(Int64, length(spec)), n_orbs) # Please see similar line earlier in the script

                # NOTE! Important to have this append!(spec,0.0) line AFTER the rows and cols assignment lines
                append!(spec,0.0) # Please see similar line earlier in the script
                Wi = dropzeros(sparse(rows,cols,spec)) # Please see similar line earlier in the script
                Wi # Please see similar line earlier in the script
            end
        end
    else # ... if you do not use multiple cores, good luck!
        spec = calcOrbSpec(orbs[1], orbs_n_fast[1], forwardmodel; kwargs...) # Calculate the diagnostic energy spectrum for the first orbit
        rows = append!(collect(1:length(spec)),length(spec)) # # Please see similar line earlier in the script
        cols = append!(1 .*ones(Int64, length(spec)), n_orbs) # # Please see similar line earlier in the script

        # NOTE! Important to have this append!(spec,0.0) line AFTER the rows and cols assignment lines
        append!(spec,0.0)
        Wtot = dropzeros(sparse(rows,cols,spec)) # Pre-allocate a sparse matrix

        for i=2:n_orbs
            verbose && println("Calculating spectra for orbit $(i) of $(n_orbs)... ")
            local spec = calcOrbSpec(orbs[i], orbs_n_fast[i], forwardmodel; kwargs...) # Calculate the diagnostic energy spectrum for it
            local rows = append!(collect(1:length(spec)),length(spec)) # Please see similar line earlier in the script
            local cols = append!(i .*ones(Int64, length(spec)), n_orbs) # Please see similar line earlier in the script

            # NOTE! Important to have this append!(spec,0.0) line AFTER the rows and cols assignment lines
            append!(spec,0.0) # Please see similar line earlier in the script
            Wtot += dropzeros(sparse(rows,cols,spec)) # Please see similar line earlier in the script
        end
    end
    return Wtot
end