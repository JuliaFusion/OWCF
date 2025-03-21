################################################################ calcOrbSpec.jl ###################################################################

#### Description:
# This helper script is comprised of functions that take orbit(s), extract their (E,p,R,z) points, weight them and then sends them into a forward model.
# The forward model then acts as a 'black box' and computes the expected signal for that orbit. It is then returned to the functions in 
# calcOrbSpec.jl who, in turn, return them to whatever script requested it.
#
# As of this version of the OWCF, the functions of calcOrbSpec.jl can compute expected signals via the following two forward models:
# - Using Monte Carlo sampling with the DRESS code (J. Eriksson et al, Computer Physics Communications 199 (2016) 40–46)
# - Using analytic formulas developed by A. Valentini ([https://doi.org/10.1088/1741-4326/adc1df, https://doi.org/10.1088/1741-4326/ad9bc8])
#
# As of this version of the OWCF, calcOrbSpec() has to be used/included ONLY in calcOrbWeights.jl or ps2WF.jl, after the Julia and Python packages, 
# ,and scripts, have been loaded. That is, the user needs to first load all the Julia packages (with @everywhere if you use the Julia package Distributed).
# Then, load all the Python packages (with @everywhere ...). Use py""" import ... """ and make sure to include the all necessary Python scripts. Please
# see calcOrbWeights.jl for an example. Finally, load this script with '@everywhere begin include("calcOrbSpec.jl") end'

# Script written by Henrik Järleblad and Andrea Valentini. Last maintained 2025-03-07.
###################################################################################################################################################

include("../misc/temp_n_dens.jl") # Load 'analytical' thermal species temperature and thermal species density profiles
include("../misc/availReacts.jl") # Load the possibility to treat fusion reaction strings
include("../vcone.jl") # Load the possibility to model diagnostic sightlines (in Julia)
include("../forward.jl") # Load the possibility to compute synthetic signals without Monte Carlo methods

"""
    calcOrbSpec(M, o, nfast, forward, thermal_dist, Ed_bin_edges, reaction)
    calcOrbSpec(-||-; o_interp_length=500, thermal_temp=3.0, thermal_dens=1.0e19, analytic=false, debug=false)

Calculate the expected diagnostic spectrum for one (drift) orbit. The input arguments are as follows:
- M - An axisymmetric equilibrium from the Equilibrium.jl Julia package. It is used to utilize its magnetic field.
- o - The orbit for which to compute the expected diagnostic spectrum. It is an Orbit struct from GuidingCenterOrbits.jl/orbit.jl.
- nfast - The number of fast ions on the orbit. Set to 1, by default.
- forward - The Forward Python DRESS object from the forward.py DRESS script. Used to compute the expected diagnostic spectrum. When the 
            'reaction input variable is specified as "b" where b is a fast-ion particle species, this input variable should instead be specified 
            as a String with the filepath to an output file from the LINE21 code or the OWCF/extra/createCustomLOS.jl script.
- thermal_dist -  The Thermal/FBM Python object from the transp_dists.py script. Used to represent the thermal species distribution. Can also be nothing.
- Ed_bin_edges - The diagnostic measurement bin edges into which the synthetic measurements will be binned (in keV or m/s).
- reaction - Fusion reaction, on any of the forms described in the OWCF/misc/availReacts.jl script.
The keyword arguments are as follows:
- o_interp_length - The number of time points onto which the orbit trajectory will be interpolated. Equidistant points in time are used.
- thermal_temp - If the 'thermal_dist' input variable is nothing, the 'thermal_temp' keyword argument will be used to model the thermal particle species temperature.
                 If the 'thermal_temp' keyword argument is specified as an Int64/Float64, the default OWCF temperature profile will be used (please see OWCF/misc/default_temp_n_dens.png)
                 with the 'thermal_temp' keyword argument value used as the thermal particle species temperature value (keV assumed) on-axis. If the 'thermal_temp' keyword argument is 
                 specified as an extrapolation object, that will be used to compute the thermal particle species temperature value at all ρ_pol values of interest. 
- thermal_dens - If the 'thermal_dist' input variable is nothing, the 'thermal_dens' keyword argument will be used to model the thermal particle species density.
                 If the 'thermal_dens' keyword argument is specified as an Int64/Float64, the default OWCF density profile will be used (please see OWCF/misc/default_temp_n_dens.png)
                 with the 'thermal_dens' keyword argument value used as the thermal particle species density value (m^-3 assumed) on-axis. If the 'thermal_dens' keyword argument is 
                 specified as an extrapolation object, that will be used to compute the thermal particle species density value at all ρ_pol values of interest. 
- analytic - If set to true, the expected diagnostic spectrum will be computed using analytic equations, instead of Monte Carlo sampling. If set to true, !getReactionForm(reaction)==3 must be true (OWCF/misc/availReacts.jl)
- debug - A boolean debug input variable. If set to true, the function will run in debug-mode.
"""
function calcOrbSpec(M::AbstractEquilibrium, o::Orbit{Float64, EPRCoordinate{Float64}}, nfast::Float64, forward::Union{PyObject,AbstractString}, thermal_dist::Union{Nothing,PyObject}, Ed_bin_edges::AbstractArray, reaction::AbstractString; 
                     o_interp_length=500, thermal_temp::Union{Nothing,Float64,Int64,Interpolations.Extrapolation,Interpolations.FilledExtrapolation}=3.0, 
                     thermal_dens::Union{Nothing,Float64,Int64,Interpolations.Extrapolation,Interpolations.FilledExtrapolation}=1.0e19, analytic::Bool=false, 
                     debug::Bool=false)

    if analytic
        # If the orbit spectrum is to be computed via analytic equations, the thermal particle species temperature needs to be assumed to be zero
        thermal_temp = eps() # For numerical reasons, eps() and not 0
    end

    # Ensure that every orbit has o_interp_length number of EpRz points to serve as input to the spectrum calculator
    # We can assume that the EpRz points are equally distant in time, because this is ensured by GuidingCenterOrbits.get_orbit()
    o_spline = OrbitTomography.OrbitSpline(o, range(0.0, o.tau_p, length=length(o))) # Make a spline to enable interpolation onto o_interp_length points
    tpl = range(0.0, o.tau_p, length=o_interp_length) # The time 'grid' we want to interpolate onto
    o_energy = zeros(length(tpl)) # All the energy values of the points along the orbit. Vector.
    o_pitch = zeros(length(tpl)) # All the pitch values of the points along the orbit. Vector.
    o_R = zeros(length(tpl)) # All the R values of the points along the orbit. Vector.
    o_z = zeros(length(tpl)) # All the z values of the points along the orbit. Vector.
    for i in eachindex(tpl)
        Si = o_spline(tpl[i]) #Interpolate onto the tpl time points
        o_energy[i] = Si[1] # Store energy at point i
        o_pitch[i] = Si[2] # Store pitch at point i
        o_R[i] = Si[3] # Store R position at point i
        o_z[i] = Si[4] # Store z position at point i
    end
    o_w = eltype(o_R)[(tpl[min(i+1,o_interp_length)] - tpl[i]) for i=1:o_interp_length] # incremental time (dt) weights for all EpRz points
    o_w = o_w .* (nfast/o.tau_p) # Normalize the time weights by the poloidal transit time (and also multiply by some number of fast ions, nfast. nfast is usually set to 1.0)
    o_B = zeros(3,length(o_R)) # Pre-allocate the magnetic field vector at every point along the orbit
    vec_Rz = zip(o_R,o_z) # zip the R-values and the z-values together in a Tuple-vector


    if debug
        # WRITE CODE HERE FOR DEBUGGING
    end
    # If you have not specified a TRANSP thermal particle species distribution, and you are not computing projected velocity orbit spectra,..
    if isnothing(thermal_dist) && !(getReactionForm(reaction)==3)
        thermal_temp_interp = zeros(length(o_energy))
        thermal_dens_interp = zeros(length(o_energy))
        for (i,crs_c) in enumerate(vec_Rz) # Go through them all and their indices
            myB = Equilibrium.Bfield(M,crs_c[1],crs_c[2]) # Calculate the B-field vector at the (R,z) point
            ψ_rz = M(crs_c[1],crs_c[2])
            psi_on_axis, psi_at_bdry = psi_limits(M)
            ρ_pol_rz = sqrt(ψ_rz-psi_on_axis)/(psi_at_bdry-psi_on_axis) # The formula for the normalized flux coordinate ρ_pol = (ψ-ψ_axis)/(ψ_edge-ψ_axis)
            if debug
                # WRITE CODE HERE FOR DEBUGGING
            end
            if typeof(thermal_temp) <: Union{Interpolations.Extrapolation,Interpolations.FilledExtrapolation} # If thermal_temp is an interpolation/extrapolation object
                thermal_temp_interp[i] = thermal_temp(ρ_pol_rz) # Interpolate
            else # If it is not (and thus, it must be a Float64/Int64)
                thermal_temp_interp[i] = getAnalyticalTemp(thermal_temp,ρ_pol_rz) # Use it as the temperature on axis
            end
            if typeof(thermal_dens) <: Union{Interpolations.Extrapolation,Interpolations.FilledExtrapolation} # If thermal_dens is an interpolation/extrapolation object
                thermal_dens_interp[i] = thermal_dens(ρ_pol_rz) # Interpolate
            else # If it is not (and thus, it must be a Float64/Int64)
                thermal_dens_interp[i] = getAnalyticalDens(thermal_dens,ρ_pol_rz) # And use thermal_dens as the density on axis
            end
            o_B[1,i] = myB[1] # Store the BR
            o_B[2,i] = myB[2] # Store the Bphi
            o_B[3,i] = myB[3] # Store the Bz
        end
    else
        thermal_temp_interp = "" # If thermal_dist is defined, then we don't need to manually define thermal_temp (forward.py will take care of it automatically)
        thermal_dens_interp = "" # If thermal_dist is defined, then we don't need to manually define thermal_dens (forward.py will take care of it automatically)
        for (i,crs_c) in enumerate(vec_Rz) # Go through them all and their indices
            myB = Equilibrium.Bfield(M,crs_c[1],crs_c[2]) # Calculate the B-field vector at the (R,z) point
            o_B[1,i] = myB[1] # Store the BR
            o_B[2,i] = myB[2] # Store the Bphi
            o_B[3,i] = myB[3] # Store the Bz
        end
    end

    if debug
        # WRITE CODE TO RETURN DEBUG QUANTITIES OF INTEREST. For example thermal_temp, o_R etc
        return thermal_temp_interp, o_R
    end

    if analytic
        if isfile(forward)
            viewing_cone = ViewingCone(forward)
        else
            error("Analytical model with spherical 4*pi emission is NOT supported! Please correct and re-try.")
        end
        spec = forward_calc(viewing_cone, o_energy, clamp.(o_pitch,-1,1), o_R, o_z, o_w, Ed_bin_edges, o_B; reaction=reaction, bulk_dens=thermal_dens_interp)
        return spec
    else 
        spec = zeros(length(Ed_bin_edges)-1)
        stats = 1
        energy_state = occursin("-",reaction) ? getEmittedParticleEnergyLevel(reaction) : "GS"
        
        if !(energy_state=="GS")
            stats = 100 # Need a lot more statistics for two-step gamma-ray reactions
        end
        for i=1:stats
            py"""
            forward = $forward # Convert the Forward Python object from Julia to Python.
            # Please note, the '$' symbol is used below to convert objects from Julia to Python. Even Julia PyObjects
            spectrum = forward.calc($o_energy, $o_pitch, $o_R, $o_z, $o_w, $thermal_dist, $Ed_bin_edges, $o_B, n_repeat=50, reaction=$reaction, bulk_temp=$thermal_temp_interp, bulk_dens=$thermal_dens_interp) # Please see the forward.py script, for further specifications
            """
            spec += vec(py"spectrum")
        end
        return spec ./ stats
    end
end

"""
    calcOrbSpecs(M, og_orbs, F_os, forward, thermal_dist, Ed_bin_edges, reaction)
    calcOrbSpecs(-||-; thermal_temp=3.0, thermal_dens=1.0e19, distributed=false, visualizeProgress=false, verbose=false, kwargs... )

Like calcOrbSpec(), but with many orbits. Return all the spectra as a matrix. This function is essentially the script calcOrbWeights.jl, 
but re-written as a callable function.

Calculate the expected diagnostic spectra of the orbits in og_orbs. The inputs are as follows
- M - An axisymmetric equilibrium from the Equilibrium.jl Julia package. It is used to utilize its magnetic field.
- og_orbs - The orbits for which to compute the expected synthetic spectra. It is a vector with orbit structs from GuidingCenterOrbits.jl/orbit.jl.
- F_os - A vector of length(og_orbs) where each element is the number of fast ions on the corresponding orbit in og_orbs. Usually, the value of every element is set to 1.0, by default.
- forward - The Forward Python DRESS object from the forward.py DRESS script. Used to compute the expected synthetic spectrum.
- thermal_dist -  The Thermal/FBM Python object from the transp_dists.py script. Used to represent the thermal species distribution. Can also be nothing.
- Ed_bin_edges - The diagnostic measurement bin edges into which the synthetic measurements will be binned (in keV or m/s).
- reaction - Fusion reaction, on any of the forms described in the OWCF/misc/availReacts.jl script.
Keyword arguments include
- distributed - If true, the orbit spectra will be computed using multiple CPU/GPU processes. If false, single-threaded computing will be used
- visualizeProgress - If true, a progress bar will be displayed when computing orbit spectra
- verbose - If true, the function execution will be very talkative
"""
function calcOrbSpecs(M::AbstractEquilibrium, og_orbs::Vector{Orbit{Float64, EPRCoordinate{Float64}}}, F_os::Union{Array{Int64,1},Array{Float64,1}}, forward::PyObject, thermal_dist::Union{PyObject,Nothing}, Ed_bin_edges::AbstractArray, reaction::AbstractString; distributed::Bool=false, visualizeProgress::Bool=false, verbose::Bool=false, kwargs...)

    norbs = length(og_orbs)
    if distributed # If parallel computating is desired...
        if visualizeProgress # if you want the progress to be visualized...
            p = Progress(norbs) # Create a progress bar that is 'norbs' long
            channel = RemoteChannel(()->Channel{Bool}(norbs), 1) # Utilize a channel
            Wtot = fetch(@sync begin
                @async while take!(channel)
                    ProgressMeter.next!(p)
                end
                @async begin
                    W = @distributed (+) for i=1:norbs
                        spec = calcOrbSpec(M, og_orbs[i], F_os[i], forward, thermal_dist, Ed_bin_edges, reaction; kwargs...) # Calculate the diagnostic energy spectrum for it
                        rows = append!(collect(1:length(spec)),length(spec)) # To be able to tell the sparse framework about the real size of the weight matrix
                        cols = append!(i .*ones(Int64, length(spec)), norbs) # To be able to tell the sparse framework about the real size of the weight matrix

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
            Wtot = @distributed (+) for i=1:norbs
            spec = calcOrbSpec(M, og_orbs[i], F_os[i], forward, thermal_dist, Ed_bin_edges, reaction; kwargs...) # Calculate the diagnostic energy spectrum for it
                rows = append!(collect(1:length(spec)),length(spec)) # Please see similar line earlier in the script
                cols = append!(i .*ones(Int64, length(spec)), norbs) # Please see similar line earlier in the script

                # NOTE! Important to have this append!(spec,0.0) line AFTER the rows and cols assignment lines
                append!(spec,0.0) # Please see similar line earlier in the script
                Wi = dropzeros(sparse(rows,cols,spec)) # Please see similar line earlier in the script
                Wi # Please see similar line earlier in the script
            end
        end
    else # ... if you do not use multiple cores, good luck!
        spec = calcOrbSpec(M, og_orbs[1], F_os[1], forward, thermal_dist, Ed_bin_edges, reaction; kwargs...) # Calculate the diagnostic energy spectrum for the first orbit
        rows = append!(collect(1:length(spec)),length(spec)) # # Please see similar line earlier in the script
        cols = append!(1 .*ones(Int64, length(spec)), norbs) # # Please see similar line earlier in the script

        # NOTE! Important to have this append!(spec,0.0) line AFTER the rows and cols assignment lines
        append!(spec,0.0)
        Wtot = dropzeros(sparse(rows,cols,spec)) # Pre-allocate a sparse matrix

        for i=2:norbs
            verbose && println("Calculating spectra for orbit $(i) of $(norbs)... ")
            local spec = calcOrbSpec(M, og_orbs[i], F_os[i], forward, thermal_dist, Ed_bin_edges, reaction; kwargs...) # Calculate the diagnostic energy spectrum for it
            local rows = append!(collect(1:length(spec)),length(spec)) # Please see similar line earlier in the script
            local cols = append!(i .*ones(Int64, length(spec)), norbs) # Please see similar line earlier in the script

            # NOTE! Important to have this append!(spec,0.0) line AFTER the rows and cols assignment lines
            append!(spec,0.0) # Please see similar line earlier in the script
            Wtot += dropzeros(sparse(rows,cols,spec)) # Please see similar line earlier in the script
        end
    end
    return Wtot
end