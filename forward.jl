########################################## forward.jl ##################################################
# This code is a julia transposition of forward.py and spec.py, with the key difference that the calculations
# do not account for the thermal reactant dynamics. This is in favor of faster computations with no Monte 
# Carlo sampling needed.
# Miscellaneous comments:
#  - the a+b->1+2 notation from DRESS is swapped here with b+t->d+r, meaning that a beam particle 'b' reacts
#    with a target particle 't' and generates a detected particle 'd' and a residual particle 'r'. 
#  - 
#
# Queries for Henrik:
#  - specifying all the argument types for all functions is not necessary here, since most inputs are defined ¨
#    right here in this script. This saves on computation time.
#  - a few packages were added, and so the manifest file will be updated. Also the project one I think.
#  - for now orbit pitch values are clamped in calcOrbSpec before they are fed here; should we do the same with DRESS?
#  - 
#
# TO DO:
#  - check the jacobian for the gamma-ray in spec.py as well. It might not be there
#  - 
#
# THIS SCRIPT IS FUNCTIONALLY DONE. ONLY DOCUMENTATION REMAINS TO BE FINALIZED.
#
# Written by Andrea Valentini. Last maintained 2024-11-06
######################################################################################################

using NaNMath
using LoopVectorization
using CSV
using DataFrames
using LegendrePolynomials  # For Legendre polynomials
using LinearAlgebra

cnst = pyimport("constants")
fsrct = pyimport("fusreact")

# function get_lims(E_b::Vector{Float64}, l_b::Vector{Float64}, 
#                  m_b::Float64, m_t::Float64, m_d::Float64, m_r::Float64)
#     # ... function body ...
# end

# function Doppler_spectrum(E_nominal::Float64, E_g::Float64, phi::Vector{Float64}, 
#                          n_b::Vector{Float64}, E_b::Vector{Float64}, l_b::Vector{Float64}, 
#                          m_b::Float64, dOmega::Float64)
#     # ... function body ...
# end

# function relativistic_inverse_jacobian(γ_b::Vector{Float64}, cosLAB::Vector{Float64}, 
#                                      E_b::Vector{Float64}, m_b::Float64, m_t::Float64, 
#                                      m_d::Float64, m_r::Float64, A::Vector{Float64}, 
#                                      B::Vector{Float64}, D::Vector{Float64})
#     # ... function body ...
# end

# function dn_dtdγ(γ_b::Vector{Float64}, l_d::Vector{Float64}, m_d::Float64, 
#                 n_b::Vector{Float64}, E_b::Vector{Float64}, l_b::Vector{Float64}, 
#                 m_b::Float64, n_t::Vector{Float64}, m_t::Float64, m_r::Float64, 
#                 dSigma_dcosCM::Function)
#     # ... function body ...
# end

# function dn_dtdE(E_d::Vector{Float64}, l_d::Vector{Float64}, m_d::Float64, 
#                 n_b::Vector{Float64}, E_b::Vector{Float64}, l_b::Vector{Float64}, 
#                 m_b::Float64, n_t::Vector{Float64}, m_t::Float64, m_r::Float64, 
#                 dSigma_dcosCM::Function)
#     # ... function body ...
# end

# function SetReaction(reaction_name::String, product_state::String="GS")
#     # ... function body ...
# end

# function GetOrbSpec(viewing_cone::Union{Nothing,ViewingCone}, 
#                    o_E::Vector{Float64}, o_p::Vector{Float64}, 
#                    o_R::Vector{Float64}, o_z::Vector{Float64}, 
#                    o_w::Vector{Float64}, Ed_bins::Vector{Float64}, 
#                    o_B::Matrix{Float64};
#                    product_state::String="GS", 
#                    reaction_name::String="TD_neutrons", 
#                    bulk_dens::Union{Nothing,Vector{Float64}}=1.0e19)
#     # ... function body ...
# end

function get_lims(E_b, l_b, m_b, m_t, m_d, m_r)
    cosLAB = 1.0

    A = (E_b .+ m_b .+ m_t).^2 .- E_b .* (E_b .+ 2 .* m_b) .+ m_d^2 .- m_r^2
    B = E_b .+ m_b .+ m_t
    D = B .* m_d .- A ./ 2
    a = B.^2 .- E_b .* (E_b .+ 2 .* m_b) .* cosLAB.^2
    bi_half = B .* D .- E_b .* (E_b .+ 2 .* m_b) .* m_d .* cosLAB.^2
    ci = D.^2
    delta_fourth = bi_half.^2 .- a .* ci

    cosLABmin = real.(sqrt.(Complex.((B.^2 .* m_d.^2 .- A.^2 ./ 4) ./ (E_b .* (E_b .+ 2 .* m_b) .* m_d.^2))))

    return hcat(max.((-bi_half .- 1.1 .* sqrt.(delta_fourth)) ./ a, 0.001), 
                (-bi_half .+ 1.1 .* sqrt.(delta_fourth)) ./ a, 
                max.(l_b .- 1.1 .* acos.(cosLABmin), 0.001),
                min.(l_b .+ 1.1 .* acos.(cosLABmin), 0.999*π),
                max.(-π/2 .- 1.15 .* acos.(cosLABmin), -0.999*π),
                min.(-π/2 .+ 1.15 .* acos.(cosLABmin), -0.001),)
end

function Doppler_spectrum(E_nominal, E_g, phi, n_b, E_b, l_b, m_b, dOmega)
    G_b = @. E_b / m_b + 1
    G_b2 = @. G_b^2
    costanti = @. n_b / (4 * π^2) * dOmega
    v_b = @. cnst.c * sqrt((G_b2 - 1) / G_b2)
    v_d = cnst.c
    cos_l_b = cos.(l_b)
    sin_l_b = sin.(l_b)
    cos_phi = cos.(phi)
    sin_phi = sin.(phi)

    sin_gyro = @. (( (E_g/E_nominal)^-1 - G_b + sqrt(G_b2 - 1) * cos_l_b * cos_phi ) / 
                  ( sqrt(G_b2 - 1) * sin_l_b * sin_phi ))
    cosLAB = @. cos_l_b * cos_phi - sin_l_b * sin_phi * sin_gyro

    dvd_dcosLAB = 0
    dcosCM_dcosLAB = @. G_b * v_d * ((v_d * (v_d - v_b * cosLAB) + v_b * (1 - cosLAB^2) * dvd_dcosLAB) /
                     ((G_b2 - 1) * v_d^2 * cosLAB^2 - 2 * G_b2 * v_d * v_b * cosLAB + 
                      v_d^2 + G_b2 * v_b^2)^1.5)
    
    sin_gyro_sq = @. ifelse(abs(sin_gyro) > 0.99999, -Inf, sin_gyro^2)
    jac = @. inv(sqrt(1 - sin_gyro_sq)) * abs(E_nominal / E_g^2 / sqrt(G_b2 - 1) / sin_l_b / sin_phi)

    return @. costanti * jac * dcosCM_dcosLAB
end

function relativistic_inverse_jacobian(γ_b, cosLAB, E_b, m_b, m_t, m_d, m_r, A, B, D)
    s = -sign.(cos.(γ_b))

    jac = s .* ( (E_b .* (E_b .+ 2 .* m_b) .* cosLAB .* 
            (-2 .* B.^3 .* D .* m_d .+ D.^2 .* E_b .* (E_b .+ 2 .* m_b) .* cosLAB.^2 
            .- s .* 2 .* B .* D .* (s .* E_b.^2 .* m_d .* cosLAB.^2 .+ s .* 2 .* E_b .* m_b .* m_d .* cosLAB.^2 
            .+ NaNMath.sqrt.(E_b .* (E_b .+ 2 .* m_b) .* cosLAB.^2 .* (D.^2 .- 2 .* B .* D .* m_d .+ E_b .* (E_b .+ 2 .* m_b) .* m_d.^2 .* cosLAB.^2)))
            .+ B.^2 .* (D.^2 .+ 2 .* m_d .* (E_b.^2 .* m_d .* cosLAB.^2 .+ 2 .* E_b .* m_b .* m_d .* cosLAB.^2 
            .+ s .* NaNMath.sqrt.(E_b .* (E_b .+ 2 .* m_b) .* cosLAB.^2 .* (D.^2 .- 2 .* B .* D .* m_d .+ E_b .* (E_b .+ 2 .* m_b) .* m_d.^2 .* cosLAB.^2)))))) 
        ) ./ ((B.^2 .- E_b .* (E_b .+ 2 .* m_b) .* cosLAB.^2).^2 .* NaNMath.sqrt.(E_b .* (E_b .+ 2 .* m_b) .* cosLAB.^2 .* (D.^2 .- 2 .* B .* D .* m_d .+ E_b .* (E_b .+ 2 .* m_b) .* m_d.^2 .* cosLAB.^2)))

    return jac
end

function dn_dtdγ(γ_b, l_d, m_d, n_b, E_b, l_b, m_b, n_t, m_t, m_r, dSigma_dcosCM)
    A = (E_b .+ m_b .+ m_t).^2 .- E_b .* (E_b .+ 2 .* m_b) .+ m_d^2 .- m_r^2
    B = E_b .+ m_b .+ m_t
    D = B .* m_d .- A ./ 2

    cosLAB = cos.(l_b) .* cos.(l_d) .- sin.(l_b) .* sin.(l_d) .* sin.(γ_b)
    a = B.^2 .- E_b .* (E_b .+ 2 .* m_b) .* cosLAB.^2
    bi_half = B .* D .- E_b .* (E_b .+ 2 .* m_b) .* m_d .* cosLAB.^2
    ci = D.^2
    E_d = (-1 .* bi_half .- sign.(cos.(γ_b)) .* NaNMath.sqrt.(bi_half.^2 .- a .* ci)) ./ a

    G_b = E_b ./ m_b .+ 1
    v_b = sqrt.((G_b.^2 .- 1) ./ G_b.^2) * cnst.c
    G_d = E_d ./ m_d .+ 1
    v_d = sqrt.((G_d.^2 .- 1) ./ G_d.^2) * cnst.c
    v_CM = sqrt.(E_b .* (E_b .+ 2 .* m_b)) ./ ((G_b .* m_b .+ m_t) / cnst.c)
    G_CM = 1 ./ sqrt.(1 .- v_CM.^2 / cnst.c^2)

    cosCM = (G_CM .* (v_d .* cosLAB .- v_CM)) ./ sqrt.((G_CM.^2 .- 1) .* v_d.^2 .* cosLAB.^2 .- 2 .* G_CM.^2 .* v_d .* v_CM .* cosLAB .+ v_d.^2 .+ G_CM.^2 .* v_CM.^2)

    dEd_dcosLAB = relativistic_inverse_jacobian(γ_b, cosLAB, E_b, m_b, m_t, m_d, m_r, A, B, D)
    dvd_dcosLAB = cnst.c .* (m_d .* G_d.^2 .* sqrt.(G_d.^2 .- 1)).^-1 .* dEd_dcosLAB
    dcosCM_dcosLAB = G_CM .* v_d .* ((v_d .* (v_d .- v_CM .* cosLAB) .+ v_CM .* (1 .- cosLAB.^2) .* dvd_dcosLAB) ./ ((G_CM.^2 .- 1) .* v_d.^2 .* cosLAB.^2 .- 2 .* G_CM.^2 .* v_d .* v_CM .* cosLAB .+ v_d.^2 .+ G_CM.^2 .* v_CM.^2).^1.5)

    costanti = n_b .* n_t .* v_b ./ π

    dSigma_dcosLAB = dSigma_dcosCM(E_b, clamp.(cosCM,-1,1)) .* abs.(dcosCM_dcosLAB)

    return costanti .* dSigma_dcosLAB, E_d
end

function dn_dtdE(E_d::Vector{Float64}, l_d::Vector{Float64}, m_d, n_b, E_b, l_b, m_b, n_t, m_t, m_r, dSigma_dcosCM::Function)
    A = (E_b .+ m_b .+ m_t).^2 .- E_b .* (E_b .+ 2 .* m_b) .+ m_d^2 .- m_r^2
    B = E_b .+ m_b .+ m_t

    cosLAB = (B .* E_d .+ B .* m_d .- A / 2) ./ sqrt.(E_b .* (E_b .+ 2 .* m_b) .* E_d .* (E_d .+ 2 .* m_d))
    cosX = (-cos.(l_b) .* cos.(l_d) .+ cosLAB) ./ (sin.(l_b) .* sin.(l_d))

    G_b = E_b ./ m_b .+ 1
    v_b = sqrt.((G_b.^2 .- 1) ./ G_b.^2) * cnst.c
    G_d = E_d ./ m_d .+ 1
    v_d = sqrt.((G_d.^2 .- 1) ./ G_d.^2) * cnst.c
    v_CM = sqrt.(E_b .* (E_b .+ 2 .* m_b)) ./ ((G_b .* m_b .+ m_t) / cnst.c)
    G_CM = 1 ./ sqrt.(1 .- v_CM.^2 / cnst.c^2)

    cosCM = (G_CM .* (v_d .* cosLAB .- v_CM)) ./ sqrt.((G_CM.^2 .- 1) .* v_d.^2 .* cosLAB.^2 .- 2 .* G_CM.^2 .* v_d .* v_CM .* cosLAB .+ v_d.^2 .+ G_CM.^2 .* v_CM.^2)

    dEd_dcosLAB = (2 .* sqrt.(E_b .* (E_b .+ 2 .* m_b)) .* (E_d .* (E_d .+ 2 .* m_d)).^1.5) ./ (A .* E_d .+ (A .- 2 .* B .* m_d) .* m_d)
    dvd_dcosLAB = cnst.c .* (m_d .* G_d.^2 .* sqrt.(G_d.^2 .- 1)).^-1 .* dEd_dcosLAB
    dcosCM_dcosLAB = G_CM .* v_d .* ((v_d .* (v_d .- v_CM .* cosLAB) .+ v_CM .* (1 .- cosLAB.^2) .* dvd_dcosLAB) ./ ((G_CM.^2 .- 1) .* v_d.^2 .* cosLAB.^2 .- 2 .* G_CM.^2 .* v_d .* v_CM .* cosLAB .+ v_d.^2 .+ G_CM.^2 .* v_CM.^2).^1.5)

    costanti = n_b .* n_t .* v_b ./ (2*π^2)
    jacobian = abs.( real.(Complex.(1 .- cosX.^2).^-0.5) .* (sin.(l_b) .* sin.(l_d)).^-1 .* abs.(dEd_dcosLAB).^-1 )

    dSigma_dcosLAB = dSigma_dcosCM(E_b, clamp.(cosCM,-1,1)) .* abs.(dcosCM_dcosLAB)

    return costanti .* dSigma_dcosLAB .* jacobian
end

function SetReaction(reaction_name, product_state="GS")
    # Beware: fast and thermal reactants are re-mapped here according to a->t, b->b, c->d, d->r (see miscellaneous comments on top)
    if reaction_name in ["TD_neutrons", "DT_neutrons", "TD_alphas", "DT_alphas"]
        mass_ratio = (reaction_name == "TD_neutrons" || reaction_name == "TD_alphas") ? cnst.mt / (cnst.md + cnst.mt) : cnst.md / (cnst.md + cnst.mt)
        # Read the data file
        df = CSV.read("fit_data/dtn4He_legendre_Drosg.txt", DataFrame; delim='\t', header=false)
        # Convert column data for interpolation
        E_col = Array{Float64}(df[!, 1])
        A_data = Matrix{Float64}(df[!, 3:end])  # Exclude first two columns
        # Interpolation in 1D
        A = [ extrapolate(interpolate((E_col*1000,), A_data[:,i],Gridded(Linear())),0) for i in 1:size(A_data, 2) ]
        # Define differential cross section anonymous function
        lmax = size(A, 1) - 1  # Assuming df.shape[1] - 2 is equivalent to number of columns in `A` matrix
        flip = (reaction_name == "TD_alphas" || reaction_name == "DT_alphas") ? -1 : 1
        # dSigma_dcosCM = (x,y) -> (fsrct.sigma_tot(mass_ratio * x, "d-t") ./ (2 .* A[1](x))) .* sum([A[l+1](x) .* Pl(flip*y,l) for l in 0:lmax])
        dSigma_dcosCM = (x, y) -> (fsrct.sigma_tot(mass_ratio .* x, "d-t") ./ (2 .* A[1].(x))) .* sum(A[l+1].(x) .* Pl.(flip .* y, l) for l in 0:lmax)

        masses = (reaction_name == "TD_neutrons") ? [cnst.md, cnst.mt, cnst.mn,   cnst.m4He] .* cnst.u_keV   :
                 (reaction_name == "DT_neutrons") ? [cnst.mt, cnst.md, cnst.mn,   cnst.m4He] .* cnst.u_keV   :
                 (reaction_name == "TD_alphas") ?   [cnst.md, cnst.mt, cnst.m4He, cnst.mn  ] .* cnst.u_keV   :
                 (reaction_name == "DT_alphas") ?   [cnst.mt, cnst.md, cnst.m4He, cnst.mn  ] .* cnst.u_keV   : error("Invalid reaction name")

        return masses..., dSigma_dcosCM
    # Beware: fast and thermal reactants are re-mapped here according to a->t, b->b, c->d, d->r (see miscellaneous comments on top)    
    elseif reaction_name in ["9Be4He_gammas", "9Be4He_neutrons"]
        if product_state == "1L"
            # Read the data file
            lines = readlines("fit_data/AlphaBe9_1L.txt")
            # 1. Energies of CM in MeV (from the second line of the file)
            energy_line = lines[2]  # Assumes energies are on the second line
            energies = parse.(Float64, split(energy_line, r"\s+"))  # Split by spaces and parse as Float64
            # 2. Cosine values of the scattering angle in the CM frame
            cosines = [parse(Float64, split(line, r"\s+")[1]) for line in lines[3:end] if !isempty(line)]  # First value in each row
            # 3. Differential cross-section values
            diff_cross = [
                parse.(Float64, split(line, r"\s+")[2:end]) for line in lines[3:end] if !isempty(line)
            ]  # Parse remaining columns
            # Define 2D interpolator
            energies = collect(energies * 10^3 * (cnst.m4He + cnst.m9Be) / cnst.m9Be)
            diff_cross = hcat(diff_cross...) * 1e-31
            if reaction_name == "9Be4He_gammas"
                diff_cross = reverse(diff_cross, dims=2)
            end
            dSigma_dcosCM = extrapolate(interpolate((energies, cosines), 2 * π * diff_cross, Gridded(Linear())), 0)
            
            masses = (reaction_name == "9Be4He_gammas") ? [cnst.m4He, cnst.m9Be, cnst.m12C+4439.82/cnst.u_keV, cnst.mn                     ] .* cnst.u_keV : 
                                                          [cnst.m4He, cnst.m9Be, cnst.mn,                      cnst.m12C+4439.82/cnst.u_keV] .* cnst.u_keV

            return masses..., 4438.94, (x,y) -> dSigma_dcosCM.(x,y)
        elseif product_state == "2L"
            # Read the data file
            lines = readlines("fit_data/AlphaBe9_2L.txt")
            # 1. Energies of CM in MeV (from the second line of the file)
            energy_line = lines[2]  # Assumes energies are on the second line
            energies = parse.(Float64, split(energy_line, r"\s+"))  # Split by spaces and parse as Float64
            # 2. Cosine values of the scattering angle in the CM frame
            cosines = [parse(Float64, split(line, r"\s+")[1]) for line in lines[3:end] if !isempty(line)]  # First value in each row
            # 3. Differential cross-section values
            diff_cross = [
                parse.(Float64, split(line, r"\s+")[2:end]) for line in lines[3:end] if !isempty(line)
            ]  # Parse remaining columns
            # Define 2D interpolator
            energies = collect(energies * 10^3 * (cnst.m4He + cnst.m9Be) / cnst.m9Be)
            diff_cross = hcat(diff_cross...) * 1e-31
            if reaction_name == "9Be4He_gammas"
                diff_cross = reverse(diff_cross, dims=2)
            end
            dSigma_dcosCM = extrapolate(interpolate((energies, cosines), 2 * π * diff_cross, Gridded(Linear())), 0)
            
            masses = (reaction_name == "9Be4He_gammas") ? [cnst.m4He, cnst.m9Be, cnst.m12C+7654.07/cnst.u_keV, cnst.mn                     ] .* cnst.u_keV :      
                                                          [cnst.m4He, cnst.m9Be, cnst.mn,                      cnst.m12C+7654.07/cnst.u_keV] .* cnst.u_keV

            return masses..., 3213.79, (x,y) -> dSigma_dcosCM.(x,y)
        end
    else 
        error("Reaction is not a valid fusion reaction available currently")
    end
end

function GetOrbSpec(viewing_cone, o_E, o_p, o_R, o_z, o_w, Ed_bins, o_B; product_state="GS", reaction_name="TD_neutrons", bulk_dens=1.0e19)
    # Intersect calculated orbits with LOS, if needed
    if isnothing(viewing_cone)
        E_b_vc = o_E # No viewing cone specified. All energy points are accepted (/inside the universe)
        p_b_vc = o_p # No viewing cone specified. All pitch points are accepted (/inside the universe)
        R_vc = o_R # No viewing cone specified. All R points are accepted (/inside the universe)
        z_vc = o_z # No viewing cone specified. All z points are accepted (/inside the universe)
        bulk_dens_vc = bulk_dens # No viewing cone specified. All bulk densities are accepted
        weights_vc = o_w # No viewing cone specified. All weights are accepted
        B_vc = o_B # No viewing cone specified. All magnetic field vectors are accepted
        omega = 4*np.pi # Spherical emission, in all directions
        dphi = 2*np.pi # All toroidal angles are included
    else
        i_voxels, i_points = map_points_voxels(viewing_cone, o_R, o_z) # Map the R,z sample points to the viewing cone voxels
        if length(i_points) == 0
            # println("No points inside viewing cone!")
            return zeros(length(Ed_bins) - 1)
        end
        E_b_vc = o_E[i_points] # Extract energy points within diagnostic viewing cone
        p_b_vc = o_p[i_points] # Extract pitch points within diagnostic viewing cone
        R_vc = o_R[i_points] # Extract R points within diagnostic viewing cone
        z_vc = o_z[i_points] # Extract z points within diagnostic viewing cone
        bulk_dens_vc = bulk_dens[i_points]
        weights_vc = o_w[i_points] # Extract weights corresponding to R,z points within diagnostic viewing cone
        B_vc = o_B[:,i_points] # Extract magnetic field vectors corresponding to R,z points within diagnostic viewing cone
        omega = viewing_cone.OMEGA[i_voxels] # Only specific solid angles, in steradian
        dphi = viewing_cone.dPHI # What is the incremental toroidal angle that one diagnostic viewing cone voxel occupies?
        ud_norm = sqrt.(sum(vc.U[:,i_voxels].^2, dims=1))
        ud = vc.U[:,i_voxels] ./ ud_norm # What is the emission direction of the viewing cone?
    end
    
    mag_B_vc = sqrt.(sum(B_vc.^2, dims=1))
    lambda_d = vec( acos.( clamp.( sum(ud .* B_vc, dims=1)./mag_B_vc, -1, 1 ) ) )

    Delta_Ed = Ed_bins[2] - Ed_bins[1]
    Ed_vals = Ed_bins[1:end-1] .+ Delta_Ed / 2
    orbit_spectrum = zeros(length(Ed_vals))

    if product_state != "GS"
        m_b,m_t,m_d,m_r,g_ray,dSigma_dcosCM = SetReaction(reaction_name, product_state)

        grid_lims = get_lims(E_b_vc, acos.(p_b_vc), m_b, m_t, m_d, m_r)
        # γ_b = LinRange(minimum(grid_lims[:,5]), maximum(grid_lims[:,6]), 50)
        # dγ = γ_b[2] - γ_b[1]
        p_d = LinRange(cos(maximum(grid_lims[:,4])), cos(minimum(grid_lims[:,3])), 50)
        dp = p_d[2] - p_d[1]
        
        E_b_vc = E_b_vc[1]
        @inbounds for i in eachindex(p_b_vc)
            γ_b = LinRange(grid_lims[i,5], grid_lims[i,6], 50)
            dγ = γ_b[2] - γ_b[1]
            # p_d = LinRange(cos(grid_lims[i,4]), cos(grid_lims[i,3]), 50)
            # dp = p_d[2] - p_d[1]
        
            p_list = vec( ones(length(γ_b))' .* p_d )
            γ_list = vec( γ_b' .* ones(length(p_d)) )
        
            w_list,e_list = dn_dtdγ(    γ_list, acos.(p_list), m_d, 
                                        10. ^ 0 .* weights_vc[i] , E_b_vc[i], acos.(p_b_vc[i]), m_b, 
                                        bulk_dens_vc[i], m_t, 
                                        m_r, 
                                        dSigma_dcosCM)
            w_list .*= dγ .* dp
            # replace!(w_list, NaN => 0.0)
            # replace!(e_list, NaN => 1.e-6)
            deleteat!(e_list, isnan.(w_list))
            deleteat!(p_list, isnan.(w_list))
            filter!(!isnan, w_list)
        
            @turbo for j in eachindex(Ed_vals)
                orbit_spectrum[j] += sum(Doppler_spectrum( g_ray, Ed_vals[j], lambda_d[i], 
                                                           w_list,  e_list, 
                                                           acos.(p_list), m_d, omega[i]) ) * (dphi / (2π))
            end
        end
            else 
        m_b,m_t,m_d,m_r,dSigma_dcosCM = SetReaction(reaction_name)

        for i in eachindex(Ed_vals)
            orbit_spectrum[i] += sum(dn_dtdE(   Ed_vals[i], lambda_d, m_d, 
                                                10. ^ 0 .* weights_vc .* omega , E_b_vc, acos.(p_b_vc), m_b, 
                                                bulk_dens_vc, m_t, 
                                                m_r, 
                                                dSigma_dcosCM)) * (dphi/(2*π))
        end
    end


    return orbit_spectrum
end

