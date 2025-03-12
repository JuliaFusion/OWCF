########################################## forward.jl ##################################################
# This code is a Julia transposition of forward.py and spec.py, with the key difference that the calculations
# do not account for the thermal reactant dynamics. This is in favor of faster computations with no Monte 
# Carlo sampling needed.
#
# Miscellaneous comments:
#  - the a+b->1+2 notation from DRESS is swapped here with b+t->d+r, meaning that a beam particle 'b' reacts
#    with a target particle 't' and generates a detected particle 'd' and a residual particle 'r'. 
#  - 
#
# TO DO:
#  - check the jacobian for the gamma-ray in spec.py as well. It might not be there
#  - 
#
# THIS SCRIPT IS FUNCTIONALLY DONE. ONLY DOCUMENTATION REMAINS TO BE FINALIZED.
#
# Written by Andrea Valentini. Last maintained 2025-01-20
######################################################################################################

using NaNMath
using LoopVectorization
using DelimitedFiles
using DataFrames
using LegendrePolynomials  # For Legendre polynomials
using LinearAlgebra
using PyCall
include("misc/availReacts.jl") # To be able to handle fusion reaction strings
pushfirst!(PyVector(pyimport("sys")."path"), "") # To add constants and fusreact modules

cnst = pyimport("constants")
fsrct = pyimport("fusreact")

"""
get_lims(E_b::Vector{Float64}, l_b::Vector{Float64}, m_b::Float64, m_t::Float64, m_d::Float64, m_r::Float64)

Method to calculate detected particle energy, pitch-angle limits (first 4 outputs) given a specific beam particle cold ring (E_b,l_b). 
Outputs 5 and 6 are the corresponding beam particle gyro-angle limits (since, for every E_d, one and only one γ_b corresponds). 

This helps because, once the arrays lenghts are fixed, smaller grids provide higher resolution.
The inputs are as follows
- E_b - beam particle energy in keV
- l_b - beam particle pitch-angle in rad
- m_b, m_t, m_d, m_r - masses in keV for beam, target, detected and residual particles
"""
function get_lims(E_b, l_b, m_b, m_t, m_d, m_r)
    cosLAB = 1.0 # grid limit in energy is found at maximum value for the cosine of the emission angle

    # Solving the relativistic 4-momentum conservation equations 
    A = (E_b .+ m_b .+ m_t).^2 .- E_b .* (E_b .+ 2 .* m_b) .+ m_d^2 .- m_r^2
    B = E_b .+ m_b .+ m_t
    D = B .* m_d .- A ./ 2
    a = B.^2 .- E_b .* (E_b .+ 2 .* m_b) .* cosLAB.^2
    bi_half = B .* D .- E_b .* (E_b .+ 2 .* m_b) .* m_d .* cosLAB.^2
    ci = D.^2
    delta_fourth = bi_half.^2 .- a .* ci

    # If the radicand is negative, then no minimum emission angle is defined. We throw zeros here for simplicity, but some tweaks in the future might be needed
    cosLABmin = real.(sqrt.(Complex.((B.^2 .* m_d.^2 .- A.^2 ./ 4) ./ (E_b .* (E_b .+ 2 .* m_b) .* m_d.^2))))

    # Return the lower and upper limits for E_d, l_d and γ_b; the factors are intended to avoid singularities at the domain edges (by some fairly arbitrary 'safe margin')
    return hcat(max.((-bi_half .- 1.1 .* sqrt.(delta_fourth)) ./ a, 0.001), # minimum E_d
                (-bi_half .+ 1.1 .* sqrt.(delta_fourth)) ./ a, # maximum E_d
                max.(l_b .- 1.1 .* acos.(cosLABmin), 0.001), # minimum l_d
                min.(l_b .+ 1.1 .* acos.(cosLABmin), 0.999*π), # maximum l_d
                max.(-π/2 .- 1.15 .* acos.(cosLABmin), -0.999*π), # minimum γ_b corresponding to minimum E_d (highly dependent on coordinate system at use, please note)
                min.(-π/2 .+ 1.15 .* acos.(cosLABmin), -0.001),) # maximum γ_b corresponding to maximum E_d (highly dependent on coordinate system at use, please note)
end

"""
Doppler_spectrum(E_nominal::Float64, E_g::Float64, phi::Vector{Float64}, n_b::Vector{Float64}, E_b::Vector{Float64}, l_b::Vector{Float64}, m_b::Float64, dOmega::Float64)

Method to calculate the relativistic-Doppler-shifted gamma-ray spectrum at energy E_g, given an array of beam particle cold rings (E_b,l_b).
Please note, the beam now is the excited product, eg a carbon-12, and not the same beam particle of the 1st step of the reaction!

The inputs are as follows
- E_nominal - gamma-ray energy in the frame of the moving particle
- E_g - gamma-ray energy at a specific point of the spectrum, corresponding to the middle-point of an energy bin
- phi - emission angle of the gamma-ray with respect to the local magnetic field
- n_b, E_b, l_b, m_b - density in m^-3, energy in keV, pitch-angle in rad and mass in keV of the inputted beam particle cold rings
- dOmega - fraction of solid angle seen by the detector
"""
function Doppler_spectrum(E_nominal, E_g, phi, n_b, E_b, l_b, m_b, dOmega)
    # Calculating the necessary variables
    G_b = @. E_b / m_b + 1
    G_b2 = @. G_b^2
    costanti = @. n_b / (4 * π^2) * dOmega
    v_b = @. cnst.c * sqrt((G_b2 - 1) / G_b2)
    v_d = cnst.c
    cos_l_b = cos.(l_b)
    sin_l_b = sin.(l_b)
    cos_phi = cos.(phi)
    sin_phi = sin.(phi)

    # The sine of the beam particle gyro-angle and the cosine of the emission angle wrt the beam particle are calculated
    sin_gyro = @. (( (E_g/E_nominal)^-1 - G_b + sqrt(G_b2 - 1) * cos_l_b * cos_phi ) / 
                  ( sqrt(G_b2 - 1) * sin_l_b * sin_phi ))
    cosLAB = @. cos_l_b * cos_phi - sin_l_b * sin_phi * sin_gyro

    # The needed Jacobians are calculated, according to ...paperDOI...
    dvd_dcosLAB = 0. # the product, if it's a gamma-ray, does not change its velocity with the emission angle
    dcosCM_dcosLAB = @. G_b * v_d * ((v_d * (v_d - v_b * cosLAB) + v_b * (1 - cosLAB^2) * dvd_dcosLAB) /
                     ((G_b2 - 1) * v_d^2 * cosLAB^2 - 2 * G_b2 * v_d * v_b * cosLAB + 
                      v_d^2 + G_b2 * v_b^2)^1.5)
    
    sin_gyro_sq = @. ifelse(abs(sin_gyro) > 0.99999, -Inf, sin_gyro^2)
    jac = @. inv(sqrt(1 - sin_gyro_sq)) * abs(E_nominal / E_g^2 / sqrt(G_b2 - 1) / sin_l_b / sin_phi)

    # Note that no differential cross section is included, since the emission is isotropic in the excited nucleus rest frame
    return @. costanti * jac * dcosCM_dcosLAB
end

"""
relativistic_inverse_jacobian(γ_b::Vector{Float64}, cosLAB::Vector{Float64}, E_b::Vector{Float64}, m_b::Float64, m_t::Float64, m_d::Float64, m_r::Float64, A::Vector{Float64}, B::Vector{Float64}, D::Vector{Float64})

Method to calculate the relativistic inverse jacobian for the cosLAB->E_d transformation. Obtained by differentiating 4-momentum conservation laws.
Please look at: ...paperDOI... for thorough explanation on this.

The inputs are as follows
- γ_b - beam particle gyro-angle in rad, which in case of beam-target reactions fully specifies E_d as well (but has better resolution for same grid size, hence why it is used)
- cosLAB - cosine of emission angle with respect to the beam particle velocity
- E_b, m_b - energy in keV and mass in keV for beam particle cold ring
- m_t, m_d, m_r - masses in keV for target, detected and residual particles
- A, B, D - pre-computed coefficients in keV^2, keV and keV^2 respectively; they depend on energies and masses of the reactants. See ...paperDOI...
"""
function relativistic_inverse_jacobian(γ_b, cosLAB, E_b, m_b, m_t, m_d, m_r, A, B, D)
    # We have to account for the sign in the quadratic solution, before differentiating
    s = -sign.(cos.(γ_b))

    jac = s .* ( (E_b .* (E_b .+ 2 * m_b) .* cosLAB .* 
            (-2 .* B.^3 .* D .* m_d .+ D.^2 .* E_b .* (E_b .+ 2 * m_b) .* cosLAB.^2 
            .- s .* 2 .* B .* D .* (s .* E_b.^2 .* m_d .* cosLAB.^2 .+ s .* 2 .* E_b .* m_b .* m_d .* cosLAB.^2 
            .+ NaNMath.sqrt.(E_b .* (E_b .+ 2 * m_b) .* cosLAB.^2 .* (D.^2 .- 2 .* B .* D .* m_d .+ E_b .* (E_b .+ 2 * m_b) .* m_d.^2 .* cosLAB.^2)))
            .+ B.^2 .* (D.^2 .+ 2 .* m_d .* (E_b.^2 .* m_d .* cosLAB.^2 .+ 2 .* E_b .* m_b .* m_d .* cosLAB.^2 
            .+ s .* NaNMath.sqrt.(E_b .* (E_b .+ 2 * m_b) .* cosLAB.^2 .* (D.^2 .- 2 .* B .* D .* m_d .+ E_b .* (E_b .+ 2 * m_b) .* m_d.^2 .* cosLAB.^2)))))) 
        ) ./ ((B.^2 .- E_b .* (E_b .+ 2 * m_b) .* cosLAB.^2).^2 .* NaNMath.sqrt.(E_b .* (E_b .+ 2 * m_b) .* cosLAB.^2 .* (D.^2 .- 2 .* B .* D .* m_d .+ E_b .* (E_b .+ 2 * m_b) .* m_d.^2 .* cosLAB.^2)))

    return jac
end

"""
dn_dtdγ(γ_b::Vector{Float64}, l_d::Vector{Float64}, m_d::Float64, n_b::Vector{Float64}, E_b::Vector{Float64}, l_b::Vector{Float64}, m_b::Float64, n_t::Vector{Float64}, m_t::Float64, m_r::Float64, dSigma_dcosCM::Function)

Method to calculate differential rate of excited product emitted per unit beam particle gyro angle dγ_b (subscript omitted for brevity).
This is a new algorithm developed just for an even faster implementation of 2step reaction spectra. Most formulas can be still derived from ...paperDOI...

The inputs are as follows
- γ_b - beam particle gyro-angle in rad, which in case of beam-target reactions fully specifies E_d as well (but has better resolution for same grid size, hence why it is used)
- l_d, m_d - pitch-angle in rad and mass in keV for the detected particle
- n_b, E_b, l_b, m_b - density in m^-3, energy in keV, pitch-angle in rad and mass in keV of the inputted beam particle cold rings
- n_t, m_t - density in m^-3 and mass in keV of the target particles
- m_r - mass in keV of the residual particle
- dSigma_dcosCM - differential cross section in m^2 for specified reaction. NB it is not per unit steradian!
"""
function dn_dtdγ(γ_b, l_d, m_d, n_b, E_b, l_b, m_b, n_t, m_t, m_r, dSigma_dcosCM)
    # Defining particle-dependent coefficients
    A = (E_b .+ m_b .+ m_t).^2 .- E_b .* (E_b .+ 2 .* m_b) .+ m_d^2 .- m_r^2
    B = E_b .+ m_b .+ m_t
    D = B .* m_d .- A ./ 2

    # Solving the quadratic equation (see ...paperDOI...)
    cosLAB = cos.(l_b) .* cos.(l_d) .- sin.(l_b) .* sin.(l_d) .* sin.(γ_b)
    a = B.^2 .- E_b .* (E_b .+ 2 .* m_b) .* cosLAB.^2
    bi_half = B .* D .- E_b .* (E_b .+ 2 .* m_b) .* m_d .* cosLAB.^2
    ci = D.^2
    E_d = (-1 .* bi_half .- sign.(cos.(γ_b)) .* NaNMath.sqrt.(bi_half.^2 .- a .* ci)) ./ a

    # Pre-computing velocities and Lorentz factors for the particles and center-of-mass
    G_b = E_b ./ m_b .+ 1
    v_b = sqrt.((G_b.^2 .- 1) ./ G_b.^2) * cnst.c
    G_d = E_d ./ m_d .+ 1
    v_d = sqrt.((G_d.^2 .- 1) ./ G_d.^2) * cnst.c
    v_CM = sqrt.(E_b .* (E_b .+ 2 .* m_b)) ./ ((G_b .* m_b .+ m_t) / cnst.c)
    G_CM = 1 ./ sqrt.(1 .- v_CM.^2 / cnst.c^2)

    # Calculating emission angle wrt beam particle in the center-of-mass rest frame
    cosCM = (G_CM .* (v_d .* cosLAB .- v_CM)) ./ sqrt.((G_CM.^2 .- 1) .* v_d.^2 .* cosLAB.^2 .- 2 .* G_CM.^2 .* v_d .* v_CM .* cosLAB .+ v_d.^2 .+ G_CM.^2 .* v_CM.^2)

    # Needed jacobians
    dEd_dcosLAB = relativistic_inverse_jacobian(γ_b, cosLAB, E_b, m_b, m_t, m_d, m_r, A, B, D)
    dvd_dcosLAB = cnst.c .* (m_d .* G_d.^2 .* sqrt.(G_d.^2 .- 1)).^-1 .* dEd_dcosLAB
    dcosCM_dcosLAB = G_CM .* v_d .* ((v_d .* (v_d .- v_CM .* cosLAB) .+ v_CM .* (1 .- cosLAB.^2) .* dvd_dcosLAB) ./ ((G_CM.^2 .- 1) .* v_d.^2 .* cosLAB.^2 .- 2 .* G_CM.^2 .* v_d .* v_CM .* cosLAB .+ v_d.^2 .+ G_CM.^2 .* v_CM.^2).^1.5)

    costanti = n_b .* n_t .* v_b ./ π

    dSigma_dcosLAB = dSigma_dcosCM(E_b, clamp.(cosCM,-1,1)) .* abs.(dcosCM_dcosLAB)

    # Return the differential rate AND the corresponding energies of the detected particles (remember γ_b <-> E_d)
    return costanti .* dSigma_dcosLAB, E_d
end

"""
dn_dtdE(E_d::Vector{Float64}, l_d::Vector{Float64}, m_d::Float64, n_b::Vector{Float64}, E_b::Vector{Float64}, l_b::Vector{Float64}, m_b::Float64, n_t::Vector{Float64}, m_t::Float64, m_r::Float64, dSigma_dcosCM::Function)

Method to calculate differential rate of excited product emitted per unit detected particle energy E_d (subscript omitted for brevity).
See https://doi.org/10.1088/1741-4326/ad9bc8 and https://doi.org/10.1063/5.0216680 for thorough explanation on this.

The inputs are as follows
- E_d - energy in keV for the emitted particle 
- l_d, m_d - pitch-angle in rad and mass in keV for the detected particle
- n_b, E_b, l_b, m_b - density in m^-3, energy in keV, pitch-angle in rad and mass in keV of the inputted beam particle cold rings
- n_t, m_t - density in m^-3 and mass in keV of the target particles
- m_r - mass in keV of the residual particle
- dSigma_dcosCM - differential cross section in m^2 for specified reaction. NB it is not per unit steradian!
"""
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

"""
SetReaction(reaction::String)

Set reactants, products and, where needed, gamma-ray nominal energy, according to reaction string with format a(b,c)d-l.

"""
function SetReaction(reaction)
    if !reactionIsAvailableAnalytically(reaction)
        error("Expected spectra from fusion reaction $(reaction) is currently not available for computation via analytic equations. Currently analytically available fusion reactions include: $(OWCF_AVAILABLE_FUSION_REACTIONS_FOR_ANALYTIC_COMPUTATION). Please correct and re-try.")
    end
    reactants = getFusionReactants(reaction)
    products = getFusionProducts(reaction)
    nuclear_state = getEmittedParticleEnergyLevel(reaction)

    # Beware: fast and thermal reactants are re-mapped here according to a->t, b->b, c->d, d->r (see miscellaneous comments on top)
    if lowercase.(reactants) == ["d","t"] || lowercase.(reactants) == ["t","d"]
        mass_ratio = (lowercase(reactants[2]) == "d") ? cnst.mt / (cnst.md + cnst.mt) : cnst.md / (cnst.md + cnst.mt)
        # Read the data file
        df = readdlm("fit_data/dtn4He_legendre_Drosg.txt", '\t', Float64, '\n')  
        # Convert column data for interpolation
        E_col = Array{Float64}(df[!, 1])
        A_data = Matrix{Float64}(df[!, 3:end])  # Exclude first two columns
        # Interpolation in 1D
        A = [ extrapolate(interpolate((E_col*1000,), A_data[:,i],Gridded(Linear())),0) for i in 1:size(A_data, 2) ]
        # Define differential cross section anonymous function
        lmax = size(A, 1) - 1  # Assuming df.shape[1] - 2 is equivalent to number of columns in `A` matrix
        flip = (lowercase(products[1]) == "4he") ? -1 : 1
        # dSigma_dcosCM = (x,y) -> (fsrct.sigma_tot(mass_ratio * x, "d-t") ./ (2 .* A[1](x))) .* sum([A[l+1](x) .* Pl(flip*y,l) for l in 0:lmax])
        dSigma_dcosCM = (x, y) -> (fsrct.sigma_tot(mass_ratio .* x, "d-t") ./ (2 .* A[1].(x))) .* sum(A[l+1].(x) .* Pl.(flip .* y, l) for l in 0:lmax)

        masses = (lowercase.(reactants) == ["t","d"] && lowercase.(products) == ["n","4he"]) ? [cnst.md, cnst.mt, cnst.mn,   cnst.m4He] .* cnst.u_keV   :
                 (lowercase.(reactants) == ["d","t"] && lowercase.(products) == ["n","4he"]) ? [cnst.mt, cnst.md, cnst.mn,   cnst.m4He] .* cnst.u_keV   :
                 (lowercase.(reactants) == ["t","d"] && lowercase.(products) == ["4he","n"]) ?   [cnst.md, cnst.mt, cnst.m4He, cnst.mn  ] .* cnst.u_keV   :
                 (lowercase.(reactants) == ["d","t"] && lowercase.(products) == ["4he","n"]) ?   [cnst.mt, cnst.md, cnst.m4He, cnst.mn  ] .* cnst.u_keV   : error("Invalid reaction name")

        return masses..., dSigma_dcosCM
    # Beware: fast and thermal reactants are re-mapped here according to a->t, b->b, c->d, d->r (see miscellaneous comments on top)    
    elseif lowercase.(reactants) == ["4he","9be"] || lowercase.(reactants) == ["9be","4he"]
        if nuclear_state == "1L"
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
            if lowercase(products[1]) == "12c"
                diff_cross = reverse(diff_cross, dims=2)
            end
            dSigma_dcosCM = extrapolate(interpolate((energies, cosines), 2 * π * diff_cross, Gridded(Linear())), 0)
            
            masses = (lowercase(products[1]) == "12c") ?    [cnst.m4He, cnst.m9Be, cnst.m12C+4439.82/cnst.u_keV, cnst.mn                     ] .* cnst.u_keV : 
                                                            [cnst.m4He, cnst.m9Be, cnst.mn,                      cnst.m12C+4439.82/cnst.u_keV] .* cnst.u_keV

            return masses..., 4438.94, (x,y) -> dSigma_dcosCM.(x,y)
        elseif nuclear_state == "2L"
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
            if lowercase(products[1]) == "12c"
                diff_cross = reverse(diff_cross, dims=2)
            end
            dSigma_dcosCM = extrapolate(interpolate((energies, cosines), 2 * π * diff_cross, Gridded(Linear())), 0)
            
            masses = (lowercase(products[1]) == "12c") ?    [cnst.m4He, cnst.m9Be, cnst.m12C+7654.07/cnst.u_keV, cnst.mn                     ] .* cnst.u_keV :      
                                                            [cnst.m4He, cnst.m9Be, cnst.mn,                      cnst.m12C+7654.07/cnst.u_keV] .* cnst.u_keV

            return masses..., 3213.79, (x,y) -> dSigma_dcosCM.(x,y)
        else
            error("The $(reaction) fusion reaction needs to have the energy state 'l' (a(b,c)d-l) in either 1L or 2L. Please correct and re-try.")
        end
    else 
        error("The $(reaction) fusion reaction is currently not available for analytic computation of expected diagnostic spectrum. The following fusion reaction are available: $(OWCF_AVAILABLE_FUSION_REACTIONS_FOR_ANALYTIC_COMPUTATION). Please correct and re-try.")
    end
end

"""
forward_calc(viewing_cone::ViewingCone, o_E::Vector{Float64}, o_p::Vector{Float64}, o_R::Vector{Float64}, o_z::Vector{Float64}, o_w::Vector{Float64}, Ed_bins::Vector{Float64}, o_B::Matrix{Float64}; 
            reaction::String="T(D,n)4He-GS", bulk_dens::Union{Nothing,Vector{Float64}}=1.0e19)

For the input (E,p,R,z) points and weights (please see inputs explanation below), calculate expected diagnostic spectrum for 1- and 2-step fusion reactions using the analytic equations obtained when making the 
beam-target approximation. Please see [CITE A. VALENTINI PUBLICATION]. Since the approach is straight forward (no Monte Carlo computations needed), a massive speed-up of computations compared to established 
Monte Carlo codes (e.g. DRESS/GENESIS) is achieved. However, the approximation is only valid as long as the thermal (slow) fusion product particle population can be assumed to be at rest (zero temperature),
relative to the energetic (fast) fusion product particle population.

The inputs are as follows
- viewing_cone - A viewing cone object as defined in vcone.jl
- o_E, o_p, o_R, o_z, o_w - Weighted (E,p,R,z) points; fast-ion energies (o_E) in keV, fast-ion pitch (o_p) values, major radius (R) positions in meters, vertical positions (z) in meters, and geometrical weights (w).
- Ed_bins - Diagnostic measurement bin edges for product energy spectrum as defined by the user. Units in keVs
- o_B - A (3,N) array containing the (R,phi,z) components of the magnetic field vector at all (E,p,R,z) points (in Teslas)
- reaction - Fusion reaction. Reaction format 1 and 2 is supported. Please see OWCF/misc/availReacts.jl for explanation.
- bulk_dens - Number densities in m^-3 for target (thermal) particle species at the (E,p,R,z) points.
"""
function forward_calc(viewing_cone, o_E, o_p, o_R, o_z, o_w, Ed_bins, o_B; reaction="T(D,n)4He-GS", bulk_dens=ones(length(o_E)).*1.0e19) 
    
    if getReactionForm(reaction)==3
        error("Fusion reaction specified as $(reaction). Spectrum computation via projected velocities assumed. This is not compatible with the analytic forward model. Please use fusion reaction form 1 or 2. Please see OWCF/misc/availReacts.jl for explanation. Please correct and re-try.")
    end

    m_b,m_t,m_d,m_r,g_ray,dSigma_dcosCM = SetReaction(reaction) # Also checks whether reaction is analytically available in the OWCF

    # Keep only (E,p,R,z) points inside diagnostic line-of-sight (LOS)
    i_voxels, i_points = map_points_voxels(viewing_cone, o_R, o_z) # Map the R,z sample points to the viewing cone voxels. See OWCF/vcone.jl
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
    
    mag_B_vc = sqrt.(sum(B_vc.^2, dims=1))
    lambda_d = vec( acos.( clamp.( sum(ud .* B_vc, dims=1)./mag_B_vc, -1, 1 ) ) )

    Delta_Ed = Ed_bins[2] - Ed_bins[1]
    Ed_vals = Ed_bins[1:end-1] .+ Delta_Ed / 2
    orbit_spectrum = zeros(length(Ed_vals))

    if !(getEmittedParticleEnergyLevel(reaction)=="GS")
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

