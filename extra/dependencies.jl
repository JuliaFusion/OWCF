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
# Finally, please note that some functions in dependencies.jl might be under construction. Hence, the bad code.
#
# Written by H. Järleblad. Last maintained 2025-07-31.
###################################################################################################

println("Loading the Julia packages for the OWCF dependencies... ")
using Base.Iterators
using Contour
using Distributed
using Distributions
using EFIT # For calculating magn½etic equilibrium quantities
using Equilibrium # For loading flux function data, tokamak geometry data etc.
using FileIO # To write/open files in general
using ForwardDiff
using GuidingCenterOrbits # For computing guiding-center orbit trajectories
using HDF5
using Interpolations # To be able to interpolate
using JLD2 # To write/open .jld2 files (Julia data files, basically)
using LinearAlgebra
using NearestNeighbors
using NetCDF # To enable write/open .cdf files
using OrbitTomography # This is what all this is about!
using ProgressMeter # To display computational progress during parallel computations
using SparseArrays # To enable utilization of sparse matrices/vectors
using VoronoiDelaunay
include("../misc/convert_units.jl") # Some functions in dependencies.jl need units functions
include("../misc/species_func.jl") # Some functions in dependencies.jl need species functions (including OWCF/extra/constants.jl)

###### Structures needed for dependencies ######

"""
    OWCF_grid(r2d,z2d,r,z,phi,nr,nz,nphi)

A structure to represent an (R,z) grid in 2D, together with a 
grid in the toroidal phi direction. The r2d and z2d fields are mesh 
grids of the (R,z) grid points.
"""
struct OWCF_grid
    r2d::Matrix{Float64}
    z2d::Matrix{Float64}
    r::Vector{Float64}
    z::Vector{Float64}
    phi::Vector{Float64}
    nr::Int64
    nz::Int64
    nphi::Int64
end

###### Functions ######

###### Mathematics

"""
    gaussian_2D(x_grid, y_grid, x_mean, y_mean, x_std, y_std)
    gaussian_2D(-||-; floor_level=0.0, verbose=false)
"""
function gaussian_2D(x_grid::AbstractVector, y_grid::AbstractVector, x_mean::T, y_mean::T, x_std::T, y_std::T; floor_level=0.0, verbose=false) where T<:Real
    nx = length(x_grid)
    ny = length(y_grid)

    g = zeros(nx,ny)
    g_max = inv(2*pi*x_std*y_std)
    for i=1:nx, j=1:ny
        a = ((x_grid[i] - x_mean)/x_std)^2 + ((y_grid[j] - y_mean)/y_std)^2
        g[i,j] = g_max *exp(-0.5*a)
    end

    if floor_level>0.0
        g = map(gg-> gg<floor_level*g_max ? 0.0 : gg, g)
    end

    return g
end

"""
    gaussian(μ, σ)
    gaussian(-||-; mx=μ .+6 .*σ, mn=μ .-6 .*σ, n=50)

Compute the (multi-variate) Gaussian distribution with mean 'μ' and standard deviation 'σ'.
By default, the upper bounds 'mx' of the grid are found by adding 6σ to the mean μ.
The lower bounds 'mn' of the grid are found by subtracting 6σ from the mean μ. By default, the number 
of grid points in all dimensions is the same: 50. To use the function to create a D-dimensional Gaussian
distribution, where D is any integer >0, do e.g. the following. The example below is for a 3-dimensional case:
μ = [100.0, 0.6, 3.3]
σ = [25.0, 0.01, 0.04]
myGauss = gaussian(μ, σ; mx=[200.0, 1.0, 3.8], mn=[0.0, -1.0, 3.0], n=[10,101,104])

The variable 'myGauss' will then be a 10x101x104 array. In the first dimension, the lower and upper bounds will 
be 0.0 and 200.0, respectively. And so on for the other dimensions.

The floor_level keyword argument can be used to set all values smaller than floor_level*maximum(f)
to 0.0 before returning the output. f is the Gaussian distribution.

The verbose keyword argument will make the function talk more!
"""
function gaussian(μ::AbstractVector, σ::AbstractVector; mn::AbstractVector=μ .-6 .*σ, mx::AbstractVector=μ .+6 .*σ, n::Union{Int64,Vector{Int64}}=50, floor_level::Float64=0.0, verbose::Bool=false)
    DIM=length(μ) # The number of dimensions
    if !(DIM==length(σ))
        error("length(μ)=$(DIM) while length(σ)=$(length(σ)). The number of mean (μ) points must be equal to the number of standard deviation (σ) points. Please correct and re-try.")
    end
    if !(DIM==length(mx))
        error("length(μ)=$(DIM) while length(mx)=$(length(mx)). The number of upper bound (mx) points must equal the number of mean (μ) points. Please correct and re-try.")
    end
    if !(DIM==length(mn))
        error("length(μ)=$(DIM) while length(mn)=$(length(mx)). The number of lower bound (mn) points must equal the number of mean (μ) points. Please correct and re-try.")
    end
    if !(DIM==length(n))
        verbose && println("Matching length of grid points 'n' to $(DIM)-dimensional space... ")
        n = repeat(vcat(n),DIM)
    end
    verbose && println("Upper bound(s) for $(DIM)-dimensional grid: $(mx)")
    verbose && println("Lower bound(s) for $(DIM)-dimensional grid: $(mn)")
    verbose && println("Number of grid points (in each dimension): $(n)")

    if DIM==2
        D1_vector = collect(range(mn[1], stop=mx[1], length=n[1]))
        D2_vector = collect(range(mn[2], stop=mx[2], length=n[2]))
        return gaussian_2D(D1_vector, D2_vector, μ[1], μ[2], σ[1], σ[2]; floor_level=floor_level, verbose=verbose)
    end

    v = σ.^2 # Compute the variance from the standard deviation
    verbose && println("Creating $(DIM)-dimensional grid for Gaussian distribution... ")
    query_vecs_n_inds = () # A tuple to hold all query points and their indices. Structure: ((vector,indices),(vector,indices),...)
    for i in 1:DIM # For all grid dimensions... 
        query_vecs_n_inds = tuple(query_vecs_n_inds[:]...,collect(zip(collect(range(mn[i],stop=mx[i],length=n[i])),1:n[i]))) # Add the (vector,indices) pairs one by one  into a big tuple (tuples are immutable, hence the cumbersome code)
    end
    query_points_n_coords = Iterators.product(query_vecs_n_inds...) # Create a long list of all coordinate space grid points and their coordinates by computing a product between all query point-index vectors. Example structure (if 3 dimensions): [((x1_1,1),(x2_1,1),(x3_1,1)),((x1_2,2),(x2_1,1),(x3_1,1)),...]
    verbose && print("Computing Gaussian distribution with mean $(μ) and standard deviation $(σ)...")
    gauss_distr = zeros(tuple(n...)) # Pre-allocate Gaussian distribution
    for query_point_n_coord in query_points_n_coords
        point = [p for p in map(x-> x[1],query_point_n_coord)] # The point to compute the Gaussian at. E.g. (100.0,0.3) in energy (keV),pitch
        coord = map(x-> x[2],query_point_n_coord) # The coordinate of that point. E.g. (53,14)
        gauss_distr[coord...] = ((2*pi)^(-DIM/2))*inv(sqrt(det(diagm(v)))) *exp(-0.5*transpose(point - μ)*inv(diagm(v))*(point - μ))
    end
    verbose && println("Done!")

    if floor_level>0
        verbose && print("Grid points with values below $(floor_level)*maximum(gauss_distr) will be manually set to 0.0...")
        max_g = maximum(gauss_distr)
        gauss_distr = map(x-> x<floor_level*max_g ? 0.0 : x, gauss_distr)
        verbose && println("Done!")
    end
    return gauss_distr
end

"""
erf(x::Real)
erf(x; resolution::Int64 = 1000, sigma=1/sqrt(2))

This is the error function, defined as 

erf(x) = (2/(sqrt(2*π)*σ)) ∫ exp(-t^2 / (2 σ^2)) dt

where the lower integration limit is 0 and the upper integration limit is x. σ is set to 1/sqrt(2) by default, via 
the keyword argument 'sigma'. This function is a quick approximation, since a sum is used instead of integration.
The number of summation elements can be set via the 'resolution' keyword argument.
"""
function erf(x::Real; resolution::Int64=1000, sigma::Float64=1/sqrt(2))
    if x<0.0
        t_array = collect(range(x,stop=0.0,length=resolution))
    else
        t_array = collect(range(0.0,stop=x,length=resolution))
    end
    dt = x/resolution
    return clamp((2/(sqrt(2*pi)*sigma))*sum(dt .*exp.((-1/(2*sigma^2)) .*(t_array).^2)), -1, 1) # Clamp polishes insufficient resolution issues
end

"""
erfc(x::Real)
erfc(x; resolution::Int64 = 1000)

This is the complementary error function, defined as 

erfc(x) = 1 - erf(x)

This function is a quick approximation, since a sum is used instead of integration for erf(x).
The number of summation elements can be set via the 'resolution' keyword argument.
"""
function erfc(x::Real; resolution::Int64=1000, sigma::Float64=1/sqrt(2))
    return 1 - erf(x; resolution=resolution, sigma=sigma)
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
    forward_difference_matrix(space_size::Tuple, dim::Int64)
    forward_difference_matrix(-||-; verbose=false)

Compute the forward difference matrix (l1) for a coordinate space of size 'space_size' along dimension 'dim'. 
For example, if space_size=(10,30) and dim=2, then 'finite_difference_matrix(space_size,dim)' will output 
a 300x300 matrix Δ with value '-1' for all on-diagonal elements and '+1' for specific off-diagonal elements.
All other elements in Δ will be '0'. For a specific row of the Δ matrix, the element with '+1' value 
correspond to the grid point of the (10,30) grid that have ONE INDEX GREATER (forward difference) in the 
2nd dimension (since dim=2) of the (10,30) grid compared to the grid point to which the on-diagonal element 
correpond. The element order of the 300 rows and columns correspond to the order output of the function 
'CartesianIndices(space_size)'. 1<=dim<=length(space_size) must hold. Output matrix will be in sparse 
matrix format, i.e. SparseMatrixCSC (see SparseArrays.jl package).
"""
function forward_difference_matrix(space_size::Tuple, dim::Int64)
    if !(dim<=length(space_size)) || dim<1
        error("Invalid input. 1<=dim<=length(space_size) must hold. Got length(space_size)=$(length(space_size)) and dim=$(dim).")
    end

    space_indices = CartesianIndices(space_size)

    l1 = spzeros(length(space_indices),length(space_indices)) # The finite difference matrix for a specific dimension (dim)
    for (i,space_coord_i) in enumerate(space_indices) # For every space point (every row in the finite difference matrix)
        for (j,space_coord_j) in enumerate(space_indices) # -||- (for every column in the finite difference matrix)
            if sum(Tuple(space_coord_j) .- Tuple(space_coord_i))==1 && (space_coord_j[dim]-space_coord_i[dim])==1 # If the coordinates differ by one (neighbours with j coordinate larger than i coordinate) and the difference is in the right dimension (dim)
                l1[i,i] = -1 # Set value '-1' for the on-diagonal element
                l1[i,j] = 1 # Set value '+1' for the j off-diagonal element
                break # Only one such nearest neighbour for every row, so move on to the next row
            end
        end
    end

    return l1
end

"""
    forward_difference_matrix(space_size::Tuple)
    forward_difference_matrix(-||-; verbose=false)

Compute the forward difference matrix (L1) for a coordinate space of size 'space_size'. The size of
the L1 matrix will be (length(space_size) * reduce(*,space_size), reduce(*,space_size)). The element 
order of the columns corresponds to the order output of the function 'CartesianIndices(space_size)'.
Output matrix will be in sparse matrix format, i.e. SparseMatrixCSC (see SparseArrays.jl package).
Please see function 'forward_difference_matrix(space_size,dim)' above for detailed info.
The keyword arguments are:
    - verbose - If set to true, the function will be talkative! - Bool
"""
function forward_difference_matrix(space_size::Tuple; verbose=false)
    L1 = Vector{SparseMatrixCSC{Float64,Int64}}(undef,length(space_size)) # All finite difference matrices (one for each dimension) (A Vector of sparce matrices (SparseMatrixCSC) CSC means 'compressed sparse column')
    for idim=1:length(space_size) # For each dimension
        verbose && println("Computing forward difference matrix for dimension $(idim) of $(length(space_size))... ")
        L1[idim] = forward_difference_matrix(space_size, idim) # The forward difference matrix for this particular dimension
    end
    
    return sparse_vcat(L1...) # Concatenate vertically and return
end

###### Geometry

"""
rz_grid(rmin, rmax, nr, zmin, zmax, nz, phimin=0.0, phimax=0.0, nphi=1)

Creates an interpolation grid.

#### Input arguments
## rmin - Minimum radius [cm]
## rmax - Maximum radius [cm]
## nr - Number of radii
## zmin - Minimum z value [cm]
## zmax - Maximum z value [cm]
## nz - Number of z values
## phimin - Minimum Phi value [rad]
## phimax - Maximum Phi value [rad]
## nphi - Number of Phi values 

#### Return Value
## Interpolation grid dictionary

####Example Usage
## ```julia
## julia> grid = rz_grid(0,200.0,200,-100,100,200,phimin=4*np.pi/3,phimax=5*np.pi/3,nphi=5)
## ```
"""
function rz_grid(rmin::Union{Int64,Float64}, rmax::Union{Int64,Float64}, nr::Int64, zmin::Union{Int64,Float64}, zmax::Union{Int64,Float64}, nz::Int64; phimin::Union{Int64,Float64}=0.0, phimax::Union{Int64,Float64}=2*pi-eps(), nphi::Int64=360)
    dr = (rmax - rmin) / nr
    dz = (zmax - zmin) / nz
    dphi = (phimax - phimin) / nphi
    r = rmin .+ dr * (0:(nr-1))
    z = zmin .+ dz * (0:(nz-1))
    phi = phimin .+ dphi * (0:(nphi-1))

    r2d = repeat(r,1,nz)
    z2d = repeat(z,1,nr)'

    return OWCF_grid(r2d,z2d,r,z,phi,nr,nz,nphi)
end

###### Physics

"""
    debye_length(n_e::Real, T_e::Real, species_th_vec::Vector{String}, n_th:vec::Vector{T}, T_th_vec::Vector{T}) where T<:Real
    debye_length(-||-; species_2::String=species_1, n_2::Real=n_1, T_2::Real=T_1)

Compute the (exact) plasma debye length. The inputs are:
    - n_e: The electron density [m^⁻3]
    - T_e: The electron temperature [keV]
    - species_th_vec: A vector containing all thermal species string identifiers, e.g. ["D", "T", "3he"] [-]
    - n_th_vec: A vector containing all thermal species densities. n_th_vec[i] is the density of species_th_vec[i] [m^-3]
    - T_th_vec: A vector containing all thermal species temperatures. T_th_vec[i] is the temperature of species_th_vec[i] [keV]
The outputs are:
    - debye_length: The debye length [m]
"""
function debye_length(n_e::Float64, T_e::Real, species_th_vec::Vector{String}, n_th_vec::Vector{T} where {T<:Real}, T_th_vec::Vector{T} where {T<:Real})
    temp_e = T_e*1000*(GuidingCenterOrbits.e0)/OWCF_kB
    temp_th_vec = (1000*(GuidingCenterOrbits.e0)/OWCF_kB) .*T_th_vec
    Z_th_vec = getSpeciesEcu.(species_th_vec)
    denom = (n_e/temp_e) + reduce(+, Z_th_vec .*Z_th_vec .*n_th_vec ./temp_th_vec)
    return sqrt((OWCF_ϵ0*OWCF_kB/(GuidingCenterOrbits.e0)^2)/denom)
end

"""
    gyro_radius(M::AbstractEquilibrium,p::GCParticle)

Compute the (relativistic) gyro radius for guiding-center particle p, given the magnetic 
equilibrium M. Output in meters. Take relativity into account. The inputs are:
    - M: An abstract equilibrium. Either from an .eqdsk file or an output file from OWCF/extra/compSolovev.jl [-]
    - p: The guiding-center particle object for the particle species of interest, e.g. GCDeuteron [-]
The output is:
    - gyro_radius: The relativistically correct gyro radius [m]
"""
function gyro_radius(M::AbstractEquilibrium,p::GCParticle)
    γ = GuidingCenterOrbits.lorentz_factor(p)
    Babs = norm(Bfield(M, p.r, p.z)) # Magnetic field magnitude. Tesla
    m = p.m # Mass of particle
    KE = (γ-1)*p.energy # Relativistic kinetic energy. keV
    mc2 = m*GuidingCenterOrbits.c0*GuidingCenterOrbits.c0 # Rest energy. Joule
    KE_j = GuidingCenterOrbits.e0*KE*1e3 # Kinetic energy. Joule
    p_rel2 = ((KE_j + mc2)^2 - mc2^2)/(GuidingCenterOrbits.c0*GuidingCenterOrbits.c0) # Relativistic momentum
    p_perp2 = p_rel2*(1-p.pitch^2) # Square of relativistic perpendicular momentum

    return sqrt(p_perp2) / (abs(GuidingCenterOrbits.e0*p.q)*Babs*γ)
end

"""
    spitzer_slowdown_time(n_e, T_e, species_f, species_th_vec, n_th_vec, T_th_vec)
    spitzer_slowdown_time(-||-; plasma_model = :texas, returnExtra = false)

Compute the non-relativistic Spitzer slowing-down time (in seconds) for fast-ion species 'species_f', following the equation in the ITER Physics Basis (http://sites.apam.columbia.edu/fusion/IPB_Chap_5.pdf). 
Assume multiple thermal species via the vector inputs. By default, use the texas (University of Texas) model for the Coloumb logarithm.
if returnExtra, in addition to τ_s, return the coulomb logarithm as well as the Debye length. The inputs are: 
    - n_e: The electron density [m^-3]
    - T_e: The electron temperature [keV]
    - species_f: The beam injection particle species, e.g. "D", "T" etc [-]
    - species_th_vec: A vector containing all thermal species string identifiers, e.g. ["D", "T", "3he"] [-]
    - n_th_vec: A vector containing all thermal species densities. n_th_vec[i] is the density of species_th_vec[i] [m^-3]
    - T_th_vec: A vector containing all thermal species temperatures. T_th_vec[i] is the temperature of species_th_vec[i] [keV]
The keyword arguments are:
    - plasma_model: The model to use for the plasma parameter. Currently supported :salewski or :texas [-]
    - returnExtra: If true, in addition to the Spitzer slowing-down time, the Coulomb logarithm and Debye length will be returned as well [-]
"""
function spitzer_slowdown_time(n_e::Real, T_e::Real, species_f::String, species_th_vec::Union{String,Vector{String}}, n_th_vec::Union{T,Vector{T}} where {T<:Real}, T_th_vec::Union{T,Vector{T}} where {T<:Real}; plasma_model::Symbol = :texas, returnExtra::Bool = false)
    species_th_vec = vcat(species_th_vec) # In case input was a String, and not a Vector
    n_th_vec = vcat(n_th_vec) # -||-
    T_th_vec = vcat(T_th_vec) # -||-

    m_f = getSpeciesMass(species_f) # The fast-ion species mass, kg
    m_e = (GuidingCenterOrbits.e_amu)*(GuidingCenterOrbits.mass_u) # Electron mass, kg

    λ_D = debye_length(n_e, T_e, species_th_vec, n_th_vec, T_th_vec)

    if plasma_model==:salewski
        Λ = 6*pi*n_e*λ_D^3 # From M. Salewski, A.H. Nielsen, Plasma Physics: lectures notes, 2021.
        Λ_c = Λ # Actually, this is Λ_c ≈ Λ, where Λ is the plasma parameter.
    elseif plasma_model==:texas
        Z_th_vec = getSpeciesEcu.(species_th_vec) # Atomic charge number for all thermal plasma species
        q_th_vec = (GuidingCenterOrbits.e0) .*Z_th_vec # Charge (in Coulombs) for all thermal plasma species
        q2_avg = reduce(+,map(x-> x[1]*x[2],collect(Iterators.product(q_th_vec,q_th_vec))[:]))/(length(q_th_vec)^2) # The average of q_i*q_j for all species pairs (i,j)
        r_closest = q2_avg/(4*pi*OWCF_ϵ0*T_e*GuidingCenterOrbits.e0*1000) # r_closest = <q_i*q_j>/(4*pi*ϵ0*T). Mean value of closest approach, University of Texas. Assume same termperature.
        Λ_c = λ_D/r_closest
    else
        error("Currently supported models for the plasma parameter are :salewski and :texas. Please correct and re-try.")
    end
    coulomb_log = log(Λ_c)
    A_D = (n_e*((GuidingCenterOrbits.e0)^4)*coulomb_log)/(2*pi*(OWCF_ϵ0^2)*(m_f^2))
    τ_s = (3*sqrt(2*pi)*((T_e*GuidingCenterOrbits.e0*1000)^(3/2)))/(sqrt(m_e)*m_f*A_D)
    if returnExtra
        return τ_s, coulomb_log, λ_D 
    end
    return τ_s
end

"""
    gaussianBasisMatrix()
"""
function gaussianBasisMatrix(abscissas, corr_lengths)
    DIM = length(abscissas) # The dimensionality of the space, e.g. 3
    SIZE = length.(abscissas) # The size of the space, e.g. [30, 31, 32]
    MAXIMA = maximum.(abscissas) # The upper boundary in each dimension
    MINIMA = minimum.(abscissas) # The lower boundary in each dimension

    means_in_each_dimension = () # A tuple to hold the gaussian means, in each dimension
    for i in 1:DIM # For all grid dimensions... 
        means_in_each_dimension = tuple(means_in_each_dimension[:]...,collect(range(MINIMA[i], stop=MAXIMA[i], step=corr_lengths[i])))
    end
    means = Iterators.product(means_in_each_dimension...)

    basis_matrix = zeros(reduce(*, SIZE), length(means))
    for (i_col, mean) in enumerate(means)
        basis_matrix[:,i_col] = gaussian([m for m in mean], corr_lengths; mx=MAXIMA, mn=MINIMA, n=SIZE)[:]
    end

    return basis_matrix
end

"""
    basis_prior()
"""
function basis_prior(abscissas; corr_lengths=ones(length(abscissas)), w_min=0.0, type=:gaussian, verbose=false)

    if type==:gaussian
        basis_matrix = gaussianBasisMatrix(abscissas, corr_lengths)
    else
        error("Currently supported options for 'type' keyword argument: (1) :gaussian. Please correct and re-try")
    end

    w = rand(size(basis_matrix,2)) # Initialize random weight in (0,1) for every basis function
    w[findall(x-> x<w_min, w)] .= 0.0 # Set all weights below w_min to 0.0

    return reshape(basis_matrix*w, Tuple(length.(abscissas)))
end

"""
    icrf_streamlines(M, abscissas, FI_species, wave_frequency, cyclotron_harmonic)
    icrf_streamlines(-||-; toroidal_mode_number=1, coord_system="(E,p)", coord_system_order=Tuple(1:length(abscissas)), R_of_interest=R_axis, 
                           z_of_interest=z_axis, verbose=false)

Compute the phase-space streamlines that indicate the direction of (fast-ion) particle transport during electromagnetic wave heating in the ion cyclotron range 
of frequencies (ICRF). Currently, the available phase-space options include: 
    - (E,p), or (p,E), where E is the fast-ion energy in keV, and p is the pitch (vpara/v with v and vpara being the total speed and the speed component 
      parallel to the magnetic field, respectively)
    - (vpara,vperp), or (vperp,vpara), where vpara and vperp are the speed components parallel and perpendicular to the magnetic field, respectively
    - (E,mu,Pphi), or any permutations thereof, where E is the fast-ion energy in keV, mu is the magnetic moment in A*m^2 and Pphi is the toroidal canonical 
      angular momentum in kg*m^2*s^-1
    - (E,mu,Pphi,sigma), or any permutations thereof, where -||- and sigma is a binary coordinate (-1,+1)
    - (E,Lambda,Pphi_n), or any permutations thereof, where -||-, Lambda is the normalized magnetic moment (Lambda=mu*B0/E with B0 the magnetic field on-axis 
      in Teslas and E is the fast-ion energy in Joules) and Pphi_n is the normalized toroidal canonical angular momentum (Pphi_n=Pphi/(q*Ψ_w) with q being the 
      particle charge in Coulomb and Ψ_w the magnetic flux at the separatrix)
    - (E,Lambda,Pphi_n,sigma), or any permutations thereof, where -||- and -||-

The information regarding which one of these phase-space coordinate spaces is discretized via grid points, included as the 'abscissas' input, must be specified 
using the 'coord_system' keyword argument. An explanation of the input variables are as follows:
    - M - The magnetic equilibrium object, as given by the Equilibrium.jl package - AbstractEquilibrium
    - abscissas - The grid points of the abscissas of the phase space of interest - Vector{Vector}
    - FI_species - The (fast-ion) particle of interest. Please see OWCF/misc/species_func.jl for a list of available particle species - String
    - wave_frequency - The wave frequency of the ICRF wave. In Hz - Float64
    - cyclotron_harmonic - The cyclotron harmonic integer of the ICRF heating scheme, for the 'FI_species' particle species and a cyclotron frequency on-axis - Int64
The keyword arguments are as follows:
    - toroidal_mode_number - The toroidal mode number of the ICRF wave. Defaults to 1 - Int64
    - coord_system - The coordinate system of the phase space of interest (please see explanation above) - String
    - coord_system_order - If phase-space abscissas in 'abscissas' input are not given in the common order, e.g. (E,Pphi,mu) instead of (E,mu,Pphi), this need to be 
                           specified via the 'coord_system_order' keyword argument as (1,3,2) instead of the common (1,2,3) order etc. Defaults to (1:length(abscissas)) - Tuple
    - R_of_interest - If "(E,p)" or "(vpara,vperp)" are given as input to the 'coord_system' keyword argument, an (R,z) point of interest needs to be specified. In meters. - Float64
    - z_of_interest - If "(E,p)" or "(vpara,vperp)" are given as input to -||-
    - verbose - If set to true, the function will talk a lot! - Bool

The streamlines will be returned as an Array with dimensionality equal to length(abscissas). Each element of this Array will be a length(abscissas)-dimensional 
unit vector, corresponding to the tangent of the icrf streamlines for that grid point.

PLEASE NOTE! Energy grid points are assumed to be given in keV for all coordinate systems!!!
Velocity grid points are assumed to be given in m/s for all coordinate systems!!!

The icrf_streamlines() function is based on the content in M. Rud et al 2025 Nucl. Fusion 65 056008.
"""
function icrf_streamlines(M::AbstractEquilibrium, abscissas::Vector{Vector{T}} where {T<:Real}, FI_species::String, wave_frequency::Float64, cyclotron_harmonic::Int64;
                          toroidal_mode_number::Int64=1, coord_system::String="(E,p)", coord_system_order=Tuple(1:length(abscissas)), 
                          R_of_interest=magnetic_axis(M)[1], z_of_interest=magnetic_axis(M)[2], verbose::Bool=false)
    coord_system_lc = lowercase(coord_system)
    coordinates = split(coord_system_lc,",")
    coordinates[1] = coordinates[1][2:end] # Remove "("
    coordinates[end] = coordinates[end][1:end-1] # Remove ")"
    
    # Necessary space quantities
    space_dimensionality = length(abscissas) # The dimensionality of the space e.g. 3
    space_size = Tuple(length.(abscissas)) # The size of the space, e.g. (33,52,12). 
    space_indices = CartesianIndices(space_size) # All space indices. E.g. [(1,1,1), (1,1,2), ..., (N1,N2,N3)] where Ni is the number of grid points in the i:th space dimension

    # ICRF wave characteristics independent of phase space and dimensionality
    B_vec_on_axis = Equilibrium.Bfield(M,magnetic_axis(M)...)
    ω     = 2*pi*wave_frequency # The angular frequency of the wave
    ω_0   = (GuidingCenterOrbits.e0)*norm(B_vec_on_axis)/getSpeciesMass(FI_species) # The gyro-motion angular frequency on-axis
    Λ_∞ = cyclotron_harmonic*ω_0/ω # The Λ_∞ parameter (see e.g. M. Rud et al, Nucl. Fusion 2025)

    verbose && println("icrf_streamlines(): Λ_∞: $(Λ_∞)     ω_0: $(ω_0)")

    if length(coordinates)==2 && ("e" in coordinates) && ("p" in coordinates)
        iE = coord_system_order[1]; ip = coord_system_order[2]
        E_array = abscissas[iE]; p_array = abscissas[ip]
        p_res = (1-Λ_∞)*norm(Equilibrium.Bfield(M,R_of_interest,z_of_interest))/norm(B_vec_on_axis)

        verbose && println("icrf_streamlines(): p_res: $(p_res)")

        dE = diff(E_array)[1] # Assume equidistant grid points
        dp = diff(p_array)[1] # Back-up dp, in case of NaN or Inf
        # For an m x n grid, epsilon will be an m x n array, where each element is a streamlines tangent (unit) 2-element vector
        epsilon = Array{Vector{Float64},space_dimensionality}(undef,space_size...)
        for c in space_indices
            E = E_array[c[iE]]; p = p_array[c[ip]]
            dp_c = -dE*(p^2)*E/(2*(E^2)*p)+(p_res^2)*dE*E/(E^2)/(2*p)
            if dp==Inf || isnan(dp)
                dE_c = 0.0
                dp_c = dp 
            else
                dE_c = dE
                dp_c = dp_c
            end

            tangent_vector = Vector{Float64}(undef,space_dimensionality)
            tangent_vector[iE] = dE_c
            tangent_vector[ip] = dp_c
            tangent_vector /= norm(tangent_vector)
            epsilon[c] = tangent_vector
        end
    elseif length(coordinates)==2 && ("vpara" in coordinates) && ("vperp" in coordinates)
        ivpa = coord_system_order[1]; ivpe = coord_system_order[2]
        vpa_array = abscissas[ivpa]; vpe_array = abscissas[ivpe]
        dvpa = diff(vpa_array)[1] # Assume equidistant grid points
        p_res = (1-Λ_∞)*norm(Equilibrium.Bfield(M,R_of_interest,z_of_interest))/norm(B_vec_on_axis)

        verbose && println("icrf_streamlines(): p_res: $(p_res)")

        # For an m x n grid, epsilon will be an m x n array, where each element is a streamlines tangent (unit) 2-element vector
        epsilon = Array{Vector{Float64},space_dimensionality}(undef,space_size...)
        for c in space_indices
            vpa = vpa_array[c[ivpa]]; vpe = vpe_array[c[ivpe]]
            if vpe==0.0
                dvpa_c = 0.0
                if vpa==0.0
                    dvpe_c = -dvpa
                else
                    dvpe_c = dvpa*sign(vpa)
                end
            else
                dvpa_c = dvpa
                dvpe_c = dvpa*(1-p_res^2)*vpa/(p_res*p_res*vpe)
            end

            tangent_vector = Vector{Float64}(undef,space_dimensionality)
            tangent_vector[ivpa] = dvpa_c
            tangent_vector[ivpe] = dvpe_c
            tangent_vector /= norm(tangent_vector)
            epsilon[c] = tangent_vector
        end
    elseif (length(coordinates)==3 || length(coordinates)==4) && ("e" in coordinates) && ("mu" in coordinates) && ("pphi" in coordinates)
        iE = coord_system_order[1]; imu = coord_system_order[2]; ipphi = coord_system_order[3]
        if ("sigma" in coordinates)
            isigma = coord_system_order[4]
        end

        E_array = abscissas[iE]; mu_array = abscissas[imu]; Pphi_array = abscissas[ipphi] # From keV to Joule

        dE = diff(E_array)[1] # Assume equidistant grid points
        dmu = dE*cyclotron_harmonic*getSpeciesCharge(FI_species)/(getSpeciesMass(FI_species)*ω)
        dPphi = dE*toroidal_mode_number/ω

        tangent_vector = Vector{Float64}(undef,space_dimensionality)
        tangent_vector[iE] = dE
        tangent_vector[imu] = dmu
        tangent_vector[ipphi] = dPphi
        if ("sigma" in coordinates)
            tangent_vector[isigma] = 0.0
        end
        tangent_vector /= norm(tangent_vector)

        # For an m x n x l grid, epsilon will be an m x n x l array, where each element is an ICRF streamlines tangent (unit) 3-element vector
        # For an m x n x l x 2 grid, epsilon will be an m x n x l x 2 array, where each element is a streamlines tangent (unit) 4-element vector, with the sigma directional element equal to 0
        epsilon = [tangent_vector for c in space_indices]
    elseif (length(coordinates)==3 || length(coordinates)==4) && ("e" in coordinates) && ("lambda" in coordinates) && ("pphi_n" in coordinates)
        iE = coord_system_order[1]; il = coord_system_order[2]; ippn = coord_system_order[3]
        if ("sigma" in coordinates)
            isigma = coord_system_order[4]
        end

        E_array = abscissas[iE]; Lambda_array = abscissas[il]; Pphi_n_array = abscissas[ippn]

        FI_species_q = getSpeciesCharge(FI_species)
        psi_axis, psi_bdry = psi_limits(M)
        if psi_bdry==0
            @warn "The magnetic flux at the last closed flux surface (LCFS) is found to be 0 for the 'M' input to os2COM(). Pϕ_n=Pϕ/(q*|Ψ_w|) where Ψ_w=Ψ(mag. axis) is assumed instead of Ψ_w=Ψ(LCFS)."
            Ψ_w_norm = abs(psi_axis)
        else
            Ψ_w_norm = abs(psi_bdry)
        end

        dE = diff(E_array)[1] # Assume equidistant grid points
        dPphi_n = inv(FI_species_q*Ψ_w_norm) *toroidal_mode_number*(1000*GuidingCenterOrbits.e0*dE)/ω

        # For an m x n x l grid, epsilon will be an m x n x l array, where each element is a streamlines tangent (unit) 3-element vector
        # For an m x n x l x 2 grid, epsilon will be an m x n x l x 2 array, where each element is a streamlines tangent (unit) 4-element vector, with the sigma directional element equal to 0
        epsilon = Array{Vector{Float64},space_dimensionality}(undef,space_size...)
        for c in space_indices
            L_diff = (Λ_∞ - Lambda_array[c[il]])
            if L_diff<0
                epsilon[c] = zeros(space_dimensionality)
                continue
            end
            tangent_vector = Vector{Float64}(undef,space_dimensionality)
            tangent_vector[iE] = dE
            tangent_vector[il] = L_diff*dE/E_array[c[iE]]
            tangent_vector[ippn] = dPphi_n
            if ("sigma" in coordinates)
                tangent_vector[isigma] = 0.0
            end
            tangent_vector /= norm(tangent_vector)
            epsilon[c] = tangent_vector
        end
    else
        error("Coordinate space $(coord_system) is not supported. Currently supported options include $(OWCF_ICRF_STREAMLINES_SUPPORTED_COORD_SYSTEMS) or any internal permutations thereof. Please correct and re-try.")
    end

    return epsilon
end

###### Data manipulation

"""
    add_noise(S, b)
    add_noise(-||-; k=0.1)

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

A function that estimates the noise of a signal S_tot. With the 'bad_inds' keyword argument,
estimate the noise of the S_tot-elements with the specific array indices in 'bad_inds'. 
Assume that the true signal is a smooth function without the noise.
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
        @warn "Instrumental response matrix input completely outside Ed_array input range. No diagnostic response will be applied."
        return S
    else
        if lo==1
            @warn "Lower bound of instrumental response matrix input might not be low enough to cover Ed_array input range. Diagnostic response representation might be inaccurate."
        end
        if hi==length(instrumental_response_input)
            @warn "Upper bound of instrumental response matrix input might not be high enough to cover Ed_array input range. Diagnostic response representation might be inaccurate."
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

Apply the instrumental response given by instrumental_response_matrix to each column of Q. Assume the diagnostic spectrum is accessed 
via the first dimension of Q. That is, Q[1,...] corresponds to the first weight function, Q[2,...] corresponds to the second weight function
and so on.
"""
function apply_instrumental_response(Q::Array{Float64,2}, Ed_array::AbstractVector, instrumental_response_input::AbstractVector, instrumental_response_output::AbstractVector, instrumental_response_matrix::AbstractMatrix; lo::Union{Int64,Nothing}=nothing, hi::Union{Int64,Nothing}=nothing)
    if isnothing(lo)
        lo = findfirst(x-> x>minimum(Ed_array),instrumental_response_input)
    end
    if isnothing(hi)
        hi = findlast(x-> x<maximum(Ed_array),instrumental_response_input)
    end
    
    if isnothing(lo) || isnothing(hi)
        @warn "Instrumental response matrix input completely outside Ed_array input range. No diagnostic response will be applied."
        return Q
    else
        Q_out = zeros(length(instrumental_response_output), size(Q,2))
        if lo==1
            @warn "Lower bound of instrumental response matrix input might not be low enough to cover Ed_array input range. Diagnostic response representation might be inaccurate."
        end
        if hi==length(instrumental_response_input)
            @warn "Upper bound of instrumental response matrix input might not be high enough to cover Ed_array input range. Diagnostic response representation might be inaccurate."
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

Apply the instrumental response given by instrumental_response_matrix to each matrix of Q. Assume the diagnostic spectrum is accessed 
via the first dimension of Q. That is, Q[1,...] corresponds to the first weight function, Q[2,...] corresponds to the second weight function
and so on.
"""
function apply_instrumental_response(Q::Array{Float64,3}, Ed_array::AbstractVector, instrumental_response_input::AbstractVector, instrumental_response_output::AbstractVector, instrumental_response_matrix::AbstractMatrix; lo::Union{Int64,Nothing}=nothing, hi::Union{Int64,Nothing}=nothing)
    if isnothing(lo)
        lo = findfirst(x-> x>minimum(Ed_array),instrumental_response_input)
    end
    if isnothing(hi)
        hi = findlast(x-> x<maximum(Ed_array),instrumental_response_input)
    end
    
    if isnothing(lo) || isnothing(hi)
        @warn "Instrumental response matrix input completely outside Ed_array input range. No diagnostic response will be applied."
        return Q
    else
        Q_out = zeros(length(instrumental_response_output), size(Q,2), size(Q,3))
        if lo==1
            @warn "Lower bound of instrumental response matrix input might not be low enough to cover Ed_array input range. Diagnostic response representation might be inaccurate."
        end
        if hi==length(instrumental_response_input)
            @warn "Upper bound of instrumental response matrix input might not be high enough to cover Ed_array input range. Diagnostic response representation might be inaccurate."
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

Apply the instrumental response given by instrumental_response_matrix to each 3D array of Q. Assume the diagnostic spectrum is accessed 
via the first dimension of Q. That is, Q[1,...] corresponds to the first weight function, Q[2,...] corresponds to the second weight function
and so on.
"""
function apply_instrumental_response(Q::Array{Float64,4}, Ed_array::AbstractVector, instrumental_response_input::AbstractVector, instrumental_response_output::AbstractVector, instrumental_response_matrix::AbstractMatrix; lo::Union{Int64,Nothing}=nothing, hi::Union{Int64,Nothing}=nothing)
    if isnothing(lo)
        lo = findfirst(x-> x>minimum(Ed_array),instrumental_response_input)
    end
    if isnothing(hi)
        hi = findlast(x-> x<maximum(Ed_array),instrumental_response_input)
    end
    
    if isnothing(lo) || isnothing(hi)
        @warn "Instrumental response matrix input completely outside Ed_array input range. No diagnostic response will be applied."
        return Q
    else
        Q_out = zeros(length(instrumental_response_output), size(Q,2), size(Q,3), size(Q,4))
        if lo==1
            @warn "Lower bound of instrumental response matrix input might not be low enough to cover Ed_array input range. Diagnostic response representation might be inaccurate."
        end
        if hi==length(instrumental_response_input)
            @warn "Upper bound of instrumental response matrix input might not be high enough to cover Ed_array input range. Diagnostic response representation might be inaccurate."
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

Apply the instrumental response given by instrumental_response_matrix to each 4D array of Q. Assume the diagnostic spectrum is accessed 
via the first dimension of Q. That is, Q[1,...] corresponds to the first weight function, Q[2,...] corresponds to the second weight function
and so on.
"""
function apply_instrumental_response(Q::Array{Float64,5}, Ed_array::AbstractVector, instrumental_response_input::AbstractVector, instrumental_response_output::AbstractVector, instrumental_response_matrix::AbstractMatrix; lo::Union{Int64,Nothing}=nothing, hi::Union{Int64,Nothing}=nothing)
    if isnothing(lo)
        lo = findfirst(x-> x>minimum(Ed_array),instrumental_response_input)
    end
    if isnothing(hi)
        hi = findlast(x-> x<maximum(Ed_array),instrumental_response_input)
    end
    
    if isnothing(lo) || isnothing(hi)
        @warn "Instrumental response matrix input completely outside Ed_array input range. No diagnostic response will be applied."
        return Q
    else
        Q_out = zeros(length(instrumental_response_output), size(Q,2), size(Q,3), size(Q,4), size(Q,5))
        if lo==1
            @warn "Lower bound of instrumental response matrix input might not be low enough to cover Ed_array input range. Diagnostic response representation might be inaccurate."
        end
        if hi==length(instrumental_response_input)
            @warn "Upper bound of instrumental response matrix input might not be high enough to cover Ed_array input range. Diagnostic response representation might be inaccurate."
        end
        for i=1:size(Q,5) # Apply instrumental response to each 4D array
            Q_out[:,:,:,:,i] = apply_instrumental_response(Q[:,:,:,:,i], Ed_array, instrumental_response_input, instrumental_response_output, instrumental_response_matrix; lo=lo, hi=hi)
        end
        return Q_out
    end
end

###### (Drift) Orbit-related functions

"""
    mu_func(energy, B, Pϕ, Ψ, RBϕ)
    mu_func(-||-; energy_in_keV=true, m=GuidingCenterOrbits.H2_amu*GuidingCenterOrbits.mass_u, q=1*GuidingCenterOrbits.e0)

Compute the magnetic moment mu, given the fast-ion energy, magnetic field B, toroidal canonical angular momentum Pϕ, magnetic flux Ψ and flux function RBϕ.
Use a 0-cap, meaning that all negative values are set to 0. The keyword arguments are:
- energy_in_keV - If true, the function will assume that the 'energy' input argument is given in kiloelectronvolt (keV). If false, assume Joule (J) - Bool
- m - The particle mass (kg) - Float64
- q - The particle charge (Coloumb) - Float64
"""
function mu_func(energy::T, B::T, Pϕ::T, Ψ::T, RBϕ::T; energy_in_keV::Bool=true, m::Float64=GuidingCenterOrbits.H2_amu*GuidingCenterOrbits.mass_u, q::Float64=1*GuidingCenterOrbits.e0) where {T}
    E = energy_in_keV ? energy*1000*GuidingCenterOrbits.e0 : energy # Else, assume energy is input in Joule
    res = E/B - (B/(2*m)) * ((Pϕ-q*Ψ)/RBϕ)^2
    return (res > zero(E)) ? res : zero(E)
end

"""
    get_orbel_volume(og)

Get the volume elements of the orbit-space grid voxels. 
Return a 3D array with all the volumes elements.
"""
function get_orbel_volume(og::OrbitGrid)
    return get3DVols(og.energy, og.pitch, og.r)
end
get_orbel_volume(og::OrbitGrid, dummy_bool::Bool) = get_orbel_volume(og) # For backwards compatibility

###### Inverse problem (fast-ion tomography) solving related functions

"""
    is_energy_pitch(abscissas,abscissas_units)
    is_energy_pitch(-||-; verbose=false, returnExtra=false)

Check if the reconstruction space discretized via the grid points in the abscissas in the Vector{Vector}
input variable 'abscissas' is the (E,p) space, where E is the energy and p=v_||/v is the pitch of the 
fast ion. The units of the abscissas are provided as elements in a Vector input via the 'abscissas_units' 
input variable. See OWCF/misc/convert_units.jl for how to specify units. If the reconstruction space is 
found to be (E,p), return true. If not, return false.

The keyword arguments are:
- verbose - If set to true, the function will talk a lot!
- returnExtra - If set to true, instead of a Bool, a Tuple will be returned, with the form (b,iE,ip).
                'b' is the Bool returned when returnExtra=false. 'iE' gives access to the E grid points 
                as 'abscissas[iE]'. 'ip' gives access to the p grid points as 'abscissas[ip]'.
"""
function is_energy_pitch(abscissas::Vector{Vector{T}} where {T<:Real}, abscissas_units::Vector{String}; verbose=false, returnExtra=false)
    if length(abscissas_units)!=2
        verbose && println("is_energy_pitch(): Number of reconstruction space abscissas is not equal to 2. Returning false... ")
        return (returnExtra ? (false, 0, 0) : false)
    end

    verbose && println("is_energy_pitch(): abscissas_units: $(abscissas_units)")

    units_1 = abscissas_units[1]
    units_2 = abscissas_units[2]
    units_tot = vcat(units_1, units_2)
    
    energy_ind = findall(x-> x in keys(ENERGY_UNITS) || x in keys(ENERGY_UNITS_LONG), units_tot)
    pitch_ind = findall(x-> x in keys(DIMENSIONLESS_UNITS) || x in keys(DIMENSIONLESS_UNITS_LONG), units_tot)

    if !(length(energy_ind)==1 && length(pitch_ind)==1)
        verbose && println("is_energy_pitch(): (E,p) coordinates not found (length(energy_ind)=$(length(energy_ind)))(length(pitch_ind)=$(length(pitch_ind))). Returning false... ")
        return (returnExtra ? (false, 0, 0) : false)
    end

    verbose && println("is_energy_pitch(): (E,p) coordinates confirmed! Returning true... ")
    return (returnExtra ? (true, energy_ind[1], pitch_ind[1]) : true) # [1] because we know there should be only 1 element returned by findall()
end

"""
    is_vpara_vperp(abscissas, abscissas_units)
    is_vpara_vperp(-||-; verbose=false, returnExtra=false)

Check if the reconstruction space discretized via the grid points in the abscissas in the Vector{Vector}
input variable 'abscissas' is the (vpara,vperp) space, where vpara and vperp are the particle velocities
parallel and perpendicular to the magnetic field, respectively. The units of the abscissas are provided 
as elements in a Vector input via the 'abscissas_units' input variable. See OWCF/misc/convert_units.jl 
for how to specify units. If the reconstruction space is found to be (vpara,vperp), return true. If not, 
return false.

The keyword arguments are:
- verbose - If set to true, the function will talk a lot!
- returnExtra - If set to true, instead of a Bool, a Tuple will be returned, with the form (b,ivpa,ivpe).
                'b' is the Bool returned when returnExtra=false. 'ivpa' gives access to the vpara grid points 
                as 'abscissas[ivpa]'. 'ivpe' gives access to the vperp grid points as 'abscissas[ivpe]'.
"""
function is_vpara_vperp(abscissas::Vector{Vector{T}} where {T<:Real}, abscissas_units::Vector{String}; verbose=false, returnExtra=false)
    if length(abscissas_units)!=2
        verbose && println("is_vpara_vperp(): Number of reconstruction space abscissas is not equal to 2. Returning false... ")
        return (returnExtra ? (false, 0, 0) : false)
    end

    units_1 = abscissas_units[1]
    units_2 = abscissas_units[2]
    units_tot = vcat(units_1, units_2)
    w_speed_inds = findall(x-> units_are_speed(x), units_tot)

    if !(length(w_speed_inds)==2)
        verbose && println("is_vpara_vperp(): (vpara,vperp) coordinates not found. Returning false... ")
        return (returnExtra ? (false, 0, 0) : false)
    else
        verbose && print("is_vpara_vperp(): (vpara,vperp) coordinates found! Distinguishing (vpara, vperp) arrays...")
        w_vel_arrays = abscissas[w_speed_inds]
        if minimum(w_vel_arrays[1])<0 && minimum(w_vel_arrays[2])>=0 # vpara can be negative. vperp cannot
            verbose && println("ok!")
            w_vpara_ind = w_speed_inds[1] # Order of abscissas in w_vel_arrays is the same as the order in w_speed_inds
            w_vperp_ind = w_speed_inds[2] # Order of abscissas in w_vel_arrays is the same as the order in w_speed_inds
        elseif minimum(w_vel_arrays[2])<0 && minimum(w_vel_arrays[1])>=0 # vpara can be negative. vperp cannot
            verbose && println("ok!")
            w_vpara_ind = w_speed_inds[2] # Order of abscissas in w_vel_arrays is the same as the order in w_speed_inds
            w_vperp_ind = w_speed_inds[1] # Order of abscissas in w_vel_arrays is the same as the order in w_speed_inds
        else
            verbose && println("")
            @warn "Could not distinguish (vpara,vperp) arrays from weight function abscissas. Assuming abscissa with index $(w_speed_inds[1]) to be vpara and abscissa with index $(w_speed_inds[2]) to be vperp."
            w_vpara_ind = w_speed_inds[1] # Order of abscissas in w_vel_arrays is the same as the order in w_speed_inds
            w_vperp_ind = w_speed_inds[2] # Order of abscissas in w_vel_arrays is the same as the order in w_speed_inds
        end
        verbose && println("is_vpara_vperp(): Returning true... ")
        return (returnExtra ? (true, w_vpara_ind, w_vperp_ind) : true)
    end
end

"""
    is_COM(abscissas, abscissas_units)
    is_COM(-||-; verbose=false, returnExtra=false)

Check if the reconstruction space discretized via the grid points in the abscissas in the Vector{Vector}
input variable 'abscissas' is the constants-of-motion (COM) space, i.e. (E,mu,Pphi) or (E,mu,Pphi,sigma),
where E is the particle energy, mu is the magnetic moment and Pphi is the toroidal canonical angular 
momentum. sigma is a binary coordinate, included to distinguish between co- and counter-going orbits.
The units of the abscissas are provided as elements in a Vector input via the 'abscissas_units' input variable. 
See OWCF/misc/convert_units.jl for how to specify units. If the reconstruction space is found to be (E,mu,Pphi)
or (E,mu,Pphi,sigma), return true. If not, return false.

The keyword arguments are:
- verbose - If set to true, the function will talk a lot!
- returnExtra - If set to true, instead of a Bool, a Tuple will be returned, with the form (b,iE,imu,iPphi) or 
                (b,iE,imu,iPphi,isigma), depending on the length of input 'abscissas'. 'b' is the Bool returned 
                when returnExtra=false. 'iE' gives access to the E grid points as 'abscissas[iE]'. Similar for 
                imu, iPphi and isigma.
"""
function is_COM(abscissas::Vector{Vector{T}} where {T<:Real}, abscissas_units::Vector{String}; verbose=false, returnExtra=false)
    if length(abscissas_units)!=3 && length(abscissas_units)!=4
        verbose && println("is_COM(): Number of reconstruction space abscissas is not equal to 3, nor 4. Returning false... ")
        return (returnExtra ? (false, zeros(Int64,length(abscissas_units))...) : false)
    end

    units_1 = abscissas_units[1]
    units_2 = abscissas_units[2]
    units_3 = abscissas_units[3]
    units_4 = length(abscissas_units)==4 ? abscissas_units[4] : "dimensionless"

    units_tot = vcat(units_1, units_2, units_3, units_4)
    w_energy_ind = findall(x-> x in keys(ENERGY_UNITS) || x in keys(ENERGY_UNITS_LONG), units_tot)
    w_mu_ind = findall(x-> units_are_equal_base(x,"m^2_A"), units_tot) # SI units of magnetic moment
    w_Pphi_ind = findall(x-> units_are_equal_base(x,"kg_m^2_s^-1"), units_tot) # SI units of toroidal angular canonical momentum
    w_sigma_ind = findall(x-> x in keys(DIMENSIONLESS_UNITS) || x in keys(DIMENSIONLESS_UNITS_LONG), units_tot)

    if !(length(w_energy_ind)==1 && length(w_mu_ind)==1 && length(w_Pphi_ind)==1 && length(w_sigma_ind)==1)
        verbose && println("is_COM(): (E,mu,Pphi) coordinates not found. Returning false... ")
        return (returnExtra ? (false, zeros(Int64,length(abscissas_units))...) : false)
    end

    verbose && println("is_COM(): (E,mu,Pphi) coordinates confirmed! Returning true... ")
    if !returnExtra
        return true
    end
    if length(abscissas_units)==4
        # [1] because we know there should be only 1 element returned by findall()
        output_tuple = (true, w_energy_ind[1], w_mu_ind[1], w_Pphi_ind[1], w_sigma_ind[1])
    else # Must be 3 (see above)
        # [1] because we know there should be only 1 element returned by findall()
        output_tuple = (true, w_energy_ind[1], w_mu_ind[1], w_Pphi_ind[1])
    end
    return (returnExtra ? output_tuple : true)
end

"""
    is_normalized_COM(abscissas, abscissas_units)
    is_normalized_COM(-||-; verbose=false, returnExtra=false)

Check if the reconstruction space discretized via the grid points in the abscissas in the Vector{Vector}
input variable 'abscissas' is the normalized constants-of-motion (COM) space, i.e. (E,Lambda,Pphi_n) or 
(E,Lambda,Pphi_n,sigma), where E is the particle energy, Lambda is the normalized magnetic moment and 
Pphi_n is the toroidal canonical angular momentum. sigma is a binary coordinate, included to distinguish 
between co- and counter-going orbits. The units of the abscissas are provided as elements in a Vector input 
via the 'abscissas_units' input variable. See OWCF/misc/convert_units.jl for how to specify units. If the 
reconstruction space is found to be (E,Lambda,Pphi_n) or (E,Lambda,Pphi_n,sigma), return true. If not, 
return false.

The keyword arguments are:
- verbose - If set to true, the function will talk a lot!
- returnExtra - If set to true, instead of a Bool, a Tuple will be returned, with the form (b,iE,iL,iPpn) or 
                (b,iE,iL,iPpn,isigma), depending on the length of input 'abscissas'. 'b' is the Bool returned 
                when returnExtra=false. 'iE' gives access to the E grid points as 'abscissas[iE]'. Similar for 
                iL, iPpn and isigma.
"""
function is_normalized_COM(abscissas::Vector{Vector{T}} where {T<:Real}, abscissas_units::Vector{String}; verbose=false, returnExtra=false)
    rec_space_DIM = length(abscissas_units) # Reconstruction space dimensionality is the number of weight function abscissas

    if rec_space_DIM!=3 && rec_space_DIM!=4
        verbose && println("is_normalized_COM(1): Number of reconstruction space abscissas is not equal to 3, nor 4. Returning false... ")
        return (returnExtra ? (false, zeros(Int64,rec_space_DIM)...) : false)
    end

    sigma_included = false # By default, assume that the binary coordinate sigma is not included as a reconstruction space coordinate
    if rec_space_DIM==4
        sigma_included = true # If the dimensionality is 4, sigma is assumed to be one of the reconstruction space coordinates
    end

    units_1 = abscissas_units[1]
    units_2 = abscissas_units[2]
    units_3 = abscissas_units[3]
    units_4 = sigma_included ? abscissas_units[4] : "dimensionless"

    units_tot = vcat(units_1, units_2, units_3, units_4)
    w_energy_ind = findall(x-> x in keys(ENERGY_UNITS) || x in keys(ENERGY_UNITS_LONG), units_tot)
    w_dimensionless_inds = findall(x-> x in keys(DIMENSIONLESS_UNITS) || x in keys(DIMENSIONLESS_UNITS_LONG), units_tot)

    if !(length(w_energy_ind)==1 && length(w_dimensionless_inds)==3)
        verbose && println("is_normalized_COM(2): (E,Lambda,Pphi_n) coordinates not found. Returning false... ")
        return (returnExtra ? (false, zeros(Int64,rec_space_DIM)...) : false)
    end

    if !sigma_included # If sigma was NOT included as a reconstruction space abscissa
        filter!(x-> x!=4,w_dimensionless_inds) # Remove the dummy index with value 4. We know it will be present in w_dimensionless_inds, because we added it ourselves (see above)
    end

    verbose && println("is_normalized_COM(3): (E,Lambda,Pphi_n) coordinates assumed. Attempting to deduce coordinate order... ")
    w_Lambda_ind = nothing
    w_Pphi_n_ind = nothing
    w_sigma_ind = sigma_included ? nothing : -1 # If reconstruction space abscissas do not include sigma, initialize to -1

    for ind in w_dimensionless_inds # For each index ind corresponding to an abscissa for one of the dimensionless coordinates Lambda, Pphi_n or sigma (we don't know which one, we try to deduce it below)
        abscissa = abscissas[ind]
        if minimum(abscissa)<0 && maximum(abscissa)<=0 # If all grid points of the abscissa are placed at zero and more negative values
            if isnothing(w_Pphi_n_ind) # and the Pphi_n index has not yet been set
                w_Pphi_n_ind = ind # We know ind MUST be the Pphi_n index, since this is NOT possible for the other coordintes (mu and sigma)
                continue # Continue to next ind in w_dimensionless_inds
            else # If the Pphi_n index has already been set, the assumption of (E,Lambda,Pphi_n) coordinates cannot be correct
                verbose && println("is_normalized_COM(4): The assumption of (E,Lambda,Pphi_n) coordinates is found to be incorrect. Returning false... ")
                return (returnExtra ? (false, zeros(Int64,rec_space_DIM)...) : false)
            end
        end
        if minimum(abscissa)>=0 && maximum(abscissa)>0 # If the grid points of the abscissa are placed at zero and greater positive values
            if isnothing(w_Lambda_ind) # and the Lambda index has not yet been set
                w_Lambda_ind = ind # Assume ind is the Lambda index
                continue # Continue to next ind in w_dimensionless_inds
            else # If the Lambda index has already been set, the deduction becomes practically impossible, since then the abscissas do not span the whole of COM-space
                @warn "is_normalized_COM(5): Unable to deduce (E,Lambda_Pphi_n,sigma) coordinate order. Assuming normal order $(sigma_included ? "(i_E,i_Lambda,i_Pphi_n;i_sigma)=(1,2,3;4)" : "(i_E,i_Lambda,i_Pphi_n)=(1,2,3)")"
                return (returnExtra ? (true, (1:rec_space_DIM)...) : true) # Start at index 2, to align with weight function abscissa indexing (index 1 is always the diagnostic measurement bin centers)
            end
        end
        if iszero(abscissa) # If all grid points are 0, we can do nothing, this is some error-level sh*t!
            @warn "is_normalized_COM(6): Unable to deduce (E,Lambda_Pphi_n,sigma) coordinate order. Assuming normal order $(sigma_included ? "(i_E,i_Lambda,i_Pphi_n;i_sigma)=(1,2,3;4)" : "(i_E,i_Lambda,i_Pphi_n)=(1,2,3)")"
            return (returnExtra ? (true, (1:rec_space_DIM)...) : true) # Start at index 2, to align with weight function abscissa indexing (index 1 is always the diagnostic measurement bin centers)
        end
        if length(abscissa)==2 # If the abscissa has exactly two grid points (one positive and one negative)
            if isnothing(w_sigma_ind) # and the sigma index has not yet been set
                w_sigma_ind = ind # (It's extremely reasonable to) assume ind is the sigma index
                continue
            elseif w_sigma_ind==-1 # If the reconstruction space dimensionality is 3
                # Don't do anything, this must be the w_Pphi_n_ind, and will be set in the if statement below
            else # If the sigma index has ALREADY been set (and the reconstruction space dimensionality is NOT 3, but 4)
                # Then, the Pphi_n and sigma abscissas have exactly the same characteristics (2 grid points, one negative and one positive). Distinguishing between the two is impossible
                @warn "is_normalized_COM(7): Unable to deduce (E,Lambda_Pphi_n,sigma) coordinate order. Assuming normal order $(sigma_included ? "(i_E,i_Lambda,i_Pphi_n;i_sigma)=(1,2,3;4)" : "(i_E,i_Lambda,i_Pphi_n)=(1,2,3)")"
                return (returnExtra ? (true, (1:rec_space_DIM)...) : true) # Start at index 2, to align with weight function abscissa indexing (index 1 is always the diagnostic measurement bin centers)
            end
        end
        # If the grid points of the abscissa are placed at both negative and positive values, and there are not exactly two grid points...
        if isnothing(w_Pphi_n_ind) # and the Pphi_n index has not yet been set
            w_Pphi_n_ind = ind # We know ind MUST be the Pphi_n index, since this is NOT possible for the other coordintes (mu and sigma)
            continue # Continue to next ind in w_dimensionless_inds
        else # If the Pphi_n index has already been set, the assumption of (E,Lambda,Pphi_n) coordinates cannot be correct
            verbose && println("is_normalized_COM(8): The assumption of (E,Lambda,Pphi_n) coordinates is found to be incorrect. Returning false... ")
            return (returnExtra ? (false, zeros(Int64,rec_space_DIM)...) : false)
        end
    end

    if sigma_included && sum(isnothing.([w_Lambda_ind, w_Pphi_n_ind, w_sigma_ind]))==0
        verbose && println("is_normalized_COM(9): Coordinate order deduced. Returning true... ")
        output_tuple = (true, w_energy_ind[1], w_Lambda_ind[1], w_Pphi_n_ind[1], w_sigma_ind[1])
    elseif !sigma_included && sum(isnothing.([w_Lambda_ind, w_Pphi_n_ind]))==0
        verbose && println("is_normalized_COM(10): Coordinate order deduced. Returning true... ")
        output_tuple = (true, w_energy_ind[1], w_Lambda_ind[1], w_Pphi_n_ind[1])
    else # There must have been some error. This should not be possible.
        error("This error should be logically impossible to reach. Please post an issue at https://github.com/juliaFusion/owCF/issues or try to directly contact henrikj@dtu.dk or anvalen@dtu.dk.")
    end
    return (returnExtra ? output_tuple : true)
end

"""
    is_EpRz(abscissas, abscissas_units)
    is_EpRz(-||-; verbose=false, returnExtra=false)

Check if the reconstruction space discretized via the grid points in the abscissas in the Vector{Vector}
input variable 'abscissas' is the (E,p,R,z) space, where E is the particle energy, p is the particle pitch 
(pitch=v_||/v), R is the major radius and z is the vertical coordinate. The units of the abscissas are 
provided as elements in a Vector input via the 'abscissas_units' input variable. See OWCF/misc/convert_units.jl 
for how to specify units. If the reconstruction space is found to be (E,p,R,z), return true. If not, return false.

The keyword arguments are:
- verbose - If set to true, the function will talk a lot!
- returnExtra - If set to true, instead of a Bool, a Tuple will be returned, with the form (b,iE,ip,iR,iz). 
                'b' is the Bool returned when returnExtra=false. 'iE' gives access to the E grid points as 
                'abscissas[iE]'. Similar for ip, iR and iz.
"""
function is_EpRz(abscissas::Vector{Vector{T}} where {T<:Real}, abscissas_units::Vector{String}; verbose=false, returnExtra=false)
    if length(abscissas_units)!=4
        verbose && println("is_EpRz(): Number of reconstruction space abscissas is not equal to 4. Returning false... ")
        return (returnExtra ? (false, 0, 0, 0, 0, [], []) : false)
    end
    units_1 = abscissas_units[1]
    units_2 = abscissas_units[2]
    units_3 = abscissas_units[3]
    units_4 = abscissas_units[4]

    units_tot = vcat(units_1, units_2, units_3, units_4)
    w_energy_ind = findall(x-> x in keys(ENERGY_UNITS) || x in keys(ENERGY_UNITS_LONG), units_tot)
    w_pitch_ind = findall(x-> x in keys(DIMENSIONLESS_UNITS) || x in keys(DIMENSIONLESS_UNITS_LONG), units_tot)
    w_Rz_inds = findall(x-> x in keys(LENGTH_UNITS) || x in keys(LENGTH_UNITS_LONG), units_tot)

    if !(length(w_energy_ind)==1 && length(w_pitch_ind)==1 && length(w_Rz_inds)==2)
        verbose && println("is_EpRz(): (E,p,R,z) coordinates not found. Returning false... ")
        return (returnExtra ? (false, 0, 0, 0, 0, [], []) : false)
    end

    verbose && print("is_EpRz(): (E,p,R,z) coordinates found! Distinguishing (R,z) arrays... ")
    w_RnZ_arrays = abscissas[w_Rz_inds]
    if minimum(w_RnZ_arrays[2])<0 && minimum(w_RnZ_arrays[1])>0 # If the second LENGTH_UNITS abscissa has negative elements, and the first one does not..
        verbose && println("ok!")
        w_R_ind = w_Rz_inds[1] # Order of abscissas in w_RnZ_arrays is the same as the order in w_Rz_inds
        w_z_ind = w_Rz_inds[2] # Order of abscissas in w_RnZ_arrays is the same as the order in w_Rz_inds
        R_of_interests = w_RnZ_arrays[1] # The first LENGTH_UNITS abscissa is very likely to be the R grid points...
        z_of_interests = w_RnZ_arrays[2] # ...and the second LENGTH_UNITS abscissa is very likely to be the z grid points
    elseif minimum(w_RnZ_arrays[1])<0 && minimum(w_RnZ_arrays[2])>0 # If it's the other way around...
        verbose && println("ok!")
        w_R_ind = w_Rz_inds[2] # Order of abscissas in w_RnZ_arrays is the same as the order in w_Rz_inds
        w_z_ind = w_Rz_inds[1] # Order of abscissas in w_RnZ_arrays is the same as the order in w_Rz_inds
        R_of_interests = w_RnZ_arrays[2] # ...it's very likely to be the other way around.
        z_of_interests = w_RnZ_arrays[1] # ...it's very likely to be the other way around.
    else
        verbose && println("")
        @warn "Could not deduce (R,z) arrays from weight function abscissas. Assuming abscissa with index $(w_Rz_inds[1]) to be R and abscissa with index $(w_Rz_inds[2]) to be z."
        w_R_ind = w_Rz_inds[1] # Order of abscissas in w_RnZ_arrays is the same as the order in w_Rz_inds
        w_z_ind = w_Rz_inds[2] # Order of abscissas in w_RnZ_arrays is the same as the order in w_Rz_inds
        R_of_interests = w_RnZ_arrays[1]
        z_of_interests = w_RnZ_arrays[2]
    end
    verbose && println("is_EpRz(): Returning true... ")
    return (returnExtra ? (true, w_energy_ind[1], w_pitch_ind[1], w_R_ind, w_z_ind, R_of_interests, z_of_interests) : true)
end

"""
    get_energy_abscissa(abscissas, abscissas_units)
    get_energy_abscissa(-||-; verbose=false, returnExtra=false)

If there is an energy abscissa in the collection of abscissas input 'abscissas', return it. Otherwise, return nothing.

The keyword arguments are:
    - verbose - If true, the function will talk a lot!
    - returnExtra - If set to true, instead of a Bool, a Tuple will be returned, with the form (b,iE). 
                'b' is the Bool returned when returnExtra=false. 'iE' gives access to the energy grid points as 
                'abscissas[iE]'.
"""
function get_energy_abscissa(abscissas::Vector{Vector{T}} where {T<:Real}, abscissas_units::Vector{String}; verbose=false, returnExtra=false)
    # Assume there is only one energy abscissa
    w_energy_ind = findfirst(x-> x in keys(ENERGY_UNITS) || x in keys(ENERGY_UNITS_LONG), abscissas_units)
    if isnothing(w_energy_ind)
        verbose && println("No energy abscissa found in input abscissas! Returning nothing... ")
        return returnExtra ? (nothing,nothing) : nothing
    end

    return returnExtra ? (abscissas[w_energy_ind], w_energy_ind) : abscissas[w_energy_ind]
end

"""
    get_pitch_abscissa(abscissas, abscissas_units)
    get_pitch_abscissa(-||-; verbose=false, returnExtra=false)

If there is a pitch abscissa in the collection of abscissas input 'abscissas', return it. Otherwise, return nothing.

The keyword arguments are:
    - verbose - If true, the function will talk a lot!
    - returnExtra - If set to true, instead of a Bool, a Tuple will be returned, with the form (b,ip). 
                'b' is the Bool returned when returnExtra=false. 'ip' gives access to the pitch grid points as 
                'abscissas[ip]'.
"""
function get_pitch_abscissa(abscissas::Vector{Vector{T}} where {T<:Real}, abscissas_units::Vector{String}; verbose=false, returnExtra=false)
    w_dimensionless_inds = findall(x-> x in keys(DIMENSIONLESS_UNITS) || x in keys(DIMENSIONLESS_UNITS_LONG), abscissas_units)

    if isempty(w_dimensionless_inds)
        verbose && println("No dimensionless abscissa found in input abscissas! Returning nothing...")
        return returnExtra ? (nothing,nothing) : nothing
    end

    w_pitch_ind = 0
    for w_dimensionless_ind in w_dimensionless_inds
        w_abscissa = abscissas[w_dimensionless_ind]
        if maximum(w_abscissa)<=1 && minimum(w_abscissa)>=-1
            w_pitch_ind = w_dimensionless_ind
            break # Assume there is only one pitch-like coordinate
        end
    end

    if w_pitch_ind==0
        verbose && println("No pitch-like dimensionless abscissa found in input abscissas! Returning nothing...")
        return returnExtra ? (nothing,nothing) : nothing
    end

    return returnExtra ? (abscissas[w_pitch_ind], w_pitch_ind) : abscissas[w_pitch_ind]
end

"""
    get_vpara_abscissa(abscissas, abscissas_units)
    get_vpara_abscissa(-||-; verbose=false, returnExtra=false)

If there is a vpara abscissa in the collection of abscissas input 'abscissas', return it. Otherwise, return nothing.

The keyword arguments are:
    - verbose - If true, the function will talk a lot!
    - returnExtra - If set to true, instead of a Bool, a Tuple will be returned, with the form (b,ivpa). 
                'b' is the Bool returned when returnExtra=false. 'ivpa' gives access to the vpara grid points as 
                'abscissas[ivpa]'.
"""
function get_vpara_abscissa(abscissas::Vector{Vector{T}} where {T<:Real}, abscissas_units::Vector{String}; verbose=false, returnExtra=false)
    w_speed_inds = findall(x-> units_are_speed(x), abscissas_units)
    if isempty(w_speed_inds)
        verbose && println("No speed abscissa found in input abscissas! Returning nothing...")
        return returnExtra ? (nothing,nothing) : nothing
    end
    if length(w_speed_inds)>2
        @warn "Impossible to determine vpara since more than two abscissas with unit of measurement 'speed' was found. Returning nothing... "
        return returnExtra ? (nothing,nothing) : nothing
    end

    if length(w_speed_inds)==1
        verbose && println("One abscissa with unit of measurement 'speed' was found. Assuming it is vpara. Returning... ")
        return returnExtra ? (abscissas[w_speed_inds], w_speed_inds[1]) : abscissas[w_speed_inds]
    end

    # w_speed_inds must be length 2 then
    w_speed_abscissa_1 = abscissas[w_speed_inds[1]]
    w_speed_abscissa_2 = abscissas[w_speed_inds[2]]
    if minimum(w_speed_abscissa_1)<0 # Must be vpara since vperp>= always holds
        return returnExtra ? (w_speed_abscissa_1, w_speed_inds[1]) : w_speed_abscissa_1
    elseif minimum(w_speed_abscissa_2)<0 # -||-
        return returnExtra ? (w_speed_abscissa_2, w_speed_inds[2]) : w_speed_abscissa_2
    else # Ambiguous
        @warn "Cannot deduce vpara from abscissas and abscissas_units. Returning abscissas with smallest index as vpara."
        return returnExtra ? (w_speed_abscissa_1, w_speed_inds[1]) : w_speed_abscissa_1
    end
end

"""
    get_vperp_abscissa(abscissas, abscissas_units)
    get_vperp_abscissa(-||-; verbose=false, returnExtra=false)

If there is a vperp abscissa in the collection of abscissas input 'abscissas', return it. Otherwise, return nothing.

The keyword arguments are:
    - verbose - If true, the function will talk a lot!
    - returnExtra - If set to true, instead of a Bool, a Tuple will be returned, with the form (b,ivpe). 
                'b' is the Bool returned when returnExtra=false. 'ivpe' gives access to the vperp grid points as 
                'abscissas[ivpe]'.
"""
function get_vperp_abscissa(abscissas::Vector{Vector{T}} where {T<:Real}, abscissas_units::Vector{String}; verbose=false, returnExtra=false)
    w_speed_inds = findall(x-> units_are_speed(x), abscissas_units)
    if isempty(w_speed_inds)
        verbose && println("No speed abscissa found in input abscissas! Returning nothing...")
        return returnExtra ? (nothing,nothing) : nothing
    end
    if length(w_speed_inds)>2
        @warn "Impossible to determine vperp since more than two abscissas with unit of measurement 'speed' was found. Returning nothing... "
        return returnExtra ? (nothing,nothing) : nothing
    end

    if length(w_speed_inds)==1
        verbose && println("One abscissa with unit of measurement 'speed' was found. Assuming it is vperp. Returning... ")
        return returnExtra ? (abscissas[w_speed_inds], w_speed_inds[1]) : abscissas[w_speed_inds]
    end

    # w_speed_inds must be length 2 then
    w_vpara_abscissa, w_vpara_ind = get_vpara_abscissa(abscissas, abscissas_units; verbose=verbose, returnExtra=returnExtra)
    w_vperp_ind = filter(x-> x!=w_vpara_ind,w_speed_inds)[1]

    return returnExtra ? (abscissas[w_vperp_ind], w_vperp_ind) : abscissas[w_vperp_ind]
end

"""
    icrf_regularization_matrix(M, abscissas, FI_species, wave_frequency, cyclotron_harmonic)
    icrf_regularization_matrix(-||-; toroidal_mode_number=1, coord_system="(E,p)", coord_system_order=Tuple(1:length(abscissas)), R_of_interest=R_axis, 
                           z_of_interest=z_axis, L1=forward_difference_matrix(Tuple(length.(abscissas))), verbose=false)

Compute a regularization matrix that enables a solution to an inverse problem to be regularized with respect to the physics of electromagnetic wave heating in 
the ion cyclotron range of frequencies (ICRF). Currently, the available phase-space options include: 
    - (E,p), or (p,E), where E is the fast-ion energy in keV, and p is the pitch (vpara/v with v and vpara being the total speed and the speed component 
      parallel to the magnetic field, respectively)
    - (vpara,vperp), or (vperp,vpara), where vpara and vperp are the speed components parallel and perpendicular to the magnetic field, respectively
    - (E,mu,Pphi), or any permutations thereof, where E is the fast-ion energy in keV, mu is the magnetic moment in A*m^2 and Pphi is the toroidal canonical 
      angular momentum in kg*m^2*s^-1
    - (E,mu,Pphi,sigma), or any permutations thereof, where -||- and sigma is a binary coordinate (-1,+1)
    - (E,Lambda,Pphi_n), or any permutations thereof, where -||-, Lambda is the normalized magnetic moment (Lambda=mu*B0/E with B0 the magnetic field on-axis 
      in Teslas and E is the fast-ion energy in Joules) and Pphi_n is the normalized toroidal canonical angular momentum (Pphi_n=Pphi/(q*Ψ_w) with q being the 
      particle charge in Coulomb and Ψ_w the magnetic flux at the separatrix)
    - (E,Lambda,Pphi_n,sigma), or any permutations thereof, where -||- and -||-

The information regarding which one of these phase-space coordinate spaces is discretized via grid points, included as the 'abscissas' input, must be specified 
using the 'coord_system' keyword argument. An explanation of the input variables are as follows:
    - M - The magnetic equilibrium object, as given by the Equilibrium.jl package - AbstractEquilibrium
    - abscissas - The grid points of the abscissas of the phase space of interest - Vector{Vector}
    - FI_species - The (fast-ion) particle of interest. Please see OWCF/misc/species_func.jl for a list of available particle species - String
    - wave_frequency - The wave frequency of the ICRF wave. In Hz - Float64
    - cyclotron_harmonic - The cyclotron harmonic integer of the ICRF heating scheme, for the 'FI_species' particle species and a cyclotron frequency on-axis - Int64
The keyword arguments are as follows:
    - toroidal_mode_number - The toroidal mode number of the ICRF wave. Defaults to 1 - Int64
    - coord_system - The coordinate system of the phase space of interest (please see explanation above) - String
    - coord_system_order - If phase-space abscissas in 'abscissas' input are not given in the common order, e.g. (E,Pphi,mu) instead of (E,mu,Pphi), this needs to be 
                           specified via the 'coord_system_order' keyword argument as e.g. (1,3,2) instead of the common (1,2,3) order etc. Defaults to (1:length(abscissas)) - Tuple
    - R_of_interest - If "(E,p)", "(p,E)", "(vpara,vperp)" or "(vperp,vpara)" are given as input to the 'coord_system' keyword argument, an (R,z) point of 
      interest needs to be specified. In meters. - Float64
    - z_of_interest - If "(E,p)", "(p,E)", "(vpara,vperp)" or "(vperp,vpara)" are given as input to -||-. In meters. - Float64
    - L1 - The finite difference matrix for the coordinate system space, as returned by e.g. the forward_difference_matrix() OWCF function (see OWCF/extra/dependencies.jl) - Matrix
    - verbose - If set to true, the function will talk a lot! - Bool

The ICRF regularization matrix will be returned as a Matrix{Float64} with nCols=reduce(*,Tuple(length.(abscissas))) number of columns and 
length(abscissas)*nCols number of rows.

PLEASE NOTE! Energy grid points are assumed to be given in keV for all coordinate systems!!!
Velocity grid points are assumed to be given in m/s for all coordinate systems!!!

The icrf_regularization_matrix() function is based on the content in M. Rud et al 2025 Nucl. Fusion 65 056008.
"""
function icrf_regularization_matrix(M::AbstractEquilibrium, abscissas::Vector{Vector{T}} where {T<:Real}, FI_species::String, wave_frequency::Float64, cyclotron_harmonic::Int64;
                           toroidal_mode_number::Int64=1, coord_system::String="(E,p)", coord_system_order=Tuple(1:length(abscissas)), 
                           R_of_interest=magnetic_axis(M)[1], z_of_interest=magnetic_axis(M)[2], L1::AbstractMatrix=forward_difference_matrix(Tuple(length.(abscissas))),
                           verbose::Bool=false)
    epsilon = icrf_streamlines(M, abscissas, FI_species, wave_frequency, cyclotron_harmonic; 
                               toroidal_mode_number=toroidal_mode_number, coord_system=coord_system, coord_system_order=coord_system_order, 
                               R_of_interest=R_of_interest, z_of_interest=z_of_interest, verbose=verbose)

    # Necessary space quantities
    space_dimensionality = length(abscissas) # The space dimensionality, e.g. 3
    space_size = Tuple(length.(abscissas)) # The size of the space, e.g. (33,52,12). 
    space_indices = CartesianIndices(space_size) # All space indices. E.g. [(1,1,1), (1,1,2), ..., (N1,N2,N3)] where Ni is the number of grid points in the i:th space dimension
    numOCoords = length(space_indices) # The number of grid points

    verbose && println("icrf_regularization_matrix() function input info:")
    verbose && println("---> Coordinate system input: $(coord_system)")
    verbose && println("---> Space dimensionality: $(space_dimensionality)")
    verbose && println("---> Space size: $(space_size)")
    verbose && println("---> Number of grid points: $(numOCoords)")
    verbose && println("---> L1 (finite differences) matrix sparsity: $(round(100*(1-length(L1.nzval)/length(L1)),sigdigits=5)) %")

    L1_ICRF = spzeros(size(L1))

    for (i,c) in enumerate(space_indices)
        epsilon_i = epsilon[c]
        finite_difference = Vector{Vector{Float64}}(undef,space_dimensionality)
        for d in coord_system_order
            finite_difference[d] = L1[i+(d-1)*numOCoords,:]
        end
        # The projection equation. 'epsilon*transpose(epsilon)' will be a space_dimensionality x space_dimensionality matrix. finite_difference will be a space_dimensionality x 1 Vector.
        icrf_difference = (epsilon_i * transpose(epsilon_i)) * finite_difference 
        for d in coord_system_order
            L1_ICRF[i+(d-1)*numOCoords,:] = icrf_difference[d] # Put each row (icrf_difference) where it belongs (the right coordinate 'i' and the right gradient component 'd')
        end
    end

    verbose && println("icrf_regularization_matrix() function output info:")
    verbose && println("---> L1_ICRF (L1 projected onto ICRF streamlines) matrix sparsity: $(round(100*(1-length(L1_ICRF.nzval)/length(L1_ICRF)),sigdigits=5)) %")

    return L1_ICRF
end

###### UNLABELLED FUNCTIONS BELOW

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
    for Ei in eachindex(Eq)
        for pi in eachindex(pq)
            for Ri in eachindex(Rq)
                for zi in eachindex(zq)
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
    do4Dto2D(-||-, returnComplement=false)

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
    for c in axes(W4D,1)
        for Ei in axes(W4D,2)
            for pmi in axes(W4D,3)
                for Rmi in axes(W4D,4)
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
    flip_along_pm(W, with_channels=false)

Flip the elements of an orbit-space function along the pm axis. If with_channels=true,
then a 4D orbit-space function on the form (channel,E,pm,Rm) or (channel,Rm,pm,E) is assumed.
Otherwise a form of (E,pm,Rm) is assumed.
"""
function flip_along_pm(W::AbstractArray; with_channels::Bool=false)

    Wres = zeros(size(W))
    if with_channels
        for pmi in axes(W,3)
            Wres[:,:,(end+1)-pmi,:] .= W[:,:,pmi,:]
        end
    else
        for pmi in axes(W,2)
            Wres[:,(end+1)-pmi,:] .= W[:,pmi,:]
        end
    end

    return Wres
end

"""
    h5to4D(filepath_distr)
    h5to4D(filepath_distr, rowmajor=false, verbose=false)

Load and read a .h5/.hdf5 file containing the data necessary to construct a 4D fast-ion distribution, with dimensions (E,p,R,z).
The keyword arguments include:
    rowmajor - If the .h5/.hdf5 file was created and saved using a row-major ordering programming language (e.g. C/C++/NumPy in Python,
               for more info please study https://en.wikipedia.org/wiki/Row-_and_column-major_order), the 'rowmajor' keyword argument 
               should be set to true.
    verbose - If set to true, the function will talk a lot!
"""
function h5to4D(filepath_distr::AbstractString; rowmajor::Bool = false, verbose::Bool = false)

    if verbose
        println("Loading the 4D distribution from .hdf5 file... ")
    end
    myfile = h5open(filepath_distr,"r")
    if haskey(myfile,"F_ps")
        verbose && println("Found F_ps key in .hdf5 file.")
        f_data = read(myfile["F_ps"])
    elseif haskey(myfile,"f")
        verbose && println("Found f key in .hdf5 file.")
        f_data = read(myfile["f"])
    elseif haskey(myfile,"F_EpRz")
        verbose && println("Found F_EpRz key in .hdf5 file.")
        f_data = read(myfile["F_EpRz"])
    else
        error("Fast-ion distribution .hdf5 file did not have any expected keys for the distribution (F_ps, f or F_EpRz).")
    end
    if haskey(myfile,"energy")
        energy = read(myfile["energy"])
    elseif haskey(myfile,"E_array")
        energy = read(myfile["E_array"])
    else
        error("Fast-ion distribution .hdf5 file did not have any expected keys for the energy grid points (energy or E_array).")
    end
    if haskey(myfile,"pitch")
        pitch = read(myfile["pitch"])
    elseif haskey(myfile,"p_array")
        pitch = read(myfile["p_array"])
    else
        error("Fast-ion distribution .hdf5 file did not have any expected keys for the energy grid points (pitch or p_array).")
    end
    if haskey(myfile,"R")
        R = read(myfile["R"])
    elseif haskey(myfile,"R_array")
        R = read(myfile["R_array"])
    else
        error("Fast-ion distribution .hdf5 file did not have any expected keys for the R grid points (R or R_array).")
    end
    if haskey(myfile,"z")
        verbose && println("Found z key in .hdf5 file.")
        z = read(myfile["z"])
    elseif haskey(myfile,"Z")
        verbose && println("Found Z key in .hdf5 file.")
        z = read(myfile["Z"])
    elseif haskey(myfile,"z_array")
        z = read(myfile["z_array"])
    elseif haskey(myfile,"Z_array")
        z = read(myfile["Z_array"])
    else
        error("Fast-ion distribution .hdf5 file did not have any expected keys for the z grid points (z, Z, z_array or Z_array).")
    end
    close(myfile)

    if rowmajor
        verbose && println("Fast-ion distribution data is permutated. Solving... ")
        f = zeros(size(f_data,4),size(f_data,3),size(f_data,2),size(f_data,1))
        for i in axes(f,1)
            for j in axes(f,2)
                for k in axes(f,3)
                    for l in axes(f,4)
                        f[i,j,k,l] = f_data[l,k,j,i]
                    end
                end
            end
        end
    else
        f = f_data # If .h5 file was not created with Python, the 4D array is not permutated
    end

    return f, energy, pitch, R, z
end
hdf5to4D = h5to4D # Function name synonym

"""
    jld2tohdf5(filepath_jld2)
    jld2tohdf5(filepath_jld2; verbose=false)

Convert a .jld2 file to a .hdf5 file.
"""
function jld2tohdf5(filepath_jld2::String; verbose::Bool=false)
    filepath_hdf5 = reduce(*, split(filepath_jld2,".")[1:end-1])*".hdf5"
    myfile_jld2 = jldopen(filepath_jld2,false,false,false,IOStream)
    myfile_hdf5 = h5open(filepath_hdf5,"w")

    for key in keys(myfile_jld2)
        verbose && println("Writing data: "*key)
        data = myfile_jld2[key]
        verbose && print("          ↳ Type: $(typeof(data))")
        if typeof(data) <: Dict # If the data is a dictionary
            verbose && println("   Not ok!")
            @warn "Part of the data in "*filepath_jld2*" is a dictionary. The OWCF currently does not support saving dictionaries in .hdf5 file format. The "*key*" data has therefore been omitted."
        elseif typeof(data) <: StepRangeLen
            verbose && println("   Ok! (Type StepRangeLen converted to type Array via collect())")
            data = collect(data)
            write(myfile_hdf5,key,data)
        else
            verbose && println("   Ok!")
            write(myfile_hdf5,key,data)
        end
    end
    close(myfile_jld2)
    close(myfile_hdf5)
    verbose && println("The .jld2 file has been re-written as a .hdf5 file at "*filepath_hdf5)
    return filepath_hdf5
end

"""
    hdf5tojld2(filepath_hdf5)
    hdf5tojld2(filepath_hdf5; verbose=false)

Convert a .hdf5 file to a .jld2 file.
"""
function hdf5tojld2(filepath_hdf5::String; verbose::Bool=false)
    filepath_jld2 = reduce(*, split(filepath_hdf5,".")[1:end-1])*".jld2"
    myfile_hdf5 = h5open(filepath_hdf5,"r")
    myfile_jld2 = jldopen(filepath_jld2,true,true,false,IOStream)

    for key in keys(myfile_hdf5)
        verbose && println("Writing data: "*key)
        data = myfile_hdf5[key]
        verbose && print("          ↳ Type: $(typeof(data))")
        verbose && println("   Ok!")
        write(myfile_jld2,key,data)
    end
    close(myfile_jld2)
    close(myfile_hdf5)
    verbose && println("The .hdf5 file has been re-written as a .jld2 file at "*filepath_jld2)
    return filepath_jld2
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
    if haskey(myfile,"energy")
        E_array = myfile["energy"]
    elseif haskey(myfile,"E_array")
        E_array = myfile["E_array"]
    else
        error("Fast-ion distribution .jld2 file did not have any expected keys for the energy grid points (energy or E_array).")
    end
    if haskey(myfile,"pitch")
        p_array = myfile["pitch"]
    elseif haskey(myfile,"p_array")
        p_array = myfile["p_array"]
    else
        error("Fast-ion distribution .jld2 file did not have any expected keys for the energy grid points (pitch or p_array).")
    end
    if haskey(myfile,"R")
        R_array = myfile["R"]
    elseif haskey(myfile,"R_array")
        R_array = myfile["R_array"]
    else
        error("Fast-ion distribution .jld2 file did not have any expected keys for the R grid points (R or R_array).")
    end
    if haskey(myfile,"z")
        verbose && println("Found z key in .jld2 file.")
        z_array = myfile["z"]
    elseif haskey(myfile,"Z")
        verbose && println("Found Z key in .jld2 file.")
        z_array = myfile["Z"]
    elseif haskey(myfile,"z_array")
        z_array = myfile["z_array"]
    elseif haskey(myfile,"Z_array")
        z_array = myfile["Z_array"]
    else
        error("Fast-ion distribution .jld2 file did not have any expected keys for the z grid points (z, Z, z_array or Z_array).")
    end
    close(myfile)

    return F_EpRz, E_array, p_array, R_array, z_array
end

"""
read_ncdf(filepath; wanted_keys=nothing)

Open and read the contents of a file in .cdf file format. The 'wanted_keys' keyword argument is a Vector of Strings that 
specifies the file keys to load. If no file keys have been specified via the 'wanted_keys' keyword argument, load all keys 
of the .cdf file, and return a dictionary with the .cdf file contents. The keys of the dictionary are the same as the keys 
of the .cdf file.
"""
function read_ncdf(filepath::String; wanted_keys=nothing)

    d = Dict()
    d["err"] = 1
    if isfile(filepath)
        d["err"] = 0
        NetCDF.open(filepath) do nc
            cdf_variables = nc.vars
            if !isnothing(wanted_keys)
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

"""
CDFto4D(filepath_distr, R_array, z_array)
CDFto4D(-||-; E_range=nothing, p_range=nothing, btipsign=-1, species=1, verbose=false, vverbose=false)

Read the TRANSP NUBEAM-generated data in the .cdf-file 'filepath_distr'. Use 'R_array' and 'z_array' as the (R,z) grid to interpolate the fast-ion data onto.
NUBEAM-generated fast-ion distributions are given on an irrregular grid spiraling outwards from the magnetic axis. CDFto4D() will 
interpolate that data onto a regular (E,p,R,z) grid.

- 'filepath_distr' is the filepath to the TRANSP NUBEAM-generated output file in .cdf file format 
- 'R_array' are the major radius (R) grid points onto which the fast-ion distribution will be interpolated. In meters
- 'z_array' are the vertical coordinate (z) grid points onto which the fast-ion distribution will be interpolated. In meters
- If 'E_range' and/or 'p_range' are not specified, use TRANSP-data ranges by default. E_range should be specified in keV. If specified, only use TRANSP-data 
within the 'E_range'/'p_range' ranges to create the fast-ion distribution f(E,p,R,z). - Union{Nothing,Tuple{Float64,Float64}}
- If 'btipsign' is not specified, assume 1 by default. 'btipsign' is the sign of the dot-product between the magnetic field and the plasma current.
- If 'species' is not specified, assume 1 by default. 'species' is the integer index of the fast-ion species in the TRANSP .cdf-file. Usually the first
fast-ion species '1' is wanted.
- If 'verbose' is not specified, assume false by default. If 'verbose' is set to true, read_nubeam() will talk a lot!
- If 'vverbose' is not specified, assume false by default. If 'vverbose' it set to true, read_nubeam() will talk a lot more (including printing for-loop information etc.)!

Return f(E,p,R,z) in keV^-1_m^-3 and R- and z-arrays in meters.
Original function written in Python by Luke Stagner as part of the FIDASIM code (https://github.com/D3DEnergetic/FIDASIM).
"""
function CDFto4D(filepath_distr::String, R_array::Union{T,Vector{T}} where {T<:Real}, z_array::Union{T,Vector{T}} where {T<:Real}; E_range::Union{Nothing,Tuple{T,T}} where {T<:Real}=nothing, p_range::Union{Nothing,Tuple{Float64,Float64}}=nothing, btipsign::Int64=1, species::Int64=1, verbose::Bool = false, vverbose::Bool = false)
    
    R_array = vcat(R_array) # Convert to Vector, if Real. Do nothing, if Vector
    z_array = vcat(z_array) # Convert to Vector, if Real. Do nothing, if Vector

    verbose && println("Converting (R,z) points from meters to centimeters (TRANSP default)... ")
    R_array = 100*R_array
    z_array = 100*z_array
    
    verbose && println("Creating (R,z) grid from R,z input... ")
    mygrid = rz_grid(minimum(R_array), maximum(R_array), length(R_array), minimum(z_array), maximum(z_array), length(z_array))

    verbose && println("Loading data from TRANSP .cdf file... ")
    species_var = "SPECIES_$(species)" # The fast-ion species identification pattern for TRANSP data
    sstr = read_ncdf(filepath_distr, wanted_keys=[species_var])[species_var] # Will return a string with the fast-ion species of the TRANSP shot
    TRANSP_data = read_ncdf(filepath_distr, wanted_keys=["TIME","R2D","Z2D","E_"*sstr,"A_"*sstr,"F_"*sstr,"RSURF","ZSURF","BMVOL"]) # Will return a struct with the requested TRANSP data

    ngrid = length(TRANSP_data["R2D"]) # The number of (spiral) grid points on which the TRANSP data is given

    if length(TRANSP_data["TIME"])>1
        TRANSP_timepoint = TRANSP_data["TIME"][1] # The tokamak shot timepoint of the TRANSP-data. In seconds
    else
        TRANSP_timepoint = TRANSP_data["TIME"] # The tokamak shot timepoint of the TRANSP-data. In seconds
    end

    TRANSP_R_vector = TRANSP_data["R2D"] # 1D array. R values in cm 
    TRANSP_z_vector = TRANSP_data["Z2D"] # 1D array. z values in cm
    TRANSP_R_mesh = TRANSP_data["RSURF"] # 2D array. R values in cm
    TRANSP_z_mesh = TRANSP_data["ZSURF"] # 2D array. z values in cm
    bmvol = TRANSP_data["BMVOL"] # Volume of elements at (spiral) grid points. In cm^-3
    p_vector = TRANSP_data["A_"*sstr] # 1D array
    E_vector = TRANSP_data["E_"*sstr]*(1.0e-3) # 1D_array. Energy values in keV (after multiplication with (1.0e-3) factor)
    F_Epbm = (TRANSP_data["F_"*sstr])*(1.0e3) # 3D array. size(F_Epbm)=(length(E_vector), length(p_vector), length(bmvol)). Distribution values in (keV)^-1 cm^-3
    F_Epbm = map(x-> x > 0.0 ? 0.5*x : 0.0, F_Epbm) #Multiply elements by 0.5 to convert to pitch instead of solid angle d_omega/4pi

    if btipsign < 0
        verbose && println("Flipping fast-ion distribution in pitch since sign(J ⋅ B)< 0... ")
        F_Epbm = reverse(F_Epbm,dims=2) # Reverse all elements in dimension of pitch, to account for plasma current and magnetic field pointing in different directions
    end
    if E_range==nothing # If an energy range has not been specified as input with the keyword arguments...
        E_range = (minimum(E_vector), maximum(E_vector))
    end
    if p_range==nothing # If a pitch range has not been specified as input with the keyword arguments...
        p_range = (minimum(p_vector), maximum(p_vector))
    end

    we = findall(x-> x >= E_range[1] && x <= E_range[2], E_vector) # Find the indices of all energies within specified energy range e_range
    wp = findall(x-> x >= p_range[1] && x <= p_range[2], p_vector) # Find the indices of all pitches within specified pitch range p_range
    E_vector = E_vector[we] # Trim energy vector accordingly
    nenergy = length(E_vector)
    p_vector = p_vector[wp] # Trim pitch vector accordingly
    npitch = length(p_vector)
    F_Epbm = F_Epbm[we,:,:] # Trim fast-ion distribution accordingly
    F_Epbm = F_Epbm[:,wp,:]
    dE = abs(E_vector[2] - E_vector[1]) # Assume uniform energy grid...
    dp = abs(p_vector[2] - p_vector[1]) # Assume uniform pitch grid
    emin, emax = clamp(minimum(E_vector) - 0.5*dE, 0.0, minimum(E_vector)), maximum(E_vector) + 0.5*dE # Make sure that lower energy bound is larger than zero
    pmin, pmax = clamp(minimum(p_vector) - 0.5*dp, -1.0, minimum(p_vector)), clamp(maximum(p_vector)+0.5*dp, maximum(p_vector), 1.0) # Make sure pitch boundaries are within [-1.0, 1.0]

    verbose && println("Energy min/max: $((emin, emax))")
    verbose && println("Pitch min/max: $((pmin, pmax))")

    # Query (R,z) grid data on which to interpolate on
    query_nR = mygrid.nr # Number of query R points
    query_nz = mygrid.nz
    query_R = mygrid.r
    query_z = mygrid.z
    query_R_mesh = mygrid.r2d
    query_z_mesh = mygrid.z2d
    dR = abs(query_R[2] - query_R[1])
    dz = abs(query_z[2] - query_z[1])

    F_bm = dropdims(sum(F_Epbm,dims=(1,2))*dE*dp,dims=(1,2)) # The fast-ion distribution as a function of the bm volume elements (units: cm^-3)
    ntot = sum(F_bm .*bmvol) # Multiply F_bm with the bm volume elements to get the total number of fast ions. Units: (cm^-3)*(cm^3) = #fast ions
    verbose && println("Ntotal in phase space: $(ntot)")

    verbose && println("Creating the Delaunay tessellation... ")
    tri = DelaunayTessellation(length(TRANSP_R_vector)) # Instantiate a Delaunay tessellation of approximately length(TRANSP_R_vector)
    mnr = minimum([minimum(TRANSP_R_vector),minimum(query_R)]) # The definite minimum of every possible R coordinate in the algorithm
    mxr = maximum([maximum(TRANSP_R_vector),maximum(query_R)]) # The definite maximum of every possible R coordinate in the algorithm
    mnz = minimum([minimum(TRANSP_z_vector),minimum(query_z)]) # The definite minimum of every possible z coordinate in the algorithm
    mxz = maximum([maximum(TRANSP_z_vector),maximum(query_z)]) # The definite maximum of every possible z coordinate in the algorithm
    r_tri(x) = min_coord + (max_coord-min_coord)*(x-mnr)/(mxr-mnr) # DelaunayTesselation expects all points to be within [min_coord,max_coord] (given by VoronoiDelaunay.min_coord etc). We therefore need a conversion function.
    z_tri(y) = min_coord + (max_coord-min_coord)*(y-mnz)/(mxz-mnz) # DelaunayTesselation expects all points to be within [min_coord,max_coord] (given by VoronoiDelaunay.min_coord etc). We therefore need a conversion function.
    r_tri_inv(x_tri) = mnr + (mxr-mnr)*(x_tri-min_coord)/(max_coord-min_coord) # We need a function to convert back to normal R values
    z_tri_inv(y_tri) = mnz + (mxz-mnz)*(y_tri-min_coord)/(max_coord-min_coord) # We need a function to convert back to normal z values

    TRANSP_spiral_rz_iterator_tri = zip(r_tri.(TRANSP_R_vector), z_tri.(TRANSP_z_vector)) # Create an iterator with all the TRANSP (R,z) coordinates as tuples
    a = Point2D[Point(rr_tri, zz_tri) for (rr_tri,zz_tri) in TRANSP_spiral_rz_iterator_tri] # Put them all in a Point() object within a Point2D array
    push!(tri, a) # Feed all points to the Delaunay tessellation

    verbose && println("Mapping all triangles to original vertex indices... ")
    vertices = Dict{String,Tuple{Int64,Int64,Int64}}() # A hashmap. The hashes are triangle IDs and the values are tuples with vertices indices in the same order as the a,b and c fields of the DelaunayTriangle object returned by locate(tri, point). The DelaunayTesselation struct provided no way to keep track of what indices of the TRANSP_R_vector, TRANSP_z_vector inputs that make up the vertices of the triangles in the tessellation... So I had to do it myself!
    count = 1
    for t in tri # For triangles in the tessellation
        vverbose && println("Triangle $(count)")
        ind_a = findall(x-> x==(getx(geta(t)),gety(geta(t))), collect(TRANSP_spiral_rz_iterator_tri)) # Find the (R,z) data index of the first vertex to triangle t
        ind_b = findall(x-> x==(getx(getb(t)),gety(getb(t))), collect(TRANSP_spiral_rz_iterator_tri)) # Find the (R,z) data index of the second vertex to triangle t
        ind_c = findall(x-> x==(getx(getc(t)),gety(getc(t))), collect(TRANSP_spiral_rz_iterator_tri)) # Find the (R,z) data index of the third vertex to triangle t

        vverbose && println("ind_a: $(ind_a)    ind_b: $(ind_b)     ind_c: $(ind_c)")

        triangle_ID = "$(t._neighbour_a)_$(t._neighbour_b)_$(t._neighbour_c)" # No other triangle in the tessellation will have this ID. Thus making it easily look-up-able later
        vertices[triangle_ID] = (first(ind_a), first(ind_b), first(ind_c))
        vverbose && println("Triangle ID: "*triangle_ID)
        vverbose && println()
        count += 1
    end

    # Test hashmap on test point
    verbose && println("Creating nearest neighbour (NN) tree for extrapolation of points outside the tessellation... ")
    brutetree = BruteTree(hcat(TRANSP_R_vector,TRANSP_z_vector)') # There are also other tree types. However, BruteTree was deemed good enough

    F_Rz = zeros(query_nR,query_nz)
    F_EpRz = zeros(nenergy,npitch,query_nR,query_nz)
    count = 1
    verbose && println("Interpolating fast-ion distribution onto all query points... ")
    for (ir,rr) in enumerate(query_R), (iz,zz) in enumerate(query_z)
        vverbose && println("$(count)/$(length(query_R)*length(query_z)): (R,z)=($(rr),$(zz))")
        t = locate(tri, Point(r_tri(rr), z_tri(zz))) # Try to locate the query point inside of the tessellation

        if isexternal(t) == true # If the query point was outside of the tessellation...
            idx, dist = nn(brutetree, [rr, zz]) # Use brute-force nearest-neighbour extrapolation instead! Notice that here we don't need to r_tri() and z_tri() functions
            F_Rz[ir,iz] = F_bm[idx]
            F_EpRz[:,:,ir,iz] = F_Epbm[:,:,idx]
        else
            tID = "$(t._neighbour_a)_$(t._neighbour_b)_$(t._neighbour_c)" # Contruct the triangle ID from the order of the neighbours
            ia, ib, ic = vertices[tID] # Look up the original indices of the vertices
            x = [r_tri(rr), z_tri(zz)] # Query point in tri coordinates
            xa = [getx(geta(t)), gety(geta(t))] # Vertex a of triangle t
            xb = [getx(getb(t)), gety(getb(t))] # Vertex b of triangle t
            xc = [getx(getc(t)), gety(getc(t))] # Vertex c of triangle t

            μ = [(1/sum(abs2,x-xa)),(1/sum(abs2,x-xb)),(1/sum(abs2,x-xc))] 
            μ = (1/sum(μ)) .* μ # Barycentric weights. sum(μ)=1 must hold

            # Barycentric interpolation. All components of μ sum to one. 'Barycentric' simply means 'center-of-mass-like'. https://en.wikipedia.org/wiki/Barycentric_coordinate_system 
            F_Rz[ir,iz] = μ[1]*F_bm[ia]+μ[2]*F_bm[ib]+μ[3]*F_bm[ic]
            F_EpRz[:,:,ir,iz] = μ[1]*F_Epbm[:,:,ia]+μ[2]*F_Epbm[:,:,ib]+μ[3]*F_Epbm[:,:,ic]
        end
        count += 1
    end 

    F_Rz = map(x-> x<0.0 ? 0.0 : x, F_Rz) # Any negative values? Set them to zero instead!

    # Correct for points outside of seperatrix. They should all be zero
    verbose && println("Identifying query points outside of the separatrix... ")
    R_maxis = mean(TRANSP_R_mesh[:,1]) # The R coordinate of the magnetic axis
    z_maxis = mean(TRANSP_z_mesh[:,1])
    R_sep = TRANSP_R_mesh[:,end] # The R coordintes of the separatrix
    z_sep = TRANSP_z_mesh[:,end]

    R_bdry = R_sep .- R_maxis # The R coordinates of the separatrix, relative to the magnetic axis
    z_bdry = z_sep .- z_maxis

    r_bdry = sqrt.(R_bdry.^2 .+ z_bdry.^2) # The r coordinates (local toroidal coordinate system. http://fusionwiki.ciemat.es/wiki/Toroidal_coordinates)
    theta_bdry = atan.(z_bdry,R_bdry) # The θ coordinates from arctan of (z-z_axis)/(R-R_axis) values
    theta_bdry = map(x-> x < 0.0 ? (x + 2*pi) : x,theta_bdry) # Ensure θ coordinates within [0,2π]
    w = sortperm(theta_bdry) # Sort theta_bdry but return only the indices that would give the sorted order
    theta_bdry = theta_bdry[w]
    r_bdry = r_bdry[w]
    w = unique(i -> theta_bdry[i], 1:length(theta_bdry)) # Keep only unique elements
    theta_bdry = theta_bdry[w]
    r_bdry = r_bdry[w]

    itp = Interpolations.interpolate((theta_bdry,), r_bdry, Gridded(Linear())) # Create an interpolation object, to be able to interpolate new values of r from input theta values
    etpf = extrapolate(itp, Interpolations.Flat()) # If outside of domain, just return flat (constant) values (flat=last known data value before extrapolation domain)


    R_query_pts = query_R_mesh .- R_maxis # Same procedure as for data points, but for query points
    z_query_pts = query_z_mesh .- z_maxis
    r_query_pts = sqrt.(R_query_pts.^2 .+ z_query_pts.^2)
    theta_query_pts = atan.(z_query_pts,R_query_pts)
    theta_query_pts = map(x-> x < 0.0 ? (x + 2*pi) : x, theta_query_pts) # [0,2pi]

    # For all query theta angles, find the cooresponding r points on the separatrix
    r_query_bdry = zeros(size(theta_query_pts))
    for (i,theta_query_pt) in enumerate(theta_query_pts)
        r_query_bdry[i] = etpf(theta_query_pt)
    end

    w = findall(x-> x >= 0.0, (r_query_pts - r_query_bdry) .- 2) # If the query points are more than 2 cm outside of the separatrix...
    F_Rz[w] .= 0.0 # ...make sure that the fast-ion distribution is zero there and outwards (> 2cm outside of the separatrix)
    F_EpRz[:,:,w] .= 0.0 # ...make sure that the fast-ion distribution is zero there and outwards (> 2cm outside of the separatrix)

    # Enforce correct normalization (sum(distribution) = #total number of fast ions) for interpolated fast-ion distributions
    verbose && println("Enforcing correct normalization for interpolated fast-ion distribution (sum(f)=#total number of fast ions)... ")
    ntot_denf = 2*pi*dR*dz*sum(reshape(query_R,(length(query_R),1)) .*F_Rz)
    F_Rz = (ntot/ntot_denf) .* F_Rz 
    ntot_fbm = (2*pi*dE*dp*dR*dz)*sum(reshape(query_R,(1,1,length(query_R),1)) .*F_EpRz)
    F_EpRz = (ntot/ntot_fbm) .*F_EpRz

    verbose && println("Returning (E,p,R,z) fast-ion distribution in units [keV^-1 m^-3]")
    return (1.0e6) .*F_EpRz, E_vector, p_vector, query_R*0.01, query_z*0.01
end

"""
    h5toSampleReady(filepath_distr, interp)
    h5toSampleReady(-||-, verbose=false)

Load the fast-ion distribution from a .h5/.hdf5 file and forward to function toSampleReady().

--- Input:
filepath_distr - The fast-ion distribution in .h5 file format - String
The keyword arguments include: 
    rowmajor - If the .h5/.hdf5 file was created and saved using a row-major ordering programming language (e.g. C/C++/NumPy in Python,
               for more info please study https://en.wikipedia.org/wiki/Row-_and_column-major_order), the 'rowmajor' keyword argument 
               should be set to true.
    verbose - If set to true, you will get a talkative function! - Bool
"""
function h5toSampleReady(filepath_distr::AbstractString; rowmajor::Bool = false, verbose::Bool = false, kwargs...)
    verbose && println("Loading fast-ion distribution from .h5 file... ")
    F_EpRz, E_array, p_array, R_array, z_array = h5to4D(filepath_distr, rowmajor=rowmajor, verbose=verbose)

    return toSampleReady(F_EpRz, E_array, p_array, R_array, z_array; verbose=verbose, kwargs...)
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
    toSampleReady(-||-, nE_ps=202, np_ps=57, nR_ps=62, nz_ps=64, slices_of_interest=[nothing,nothing,nothing,nothing], verbose=true)

Interpolate the fast-ion distribution and grid arrays if specified, extract the E-, p-, R- and/or z-planes of interest 
if specified, then compute the necessary quantities needed to be able to sample from the fast-ion distribution with the spaghettification
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
slices_of_interest - Specify Float64 instead of nothing to extract fast-ion distribution for specific E-, p-, R- and/or z-planes - Vector{Union{Nothing,Float64}}
verbose - If set to true, you will get a talkative function indeed - Bool
"""
function toSampleReady(F_EpRz::Array{Float64,4}, E_array::Array{Float64,1}, p_array::Array{Float64,1}, R_array::Array{Float64,1}, z_array::Array{Float64,1}; interp=false, nE_ps::Int64=100, np_ps::Int64=51, nR_ps::Int64=60, nz_ps::Int64=61, slices_of_interest=[nothing,nothing,nothing,nothing], verbose::Bool = false)
    if interp
        energyq = range(E_array[1], E_array[end], length=nE_ps) # Set the energy query points
        pitchq = range(p_array[1], p_array[end], length=np_ps) # Etc
        Rq = range(R_array[1], R_array[end], length=nR_ps) # Etc
        zq = range(z_array[1], z_array[end], length=nz_ps) # Etc
        verbose && println("toSampleReady(): Interpolating Julia fast-ion distribution... ")
        F_EpRz = interpFps(F_EpRz, E_array, p_array, R_array, z_array, energyq, pitchq, Rq, zq) # Interpolate
        E_array = energyq # The energy query points are now the energy points
        p_array = pitchq # Etc
        R_array = Rq
        z_array = zq
    end

    verbose && println("toSampleReady(): Checking for E-, p-, R- and/or z-slices of interest... ")
    super_array = [E_array, p_array, R_array, z_array] # Put all (E,p,R,z) vectors in a supervector
    string_array = ["energy","pitch","R","z"]
    inds_of_interest = findall(x-> x!==nothing,slices_of_interest) # Are there any specific E-, p-, R- and/or z-slices of interest?
    if length(inds_of_interest)>0
        for ind in inds_of_interest
            verbose && println("toSampleReady(): Found "*string_array[ind]*" slice of interest. Processing... ")
            s = slices_of_interest[ind] # One Float64 value
            s_array = super_array[ind] # One array

            if s>maximum(s_array)
                verbose && (@warn "toSampleReady(): Requested "*string_array[ind]*" value ($(s)) is larger than maximum("*string_array[ind]*"_array) ($(maximum(s_array))). maximum("*string_array[ind]*"_array) will be used instead (extrapolation not possible)")
                s = maximum(s_array)
            elseif s<minimum(s_array)
                verbose && (@warn "toSampleReady(): Requested "*string_array[ind]*" value ($(s)) is smaller than minimum("*string_array[ind]*"_array) ($(minimum(s_array))). minimum("*string_array[ind]*"_array) will be used instead (extrapolation not possible)")
                s = minimum(s_array)
            else
                verbose && println("toSampleReady(): Requested "*string_array[ind]*" value ($(s)) inside of fast-ion distribution grid $(extrema(s_array)).")
            end
        end
        Eq = (1 in inds_of_interest) ? [slices_of_interest[1]] : E_array
        pq = (2 in inds_of_interest) ? [slices_of_interest[2]] : p_array
        Rq = (3 in inds_of_interest) ? [slices_of_interest[3]] : R_array
        zq = (4 in inds_of_interest) ? [slices_of_interest[4]] : z_array
        F_o_interest = interpFps(F_EpRz,E_array,p_array,R_array,z_array,Eq,pq,Rq,zq)
        E_array = Eq # The energy query points are now the energy points
        p_array = pq # Etc
        R_array = Rq
        z_array = zq
    else
        F_o_interest = F_EpRz
    end
    fr = F_o_interest.*reshape(R_array,(1,1,length(R_array),1))
    verbose && println("toSampleReady(): Computing 4D vols... ")
    dvols = get4DVols(E_array, p_array, R_array, z_array)
    nfast = sum(fr.*dvols)*(2*pi)

    # Checking if units of R,z is meters or centimeters. Correct to meters if units is centimeters
    # PLEASE NOTE! IT IS VERY IMPORTANT TO DO THIS STEP AFTER YOU HAVE COMPUTED fr, dvols and nfast
    # OTHERWISE, THE TOTAL NUMBER OF FAST IONS WILL BE WRONG, BECAUSE THE ORIGINAL F_EpRz DATA WAS
    # IN UNITS OF PER CENTIMETER
    if maximum(R_array)>100.0 # Assume no tokamak has a major radius larger than 100 meters...
        verbose && println("toSampleReady(): Converting R- and z-arrays from centimeters to meters... ")
        R_array = R_array./100.0
        z_array = z_array./100.0
    end
    dE_array = (length(E_array)>1) ? vcat(abs.(diff(E_array)),abs(E_array[end]-E_array[end-1])) : [0.0]
    dp_array = (length(p_array)>1) ? vcat(abs.(diff(p_array)),abs(p_array[end]-p_array[end-1])) : [0.0]
    dR_array = (length(R_array)>1) ? vcat(abs.(diff(R_array)),abs(R_array[end]-R_array[end-1])) : [0.0]
    dz_array = (length(z_array)>1) ? vcat(abs.(diff(z_array)),abs(z_array[end]-z_array[end-1])) : [0.0]

    dims = size(fr) # Tuple
    verbose && println("toSampleReady(): Computing cumulative sum vector... ")
    frdvols_cumsum_vector = cumsum(vec(fr .*dvols)) # Vector. Reaaally long vector
    subs = CartesianIndices(dims) # 4D matrix

    return frdvols_cumsum_vector, subs, E_array, p_array, R_array, z_array, dE_array, dp_array, dR_array, dz_array, nfast

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

Calculate and return 4D-arrays of the diffs of particle-space coordinates. If length>2, assume edge diff to be same as next-to-edge diff.
If length==1, use diff=1. If length==0, throw error.
"""
function get4DDiffs(E::AbstractVector, p::AbstractVector, R::AbstractVector, z::AbstractVector)
    if length(R)>1
        dR = vcat(abs.(diff(R)),abs(R[end]-R[end-1]))
        dR = reshape(dR,1,1,length(dR),1)
    elseif length(R)==1
        dR = reshape([0.0],(1,1,1,1))
    else
        error("R input was of 0 length. Please correct and re-try.")
    end
    dR4D = repeat(dR,length(E),length(p),1,length(z))

    if length(z)>1
        dz = vcat(abs.(diff(z)),abs(z[end]-z[end-1]))
        dz = reshape(dz,1,1,1,length(dz))
    elseif length(z)==1
        dz = reshape([0.0],(1,1,1,1))
    else
        error("z input was of 0 length. Please correct and re-try.")
    end
    dz4D = repeat(dz,length(E),length(p),length(R),1)

    if length(E)>1
        dE = vcat(abs.(diff(E)),abs(E[end]-E[end-1]))
        dE = reshape(dE,length(dE),1,1,1)
    elseif length(E)==1
        dE = reshape([0.0],(1,1,1,1))
    else
        error("E input was of 0 length. Please correct and re-try.")
    end
    dE4D = repeat(dE,1,length(p),length(R),length(z))

    if length(p)>1
        dp = vcat(abs.(diff(p)),abs(p[end]-p[end-1]))
        dp = reshape(dp,1,length(dp),1,1)
    elseif length(p)==1
        dp = reshape([0.0],(1,1,1,1))
    else
        error("p input was of 0 length. Please correct and re-try.")
    end
    dp4D = repeat(dp,length(E),1,length(R),length(z))

    return dE4D, dp4D, dR4D, dz4D
end

"""
    get4DVols(E, p, R, z)

This function will calculate the volume of all the hyper-voxels pertaining to the 4D grid. It assumes the hyper-voxels pertaining to the
upper-end (edge) grid-points have the same volumes as the hyper-voxels just inside of them. It return a 4D array, containing all the hyper-voxel
volumes. The 4D array will have size()=(length(E), length(p), length(R), length(z)). The function assumes a rectangular 4D grid.
If any of the dimensions only has a single grid point, compute the resulting lower-dimensional voxel volume (3D: volume. 2D: area. 1D: length). 
"""
function get4DVols(E::AbstractVector, p::AbstractVector, R::AbstractVector, z::AbstractVector)

    # Safety-check to ensure vectors
    if !(1==length(size(E))==length(size(p))==length(size(R))==length(size(z)))
        throw(ArgumentError("E, p, R, z inputs are not all vectors. Please correct and re-try."))
    end

    dE4D, dp4D, dR4D, dz4D = get4DDiffs(E, p, R, z)
    if sum(dE4D)==0.0
        dE = reshape([1.0],(1,1,1,1))
        dE4D = repeat(dE,1,length(p),length(R),length(z))
    end
    if sum(dp4D)==0.0
        dp = reshape([1.0],(1,1,1,1))
        dp4D = repeat(dp,length(E),1,length(R),length(z))
    end
    if sum(dR4D)==0.0
        dR = reshape([1.0],(1,1,1,1))
        dR4D = repeat(dR,length(E),length(p),1,length(z))
    end
    if sum(dz4D)==0.0
        dz = reshape([1.0],(1,1,1,1))
        dz4D = repeat(dz,length(E),length(p),length(R),1)
    end
    dvols = dE4D .*dp4D .*dR4D .*dz4D
    return dvols
end

"""
    OWCF_sample_f(fr,dvols,energy,pitch,R,z)
    OWCF_sample_f(-||-;n=100_000)

This function is the OWCF version of the original function sample_f which already exists in the OrbitTomography.jl package.
This version has the ability to take non-equidistant 4D grid-points into consideration. w is energy, x is pitch, y is R and z is z.
"""
function OWCF_sample_f(fr::Array{T,N}, dvols, w::AbstractVector, x::AbstractVector, y::AbstractVector, z::AbstractVector; n::Int64=100_000) where {T,N}

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
    class2int(M,o)
    class2int(-||-; sigma=0, plot=true, distinguishLost=false, distinguishIncomplete=false)

Take an orbit, examine its class (stagnation, trapped, co-passing etc) and return 
the appropriate integer. The integers are as follows
1 - stagnation
2 - trapped
3 - co-passing
4 - counter-passing
5 - potato
9 - invalid

Depending on the keyword arguments, give orbit appropriate integer
- If sigma==0, all integers are available (default). If sigma==1 and if orbit is counter-current,
  give 9. If sigma==-1 and if orbit is co-current, give 9.
- If plot, give counter-stagnation orbits integer 8. Else, give 6
- If distinguishLost, give lost orbits integer 7. Else, give 9
- If distinguishIncomplete and plot, give incomplete orbits integer 6. Otherwise, give 9

Keyword argument 'plot' is meant to ensure correct color map when using output for OWCF apps.
If 'plot' is set to false, then you can simply feed all integers into a 6 element array and 
have the valid orbits (stagnation, trapped, co-passing, counter-passing, potato and counter-stagnation) in bins 1 to 6.
"""
function class2int(M::AbstractEquilibrium, o::GuidingCenterOrbits.Orbit; sigma::Int64=0, plot::Bool=true, distinguishLost::Bool=false, distinguishIncomplete::Bool=false)
    if !(sigma==0 || sigma==1 || sigma==-1)
        error("Keyword argument 'sigma' not equal to 0, 1 or -1. Please correct and re-try.")
    end
    if hasproperty(o.coordinate,:r) # If the orbit has an (E,pm,Rm) coordinate
        Rm = o.coordinate.r
    else # It must have an (E,mu,Pphi) coordinate
        if isempty(o.path.r)
            @warn "Orbit object input does not have an attached OrbitPath object with non-empty (R,z) path. Cannot uniquely infer integer identification number. For stagnation orbits, output might be erronous."
            Rm = magnetic_axis(M)[1] # Just put it to magnetic axis, since we cannot know without computing the mu-contour
        end
        Rm = maximum(o.path.r)
    end
    if (o.class == :lost) && distinguishLost
        return 7
    elseif (o.class == :incomplete) && distinguishIncomplete && plot
        return 6
    elseif (o.class == :trapped) && !(sigma==-1)
        return 2
    elseif (o.class == :co_passing) && !(sigma==-1)
        return 3
    elseif (o.class == :stagnation) && (Rm>=magnetic_axis(M)[1]) && !(sigma==-1) # Regular stagnation orbit
        return 1
    elseif (o.class == :stagnation) && (Rm<magnetic_axis(M)[1]) && plot && !(sigma==1) # Counterstagnation orbit
        return 8
    elseif (o.class == :stagnation) && (Rm<magnetic_axis(M)[1]) && !plot && !(sigma==1) # Counterstagnation orbit, but treat as normal stagnation orbit
        return 6
    elseif (o.class == :potato) && !(sigma==-1)
        return 5
    elseif (o.class == :ctr_passing) && !(sigma==1)
        return 4
    else
        return 9
    end
end

"""
    OWCF_map_orbits(og, f)
    OWCF_map_orbits(-||-;weights=false)

This function is the OWCF version of the original function map_orbits which already exists in the OrbitTomography.jl package.
It takes non-equidistant 3D grid-points into account, as well as weights (then do not divide by orbit-space volume).
"""
function OWCF_map_orbits(og::OrbitGrid, f::Vector; weights::Bool=false)
    if length(og.counts) != length(f)
        throw(ArgumentError("Incompatible sizes for input og and f"))
    end

    if !(weights)
        dorb = get_orbel_volume(og) # Get the volumes of the voxels of the orbit grid
        return [i == 0 ? zero(f[1]) : f[i]/(og.counts[i]*dorb[i]) for i in og.orbit_index] # This line is complicated...
    else
        return [i == 0 ? zero(f[1]) : f[i]/(og.counts[i]*1.0) for i in og.orbit_index] # This line is complicated...
    end
end

OWCF_map_orbits(og::OrbitGrid, f::Vector, dummy_bool::Bool; kwargs...) = OWCF_map_orbits(og, f; kwargs...) # For backwards compatibility

"""
    ps2os_streamlined(filepath_distr, filepath_equil, og; numOsamples)
    ps2os_streamlined(-||-; numOsamples, verbose=false, kwargs...)

This function is a streamlined version of sampling particle-space (PS) and converting those samples to orbit-space (OS).
It allows transformation from (E,p,R,z) to (E,pm,Rm) directly from filepath_distr, using a magnetic equilibrium found in the filepath_equil file.
numOsamples sets the number of Monte-Carlo samples for the transformation. Verbose switches function print statements on/off.
"""
function ps2os_streamlined(filepath_distr::AbstractString, filepath_equil::AbstractString, og::OrbitGrid; rowmajor::Bool=false, numOsamples::Int64, verbose::Bool=false, kwargs...)

    if verbose
        println("")
        println("Loading f, energy, pitch, R and z from file... ")
    end
    f, energy, pitch, R, z = h5to4D(filepath_distr; rowmajor=rowmajor, verbose = verbose) # Load fast-ion distribution

    return ps2os_streamlined(f,energy,pitch,R,z,filepath_equil, og; numOsamples=numOsamples, verbose=verbose, kwargs...)
end

"""
ps2os_streamlined(F_EpRz, energy, pitch, R, z, filepath_equil, og)
ps2os_streamlined(-||-;numOsamples, verbose=false, distr_dim = [], sign_o_pitch_wrt_B=false, FI_species = "D", distributed=true, nbatch = 1_000_000, clockwise_phi, kwargs...)

Continuation of the ps2os_streamlined function above. F_EpRz is f(E,p,R,z). Energy is a vector with the energy grid points (keV). Pitch is a vector with the pitch grid points.
R is a vector with the major radius grid points (meters). z is a vector with the vertical coordinate grid points (meters). filepath_equil is a string pointing to the filepath 
of a magnetic equilibrium file (either .eqdsk file or output from extra/compSolovev.jl). og is an orbit grid, computed with the orbit_grid() function (OrbitTomography.jl Julia package). 
numOsamples is a necessary keyword argument that sets the number of Monte-Carlo samples for the transformation. verbose is a bool switch that turns function print statements on/off.
distr_dim is a list of exactly 4 elements. If specified, the fast-ion distribution will be interpolated to have this size in each of the 4 dimensions. For example, 
distr_dim=[50,40,30,35] will lead to that f(E,p,R,z) has size 50x40x30x35. sign_o_pitch_wrt_B is a bool switch that f(E,p,R,z)->f(E,-p,R,z). clockwise_phi is a bool switch for magnetic 
equilibria. Most magnetic equilibria have their cylindrical coordinate phi-direction in the anti-clockwise direction, tokamak viewed from above. However, some don't. For these, 
set clockwise_phi to true.
"""
function ps2os_streamlined(F_EpRz::Array{Float64,4}, energy::AbstractVector, pitch::AbstractVector, R::AbstractVector, z::AbstractVector, filepath_equil::AbstractString, 
                           og::OrbitGrid; 
                           numOsamples::Int64, verbose::Bool=false, distr_dim = [], sign_o_pitch_wrt_B::Bool=false, clockwise_phi::Bool=false, kwargs...)

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

    # This operation will be necessary if the loaded fast-ion distribution defines the sign of the pitch(p) w.r.t. the magnetic field (B)
    # and NOT the plasma current (J). That is, if the following three statements can be simultaneously true:
    #   1. p > 0
    #   2. sign(dot(v,B))>0 
    #   3. sign(dot(v,J))<0
    # or, equivalently if the following three statements can be simultaneously true:
    #   4. p < 0
    #   5. sign(dot(v,B))<0 
    #   6. sign(dot(v,J))>0.
    # This is because the OWCF defines the pitch to be positive w.r.t. the plasma current. I.e.:
    #   7. p > 0
    #   8. sign(dot(v,B))<0 
    #   9. sign(dot(v,J))>0
    # or
    #   10. p < 0
    #   11. sign(dot(v,B))>0 
    #   12. sign(dot(v,J))<0.
    # In short, some codes define the sign of the pitch via equalities 1-6, while the OWCF uses equalities 7-12.
    # This convention is not possible to deduce from the data structures alone, hence the information needs to be provided manually.
    if sign_o_pitch_wrt_B && (sign(jdotb)<0)
        verbose && println("sign_o_pitch_wrt_B is true and sign(dot(J,B))<0. Flipping pitch direction...")
        F_EpRz = reverse(F_EpRz,dims=2)
        pitch = -reverse(pitch) # 
    end

    return ps2os(M, wall, F_EpRz, energy, pitch, R, z, og; numOsamples=numOsamples, verbose=verbose, kwargs...)
end

"""
    ps2os(M, wall, F_EpRz, energy, pitch, R, z, og)
    ps2os(-||-; distributed=true, FI_species = "D", nbatch=100_000, numOsamples=1_000_000, performance=true, progress_file_name="ps2os", saveProgress=true, verbose=false, kwargs...)

Take the (E,p,R,z) fast-ion distribution and use Monte-Carlo sampling to transform it to (E,pm,Rm) orbit space. The keyword arguments are as follows:
- distributed - If true, multi-core processing will be used when MC sampling - Bool
- FI_species - The fast-ion particle species. Please see OWCF/misc/species_func.jl for more info - String
- nbatch - The size of the batch in which to perform MC sampling. Smaller batch, less RAM required but (likely) longer computational time - Int64
- numOsamples - The number of Monte-Carlo samples - Int64
- performance - If true, performance sampling will be used. Normal sampling will be deprecated in newer version of the OWCF - Bool
- progress_file_name - If saveProgress is set to true (see below), a progress file named progress_[progress_file_name].jld2 is saved, to ensure a warm start if 
                       the ps2os() progress is terminated prematurely. Defaults to "ps2os" - String
- save_progress - If true, a file named "progress_[progress_file_name].jld2" will be saved to enable a warm start of the ps2os() progress, if the process is 
                  terminated prematurely - Bool
- verbose - If set to true, the function will talk a lot! - Bool
"""
function ps2os(M::AbstractEquilibrium, wall::Boundary, F_EpRz::Array{Float64,4}, energy::AbstractVector, pitch::AbstractVector, R::AbstractVector, z::AbstractVector, og::OrbitGrid; 
               distributed::Bool=true, FI_species = "D", nbatch::Int64 = 100_000, numOsamples::Int64=1_000_000, performance::Bool=true, progress_file_name::String="ps2os",
               saveProgress::Bool=true, verbose::Bool=false, kwargs...)

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
    # Handle re-start of ps2os-transformation process, if terminated prematurely
    if isfile("progress_$(progress_file_name).jld2")
        myfile = jldopen("progress_$(progress_file_name).jld2",false,false,false,IOStream)
        numOsamples_sofar = deepcopy(myfile["numOsamples"])
        result_sofar = deepcopy(myfile["F_os"])
        class_distr_sofar = deepcopy(myfile["class_distr"])
        numOsamples = numOsamples - numOsamples_sofar

        if length(og.counts) != length(result_sofar)
            error("Loaded orbit-space fast-ion distribution from progress-file is non-compatible with specified orbit grid. Please correct orbit grid, provide valid progress-file or remove progress-file from work directory.")
        end
        close(myfile)
    else
        numOsamples_sofar = 0
        result_sofar = zeros(size(og.counts))
        class_distr_sofar = zeros(9)
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
        return ps2os_performance(M, wall, frdvols_cumsum_vector, subs, nfast, energy, pitch, R, z, og; numOsamples=numOsamples, numOsamples_sofar=numOsamples_sofar, 
                                 result_sofar=result_sofar, class_distr_sofar=class_distr_sofar, distributed=distributed, FI_species=FI_species, saveProgress=saveProgress, 
                                 progress_file_name=progress_file_name, verbose=verbose, kwargs...)
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

            result_p, class_distr_p = sample_helper(M, nbatch, fr, dvols, energy, pitch, R, z, og; wall=wall, FI_species=FI_species, kwargs...)
            result_sofar .+= result_p
            class_distr_sofar .+= class_distr_p
            numOsamples_sofar += nbatch
            if saveProgress
                rm("progress_$(progress_file_name).jld2", force=true) #clear the previous file
                myfile = jldopen("progress_$(progress_file_name).jld2",true,true,false,IOStream)
                write(myfile,"F_os",result_sofar)
                write(myfile,"class_distr",class_distr_sofar)
                write(myfile,"numOsamples",numOsamples_sofar)
                close(myfile)
            end
        end
        if verbose
            println("(Rest) Samples left: $(numOsamples)")
        end
        result_rest, class_distr_rest = sample_helper(M, numOsamples, fr, dvols, energy, pitch, R, z, og; wall=wall, FI_species=FI_species, kwargs...)
        numOsamples_rest = numOsamples

        if subdivide
            result = result_sofar + result_rest
            class_distr = class_distr_sofar + class_distr_rest
            numOsamples = numOsamples_sofar + numOsamples_rest
        else
            result = result_rest
            class_distr = class_distr_rest
            numOsamples = numOsamples_rest
        end
    else # ...if not parallel computing... I wish you good luck.
        dE_os_end = abs((og.energy)[end]-(og.energy)[end-1])
        dE_os_1 = abs((og.energy)[2]-(og.energy)[1])
        
        for i=numOsamples_sofar+1:numOsamples
            if verbose
                println("Sample number: $(i)")
            end
            E_sample, p_sample, R_sample, z_sample = OWCF_sample_f(fr, dvols, energy, pitch, R, z, n=1) # Returns 4 1-element arrays

            # CHECK IF IT'S A GOOD SAMPLE
            good_sample = checkIfGoodSample(E_sample[1], p_sample[1], R_sample[1], z_sample[1], energy, pitch, R, z)

            if good_sample # getGCP function is from OWCF/misc/species_func.jl
                o = get_orbit(M,getGCP(FI_species; E=E_sample[1], p=p_sample[1], R=R_sample[1], z=z_sample[1]); store_path=false, wall=wall, kwargs...)
                if (o.coordinate.energy <= (og.energy[end]+dE_os_end/2) && o.coordinate.energy >= (og.energy[1]-dE_os_1/2)) # Make sure it's within the energy bounds (+one half grid cell)
                    F_os_i = bin_orbits(og,Vector([o.coordinate]),weights=Vector([1.0]))
                    class_distr_i = zeros(9)
                    class_distr_i[class2int(M,o; plot=true, distinguishLost=true, distinguishIncomplete=true)] = 1
                else
                    F_os_i = zeros(length(og.counts))
                    class_distr_i = zeros(9)
                end
            else
                F_os_i = zeros(length(og.counts))
                class_distr_i = zeros(9)
            end

            result_sofar .+= F_os_i
            class_distr_sofar .+= class_distr_i

            if (i%nbatch)==0 && saveProgress # Every nbatch sample, save
                rm("progress_$(progress_file_name).jld2", force=true) #clear the previous file
                myfile = jldopen("progress_$(progress_file_name).jld2",true,true,false,IOStream)
                write(myfile,"F_os",result_sofar)
                write(myfile,"class_distr",class_distr_sofar)
                write(myfile,"numOsamples",i)
                close(myfile)
            end
        end
        result = result_sofar
        class_distr = class_distr_sofar
    end
    if verbose
        println("Number of good samples/All samples: $(sum(result)/numOsamples)")
    end

    rm("progress_$(progress_file_name).jld2", force=true) # Finally, remove the file that is no longer needed

    return result, class_distr, nfast
end

"""
    ps2os_performance(M, wall, fr, dvols, nfast, energy, pitch, R, z, og)
    ps2os_performance(-||-; numOsamples, numOsamples_sofar=0, result_sofar=zeros(size(og.counts))), distributed=true, FI_species="D", 
                            progress_file_name="ps2os", saveProgess=true, nbatch=1_000_000, verbose=false, kwargs...)

The performance version of part of ps2os(). This function will likely completely replace ps2os() in the near future. It computes necessary quantities
once instead of for every sample (as ps2os() does). Such as the element-wise product of fr and dvols, and its cumulative sum.
The input arguments are:
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
The keyword arguments are:
    numOsamples - The total number of Monte-Carlo samples (includes numOsamples_sofar)
    numOsamples_sofar - The number of Monte-Carlo samples sampled so far
    result_sofar - The (E,pm,Rm) fast-ion distribution in 1D compressed vector format. So far, having completed numOsamples_sofar number of samples
    distributed - If true, multi-core processing will be used
    FI_species - The fast-ion species. Please see OWCF/misc/species_func.jl
    progress_file_name - The name of the progress file will be "progress_[progress_file_name].jld2", if saveProgress is set to true
    saveProgress - If true, Monte-Carlo sampling process will be saved every nbatch samples
    nbatch - Compute the Monte-Carlo samples in batches, to optimize computation efficiency and enable subsequent progress saving
    verbose - If true, the function will talk a lot
    visualizeProgress - If true, a progress bar (or equivalent in the single-CPU core case) will be displayed during computations
"""
function ps2os_performance(M::AbstractEquilibrium, wall::Boundary, frdvols_cumsum_vector::AbstractVector, subs::CartesianIndices{4,NTuple{4,Base.OneTo{Int64}}},
                            nfast::Real, energy::AbstractVector, pitch::AbstractVector, R::AbstractVector, z::AbstractVector, og::OrbitGrid;
                            numOsamples::Int64, numOsamples_sofar::Int64=0, result_sofar=zeros(size(og.counts)), class_distr_sofar=zeros(9), distributed::Bool=true,
                            FI_species="D", progress_file_name="ps2os", saveProgress::Bool=true, nbatch::Int64 = 1_000_000, verbose::Bool=false, visualizeProgress::Bool=false,
                            kwargs...)
    verbose && println("Pre-computing difference vectors... ")
    dE_vector = vcat(abs.(diff(energy)),abs(energy[end]-energy[end-1]))
    dp_vector = vcat(abs.(diff(pitch)),abs(pitch[end]-pitch[end-1]))
    dR_vector = vcat(abs.(diff(R)),abs(R[end]-R[end-1]))
    dz_vector = vcat(abs.(diff(z)),abs(z[end]-z[end-1]))

    verbose && println("Starting Monte-Carlo computations... ")
    if distributed
        subdivide = false
        while numOsamples > nbatch
            subdivide = true
            verbose && println("Samples left: $(numOsamples)")
            numOsamples = numOsamples - nbatch
            result_p, class_distr_p = performance_helper(M, nbatch, frdvols_cumsum_vector, subs, dE_vector, dp_vector, dR_vector, dz_vector, energy, pitch, R, z, og; wall=wall, FI_species=FI_species, kwargs...)
            result_sofar .+= result_p
            class_distr_sofar .+= class_distr_p
            numOsamples_sofar += nbatch
            if saveProgress
                rm("progress_$(progress_file_name).jld2", force=true) #clear the previous file
                myfile = jldopen("progress_$(progress_file_name).jld2",true,true,false,IOStream)
                write(myfile,"F_os",result_sofar)
                write(myfile,"class_distr", class_distr_sofar)
                write(myfile,"numOsamples",numOsamples_sofar)
                close(myfile)
            end
        end
        verbose && println("(Rest) Samples left: $(numOsamples)")
        result_rest, class_distr_rest = performance_helper(M, numOsamples, frdvols_cumsum_vector, subs, dE_vector, dp_vector, dR_vector, dz_vector, energy, pitch, R, z, og; wall=wall, FI_species=FI_species, kwargs...)
        numOsamples_rest = numOsamples

        if subdivide
            result = result_sofar + result_rest
            class_distr = class_distr_sofar + class_distr_rest
            numOsamples = numOsamples_sofar + numOsamples_rest
        else
            result = result_rest
            class_distr = class_distr_rest
            numOsamples = numOsamples_rest
        end
    else # ...if not parallel computing... I wish you good luck.
        dE_os_end = abs((og.energy)[end]-(og.energy)[end-1])
        dE_os_1 = abs((og.energy)[2]-(og.energy)[1])

        for i=numOsamples_sofar+1:numOsamples
            visualizeProgress && print("")
            visualizeProgress && print("Sample number: $(i)")
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
            visualizeProgress && print("   ($(round(E_sample,sigdigits=5)),$(round(p_sample,sigdigits=3)),$(round(R_sample,sigdigits=4)),$(round(z_sample,sigdigits=4)))")
            good_sample = checkIfGoodSample(E_sample, p_sample, R_sample, z_sample, energy, pitch, R, z)
            visualizeProgress && print("   $(good_sample)")
            if good_sample
                o = get_orbit(M,getGCP(FI_species; E=E_sample,p=p_sample,R=R_sample,z=z_sample); store_path=false, wall=wall, kwargs...)
                visualizeProgress && print("   $(o.class)")
                visualizeProgress && print("   $(o.coordinate.energy)")
                if (o.coordinate.energy <= (maximum(og.energy)+dE_os_end/2) && o.coordinate.energy >= (minimum(og.energy)-dE_os_1/2)) # Make sure it's within the energy bounds (+one half grid cell)
                    F_os_i = bin_orbits(og,Vector([o.coordinate]),weights=Vector([1.0]))
                    class_distr_i = zeros(9)
                    class_distr_i[class2int(M,o; plot=true, distinguishLost=true, distinguishIncomplete=true)] = 1
                else
                    F_os_i = zeros(length(og.counts))
                    class_distr_i = zeros(9) # No orbit class
                end
            else
                visualizeProgress && print("   $(round(energy[1],sigdigits=5))<E<$(round(energy[end],sigdigits=5))")
                visualizeProgress && print("   $(round(pitch[1],sigdigits=3))<p<$(round(pitch[end],sigdigits=3))")
                visualizeProgress && print("   $(round(R[1],sigdigits=4))<R<$(round(R[end],sigdigits=4))")
                visualizeProgress && print("   $(round(z[1],sigdigits=4))<z<$(round(z[end],sigdigits=4))")
                F_os_i = zeros(length(og.counts))
                class_distr_i = zeros(9) # No orbit class
            end

            result_sofar .+= F_os_i
            class_distr_sofar .+= class_distr_i

            if (i%nbatch)==0 && saveProgress # Every nbatch sample, save
                rm("progress_$(progress_file_name).jld2", force=true) #clear the previous file
                myfile = jldopen("progress_$(progress_file_name).jld2",true,true,false,IOStream)
                write(myfile,"F_os",result_sofar)
                write(myfile,"class_distr", class_distr_sofar)
                write(myfile,"numOsamples",i)
                close(myfile)
            end
            visualizeProgress && println("")
        end
        result = result_sofar
        class_distr = class_distr_sofar
    end

    verbose && println("Number of good samples/All samples: $(sum(result)/numOsamples)")
    rm("progress_$(progress_file_name).jld2", force=true) # As in ps2os(), remove the progress file that is no longer needed

    return result, class_distr, nfast
end

"""
    sample_helper(M, numOsamples, fr, dvols, energy, pitch, R, z, og)

Help the function ps2os() with acquiring orbit samples when parallel computations are desired.
This is to enable the sampling process to be saved regularly when calculating a large number of samples.
If the sampling process is not saved, then progress will be lost when the super-user of the HPC terminates
the sampling process early, due to misinterpretation of Julia's way of distributed computing.
"""
function sample_helper(M::AbstractEquilibrium, numOsamples::Int64, fr::AbstractArray, dvols::AbstractArray, energy::AbstractVector, pitch::AbstractVector, R::AbstractVector, z::AbstractVector, og::OrbitGrid; wall::Union{Nothing,Boundary}, FI_species="D", visualizeProgress::Bool=false, kwargs...)

    dE_os_end = abs((og.energy)[end]-(og.energy)[end-1])
    dE_os_1 = abs((og.energy)[2]-(og.energy)[1])

    if numOsamples>0.0 # If there are actually a non-zero number of samples left to sample...
        if visualizeProgress
            prog = ProgressMeter.Progress(numOsamples) # Define the progress bar
            channel = RemoteChannel(()->Channel{Bool}(numOsamples),1) # Define the channel from which the progress bar draws data
            result = fetch(@sync begin # Start the distributed computational process, fetch result when done
                @async while take!(channel) # An asynchronous process, with no need for sync, since it simply displays the progress bar
                    ProgressMeter.next!(prog)
                end

                @async begin # No internal syncronization needed here either, only external sync needed
                    result_async = @distributed (+) for i=1:numOsamples # Compute one result, and reduce (add) it to a resulting vector F_os

                        # Sample
                        E_sample, p_sample, R_sample, z_sample = OWCF_sample_f(fr, dvols, energy, pitch, R, z, n=1) # Returns 4 1-element arrays
                        # CHECK IF IT'S A GOOD SAMPLE
                        good_sample = checkIfGoodSample(E_sample[1], p_sample[1], R_sample[1], z_sample[1], energy, pitch, R, z)
                        if good_sample # getGCP function is from OWCF/misc/species_func.jl
                            o = get_orbit(M,getGCP(FI_species; E=E_sample[1], p=p_sample[1], R=R_sample[1], z=z_sample[1]); store_path=false, wall=wall, kwargs...) # Calculate the orbit
                            if (o.coordinate.energy <= (maximum(og.energy)+dE_os_end/2) && o.coordinate.energy >= (minimum(og.energy)-dE_os_1/2)) # Make sure it's within the energy bounds (+one half grid cell)
                                F_os_i = bin_orbits(og,Vector([o.coordinate]),weights=Vector([1.0])) # Bin to the orbit grid
                                class_distr_i = zeros(9)
                                class_distr_i[class2int(M,o; plot=true, distinguishLost=true, distinguishIncomplete=true)] = 1
                            else
                                F_os_i = zeros(length(og.counts)) # Otherwise, zero
                                class_distr_i = zeros(9) # No orbit class
                            end
                        else
                            F_os_i = zeros(length(og.counts)) # Otherwise, zero
                            class_distr_i = zeros(9) # No orbit class
                        end
                        put!(channel,true) # Update the progress bar
                        result_i = [F_os_i, class_distr_i]
                        result_i # Declare result_i as a result to add to result_async
                    end
                    put!(channel,false) # Update progress bar
                    result_async # Delcare result_async_os as done/result, so it can be fetched
                end
            end)
        else
            result = @distributed (+) for i=1:numOsamples
                E_sample, p_sample, R_sample, z_sample = OWCF_sample_f(fr, dvols, energy, pitch, R, z, n=1) # Returns 4 1-element arrays
                good_sample = checkIfGoodSample(E_sample[1], p_sample[1], R_sample[1], z_sample[1], energy, pitch, R, z)
                if good_sample # getGCP function is from OWCF/misc/species_func.jl
                    o = get_orbit(M,getGCP(FI_species; E=E_sample[1],p=p_sample[1],R=R_sample[1],z=z_sample[1]); store_path=false, wall=wall, kwargs...) # Calculate the orbit
                    if (o.coordinate.energy <= (maximum(og.energy)+dE_os_end/2) && o.coordinate.energy >= (minimum(og.energy)-dE_os_1/2)) # Make sure it's within the energy bounds (+one half grid cell)
                        F_os_i = bin_orbits(og,Vector([o.coordinate]),weights=Vector([1.0])) # Bin to the orbit grid
                        class_distr_i = zeros(9)
                        class_distr_i[class2int(M,o; plot=true, distinguishLost=true, distinguishIncomplete=true)] = 1
                    else
                        F_os_i = zeros(length(og.counts)) # Otherwise, zero
                        class_distr_i = zeros(9) # No orbit class
                    end
                else
                    F_os_i = zeros(length(og.counts)) # Otherwise, zero
                    class_distr_i = zeros(9) # No orbit class
                end
                result_i = [F_os_i, class_distr_i]
                result_i # Declare result_i as a result to add to result_async
            end
        end

        class_distr = result[2] # Extract the (non-interpolated) distribution of orbit classes (including lost, invalid and incomplete orbits)
        result = result[1] # A vector representing the (interpolated) fast-ion distribution, binned onto a (finite) orbit grid
        return result, class_distr
    else # ...otherwise just return a zero vector with the length of the number of valid orbits
        return zeros(length(og.counts)), zeros(9)
    end
end

"""
    performance_helper(M, nbatch, frdvols_cumsum_vector, subs, dE_vector, dp_vector, dR_vector, dz_vector, energy, pitch, R, z, og; wall=wall, FI_species="D", kwargs...)

Help the function ps2os_performance() with acquiring orbit samples when parallel computations are desired.
This is to enable the sampling process to be saved regularly when calculating a large number of samples.
If the sampling process is not saved, then progress will be lost when the super-user of the HPC terminates
the sampling process early, due to misinterpretation of Julia's way of distributed computing.
"""
function performance_helper(M::AbstractEquilibrium, numOsamples::Int64, frdvols_cumsum_vector::AbstractVector, subs::CartesianIndices{4,NTuple{4,Base.OneTo{Int64}}}, 
                            dE_vector::AbstractVector, dp_vector::AbstractVector, dR_vector::AbstractVector, dz_vector::AbstractVector, energy::AbstractVector, 
                            pitch::AbstractVector, R::AbstractVector, z::AbstractVector, og::OrbitGrid; 
                            wall::Union{Nothing,Boundary{Float64}}, FI_species="D", visualizeProgress::Bool=false, kwargs...)

    dE_os_end = abs((og.energy)[end]-(og.energy)[end-1])
    dE_os_1 = abs((og.energy)[2]-(og.energy)[1])

    if numOsamples>0 # If there are actually a non-zero number of samples left to sample...
        if visualizeProgress
            prog = ProgressMeter.Progress(numOsamples) # Define the progress bar
            channel = RemoteChannel(()->Channel{Bool}(numOsamples),1) # Define the channel from which the progress bar draws data
            result = fetch(@sync begin # Start the distributed computational process, fetch result when done
                @async while take!(channel) # An asynchronous process, with no need for sync, since it simply displays the progress bar
                    ProgressMeter.next!(prog)
                end

                @async begin # No internal syncronization needed here either, only external sync needed
                    result_async = @distributed (+) for i=1:numOsamples # Compute one result, and reduce (add) it to a resulting vector F_os

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
                        if good_sample # getGCP function is from OWCF/misc/species_func.jl
                            o = get_orbit(M,getGCP(FI_species; E=E_sample,p=p_sample,R=R_sample,z=z_sample); store_path=false, wall=wall, kwargs...) # Calculate the orbit
                            if (o.coordinate.energy <= (maximum(og.energy)+dE_os_end/2) && o.coordinate.energy >= (minimum(og.energy)-dE_os_1/2)) # Make sure it's within the energy bounds (+one half grid cell)
                                F_os_i = bin_orbits(og,Vector([o.coordinate]),weights=Vector([1.0])) # Bin to the orbit grid
                                class_distr_i = zeros(9)
                                class_distr_i[class2int(M,o; plot=true, distinguishLost=true, distinguishIncomplete=true)] = 1
                            else
                                F_os_i = zeros(length(og.counts)) # Otherwise, zero
                                class_distr_i = zeros(9) # No orbit class
                            end
                        else
                            F_os_i = zeros(length(og.counts)) # Otherwise, zero
                            class_distr_i = zeros(9) # No orbit class
                        end
                        put!(channel,true) # Update the progress bar
                        result_i = [F_os_i, class_distr_i]
                        result_i # Declare result_i as a result to add to result_async
                    end
                    put!(channel,false) # Update progress bar
                    result_async # Delcare result as done/result, so it can be fetched
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
                if good_sample # getGCP function is from OWCF/misc/species_func.jl
                    o = get_orbit(M,getGCP(FI_species; E=E_sample, p=p_sample, R=R_sample, z=z_sample); store_path=false, wall=wall, kwargs...) # Calculate the orbit
                    if (o.coordinate.energy <= (maximum(og.energy)+dE_os_end/2) && o.coordinate.energy >= (minimum(og.energy)-dE_os_1/2)) # Make sure it's within the energy bounds (+/- one half grid cell)
                        F_os_i = bin_orbits(og,Vector([o.coordinate]),weights=Vector([1.0])) # Bin to the orbit grid
                        class_distr_i = zeros(9)
                        class_distr_i[class2int(M,o; plot=true, distinguishLost=true, distinguishIncomplete=true)] = 1
                    else
                        F_os_i = zeros(length(og.counts)) # Otherwise, zero
                        class_distr_i = zeros(9) # No orbit class
                    end
                else
                    F_os_i = zeros(length(og.counts)) # Otherwise, zero
                    class_distr_i = zeros(9) # No orbit class
                end
                result_i = [F_os_i, class_distr_i]
                result_i # Declare result_i as a result to add to result
            end
        end

        class_distr = result[2] # Extract the (non-interpolated) distribution of orbit classes (including lost, invalid and incomplete orbits)
        result = result[1] # A vector representing the (interpolated) fast-ion distribution, binned onto a (finite) orbit grid
        return result, class_distr
    else # ...otherwise just return a zero vector with the length of the number of valid orbits
        return zeros(length(og.counts)), zeros(9) # USE class2int(M, o; plot=true, distinguishLost=true, distinguishIncomplete=true) to bin orbit classes to 9-element Int vector. Distributed can reduce vector of vectors of different lengths
    end
end

"""

checkIfGoodSample(E_sample, p_sample, R_sample, z_sample, E_array, p_array, R_array, z_array)

The function checks if a sample is within bounds. Returns true if that is the case. Otherwise false.
"""
function checkIfGoodSample(E_sample, p_sample, R_sample, z_sample, E_array, p_array, R_array, z_array)
    
    if (E_sample <= minimum(E_array) || E_sample >= maximum(E_array)  || 
        p_sample <= minimum(p_array) || p_sample >= maximum(p_array)  || 
        R_sample <= minimum(R_array) || R_sample >= maximum(R_array)  || 
        z_sample <= minimum(z_array) || z_sample >= maximum(z_array))
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

    dorbs = get_orbel_volume(og)
    if verbose
        println("Unmapping the orbits... ")
    end
    F_os = zeros(length(og.counts))
    for i in eachindex(F_os_3D)
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
    FI_species="D", wall, verbose=false, kwargs...)

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
        lorbs = reshape([get_orbit(M, getGCP(FI_species; E=energy[k], p=pitch[l], R=rr, z=zz); wall=wall, kwargs...) for k=1:nenergy,l=1:npitch],nenergy*npitch)

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
        end
        last_ind = I
    end

    return EPRZDensity(f_eprz,energy,pitch,r,z) # Return as EPRZDensity object
end

"""
    getOSTopoMap(M, E, pmRm_inds, pm_array, Rm_array)
    getOSTopoMap(-||-; FI_species="D", wall=nothing, distinguishLost=false, distinguishIncomplete=false, extra_kw_args=Dict(:toa => true))

Compute the orbit-sapce topological map for the energy slice specified by E (keV). Use pmRm_inds to iterate through the pm and Rm coordinates. 
"""
function getOSTopoMap(M::AbstractEquilibrium, E::Float64, pmRm_inds_array::Vector{Tuple{CartesianIndex{1}, CartesianIndex{1}}}, pm_array::AbstractVector, 
                      Rm_array::AbstractVector; FI_species::String="D", wall::Union{Nothing,Boundary{Float64}}=nothing, distinguishLost::Bool=false, 
                      distinguishIncomplete::Bool=false, show_progress::Bool=false, extra_kw_args=Dict(:toa => true))
    topoMap = @showprogress enabled=show_progress dt=1 desc="Computing topoMap for E=$(E) keV... " @distributed (+) for i in eachindex(pmRm_inds_array)
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
    getOSTopoMap(-||- ; kwargs... )

Compute the orbit-space topological map, given a grid spanned by the inputs E_array, pm_array and Rm_array.
Output the topological map as a 3D Array of Float64. Float64 type is advantageous for compatibility reasons (with other OWCF functions/scripts)
"""
function getOSTopoMap(M::AbstractEquilibrium, E_array::AbstractVector, pm_array::AbstractVector, Rm_array::AbstractVector; kwargs... )
    pm_inds_rep = repeat(CartesianIndices(pm_array),inner=length(Rm_array)) # To get all points
    Rm_inds_rep = repeat(CartesianIndices(Rm_array),outer=length(pm_array)) # To get all points
    topoMap = zeros(length(E_array),length(pm_array),length(Rm_array)) # The total 3D topological map
    pmRm_inds_array = collect(zip(pm_inds_rep,Rm_inds_rep)) # All points in an energy slice, in one long vector
    for iE in eachindex(E_array)
        topoMap[iE,:,:] = getOSTopoMap(M, E_array[iE], pmRm_inds_array, pm_array, Rm_array; kwargs... )
    end
    return Float64.(topoMap)
end

"""
    getCOMTopoMap(M, E, Λ_array, Pϕ_n_array)
    getCOMTopoMap(-||-; sigma=1, kwargs...)

Compute the constants-of-motion (COM) space topological map for a specific energy slice given by the energy 'E' (keV).
Output the topological map as a 2D Array of Int64. The input arguments are:
- M - The abstract equilibrium Struct (Equilibrium.jl) containing all information about the magnetic equilibrium - AbstractEquilibrium
- E - The energy (keV) of the topological map - Real
- Λ_array - The normalized magnetic moment (Λ=μ*B0/E with B0=B(mag. axis)) grid points of the topological map - AbstractVector
- Pϕ_n_array - The normalized toroidal canonical angular momentum (Pϕ_n=Pϕ/(q*|Ψ_w|) where q is the particle charge in Coulomb and Ψ_w is the poloidal
                 flux at the last closed flux surface (LCFS)) grid points of the topological map. Please note! If Ψ(LCFS)==0 for some reason (e.g. by some
                 convention), then Ψ_w=Ψ(mag. axis) is assumed - AbstractVector
The keyword arguments are:
- sigma - Binary index to identify co-passing (1) or counter-passing (-1) orbit, when (E,Λ,Pϕ_n) is a degenerate triplet - Int64
- R_array - An array with major radius (R) points to be included as R grid points for the mu-contours - AbstractVector
- z_array - An array with vertical coordinate (z) points to be included as z grid points for the mu-contours - AbstractVector
- nR - If R_array is unspecified, nR major radius grid points between R_LCFS_HFS and R_LCFS_LFS will be used. Defaults to 100 - Int64
- nz - If z_array is unspecified, nz major radius grid points between z_LCFS_HFS and z_LCFS_LFS will be used. Defaults to 120 - Int64
- B_abs_Rz - The norm of the magnetic field at all (R,z) points, i.e. |B(R,z)|. Computed internally, if not specified - AbstractMatrix
- psi_Rz - The poloidal flux at all (R,z) points, i.e. ψ(R,z). Computed internally, if not specified - AbstractMatrix
- RBT_Rz - The poloidal current at all (R,z) points, i.e. R(R,z) .*Bϕ(R,z). Computed internally, if not specified - AbstractMatrix
- FI_species - The particle species of the ion. Please see OWCF/misc/species_func.jl for a list of available particle species - String
- show_progress - If set to true, a progress meter will be shown to indicate the computational progress - Bool
- verbose - If set to true, the function will talk a lot! - Bool
- vverbose - If set to true, the function will talk even more, and plot too! - Bool
"""
function getCOMTopoMap(M::AbstractEquilibrium, E::Real, Λ_array::AbstractVector, Pϕ_n_array::AbstractVector; 
                       sigma::Int64=1, R_array::AbstractVector=nothing, z_array::AbstractVector=nothing, 
                       nR::Int64=100, nz::Int64=120, B_abs_Rz::Union{AbstractMatrix,Nothing}=nothing,
                       psi_Rz::Union{AbstractMatrix,Nothing}=nothing,RBT_Rz::Union{AbstractMatrix,Nothing}=nothing,
                       FI_species::String="D", show_progress::Bool=false, verbose::Bool=false, vverbose::Bool=false)
    if !(sigma==1 || sigma==-1)
        error("Keyword 'sigma' must be specified as either +1 or -1. Please correct and re-try.")
    end
    orb_select_func = sigma==1 ? argmax : argmin # If co-current (sigma==1) orbits are of interest, use argmax function. Otherwise (if sigma==-1), use argmin
    E_joule = E*(GuidingCenterOrbits.e0)*1000 # from keV to Joule
    m_FI = getSpeciesMass(FI_species) # kg
    q_FI = getSpeciesCharge(FI_species) # Coloumb
    psi_axis, psi_bdry = psi_limits(M)
    if psi_bdry==0
        @warn "The magnetic flux at the last closed flux surface (LCFS) is found to be 0 for the 'M' input to getCOMTopoMap(). Pϕ_n=Pϕ/(q*|Ψ_w|) where Ψ_w=Ψ(mag. axis) is assumed instead of Ψ_w=Ψ(LCFS)."
        Ψ_w_norm = abs(psi_axis)
    else
        Ψ_w_norm = abs(psi_bdry)
    end
    B0 = norm(Equilibrium.Bfield(M,magnetic_axis(M)...))
    mu_array = (E_joule/B0) .*Λ_array # Transform to μ, for easy contour computation
    Pphi_array = (q_FI*Ψ_w_norm) .*Pϕ_n_array # Transform to Pϕ, for easy contour computation

    if vverbose || isnothing(R_array) || isnothing(z_array)
        if :r in fieldnames(typeof(M)) # If the magnetic equilibrium is an .eqdsk-equilibrium...
            LCFS = Equilibrium.boundary(M,psi_bdry) # Default 0.01 m (R,z) resolution to identify LCFS
        else # ...else, it must be a Solov'ev equilibrium (currently the only two types supported by the package Equilibrium.jl)
            LCFS = Equilibrium.boundary(M;n=500) # Default 500 (R,z) points to represent LCFS
        end
    end
    if isnothing(R_array) || isnothing(z_array)
        verbose && println("getCOMTopoMap(): R_array/z_array keyword argument not specified. Computing internally... ")
        if isnothing(R_array)
            R_array = collect(range(minimum(LCFS.r),stop=maximum(LCFS.r),length=nR))
        end
        if isnothing(z_array)
            z_array = collect(range(minimum(LCFS.z),stop=maximum(LCFS.z),length=nz))
        end
    end
    if isnothing(B_abs_Rz)
        verbose && println("getCOMTopoMap(): B_abs_Rz keyword argument not specified. Computing internally... ")
        B_abs_Rz = [norm(Equilibrium.Bfield(M,R,z)) for R in R_array, z in z_array]
    end
    if isnothing(psi_Rz)
        verbose && println("getCOMTopoMap(): psi_Rz keyword argument not specified. Computing internally... ")
        psi_Rz = [M(R,z) for R in R_array, z in z_array]
    end
    if isnothing(RBT_Rz)
        verbose && println("getCOMTopoMap(): RBT_Rz keyword argument not specified. Computing internally... ")
        RBT_Rz = [poloidal_current(M,M(R,z)) for R in R_array, z in z_array] # RBT is short for R*Bϕ
    end
    verbose && println("getCOMTopoMap(): Defining magnetic moment (mu) as a function of toroidal canonical angular momentum (Pphi) input... ")
    function _mu_func(Pphi::Real)
        res = E_joule ./B_abs_Rz .- (B_abs_Rz ./(2*m_FI)) .* ((Pphi .- q_FI .*psi_Rz) ./RBT_Rz).^2
        return map(x-> (x > 0.0) ? x : 0.0, res)
    end

    topoMap = @showprogress enabled=show_progress dt=1 desc="Computing topoMap for E=$(E) keV... " @distributed (+) for iPphi in eachindex(Pphi_array)
        Pphi = Pphi_array[iPphi]
        topoMap_i = zeros(Int64,(length(mu_array),length(Pphi_array))) # Initially assume all (E,mu,Pphi) triplets are invalid (9)
        mu_Rz = _mu_func(Pphi) # Mu as a function of (R,z), for a specific Pphi
        mu_max = maximum(mu_Rz)
        valid_mu_indices = findall(x-> x<=mu_max && x>0,mu_array)
        cls = Contour.contours(R_array,z_array,mu_Rz,mu_array[valid_mu_indices])
        cls = Contour.levels(cls) # Turn into iterable
        for j in eachindex(cls) # For each mu-contour index
            lns = Contour.lines(cls[j]) # Get the lines of that contour level (might be several lines)
            closed_lines = []
            Rm_of_closed_lines = []
            for ln in lns # For line in lines
                Rl, zl = Contour.coordinates(ln) # Get the (R,z) points of the contour line
                if Rl[1]==Rl[end] && zl[1]==zl[end] # Is it a closed mu-contour? If so, Contour.jl puts the first point equal to the last
                    Rm = Rl[argmax(Rl)] # The maximum major radius (Rm) coordinate of the closed mu-contour
                    push!(closed_lines,ln)
                    push!(Rm_of_closed_lines,Rm)
                end
            end
            if !isempty(closed_lines) # If there are valid orbits to consider
                mu = mu_array[valid_mu_indices[j]] # The mu value for this particular orbit (closed line)
                oi = orb_select_func(Rm_of_closed_lines) # Select the relevant one (co- or counter-current)
                Ro, zo = Contour.coordinates(closed_lines[oi]) # Get the (R,z) points of the orbit
                po = [sqrt(getSpeciesMass(FI_species)/(2*E_joule))*(Pphi-getSpeciesCharge(FI_species)*M(Ro[k],zo[k]))*norm(vec(Equilibrium.Bfield(M,Ro[k],zo[k])))/(getSpeciesMass(FI_species)*Ro[k]*vec(Equilibrium.Bfield(M,Ro[k],zo[k]))[2]) for k in eachindex(Ro)] # Compute the (classical) pitch (v_||/v) from Pphi and E. v_||/v = sqrt(m/(2*E))*(Pphi-q*ψ)*B/(m*R*B_phi)
                po = clamp.(po,-1.0,1.0) # Clamp the pitch between -1.0 and 1.0, to adress numerical inaccuracies
                o_class = GuidingCenterOrbits.classify(Ro,zo,po,magnetic_axis(M)) # Deduce the orbit class from the (R,z), pitch and magnetic axis points
                if vverbose 
                    # To use this, include 'using Plots' at the top of dependencies.jl, and un-comment the lines below
                    #plt_crs = Plots.plot(Ro,zo,title="$(o_class)",aspect_ratio=:equal,label="")
                    #plt_crs = Plots.plot!(plt_crs,LCFS.r,LCFS.z,label="LCFS")
                    #plt_p = Plots.plot(po,title="Pitch evolution (v_||/v)")
                    #Plots.display(Plots.plot(plt_crs,plt_p,layout=(1,2)))
                end
                o_path = OrbitPath(false,true,E*ones(length(Ro)),po,Ro,zo,zero(zo),zero(zo)) # Create an OrbitPath object (struct). Don't assume vacuum (false), include drift effects (true), and we don't know anything about phi and dt (zero(zo))
                o = Orbit(HamiltonianCoordinate(E, mu, Pphi; amu=getSpeciesAmu(FI_species), q=getSpeciesEcu(FI_species)), o_class, zero(E), zero(E), o_path, false) # Create an Orbit struct. We don't know tau_p and tau_t (zero(E)) and we don't know if guiding-center equation of motion are valid (false)
                topoMap_i[valid_mu_indices[j],iPphi] = class2int(M, o; sigma=sigma) # Give the orbit an integer identification number, according to OWCF/calcTopoMap.jl rules
            end
        end
        topoMap_i # Declare for @distribution reduction (+)
    end
    topoMap[topoMap .== 0.0] .= 9 # Set all invalid orbits to integer 9
    return topoMap
end

"""
    getCOMTopoMap(M, E_array, Λ_array, Pϕ_n_array)
    getCOMTopoMap(-||- ; sigma=1, R_array=nothing, z_array=nothing, nR=100, nz=120, B_abs_Rz=nothing, psi_Rz=nothing, RBT_Rz=nothing, verbose=false, kwargs...)

Compute the constants-of-motion space topological map, given a grid spanned by the inputs E_array, Λ_array and Pϕ_n_array.
Output the topological map as a 3D Array of Int64. The input arguments are:
- M - The abstract equilibrium Struct (Equilibrium.jl) containing all information on the magnetic equilibrium - AbstractEquilibrium
- E_array - The energy (keV) grid points of the topological map - AbstractVector
- Λ_array - The normalized magnetic moment (μ) grid points of the topological map. Λ=μ*B0/E where B0=B(mag. axis) - AbstractVector
- Pϕ_n_array - The normalized toroidal canonical angular momentum (Pϕ) grid points of the topological map. Pϕ_n=Pϕ/(q*|Ψ_w|) where
               q is the particle charge in Coulomb and Ψ_w=Ψ(LCFS). Ψ is the poloidal magnetic flux and LCFS means 'last closed flux
               surface'. If Ψ(LCFS)==0 for some reason (e.g. because of some convention), Ψ_w=Ψ(mag. axis) is assumed - AbstractVector
The keyword arguments are:
- sigma - Binary index to identify co-current (1) or counter-current (-1) orbit, when (E,Λ,Pϕ) is a degenerate triplet - Int64
- R_array - An array with major radius (R) points to be included as R grid points for the mu-contours - AbstractVector
- z_array - An array with vertical coordinate (z) points to be included as z grid points for the mu-contours - AbstractVector
- nR - If R_array is unspecified, nR major radius grid points between R_LCFS_HFS and R_LCFS_LFS will be used. Defaults to 100 - Int64
- nz - If z_array is unspecified, nz major radius grid points between z_LCFS_HFS and z_LCFS_LFS will be used. Defaults to 120 - Int64
- B_abs_Rz - The norm of the magnetic field at all (R,z) points, i.e. ||B(R,z)|| - AbstractMatrix
- psi_Rz - The poloidal flux at all (R,z) points, i.e. ψ(R,z) - AbstractMatrix
- RBT_Rz - The poloidal current at all (R,z) points, i.e. R(R,z) .*Bϕ(R,z) - AbstractMatrix
- verbose - If set to true, the function will talk a lot! - Bool
"""
function getCOMTopoMap(M::AbstractEquilibrium, E_array::AbstractVector, Λ_array::AbstractVector, 
                       Pϕ_n_array::AbstractVector; sigma::Int64=1, R_array::Union{AbstractVector,Nothing}=nothing, 
                       z_array::Union{AbstractVector,Nothing}=nothing, nR::Int64=100, nz::Int64=120, 
                       B_abs_Rz::Union{AbstractMatrix,Nothing}=nothing, psi_Rz::Union{AbstractMatrix,Nothing}=nothing,
                       RBT_Rz::Union{AbstractMatrix,Nothing}=nothing, verbose::Bool=false, kwargs...)
    if !(sigma==1 || sigma==-1)
        error("Keyword 'sigma' must be specified as either +1 or -1. Please correct and re-try.")
    end

    if isnothing(R_array) || isnothing(z_array)
        verbose && println("getCOMTopoMap(): R_array/z_array keyword argument not specified. Computing internally... ")
        if :r in fieldnames(typeof(M)) # If the magnetic equilibrium is an .eqdsk-equilibrium...
            LCFS = Equilibrium.boundary(M,psi_limits(M)[2]) # Default 0.01 m (R,z) resolution to identify LCFS
        else # ...else, it must be a Solov'ev equilibrium (currently the only two types supported by the package Equilibrium.jl)
            LCFS = Equilibrium.boundary(M;n=500) # Default 500 (R,z) points to represent LCFS
        end
        if isnothing(R_array)
            R_array = collect(range(minimum(LCFS.r),stop=maximum(LCFS.r),length=nR))
        end
        if isnothing(z_array)
            z_array = collect(range(minimum(LCFS.z),stop=maximum(LCFS.z),length=nz))
        end
    end
    if isnothing(B_abs_Rz)
        verbose && println("getCOMTopoMap(): B_abs_Rz keyword argument not specified. Computing internally... ")
        B_abs_Rz = [norm(Equilibrium.Bfield(M,R,z)) for R in R_array, z in z_array]
    end
    if isnothing(psi_Rz)
        verbose && println("getCOMTopoMap(): psi_Rz keyword argument not specified. Computing internally... ")
        psi_Rz = [M(R,z) for R in R_array, z in z_array]
    end
    if isnothing(RBT_Rz)
        verbose && println("getCOMTopoMap(): RBT_Rz keyword argument not specified. Computing internally... ")
        RBT_Rz = [poloidal_current(M,M(R,z)) for R in R_array, z in z_array] # RBT is short for R*Bϕ
    end

    topoMap = 9 .*ones(Int64, (length(E_array),length(Λ_array),length(Pϕ_n_array))) # Initially assume all (E,Λ,Pϕ) triplets are invalid (9)
    for (iE,E) in enumerate(E_array)
        topoMap[iE,:,:] = getCOMTopoMap(M, E, Λ_array, Pϕ_n_array; sigma=sigma, R_array=R_array, z_array=z_array, B_abs_Rz=B_abs_Rz, psi_Rz=psi_Rz, RBT_Rz=RBT_Rz, verbose=verbose, kwargs...)
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
function pmRm_2_μPϕ(M::AbstractEquilibrium, good_coords_pmRm::Vector{CartesianIndex{2}}, data::Array{Float64,2}, E::Real, pm_array::AbstractVector, Rm_array::AbstractVector, 
                    FI_species::AbstractString; nμ::Int64=length(pm_array), nPϕ::Int64=length(Rm_array), isTopoMap::Bool=false, needJac::Bool=false, 
                    transform = x -> x, verbose::Bool=false, vverbose::Bool=false, debug::Bool=false)

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
            # Write your own debug code here
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
        verbose && println("Returning debug quantities (pmRm_2_μPϕ())... ")
        return μ_values_cogo, μ_values_ctgo, Pϕ_values_cogo, Pϕ_values_ctgo, μ_array, Pϕ_array
    end

    return map(x-> isnan(x) ? 0.0 : x, data_COM), E, μ_array, Pϕ_array # Make sure to remove all NaNs. Sometimes, NaNs arise. It's inevitable.
end


"""
    os2COM(M, good_coords, data, E_array, pm_array, Rm_array, FI_species)
    os2COM(M, good_coords, data, E_array, pm_array, Rm_array, FI_species; nl=length(pm_array), npp=length(Rm_array), verbose=false)

This function maps an orbit-space (E,pm,Rm) quantity into COM-space (E, Λ, Pϕ_n; σ). E is the energy (keV), pm is the pitch (v_||/v) at the 
maximum major radius 'Rm' position of the orbit. Λ=μ*B_0/E where μ is the magnetic moment and B_0 is the norm of the magnetic field vector (T) at the 
magnetic axis. Pϕ_n=Pϕ/(q*|Ψ_w|) where Pϕ is the toroidal canonical angular momentum, q is the particle charge (Coulomb) and Ψ_w is the value of the (poloidal)
magnetic flux at the last closed flux surface (unless it is 0 due to some convention, then the value of the poloidal magnetic flux at the magnetic axis should be used).
The os2COM() function assumes the quantity 'data' is given a 3D array and that the (E,pm,Rm) coordinates suitable for mapping are given as good_coords, 
a Vector{CartesianIndex{3}}. The output data will be 4D. The last dimension will correspond to the binary σ coordinate (+1 or -1) that keeps track of co- and 
counter-passing orbits with the same (E, Λ, Pϕ_n) coordinate.

The input arguments are:
M - The magnetic equilibrium object
good_coords - The (E,pm,Rm) coordinates suitable (valid) for mapping
data - The data (a function of (E,pm,Rm)) to be mapped to (E, Λ, Pϕ_n; σ)
E_array - The energy grid points. In keV
pm_array - The pitch maximum grid points
Rm_array - The maximum major radius grid points. In meters
FI_species - The particle species (of the fast ion). Please see OWCF/misc/species_func.jl for more info

The keyword arguments are:
isTopoMap - If set to true, the function will assume that the input 'data' is a topological map. Defaults to false
Lambda_array - The Λ grid points to be mapped onto. Defaults to nothing
Pphi_n_array - The Pϕ_n grid points to be mapped onto. Please note! Pphi_n=Pphi/(q*|Ψ_w|) where Ψ_w=Ψ(LCFS). If Ψ(LCFS)==0, Ψ_w=Ψ(mag. axis) should be used instead. Defaults to nothing
nl - If Lambda_array is not specified, the number of Λ grid points to be mapped onto (max/min Λ will be automatically determined). Defaults to length(pm_array)
npp - If Pphi_n_array is not specified, the number of Pϕ_n grid points to be mapped onto (max/min Pϕ_n will be automatically determined). Defaults to length(Rm_array)
needJac - If set to true, the algorithm will assume that jacobian is needed when mapping between spaces (e.g. for distributions)
verbose - If set to true, the function will talk a lot!

Please note! This version of the function is NOT the most efficient. It could be improved so as to avoid having to run pm2Rm_2_μPϕ() when the Lambda_array
and Pphi_n_array keyword arguments are specified. This might be done in future versions of the OWCF.
"""
function os2COM(M::AbstractEquilibrium, good_coords::Vector{CartesianIndex{3}}, data::Array{Float64, 3}, E_array::AbstractVector, pm_array::AbstractVector, 
                Rm_array::AbstractVector, FI_species::AbstractString; isTopoMap::Bool=false, Lambda_array::Union{Nothing,AbstractVector}=nothing, 
                Pphi_n_array::Union{Nothing,AbstractVector}=nothing, nl::Int64=length(pm_array), npp::Int64=length(Rm_array), needJac::Bool=false, 
                verbose::Bool=false, kwargs...)

    data_COM_raw = zeros(length(E_array),nl,npp,2)
    μ_matrix = zeros(length(E_array),nl) # Create matrix, because we need to save all possible μ-values for all possible energies (the min/max values of μ scale with the energy)
    Pϕ_matrix = zeros(length(E_array),npp) # Create matrix, because we need to save all possible Pϕ-values for all possible energies (the min/max values of Pϕ scale with the energy)
    for (iE, E) in enumerate(E_array)
        verbose && println("Mapping (E,pm,Rm)->(E,μ,Pϕ;σ) for energy slice $(iE) of $(length(E_array))... ")
        good_coords_iE = findall(x-> x[1]==iE, good_coords) # Returns a 1D vector with Integer elements
        good_coords_pmRm = map(x-> CartesianIndex(x[2],x[3]), good_coords[good_coords_iE]) # Returns a vector with CartesianIndex{2} elements
        data_COM_raw[iE,:,:,:], E, μ_array, Pϕ_array = pmRm_2_μPϕ(M, good_coords_pmRm, data[iE,:,:], E, pm_array, Rm_array, FI_species; isTopoMap=isTopoMap, nμ=nl, nPϕ=npp, needJac=needJac, verbose=verbose, kwargs...)
        μ_matrix[iE,:] = μ_array # Save the array of possible μ-values FOR THE CURRENT ENERGY E in a matrix at row iE
        Pϕ_matrix[iE,:] = Pϕ_array # Save the array of possible Pϕ-values FOR THE CURRENT ENERGY E in a matrix at row iE
    end

    B0 = norm(Equilibrium.Bfield(M,magnetic_axis(M)...))
    if isnothing(Lambda_array)
        verbose && println("'Lambda_array' input to os2COM() not specified. Computing from μ by using B0=$(B0) T... ")
        Λ_matrix = (B0 ./((1000*GuidingCenterOrbits.e0) .*E_array)) .*μ_matrix
        min_Lambda, max_Lambda = extrema(Λ_matrix)
        Lambda_array = collect(range(min_Lambda,stop=max_Lambda,length=nl))
    end

    q = getSpeciesCharge(FI_species)
    psi_axis, psi_bdry = psi_limits(M)
    if psi_bdry==0
        @warn "The magnetic flux at the last closed flux surface (LCFS) is found to be 0 for the 'M' input to os2COM(). Pϕ_n=Pϕ/(q*|Ψ_w|) where Ψ_w=Ψ(mag. axis) is assumed instead of Ψ_w=Ψ(LCFS)."
        Ψ_w_norm = abs(psi_axis)
    else
        Ψ_w_norm = abs(psi_bdry)
    end

    if isnothing(Pphi_n_array)
        verbose && println("'Pphi_n_array' input to os2COM() not specified. Computing from Pϕ by using q=$(q) Coulomb and |Ψ_w|=$(Ψ_w_norm)... ")
        Pϕ_n_matrix = inv(q*Ψ_w_norm) .*Pϕ_matrix
        min_Pphi_n, max_Pphi_n = extrema(Pϕ_n_matrix)
        Pphi_n_array = collect(range(min_Pphi_n,stop=max_Pphi_n,length=npp))
    end

    interp_func = isTopoMap ? Interpolations.Constant() : Interpolations.Linear() # For topological maps, nearest neighbour method should be used
    data_COM = zeros(length(E_array),length(Lambda_array),length(Pphi_n_array),2)
    for (iE,E) in enumerate(E_array)
        verbose && println("Mapping (E,μ,Pϕ;σ)->(E,Λ,Pϕ_n;σ) for energy slice $(iE) of $(length(E_array))... ")
        E_joule = (1000*GuidingCenterOrbits.e0)*E
        nodes = ((B0/E_joule) .*μ_matrix[iE,:], (1/(q*Ψ_w_norm)) .*Pϕ_matrix[iE,:])
        jac = needJac ? (E_joule*q*Ψ_w_norm)/B0 : 1.0
        data_COM_ctgo_raw_itp = Interpolations.interpolate(nodes,data_COM_raw[iE,:,:,1],Gridded(interp_func))
        data_COM_ctgo_raw_etp = Interpolations.extrapolate(data_COM_ctgo_raw_itp,0.0)
        data_COM_cogo_raw_itp = Interpolations.interpolate(nodes,data_COM_raw[iE,:,:,2],Gridded(interp_func))
        data_COM_cogo_raw_etp = Interpolations.extrapolate(data_COM_cogo_raw_itp,0.0)
        for (il,Λ) in enumerate(Lambda_array)
            for (ip,Pϕ_n) in enumerate(Pphi_n_array)
                data_COM[iE,il,ip,1] = jac*data_COM_ctgo_raw_etp(Λ,Pϕ_n)
                data_COM[iE,il,ip,2] = jac*data_COM_cogo_raw_etp(Λ,Pϕ_n)
            end
        end
    end

    if isTopoMap
        data_COM[data_COM .== 0] .= 9 # Invalid orbit value is 9, not 0
        data_COM = Int64.(data_COM)
    end

    return data_COM, E_array, Lambda_array, Pphi_n_array
end

"""
    os2COM(M, data, E_array, pm_array, Rm_array, FI_species)
    os2COM(M, data, E_array, pm_array, Rm_array, FI_species; nl=length(pm_array), npp=length(Rm_array), isTopoMap=false, needJac=false, verbose=false, 
                                                             vverbose=false, good_coords=nothing, wall=nothing, extra_kw_args=Dict(:toa => true, :limit_phi => true, :max_tries => 0), 
                                                             kwargs...)

This function maps an orbit-space (E,pm,Rm) quantity into COM-space (E,Λ,Pϕ_n;σ). E is the energy (keV), pm is the pitch (v_||/v) at the 
maximum major radius 'Rm' position of the orbit. Λ=μ*B_0/E where μ is the magnetic moment and B_0 is the norm of the magnetic field vector (T) at the 
magnetic axis. Pϕ_n=Pϕ/(q*Ψ_w) where Pϕ is the toroidal canonical angular momentum, q is the particle charge (Coulomb) and Ψ_w is the value of the (poloidal)
magnetic flux at the last closed flux surface (unless it is 0 due to some convention, then the value of the poloidal magnetic flux at the magnetic axis should be used).
If a 4D quantity is given as input 'data', the function will assume it's a collection of 3D (E,pm,Rm) quantities with the last three dimensions corresponding to E, pm and Rm. 
Return transformed data in COM, along with E, Λ, Pϕ_n arrays. The output data will thus be either 4D or 5D. The last dimension of the output data will correspond to 
the binary σ coordinate (-1 or +1) that keeps track of counter- and co-going orbits with the same (E, Λ, Pϕ_n) coordinate. The first Julia index (1) corresponds 
to σ=-1 and the second Julia index (2) corresponds to σ=+1.

If the input data is not a topological map and the 'good_coords' keyword argument has not been specified, the function will identify the (E,pm,Rm) coordinates suitable for 
mapping by first computing the corresponding orbit grid.

The input arguments are:
M - The magnetic equilibrium object
data - The data (a function, or several functions, of (E,pm,Rm)) to be mapped to (E,Λ,Pϕ_n;σ)
E_array - The energy grid points. In keV
pm_array - The pitch maximum grid points
Rm_array - The maximum major radius grid points. In meters
FI_species - The particle species (of the fast ion). Please see OWCF/misc/species_func.jl for more info

The keyword arguments are:
nl - The number of Λ grid points. Defaults to length(pm_array)
npp - The number of Pϕ_n grid points. Defaults to length(Rm_array)
isTopoMap - If set to true, the os2COM() function will assume the input argument 'data' is a topological map. Defaults to false
verbose - If set to true, the function will talk a lot!
vverbose - If set to true, the function will REALLY talk a lot!
good_coords - The (E,pm,Rm) coordinates suitable (valid) for mapping. Defaults to nothing (unknown)
wall - The wall/boundary object to be included when integrating the equations of motion to compute the orbit grid, if needed. Defaults to nothing (unknown)
extra_kw_args - Extra keyword arguments to be included when integrating the -||-.

PLEASE NOTE! If you want to map a function that needs a Jacobian (e.g. a fast-ion distribution), you need to supply a keyword argument called 'needJac' and 
specify it to true.
"""
function os2COM(M::AbstractEquilibrium, data::Union{Array{Float64, 3},Array{Float64, 4}}, E_array::AbstractVector, pm_array::AbstractVector, Rm_array::AbstractVector, 
                FI_species::AbstractString; nl::Int64=length(pm_array), npp::Int64=length(Rm_array), isTopoMap::Bool=false, verbose::Bool=false, vverbose::Bool=false, 
                good_coords=nothing, wall::Union{Nothing,Boundary{Float64}}=nothing, extra_kw_args=Dict(:toa => true, :limit_phi => true, :max_tries => 0), kwargs...)
    if !isTopoMap && !(typeof(good_coords)==Vector{CartesianIndex{3}})
        verbose && println("Input data is not a topological map, and keyword argument 'good_coords' has not been (correctly?) provided... ")
        verbose && println("--------> Computing orbit grid for (E,pm,Rm)-values to be able to deduce valid orbits for (E,pm,Rm)-grid... ")
        if !verbose
            oldstd = stdout
            redirect_stdout(devnull) # Redirect prints to null, to avoid orbit_grid progress bar messing up log files
        end
        orbs, og = OrbitTomography.orbit_grid(M, E_array, pm_array, Rm_array; q=getSpeciesEcu(FI_species), amu=getSpeciesAmu(FI_species), wall=wall, extra_kw_args...)
        if !verbose
            redirect_stdout(oldstd) # Start printing outputs again
        end
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

    data_COM = zeros(size(data,1),length(E_array),nl,npp,2) # Pre-allocate the output data
    verbose && println("Mapping (E,pm,Rm) -> (E,Λ,Pϕ_n;σ) for 3D quantity 1 of $(size(data,1))... ")
    data_COM[1, :, :, :, :], E_array, Λ_array, Pϕ_n_array = os2COM(M, good_coords, data[1,:,:,:], E_array, pm_array, Rm_array, FI_species; nl=nl, npp=npp, isTopoMap=isTopoMap, verbose=verbose, vverbose=vverbose, kwargs...)
    if size(data,1)>1 # To avoid distributed errors
        data_COM_mp = @showprogress 1 "Mapping (E,pm,Rm) -> (E,Λ,Pϕ_n;σ)... " @distributed (+) for iEd in collect(2:size(data,1))
            #verbose && println("Mapping (E,pm,Rm) -> (E,Λ,Pϕ_n;σ) for 3D quantity $(iEd) of $(size(data,1))... ")
            data_COM_part = zeros(size(data,1),length(E_array),nl,npp,2)
            # b,c,d do not matter. But they need to be returned by os2COM.
            data_COM_part[iEd, :, :, :, :], b, c, d = os2COM(M, good_coords, data[iEd,:,:,:], E_array, pm_array, Rm_array, FI_species; Lambda_array=Λ_array, Pphi_n_array=Pϕ_n_array, nl=nl, npp=npp, isTopoMap=isTopoMap, verbose=false, vverbose=false, kwargs...)
            data_COM_part # Declare ready for reduction (data_COM += data_COM_part but distributed over several CPU cores)
        end

        data_COM = data_COM + data_COM_mp # Maybe inefficient/non-optimal coding, but I don't care. It works!
    end

    if !input_4D
        verbose && println("Reshaping output data from 5D to 4D... ")
        return dropdims(data_COM,dims=1), E_array, Λ_array, Pϕ_n_array 
    end

    if isTopoMap
        data_COM = Int64.(data_COM) # Convert to Int64, just in case
    end
    return data_COM, E_array, Λ_array, Pϕ_n_array
end

"""
    com2OS()

PLEASE NOTE! THIS FUNCTION IS UNDER ACTIVE DEVELOPMENT!!!
PLEASE NOTE! THIS FUNCTION IS UNDER ACTIVE DEVELOPMENT!!!
PLEASE NOTE! THIS FUNCTION IS UNDER ACTIVE DEVELOPMENT!!!
"""
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

    # If data is NOT a topological map, we need to extract valid (E,pm,Rm) coordinates to map to (E,mu,Pphi,sigma)
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
        @warn "The interpolation algorithm failed to interpolate $(numObadInds) (E,mu,Pphi,sigma) points."
    end

    if debug
        return data_OS, E_array, pm_array, Rm_array, Hc_OS, valid_cogo_OScoords
    end

    return data_OS, E_array, pm_array, Rm_array
end

"""
    com2OS(M, data, E_array, mu_matrix, Pphi_matrix, FI_species)
    com2OS(M, data, E_array, mu_matrix, Pphi_matrix, FI_species; npm=length(mu_matrix[1,:]), nRm=length(Pphi_matrix[1,:]), wall=nothing, dataIsTopoMap=false, topoMap=nothing, verbose=false, kwargs...)

PLEASE NOTE! THIS FUNCTION IS UNDER ACTIVE DEVELOPMENT!!!
PLEASE NOTE! THIS FUNCTION IS UNDER ACTIVE DEVELOPMENT!!!
PLEASE NOTE! THIS FUNCTION IS UNDER ACTIVE DEVELOPMENT!!!

This function maps a constants-of-motion (E,mu,Pphi,sigma) quantity into orbit space (E,pm,Rm). The input data is given as a 5D quantity, and the function will assume it's a collection of 4D (E,mu,Pphi,sigma) quantities with
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
the number of invalid/impossible (E,mu,Pphi,sigma) coordinates.

kwargs... => dataIsTopoMap, needJac, vverbose, debug, extra_kw_args

dataIsTopoMap - Should be set to true if 'data' is a topological map.
needJac - Should be set to true if mapping a quantity that requires a Jacobian e.g. a fast-ion distribution.
vverbose - Should be set to true if you want the algorithm to be VERY talkative!
debug - If set to true, debug mode will be active.
extra_kw_args - A dictionary with extra keyword arguments for the orbit integration algorithm found in GuidingCenterOrbits.jl/orbit.jl, used to finalize (E,mu,Phi,sigma) -> (E,pm,Rm) transform.
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
    Ep2VparaVperp(E_array, p_array, Q; my_gcp=GCDeuteron, needJac=false, isTopoMap=false, verbose=false, returnAbscissas=false)

Transform a 2D (E,p) quantity Q into (v_para,v_perp) space. If needJac, assume a jacobian is needed (e.g distribution).
If isTopoMap, assume Q is a topological map, implying special care is needed for successful transform.
If returnAbscissas, then the vpara, vperp grid points will be returned as separate output variables. That is:

    Q_vel, vpara_vector, vperp_vector = Ep2VparaVperp(E_array, p_array, Q; returnAbscissas=true)

Otherwise, please use as:

    Q_vel = Ep2VparaVperp(E_array, p_array, Q)

PLEASE NOTE! E_array MUST be given in keV.
"""
function Ep2VparaVperp(E_array::Vector, p_array::Vector, Q::Matrix{Float64}; my_gcp::AbstractParticle{T}=GCDeuteron(0.0,0.0,0.0,0.0), needJac::Bool=false, isTopoMap::Bool=false, verbose::Bool=false, returnAbscissas::Bool=false) where {T}
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

    #brutetree = BruteTree(hcat(E_rep,p_rep)')
    for (ivpara, vpara) in enumerate(vpara_array), (ivperp, vperp) in enumerate(vperp_array)
        v = sqrt(vpara^2 + vperp^2)
        Eq = ((my_gcp.m)*(GuidingCenterOrbits.c0)^2)*(sqrt(1/(1-(v/(GuidingCenterOrbits.c0))^(2)))-1) # Relativistic energy
        pq = vpara/v
        
        #if (Eq < min_keV_E_array) || (Eq > max_keV_E_array) || (pq < min_p_array) || (pq > max_p_array) # If outside of bounds
            #verbose && println("Point out-of-bounds found!")
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

    if returnAbscissas
        return Q_VEL, vpara_array, vperp_array
    else
        return Q_VEL
    end
end

"""
_slowing_down_function_core(v_0, p_0, v_array, p_array, n_e, T_e, species_f, species_th_vec, n_th_vec, T_th_vec)
_slowing_down_function_core(-||-; S0=1.0, g=(v-> v), dampen=false, damp_type=:erfc, sigma=nothing, v_damp=nothing, returnExtra=false)

Compute a neutral beam injection slowing-down distribution function, using the formulas in W.G.F. Core, Nucl. Fusion 33, 829 (1993)
and J.D. Gaffey, J. Plasma Physics 16, 149-169, (1976).
Function arguments:
- v_0 - Injection speed [m/s]
- p_0 - Injection pitch 
- v_array - Speed grid points [m/s]
- p_array - Pitch grid points
- n_e - Electron density [m^-3]
- T_e: The electron temperature [keV]
- species_f: The (fast) beam particle species [-]
- species_th_vec: A vector with all thermal particle species, e.g. ["D", "T"] [-]
- n_th_vec: A vector with all thermal particle densities, e.g. [1.0e20, 2.0e20] [m^-3]
- T_th_vec: A vector with all thermal particle temperatures, e.g. [7.0, 5.0] [keV]
Keyword arguments:
- S0 - Source amplitude at injection point
- g - Function of speed (v). exp(-g(v)) governs speed diffusion above injection speed
- dampen - If true, the slowing-down function will be damped below v_L (please see final eq. in W.G.F. Core, Nucl. Fusion 33, 829 (1993)) using the error function complement
- damp_type - If 'dampen' is set to true, damp_type specifies the type of damping. :linear and :erfc (complementory error function) currently supported
- sigma - If 'damp_type' is set to :erfc, sigma is the sigma in the Gaussian distribution
- v_damp - If specified, the slowing-down function (SDF) will be dampened below the speed of v_damp (in m/s). If not specified, 'dampen' is set to true and 'damp_type' is set to :erfc, the SDF will be dampened below v_L (see function code).
- returnExtra - If true, (in addition to f_SD) v_c, v_L, τ_s and α_0 will be returned

This function should not be used directly, but always via the slowing_down_function() function, to avoid confusing between 
energy (E) and speed (v) grid points input.
"""
function _slowing_down_function_core(v_0::Real, p_0::Real, v_array::Vector{T} where {T<:Real}, p_array::Vector{T} where {T<:Real}, n_e::Real, T_e::Real, 
                                     species_f::String, species_th_vec::Vector{String}, n_th_vec::Vector{T} where {T<:Real}, T_th_vec::Vector{T} where {T<:Real}; 
                                     S0::Float64=1.0, g=(v-> v), dampen::Bool=false, damp_type::Symbol=:erfc, sigma::Union{Nothing,Real}=nothing, 
                                     v_damp::Union{Nothing,Real}=nothing, returnExtra::Bool=false)
    p_0 = clamp(p_0, -0.999, 0.999) # Clamp the injection pitch, to avoid singularities at (-1.0, 1.0)
    m_e = (GuidingCenterOrbits.e_amu)*(GuidingCenterOrbits.mass_u) # Electron mass, kg
    τ_s = spitzer_slowdown_time(n_e, T_e, species_f, species_th_vec, n_th_vec, T_th_vec) # Spitzer slowing-down time, s
    Z_th_vec = getSpeciesEcu.(species_th_vec) # Atomic charge number for all thermal plasma ion species
    M_th_vec = getSpeciesMass.(species_th_vec) # Mass for all thermal plasma ion species, kg
    Z_1 = (getSpeciesMass(species_f)/n_e)*reduce(+, n_th_vec .*Z_th_vec .*Z_th_vec ./M_th_vec) # Slowing-down constant. 
    E_e = T_e*1000*GuidingCenterOrbits.e0 # Thermal electron energy, keV -> J
    v_e = (GuidingCenterOrbits.c0)*sqrt(1-(1/((E_e/(m_e*(GuidingCenterOrbits.c0)^2))+1)^2))# Thermal electron (relativistic) speed, m/s
    v_c = (((3*sqrt(pi)*m_e*Z_1)/(4*getSpeciesMass(species_f)))^(1/3))*v_e # Cross-over speed, m/s
    Z_2 = (1/n_e)*(1/Z_1)*reduce(+, n_th_vec .*Z_th_vec .*Z_th_vec) # Slowing-down constant
    β = Z_2/2 # A Fokker-planck plasma parameter
    α_0 = β*(1-p_0^2)/3 # W.G.F Core parameter
    f_SD = zeros(length(v_array),length(p_array)) # Pre-allocate the slowing-down distribution function
    for (iv,v) in enumerate(v_array), (ip,p) in enumerate(p_array)
        α = (v > v_0 ? -1.0 : 1.0) * α_0 * log((1+(v_c/v)^3)/(1+(v_c/v_0)^3))
        F_SD = (S0*τ_s/(v^3 + v_c^3))*exp(-(p-p_0)^2 /(4*α))/(2*sqrt(pi*α)) 
        f_SD[iv,ip] = F_SD * exp(-(v > v_0 ? g(v) : 0)) # Heavyside function with g(v)
    end
    if returnExtra || dampen
        v_L = v_c * ((1+(v_c/v_0)^3)*exp((3/(4*β))*((1-abs(p_0))/(1+abs(p_0))))-1)^(-1/3)
    end
    if dampen
        v_damp = isnothing(v_damp) ? v_L : v_damp # If v_damp is specified, use v_damp. If not, use v_L
        nv = length(v_array)
        if damp_type==:erfc && (v_damp>0.0)
            sigma = !isnothing(sigma) ? sigma : (v_damp/10)/(2*sqrt(2*log(2))) # If sigma was not specified, use formula based on the formula for FWHM for Gaussians
            damp_factors = reshape(erfc.((-1) .*(v_array .- v_damp); sigma=sigma) ./2,(nv,1)) # Factors range between 0 and 1
            f_SD = damp_factors .*f_SD
        elseif damp_type==:linear
            damp_factors = reshape(collect(range(eps(),stop=1,length=nv)),(nv,1))
            f_SD = damp_factors .*f_SD
        else
            if !(damp_type in [:erfc, :linear])
                error("Unknown damp_type. Currently valid types are :erfc and :linear. Please correct and re-try.")
            end
        end
    end
    if returnExtra
        return f_SD, v_c, v_L, τ_s, α_0
    end
    return f_SD
end

"""
slowing_down_function(E_0, p_0, E_array, p_array, n_e, T_e, species_f, species_th_vec, n_th_vec, T_th_vec)
slowing_down_function(-||-; type=:core, returnExtra=false, E_tail_length=nothing, kwargs...)

Compute a neutral beam injection slowing-down distribution function. The inputs are 
    - E_0: The beam injection energy [keV]
    - p_0: The beam pitch [-]
    - E_array: The energy grid points [keV]
    - p_array: The pitch grid points [-]
    - n_e: The electron density [m^-3]
    - T_e: The electron temperature [keV]
    - species_f: The (fast) beam particle species [-]
    - species_th_vec: A vector with all thermal particle species, e.g. ["D", "T"] [-]
    - n_th_vec: A vector with all thermal particle densities, e.g. [1.0e20, 2.0e20] [m^-3]
    - T_th_vec: A vector with all thermal particle temperatures, e.g. [7.0, 5.0] [keV]
The keyword arguments are:
    - type: The type of slowing-down function you want to compute (currently, only :core is supported) [-]
    - returnExtra: If true, extra outputs will be returned (please see below)
    - E_tail_length: If specified (Float or Int), the slowing-down function will be dampened below (E0-E_tail_length) [keV]
    - sigma: If specified, (Float or Int), the slowing-down function will be dampened this quickly [keV]
"""
function slowing_down_function(E_0::Real, p_0::Real, E_array::Vector{T} where {T<:Real}, p_array::Vector{T} where {T<:Real}, n_e::Real, T_e::Real, 
                               species_f::String, species_th_vec::Vector{String}, n_th_vec::Vector{T} where {T<:Real}, T_th_vec::Vector{T} where {T<:Real}; 
                               type::Symbol=:core, returnExtra::Bool=false, E_tail_length::Union{Nothing,Real}=nothing, sigma::Union{Nothing,Real}=nothing, 
                               kwargs...)
    m_f = getSpeciesMass(species_f) # Beam (fast) particle mass, kg
    E2v_rel = (E-> (GuidingCenterOrbits.c0)*sqrt(1-(1/(((E*1000*GuidingCenterOrbits.e0)/(m_f*(GuidingCenterOrbits.c0)^2))+1)^2))) # A one-line function to transform from energy (keV) to relativistic speed (m/s)
    v2E_rel = (v-> inv(1000*GuidingCenterOrbits.e0)*m_f*((GuidingCenterOrbits.c0)^2)*(sqrt(1/(1-(v/(GuidingCenterOrbits.c0))^(2)))-1)) # A one-line function to transform from relativistic speed (m/s) to energy (keV)

    v_damp = nothing
    dampen = false
    if !isnothing(E_tail_length)
        v_damp = E2v_rel(clamp(E_0 - E_tail_length,0,Inf))
        dampen = true
    end

    v_sigma = nothing
    if !isnothing(sigma)
        v_sigma = E2v_rel(clamp(sigma,0,Inf))
    end

    # Use the specified type of slowing-down function
    if type==:core
        my_SD_func = _slowing_down_function_core
    else
        error("slowing_down_function currently only supports W.G.F. Core type. Please correct and re-try.")
    end

    # If the source point is right on a grid point...
    E_array_orig = deepcopy(E_array)
    p_array_orig = deepcopy(p_array)
    if E_0 in E_array_orig || p_0 in p_array_orig
        if E_0 in E_array_orig
            E_array = collect(range(minimum(E_array_orig),stop=maximum(E_array_orig), length=(length(E_array_orig)+1))) # Increase number of energy grid points by 1, to avoid having grid on E_0
        end
        if p_0 in p_array_orig
            p_array = collect(range(minimum(p_array_orig),stop=maximum(p_array_orig), length=(length(p_array_orig)+1))) # Increase number of pitch grid points by 1, to avoid having grid on p_0
        end
    end
    v_0 =  E2v_rel(E_0) # Convert energy source point to velocity (relativistically)
    v_array = E2v_rel.(E_array) # Convert energy grid to velocity (relativistically)
    if returnExtra # If extra output is desired... 
        f_SD_vp, v_c, v_L, τ_s, α_0 = my_SD_func(v_0, p_0, v_array, p_array, n_e, T_e, species_f, species_th_vec, n_th_vec, T_th_vec; returnExtra=returnExtra, dampen=dampen, v_damp=v_damp, sigma=v_sigma, kwargs...)
    else # If not...
        f_SD_vp = my_SD_func(v_0, p_0, v_array, p_array, n_e, T_e, species_f, species_th_vec, n_th_vec, T_th_vec; dampen=dampen, v_damp=v_damp, sigma=v_sigma, kwargs...)
    end
    dvdE_array = [(1/sqrt(2*m_f*v2E_rel(v))) for v in v_array] # Compute Jacobian from v to E
    f_SD = reshape(dvdE_array,(length(dvdE_array),1)) .*f_SD_vp # Transform f(v,p) to f(E,p)
    if E_0 in E_array_orig || p_0 in p_array_orig # If there were any source points on grid
        nodes = (E_array,p_array)
        itp = Interpolations.interpolate(nodes,f_SD,Gridded(Linear()))
        f_SD_orig = [itp(E,p) for E in E_array_orig, p in p_array_orig] # Interpolate onto original (E,p) grid (having avoided singularities when computing)
    else
        f_SD_orig = f_SD
    end
    if returnExtra
        return f_SD_orig, v_c, v_L, τ_s, α_0 # Return critical velocity, lower-threshold velocity, Spitzer slowing-down time and α_0 parameter
    else
        return f_SD_orig # Return only slowing-down function in (E,p) coordinates
    end
end