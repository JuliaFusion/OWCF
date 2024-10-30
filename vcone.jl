########################################## vcone.jl ##################################################
# This script contains structures and functions to be able to load and represent a viewing cone 
# (aka line-of-sight, aka sightline) stored in a file produced by the LINE21 code (or similar).
# It is a direct Julia re-implementation of the vcone.py file in the DRESS code (J. Eriksson et al, 
# Computer Physics Communications 199 (2016) 40–46).
#
# THIS SCRIPT IS FUNCTIONALLY DONE. ONLY DOCUMENTATION REMAINS TO BE FINALIZED.
#
# Written by Henrik Järleblad. Last maintained 2024-10-30
######################################################################################################

using DelimitedFiles
using LinearAlgebra

"""
    ViewingCone

STRUCTURE DESCRIPTION TO BE WRITTEN HERE.
"""
struct ViewingCone
    file::String
    X::Vector{Float64}
    Y::Vector{Float64}
    Z::Vector{Float64}
    C::Vector{Float64}
    V::Vector{Float64}
    UX::Vector{Float64}
    UY::Vector{Float64}
    UZ::Vector{Float64}
    R::Vector{Float64}
    PHI::Vector{Float64}
    OMEGA::Vector{Float64}
    Nvoxels::Int64
    Vtot::Float64
    UR::Vector{Float64}
    Uphi::Vector{Float64}
    U::Matrix{Float64}
    Rvals::Vector{Float64}
    Zvals::Vector{Float64}
    PHIvals::Vector{Float64}
    NR::Int64
    NZ::Int64
    NP::Int64
    NPHI::Int64
    dR::Float64
    dZ::Float64
    dPHI::Float64
    IR::Vector{Int64}
    IZ::Vector{Int64}
    IP::Vector{Int64}
    IPvals::Vector{Int64}
    IVOX::Vector{Vector{Int64}}
    NVOX::Vector{Int64}
    RP::Vector{Float64}
    ZP::Vector{Float64}
    CP::Vector{Float64}
end

"""
    ViewingCone(name)

FUNCTION DESCRIPTION TO BE WRITTEN HERE
"""
function ViewingCone(name::String)
    # Determine file path based on name
    file_path = (lowercase(name) == "km11" || lowercase(name) == "tofor") ? "vc_data/TOFOR/TOFOR.vc" :
                (lowercase(name) == "ab" ? "vc_data/AB/AB.vc" : name)

    # Load data from file
    VC = readdlm(file_path)

    # Assign columns to respective variables
    X, Y, Z = VC[:, 1], VC[:, 2], VC[:, 3]
    C, V = VC[:, 4], VC[:, 5]
    UX, UY, UZ = VC[:, 6], VC[:, 7], VC[:, 8]
    R, PHI = VC[:, 9], VC[:, 10]

    # Calculate OMEGA if not present
    OMEGA = size(VC, 2) == 11 ? VC[:, 11] : ((4*π .*C) ./ V)

    # Basic attributes
    Nvoxels = length(C)
    Vtot = sum(V)

    # Convert detector direction vectors to cylindrical coordinates
    UR = UX .* cos.(PHI) .+ UY .* sin.(PHI)
    Uphi = -UX .* sin.(PHI) .+ UY .* cos.(PHI)
    U = hcat(UR, Uphi, UZ)'

    # Reconstruct grid used by LINE21 and unique values
    Rvals, Zvals, PHIvals = sort(unique(R)), sort(unique(Z)), sort(unique(PHI))
    NR, NZ, NPHI = length(Rvals), length(Zvals), length(PHIvals)
    NP = NR * NZ # Number of poloidal grid points

    # Grid spacing (assumes uniform)
    dR, dZ, dPHI = Rvals[2] - Rvals[1], Zvals[2] - Zvals[1], PHIvals[2] - PHIvals[1]

    # Map each voxel to the corrresponding poloidal grid point
    IR = round.(Int64, (R .- Rvals[1]) ./ dR) # Indices starting at 0, due to easy math operations
    IZ = round.(Int64, (Z .- Zvals[1]) ./ dZ) # Indices starting at 0, due to easy math operations
    IP = IR + NR .*IZ # Indices starting at 0, due to easy math operations
    IPvals = sort(unique(IP)) # Indices starting at 0, due to easy math operations

    # Map each poloidal grid point to the corresponding voxels
    IVOX = [Int64[] for _ in 1:NP] # A Vector of (empty) Vectors of Int64 elements
    for (ivox, ip) in enumerate(IP)
        push!(IVOX[ip + 1], ivox)
    end

    # Count number of voxels at each poloidal grid point
    NVOX = [length(IVOX[i]) for i in 1:NP]

    # Generate poloidal projection of the viewing cone
    IRP, IZP = IPvals .% NR, div.(IPvals, NR) # The remainder is the R index. The quotient is the z index
    RP, ZP = Rvals[IRP .+ 1], Zvals[IZP .+ 1]
    CP = [sum(C[IVOX[ip + 1]]) for ip in IPvals] # Total weight of all voxels at the given poloidal location

    # Correct for Julia array index convention (start at index 1 instead of 0)
    IR .+= 1
    IZ .+= 1
    IP .+= 1
    IPvals .+= 1

    return ViewingCone(file_path, X, Y, Z, C, V, UX, UY, UZ, R, PHI, OMEGA, Nvoxels, Vtot, UR, Uphi, U, Rvals, Zvals, PHIvals,
                       NR, NZ, NP, NPHI, dR, dZ, dPHI, IR, IZ, IP, IPvals, IVOX, NVOX, RP, ZP, CP)
end

"""
get_poloidal_index(vc,r,z)

Method to check poloidal index. 
Check which poloidal bin the input (R,z) is in.
Requires that the poloidal grid spacing is uniform.

WRITE EXPLANATION FOR INPUTS AND OUTPUT
"""
function get_poloidal_index(vc::ViewingCone, R::Union{Real,T}, z::Union{Real,T}) where {T<:AbstractArray}
    R, z = R |> vcat, z |> vcat  # Convert scalars to 1-element Vectors
    ir = Int64.(floor.((R .- vc.Rvals[1] .+ (vc.dR / 2)) ./ vc.dR)) # Indices starting at 0, due to easy math operations
    iz = Int64.(floor.((z .- vc.Zvals[1] .+ (vc.dZ / 2)) ./ vc.dZ)) # Indices starting at 0, due to easy math operations
    ip = ir .+ vc.NR .*iz # Indices starting at 0, due to easy math operations

    # Mark out-of-bounds indices
    ip[ir .< 0 .|| ir .>= vc.NR .|| iz .< 0 .|| iz .>= vc.NZ] .= -1 # Set all poloidal indices with (R,z) points outside the grid to -1
    ip = map(x-> x==-1 ? x : x+1, ip) # Return to using 1 to refer to index of first element in array, but keep -1 elements
    outside = .!(in.(ip, Ref(vc.IP))) # If there somehow were indices created outside of the poloidal index range...
    ip[outside] .= -1 # ...set them to -1 as well
    return ip
end

"""
is_inside(vc,R,z)

Method to check if points are inside the viewing cone

WRITE EXPLANATION FOR INPUTS AND OUTPUT
"""
function is_inside(vc::ViewingCone, R::Union{Real,T}, z::Union{Real,T}) where {T<:AbstractArray}
    ip = get_poloidal_index(vc, R, z)
    return in.(ip, Ref(vc.IP)) # Return an array where 1 means inside the vc, 0 means outside
end

"""
get_voxels(vc,R,z)

Method to get all voxels including a point in the poloidal plane

WRITE EXPLANATION FOR INPUTS AND OUTPUT
"""
function get_voxels(vc::ViewingCone, R::Union{Real,T}, z::Union{Real,T}) where {T<:AbstractArray}
    ip = get_poloidal_index(vc, R, z)
    return map(x-> x == -1 ? -1 : vc.IVOX[x],ip)
end

"""

map_points_voxels(vc,R,z)

Map given (R,z) points to each voxel that encloses them, under the
assumption of toroidal symmetry. 

Inputs
-------
WRITE EXPLANATION FOR INPUTS

Outputs
-------
i_voxels: array of voxel indices. Each element correponds to 
    one point in the voxel with that index.

i_points: array of point indices (same length as i_voxels). 
    For instance, i_point[j] tells us that that point is
    inside the voxel with index i_voxels[j].
    (A point with given R,z coordinates might be inside 
    multiple voxels, since toroidal symmetry is assumed).
"""
function map_points_voxels(vc::ViewingCone, R::Union{Real,T}, z::Union{Real,T}) where {T<:AbstractArray}
    ip = get_poloidal_index(vc, R, z)
    i_all_points = collect(1:length(ip))

    # Filter points inside the viewing cone
    inside = findall(x-> x != -1, ip)
    ip, i_points_inside = ip[inside], i_all_points[inside]

    if isempty(ip)
        # No points inside any voxel
        return Int64[], Int64[]
    end

    # Map points to voxel indices
    ivox = vc.IVOX[ip]
    nvox = vc.NVOX[ip]
    
    i_voxels = vcat(ivox) # Reshape all voxel 1-element arrays into one long vector
    i_points = similar(i_voxels) # Pre-allocate i_points

    # Assign points to voxel indices
    i0 = 1
    for (i, n) in enumerate(nvox)
        i1 = i0 + n - 1
        i_points[i0:i1] .= i_points_inside[i]
        i0 = i1 + 1
    end

    return i_voxels, i_points
end
