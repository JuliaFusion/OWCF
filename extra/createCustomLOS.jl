###################################### createCustomLOS.jl ##########################################

#### Description:
# This script allows the user of the OWCF to create a custom line-of-sight (LOS) to use as a 
# diagnostic for computing weight functions. The LOS will be saved as a .vc file that can be used by 
# other OWCF script, often as the 'filepath_diagnostic' input variable. The LOS can be constructed in 
# any of the following 5 ways:
# Case 1: The LOS is constructed using a Cartesian (x,y,z) vector as the 'LOS_vec' input variable
#         (further explained below) and a Cartesian point for the 'detector_location' input variable.
#         The LOS_vec is then a vector (does not need to be a unit vector) that runs along the LOS,
#         TOWARDS the detector location.
# Case 2: The LOS is constructed using a cylindrical (R,phi,z) vector as the 'LOS_vec' input variable
#         and a cylindrical point for the 'detector_location' input variable. The LOS_vec is then a vector 
#         (does not need to be a unit vector) that runs along the LOS, TOWARDS the detector location. 
# Case 3: The LOS is constructed using two angles (θ_x, θ_z) as the 'LOS_vec' input variable and a 
#         Cartesian point for the 'detector_location' input variable. The θ_x angle is the angle 
#         in degrees between the (projection of the) LOS vector and the x-axis in the (x,y) plane. 
#         The θ_z angle is the angle in degrees between the LOS vector and the z-axis. The LOS_vec is 
#         then a vector (does not need to be a unit vector) that runs along the LOS, TOWARDS the detector 
#         location.
# Case 4: The LOS is constructed using two angles (θ_R, phi) as the 'LOS_vec' input variable and a 
#         cylindrical point for the 'detector_location' input variable. The θ_R angle is the angle 
#         in degrees between the LOS vector and the R-axis in the (R,z) plane. The phi angle is the 
#         standard cylindrical coordinate angle. The LOS vector is assumed to run along the LOS, 
#         TOWARDS the detector location.
# Case 5: The LOS is constructed using a single angle θ_u as the 'LOS_vec' (a 1-element vector). θ_u
#         is the angle between the LOS vector and the magnetic field at an (R,z) point of interest.
#         This (R,z) point has to be specified (further explanation below). The 'detector_location'
#         is inferred automatically in this case, and can be left unspecified.

#### Input variables:
# The input variables are defined in e.g. the start_createCustomLOS_template.jl file. Please see 
# OWCF/templates/start_createCustomLOS_template.jl for explanation.
#
# You run the createCustomLOS.jl script by making a copy of the start_createCustomLOS_template.jl 
# file, moving it to the OWCF/ folder, specifying the inputs and executing it by executing 
# 'julia start_createCustomLOS_template.jl' or equivalent.

#### Outputs:
# A .vc file modelling the LOS. The .vc file follows the structure set originally by the LINE21 code.
# The rows of the file correspond to LOS voxels. The columns of the file correspond to different 
# LOS quantities for each voxel. They are as follows.
# - The first column corresponds to Cartesian x points
# - The second column corresponds to Cartesian y points
# - The third column corresponds to Cartesian/Cylindrical z points
# - The fourth column corresponds to redundant weights C
# - The fifth column corresponds to voxel volumes V
# - The sixth column corresponds to x components of the vector U pointing towards the detector
# - The seventh column corresponds to y components of the vector U pointing towards the detector
# - The eigth column corresponds to z components of the vector U pointing towards the detector
# - The ninth column corresponds to major radius R points
# - The tenth column corresponds to toroidal phi angle points
# - The eleventh column corresponds to solid angles OMEGA

# Script written by Henrik Järleblad. Last maintained 2025-06-03.
###############################################################################################################

## ---------------------------------------------------------------------------------------------
println("Loading Julia packages... ")
@everywhere begin
    using Dates
    using FileIO
    using JLD2
    using LinearAlgebra
    using Equilibrium
    plot_LOS && (using Plots)
    debug = debug # This is always set to false, except when OWCF developers are debugging
    debug && (using Plots)
end

## ---------------------------------------------------------------------------------------------
# Determine which case it is
verbose && print("Determining custom LOS build case... ")
case = :UNKNOWN
if collimated
    possible_cases = Array{Int64,1}()
    if coordinate_system==:cartesian
        push!(possible_cases,1)
        push!(possible_cases,3)
    elseif coordinate_system==:cylindrical
        push!(possible_cases,2)
        push!(possible_cases,4)
    elseif coordinate_system==:magrel
        verbose && println("It's 5!")
        case = 5
    else
        error("Unknown coordinate system. Please correct and re-try. Use :cartesian, :cylindrical or :magrel.")
    end
    if case==:UNKNOWN
        numOcoords = length(LOS_vec)
        if numOcoords==2
            push!(possible_cases,3)
            push!(possible_cases,4)
        elseif numOcoords==3
            push!(possible_cases,1)
            push!(possible_cases,2)
        else
            error("Unknown number of coordinates. Please correct and re-try. If coordinate_system is not set to :magrel, please use 2 (angles) or 3 (coordinates) to specify 'LOS_vec'.")
        end
        # Find which case it is. 1, 2, 3 or 4?
        case = argmax([length(findall(x-> x==0, possible_cases .- 1)),length(findall(x-> x==0, possible_cases .- 2)),length(findall(x-> x==0, possible_cases .- 3)),length(findall(x-> x==0, possible_cases .- 4))])
    end
else
    if !(coordinate_system in [:cartesian, :cylindrical])
        error("Input variable 'collimated' was set to $(collimated) but input variable 'coordinate_system' was set to $(coordinate_system). When 'collimated=false', currently supported options for 'coordinate_system' include :cartesian and :cylindrical. Please correct and re-try.")
    end
    if coordinate_system==:cartesian
        case = 1
    end
    if coordinate_system==:cylindrical
        case = 2
    end
end
verbose && println("It's $(case)!")
## ---------------------------------------------------------------------------------------------
# Loading tokamak equilibrium
verbose && println("Loading tokamak equilibrium... ")
if ((split(filepath_equil,"."))[end] == "eqdsk") || ((split(filepath_equil,"."))[end] == "geqdsk")
    M, wall = read_geqdsk(filepath_equil,clockwise_phi=false) # Assume counter-clockwise phi-direction
else # Otherwise, assume magnetic equilibrium is a saved .jld2 file
    myfile = jldopen(filepath_equil,false,false,false,IOStream)
    M, wall = myfile["S"], myfile["wall"]
    close(myfile)
end
wall_dR = maximum(wall.r)-minimum(wall.r)
wall_dz = maximum(wall.z)-minimum(wall.z)

## ---------------------------------------------------------------------------------------------
# Define the function that returns a rotation matrix that can rotate a 3D vector an angle 'theta' about an arbitrary axis (defined by the unit vector 'u')
verbose && println("Defining rotation matrix helper function... ")
"""
    getRotationMatrix(u,theta)
    getRotationMatrix(-||-; verbose=false)

Return a rotation matrix that can rotate a 3D vector an angle 'theta' about an arbitrary axis (defined by the unit vector 'u'
"""
function getRotationMatrix(u::Vector{Float64},theta::Float64; verbose=false)
    if length(u)!=3
        @warn "getRotationMatrix() not supported for other dimensionalities than 3 at the moment. Returning zero matrix."
        return zeros(length(u),length(u))
    end
    R11 = cos(theta)+u[1]*u[1]*(1-cos(theta))
    R12 = u[1]*u[2]*(1-cos(theta))-u[3]*sin(theta)
    R13 = u[1]*u[3]*(1-cos(theta))+u[2]*sin(theta)
    R21 = u[2]*u[1]*(1-cos(theta))+u[3]*sin(theta)
    R22 = cos(theta)+u[2]*u[2]*(1-cos(theta))
    R23 = u[2]*u[3]*(1-cos(theta))-u[1]*sin(theta)
    R31 = u[3]*u[1]*(1-cos(theta))-u[2]*sin(theta)
    R32 = u[3]*u[2]*(1-cos(theta))+u[1]*sin(theta)
    R33 = cos(theta)+u[3]*u[3]*(1-cos(theta))

    R = [[R11,R21,R31] [R12,R22,R32] [R13,R23,R33]]
    return R
end

## ---------------------------------------------------------------------------------------------
# Define functions that help to model the diagnostic LOS
verbose && println("Defining helper functions... ")

"""
    getPhiFromXandY(x,y)

Compute the cylindrical coordinate system coordinate phi (in radians), from the Cartesian coordinate system coordinates x and y.
"""
function getPhiFromXandY(x::Number, y::Number)
    if x==0 && y==0
        return 0.0 # Just define atan(0/0)=:0, for functionality's sake
    elseif x==0 && !(y==0)
        return (pi/2)*y/abs(y)
    elseif x>0.0
        return atan(y/x)
    elseif x<0.0 && y>=0.0
        return atan(y/x)+pi
    elseif x<0.0 && y<0.0
        return atan(y/x)-pi
    else
        error("This should not be possible! Fatal error! Re-start Julia please.")
    end
end

"""
    visible_by_detector(wall::Boundary, point::AbstractVector, detector_location::AbstractVector)

Given the tokamak wall 'wall', deduce if the Cartesian point 'point' is visible by the detector 
at the Cartesian location 'detector_location'. If it is, return true. Otherwise, return false.
The point is deemed not visible if:
    - The point is outside of the tokamak wall
    - The point is inside  of the tokamak wall and a straight line 's' from the point to the detector 
      location passes through a region outside of the tokamak wall where ∂s/∂R<0, and the region has a 
      major radius value smaller than the detector location and the minimum of the tokamak wall

Keyword arguments include:
    - verbose - If set to true, the function will talk a lot! - Bool

The AbstractEquilibrium struct is the AbstractEquilibrium struct in the Equilibrium.jl package.
The Boundary struct is the Boundary struct in the Equilibrium.jl package.
"""
function visible_by_detector(wall::Boundary, point::AbstractVector, detector_location::AbstractVector; verbose::Bool=false)

    R_point = sqrt(point[1]^2 + point[2]^2) # R of point
    if !in_boundary(wall,R_point,point[3]) # If the point is outside of the tokamak wall, it is deemed not visible
        verbose && println("visible_by_detector(): Point $(round.([R_point,point[3]],digits=2)) is outside of the tokamak first wall. Returning false... ")
        return false
    end

    min_R_wall = minimum(wall.r) # The minimum of the tokamak first wall
    R_detector = sqrt(detector_location[1]^2 + detector_location[2]^2) # R of detector location

    U = detector_location - point # A vector pointing from the point to the detector
    distance_to_detector = norm(U) # The distance (in meters) from the point to the detector location
    u = U ./distance_to_detector # Unit vector pointing from voxel to detector
    
    prev_R = -Inf
    for l in round.(cumsum(0.1 .*ones(length(0.0:0.1:distance_to_detector)-1)),digits=3) # For every decimeter (except for the first one) from the point to the detector location...
        point_to_check = point + (l .*u)
        R = sqrt(point_to_check[1]*point_to_check[1] + point_to_check[2]*point_to_check[2])
        if !in_boundary(wall, R, point_to_check[3]) && (R-prev_R)<0 && R<R_detector && R<min_R_wall # This is an approximation. It should really be R_wall(z) where R_wall is the R value as a function of z for R<argmax(z(R))
            verbose && println("visible_by_detector(): Point $(round.([R,point_to_check[3]],digits=2)) is outside of the tokamak wall and has ∂s/∂R<0 (R: $(round(R,digits=2)), prev_R: $(round(prev_R,digits=2))). Returning false... ")
            return false
        end
        prev_R = R
    end
    return true
end
## ---------------------------------------------------------------------------------------------
verbose && println("Creating uniform (R,phi,z) grid for storing LOS data... ")
R_grid = collect(range(minimum(wall.r),stop=maximum(wall.r),length=nR)) # Major radius in meters
phi_grid = (pi/180) .*collect(range(-180,stop=180,length=nphi))[1:end-1] # Toroidal angle in radians. Exclude phi=pi (since it is equivalent to phi=-pi)
z_grid = collect(range(minimum(wall.z),stop=maximum(wall.z),length=nz)) # Vertical coordinate in meters
ΔR = diff(R_grid)[1]; Δphi = diff(phi_grid)[1]; Δz = diff(z_grid)[1] # The grid spacings
voxel_volumes = map(c-> ΔR*Δz*R_grid[c[1]]*Δphi, CartesianIndices((nR,nphi-1,nz))) # The volume of each voxel
omega_3D_array = zeros(nR,nphi-1,nz) # A 3D array to store all solid angle values for each point of the uniform (R,phi,z) grid
if case==1
    verbose && println("LOS vector and detector location specified in Cartesian coordinates. Taken as is...")
    LOS_vec = LOS_vec
    detector_location = detector_location
elseif case==2
    verbose && println("LOS vector and detector location given in (R,phi,z). Transforming to (x,y,z) for easy computation... ")
    LOS_vec = !collimated ? [1.0,1.0,1.0] : LOS_vec # If !collimated, LOS_vec does not matter
    LOS_vec_x = LOS_vec[1]*cos(LOS_vec[2]*pi/180)
    LOS_vec_y = LOS_vec[1]*sin(LOS_vec[2]*pi/180)
    detector_location_x = detector_location[1]*cos(detector_location[2]*pi/180)
    detector_location_y = detector_location[1]*sin(detector_location[2]*pi/180)
    LOS_vec = [LOS_vec_x, LOS_vec_y, LOS_vec[3]] # Cartesian
    detector_location = [detector_location_x, detector_location_y, detector_location[3]] # Cartesian
elseif case==3
    verbose && println("LOS specified as angles w.r.t. the x- and z-axes. Computing (x,y,z) LOS vector... ")
    LOS_vec_x = sin(LOS_vec[2]*pi/180)*cos(LOS_vec[1]*pi/180)
    LOS_vec_y = sin(LOS_vec[2]*pi/180)*sin(LOS_vec[1]*pi/180)
    LOS_vec_z = cos(LOS_vec[2]*pi/180)
    LOS_vec = [LOS_vec_x, LOS_vec_y, LOS_vec_z] # Cartesian
    detector_location = detector_location
elseif case==4
    verbose && println("LOS specified as phi-angle and angle w.r.t. the R-axis. Computing (x,y,z) LOS vector... ")
    LOS_vec_x = cos(LOS_vec[1]*pi/180) * cos(LOS_vec[2]*pi/180)
    LOS_vec_y = cos(LOS_vec[1]*pi/180) * sin(LOS_vec[2]*pi/180)
    LOS_vec_z = sin(LOS_vec[1]*pi/180)
    LOS_vec = [LOS_vec_x, LOS_vec_y, LOS_vec_z] # Cartesian
    detector_location_x = detector_location[1]*cos(detector_location[2]*pi/180)
    detector_location_y = detector_location[1]*sin(detector_location[2]*pi/180)
    detector_location = [detector_location_x, detector_location_y, detector_location[3]] # Cartesian
else # case must be 5
    verbose && println("LOS specified as angle w.r.t. the magnetic field. Computing (x,y,z) LOS vecotrs... ")
    # The most complicated one...
    if R_of_interest==:mag_axis
        R_of_interest = magnetic_axis(M)[1]
    end
    if !(typeof(R_of_interest)==Float64 || typeof(R_of_interest)==Int64)
        error("R_of_interest specified incorrectly! Please correct and re-try.")
    end
    if z_of_interest==:mag_axis
        z_of_interest = magnetic_axis(M)[2]
    end
    if !(typeof(z_of_interest)==Float64 || typeof(z_of_interest)==Int64)
        error("z_of_interest specified incorrectly! Please correct and re-try.")
    end
    B_vec = Bfield(M, R_of_interest, z_of_interest)
    R = getRotationMatrix([0.0,1.0,0.0], LOS_vec[1] *(pi/180)) # Get a rotation matrix that can rotate any vector θ_u degrees
    b_vec = B_vec ./ norm(B_vec)
    LOS_vec = R*b_vec # Rotate B-field unit vector θ_u degrees around the y_axis to get LOS_vec with θ_u degrees relative to B-field

    # Now, assume we are at phi=0. Because of toroidal symmetry, it does not matter where we are.
    # At phi=0, the (x,y,z) and (R,phi,z) coordinate systems overlap.
    # We can thus perform the following.
    LOS_vec_x = LOS_vec[1]
    LOS_vec_y = LOS_vec[2]
    LOS_vec_z = LOS_vec[3]
    LOS_vec = [LOS_vec_x, LOS_vec_y, LOS_vec_z]
    x_of_interest = R_of_interest * cos(0.0)
    y_of_interest = R_of_interest * sin(0.0)
    p_of_interest = [x_of_interest, y_of_interest, z_of_interest]
    detector_location = p_of_interest + (LOS_length/2)*LOS_vec
end
detector_location_R = sqrt(detector_location[1]^2 + detector_location[2]^2) # For plotting purposes

## ---------------------------------------------------------------------------------------------
# Check how parallel the line-of-sight vector is to one of the cylindrical coordinate directions.
# If it is too parallel, we need special care to successfully map out the LOS voxels
if collimated
    LOS_vec_R = sqrt(LOS_vec[1]*LOS_vec[1]+LOS_vec[2]*LOS_vec[2])
    LOS_vec_phi = round(getPhiFromXandY(LOS_vec[1],LOS_vec[2])*180/pi)
    LOS_vec_cyl = [LOS_vec_R, LOS_vec_phi, LOS_vec[3]]
    dot_R = dot(LOS_vec_cyl,[1.0,0.0,0.0])/norm(LOS_vec_cyl)
    dot_phi = dot(LOS_vec_cyl,[0.0,1.0,0.0])/norm(LOS_vec_cyl)
    dot_z = dot(LOS_vec_cyl,[0.0,0.0,1.0])/norm(LOS_vec_cyl)

    verbose && println("---> LOS vector (towards detector, cartesian   coordinates): $(round.(LOS_vec,digits=3))")
    verbose && println("---> LOS vector (towards detector, cylindrical coordinates (phi in degrees)): $(LOS_vec_cyl)")

    LOS_parallel_to_R = false
    if isapprox(dot_R,1.0)
        verbose && println("LOS parallel to R-axis! Special algorithm will be utilized to successfully compute the LOS.")
        LOS_parallel_to_R = true
        prev_R = sqrt(detector_location[1]*detector_location[1]+detector_location[2]*detector_location[2]) # Will be needed if LOS is too parallel to the R axis
    end
    LOS_parallel_to_phi = false
    if isapprox(dot_phi,1.0)
        verbose && println("LOS parallel to phi-axis! Special algorithm will be utilized to successfully compute the LOS.")
        LOS_parallel_to_phi = true
        prev_phi = getPhiFromXandY(detector_location[1],detector_location[2]) # Will be needed if LOS is too parallel to the phi axis
    end
    LOS_parallel_to_z = false
    if isapprox(dot_z,1.0)
        verbose && println("LOS parallel to z-axis! Special algorithm will be utilized to successfully compute the LOS.")
        LOS_parallel_to_z = true
        prev_z = detector_location[3] # Will be needed if LOS is too parallel to the z axis
    end
end
## ---------------------------------------------------------------------------------------------
# Set the resolution of the LOS, based on the smallest possible distance resolvable by the (R,phi,z) grid
if collimated
    δL = minimum([ΔR, minimum(R_grid)*Δphi, Δz]) # Smallest possible distance resolvable by the (R,phi,z) grid
    LOS_length_res = Int64(round(LOS_length/δL))
    if LOS_length_res<1
        error("The LOS was found to be shorter than the smallest possible distance resolvable by the (R,phi,z) grid. Please increase the values of the nR, nz and nphi input variables to be able to create the LOS.")
    end
    LOS_azimuthal_res = Int64(round(pi*LOS_width/δL))
    if LOS_azimuthal_res<1
        error("The LOS was found to be smaller than the smallest possible distance resolvable by the (R,phi,z) grid. Please increase the values of the nR, nz and nphi input variables to be able to create the LOS.")
    end
    LOS_radial_res = Int64(round(0.5*LOS_width/δL))
    if LOS_radial_res<1
        error("The LOS was found to be thinner than the smallest possible distance resolvable by the (R,phi,z) grid. Please increase the values of the nR, nz and nphi input variables to be able to create the LOS.")
    end
    verbose && println("LOS will be resolved in... ")
    verbose && println("---> $(LOS_length_res) points along the LOS")
    verbose && println("---> $(LOS_azimuthal_res) points around the LOS")
    verbose && println("---> $(2*LOS_radial_res) points across the (width of the) LOS")
end

## ---------------------------------------------------------------------------------------------
# Create LOS via 3D array storing solid angles, which represents a uniform (R,phi,z) grid
# If user requested plot of the LOS, or debugging, set the plot font to Computer modern
if plot_LOS || debug
    plot_font = "Computer Modern"
    my_plt = Plots.default(fontfamily=plot_font)
end
# If debugging, pre-allocate a Vector for storing plot frames to create debug .gif of LOS creation process
if debug
    plt_frames = Vector{Plots.Plot}(undef,101)
    plt_ind = 1
end

verbose && println("Creating LOS... ")
if collimated
    # Compute the unit vector in the direction of the detector. Using that vector, go stepwise from the detector along 
    # the length of the LOS. Transform the points of the LOS to (R,phi,z) to find LOS on uniform (R,phi,z) grid (LINE21 format requirement)
    # Check which points are inside the tokamak wall/plasma boundary, and keep only those
    y_hat = [0.0,1.0,0.0] # The Cartesian y unit vector. Arbitrary direction for Gram-Schmidt process.
    PERP2LOS_vec = y_hat - (dot(y_hat,LOS_vec)/dot(LOS_vec,LOS_vec))*LOS_vec # Gram-Schmidt process. To acquire vector perpendicular to line-of-sight
    perp2los_vec = PERP2LOS_vec ./norm(PERP2LOS_vec) # Get the unit vector perpendicular to the line-of-sight
    los_vec = LOS_vec ./norm(LOS_vec) # Get the normalized line-of-sight vector
    if LOS_azimuthal_res==1
        @warn "'LOS_azimuthal_res' automatically set to 1. LOS might be badly resolved."
        EDGE_ofLOS_unitvecs = [perp2los_vec]
    else
        LOS_circ_angles = collect(range(0.0,stop=2*pi,length=LOS_azimuthal_res)) # Angles to rotate around spine of LOS
        LOS_circ_angles = LOS_circ_angles[1:end-1] # Remove redundant 2*pi angle
        EDGE_ofLOS_unitvecs = Vector{Vector{Float64}}(undef,LOS_azimuthal_res-1)
        for (ia,angle) in enumerate(LOS_circ_angles)
            local R
            R = getRotationMatrix(los_vec,angle)
            EDGE_ofLOS_unitvec = R*perp2los_vec # Get the unit vectors pointing outwards to the edge of the cylindrical LOS, relative to the radial center (spine) of the LOS
            EDGE_ofLOS_unitvecs[ia] = inv(norm(EDGE_ofLOS_unitvec)) .*EDGE_ofLOS_unitvec # Should already be unit vector, but just in case of numerical inaccuracies
        end
    end

    s = collect(range(0.0,stop=LOS_length,length=LOS_length_res)) # Create increments for walking along the LOS using the normalized LOS vector
    ds = s[2] # Since s[1] is zero
    prog_proc = []
    for is=2:LOS_length_res # Skip first increment, since that would result in LOS voxels behind the detector
        debug_plot = false

        ss = s[is] # The length increment

        # The "spine" is the center of the LOS cylinder (radially)
        LOS_spine_point = detector_location - (ss-ds)*los_vec # LOS_vec is pointing towards detector. To go away from detector, use -, not +.
        spine_x = LOS_spine_point[1]; spine_y = LOS_spine_point[2]; spine_z = LOS_spine_point[3] # (x,y,z) of LOS spine point
        spine_R = sqrt(spine_x*spine_x+spine_y*spine_y) # R = sqrt(x^2 + y^2)

        # Progress bar (1 % to 100 %) and debug plot bool (every 1 %)
        if !(floor(100*is/LOS_length_res) in prog_proc)
            append!(prog_proc,floor(100*is/LOS_length_res))
            verbose && println("Creating diagnostic line-of-sight $(prog_proc[end]) %... ")
            debug && (debug_plot = true)
        end

        # Is the spine point not visible by the detector (outside the tokamak wall or obscured by the central solenoid)?
        if !visible_by_detector(wall,[spine_x, spine_y, spine_z],detector_location; verbose=debug)
            continue # Skip this increment ('is') 
        end

        # Spine point phi (in radians) and (iR,iphi,iz) index
        spine_phi = getPhiFromXandY(spine_x, spine_y) # In radians
        spine_iR = Int64(round(inv(ΔR)*(spine_R-R_grid[1])+1))
        spine_iphi = Int64(round(inv(Δphi)*(spine_phi-phi_grid[1])+1))
        spine_iz = Int64(round(inv(Δz)*(spine_z-z_grid[1])+1))

        spine_distance_to_detector = ss-ds # Distance (in meters) from the spine point to the detector
        if (ss-ds)==0 # If it's the first point along the LOS (should never be, since the detector would be then be at the edge of the machine/tokamak first wall)
            spine_distance_to_detector = eps() # Avoid divide by 0
        end
        spine_omega = inv(spine_distance_to_detector^2)*(voxel_volumes[spine_iR, spine_iphi, spine_iz])^(2/3) # Ω = A/r^2 where A is the area of the voxel (approximate by volume^(2/3)) and r is the distance to the detector
        omega_3D_array[spine_iR, spine_iphi, spine_iz] = (spine_omega > 2*pi ? 2*pi : spine_omega) # Cap the solid angle value to 2*pi. However, this should almost never be necessary

        for r_radial in range(0.0,stop=LOS_width/2,length=LOS_radial_res)[2:end] # Radially discretize the LOS into LOS_radial_res number of points
            for EDGE_ofLOS_unitvec in EDGE_ofLOS_unitvecs # Azimuthally discretize the LOS into LOS_azimuthal_res number of points
                LOS_point = LOS_spine_point + r_radial .*EDGE_ofLOS_unitvec # An (x,y,z) point inside the LOS

                local x; local y; local z; local R # Declare local scope, to avoid annoying warnings
                x = LOS_point[1]; y = LOS_point[2]; z = LOS_point[3]; R = sqrt(x^2 + y^2)
                if !in_boundary(wall,R,z) || !visible_by_detector(wall, [x,y,z], detector_location) # Is the point outside the tokamak wall/boundary? Or is it not visible by the detector (obscured by the central solenoid)?
                    continue # Ignore this LOS point
                end
                local phi; # Declare local scope, to avoid annoying warnings
                phi = getPhiFromXandY(x,y)

                iR = Int64(round(inv(ΔR)*(R-R_grid[1])+1))
                iphi = Int64(round(inv(Δphi)*(phi-phi_grid[1])+1))
                iz = Int64(round(inv(Δz)*(z-z_grid[1])+1))

                distance_to_detector = norm(LOS_point - detector_location)

                omega = inv(distance_to_detector^2)*(voxel_volumes[iR, iphi, iz])^(2/3)

                omega_3D_array[iR, iphi, iz] = (omega > 2*pi ? 2*pi : omega) # Cap the solid angle value to 2*pi. However, this should almost never be necessary
            end
        end

        # Process plotting, for debug purposes
        if debug_plot
            global my_plt; global plt_frames; global plt_ind
            LOS_Rphi = dropdims(sum(omega_3D_array,dims=3),dims=3)
            LOS_Rz = dropdims(sum(omega_3D_array,dims=2),dims=2)
            LOS_phiz = dropdims(sum(omega_3D_array,dims=1),dims=1)

            plt_Rphi = Plots.heatmap(R_grid, (180/pi) .*phi_grid, transpose(LOS_Rphi), colorbar=false, label="", xlabel="R [m]", ylabel="phi [degrees]", fillcolor=:blues)
            plt_Rz = Plots.heatmap(R_grid, z_grid, transpose(LOS_Rz), colorbar=false, label="", xlabel="R [m]", ylabel="z [-]", fillcolor=:blues)
            #plt_Rz = Plots.scatter!([detector_location_R],[detector_location[3]], label="Detector location", markershape=:star, markercolor=:purple, markerstrokewidth=2)
            plt_Rz = Plots.plot!(wall.r,wall.z,label="Tokamak first wall",linewidth=2.5,color=:gray)
            plt_Rz = Plots.scatter!([magnetic_axis(M)[1]],[magnetic_axis(M)[2]],label="Mag. axis",markershape=:xcross,markercolor=:red,markerstrokewidth=4)
            plt_Rz = Plots.plot!(title="Creating LOS... ", aspect_ratio=:equal, xlims=(minimum(wall.r)-0.1*wall_dR,maximum(wall.r)+wall_dR), ylims=extrema(wall.z))
            plt_phiz = Plots.heatmap(z_grid, (180/pi) .*phi_grid, LOS_phiz, colorbar=false, label="", xlabel="z [m]", ylabel="phi [degrees]", fillcolor=:blues)
            plt2 = Plots.plot(plt_Rphi, plt_phiz)
            my_plt = Plots.plot(plt_Rz, plt2, layout=(1,2))
            plt_frames[plt_ind] = my_plt; plt_ind += 1
        end
    end
else
    verbose && println("LOS not collimated. Creating LOS by checking which (R,phi,z) voxels are obscured by the central solenoid... ")
    mycoords = CartesianIndices((nR,nphi-1,nz))
    prog_proc = []
    for (icoord,coord) in enumerate(mycoords)
        debug_plot=false

        # Progress bar (1 % to 100 %) for ever 1 %
        if !(floor(100*icoord/length(mycoords)) in prog_proc)
            append!(prog_proc,floor(100*icoord/length(mycoords)))
            verbose && println("Creating diagnostic line-of-sight $(prog_proc[end]) %... ")
            debug && (debug_plot=true)
        end
        
        local x; local y; local z; local R; local phi
        iR, iphi, iz = coord[1], coord[2], coord[3]
        R, phi, z = R_grid[iR], phi_grid[iphi], z_grid[iz]
        x, y = R*cos(phi), R*sin(phi)

        LOS_point = [x,y,z]

        if visible_by_detector(wall, LOS_point, detector_location)
            distance_to_detector = norm(LOS_point - detector_location)
            omega = inv(distance_to_detector^2)*(voxel_volumes[iR, iphi, iz])^(2/3)
            omega_3D_array[iR,iphi,iz] = (omega > 2*pi ? 2*pi : omega) # Cap the solid angle value to 2*pi. However, this should almost never be necessary
        end

        # Process plotting, for debug purposes
        if debug_plot
            global my_plt; global plt_frames; global plt_ind
            LOS_Rphi = dropdims(sum(omega_3D_array,dims=3),dims=3)
            LOS_Rz = dropdims(sum(omega_3D_array,dims=2),dims=2)
            LOS_phiz = dropdims(sum(omega_3D_array,dims=1),dims=1)

            plt_Rphi = Plots.heatmap(R_grid, (180/pi) .*phi_grid, transpose(LOS_Rphi), colorbar=false, label="", xlabel="R [m]", ylabel="phi [degrees]", fillcolor=:blues)
            plt_Rz = Plots.heatmap(R_grid, z_grid, transpose(LOS_Rz), colorbar=false, label="", xlabel="R [m]", ylabel="z [-]", fillcolor=:blues)
            #plt_Rz = Plots.scatter!([detector_location_R],[detector_location[3]], label="Detector location", markershape=:star, markercolor=:purple, markerstrokewidth=2)
            plt_Rz = Plots.plot!(wall.r,wall.z,label="Tokamak first wall",linewidth=2.5,color=:gray)
            plt_Rz = Plots.scatter!([magnetic_axis(M)[1]],[magnetic_axis(M)[2]],label="Mag. axis",markershape=:xcross,markercolor=:red,markerstrokewidth=4)
            plt_Rz = Plots.plot!(title="Creating LOS... ", aspect_ratio=:equal, xlims=(minimum(wall.r)-0.1*wall_dR,maximum(wall.r)+wall_dR), ylims=extrema(wall.z))
            plt_phiz = Plots.heatmap(z_grid, (180/pi) .*phi_grid, LOS_phiz, colorbar=false, label="", xlabel="z [m]", ylabel="phi [degrees]", fillcolor=:blues)
            plt2 = Plots.plot(plt_Rphi, plt_phiz)
            my_plt = Plots.plot(plt_Rz, plt2, layout=(1,2))
            plt_frames[plt_ind] = my_plt; plt_ind += 1
        end
    end
end

if debug
    verbose && println("--- DEBUG ---> Saving .gif file of LOS creation process as $(folderpath_o)$(LOS_name)_debug.gif... ")
    anim = @animate for frame in plt_frames[1:plt_ind-1]
        anim_plt = Plots.plot(frame, dpi=100)
        anim_plt
    end
    gif(anim, folderpath_o*"$(LOS_name)_debug.gif",fps=5)
end

## ---------------------------------------------------------------------------------------------
# Use non-zero elements of omega_3D_array to map out LOS on uniform (R,phi,z) grid
# Compute x, y, V, UX, UY, and UZ for those points
# Together with R, phi and OMEGA (and C, redundant), put in data structures to be saved in LINE21-like file
# Filter out points that are inside the tokamak wall/plasma boundary and/or not obscured by the central solenoid
verbose && println("Mapping out LOS data on uniform (R,phi,z) grid... ")

# Find all non-zero omega_3D_array indices
nz_coords = findall(x-> x>0.0, omega_3D_array)
verbose && println("---> Number of LOS voxels with eligible (R,phi,z) coordinates (inside grid boundaries): $(length(nz_coords))")

# Instantiate arrays for x, y, z, ... etc
LOS_x = zeros(length(nz_coords))
LOS_y = zeros(length(nz_coords))
LOS_z = zeros(length(nz_coords))
LOS_C = ones(length(nz_coords)) # Does not matter
LOS_V = zeros(length(nz_coords))
LOS_UX = zeros(length(nz_coords))
LOS_UY = zeros(length(nz_coords))
LOS_UZ = zeros(length(nz_coords))
LOS_R = zeros(length(nz_coords))
LOS_phi = zeros(length(nz_coords))
LOS_OMEGA = zeros(length(nz_coords))
LOS_in = Array{Bool,1}(undef,length(nz_coords)) # To keep track of points inside tokamak wall/plasma boundary

# Loop through non-zero indices
prog_proc = []
for (inz, nz_coord) in enumerate(nz_coords)
    local R; local phi; local z; local x; local y

    iR, iphi, iz = nz_coord[1], nz_coord[2], nz_coord[3]
    R, phi, z = R_grid[iR], phi_grid[iphi], z_grid[iz]

    x, y = R*cos(phi), R*sin(phi)
    V = voxel_volumes[iR,iphi,iz] # Volume of cylindrical voxel

    if !CTS_like
        U = detector_location - [x,y,z]
        U = U ./norm(U) # Unit vector pointing from voxel to detector
    else
        U = los_vec # Approximate all voxels as having the same vector pointing towards the detector
    end
    OMEGA = omega_3D_array[nz_coord]

    LOS_x[inz] = x
    LOS_y[inz] = y
    LOS_z[inz] = z
    LOS_V[inz] = V
    LOS_UX[inz] = U[1]
    LOS_UY[inz] = U[2]
    LOS_UZ[inz] = U[3]
    LOS_R[inz] = R
    LOS_phi[inz] = phi
    LOS_OMEGA[inz] = OMEGA
    LOS_in[inz] = visible_by_detector(wall, [x,y,z], detector_location)

    if !(floor(100*inz/length(nz_coords)) in prog_proc)
        append!(prog_proc,floor(100*inz/length(nz_coords)))
        verbose && println("Removing LOS voxels outside tokamak first wall and/or obscured by central solenoid $(prog_proc[end]) %... ")
        visible_by_detector(wall, [x,y,z], detector_location; verbose=debug)
    end
end
verbose && println("---> Number of valid LOS voxels (not outside tokamak first wall and not obscured by central solenoid): $(sum(LOS_in))")

## ---------------------------------------------------------------------------------------------
# Remove voxels outside of the tokamak wall/plasma boundary. Round to 9 significant digits (like the LINE21 code)
LOS_x = round.(LOS_x[LOS_in],sigdigits=9)
LOS_y = round.(LOS_y[LOS_in],sigdigits=9)
LOS_z = round.(LOS_z[LOS_in],sigdigits=9)
LOS_C = round.(LOS_C[LOS_in],sigdigits=9)
LOS_V = round.(LOS_V[LOS_in],sigdigits=9)
LOS_UX = round.(LOS_UX[LOS_in],sigdigits=9)
LOS_UY = round.(LOS_UY[LOS_in],sigdigits=9)
LOS_UZ = round.(LOS_UZ[LOS_in],sigdigits=9)
LOS_R = round.(LOS_R[LOS_in],sigdigits=9)
LOS_phi = round.(LOS_phi[LOS_in],sigdigits=9)
LOS_OMEGA = round.(LOS_OMEGA[LOS_in],sigdigits=9)

if debug
    for i in eachindex(LOS_x)
        if !in_boundary(wall,LOS_R[i],LOS_z[i])
            println("This (R,z) point is outside of the tokamak first wall: ($(round(LOS_R[i],sigdigits=3)),$(round(LOS_z[i],sigdigits=3))))")
        end
    end
end

## ---------------------------------------------------------------------------------------------
# Determining output file name
filepath_output_orig = folderpath_o*LOS_name
filepath_output = deepcopy(filepath_output_orig)
count = 1
while isfile(filepath_output*".vc") # To take care of not overwriting files. Add _(1), _(2) etc
    global filepath_output; global count
    filepath_output = filepath_output_orig*"_($(Int64(count)))"
    count += 1 # global scope, to surpress warnings
end

## ---------------------------------------------------------------------------------------------
# Plot LOS, if requested
date_and_time = split("$(Dates.now())","T")[1]*"at"*split("$(Dates.now())","T")[2][1:5] # Might be needed immidiately below, and definitely when saving output data
if plot_LOS
    verbose && println("Plotting top view (x,y) and poloidal projection (R,z) of LOS... ")
    LOS_heatmap_res = 100
    flux_R = range(extrema(wall.r)...,length=LOS_heatmap_res); dR = diff(flux_R)[1]
    flux_z = range(extrema(wall.z)...,length=LOS_heatmap_res); dz = diff(flux_z)[1]
    inds = CartesianIndices((length(flux_R),length(flux_z)))
    psi_rz = [M(flux_R[ind[1]], flux_z[ind[2]]) for ind in inds]
    psi_mag, psi_bdry = psi_limits(M)

    R_hfs = minimum(wall.r) # R-coord of high-field side wall
    R_lfs = maximum(wall.r) # R-coord of low-field side wall
    wall_phi = collect(0:1:359).*(pi/180.0) # Toroidal angle
    topview_R_hfs_x = (R_hfs).*cos.(wall_phi)
    topview_R_hfs_y = (R_hfs).*sin.(wall_phi)
    topview_R_lfs_x = (R_lfs).*cos.(wall_phi)
    topview_R_lfs_y = (R_lfs).*sin.(wall_phi)

    flux_x = range(extrema(topview_R_lfs_x)...,length=LOS_heatmap_res); dx = diff(flux_x)[1]
    flux_y = range(extrema(topview_R_lfs_y)...,length=LOS_heatmap_res); dy = diff(flux_y)[1]

    LOS_Rz_proj = psi_mag .*ones(LOS_heatmap_res,LOS_heatmap_res) # Multiplying with psi_mag, psy_bdry purely for plotting reasons (to be able to use heatmap and contour in same figure)
    LOS_xy_proj = psi_mag .*ones(LOS_heatmap_res,LOS_heatmap_res)
    for i in eachindex(LOS_R)
        iR = Int64(round(inv(dR)*(LOS_R[i]-flux_R[1])+1))
        iz = Int64(round(inv(dz)*(LOS_z[i]-flux_z[1])+1))
        ix = Int64(round(inv(dx)*(LOS_x[i]-flux_x[1])+1))
        iy = Int64(round(inv(dy)*(LOS_y[i]-flux_y[1])+1))
        if iR>0 && iR<=LOS_heatmap_res && iz>0 && iz<=LOS_heatmap_res
            LOS_Rz_proj[iR,iz] = psi_bdry
        end
        if ix>0 && ix<=LOS_heatmap_res && iy>0 && iy<=LOS_heatmap_res
            LOS_xy_proj[ix,iy] = psi_bdry
        end
    end

    # Use correct color ordering for LOS visualization
    color_array = [:white, :green]
    dpsi = abs(psi_mag - psi_bdry)
    clims = (psi_mag-dpsi*0.01, psi_bdry+dpsi*0.01)
    if (psi_bdry-psi_mag)<0 # If COCOS ID 3, 4, 7 or 8 (O. Sauter and S. Yu. Medvedev, Computer Physics Communications 184 (2013) 293)
        color_array = [:green, :white]
        clims = (psi_bdry-dpsi*0.01, psi_mag+dpsi*0.01)
    end

    plt_crs = Plots.heatmap(flux_R, flux_z, transpose(LOS_Rz_proj), title="Poloidal proj. $(LOS_name) LOS", fillcolor=cgrad(color_array, categorical=true), colorbar=false)
    plt_crs = Plots.contour!(flux_R, flux_z, transpose(psi_rz), levels=collect(range(psi_mag, stop=psi_bdry,length=7)), color=:lightgray, clims=clims, linewidth=2.5, label="", colorbar=false)
    plt_crs = Plots.plot!(wall.r,wall.z,label="Tokamak first wall",linewidth=2.5,color=:black)
    #plt_crs = Plots.scatter!([detector_location_R],[detector_location[3]], label="Detector location", markershape=:star, markercolor=:purple, markerstrokewidth=2)
    plt_crs = Plots.scatter!([magnetic_axis(M)[1]],[magnetic_axis(M)[2]],label="Mag. axis",markershape=:xcross,markercolor=:red,markerstrokewidth=4)
    plt_crs = Plots.plot!(aspect_ratio=:equal,xlabel="R [m]",ylabel="z [m]", xlims=(minimum(wall.r)-0.1*wall_dR,maximum(wall.r)+wall_dR), ylims=extrema(wall.z))
    plt_crs = Plots.plot!(xtickfontsize=14,ytickfontsize=14,xguidefontsize=16,yguidefontsize=16)
    plt_crs = Plots.plot!(legend=:bottomright,legendfontsize=13)

    plt_top = Plots.heatmap(flux_x, flux_y, transpose(LOS_xy_proj), title="Top view $(date_and_time)", fillcolor=cgrad(color_array, categorical=true), colorbar=false)
    plt_top = Plots.scatter!([detector_location[1]],[detector_location[2]], label="Detector location", markershape=:star, markercolor=:purple, markerstrokewidth=2)
    plt_top = Plots.plot!(topview_R_hfs_x,topview_R_hfs_y,linewidth=2.5,color=:black,label="")
    plt_top = Plots.plot!(topview_R_lfs_x,topview_R_lfs_y,linewidth=2.5,color=:black,label="")
    plt_top = Plots.plot!(aspect_ratio=:equal,xlabel="x [m]",ylabel="y [m]")
    plt_top = Plots.plot!(xtickfontsize=14,ytickfontsize=14,xguidefontsize=16,yguidefontsize=16)

    myplt = Plots.plot(plt_top,plt_crs,layout=(1,2),left_margin=5Plots.mm,bottom_margin=5Plots.mm, size=(1000,500))
    myplt = Plots.plot!(dpi=100)
    display(myplt)

    if save_LOS_plot
        verbose && println("Saving LOS plot in .png file format... ")
        png(myplt,filepath_output)
    end
end

## ---------------------------------------------------------------------------------------------
# Save the computed custom LOS data as a .vc file (same format as the output of the LINE21 code)
verbose && println("Saving results... ")
global filepath_output = filepath_output*".vc"
myfile = open(filepath_output, "w")
for i=1:length(LOS_x)
    write(myfile,"   $(LOS_x[i]) $(LOS_y[i]) $(LOS_z[i]) $(LOS_C[i]) $(LOS_V[i]) $(LOS_UX[i]) $(LOS_UY[i]) $(LOS_UZ[i]) $(LOS_R[i]) $(LOS_phi[i]) $(LOS_OMEGA[i])\n")
end
close(myfile)

## ---------------------------------------------------------------------------------------------
verbose && println("Your output file can be found at $(filepath_output)")
verbose && println("~~~createCustomLOS.jl finished successfully!~~~")