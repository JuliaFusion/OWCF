###################################### createCustomLOS.jl ##########################################

#### Description:
# This script allows the user of the OWCF to create a custom line-of-sight (LOS) to use as a 
# diagnostic for computing weight functions. The LOS will be saved as a .vc file that can be used by 
# other OWCF script, often as the 'filepath_diagnostic' input variable. The LOS can be constructed in 
# either of the following 5 ways:
# Case 1: The LOS is constructed using a Cartesian (x,y,z) vector as the 'los_vec' input variable
#         (further explained below) and a Cartesian point for the 'detector_location' input variable.
#         The los_vec is then a vector (does not need to be a unit vector) than runs along the LOS,
#         TOWARDS the detector location.
# Case 2: The LOS is constructed using a cylindrical (R,phi,z) vector as the 'los_vec' input variable
#         and a cylindrical point for the 'detector_location' input variable. The los_vec... -||-.
# Case 3: The LOS is constructed using two angles (θ_x, θ_z) as the 'los_vec' input variable and a 
#         Cartesian point for the 'detector_location' input variable. The θ_x angle is the angle 
#         in degrees between the LOS vector and the x-axis in the (x,y) plane. The θ_z angle is the
#         angle in degrees between the LOS vector and the z-axis in the (x,z) plane.
# Case 4: The LOS is constructed using two angles (θ_R, θ_phi) as the 'los_vec' input variable and a 
#         cylindrical point for the 'detector_location' input variable. The θ_R angle is the angle 
#         in degrees between the LOS vector and the R-axis in the (R,z) plane. The θ_phi angle is the 
#         angle in degrees between the LOS vector and the phi-axis in the (phi,z) plane.
# Case 5: The LOS is constructed using a single angle θ_u as the 'los_vec' (a 1-element vector). θ_u
#         is the angle between the LOS vector and the magnetic field at an (R,z) point of interest.
#         This (R,z) point has to be specified (further explanation below). The 'detector_location'
#         is inferred automatically in this case, and can be left unspecified.

#### Input variables:
# folderpath_OWCF - The file path to the OWCF folder - String 
#
# coordinate_system - A symbol input to tell the script what type of coordinate system the 'los_vec'
#                     and 'detector_location' input variables are using. The accepted values are 
#                     :cartesian, :cylindrical or :magrel. The :cartesian symbol is used in cases 1
#                     and 3 (see above), the :cylindrical symbol is used in cases 2 and 4 (see above)
#                     and the :magrel symbol is used in case 5 (see above). - Symbol
# R_of_interest - The R coordinate of interest, in case the script is run as case 5. Meters. If :mag_axis, the magnetic axis is used. - Float64
# z_of_interest - The z coordinate of interest, in case the script is run as case 5. Meters. If :mag_axis, the magnetic axis is used. - Float64
# detector_location - The location of the detector at the end of the LOS. In cases 1 and 3, this 
#                     is specified as a 3-element Cartesian vector. In cases 2 and 4, this is
#                     specified as a 3-element cylindrical vector. In case 5, this can be left 
#                     unspecified. - Vector{Float64}
# filepath_equil - The file path to the magnetic equilibrium and tokamak wall/plasma boundary. This 
#                  is needed to filter out LOS points outside of the tokamak wall/plasma boundary. 
#                  In case 5, it's also needed to get the magnetic field vector and the (R,z) point 
#                  of interest. - String
# folderpath_o - The folder path to the folder where you want the output to be saved. Remember to 
#                finish folder paths with '/'. - String
# LOS_length - The length of the LOS. The algorithm will start at the detector_location and move 
#              along the LOS a distance equal to LOS_length. After, it will assume there are no more 
#              voxels relevant for the LOS and therefore stop. Specified in meters. - Float64
# LOS_name - The name of the LOS/detector. Purely for aesthetic reasons. - String
# LOS_vec - The vector pointing TOWARDS the detector. Does NOT need to be a unit vector. LOS_vec can
#           be specified in several ways. For case 1, specify as a (x,y,z) vector pointing towards 
#           the detector along the center of the LOS. For case 2, specify as a (R,phi,z) vector 
#           poiting towards the detector along the center of the LOS. For case 3, specify as an array
#           of 2 elements, θ_x and θ_z (explained above). For case 4, specify as an array of 2
#           elements, θ_R and θ_phi (explained above). For case 5, specify as an array of 1 element,
#           θ_u (explained above). - Vector{Float64}
# LOS_width - The width of the LOS. In meters. - Float64
# verbose - If set to true, the script will talk a lot - Bool
 
#### Advanced inputs (should not be changed from default values, unless you know what you're doing):
# LOS_circ_res - The number of points by which to resolve the LOS azimuthally (around the cylinder
#                LOS). - Int64
# LOS_length_res - The number of points by which to resolve the LOS along the LOS - Int64
# nR - The number of major radius grid points for the uniform (R,phi,z) grid onto which to map the LOS - Int64
# nz - The number of vertical grid points for the uniform (R,phi,z) grid onto which to map the LOS - Int64
# nphi - The number of toroidal grid points for the uniform (R,phi,z) grid onto which to map the LOS - Int64

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

# Script written by Henrik Järleblad. Last maintained 2023-09-21.
###################################################################################################

folderpath_OWCF = ""
cd(folderpath_OWCF)
using Pkg
Pkg.activate(".")

## ---------------------------------------------------------------------------------------------
# Inputs
coordinate_system = :magrel
(coordinate_system==:magrel) && ((R_of_interest, z_of_interest)=(:mag_axis, :mag_axis)) # If you want to specify the LOS as an angle relative to the magnetic field, you have to specify an (R,z) point of interest. Accepted values are Float64 and :mag_axis (magnetic axis). 
detector_location = [] # TOFOR (proxy): [2.97, 0.0, 18.8] /// NE213 (proxy): [8.35,2.1,-0.27]. When coordinate_system==:magrel, detector_location can be set to whatever (the script will figure it out automatically from the angle and LOS_length)
filepath_equil = "" # To load tokamak wall/plasma boundary and magnetic equilibrium (for coordinate_system==:magrel)
folderpath_o = ""
LOS_length = 0.0 # Length of line-of-sight. From detector to desired end. Meters. When coordinate_system==:magrel, assume (R,z) point of interest is at half the LOS length
LOS_name = "test_Detector"
LOS_vec = [0.0,0.0,0.0] # TOFOR (proxy) cartesian: [0.005,-0.001,0.9999] /// NE213 (proxy) cartesian: [0.995902688, 0.293844478e-2, -0.903836320e-1]
LOS_width = 1.0 # TOFOR (proxy): 0.27 /// NE213 (proxy): 0.25
verbose = true

# Advanced inputs
LOS_circ_res = 100 # The cylindrical shape that is the line-of-sight will be resolved azimuthally into these many points
LOS_length_res = 20000 # The cylindrical shape that is the line-of-sight will be resolved parallel to LOS_vec into these many points
nR = 200 # The number of R grid points
nz = 200 # The number of z grid points
nphi = 720 # The number of phi grid points

## ---------------------------------------------------------------------------------------------
println("Loading Julia packages... ")
using FileIO
using LinearAlgebra
using Equilibrium

## ---------------------------------------------------------------------------------------------
# Determine which case it is
verbose && print("Determining custom LOS build case... ")
case = :UNKNOWN
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
    verbose && println("It's $(case)!")
end

## ---------------------------------------------------------------------------------------------
# Loading tokamak equilibrium
verbose && println("Loading tokamak equilibrium... ")
if ((split(filepath_equil,"."))[end] == "eqdsk") || ((split(filepath_equil,"."))[end] == "geqdsk")
    M, wall = read_geqdsk(filepath_equil,clockwise_phi=false) # Assume counter-clockwise phi-direction
    jdotb = M.sigma # The sign of the dot product between the plasma current and the magnetic field

    # Extract timepoint information from .eqdsk/.geqdsk file
    eqdsk_array = split(filepath_equil,".")
    XX = (split(eqdsk_array[end-2],"-"))[end] # Assume format ...-XX.YYYY.eqdsk where XX are the seconds and YYYY are the decimals
    YYYY = eqdsk_array[end-1] # Assume format ...-XX.YYYY.eqdsk where XX are the seconds and YYYY are the decimals
    timepoint = XX*","*YYYY # Format XX,YYYY to avoid "." when including in filename of saved output
else # Otherwise, assume magnetic equilibrium is a saved .jld2 file
    myfile = jldopen(filepath_equil,false,false,false,IOStream)
    M = myfile["S"]
    wall = myfile["wall"]
    close(myfile)
    jdotb = (M.sigma_B0)*(M.sigma_Ip)

    if typeof(timepoint)==String && length(split(timepoint,","))==2
        timepoint = timepoint
    else
        timepoint = "TIMELESS" # Unknown timepoint for magnetic equilibrium
    end
end

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
# Define a function that computes phi from x and y
verbose && println("Defining (x,y)-(phi) helper function... ")
"""
    getPhiFromXandY(x,y)

Compute the cylindrical coordinate system coordinate phi, from the Cartesian coordinate system coordinates x and y.
"""
function getPhiFromXandY(x::Number,y::Number)
    if x==0 && y==0
        return NaN
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

## ---------------------------------------------------------------------------------------------
verbose && println("Creating uniform (R,phi,z) grid for storing LOS data... ")
R_grid = collect(range(minimum(wall.r),stop=maximum(wall.r),length=nR))
phi_grid = collect(range(-179,stop=180,length=nphi)).*(pi/180.0) # Toroidal angle
z_grid = collect(range(minimum(wall.z),stop=maximum(wall.z),length=nz))
omega_3D_array = zeros(nR,nphi,nz) # A 3D array to store all solid angle values for each point on the uniform (R,phi,z) grid
if case==1
    verbose && println("LOS vector and detector location specified in Cartesian coordinates. Taken as is...")
    LOS_vec = LOS_vec
    detector_location = detector_location
elseif case==2
    verbose && println("LOS vector and detector location given in (R,phi,z). Transforming to (x,y,z) for easy computation... ")
    LOS_vec_x = LOS_vec[1]*cos(LOS_vec[2]*pi/180)
    LOS_vec_y = LOS_vec[1]*sin(LOS_vec[2]*pi/180)
    detector_location_x = detector_location[1]*cos(detector_location[2]*pi/180)
    detector_location_y = detector_location[1]*sin(detector_location[2]*pi/180)
    LOS_vec = [LOS_vec_x, LOS_vec_y, LOS_vec[3]] # Cartesian
    detector_location = [detector_location_x, detector_location_y, detector_location[3]] # Cartesian
elseif case==3
    verbose && println("LOS specified as angles w.r.t. the x- and z-axes. Computing (x,y,z) LOS vector... ")
    LOS_vec_x = cos(LOS_vec[1]*pi/180)
    LOS_vec_y = sin(LOS_vec[1]*pi/180)
    LOS_vec_z = (LOS_vec[2]==0.0 ? 1.6e16 : LOS_vec_x/tan(LOS_vec[2]*pi/180)) # Avoid dividing by 0
    LOS_vec = [LOS_vec_x, LOS_vec_y, LOS_vec_z] # Cartesian
    detector_location = detector_location
elseif case==4
    verbose && println("LOS specified as angles w.r.t. the R- and phi-axes. Computing (x,y,z) LOS vector... ")
    LOS_vec_R = cos(LOS_vec[1]*pi/180)
    LOS_vec_z = sin(LOS_vec[1]*pi/180)
    LOS_vec_phi = (LOS_vec[2]==0.0 ? 1.6e16 : LOS_vec_z/tan(LOS_vec[2]*pi/180)) # Avoid dividing by 0
    LOS_vec_x = LOS_vec_R * cos(LOS_vec_phi)
    LOS_vec_y = LOS_vec_R * sin(LOS_vec_phi)
    LOS_vec_z = LOS_vec_z
    LOS_vec = [LOS_vec_x, LOS_vec_y, LOS_vec_z] # Cartesian
    detector_location_x = detector_location[1]*cos(detector_location[2]*pi/180)
    detector_location_y = detector_location[1]*sin(detector_location[2]*pi/180)
    detector_location = [detector_location_x, detector_location_y, detector_location[3]] # Cartesian
else # case must be 5
    # The most complicated one...
    if R_of_interest==:mag_axis
        R_of_interest = M.axis[1]
    end
    if !(typeof(R_of_interest)==Float64 || typeof(R_of_interest)==Int64)
        error("R_of_interest specified incorrectly! Please correct and re-try.")
    end
    if z_of_interest==:mag_axis
        z_of_interest = M.axis[2]
    end
    if !(typeof(z_of_interest)==Float64 || typeof(z_of_interest)==Int64)
        error("z_of_interest specified incorrectly! Please correct and re-try.")
    end
    B_vec = Bfield(M, R_of_interest, z_of_interest)
    R = getRotationMatrix([1.0,0.0,0.0], LOS_vec[1] *(pi/180)) # Get a rotation matrix that can rotate any vector θ_u degrees
    LOS_vec = B_vec ./ norm(B_vec)
    LOS_vec = R*LOS_vec # Rotate θ_u degrees around the R_axis to get los_vec with θ_u degrees relative to B-field

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

## ---------------------------------------------------------------------------------------------
# Check how parallel the line-of-sight vector is to one of the cylindrical coordinate directions.
# If it is too parallel, we need special care to successfully map out the LOS voxels
LOS_vec_R = sqrt(LOS_vec[1]*LOS_vec[1]+LOS_vec[2]*LOS_vec[2])
LOS_vec_phi = getPhiFromXandY(LOS_vec[1],LOS_vec[2])
LOS_vec_cyl = [LOS_vec_R, LOS_vec_phi, LOS_vec[3]]
dot_R = dot(LOS_vec_cyl,[1.0,0.0,0.0])/norm(LOS_vec_cyl)
dot_phi = dot(LOS_vec_cyl,[0.0,1.0,0.0])/norm(LOS_vec_cyl)
dot_z = dot(LOS_vec_cyl,[0.0,0.0,1.0])/norm(LOS_vec_cyl)
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

## ---------------------------------------------------------------------------------------------
# Create LOS via 3D array storing solid angles, which represents a uniform (R,phi,z) grid
verbose && println("Creating LOS... ")
# Compute the unit vector in the direction of the detector. Using that vector, go stepwise from the detector along 
# the length of the LOS. Transform points on edge of LOS to (R,phi,z) to find LOS on uniform (R,phi,z) grid (LINE21 format requirement)
# In the end, check which points are inside the tokamak wall/plasma boundary, and keep only those
y_hat = [0.0,1.0,0.0] # The Cartesian y unit vector. Arbitrary direction for Gram-Schmidt process.
PERP2LOS_vec = y_hat - (dot(y_hat,LOS_vec)/dot(LOS_vec,LOS_vec))*LOS_vec # Gram-Schmidt process. To acquire vector perpendicular to line-of-sight
perp2los_vec = PERP2LOS_vec ./norm(PERP2LOS_vec) # Get the unit vector perpendicular to the line-of-sight
EDGEofLOS_vec = (LOS_width/2)*perp2los_vec # Re-scale it to acquire vector that will show edge of LOS, relative to center (spine) of LOS
los_vec = LOS_vec ./norm(LOS_vec) # Get the normalized line-of-sight vector
if LOS_circ_res<=1
    @warn "'LOS_circ_res' set to 1 or less. LOS might be badly resolved."
    PERP2LOS_vecs = [[0.0,0.0,0.0]]
else
    LOS_circ_angles = collect(range(0.0,stop=2*pi,length=Int64(LOS_circ_res))) # Angles to rotate around spine of LOS
    LOS_circ_angles = LOS_circ_angles[1:end-1] # Remove redundant 2*pi angle
    EDGEofLOS_vecs = Vector{Vector{Float64}}(undef,LOS_circ_res-1)
    for (ia,angle) in enumerate(LOS_circ_angles)
        R = getRotationMatrix(los_vec,angle)
        EDGEofLOS_vecs[ia] = R*EDGEofLOS_vec # Get the vectors pointing to the edge of the cylindrical LOS, relative to the center (spine) of the LOS
    end
end

s = collect(range(0.0,stop=LOS_length,length=Int64(LOS_length_res))) # Create increments for walking the normalized LOS vector along LOS
ds = s[2] # Since s[1] is zero
for is=2:Int64(LOS_length_res) # Skip first element, since that would result in an LOS element behind the detector
    verbose && println("Creating diagnostic line-of-sight part $(is-1)/$(Int64(LOS_length_res)-1)... ")
    ss = s[is]
    LOS_spine_point = detector_location - (ss-ds)*los_vec # LOS_vec is pointing towards detector. To go away from detector, use -, not +.
    rminush = norm(s[is-1]*los_vec) # The distance from the detector to the center of the 'lid' of the voxel, where 'lid' is the voxel surface closest to the detector
    r = norm(detector_location - (LOS_spine_point+EDGEofLOS_vecs[1])) # The distance from the detector to the edge of the 'lid' of the voxel. Don't care which EDGEofLOS_vec. Just use the first.
    OMEGA_LOS_crs = 2*pi*(1-rminush/r) # See https://en.wikipedia.org/wiki/Solid_angle#Solid_angles_for_common_objects and pertaining figure for explanation

    if LOS_parallel_to_R
        # If the LOS is parallel to the R axis, we need to map all R points on the uniform (R,phi,z) grid via a different mechanism than the edge-cylinder approach
        LOS_spine_point_R = sqrt(LOS_spine_point[1]*LOS_spine_point[1]+LOS_spine_point[2]*LOS_spine_point[2])
        LOS_edge_Rmax = max(prev_R,LOS_spine_point_R) # To find boundaries of LOS on uniform (R,phi,z) grid (required to adhere to LINE21 format)
        LOS_edge_Rmin = min(prev_R,LOS_spine_point_R) # To find boundaries of LOS on uniform (R,phi,z) grid (required to adhere to LINE21 format)
        prev_R = LOS_spine_point_R
    else
        LOS_edge_Rmax = -Inf # To find boundaries of LOS on uniform (R,phi,z) grid (required to adhere to LINE21 format)
        LOS_edge_Rmin = Inf # To find boundaries of LOS on uniform (R,phi,z) grid (required to adhere to LINE21 format)
    end
    if LOS_parallel_to_phi
        # If the LOS is parallel to the phi axis, we need to map all R points on the uniform (R,phi,z) grid via a different mechanism than the edge-cylinder approach
        LOS_spine_point_phi = getPhiFromXandY(LOS_spine_point[1],LOS_spine_point[2])
        LOS_edge_phimax = max(prev_phi,LOS_spine_point_phi) # To find boundaries of LOS on uniform (R,phi,z) grid (required to adhere to LINE21 format)
        LOS_edge_phimin = min(prev_phi,LOS_spine_point_phi) # To find boundaries of LOS on uniform (R,phi,z) grid (required to adhere to LINE21 format)
        prev_phi = LOS_spine_point_phi
    else
        LOS_edge_phimax = -Inf # To find boundaries of LOS on uniform (R,phi,z) grid (required to adhere to LINE21 format)
        LOS_edge_phimin = Inf # To find boundaries of LOS on uniform (R,phi,z) grid (required to adhere to LINE21 format)
    end
    if LOS_parallel_to_z
        # If the LOS is parallel to the z axis, we need to map all R points on the uniform (R,phi,z) grid via a different mechanism than the edge-cylinder approach
        LOS_spine_point_z = LOS_spine_point[3]
        LOS_edge_zmax = max(prev_z,LOS_spine_point_z) # To find boundaries of LOS on uniform (R,phi,z) grid (required to adhere to LINE21 format)
        LOS_edge_zmin = min(prev_z,LOS_spine_point_z) # To find boundaries of LOS on uniform (R,phi,z) grid (required to adhere to LINE21 format)
        prev_z = LOS_spine_point_z
    else
        LOS_edge_zmax = -Inf # To find boundaries of LOS on uniform (R,phi,z) grid (required to adhere to LINE21 format)
        LOS_edge_zmin = Inf # To find boundaries of LOS on uniform (R,phi,z) grid (required to adhere to LINE21 format)
    end
    for edgeoflos_vec in EDGEofLOS_vecs
        LOS_edge_point = LOS_spine_point + edgeoflos_vec

        R_edge = sqrt(LOS_edge_point[1]*LOS_edge_point[1]+LOS_edge_point[2]*LOS_edge_point[2])
        phi_edge = getPhiFromXandY(LOS_edge_point[1],LOS_edge_point[2])
        if R_edge>LOS_edge_Rmax 
            LOS_edge_Rmax = R_edge
        end
        if R_edge<LOS_edge_Rmin
            LOS_edge_Rmin = R_edge
        end
        if phi_edge>LOS_edge_phimax
            LOS_edge_phimax = phi_edge
        end
        if phi_edge<LOS_edge_phimin
            LOS_edge_phimin = phi_edge
        end
        if LOS_edge_point[3]>LOS_edge_zmax
            LOS_edge_zmax = LOS_edge_point[3]
        end
        if LOS_edge_point[3]<LOS_edge_zmin
            LOS_edge_zmin = LOS_edge_point[3]
        end
    end
    
    edge_Rinds = findall(x-> x>=LOS_edge_Rmin && x<=LOS_edge_Rmax,R_grid)
    edge_phiinds = findall(x-> x>=LOS_edge_phimin && x<=LOS_edge_phimax,phi_grid)
    edge_zinds = findall(x-> x>=LOS_edge_zmin && x<=LOS_edge_zmax,z_grid)
    OMEGA_per_R_phi_z_voxel = OMEGA_LOS_crs/(length(edge_Rinds)*length(edge_phiinds)*length(edge_zinds)) # The solid angle of a single (R,phi,z) voxel is approximately the total solid angle for the whole circular cross section of the cylinder LOS, divided by the number of (R,phi,z) points within that circular cross section. This is an approximation, but hopefully good enough, if the resolution of R_grid, phi_grid and z_grid is high enough.
    omega_3D_array[edge_Rinds,edge_phiinds,edge_zinds] .= OMEGA_per_R_phi_z_voxel
end

## ---------------------------------------------------------------------------------------------
# Use non-zero elements of omega_3D_array to map out LOS on uniform (R,phi,z) grid
# Compute x, y, V, UX, UY, and UZ for those points
# Together with R, phi and OMEGA (and C, redundant), put in data structures to be saved in LINE21-like file
# Filter out points that are inside the tokamak wall/plasma boundary
verbose && println("Mapping out LOS data on uniform (R,phi,z) grid... ")

# Find all non-zero omega_3D_array indices
nz_coords = findall(x-> x>0.0, omega_3D_array)

# Instantiate arrays for x, y, z, ... etc
dR = abs(R_grid[2]-R_grid[1]) # Need to be uniform grid because of LINE21
dphi = abs(phi_grid[2]-phi_grid[1]) # -||-
dz = abs(z_grid[2]-z_grid[1]) # -||-
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
for (inz, nz_coord) in enumerate(nz_coords)
    iR = nz_coord[1]
    iphi = nz_coord[2]
    iz = nz_coord[3]
    R = R_grid[iR]
    phi = phi_grid[iphi]
    z = z_grid[iz]

    x = R * cos(phi)
    y = R * sin(phi)
    V = R*dphi*dR*dz # Volume of cylindrical voxel
    U = detector_location - [x,y,z]
    U = U ./norm(U) # Unit vector pointing from voxel to detector
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
    LOS_in[inz] = in_boundary(wall,R,z) # Is the point inside the tokamak wall/boundary? If so, set bool element to true
end

## ---------------------------------------------------------------------------------------------
# Remove voxels outside of the tokamak wall/plasma boundary. Round to 9 significant digits (like the LINE21 code)
verbose && println("Filtering out points outside of the tokamak wall/plasma boundary... ")
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

## ---------------------------------------------------------------------------------------------
# Save the computed custom LOS data as a .vc file (same format as the output of the LINE21 code)
verbose && println("Saving results... ")
myfile = open(folderpath_o*LOS_name*".vc", "w")
for i=1:length(LOS_x)
    write(myfile,"   $(LOS_x[i]) $(LOS_y[i]) $(LOS_z[i]) $(LOS_C[i]) $(LOS_V[i]) $(LOS_UX[i]) $(LOS_UY[i]) $(LOS_UZ[i]) $(LOS_R[i]) $(LOS_phi[i]) $(LOS_OMEGA[i])\n")
end
close(myfile)

## ---------------------------------------------------------------------------------------------
verbose && println("Your output file can be found at "*folderpath_o*LOS_name*".vc")
verbose && println("~~~createCustomLOS.jl finished successfully!~~~")