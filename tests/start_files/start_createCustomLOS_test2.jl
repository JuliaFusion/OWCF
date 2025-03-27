################################ start_createCustomLOS_test1.jl #########################################
# This file contains all the inputs that the script createCustomLOS.jl needs to create a custom diagnostic 
# line-of-sight (LOS) for synthetic diagnostics. This file also executes the script createCustomLOS.jl after the 
# inputs are defined. PLEASE NOTE! THIS FILE IS INCOMPLETE COMPARED TO THE TEMPLATE FILE
# OWCF/templates/start_createCustomLOS_template.jl. THIS IS TO MAKE IT FIT INTO THE TEST FLOW BY THE 
# OWCF/tests/run_tests.jl script.
#
# createCustomLOS.jl allows the user of the OWCF to create a custom line-of-sight (LOS) to use as a 
# diagnostic for computing weight functions. The LOS will be saved as a .vc file that can be used by 
# other OWCF script, often as the 'filepath_diagnostic' input variable. The LOS can be constructed in 
# either of the following 5 ways:
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

#### The inputs are as follows:
# folderpath_OWCF - The file path to the OWCF folder - String
# coordinate_system - A symbol input to tell the script what type of coordinate system the 'LOS_vec'
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
 
### Advanced inputs (should NOT be changed from default values, unless you know what you're doing):
# LOS_circ_res - The number of points by which to resolve the LOS azimuthally (around the cylinder
#                LOS). - Int64
# LOS_length_res - The number of points by which to resolve the LOS along the LOS - Int64
# nR - The number of major radius grid points for the uniform (R,phi,z) grid onto which to map the LOS - Int64
# nz - The number of vertical grid points for the uniform (R,phi,z) grid onto which to map the LOS - Int64
# nphi - The number of toroidal grid points for the uniform (R,phi,z) grid onto which to map the LOS - Int64

#### Other
# 

# Script written by Henrik Järleblad. Last maintained 2025-03-26.
######################################################################################################

## First you have to set the system specifications
using Distributed # Needed to be loaded, even though multi-core computations are not needed for createCustomLOS.jl.

## -----------------------------------------------------------------------------
@everywhere begin
    coordinate_system = :cylindrical # :cartesian, :cylindrical or :magrel
    (coordinate_system==:magrel) && ((R_of_interest, z_of_interest)=(:mag_axis, :mag_axis)) # If you want to specify the LOS as an angle relative to the magnetic field, you have to specify an (R,z) point of interest. Accepted values are Float64 and :mag_axis (magnetic axis). 
    detector_location = [2.97, 0.0, 18.8] # TOFOR (proxy): [2.97, 0.0, 18.8] /// NE213 (proxy): [8.35,2.1,-0.27]. When coordinate_system==:magrel, detector_location can be set to whatever (the script will figure it out automatically from the angle and LOS_length)
    filepath_equil = folderpath_OWCF*"equilibrium/JET/g96100/g96100_0-53.0012.eqdsk" # To load tokamak wall/plasma boundary and magnetic equilibrium (for coordinate_system==:magrel)
    folderpath_o = folderpath_OWCF*"tests/outputs/"
    LOS_length = 22.0 # Length of line-of-sight. From detector to desired end. Meters. When coordinate_system==:magrel, assume (R,z) point of interest is at half the LOS length
    LOS_name = "test_detector_2"
    LOS_vec = [0.005,0.0,0.9999] # TOFOR (proxy) cartesian: [0.005,-0.001,0.9999] /// NE213 (proxy) cartesian: [0.995902688, 0.293844478e-2, -0.903836320e-1]
    LOS_width = 0.27 # TOFOR (proxy): 0.27 /// NE213 (proxy): 0.25
    plot_LOS = true # If set to true, the LOS will be plotted (top view and poloidal projection) after creation. For validation purposes.
    plot_LOS && (save_LOS_plot = true) # If set to true, the LOS plot will be saved in .png file format
    verbose = true
    
    # Advanced inputs
    LOS_circ_res = 100 # The cylindrical shape that is the line-of-sight will be resolved azimuthally into these many points
    LOS_length_res = 20000 # The cylindrical shape that is the line-of-sight will be resolved parallel to LOS_vec into these many points
    nR = 200 # The number of R grid points
    nz = 200 # The number of z grid points
    nphi = 720 # The number of phi grid points
end

## -----------------------------------------------------------------------------
# Then you execute the script
include(folderpath_OWCF*"extra/createCustomLOS.jl")