################################ start_createCustomFIDistr_template.jl #########################################
# This file contains all the inputs that the script OWCF/extra/createCustomFIDistr.jl needs to create a custom 
# fast-ion distribution. This file also executes the script OWCF/extra/createCustomFIDistr.jl after the inputs 
# are defined.
#
# OWCF/extra/createCustomFIDistr.jl allows the user of the OWCF to create a custom fast-ion distribution. The 
# output file can be used as an input file in other OWCF scripts, e.g. as the 'filepath_FI_distr' input variable 
# in the calcSpec.jl script. As of the current OWCF version, the custom FI distribution can be created in the 
# following ways:
#   - Gaussian. By specifying the 'distribution_type' input variable as :gaussian, a custom Gaussian will be 
#               created and saved as the custom fast-ion distribution. Please see specific input variables below.
#   - Slowing-down. By specifying the 'distribution_type' input variable as :collisional, a custom slowing-down 
#                   distribution function will be created and saved as the -|-
#   - Custom. By specifying the 'distribution_type' input variable as :custom, a custom custom fast-ion distribution 
#             will be created and saved as the custom fast-ion distribution. PLEASE NOTE! This requires manual 
#             specifications of code in the OWCF/extra/createCustomFIDistr.jl script, as well as custom input 
#             variables to be defined.
# The fast-ion distribution can also be easily made to be constant in position space, i.e. in (R,z) space. This 
# is done via the 'constant_Rz' input variable. For further information about the input variables, please see below.

#### The inputs are as follows:
# folderpath_OWCF - The file path to the OWCF folder - String
#
# constant_Rz - Set to true, if the FI distribution should be computed as constant in (R,z) position space - Bool
# distribution_type - The distribution type of the to-be-created custom fast-ion distribution. As of the current 
#                     OWCF version, the following options are supported:
#                       - :gaussian - The fast-ion distribution will be created as a Gaussian distribution.
#                       - :collisional - The fast-ion distribution will be created as a slowing-down distribution 
#                                        function. The injection energy and characteristics of the sllowing-down
#                                        distribution function are specified via the specific input variables 
#                                        (please see below)
#                       - :custom - The fast-ion distribution will be created in a custom way. PLEASE NOTE! This
#                                   requires manual specifications of custom input variables, as well as manual 
#                                   addition of code in the OWCF/extra/createCustomFIDistr.jl script.
#                       The 'distribution_type' input variables is of the Symbol type - Symbol
# filename_out - The name of the fast-ion distribution output data .jld2 file. By default, it is unspecified ("") and the default 
#                filename format is used to name the output data .jld2 file (see OWCF/extra/createCustomFIDistr.jl for info).
#                PLEASE NOTE! Do NOT include the .jld2 filename extension in the 'filename_out' input variable - String
# folderpath_out - The path to the folder in which you want the custom fast-ion distribution data to be saved as 
#                  a .jld2 file. Remember to finish folder paths with '/' - String
# max_E - The maximum of the energy grid points of the to-be-created custom fast-ion distribution. In keV - Float64
# max_p - The maximum of the pitch grid points of the -||-. Dimensionless - Float64
# max_R - The maximum of the major radius (R) grid points of the -||-. In meters - Float64
# max_z - The maximum of the vertical (z) grid points of the -||-. In meters - Float64
# min_E - The minimum of the energy grid points of the -||-. In keV - Float64
# min_p - The minimum of the pitch grid points of the -||-. Dimensionless - Float64
# min_R - The minimum of the major radius (R) grid points of the -||-. In meters - Float64
# min_z - The minimum of the vertical (z) grid points of the -||-. In meters - Float64
# nE - The number of grid points in energy (E) - Int64
# np - The number of grid points in pitch (p) - Int64
# nR - The number of grid points in major radius (R). PLEASE NOTE! Even if constant_Rz is set to true, this needs to be specified! - Int64
# nz - The number of grid points in vertical coordinate (z). PLEASE NOTE! Even if constant_Rz is set to true, this needs to be specified! - Int64
# save_plots - If true, (E,p) and (R,z) plots of the computed FI distribution will be saved in .png file format. For post-analysis - Bool
# tot_N_FI - The total number of fast ions in the plasma. Obtained when integrating over all of (E,p,R,z) space - Int64/Float64
# verbose - If true, the script will talk a lot! - Bool
# 
# Then, depending on how the 'distribution_type' input variable is specified, different sets of specific inputs need to be specified.
# These specific inputs are as follows:
# If the 'distribution_type' input variable is set to :gaussian, the following input variables need to be specified:
#   - floor_level - Below floor_level*maximum(f_Gaussian), the FI distribution will be set to 0. For stability reasons - Float64
#   - peak_E - The energy (E) coordinate of the peak of the gaussian. In keV - Float64
#   - peak_p - The pitch (E) coordinate of the peak of the gaussian. Dimensionless - Float64
#   - peak_R - The major radius (R) coordinate of the peak of the gaussian. In meters. If constant_Rz==true, this does not matter - Float64
#   - peak_z - The vertical coordinate (z) coordinate of the peak of the gaussian. In meters. If constant_Rz==true, this does not matter - Float64
#   - Sigma - The values of Sigma (correlation lengths) for the Gaussian distribution, i.e. 
#             (1/(4*pi^2)) * (det(Sigma)^(-1/2)) * exp(-0.5*transpose(X-X_peak)*(Sigma^(-1))*(X-X_peak)) - Array{Float64}
# If the 'distribution_type' input variable is set to :collisional, the following input variables need to be specified:
#   - assume_Ti_equalTo_Te - Should the thermal ion temperature be assumed to be equal to the (thermal) electron temperature (Ti=Te)? If so,
#                            set this input variable to true - Bool
#   - assume_ni_equalTo_ne - Should the thermal ion density be assumed to be equal to the (thermal) electron density (ni=ne)? If so,
#                            set this input variable to true - Bool
#   - dampen - Set this input variable to true, if the slowing-down function should be damped below the critical energy (where ion drag equals 
#             electron drag) - Bool
#   - FI_species - The fast-ion particle species. A valid list of species can be found in OWCF/misc/species_func.jl - String
#   - filepath_equil - The file path to a file containing magnetic equilibrium data. Can be an .eqdsk file or an output file of the 
#                   OWCF/extra/createCustomMagneticEquilibrium.jl script. E.g. "g94701_0-50.7932.eqdsk" or "solovev_equilibrium_2024-09-30.jld2" - String
#   - filepath_thermal_profiles - The file path to a file containing thermal plasma temperature and density data. Can be a TRANSP output file, an output 
#                                 file of the OWCF/extra/createCustomThermalProfiles.jl or "". E.g. "96100J01.cdf", "my_custom_thermal_profiles.jld2" or ""
#                                 - String
#   - inj_E - The injection energy (E) coordinate of the collision-physics function. In keV - Float64
#   - inj_p - The injection pitch (p) coordinate of the collision-physics function. Dimensionless - Float64
#   - inj_R - The injection major radius (R) coordinate of the collision-physics function. In meters - Float64
#   - inj_z - The injection vertical (z) coordinate of the collision-physics function. In meters - Float64
#   - R_of_temp_n_dens - The major radius (R) coordinate of the thermal species temperature+density values. By default, the inj_R input variable is 
#                        used to determine the major radius (R) coordinate of the thermal temperature+density profiles. In meters - Float64
#   - z_of_temp_n_dens - The vertical (z) coordinate of the thermal species temperature+density values. By default, the inj_z input variable is 
#                        used to determine the vertical (z) coordinate of the thermal temperature+density profiles. In meters - Float64
#   - thermal_species - The thermal ion species of interest for the temp+dens profiles. Valid list of species can be found in misc/species_func.jl.
#                       PLEASE NOTE! Currently, only one thermal species is supported (multi-species thermal plasmas are not supported) - String
#   - filepath_thermal_electron_profiles - If you haven't specified the (thermal) electron temperature in any other way, please use this input 
#                                          variable to specify the file path to an output file of the OWCF/extra/createCustomThermalProfiles.jl
#                                          script. The temperature and density profiles data in this file will be used as the electron temperature 
#                                          and density - String
#   If the 'dampen' input variable is set to true, please specify...
#       - damp_type - The type of damping. Can be either :erfc (complementory error-function) or :linear - Symbol
#       - E_tail_length - The length (in keV) of the tail of the collision-physics functions below the injection energy. in keV - Float64 or Nothing
#        /\ If left as nothing, automatic damping below the v_L speed (W.G.F. Core, Nucl. Fusion 33, 829 (1993)) will be used
#   If the 'filepath_thermal_profiles' is specified as a TRANSP output file, please specify...
#       - timepoint - The timepoint of interest for the temperature+density profiles to be extracted from the TRANSP .cdf format file - Float64
#   If the 'filepath_thermal_profiles' is specified as "", please specify...
#       - dens_on_axis - The thermal ion density value on-axis. In m^-3 - Float64
#       - temp_on_axis - The thermal temperature value on-axis. In keV - Float64
# If the 'distribution_type' input variable is set to :custom, the following input variables need to be specified:
#   - custom_FI_distribution(E,p,R,z) - A function that takes an (E,p,R,z) point as input, and returns a distribution value for that point - Function
#   If the 'constant_Rz- input variable is set to true, please specify...
#       - custom_FI_Ep_distribution(E,p) - A function that takes an (E,p) point as input, and return -||- - Function

#### Other
# 

# Script written by Henrik Järleblad. Last maintained 2025-04-24.
################################################################################################################

## First you have to set the system specifications
using Distributed # Needed to be loaded, even though multi-core computations are not needed for createCustomLOS.jl.
folderpath_OWCF = "" # OWCF folder path. Finish with '/'

## Navigate to the OWCF folder and activate the OWCF environment
cd(folderpath_OWCF)
using Pkg
Pkg.activate(".")

## -----------------------------------------------------------------------------
# Input variables
@everywhere begin
    constant_Rz = false # By default, assume that the fast-ion distribution is not constant in (R,z)
    distribution_type = :gaussian # The type of the fast-ion distribution to be created. Currently supported options are :gaussian, :collisional and :custom
    filename_out = "" # The name of the output data .jld2 file. Please DO NOT include the .jlde file extension!
    folderpath_out = "" # The output folder in which to save the custom fast-ion distribution data
    max_E = 0.0 # keV
    max_p = 0.0 # Dimensionless
    max_R = 0.0 # Meters
    max_z = 0.0 # Meters
    min_E = 0.0 # keV
    min_p = 0.0 # Dimensionless
    min_R = 0.0 # Meters
    min_z = 0.0 # Meters
    nE = 0 # Number of energy grid points
    np = 0 # Number of pitch grid points
    nR = 0 # Number of major radius grid points
    nz = 0 # Number of vertical grid points
    save_plots = false
    tot_N_FI = 0 # Total number of fast-ions, i.e. ∫ f(E,p,R,z) 2*pi*R dEdpdRdz
    verbose = false

    if distribution_type==:gaussian
        floor_level = 0.0
        peak_E = 0.0 # keV
        peak_p = 0.0 # Dimensionless
        peak_R = 0.0 # Meters
        peak_z = 0.0 # Meters
        Sigma = [0.0,0.0,0.0,0.0] # In [keV,Dimensionless,m,m]
    elseif distribution_type==:collisional
        assume_Ti_equalTo_Te = false # By default, do not assume that Ti==Te
        assume_ni_equalTo_ne = true # By default, assume that ni==ne
        dampen = false
        FI_species = "D" # Please see OWCF/misc/species_func.jl for a list of available particle species
        filepath_equil = "" # E.g. "g94701_0-50.7932.eqdsk" or "solovev_equilibrium_2024-09-30.jld2"
        filepath_thermal_profiles = "" # E.g. "96100J01.cdf", "my_custom_thermal_profiles.jld2" or ""
        inj_E = 0.0 # keV
        inj_p = 0.0 # Dimensionless
        inj_R = 0.0 # Meters
        inj_z = 0.0 # Meters
        R_of_temp_n_dens = inj_R # By default, assume that the thermal ion temp+dens R point of interest is at the R injection point
        z_of_temp_n_dens = inj_z # By default, assume that the thermal ion temp+dens z point of interest is at the z injection point
        thermal_species = "" # Please see OWCF/misc/species_func.jl for a list of available particle species
        filepath_thermal_electron_profiles = "" # PLEASE NOTE! Only specify this if you have not specified the electron temp+dens in any other way
        if dampen # If you specify the 'dampen' input variable to true, please specify...
            damp_type = :erfc # The type of damping, currently available :erfc and :linear
            E_tail_length = nothing # The length of the energy tail below the damping energy. In keV
        end
        timepoint = 00.000 # If the 'filepath_thermal_profiles' input variable is specified as a TRANSP output file in .cdf file format, please specify a shot timepoint. In seconds
        if filepath_thermal_profiles=="" # If the 'filepath_thermal_profiles' input variable is left unspecified please specify... 
            dens_on_axis = 0.0 # Thermal ion density on-axis. In m^-3
            temp_on_axis = 0.0 # Thermal ion temperature on-axis. In keV
        end
    else # must be :custom
        if constant_Rz
            function custom_FI_Ep_distribution(E,p)
                # Define a function that takes an (E,p) point as input, and returns a distribution value for that point
                return f
            end
        else
            function custom_FI_distribution(E,p,R,z)
                # Define a function that takes an (E,p,R,z) point as input, and returns a distribution value for that point
                return f
            end
        end
    end
end

## -----------------------------------------------------------------------------
# Change directory to OWCF-folder on all external processes. Activate the environment there to ensure correct package versions
# as specified in the Project.toml and Manuscript.toml files.
@everywhere begin
    using Pkg
    cd(folderpath_OWCF)
    Pkg.activate(".")
end

## -----------------------------------------------------------------------------
# Then you execute the script
include(folderpath_OWCF*"extra/createCustomFIDistr.jl")