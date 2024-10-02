# This cript will compute thermal temperature and density profiles, given inputs.
# The temperature and density profiles will be computed as functions of rho_pol (normalized poloidal flux).
# This script will create an output file in .jld2 format that can be used as input file 
# to e.g. the calcSpec.jl script, to compute forward model signals in the OWCF

################################################################################
############################ Specify inputs ####################################
################################################################################
dens_profile_type = :flat # Density profile type. Only :flat works for now
dens_value_on_axis = 0.0 # Density profile value on-axis. m^-3

folderpath_OWCF = "" # /path/to/the/directory/of/the/OWCF/
filepath_out = "" # /path/to/where/you/want/your/output/file/to/be/saved.jld2

rho_grid_resolution = 100 # The number of grid points for 

temp_profile_type = :flat # Temperature profile type. Only :flat works for now
temp_value_on_axis = 0.0 # Temperature value on-axis. keV

verbose = false # If true, the script will talk a lot!
################################################################################
############################ Change to OWCF directory ##########################
############################ Load Julia packages ###############################
################################################################################
cd(folderpath_OWCF)
using Pkg
Pkg.activate(".")
using JLD2
################################################################################
############################ Main script #######################################
################################################################################
# rho_pol grid points
verbose && println("Defining rho_pol grid points... ")
rho_array = collect(range(0.0,stop=1.0,length=rho_grid_resolution))

# Temp
verbose && println("Creating temperature profile... ")
if temp_profile_type==:flat
    temp_array = temp_value_on_axis .* ones(rho_grid_resolution)
else
    error("Nothing else supported right now!")
end

# Dens
verbose && println("Creating density profile... ")
if dens_profile_type==:flat
    dens_array = dens_value_on_axis .* ones(rho_grid_resolution)
else
    error("Nothing else supported right now!")
end

################################################################################
############################ Saving ############################################
################################################################################
# Save temperature and density profiles to .jld2 file with file keys suitable for e.g. calcSpec.jl
verbose && println("Saving results to file... ")
myfile = jldopen(filepath_out,true,true,false,IOStream)
write(myfile,"rho_pol",rho_array)
write(myfile,"thermal_temp",temp_array)
write(myfile,"thermal_dens",dens_array)
close(myfile)
println("~~~createThermalProfilesForOWCF.jl done!~~~")