# This simple tutorial will show you how to compute orbits using the functions and tools provided by the OWCF.
# We will also take a look at how to visualize them.
# It is recommended that you execute this script block-by-block in VSCode.
# You do this by highlighting the code of interest, and press ctrl + enter (you can also use shift + enter, to advance line-by-line)
# Each block is separated by a line such as the one below.
# The folderpath_OWCF string you have to fill in yourself. Remember to finish with "/"
##############################################################################################################
folderpath_OWCF = "/path/to/the/OWCF/" # First, you have to specify the path to the OWCF folder, as a string. Change this to your own OWCF path
cd(folderpath_OWCF) # Then, you make the VSCode terminal change place to the OWCF folder
using Pkg # Then, you load the Pkg.jl package which is Julia's package manager
Pkg.activate(".") # Then, you activate the virtual environment of the OWCF
##############################################################################################################
using Equilibrium # Then, we load the Equilibrium.jl package, which we use to handle magnetic equilibria
using GuidingCenterOrbits # Then, we load the GuidingCenterOrbits.jl package, which we use to compute guiding-center fast-ion drift orbits
using Plots # Then, we load the Plots.jl package, to be able to plot things
include(folderpath_OWCF*"misc/species_func.jl") # Then, we load OWCF species functions to be able to easily work with different particle species
##############################################################################################################
folderpath_equil = folderpath_OWCF*"equilibrium/JET/g96100/g96100_0-53.0012.eqdsk" # Specify the filepath to a magnetic equilibrium file.
# The magnetic equilibrium file can be either an .eqdsk file or a .jld2 file obtained from the extra/compSolovev.jl script
M, wall = read_geqdsk(folderpath_equil, clockwise_phi=false) # Then, we load the magnetic equilibrium into the M variable and the 
# tokamak wall into the wall variable, via the read_geqdsk() function. We have to specify whether the coordinate system is arranged so as to 
# have the phi-coordinate in the (R,phi,z) coordinate system point in the clockwise or counter-clockwise direction, viewed from above.
# Almost all tokamak coordinate systems have phi pointing counter-clockwise, viewed from above
##############################################################################################################
E = 100.0 # Then, we specify the energy of the fast ion, in keV
pm = 0.29 # Then, we specify the pitch maximum orbit label of the fast ion
Rm = 3.4 # Then, we specify the radium maximum orbit label of the fast ion, in meters
FI_species = "D" # Then, we specify what type of particle species the fast ion is (D, T, alpha, 3he, p, e etc)
myEPRc = EPRCoordinate(M, E, pm, Rm; amu=getSpeciesAmu(FI_species), q=getSpeciesEcu(FI_species)) # Create an (E,pm,Rm)-coordinate object. Specify the magnetic equilibrium, E, pm, Rm and, optionally, the mass (in atomic mass units, amu) and the particle charge (in elemental charge units, ecu)
##############################################################################################################
o = get_orbit(M, myEPRc; wall=wall, toa=true) # Use the get_orbit() function to compute the guiding-center orbit. Specify the wall, so that the algorithm knows when the guiding-center particle is los, and set toa (which stands for 'try only adaptive') to true, since we only want to try adaptive integration of the equations-of-motion here
##############################################################################################################
Plots.plot(o.path.r, o.path.z, label="Orbit") # Plot the poloidal projection of the orbit trajectory (its R and z coordinates)
Plots.plot!(wall.r, wall.z, label="Tokamak wall", color=:black) # Plot the poloidal projection of the tokamak wall. Notice we use plot!() and not plot(). The exclamation mark ! adds more plotting on top of the already existing figure
Plots.scatter!([M.axis[1]],[M.axis[2]], markercolor=:red, label="Magnetic axis") # We can also include the magnetic axis as a scatter plot (with a single point, hence why we need to put M.axis[1] and M.axis[2] within [])
Plots.plot!(xlabel="R [m]", ylabel="z [m]", aspect_ratio=:equal) # Added information and formatting of the plot
##############################################################################################################
# That concludes this basic tutorial on how to compute and plot guiding-center orbits 
# with the tools provided by the OWCF