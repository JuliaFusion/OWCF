# This simple tutorial will show you how to compute orbits using the functions and tools provided by the OWCF.
# We will also take a look at how to visualize them.
# It is recommended that you execute this script block-by-block in VSCode.
# You do this by highlighting the code of interest, and press ctrl + enter (you can also use shift + enter, to advance line-by-line)
# Each block is separated by a line such as the one below.
# The folderpath_OWCF string you have to fill in yourself. Remember to finish with "/"
##############################################################################################################
folderpath_OWCF = "/folderpath/to/your/OWCF/" # First, you have to specify the path to the OWCF folder, as a string. Change this to your own OWCF path
cd(folderpath_OWCF) # Then, you make the VSCode terminal change place to the OWCF folder
using Pkg # Then, you load the Pkg.jl package which is Julia's package manager
Pkg.activate(".") # Then, you activate the virtual environment of the OWCF
##############################################################################################################
using Equilibrium # Then, we load the Equilibrium.jl package, which we use to handle magnetic equilibria
using GuidingCenterOrbits # Then, we load the GuidingCenterOrbits.jl package, which we use to compute fast-ion drift orbits
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
my_plot = Plots.plot(o.path.r, o.path.z, label="Orbit (equations-of-motion integration)") # Plot the poloidal projection of the orbit trajectory (its R and z coordinates)
my_plot = Plots.plot!(wall.r, wall.z, label="Tokamak wall", color=:black) # Plot the poloidal projection of the tokamak wall. Notice we use plot!() and not plot(). The exclamation mark ! adds more plotting on top of the already existing figure
my_plot = Plots.scatter!([M.axis[1]],[M.axis[2]], markercolor=:red, label="Magnetic axis") # We can also include the magnetic axis as a scatter plot (with a single point, hence why we need to put M.axis[1] and M.axis[2] within [])
my_plot = Plots.plot!(xlabel="R [m]", ylabel="z [m]", aspect_ratio=:equal) # Added information and formatting of the plot
##############################################################################################################
# If we load the extra/gui.jl script, we can investigate an animation of the guiding-center orbit
include(folderpath_OWCF*"extra/gui.jl")
plot_orbit_movie(M, myEPRc; wall=wall)
##############################################################################################################
# Then, you can also use the OWCF to compute full-orbit trajectories (where the gyro-motion is included)
# This is done as follows:
# We start by converting our (E,pm,Rm)-coordinate to a guiding-center particle (see below).
# We have to include the mass (in kg) and the charge (in elemental charge units) as arguments (the two last input arguments).
# This is so the algorithm knows what species the guiding-centre particle is.
# Please note, that even though we create a guiding-center particle object, we are still going to 
# input it to the get_full_orbit() method, which will return a full-orbit trajectory object
my_gc_particle = GCParticle(E,pm,Rm,myEPRc.z, getSpeciesMass(FI_species), getSpeciesEcu(FI_species))
cyclotron_period_gc_particle = cyclotron_period(M,my_gc_particle) # We can compute the cyclotron period for our guiding-center particle
# Please note, that even though it's a guiding-center particle, the cyclotron period is the same as for an actual particle 
# It's just a way to easily compute the cyclotron period.
##############################################################################################################
# Below, we use the get_full_orbit() method to compute the full-orbit trajectory of our (guiding-center) particle.
# We can specify the fineness of the integration via dt, which is the fraction of cyclotron periods.
# We can specify the total length of our integration via tmaxm, which is also given in terms of cyclotron periods (1000 thus mean a thousand cyclotron periods)
full_orbit_path, status = get_full_orbit(M, my_gc_particle ; dt=0.01, tmax=1000)
##############################################################################################################
# We can plot the full-orbit trajectory via the Plots.plot() function. 
# We can specify a variable N, to decide how many (R,z) points we want to plot.
N= 50000 # Specify a number N, to not plot all of the full-orbit path
# We add the full-orbit path to the equation-of-motion plot (my_plot)
Plots.plot!(my_plot, full_orbit_path.r[1:N],full_orbit_path.z[1:N],label="Orbit (full-orbit)")
##############################################################################################################
# Finally, you could use the OWCF to compute guiding-center orbits using the μ-contouring method also
# μ is the magnetic moment
# This can be done as follows
using Contour # First, load the Contour.jl package, because we need that one to make contours
using LinearAlgebra # We need the LinearAlgebra.jl package to be able to take the norm of the B-field vector
# We need to get the toroidal canonical angular momentum Pϕ and the magnetic moment μ for our guiding-center 
# orbit that we computed above. Instead of computing it manually, we can use the HamiltonianCoordinate() method 
# to easily convert our (E,pm,Rm)-coordinate object to a (E,μ,Pϕ)-coordinate object
myHc = HamiltonianCoordinate(M, myEPRc)
# The energy of myHc will be equal to the energy of myEPRc as long as we have a negligable electric field in the 
# tokamak.
# Now, we can check the value of Pϕ and μ as follows
Pϕ = myHc.p_phi
μ = myHc.mu

# Then, define the μ_func function, to be able to compute values of μ from the particle energy E,
# the magnetic field strength B, the toroidal canonical angular momentum Pϕ, the magnetic flux Ψ
# and the product between the major radius R and the toroidal magnetic field Bϕ.
# This way of computing the magnetic moment μ can be achieved by splitting the total particle energy into 
# its perpendicular and parallel parts:
# 
# E = E_⟂ + E_∥ = m*v_⟂^2 / 2 + m*v_∥^2 / 2
#
# and then using the formula for μ = m*v_⟂^2 / 2*B together with the formula for the toroidal canonical angular momentum 
# Pϕ = m*R*v_∥*B_ϕ / B + q * Ψ to isolate μ.
# Programming-wise, the way of defining functions in Julia as is done below is called 'implicit definition'.
# It's essentially a one-liner, different from the 'normal' way of defining a function.
m_kg = myHc.m # We need the mass in kilograms to be able to define the μ function 
q_Coloumb = (myHc.q)*(GuidingCenterOrbits.e0) # We need the charge in Coloumb to be able to define the μ function
μ_func(E, B, Pϕ, Ψ, RBϕ) = begin
    res = E/B - (B/(2*m_kg)) * ((Pϕ-q_Coloumb*Ψ)/RBϕ)^2
    (res > 0.0) ? res : 0.0
end
# Now, we create everything we need for the μ_func function
# The (R,z)-grid resolution of M is usually not fine enough to be able to accurately 
# compute contours of the μ function. So we create new (R,z) vectors that are more fine (at least 100x100)
R_array_fine = collect(range(minimum(M.r),stop=maximum(M.r),length=100)); nR = length(R_array_fine) # A ';' in Julia has the same function as a new line
z_array_fine = collect(range(minimum(M.z),stop=maximum(M.z),length=101)); nz = length(z_array_fine) # Never use the same grid resolution for R and z. That can lead to confusion
E_matrix = E .* ones(nR,nz)
B_matrix = zeros(nR,nz) # Pre-allocate
Ψ_matrix = zeros(nR,nz) # Pre-allocate
RBϕ_matrix = zeros(nR,nz) # Pre-allocate
for (iR,R) in enumerate(R_array_fine), (iz,z) in enumerate(z_array_fine) # enumerate() returns all the indices and values of an array. See ?enumerate()
    B_at_Rz = Bfield(M,R,z)
    B_matrix[iR,iz] = norm(B_at_Rz)
    RBϕ_matrix[iR,iz] = R*B_at_Rz[2] # The second element in the B-vector is Bϕ
    Ψ_matrix[iR,iz] = M(R,z) # Use the flux given by the M file. This could also by an S
end
Pϕ_matrix = Pϕ .* ones(size(B_matrix))
E_matrix = (E*1000*(GuidingCenterOrbits.e0)) .* ones(size(B_matrix)) # 1000*(GuidingCenterOrbits.e0) is to convert from keV to Joule
# The map() function in Julia takes any number of inputs and feeds them into the function that is the first input 
# argument to the map function(). In this case, μ_func is the first input argument, and E_matrix, B_matrix, Pϕ_matrix,
# Ψ_matrix and RBϕ_matrix are what will be fed into μ_func, via map(). The map() function then takes all the elements 
# with the same index in E_matrix, B_matrix, Pϕ_matrix, Ψ_matrix, RBϕ_matrix and feeds them into μ_func(). The output 
# of the map() function will then be a matrix with the same size as E_matrix, B_matrix, Pϕ_matrix, Ψ_matrix and RBϕ_matrix
# (they all have the same size. If they didn't, map would return an error). In short, map does this:
#
# μ_matrix_ij = μ_func(E_ij, B_ij, Pϕ_ij, Ψ_ij, RBϕ_ij)
#
# For every row element i and column element j in E_matrix, B_matrix, Pϕ_matrix, Ψ_matrix, RBϕ_matrix.
μ_matrix = map(μ_func, E_matrix, B_matrix, Pϕ_matrix, Ψ_matrix, RBϕ_matrix)
##############################################################################################################
# Now we can use the Contour.contour() function to pick out the μ-contour of interest 
# from our μ_matrix
cl = Contour.contour(R_array_fine,z_array_fine,μ_matrix,μ) # Find the contours of the μ_matrix with the value μ
l = Contour.lines(cl) # Extract the lines from those contours
indices_of_good_lines = [] # Some lines in l might be bad (not closed contours). We only want closed contours.
for li=1:length(l) # Go through all lines in l
    ll = l[li] # Check line with index li in the array of lines (l)
    Rll, zll = Contour.coordinates(ll) # Get the (R,z) coordinates of that line
    closed_curve_measure = sqrt((Rll[end]-Rll[1])^2  +(zll[end]-zll[1])^2) # Check the distance between the first and last (R,z) coordinates. If it is a closed contour, they will be the same
    closed_curve_ref = sqrt((Rll[end]-Rll[end-1])^2  +(zll[end]-zll[end-1])^2) # Use the distance between the last and second-to-last (R,z) coordinates as reference (they will also be the same, since the Contour.contour() method returns the same (R,z) point for the first, last and second-to-last (R,z) points for closed curves)
    if isapprox(closed_curve_measure, closed_curve_ref, atol=1e-1) # If the distances are approximately the same (should be close to zero)
        push!(indices_of_good_lines,li) # Label this as a good line (closed contour)
    end
end

# Extract and plot the (R,z) points of the closed μ contours
# We can have several closed μ contours (some (E,μ,Pϕ) triplets correspond to two orbits)
# Add them to the plt_EOM plot
for gi in indices_of_good_lines
    Rl_points, zl_points = Contour.coordinates(l[gi])
    my_plot = Plots.plot!(my_plot, Rl_points,zl_points,label="Orbit (mu-contouring)") # Cannot plot in a good way from inside a for-loop. So save the plot as a plot object my_plot
end
display(my_plot) # Display the plot outside of the for loop

# We can also zoom in on the area of interest (the orbit), and remove the legend
# To be able to view only the orbit
# If you change the orbit parameters (E,pm,Rm) and (E,μ,Pϕ), you have to change the 
# xlims and ylims below accordingly, to be able to see the orbit
Plots.plot!(my_plot,xlims=(3.0,3.6),ylims=(-0.3,0.7),legend=false)
# We can observe how the equations-of-motion integration, the full-orbit integration 
# and the μ-contouring result in almost the exact same (R,z) trajectory.
# At least for these default tutorial parameters.

# That concludes this tutorial.