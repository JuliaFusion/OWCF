# This simple tutorial will show you how to compute and plot the expected neutron spectrum from a neutral beam 
# injection (NBI) beam (modelled as a delta function in (E,p) energy pitch space) using the functions and tools 
# provided by the OWCF. The tutorial will show you how to compute the neutron spectrum both for a specific diagnostic 
# line-of-sight (LOS) and for 4*pi spherical emission (all neutrons in all directions).
#
# We will also take a look at how to visualize some of the input and output quantities
#.
# It is recommended that you execute this script block-by-block in VSCode.
# You do this by highlighting the code of interest, and press ctrl + enter (you can also use shift + enter, to 
# advance line-by-line). Each block is separated by a line such as the one below.
# The folderpath_OWCF string you have to fill in yourself. Remember to finish with "/"
#
# Finally, the default variable input values are set so as to model the expected neutron emission spectrum 
# of neutrons born on-axis from NBI injection in JET.
#
# Written by Henrik Järleblad. Last maintained 2025-04-03.
##############################################################################################################

folderpath_OWCF = "/path/to/the/OWCF/" # First, you have to specify the path to where the OWCF folder is on your computer, as a string. 
verbose = true # If you set this to true, the tutorial will print all the print statements :D

# Change the 'folderpath_OWCF' variable to your own OWCF path
cd(folderpath_OWCF) # Then, you make the VSCode terminal change place to the OWCF folder
using Pkg # Then, you load the Pkg.jl package which is Julia's package manager
Pkg.activate(".") # Then, you activate the virtual environment of the OWCF

##############################################################################################################
using Equilibrium # Then, we load the Equilibrium.jl package, which we use to handle magnetic equilibria
using Plots # We load the Plots.jl package, to be able to plot things
using PyCall # We load the PyCall.jl package, to be able to interface with the DRESS code (J. Eriksson et al, Comp. Phys. Comm., 2016)
include(folderpath_OWCF*"misc/availReacts.jl") # We load the OWCF fusion reaction functions to be able to
# easily work with different fusion reactions
pushfirst!(PyVector(pyimport("sys")."path"), "") # Add all DRESS code scripts, to be able to import them

##############################################################################################################
# If you want to load the magnetic field from an .eqdsk file (or from an output file of the OWCF/extra/createSolovev.jl
# script), you need to specify the path to the output file. This is done with the 'filepath_equil' variable 
# below. If you instead want to specify the magnetic field yourself, as a simple 3-element Vector, 
# simply leave the 'filepath_equil' variable unspecified (""), and change the '0.0' elements in the B_vec 
# variable to your liking below.
filepath_equil = ""
B_vec = [0.0, -2.7, 0.0] # Teslas

##############################################################################################################
# In this section, if you have specified the 'filepath_equil' variable, the magnetic equilibrium 
# is loaded and plotted.
if !(filepath_equil=="")
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
            timepoint = "00,0000" # Unknown timepoint for magnetic equilibrium
        end
    end

    verbose && println("Computing flux function on 100x100 (R,z)-grid... ")
    flux_r = range(extrema(wall.r)...,length=100)
    flux_z = range(extrema(wall.z)...,length=100)
    inds = CartesianIndices((length(flux_r),length(flux_z)))
    psi_rz = [M(flux_r[ind[1]], flux_z[ind[2]]) for ind in inds]
    psi_mag, psi_bdry = psi_limits(M)

    wall_dR = maximum(wall.r)-minimum(wall.r)
    plt_crs = Plots.plot()
    plt_crs = Plots.contour!(flux_r,flux_z,psi_rz',levels=collect(range(psi_mag,stop=psi_bdry,length=5)),color=:gray, linestyle=:dot,linewidth=1.5, label="",colorbar=false)
    plt_crs = Plots.plot!(wall.r,wall.z, label="Tokamak wall", color=:black, linewidth=1.5,xaxis=[minimum(wall.r)-wall_dR/10,maximum(wall.r)+wall_dR])
    plt_crs = Plots.scatter!([magnetic_axis(M)[1]],[magnetic_axis(M)[2]],label="Mag. axis",markershape=:xcross,mc=:red)
    plt_crs = Plots.plot!(aspect_ratio=:equal,xlabel="R [m]", ylabel="z [m]")
    display(plt_crs)
end

##############################################################################################################
# If you want to compute the expected neutron spectrum for 4*pi spherical emission, you can skip this part of
# the tutorial. If you want to compute the expected neutron spectrum for a specific diagnostic LOS, you have to 
# specify the path to a LINE21 code (S. Connory et al, private correspondance w. J. Eriksson) output file.
# OR, the path to an output file of the OWCF/extra/createCustomLOS.jl script. This is done with the 
# 'diagnostic_filepath' variable below.
diagnostic_filepath = ""
# Valid diagnostic LOS files can be found in the OWCF/vc_data/ folder

##############################################################################################################
# In this section, the Python packages are loaded. To be able to use the DRESS code interface to compute 
# expected neutron spectra
py"""
import h5py
import numpy as np

import forward
import spec
import transp_dists
import transp_output
import vcone
"""

##############################################################################################################
# In this section, the fusion reaction of interest is specified. Please see the OWCF/misc/availReacts.jl for a 
# list of available fusion reactions, and how they should be specified. The 'b' reactant is the beam (fast)
# particle species.
reaction = "T(D,n)4He"

##############################################################################################################
# In this section, the fusion reaction is checked so that it's available to use in the OWCF 
# Also, reactant data is deduced from the 'reaction' input
if !reactionIsAvailable(reaction)
    error("Fusion reaction $(reaction) is not yet available in the OWCF. The following reactions are available: $(OWCF_AVAILABLE_reactionS). For projected-velocity computations, the following particle species are available: $(OWCF_SPECIES). Please correct and re-try.")
end

thermal_species, FI_species = getFusionReactants(reaction) # Check the specified fusion reaction, and extract thermal and fast-ion species

##############################################################################################################
# In this section, the neutron spectrum calculator is initialized with the DRESS code as a numerical model
py"""
# The '$' in front of many Python variables means that the variable is defined in Julia, not in Python.
reaction = $reaction
test_thermal_particle = spec.Particle($thermal_species) # Check so that thermal species is available in DRESS code
thermal_species = $thermal_species

forwardmodel = forward.Forward($diagnostic_filepath) # Pre-initialize the forward model
"""

##############################################################################################################
# In this section, you specify an (R,z) point of origin for the neutron emission. This needs to be specified
# since (below) you will need to specify the pitch of the NBI beam (see below). The pitch of the NBI beam changes
# depending on where along the beam you are looking. If you would like to compute the expected total neutron 
# spectrum for all points along the NBI beam, please contact henrikj@dtu.dk.
R = 3.0 # meters
z = 0.3 # meters

##############################################################################################################
# In this section, you specify the energy and the pitch of the NBI beam that you want to simulate
energy = 110.0 # keV
pitch = 0.75 # - The cosine of the angle between the NBI beam and the magnetic field at your (R,z) point of
# origin for the neutrons

##############################################################################################################
# In this section, you specify a thermal species temperature and density distribution.
# In the full OWCF, this can actually be done in several different ways (via TRANSP output files, NUBEAM files,
# default profiles, intrapolation, extrapolation, flat profiles etc). However, for this simple tutorial, you 
# just specify a thermal species temperature value and density value below.
thermal_temperature = 5.0 # keV
thermal_density = 3.0e20 # m^-3

##############################################################################################################
# In this section, you specify the lower and upper bounds for your neutron spectrum energies.
# Also, you specify the width of the measurement bins. For example, if you simulate TOFOR neutron spectra,
# you should set the width to 50 keV, since that is the approximate neutron energy resolution of TOFOR.
spectrum_lower_bound = 12500.0 # keV
spectrum_upper_bound = 15500.0 # keV
spectrum_bin_width = 50.0 # keV

##############################################################################################################
# In this section, you specify
#   - How many points 'n_gyro' the beam (fast) particle gyro motion should be discretized into?
#     Please note, this is only discretization in velocity-space to obtain the up-/down-shift of the nominal
#     neutron birth energy. It is NOT a discretization in position-space (R,z). For that, please use...
#   - ...the 'include_FLR_effects' variable. If set to true, FLR (spatial) effects will be included when 
#        computing the neutron spectrum. 
#   - The plasma rotation vector 'plasma_rotation_vector', if plasma rotation should be included. If not, 
#     leave as [0.0,0.0,0.0]
n_gyro = 1000000
include_FLR_effects = false
plasma_rotation_vector = [0.0,0.0,0.0] # m/s

##############################################################################################################
# In this section, the magnetic field (B) is computed given our (R,z) point of neutron emission origin, if 
# the 'B_vec' variable was left unspecified (all elements are zero) (see above)
if iszero(B_vec)
    B_vec = Equilibrium.Bfield(M,R,z)
end

##############################################################################################################
# In this section, some variables are reshaped, to correctly match format expected by the numpy package, in Python
B_vec = collect(reshape(B_vec,3,1))
plasma_rotation_vector = collect(reshape(plasma_rotation_vector,3,1))

##############################################################################################################
# In this section, the expected neutron spectrum is computed, using the Python-Julia interface and calling the 
# DRESS code model object (see above)
py"""
# Define the diagnostic energy bins to which samples will be binned
Ed_bins = np.arange($spectrum_lower_bound,$spectrum_upper_bound,$spectrum_bin_width) # diagnostic spectrum bin edges (e.g. keV)
Ed_vals = 0.5*(Ed_bins[1:] + Ed_bins[:-1])  # bin centers (e.g. keV)
neutron_spectrum = forwardmodel.calc($[energy], $[pitch], $[R], $[z], $[1.0], "", Ed_bins, $B_vec, n_repeat=$n_gyro, reaction=$reaction, bulk_temp=$[thermal_temperature], bulk_dens=$[thermal_density], flr=$include_FLR_effects, v_rot=$plasma_rotation_vector)
"""
neutron_spectrum = Vector(py"neutron_spectrum") # Convert from Python to Julia and vectorize (just to be sure)
neutron_energies = Vector(py"Ed_vals") # -||-

##############################################################################################################
# In this section, the computed expected neutron spectrum is plotted
plt_results = Plots.scatter(neutron_energies ./(1.0e3), (1.0e3) .*neutron_spectrum, label="Expected signal", linewidth=2)
plt_results = Plots.plot!(xlabel="Neutron energy [MeV]",ylabel="Neutrons [per MeV]",legend=:topleft, dpi=600)
display(plt_results)

# Save figure as png file 
png(plt_results,"test.png")