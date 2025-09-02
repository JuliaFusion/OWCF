################################### weightsWebApp.jl ##########################################

#### Description:
# This script provides an application to visualize orbit weight functions in an interactive and
# intuitive manner. 
#
# It visualizes orbit weight functions as a function of orbit space (E,pm,Rm) coordinates.
# Since the orbit weight functions are three-dimensional, fast-ion energy (E) slices 
# are visualized. The fast-ion energy can be changed via an interactive slider to easily 
# scroll through different fast-ion energies.
#
# The weights (delta function signals) that make up the orbit weight functions can also be 
# visualized. This is done via the 'show_delta_signal' button in the web application. To be clear,
# in addition to having the web application help the user to visualize W[Ed,E,:,:] (where 'W' is 
# the orbit weight function, 'Ed' is a specific diagnostic signal measurement bin center and 'E'
# is a specific fast-ion energy), this 'show_delta_signal'button allows the user to visualize 
# W[:,E,pm,Rm] where 'pm' is a specific pitch maximum and 'Rm' is a specific major radius maximum.
#
# Prior to running this script, please make sure you have run the following scripts:
#   - calcOrbWeights.jl (or equivalent)
#   - extractTopoBounds.jl
#   - orbWeights_2Dto4D.jl (if not automatically in calcOrbWeights.jl)
# And that you have noted the paths to the outputs. You will also (partially) need the following: 
#   - .eqdsk-file containing tokamak magnetic equilibrium and geometry
#   - .hdf5-file/.jld2 with synthetic/experimental signal (please see script further down for specs)
# Note! OWCF has the following diagnostics' sightlines already built-in:
#   - TOFOR (JET neutron diagnostic). Keyword for usage: "TOFOR"
#   - Afterburner (NE213, JET neutron diagnostic). Keyword for usage: "AB"
#   - No sightline. Keyword for usage: ""
# But an arbitrary sightline can be specified as well.
#
# Continuing, the user can also choose to specify the path to a file containing data for a 
# fast-ion orbit-space distribution. The same orbit grid as the weights is required. In that case,
# the app will show a 2D slice (constant fast-ion energy) of the fast-ion distribution, and of the
# non-integrated WF density, along with everything else. Please consult the OWCF for further info.
#
# All energy slices are visualized with Rm on the x-axis and pm on the y-axis by default. This can
# be switched to the toroidal canonical angular momentum and magnetic moment via a toggle button 
# in the app.
#
# Run the script by starting the Julia command prompt, and type 'include("weightsWebApp.jl")'.
# Or 'include("apps/weightsWebApp.jl")' if standing in the OWCF/ folder (which you should).
#
# NOTE! Make sure to set the inputs in the ## ------ # Inputs section first.
# After running the above command, the script will create the web application (might take several
# minutes). When you see "Task runnable...", open a web browser of your choice and type 
# 'localhost:'+ the port number you have set with the port variable.
# The web application will then load (might take 1-2 minutes). Et voilà, have fun!
# (might be a bit slow. Especially on Windows. Future speed improvements are planned.)
#
# ALSO NOTE! It is recommended to ensure that the orbit weights and the topological boundaries
# have exactly the same dimensions and ranges, by using the orbit weights to calculate the 
# topological map in calcTopoMap.jl (set useWeightsFile to true) in the first place.

#### Inputs (units given when defined in the script):
# folderpath_OWCF - The path to the OWCF folder on your computer. - String
# port - The I/O port on which to host the web application (0001-19999) - Int64
# verbose - If true, then the web application will talk a lot! Yay! - Bool
# filepath_tb - The path to the .jld2-file containing the topological boundaries - String
# enable_COM - If set to true, it will be possible to also visualize the orbit weight functions in (E,Λ,Pϕ_n;σ) coordinates - Bool
# filepath_W_COM - If enable_COM is set to true, it is highly recommended that you provide the path to a file containing the weight matrix in (E,Λ,Pϕ_n;σ) format, e.g. from OWCF/helper/os2com.jl - String
# filepath_tm - If enable_COM is set to true, it is highly recommended that you provide the path to a file containing the pertaining topological map - String
# filepath_equil - The path to the .eqdsk-file with the tokamak magnetic equilibrium and geometry - String
# diagnostic_filepath - The path to the LINE21 file containing line-of-sight data for the diagnostic - String
# diagnostic_name - The name of the diagnostic to be visualized. For esthetic purposes. - String
# filepath_W - The path to the .jld2/.h5 weights file, containing orbit weights (4D) to be visualized - String 
# filepath_Fos3D - The path to the 3D orbit-space fast-ion distribution to be visualized as well. Enable with plot_Fos - String
# filepath_S - The path to the .hdf5/.jld2 synthetic signal spectrum for the diagnostic. Enable with plot_S - String
# specFileJLD2 - If true, the app will assume the synthetic signal file is .jld2 format - Bool
# filepath_no - The path to the 3D orbit-space fast-ion null-region boundaries. Enable with showNullOrbs - String
# filepath_WF - The path to the WF-signal to be visualized as well. Enable with plot_WF - String
# FI_species - The particle species for the orbit weight functions - String
#                     "D" for deuterium, "T" for tritium, "p" for proton, "3he" for helium-3 etc
# xlabel - The x-axis label for the WF and/or S signal to be visualized. Please specify as "X [units]", e.g. "Neutron energy [keV]" - String
# ylabel - The y-axis label for the WF and/or S signal to be visualzed. Please specify as "Y [units]", e.g. "Neutron counts [(keV*s)^-1]". If no S or WF is specified, please leave as "" or nothing - String
# verbose - If set to true, the app will talk a lot - Bool

#### Outputs
# -

#### Saved files
# - 

#### Other
# Warning! Please note! For orbit-space grids containing more than approximately 150 000 valid orbits (e.g. 20x100x100),
# you should NOT use weightsWebApp.jl (or any other interactive app). As of OWCF version 1.4, the 
# web interface simply becomes too slow. Please do instead plot the energy slices manually instead (in scratch.jl for example).
# You do this by, for example, coding
#
# folderpath_OWCF = "/the/path/to/the/OWCF/folder/"
# cd(folderpath_OWCF)
# using Pkg
# Pkg.activate(".")
# using JLD2
# using Plots
#
# myfile = jldopen(filepath_W, false, false, false, IOStream)
# W = myfile["Wtot"]
# Ed_array = myfile["Ed_array"]
# E_array = myfile["E_array"]
# pm_array = myfile["pm_array"]
# Rm_array = myfile["Rm_array"]
# close(myfile)
#
# myfile = jldopen(filepath_tb,false,false,false,IOStream)
# topoBounds = myfile["topoBounds"]
# close(myfile)
#
# Ed = 14100.0 # Example of 14.1 MeV
# E = 150.0 # Example of 150 keV
# iEd = argmin(abs.(Ed_array - Ed))
# iE = argmin(abs.(E_array - E)) # Find the closest value to E in E_array
#
# ones_carinds = findall(x-> x==1.0,topoBounds[iE,:,:])
# pm_scatvals_tb = zeros(length(ones_carinds))
# Rm_scatvals_tb = zeros(length(ones_carinds))
# for (ind,carinds) in enumerate(ones_carinds)
#     pm_scatvals_tb[ind] = pm_array[carinds[1]]
#     Rm_scatvals_tb[ind] = Rm_array[carinds[2]]
# end
#
# Plots.heatmap(Rm_array,pm_array,(W[iEd,iE,:,:])./maximum(W[iEd,iE,:,:]),colorbar=true,title="W ($(round(maximum(W[iEd,iE,:,:]),sigdigits=4)) = 1.0)", fillcolor=cgrad([:white, :darkblue, :green, :yellow, :orange, :red]))
# Plots.scatter!(Rm_scatvals_tb,pm_scatvals_tb,markersize=ms,leg=false,markercolor=:black, xlabel="Rm [m]", ylabel="pm")

# Furthermore, the filepath_W_COM input variable can be either an output of this weightsWebApp.jl script, or the os2com.jl script.

# Script written by Henrik Järleblad and Andrea Valentini. Last maintained 2025-09-01.
###############################################################################################

## --------------------------------------------------------------------------
# Please specify the OWCF folder and let the app change directory to the 
# OWCF folder when weightsWebApp.jl is executed. This is to be able to load the
# correct versions of the Julia packages as specified in the Project.toml and 
# Manifest.toml files.
folderpath_OWCF = "" # Finish with '/'
cd(folderpath_OWCF)
using Pkg
Pkg.activate(".")

## ------
# Inputs
port = 9999 # The computer connection port to connect to the web app with
verbose = true # If true, then the app will talk a lot!
filepath_tb = "" # .jld2 file with orbit-space topological boundaries (see extractTopoBounds.jl, and/or calcTopoMap.jl)
enable_COM = false # Set to true if you would like to be able to switch between (E,pm,Rm) and (E,Λ,Pϕ_n;σ). Please note! This will require extra loading time and computer resources!
if enable_COM
    filepath_W_COM = "" # If you would like to be able to switch to (E,Λ,Pϕ_n;σ), you should definitely provide the path to a file containing the orbit weight matrix mapped to COM (E,Λ,Pϕ_n;σ). This will greatly speed-up pre-app computations
    # OR
    filepath_tm = "" # If you have not computed such a file, you should at least provide the path to a topological map file. Otherwise, the orbit grid will need to be computed from scratch
end
filepath_equil = "" # The equilibrium file containing magnetic equilibrium data
diagnostic_filepath = "" # The file path to the LINE21 output file, containing viewing cone data for the diagnostic
diagnostic_name = "" # Diagnostic sightline aestethic keyword. E.g: "TOFOR", "AB" or ""
filepath_W = ""
(plot_Fos = false) && (filepath_Fos3D = "") # Enable and specify an orbit-space fast-ion distribution file (3D-format)
(plot_S = false) && (filepath_S = "") # Enable and specify a signal file to include in the plotting
specFileJLD2 = true # Assume .jld2 file format for signal file by default. Otherwise, assume .hdf5 file format. Only applicable if plot_S==true
(showNullOrbs = false) && (filepath_no = "")
(plot_WF = false) && (filepath_WF = "")
FI_species = "" # Specify with D, T, p, 3he etc
xlabel = "" # Example neutron energy: "Neutron energy [keV]". Example projected velocity: "Projected velocity [m/s]"
ylabel = "" # Example neutron count: "Neutron count [(keV*s)^-1]". Example projected velocity signal: "Projected velocity signal [m^-1]"
# The projected velocity signal units work because if we integrate the signal, we should get the number of counts per second. So [m^-1]*[m/s] = [s^-1]
verbose = true

########################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################
############### WEB APPLICATION BEGINS BELOW. DO NOT ALTER ANY CODE BELOW THIS LINE! ###################
########################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################

# Must load JLD2 package first, to be able to check filepath_tm and filepath_W for 'extra_kw_args'
using JLD2
# EXTRA KEYWORD ARGUMENTS BELOW (these will go into the orbit_grid() and get_orbit() functions from GuidingCenterOrbits.jl)
# Try to find extra_kw_args in filepath_tb or filepath_W
myfile = jldopen(filepath_tb,false,false,false,IOStream)
myfile2 = jldopen(filepath_W,false,false,false,IOStream)
if haskey(myfile,"extra_kw_args")
    extra_kw_args = myfile["extra_kw_args"]
elseif haskey(myfile2,"extra_kw_args")
    extra_kw_args = myfile2["extra_kw_args"]
else
    extra_kw_args = Dict(:toa => true, :limit_phi => true)
    # toa is 'try only adaptive'
    # limits the number of toroidal turns for orbits
end
close(myfile)
close(myfile2)

if haskey(extra_kw_args, :max_tries)
    delete!(extra_kw_args, :max_tries)
end

## ------
# Loading packages
verbose && println("Loading Julia packages... ")
using Distributed
using Interact
using EFIT
using Equilibrium
using GuidingCenterOrbits
using SparseArrays
using HDF5
using Plots
using FileIO
using Mux
using WebIO
using PyCall # For loading sightlines efficiently
include(folderpath_OWCF*"misc/species_func.jl") # To get particle masses and charges
include(folderpath_OWCF*"misc/availReacts.jl") # To check reaction availability and extract fast-ion and thermal species
include(folderpath_OWCF*"extra/dependencies.jl") # To load the (E,pm,Rm) to (E,Λ,Pϕ_n;σ) mapping function

## ------
# Loading packages on external CPU processors, if anticipated that it will be needed later
if enable_COM && !isfile(filepath_W_COM)
    verbose && println("Anticipating the need for several CPU-processors. Adding CPU-processors and loading necessary Julia packages... ")
    addprocs(2) # Two extra CPU processes should be sufficient, and available on most PCs/Macbook/Linux devices, also in terms of RAM
    @everywhere folderpath_OWCF = $folderpath_OWCF # Transfer variable to extra CPU-processes
    @everywhere begin
        cd(folderpath_OWCF)
        using Pkg
        Pkg.activate(".")
        using Equilibrium
        using EFIT
        using GuidingCenterOrbits
        include(folderpath_OWCF*"extra/dependencies.jl") # To load the (E,pm,Rm) to (E,Λ,Pϕ_n;σ) mapping function
    end
end

verbose && println("Loading Python packages... ")
py"""
import numpy as np
"""

## ------
# Extract y-axis units for S and/or WF from xlabel input variable
Ed_units = String(split(split(xlabel,"[")[2],"]")[1])

## ------
# Read the .jld2-file for displaying the topological boundaries of orbit space.
# Pre-calculated with getTopoMap.jl followed by extractTopoBounds.jl.
println("Loading topological boundaries... ")
myfile = jldopen(filepath_tb,false,false,false,IOStream)
topoBounds = myfile["topoBounds"]
close(myfile)

## ------
# Read the .jld2-file for displaying the null-region boundaries of orbit space.
if showNullOrbs
    if !isfile(filepath_no)
        error("Null-measurement boundaries set for plot (showNullOrbs=true). But filepath_no is invalid. Please correct and re-try.")
    end
    verbose && println("Loading null orbits... ")
    myfile = jldopen(filepath_no,false,false,false,IOStream)
    nullOrbs_indices = myfile["nullOrbs_indices"]
    close(myfile)
end


## ------
# Load the 3D orbit-space fast-ion distribution, the slices of which will be visualized
if plot_Fos
    if !isfile(filepath_Fos3D)
        error("Fast-ion orbit-space distribution set for plot (plot_F_os=true). But filepath_Fos3D is invalid. Please correct and re-try.")
    end
    println("Loading 3D orbit-space fast-ion distribution... ")
    myfile = jldopen(filepath_Fos3D,false,false,false,IOStream)
    if haskey(myfile,"F_os_3D")
        F_os_3D = myfile["F_os_3D"]
        F_os_4D = reshape(F_os_3D,(1,size(F_os_3D,1),size(F_os_3D,2),size(F_os_3D,3)))
    elseif haskey(myfile,"F_reconstructions_3D")
        F_os_4D = myfile["F_reconstructions_3D"]
    else
        error("No known key found for 3D orbit-space fast-ion distribution in filepath_Fos3D. Please correct and re-try.")
    end
    close(myfile)
end

## ------
# Loading tokamak equilibrium
verbose && println("Loading magnetic equilibrium... ")
M, wall, jdotb = nothing, nothing, nothing # Initialize global magnetic equilibrium variables
try
    global M; global wall; global jdotb # Declare global scope
    M, wall = read_geqdsk(filepath_equil,clockwise_phi=false) # Assume counter-clockwise phi-direction
    jdotb = M.sigma # The sign of the dot product between the plasma current and the magnetic field
catch # Otherwise, assume magnetic equilibrium is a saved .jld2 file
    global M; global wall; global jdotb; local myfile # Declare global scope
    myfile = jldopen(filepath_equil,false,false,false,IOStream)
    M = myfile["S"]
    wall = myfile["wall"]
    close(myfile)
    jdotb = (M.sigma_B0)*(M.sigma_Ip)
end

## ------
# Read the file containing the pre-computed weight function
projVel = false # Define this for its own sake
verbose && println("Loading weight function ("*filepath_W*")... ")
myfile = jldopen(filepath_W,false,false,false,IOStream)
if haskey(myfile,"W")
    W_correct = myfile["W"]
elseif haskey(myfile,"W_null")
    W_correct = myfile["W_null"]
elseif haskey(myfile,"W")
    W_correct = myfile["W"]
else
    error("No known keys in filepath_W for the orbit weight functions. Please correct and re-try.")
end
E_array = myfile["E_array"]
pm_array = myfile["pm_array"]
Rm_array = myfile["Rm_array"]
if haskey(myfile,"En_array")
    Ed_array = myfile["En_array"]
elseif haskey(myfile,"Ed_array")
    Ed_array = myfile["Ed_array"]
else
    error("weightsWebApp did not recognize known diagnostic data array type in provide 'filepath_W'. Please re-try another file.")
end
if haskey(myfile,"reaction")
    reaction = myfile["reaction"]
end
if haskey(myfile,"reaction_full")
    reaction_full = myfile["reaction_full"]
end
if haskey(myfile,"projVel")
    projVel = true
end
close(myfile)
fW_array = split(filepath_W,"_")
tokamak = fW_array[end-5]
TRANSP_id = fW_array[end-4]
timepoint = fW_array[end-3]
diagnostic = fW_array[end-2]
reaction_sscp = fW_array[end-1]

## ---------------------------------------------------------------------------------------------
## Safety check to ensure same dimensions for nullOrbits and weight functions
if showNullOrbs
    ### WRITE CODE HERE ###
end

## ---------------------------------------------------------------------------------------------
# Create null orbits matrix (if specified)
if showNullOrbs
    nullOrbs = zeros(length(Ed_array),length(E_array),length(pm_array),length(Rm_array))
    nullOrbs[nullOrbs_indices] .= 1.0
end

## ---------------------------------------------------------------------------------------------
# Determine fast-ion species from reaction
if (@isdefined reaction_full)
    thermal_species, FI_species = getFusionReactants(reaction_full)
elseif (@isdefined reaction)
    FI_species = getFusionReactants(reaction)[2]
else
    FI_species = FI_species # Already defined
end
FI_species = String(FI_species)

## ----------
# Check that the dimensions of topoBounds and F_os_4D match, and match the E-, pm- and Rm-arrays.
if plot_Fos
    println("Checking that dimensions are consistent... ")
    if !(size(F_os_4D[1,:,:,:])==size(topoBounds)) || !(size(F_os_4D,2) == length(E_array)) || !(size(F_os_4D,3) == length(pm_array)) || !(size(F_os_4D,4) == length(Rm_array))
        println("Size(F_os_4D): $(size(F_os_4D))")
        println("Size(topoBounds): $(size(topoBounds))")
        println("Length(E_array): $(length(E_array))")
        println("Length(pm_array): $(length(pm_array))")
        println("Length(Rm_array): $(length(Rm_array))")
        error("Dimensions of topoBounds and/or F_os_4D don't match given E-, pm- and/or Rm-arrays. Please correct and re-try.")
    end
end
verbose && println("Size of OWs: $(size(W_correct))")
verbose && println("Size of pm_array: $(size(pm_array))")
verbose && println("Size of Rm_array: $(size(Rm_array))")
verbose && println("Size of E_array: $(size(E_array))")
verbose && println("Size of Ed_array: $(size(Ed_array))")
if plot_Fos
    verbose && println("Size of orbit-space fast-ion distribution data: $(size(F_os_4D))")
end
if showNullOrbs
    verbose && println("Size of null orbits data: $(size(nullOrbs))")
end
verbose && println("Size of topoBounds: $(size(topoBounds))")

if showNullOrbs && plot_Fos
    if !(size(F_os_4D[1,:,:,:])==size(nullOrbs[1,:,:,:]))
        println("Size(nullOrbs): $(size(nullOrbs))")
        error("(E,pm,Rm) grid sizes of nullOrbs and F_os_4D don't match. Please correct and re-try.")
    end
end

## ----------
# Load the WF signal for plotting together with true S signal
if plot_WF
    if !isfile(filepath_WF)
        error("WF signal set for plot (plot_WF=true). But filepath_WF is invalid. Please correct and re-try.")
    end
    verbose && println("Loading WF signal ("*filepath_WF*")... ")
    myfile = jldopen(filepath_WF,false,false,false,IOStream)
    S_WF = myfile["S_WF"]
    if length(size(S_WF))==1
        S_WF = reshape(S_WF,(1,length(S_WF)))
    else
        S_WF = S_WF
    end
    if haskey(myfile,"En_array")
        Ed_array_WF = myfile["En_array"]
    elseif haskey(myfile,"Ed_array")
        Ed_array_WF = myfile["Ed_array"]
    else
        error("weightsWebApp did not recognize known diagnostic data array type in provide 'filepath_WF'. Please re-try another file.")
    end
    close(myfile)

    if length(size(S_WF))==2
    else
        S_WF = reshape(S_WF,(1,length(S_WF)))
    end
end

## ----------
# Load the true signal S for plotting and following neutron energy
if plot_S
    if !isfile(filepath_S)
        error("S signal set for plot (plot_S=true). But filepath_S is invalid. Please correct and re-try.")
    end
    verbose && println("Loading synthetic signal ("*filepath_S*")... ")
    if !specFileJLD2
        myfile = h5open(filepath_S)
        if diagnostic_name=="TOFOR"
            Ed_array_S = read(myfile["En_tofor"])
            spec = read(myfile["spec_tofor"])
        elseif diagnostic_name=="AB"
            Ed_array_S = read(myfile["En_ab"])
            spec = read(myfile["spec_ab"])
        elseif diagnostic_name==""
            error("Script does not know the keys to the signal spectrum and energy array. Please manually specify near line 410 of weightsWebApp.jl and comment out this error line.")
            Ed_array_S = read(myfile["INPUT ENERGY ARRAY KEY HERE"])
            spec = read(myfile["INPUT SPECTRUM KEY HERE"])
        else
            error("Invalid diagnostic. Please change diagnostic string in wFuncNfastionDistrWebApp.jl")
        end
    else
        myfile = jldopen(filepath_S,false,false,false,IOStream)
        if haskey(myfile,"S")
            spec = myfile["S"]
        elseif haskey(myfile,"spec")
            spec = myfile["spec"]
        else
            error("weightsWebApp did not recognize diagnostic signal data array key in provide 'filepath_S'. Please re-try another file.")
        end
        if haskey(myfile,"En_array")
            Ed_array_S = myfile["En_array"]
        elseif haskey(myfile,"Ed_array")
            Ed_array_S = myfile["Ed_array"]
        else
            error("weightsWebApp did not recognize diagnostic energy bin data array key in provide 'filepath_S'. Please re-try another file.")
        end
    end
    close(myfile)
end

## ----------
# Load sightlines
verbose && println("Loading "*diagnostic_name*" sightline ("*diagnostic_filepath*")... ")
if isfile(diagnostic_filepath)
    py"""
    VC = np.loadtxt($diagnostic_filepath)
    VC_RP = VC[:,8]
    VC_zP = VC[:,2]

    VC_X = VC[:,0]
    VC_Y = VC[:,1]
    """
    VC_RP = py"VC_RP"
    VC_zP = py"VC_zP"
    VC_X = py"VC_X"
    VC_Y = py"VC_Y"
else
    verbose && println("------> WARNING! No diagnostic viewing cone specified.")
end

E_array = vec(collect(E_array)) # Ensure type Array{Float64,1}
pm_array = vec(collect(pm_array)) # Ensure type Array{Float64,1}
Rm_array = vec(collect(Rm_array)) # Ensure type Array{Float64,1}
Ed_array = vec(collect(Ed_array)) # Ensure type Array{Float64,1}

## --------------------------------------------------------------------------
# Preparing utility for mapping (E,pm,Rm) to (E,Λ,Pϕ_n;σ)
if enable_COM
    if !(isfile(filepath_W_COM))
        if isfile(filepath_tm)
            verbose && print("Possibility to switch to (E,Λ,Pϕ_n;σ) requested. Topological map filepath specified. Attempting to load... ")
            myfile = jldopen(filepath_tm,false,false,false,IOStream)
            topoMap = myfile["topoMap"]
            E_array_tm = myfile["E_array"]
            pm_array_tm = myfile["pm_array"]
            Rm_array_tm = myfile["Rm_array"]
            close(myfile)
            verbose && println("Success!")
            verbose && print("Checking that (E,pm,Rm) grid points are consistent... ")
            if !((length(E_array)==length(E_array_tm)) && (length(pm_array)==length(pm_array_tm)) && (length(Rm_array)==length(Rm_array_tm)))
                verbose && println("")
                println("E_array:E_array_tm     $(length(E_array)):$(length(E_array_tm))")
                println("pm_array:pm_array_tm     $(length(pm_array)):$(length(pm_array_tm))")
                println("Rm_array:Rm_array_tm     $(length(Rm_array)):$(length(Rm_array_tm))")
                error("Dimensions of topological map does not match those of orbit weight functions! Please correct and re-start the app.")
            end
            verbose && println("Success!")
            verbose && print("Extracting valid orbits 3D indices from topological map... ")
            valid_orbit_indices = findall(x-> (x!=9.0) && (x!=7.0), topoMap) # 9 and 7 are the integers representing invalid and lost orbits in the calcTopoMap.jl script, respectively. We don't want them.
            verbose && println("Success!")
        else
            verbose && println("Switching to (E,Λ,Pϕ_n;σ) requested. But neither W_COM, nor topological map, filepath was not specified. Orbit grid will need to be computed (takes a long time)... ")
            valid_orbit_indices = :UNKNOWN
        end
    end
else
    verbose && println("Switching to (E,Λ,Pϕ_n;σ) will not be possible.")
end

## --------------------------------------------------------------------------
# Mapping orbit weight functions to (E,Λ,Pϕ_n;σ)
if enable_COM
    if !(isfile(filepath_W_COM))
        verbose && println(">>>>>>>>>>>>>>> Mapping orbit weight functions from (E,pm,Rm) to (E,Λ,Pϕ_n;σ)... <<<<<<<<<<<<<<<")
        W_correct_COM, E_array, Λ_array, Pϕ_n_array = os2COM(M, W_correct, E_array, pm_array, Rm_array, FI_species; nl=2*length(pm_array), npp=2*length(Rm_array), verbose=verbose, good_coords=valid_orbit_indices, wall=wall, extra_kw_args=extra_kw_args)
        verbose && println("Creating orbit weight matrix (E,Λ,Pϕ_n;σ) data... ")
        W_COM_inds = findall(x-> x>0.0, W_correct_COM)
        W_COM_inds_n_values = Array{Tuple{CartesianIndex{5},Float64}}(undef,length(W_COM_inds))
        for (ii,inds) in enumerate(W_COM_inds)
            W_COM_inds_n_values[ii] = (inds,W_correct_COM[Tuple(inds)...])
        end
        verbose && println("Saving orbit weight matrix (E,Λ,Pϕ_n;σ) data so that you do not have to re-map from (E,pm,Rm) next time... ")
        filepath_W_COM = folderpath_OWCF*"orbWeightsCOM_"*tokamak*"_"*TRANSP_id*"_at"*timepoint*"s_"*diagnostic*"_"*reaction_sscp*"_"*"$(size(W_correct_COM,1))x$(size(W_correct_COM,2))x$(size(W_correct_COM,3))x$(size(W_correct_COM,4))x2.jld2"
        verbose && println("Next time weightsWebApp.jl is run with same inputs, set filepath_W_COM = "*filepath_W_COM)
        myfile = jldopen(filepath_W_COM, true, true, false, IOStream)
        write(myfile,"W_COM_inds_n_values", W_COM_inds_n_values) # Save only non-zero indices an values, to save memory space
        write(myfile,"Ed_array",Ed_array)
        write(myfile,"E_array",E_array)
        write(myfile,"Lambda_array",Λ_array)
        write(myfile,"Pphi_n_array",Pϕ_n_array)
        close(myfile)
        # Removing extra CPU workers, if any, to minimize memory usage
        for i in workers()
            rmprocs(i)
        end

        # Release memory
        W_COM_inds = nothing
        W_COM_inds_n_values = nothing
    else
        verbose && print("Loading weight matrix in (E,Λ,Pϕ_n;σ) coordinates from filepath_W_COM... ")
        myfile = jldopen(filepath_W_COM, false, false, false, IOStream)
        if haskey(myfile,"W_COM_inds_n_values")    
            W_COM_inds_n_values = myfile["W_COM_inds_n_values"]
            Ed_array_COM = myfile["Ed_array"]
            E_array_COM = myfile["E_array"]
            Λ_array = myfile["Lambda_array"]
            Pϕ_n_array = myfile["Pphi_n_array"]
            close(myfile)
            verbose && println("Success!")

            if !(length(Ed_array)==length(Ed_array_COM))
                error("Number of diagnostic measurement bins of (E,Λ,Pϕ_n;σ) weight matrix in filepath_W_COM does not match number of diagnostic measurement bins of (E,pm,Rm) weight matrix in filepath_W. Please correct and re-start app.")
            end
            if !(length(E_array)==length(E_array_COM))
                error("Number of fast-ion energy grid points of (E,Λ,Pϕ_n;σ) weight matrix in filepath_W_COM does not match number of fast-ion energy grid points of (E,pm,Rm) weight matrix in filepath_W. Please correct and re-start app.")
            end

            verbose && println("Assembling data for orbit weight functions in (E,Λ,Pϕ_n;σ) space... ")
            W_correct_COM = zeros(length(Ed_array_COM),length(E_array_COM),length(Λ_array),length(Pϕ_n_array),2)
            for indValueTuple in W_COM_inds_n_values
                inds = indValueTuple[1]
                weight = indValueTuple[2]
                W_correct_COM[inds] = weight
            end

            # Release memory
            W_COM_inds_n_values = nothing
        else # Must be from os2com.jl
            if haskey(myfile,"Wtot")
                W_correct_COM = myfile["Wtot"]
                lk = "Wtot"
            elseif haskey(myfile,"W")
                W_correct_COM = myfile["W"]
                lk = "W"
            elseif haskey(myfile,"W4D")
                W_correct_COM = myfile["W4D"]
                lk = "W4D"
            else
                error("Unknown file key to load (E,Λ,Pϕ_n;σ) weight matrix from $(filepath_W_COM). Please correct and re-try.")
            end
            Ed_array_COM = myfile["Ed_array"]
            E_array_COM = myfile["E_array"]
            Λ_array = myfile["Lambda_array_$(lk)"]
            Pϕ_n_array = myfile["Pphi_n_array_$(lk)"]
            close(myfile)
        end
    end
end

## --------------------------------------------------------------------------
# Mapping topological boundaries to (E,Λ,Pϕ_n;σ)
if enable_COM
    B0 = norm(Equilibrium.Bfield(M,magnetic_axis(M)...)) # Tesla
    q = getSpeciesCharge(FI_species) # Coulomb
    if psi_bdry==0
        @warn "The magnetic flux at the last closed flux surface (LCFS) is found to be 0 for the magnetic equilibrium in $(filepath_equil). Pϕ_n=Pϕ/(q*|Ψ_w|) where Ψ_w=Ψ(mag. axis) is assumed instead of Ψ_w=Ψ(LCFS)."
        Ψ_w_norm = abs(psi_axis)
    else
        Ψ_w_norm = abs(psi_bdry)
    end
    verbose && println(">>>>>>>>>>>>>>> Mapping topological boundaries to (E,Λ,Pϕ_n;σ)... <<<<<<<<<<<<<<<")
    topoBounds_COM = zeros(size(W_correct_COM,2),size(W_correct_COM,3),size(W_correct_COM,4),2)
    for (iE,E) in enumerate(E_array)
        E_joule = (1000*GuidingCenterOrbits.e0)*E
        ones_carinds = findall(x-> x==1.0,topoBounds[iE,:,:])
        for carinds in ones_carinds
            pm, Rm = pm_array[carinds[1]], Rm_array[carinds[2]]
            EPRc = EPRCoordinate(M, E, pm, Rm; amu=getSpeciesAmu(FI_species), q=getSpeciesEcu(FI_species))
            myHc = HamiltonianCoordinate(M, EPRc)
            μ = myHc.mu; Λ = μ*B0/E_joule
            Pϕ = myHc.p_phi; Pϕ_n = Pϕ/(q*Ψ_w_norm)
            il, ipp = argmin(abs.(Λ_array .- Λ)), argmin(abs.(Pϕ_n_array .- Pϕ_n))
            if pm<0.0
                topoBounds_COM[iE, il, ipp, 1] = 1.0
            else
                topoBounds_COM[iE, il, ipp, 2] = 1.0
            end
        end
    end
end

## --------------------------------------------------------------------------
# Mapping null orbits to (E,Λ,Pϕ_n;σ)
if showNullOrbs && enable_COM
    verbose && println(">>>>>>>>>>>>>>> Mapping null orbits to (E,Λ,Pϕ_n;σ)... <<<<<<<<<<<<<<<")
    nullOrbs_COM = zeros(size(W_correct_COM))
    for iEd=1:size(W_correct_COM,1)
        for (iE,E) in enumerate(E_array)
            E_joule = (1000*GuidingCenterOrbits.e0)*E
            ones_carinds = findall(x-> x==1.0,nullOrbs[iEd,iE,:,:])
            for carinds in ones_carinds
                pm, Rm = pm_array[carinds[1]], Rm_array[carinds[2]]
                EPRc = EPRCoordinate(M, E, pm, Rm; amu=getSpeciesAmu(FI_species), q=getSpeciesEcu(FI_species))
                myHc = HamiltonianCoordinate(M, EPRc)
                μ = myHc.mu; Λ = μ*B0/E_joule
                Pϕ = myHc.p_phi; Pϕ_n = Pϕ/(q*Ψ_w_norm)
                il, ipp = argmin(abs.(Λ_array .- Λ)), argmin(abs.(Pϕ_n_array .- Pϕ_n))
                if pm<0.0
                    nullOrbs_COM[iEd, iE, il, ipp, 1] = 1.0
                else
                    nullOrbs_COM[iEd, iE, il, ipp, 2] = 1.0
                end
            end
        end
    end
end

## --------------------------------------------------------------------------
# Creating the array with which to switch between different reconstructions, if filepath_Fos3D was specified as a .jld2 file containing many fast-ion distributions (structured into a 4D data array)
if @isdefined F_os_4D
    Rec_array = collect(1:size(F_os_4D,1))
else
    Rec_array = [1]
end

## ------
# The web application
R_hfs = minimum(wall.r) # R-coord of high-field side wall
R_lfs = maximum(wall.r) # R-coord of low-field side wall
phi = collect(0:1:359).*(2*pi/180.0) # Toroidal angle
topview_R_hfs_x = (R_hfs).*cos.(phi)
topview_R_hfs_y = (R_hfs).*sin.(phi)
topview_R_lfs_x = (R_lfs).*cos.(phi)
topview_R_lfs_y = (R_lfs).*sin.(phi)
verbose && println("--- You can access the weightsWebApp via an internet web browser when you see 'Task (runnable)...' ")
verbose && println("--- When 'Task (runnable)...' has appeared, please visit the website localhost:$(port) ---")
verbose && println("--- Remember: It might take 1-2 minutes to load the webpage. Please be patient. ---")
function app(req)
<<<<<<< HEAD
    @manipulate for tokamak_wall = Dict("on" => true, "off" => false), include_Fos = Dict("on" => true, "off" => false), colorbar_scale = Dict("0.0-1.0" => "zero2one", "0.0-0.1" => "zero2aTenth"), phase_space = Dict("(E,μ,Pϕ;σ)" => :COM, "(E,pm,Rm)" => :OS), Ed=Ed_array, E=E_array, pm=pm_array, Rm=Rm_array, irec=Rec_array, save_plots = Dict("on" => true, "off" => false), show_coordinate = Dict("on" => true, "off" => false), show_delta_signal = Dict("on" => true, "off" => false)
=======
    @manipulate for tokamak_wall = Dict("on" => true, "off" => false), include_Fos = Dict("on" => true, "off" => false), colorbar_scale = Dict("0.0-1.0" => "zero2one", "0.0-0.1" => "zero2aTenth"), phase_space = Dict("(E,Λ,Pϕ_n;σ)" => :COM, "(E,pm,Rm)" => :OS), Ed=Ed_array, E=E_array, pm=pm_array, Rm=Rm_array, irec=Rec_array, save_plots = Dict("on" => true, "off" => false), show_coordinate = Dict("on" => true, "off" => false)
>>>>>>> main

        EPRc = EPRCoordinate(M, E, pm, Rm; amu=getSpeciesAmu(FI_species), q=getSpeciesEcu(FI_species))
        o = get_orbit(M,EPRc; wall=wall, extra_kw_args...)
        if phase_space==:COM
            E_joule = (1000*GuidingCenterOrbits.e0)*E
            myHc = HamiltonianCoordinate(M, EPRc)
            μ = myHc.mu; Λ = μ*B0/E_joule
            Pϕ = myHc.p_phi; Pϕ_n = Pϕ/(q*Ψ_w_norm)
            if (pm<0.0)
                iσ = 1 # Index 1 corresponds to σ=-1
            else
                iσ = 2 # Index 2 corresponds to σ=+1
            end
        end

        topview_o_x = cos.(o.path.phi).*(o.path.r)
        topview_o_y = sin.(o.path.phi).*(o.path.r)

        orb_color = :black
        orb_linestyle = :solid
        if o.class==:invalid
            orb_color = :gray
            orb_linestyle = :dash
        elseif o.class == :lost
            orb_color = :brown
        elseif o.class == :incomplete # If this happens, you are in trouble. Because it will likely take forever to calculate. Please just re-start the app instead.
            orb_color = :yellow
        elseif o.class == :trapped
            orb_color = :blue
        elseif o.class == :co_passing
            orb_color = :green
        elseif (o.class == :stagnation && o.coordinate.r>=magnetic_axis(M)[1]) # Regular stagnation orbit
            orb_color = :red
        elseif o.class == :potato
            orb_color = :orange
        elseif o.class == :ctr_passing
            orb_color = :purple
        elseif (o.class == :stagnation && o.coordinate.r<magnetic_axis(M)[1]) # Pinch (HFS stagnation) orbit
            orb_color = :pink
        else
            error("Something's gone terribly wrong!!! Orbit class unknown!")
        end

        ###### Cross-sectional plot ###### 
        if uppercase(diagnostic_name)=="TOFOR"
            plt_crs = Plots.plot(VC_RP,VC_zP,color=:green3, linewidth=1.2,label="")
        elseif uppercase(diagnostic_name)=="AB"
            plt_crs = Plots.plot(VC_RP,VC_zP,color=:red1, linewidth=1.2,label="")
        elseif diagnostic_name==""
            plt_crs = Plots.plot()
        else
            plt_crs = Plots.plot(VC_RP,VC_zP,color=:gray, linewidth=1.2,label="")
        end
        if (phase_space==:COM) && enable_COM
            plt_crs = Plots.plot!(title="E: $(round(E,digits=2)) keV  Λ: $(round(Λ, sigdigits=2))  Pϕ_n: $(round(Pϕ_n,sigdigits=2))")
        else # phase_space==:COM (OS = orbit space, COM = constants-of-motion)
            plt_crs = Plots.plot!(title="E: $(round(E,digits=2)) keV  pm: $(round(o.coordinate.pitch, digits=2))  Rm: $(round(o.coordinate.r,digits=2))")
        end
        plt_crs = Plots.plot!(o.path.r,o.path.z, label="$(o.class) orbit", color=orb_color, linestyle=orb_linestyle, linewidth=1.5)
        if tokamak_wall
            wall_dR = maximum(wall.r)-minimum(wall.r)
            plt_crs = Plots.plot!(wall.r,wall.z, label=tokamak*" wall", color=:black, linewidth=1.5,xaxis=[minimum(wall.r)-wall_dR/10,maximum(wall.r)+wall_dR])
        end
        plt_crs = Plots.scatter!([magnetic_axis(M)[1]],[magnetic_axis(M)[2]],label="Magnetic axis", mc=:grey, aspect_ratio=:equal, xlabel="R [m]", ylabel=" z[m]")
        plt_crs = Plots.scatter!([o.coordinate.r],[o.coordinate.z], mc=orb_color, label="(Rm,zm)")
        if save_plots
            if (phase_space==:COM) && enable_COM
                sigExt = iσ==1 ? "-1" : "+1" 
                png(plt_crs, "plt_crs_$(round(Ed, digits=2))_$(round(E, digits=2))_$(round(Λ, sigdigits=2))_$(round(Pϕ_n,sigdigits=2))_"*sigExt)
            else
                png(plt_crs, "plt_crs_$(round(Ed, digits=2))_$(round(E, digits=2))_$(round(o.coordinate.pitch, digits=2))_$(round(o.coordinate.r,digits=2))")
            end
        end

        # Find the E and Ed indices in E_array and Ed_array. For weights and fast-ion slice plot.
        Ei = (findall(x-> x==E,E_array))[1] # Should only be 1 exactly equal element
        Edi = (findall(x-> x==Ed,Ed_array))[1] # -||-

        # Extract correct topological boundaries and convert to vectors (for scatter plot)
        if (phase_space==:COM) && enable_COM
            ones_carinds = findall(x-> x==1.0, topoBounds_COM[Ei,:,:,iσ])
            Λ_scatvals_tb = zeros(length(ones_carinds))
            Pϕ_n_scatvals_tb = zeros(length(ones_carinds))
            for (ind,carinds) in enumerate(ones_carinds)
                Λ_scatvals_tb[ind] = Λ_array[carinds[1]]
                Pϕ_n_scatvals_tb[ind] = Pϕ_n_array[carinds[2]]
            end
        end
        if (phase_space==:OS) || (include_Fos && plot_Fos)
            ones_carinds = findall(x-> x==1.0,topoBounds[Ei,:,:])
            pm_scatvals_tb = zeros(length(ones_carinds))
            Rm_scatvals_tb = zeros(length(ones_carinds))
            for (ind,carinds) in enumerate(ones_carinds)
                pm_scatvals_tb[ind] = pm_array[carinds[1]]
                Rm_scatvals_tb[ind] = Rm_array[carinds[2]]
            end
        end

        # Extract correct null boundaries and convert to vectors (for scatter plot)
        if showNullOrbs
            if (phase_space==:COM) && enable_COM
                ones_carinds = findall(x-> x==1.0, nullOrbs_COM[Edi,Ei,:,:,iσ])
                Λ_scatvals_nb = zeros(length(ones_carinds))
                Pϕ_n_scatvals_nb = zeros(length(ones_carinds))
                for (ind,carinds) in enumerate(ones_carinds)
                    Λ_scatvals_nb[ind] = Λ_array[carinds[1]]
                    Pϕ_n_scatvals_nb[ind] = Pϕ_n_array[carinds[2]]
                end
            end
            if (phase_space==:OS) || (include_Fos && plot_Fos)
                ones_carinds = findall(x-> x==1.0,nullOrbs[Edi,Ei,:,:])
                pm_scatvals_nb = zeros(length(ones_carinds))
                Rm_scatvals_nb = zeros(length(ones_carinds))
                for (ind,carinds) in enumerate(ones_carinds)
                    pm_scatvals_nb[ind] = pm_array[carinds[1]]
                    Rm_scatvals_nb[ind] = Rm_array[carinds[2]]
                end
            end
        end

        if colorbar_scale == "zero2one"
            clims = (0.0,1.0)
        else
            clims = (0.0,0.1)
        end

        ###### Heatmap of the weight function slice plot ######
        if (phase_space==:COM) && enable_COM
            plt_weights = Plots.heatmap(Pϕ_n_array,Λ_array,(W_correct_COM[Edi,Ei,:,:,iσ])./maximum(W_correct_COM[Edi,Ei,:,:,iσ]),colorbar=true,title="W ($(round(maximum(W_correct_COM[Edi,Ei,:,:,iσ]),sigdigits=4)) = 1.0)", clims=clims, fillcolor=cgrad([:white, :darkblue, :green, :yellow, :orange, :red]), ylims=extrema(Λ_scatvals_tb), xlims=extrema(Pϕ_n_scatvals_tb))
        else
            plt_weights = Plots.heatmap(Rm_array,pm_array,(W_correct[Edi,Ei,:,:])./maximum(W_correct[Edi,Ei,:,:]),colorbar=true,title="W ($(round(maximum(W_correct[Edi,Ei,:,:]),sigdigits=4)) = 1.0)", clims=clims, fillcolor=cgrad([:white, :darkblue, :green, :yellow, :orange, :red]))
        end
        # Plot a scatter-plot of the topological boundaries for that specific fast-ion energy slice
        ms = save_plots ? 2.6 : 1.8
        if (phase_space==:COM) && enable_COM
            plt_weights = Plots.scatter!(Pϕ_n_scatvals_tb,Λ_scatvals_tb,markersize=ms,leg=false,markercolor=:black, xlabel="Pϕ_n", ylabel="Λ")
        else
            plt_weights = Plots.scatter!(Rm_scatvals_tb,pm_scatvals_tb,markersize=ms,leg=false,markercolor=:black, xlabel="Rm [m]", ylabel="pm")
        end
        if showNullOrbs
            # Plot a scatter-plot of the null orbits for that specific diagnostic and fast-ion energy
            if (phase_space==:COM) && enable_COM
                plt_weights = Plots.scatter!(Pϕ_n_scatvals_nb,Λ_scatvals_nb,markersize=5.0,leg=false,markercolor=:teal, markershape=:xcross, linewidth=2.5)
            else
                plt_weights = Plots.scatter!(Rm_scatvals_nb,pm_scatvals_nb,markersize=5.0,leg=false,markercolor=:teal, markershape=:xcross, linewidth=2.5)
            end
        end
        if show_coordinate
            # Plot the orbit coordinate
            if (phase_space==:COM) && enable_COM
                plt_weights = Plots.scatter!([Pϕ_n],[Λ],markershape=:circle,mc=:white,markerstrokewidth=3.0,markerstrokealpha=1.0, markerstrokecolor=orb_color,markersize=5.0) # orbit coordinate marker
            else
                plt_weights = Plots.scatter!([Rm],[pm],markershape=:circle,mc=:white,markerstrokewidth=3.0,markerstrokealpha=1.0, markerstrokecolor=orb_color,markersize=5.0) # orbit coordinate marker
            end
        end
        if save_plots
            plt_weights = Plots.plot!(dpi=600)
            if (phase_space==:COM) && enable_COM
                sigExt = iσ==1 ? "-1" : "+1" 
                png(plt_weights, "plt_weights_$(round(Ed, digits=2))_$(round(E, digits=2))_$(round(Λ, digits=2))_$(round(Pϕ_n,digits=2))_"*sigExt)
            else
                png(plt_weights, "plt_weights_$(round(Ed, digits=2))_$(round(E, digits=2))_$(round(o.coordinate.pitch, digits=2))_$(round(o.coordinate.r,digits=2))")
            end
        end

        ####### Plot of the fast-ion distribution slice, if included
        if include_Fos && plot_Fos # If you want to plot it (include_Fos) and it's been loaded (plot_Fos)
            plt_Fos = Plots.heatmap(Rm_array, pm_array, (F_os_4D[irec,Ei,:,:])./maximum(F_os_4D[irec,Ei,:,:]), colorbar=true, title="Fast-ion distribution slice  ($(round(maximum(F_os_4D[irec,Ei,:,:]), sigdigits=4)) = 1.0)", clims=clims, fillcolor=cgrad([:white, :darkblue, :green, :yellow, :orange, :red]))
            plt_Fos = Plots.scatter!(Rm_scatvals_tb,pm_scatvals_tb,markersize=1.8,leg=false,markercolor=:black, xlabel="Rm [m]", ylabel="pm")
            if showNullOrbs
                plt_Fos = Plots.scatter!(Rm_scatvals_nb,pm_scatvals_nb,markersize=2.5,leg=false,markercolor=:teal, markershape=:xcross)
            end
            if show_coordinate
                plt_Fos = Plots.scatter!([Rm],[pm],markershape=:circle,mc=:white,markerstrokewidth=1.1,markerstrokealpha=1.0, markerstrokecolor=orb_color,markersize=5.0) # orbit coordinate marker
            end
        end
        if save_plots
            plt_Fos = Plots.plot!(dpi=600)
            if (phase_space==:COM) && enable_COM
                sigExt = iσ==1 ? "-1" : "+1" 
                png(plt_Fos, "plt_Fos_$(round(Ed, digits=2))_$(round(E, digits=2))_$(round(Λ, sigdigits=2))_$(round(Pϕ_n,sigdigits=2))_"*sigExt)
            else
                png(plt_Fos, "plt_Fos_$(round(Ed, digits=2))_$(round(E, digits=2))_$(round(o.coordinate.pitch, digits=2))_$(round(o.coordinate.r,digits=2))")
            end
        end

        ####### Plot of the WF signal density, if FI distribution is included
        if include_Fos && plot_Fos
            WF_dens = (W_correct[Edi,Ei,:,:]) .* (F_os_4D[irec,Ei,:,:]) # WF, point-wise (without summation)
            plt_WFdens = Plots.heatmap(Rm_array, pm_array, (WF_dens)./maximum(WF_dens), colorbar=true, title="WF-density slice  ($(round(maximum(WF_dens), sigdigits=4)) = 1.0)", clims=clims, fillcolor=cgrad([:white, :darkblue, :green, :yellow, :orange, :red]))
            plt_WFdens = Plots.scatter!(Rm_scatvals_tb,pm_scatvals_tb,markersize=1.8,leg=false,markercolor=:black, xlabel="Rm [m]", ylabel="pm")
            if showNullOrbs
                plt_WFdens = Plots.scatter!(Rm_scatvals_nb,pm_scatvals_nb,markersize=2.5,leg=false,markercolor=:teal, markershape=:xcross)
            end
            if show_coordinate
                plt_WFdens = Plots.scatter!([Rm],[pm],markershape=:circle,mc=:white,markerstrokewidth=1.1,markerstrokealpha=1.0, markerstrokecolor=orb_color, markersize=5.0) # orbit coordinate marker
            end
        end
        if save_plots
            plt_WFdens = Plots.plot!(dpi=600)
            if (phase_space==:COM) && enable_COM
                sigExt = iσ==1 ? "-1" : "+1" 
                png(plt_WFdens, "plt_WFdens_$(round(Ed, digits=2))_$(round(E, digits=2))_$(round(Λ, sigdigits=2))_$(round(Pϕ_n,sigdigits=2))_"*sigExt)
            else
                png(plt_WFdens, "plt_WFdens_$(round(Ed, digits=2))_$(round(E, digits=2))_$(round(o.coordinate.pitch, digits=2))_$(round(o.coordinate.r,digits=2))")
            end
        end

        # Plot of the top view of the tokamak
        plt_top = Plots.plot(topview_R_lfs_x,topview_R_lfs_y, label=tokamak*" wall", color=:black, linewidth=1.5)
        plt_top = Plots.plot!(topview_R_hfs_x,topview_R_hfs_y, label="", color=:black,linewidth=1.5, aspect_ratio=:equal, title="Top view")
        if uppercase(diagnostic_name)=="TOFOR"
            plt_top = Plots.plot!(VC_X,VC_Y,color=:green3, linewidth=1.2,label="")
        elseif uppercase(diagnostic_name)=="AB"
            plt_top = Plots.plot!(VC_X,VC_Y,color=:red1, linewidth=1.2,label="")
        elseif diagnostic_name==""
        else
            plt_top = Plots.plot!(VC_X,VC_Y,color=:gray, linewidth=1.2,label="")
        end
        plt_top = Plots.plot!(topview_o_x,topview_o_y,label="$(o.class) orbit", color=orb_color, linestyle=orb_linestyle, linewidth=1.5)
        if save_plots
            if (phase_space==:COM) && enable_COM
                sigExt = iσ==1 ? "-1" : "+1" 
                png(plt_top, "plt_top_$(round(Ed, digits=2))_$(round(E, digits=2))_$(round(Λ, sigdigits=2))_$(round(Pϕ_n,sigdigits=2))_"*sigExt)
            else
                png(plt_top, "plt_top_$(round(Ed, digits=2))_$(round(E, digits=2))_$(round(o.coordinate.pitch, digits=2))_$(round(o.coordinate.r,digits=2))")
            end
        end

        # Signal plot
        if plot_S || plot_WF || show_delta_signal
            if uppercase(diagnostic_name)=="TOFOR"
                sig_color = :green3
            elseif uppercase(diagnostic_name)=="AB"
                sig_color = :red1
            elseif uppercase(diagnostic_name)=="KM6T"
                sig_color = :blue1
            else
                sig_color = :gray
            end
            plt_sig = Plots.plot()
            if plot_WF
                plt_sig = Plots.scatter!(Ed_array_WF, S_WF[irec,:] ./maximum(S_WF[irec,:]), markerstrokealpha=1.0, markerstrokecolor=sig_color, markercolor=:white, markerstrokewidth=1.5, label="WF")
                Edi_WF = (findfirst(x-> x>=Ed,Ed_array_WF))[1] # Should be perfect match with weights
                plt_sig = Plots.scatter!([Ed_array_WF[Edi_WF]],[(S_WF[irec,:])[Edi_WF] ./maximum(S_WF[irec,:])],markersize=5.0,markercolor=sig_color, label="WF: $(round((S_WF)[Edi_WF],sigdigits=4))")
                plt_sig = Plots.plot!(xlabel=xlabel,ylabel="Normalized signal [a.u.]", legend=true, title="Ed: $(round(Ed,digits=4)) "*Ed_units)
            end
            if plot_S
                plt_sig = Plots.plot!(Ed_array_S, spec ./maximum(spec), color=sig_color, linewidth=2.0, label="S")
                Edi_S = (findfirst(x-> x>=Ed,Ed_array_S))[1]
                plt_sig = Plots.scatter!([Ed_array_S[Edi_S]],[spec[Edi_S]],markersize=5.0,markercolor=sig_color, label="S: $(round((spec)[Edi_S],sigdigits=4))")
                plt_sig = Plots.plot!(xlabel=xlabel,ylabel="Normalized signal [a.u.]", legend=true, title="Ed: $(round(Ed,digits=4)) "*Ed_units)
            end
            if show_delta_signal
                pmci = argmin(abs.(pm_array .- pm)) # Find the closest match
                Rmci = argmin(abs.(Rm_array .- Rm)) # Find the closest match
                delta_signal = W_correct[:,Ei,pmci,Rmci] ./maximum(W_correct[:,Ei,pmci,Rmci])
                plt_sig = Plots.plot!(Ed_array, delta_signal, label="Orbit-pixel signal")
                plt_sig = Plots.scatter!([Ed_array[Edi]],[delta_signal[Edi]],markersize=5.0, markercolor=:black, label="Orbit-pixel sig.: $(delta_signal[Edi])", xlims=extrema(Ed_array), ylims=[-0.1,1.2])
                plt_sig = Plots.plot!(xlabel=xlabel,ylabel="Normalized signal [a.u.]", legend=true, title="Ed: $(round(Ed,digits=4)) "*Ed_units)
            end
        else
            Edi = (findfirst(x-> x>=Ed,Ed_array))[1]
            plt_sig = Plots.scatter([Ed_array[Edi]],[0.0],markersize=5.0, markercolor=:black, title="Ed: $(round(Ed,digits=4)) "*Ed_units,label="", legend=false, xlims=[minimum(Ed_array),maximum(Ed_array)], ylims=[-1.0,1.0],xlabel=xlabel,ylabel="No diagnostic signal specified [a.u.]")
        end
        if save_plots
            plt_sig = Plots.plot!(title="", legend=false)
            if (phase_space==:COM) && enable_COM
                sigExt = iσ==1 ? "-1" : "+1" 
                png(plt_top, "plt_top_$(round(Ed, digits=2))_$(round(E, digits=2))_$(round(Λ, sigdigits=2))_$(round(Pϕ_n,sigdigits=2))_"*sigExt)
            else
                png(plt_sig, "plt_sig_$(round(Ed, digits=2))_$(round(E, digits=2))_$(round(o.coordinate.pitch, digits=2))_$(round(o.coordinate.r,digits=2))")
            end
            plt_sig = Plots.plot!(title="Ed: $(round(Ed,digits=4)) keV", legend=true)
        end

        if include_Fos && plot_Fos # Three rows with two plots each (total 6)
            vbox(vskip(1em),
                md"*Please note, fast-ion distribution and WF-density will not be visualzed in (E,Λ,Pϕ_n;σ). A Jacobian going from (E,pm,Rm) is currently not implemented for the weightsWebApp.*",
                md"Also, please note that the weightsWebApp.jl gets rather slow when all the inputs have been specified.",
                md"If so, it is better to click new values on the sliders, rather than drag the sliders.",
                vskip(1em),
                hbox(Plots.plot(plt_sig),Plots.plot(plt_weights)),
                hbox(Plots.plot(plt_Fos),Plots.plot(plt_WFdens)),
                hbox(Plots.plot(plt_top),Plots.plot(plt_crs))
            )
        else
            vbox(vskip(1em), # Two rows with two plots each (total 4)
                hbox(Plots.plot(plt_sig),Plots.plot(plt_weights)),
                hbox(Plots.plot(plt_top),Plots.plot(plt_crs))
            )
        end
    end
end

webio_serve(page("/",app), port)