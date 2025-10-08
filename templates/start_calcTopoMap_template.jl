################################ start_calcTopoMap_template.jl #########################################

#### Description:
# This script contains all the inputs for the calcTopoMap.jl script. It defines the computational
# environment (batch job, cores etc), loads the inputs and then executes the calcTopoMap.jl script.
# The inputs are as follows:
#
# batch_job_SLURM - If true, the script assumes that it will be executed as part of an high-performance computational 
#                   cluster batch job within the SLURM cluster environment. The script will act accordingly. - Bool
# distributed - If true, parallel computing with multiple CPU cores will be used. - Bool
# folderpath_OWCF - The path to the OWCF folder - String
# numOcores - The number of CPU cores that will be used if distributed is set to true. - Int64
#
# distinguishIncomplete - If true, then incomplete orbits will be given integer 6 instead of 9.
# distinguishLost - If true, then lost orbits will be given integer 7 instead of 9.
# filepath_equil - The path to the .eqdsk-file (or .geqdsk/.jld2-file) with the tokamak magnetic equilibrium and geometry - String
# filepath_OG - The path to the .jld2 file containing the output of calcOrbGrid.jl - String
# filepath_W - The path to the .jld2/.h5/.hdf5 weights file to extract orbit grid info from - String
# folderpath_o - The path to the folder where the results will be saved - String
# includeExtractTopoBounds - If true, then the topological boundaries will be computed as well - Bool
# keyname_diagHDF5 - The HDF5 key identifier with which to load the weights file. if !(weightsFileJLD2). The 'E_array', 'pm_array' and 'Rm_array' key identifiers will be assumed for the orbit-space grid points - String
# saveTransitTimeMaps - If true, the output .jld2 file will also contain two 3D array, both of size(orbit grid), with poloidal and toroidal transit times for the valid orbits - Bool
# timepoint - The timepoint of the tokamak shot for the magnetic equilibrium. Format "XX,YYYY" where XX are seconds and YYYY are decimals - String
# tokamak - The tokamak for which the topological map is created - String
# TRANSP_id - The TRANSP id to extract tokamak/equilibrium data from - String
# verbose - If true, lots of information will be printed during execution - Bool
# visualizeProgress - If false, progress bar will not be displayed during computations - Bool
# weightsFileJLD2 - If true, it is assumed that the weights file is in JLD2 format - Bool

# Then, there are inputs that only need to be specified if (!isfile(filepath_W) && !isfile(filepath_OG)):
# E_array - If you want, you can manually specify a vector of fast-ion energies for which to compute the topological map - Array{Float64,1}
# FI_species - The fast-ion particle species that the topological map is to be computed for - String
# pm_array - If you want, you can manually specify a vector of pitch maximum values for which to compute the topological map - Array{Float64,1}
# Rm_array - If you want, you can manually specify a vector of radius maximum values for which to compute the topological map - Array{Float64,1}

# Then, there are inputs that only need to be specified if E_array, pm_array and/or Rm_array are set to nothing:
# Emin - The lower bound for the fast-ion energy in orbit space. if (!isfile(filepath_W) && !isfile(filepath_OG)) - Float64
# Emax - The upper bound for the fast-ion energy in orbit space. if (!isfile(filepath_W) && !isfile(filepath_OG)) - Float64
# inclPrideRockOrbs - If true, then pride rock (counter-stagnation) orbits will be included
# nE - The number of fast-ion energy grid points in orbit space. if (!isfile(filepath_W) && !isfile(filepath_OG)) - Int64
# npm - The number of fast-ion pm grid points in orbit space. if (!isfile(filepath_W) && !isfile(filepath_OG)) - Int64
# nRm - The number of fast-ion Rm grid points in orbit space. if (!isfile(filepath_W) && !isfile(filepath_OG)) - Int64
# pm_min - The lower bound for the pitch maximum in orbit space - Float64
# pm_max - The upper bound for the pitch maximum in orbit space - Float64
# Rm_min - The lower boundary for the orbit-grid Rm values - Float64
# Rm_max - The upper boundary for the orbit-grid Rm values - Float64
#
# extra_kw_args - Here you can define extra keyword arguments to be used when integrating the equations-of-motion
# for the guiding-center orbits. Please consult the OWCF manual for further info on this.

#### Outputs
# -

#### Saved files
# topoMap_[tokamak]_[TRANSP ID]_[FI species]_[timestamp]_[nE]x[npm]x[nRm].jld2
#   This saved file will have the fields:
#   topoMap - The saved topological map. Dimensions are size(orbit grid) - Array{Float64,3}
#   E_array - The fast-ion energy grid array used for orbit space - Array{Float64,1}
#   pm_array - The fast-ion pm grid array used for orbit space - Array{Float64,1}
#   Rm_array - The fast-ion Rm grid array used for orbit space - Array{Float64,1}
# And if saveTransitTimeMaps was set to true, the saved file will also have the fields:
#   polTransTimes - The poloidal transit times for all valid orbits. Dimensions are size(orbit grid) - Array{Float64,3}
#   torTransTimes - The toroidal transit times for all valid orbits. Dimensions are size(orbit grid) - Array{Float64,3}

# Script written by Henrik Järleblad. Last maintained 2025-10-07.
########################################################################################

## First you have to set the system specifications
using Distributed # Needed, even though distributed might be set to false. This is to export all inputs to all workers right away, if needed.
batch_job_SLURM = false
distributed = true
folderpath_OWCF = "" # OWCF folder path. Finish with '/'
numOcores = 4 # When executing the script as part of an HPC cluster batch job, the number of CPU cores will be detected and set automatically (i.e. the value of the numOcores variable does not matter)

## Navigate to the OWCF folder and activate the OWCF environment
cd(folderpath_OWCF)
using Pkg
Pkg.activate(".")

## If running as a batch job on a SLURM HPC cluster with CPUs
if batch_job_SLURM
    # Load the SLURM package and add the CPU cores
    using SlurmClusterManager
    addprocs(SlurmManager()) # The number of CPU cores are detected automatically
    @everywhere println("X-wing $(rand(["Red","Green","Blue"])) $(myid()) reporting for duty at $(gethostname())")
end

## If running locally and multi-threaded
if !batch_job_SLURM && distributed # Assume you are executing the script on a local laptop (/computer)
    println("Adding processes... ")
    addprocs(numOcores-(nprocs()-1)) # If you didn't execute this script as an HPC cluster job, then you need to add processors like this. Add all remaining available cores.
    # The '-(nprocs()-1)' part is simply to ensure to extra processes are added, in case script needs to be restarted on a local computer
end

## -----------------------------------------------------------------------------
@everywhere begin
    distinguishIncomplete = false
    distinguishLost = false
    filepath_equil = ""
    filepath_OG = "" # You can choose to load the orbit-space grid points from a pre-computed orbit-grid .jld2-file...
    filepath_W = "" # ...or you can choose to load the orbit-space grid points from an orbit weights .jld2- OR .hdf5-file
    folderpath_o = "../OWCF_results/template/" # Output folder path. Finish with '/'
    folderpath_OWCF = $folderpath_OWCF # Set path to OWCF folder to same as main process (hence the '$')
    includeExtractTopoBounds = false # If true, then a .jld2 file with the boundaries of the topological map will be extracted and saved as well
    keyname_diagHDF5 = "" # Needed to extract info from .hdf5-file, if weightsFileJLD2 is set to false. 
    saveTransitTimeMaps = false # If true, then poloidal and toroidal transit times for all valid orbits will be included in the results file (see description above)
    timepoint = nothing # "XX,YYY". If unknown, just leave as nothing. The algorithm will try to figure it out automatically.
    tokamak = "" # If unknown, just set ""
    TRANSP_id = "" # If unknown, just set ""
    verbose = true # If set to true, the script will be very extroverted! Yay!
    visualizeProgress = false # If set to true, a progress bar will be displayed during distributed computation
    weightsFileJLD2 = false # Set to true, if isfile(filepath_W), and filepath_W is .jld2 in file format

    # \/ Only necessary if (!isfile(filepath_W) && !isfile(filepath_OG))
    E_array = nothing # for example [1.4, 100.0, 2500.0]
    FI_species = "D" # for example "D" for deuterium, "p" for proton, "T" for tritium, "3he" for helium-3
    pm_array = nothing # for example [-0.5, 0.5]
    Rm_array = nothing # for example [3.0 3.5]

    # \/ Only necessary if (!isfile(filepath_W) && !isfile(filepath_OG)) and E_array, pm_array and/or Rm_array are set to nothing
    Emin = 0.0 # keV
    Emax = 000.0 # keV
    inclPrideRockOrbs = true # If true, then pride rock orbits will be included. Otherwise, minimum(Rm) = magnetic axis.
    nE = 0
    npm = 0
    nRm = 0
    pm_min = -1.0
    pm_max = 1.0
    Rm_min = nothing # Automatically use magnetic axis or 4/5th from HFS wall to magnetic axis. Override by specifying together with Rm_max
    Rm_max = nothing # Automatically use LFS wall. Override by specifying together with Rm_min
    # /\ Only necessary if E_array, pm_array and/or Rm_array is set to nothing
    # /\ Only necessary if (!isfile(filepath_W) && !isfile(filepath_OG))

    # EXTRA KEYWORD ARGUMENTS BELOW
    extra_kw_args = Dict(:limit_phi => true, :maxiter => 0, :toa => true)
    # limits the number of toroidal turns for orbits
    # The orbit integration algorithm will try progressively smaller timesteps these number of times
    # toa stands for "try only adaptive" orbit ingration algorithm
end

## -----------------------------------------------------------------------------
# Change directory to OWCF-folder on all external processes. Activate the environment there to ensure correct package versions
# as specified in the Project.toml and Manuscript.toml files. Do this for every 
# CPU processor (@everywhere)
@everywhere begin
    using Pkg
    cd(folderpath_OWCF)
    Pkg.activate(".")
end

## -----------------------------------------------------------------------------
# Then you execute the script
include("calcTopoMap.jl")
## -----------------------------------------------------------------------------
# Then you clean up after yourself
if batch_job_SLURM || distributed
    for i in workers()
        rmprocs(i)
    end
end