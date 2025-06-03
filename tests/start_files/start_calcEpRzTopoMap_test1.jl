################################ start_calcEpRzTopoMap_template.jl #########################################

#### Description:
# This script contains all the inputs for the calcEpRzTopoMap.jl script. It defines the computational
# environment (batch job, cores etc), loads the inputs and then executes the calcEpRzTopoMap.jl script.
# The inputs are as follows:
#
# batch_job - If true, the script assumes that it will be executed via an HPC batch job. The script will act accordingly. - Bool
# distributed - If true, parallel computing with multiple CPU cores will be used. - Bool
# folderpath_OWCF - The path to the OWCF folder - String
# numOcores - The number of CPU cores that will be used if distributed is set to true. - Int64
#
# distinguishIncomplete - If true, then incomplete orbits will be given integer 6 instead of 9 - Bool
# distinguishLost - If true, then lost orbits will be given integer 7 instead of 9 - Bool
# distrFileJLD2 - If true, it is assumed that the weights file is in JLD2 format - Bool
# FI_species - The fast-ion particle species that the topological map is to be computed for - String
# filepath_distr - The path to the .jld2/.h5 weights file to extract particle grid info from (if useDistrFile) - String
# filepath_equil - The path to the .eqdsk-file (or .geqdsk/.jld2-file) with the tokamak magnetic equilibrium and geometry - String
# folderpath_o - The path to the folder where the results will be saved - String
# saveTransitTimeMaps - If set to true, the poloidal and toroidal transit times will be saved for every valid orbit - Bool
# saveXYZJacobian - If true, the Jacobian from (x,y,z,vx,vy,vz) to (E,p,R,z) will be computed - Bool
# timepoint - The timepoint of the tokamak shot for the magnetic equilibrium. Format XX,YYYY where XX are seconds and YYYY are decimals - String
# tokamak - The tokamak for which the topological map is created - String
# TRANSP_id - The TRANSP id to extract tokamak/equilibrium data from - String
# useDistrFile - If true, then distribution file will be used to extract particle grid info from - Bool
# verbose - If true, lots of information will be printed during execution - Bool
# visualizeProgress - If false, progress bar will not be displayed during computations - Bool

# Then, there are inputs that only need to be specified if useDistrFile is set to false:
# E_array - If you want, you can manually specify a vector of fast-ion energies for which to compute the topological map - Array{Float64,1}
# p_array - If you want, you can manually specify a vector of pitch values for which to compute the topological map - Array{Float64,1}
# R_array - If you want, you can manually specify a vector of radius values for which to compute the topological map - Array{Float64,1}
# z_array - If you want, you can manually specify a vector of vertical values for which to compute the topological map - Array{Float64,1}

# Then, there are inputs that only need to be specified if E_array, p_array, R_array and/or z_array are set to nothing:
# Emin - The lower bound for the fast-ion energy in particle space. if !(useDistrFile). - Float64
# Emax - The upper bound for the fast-ion energy in particle space. if !(useDistrFile). - Float64
# nE - The number of fast-ion energy grid points in particle space. if !(useDistrFile). - Int64
# np - The number of fast-ion pitch grid points in particle space. if !(useDistrFile). - Int64
# nR - The number of fast-ion R grid points in particle space. if !(useDistrFile). - Int64
# nz - The number of fast-ion z grid points in particle space. if !(useDistrFile). - Int64
# p_min - The lower bound for the pitch in particle space - Float64
# p_max - The upper bound for the pitch in particle space - Float64
# R_min - The lower boundary for the particle-grid R values - Float64
# R_max - The upper boundary for the particle-grid R values - Float64
# z_min - The lower boundary for the particle-grid z values - Float64
# z_max - The upper boundary for the particle-grid z values - Float64
#
# extra_kw_args - Here you can define extra keyword arguments to be used when integrating the equations-of-motion
# for the guiding-center orbits. Please consult the OWCF manual for further info on this.

#### Outputs
# -

#### Saved files
# topoMap_[tokamak]_[TRANSP ID]_[FI species]_[nE]x[np]x[nR]x[nz].jld2
#   This saved file will have the fields:
#   topoMap - The saved topological map. Dimensions are size(particle space) - Array{Float64,4}
#   E_array - The fast-ion energy grid array used for particle space - Array{Float64,1}
#   p_array - The fast-ion p grid array used for particle space - Array{Float64,1}
#   R_array - The fast-ion R grid array used for particle space - Array{Float64,1}
#   z_array - The fast-ion z grid array used for particle space - Array{Float64,1}
#   FI_species - The fast-ion species. "D", "3he", "T" etc - String
#   filepath_distr - If used, the filepath to the fast-ion distribution .jld2 file - String

# Script written by Henrik Järleblad. Last maintained 2024-06-25.
########################################################################################

## First you have to set the system specifications
using Distributed # Needed, even though distributed might be set to false. This is to export all inputs to all workers right away, if needed.
batch_job = false
distributed = false
if !(@isdefined folderpath_OWCF)
    folderpath_OWCF = reduce(*,map(x-> "/"*x,split(@__DIR__,"/")[2:end-2]))*"/" # We know that the test start file is located in the OWCF/tests/start_files/ folder. Deduce the full OWCF folder path from that information
end
if !(@isdefined plot_test_results)
    plot_test_results = false
end
if !isdir(folderpath_OWCF*"tests/outputs/")
    print("The folder $(folderpath_OWCF)tests/outputs/ does not exist. Creating... ")
    mkdir(folderpath_OWCF*"tests/outputs")
    println("ok!")
end
@everywhere begin
    using Pkg
    cd(folderpath_OWCF)
    Pkg.activate(".")
end
#numOcores = 4 # When executing script via HPC cluster job, make sure you know how many cores you have requested for your batch job

## Navigate to the OWCF folder and activate the OWCF environment
#cd(folderpath_OWCF)
#using Pkg
#Pkg.activate(".")

## If running as a batch job on a SLURM CPU cluster
#if batch_job && distributed
#    # Load the SLURM CPU cores
#    using ClusterManagers
#    addprocs(SlurmManager(numOcores))
#    hosts = []
#    pids = []
#    for i in workers()
#        host, pid = fetch(@spawnat i (gethostname(), getpid()))
#        push!(hosts, host)
#        push!(pids, pid)
#    end
#    @show hosts
#end

## If running locally and multi-threaded
#if !batch_job && distributed # Assume you are executing the script on a local laptop (/computer)
#    println("Adding processes... ")
#    addprocs(numOcores-(nprocs()-1)) # If you didn't execute this script as an HPC cluster job, then you need to add processors like this. Add all remaining available cores.
#    # The '-(nprocs()-1)' part is simply to ensure to extra processes are added, in case script needs to be restarted on a local computer
#end

## -----------------------------------------------------------------------------
@everywhere begin
    folderpath_OWCF = $folderpath_OWCF

    distinguishIncomplete = true
    distinguishLost = true
    distrFileJLD2 = false # Set to true, if useDistrFile is set to true, and filepath_distr is .jld2 in file format
    FI_species = "Be" # for example "D" for deuterium, "p" for proton, "T" for tritium, "3he" for helium-3
    filepath_distr = ""
    filepath_equil = folderpath_OWCF*"equilibrium/JET/g99971/g99971_474-48.9.eqdsk"
    filename_o = "calcEpRzTopoMap_test1"
    folderpath_o = folderpath_OWCF*"tests/outputs/" # Output folder path. Finish with '/'
    plot_results = $plot_test_results
    saveTransitTimeMaps = true # If true, then calcEpRzTopoMap.jl will save 4D data of the poloidal and toroidal transit times for the valid orbits
    saveXYZJacobian = false # Set to true, if Jacobian from (x,y,z,vx,vy,vz) to (E,p,R,z) should be computed and saved
    timepoint = nothing # If unknown, just leave as nothing. The algorithm will try to figure it out automatically.
    tokamak = "JET" # If unknown, just set ""
    TRANSP_id = "" # If unknown, just set ""
    useDistrFile = false # Set to true, if you want to load the particle grid from the distribution file filepath_distr
    verbose = true # If set to true, the script will be very extroverted! Yay!
    visualizeProgress = false # If set to true, a progress bar will be displayed during distributed computation

    # \/ Only necessary if !(useDistrFile)
    E_array = nothing # for example [1.4, 100.0, 2500.0]
    p_array = nothing # for example [-0.5, 0.5]
    R_array = nothing # for example [2.5 3.5]
    z_array = nothing # for example [-1.5 1.5]

    # \/ Only necessary if E_array, p_array, R_array and/or z_array are set to nothing
    Emin = 1.0 # keV
    Emax = 10.0 # keV
    nE = 11
    np = 12
    nR = 5
    nz = 6
    p_min = -1.0
    p_max = 1.0
    R_min = nothing # Automatically HFS wall to LFS wall. Override by specifying together with R_max
    R_max = nothing # Automatically HFS wall to LFS wall. Override by specifying together with R_min
    z_min = nothing # Automatically divertor bottom to vacuum vessel roof. Override by specifying together with z_max
    z_max = nothing # Automatically divertor bottom to vacuum vessel roof. Override by specifying together with z_min
    # /\ Only necessary if E_array, p_array, R_array and/or z_array are set to nothing
    # /\ Only necessary if !(useDistrFile)

    # EXTRA KEYWORD ARGUMENTS BELOW
    extra_kw_args = Dict(:limit_phi => true, :max_tries => 0)
    # limit_phi limits the number of toroidal turns for orbits
    # The orbit integration algorithm will try progressively smaller timesteps max_tries number of times
end

## -----------------------------------------------------------------------------
# Then you execute the script
include(folderpath_OWCF*"calcEpRzTopoMap.jl")

## -----------------------------------------------------------------------------
# Then you clean up after yourself
if batch_job && distributed
    for i in workers()
        rmprocs(i)
    end
end
