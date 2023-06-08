#########################################  ps2os.jl ##############################################################################

#### Description:
# This script transforms a fast-ion distribution function from particle-space coordinates
# (E,p,R,z) to orbit-space coordinates (E,pm,Rm). As of this version of the OWCF, it can only
# employ Monte-Carlo sampling to do so.
# 
# The fast-ion distribution data can be specified in either .jld2 or .h5/.hdf5 file format.
# If the fast-ion distribution is provided as an .h5 file, it must have the keys:
# 'f' or 'F_ps' - The 4D matrix containing the fast-ion distribution.
# 'energy' - The 1D array containing the energy grid points
# 'pitch ' - The 1D array containing the pitch grid points
# 'R' - The 1D array containing the R grid points
# 'z' - The 1D array containing the z grid points
#
# If the fast-ion distribution is provided as a .jld2 file, it must have the keys:
# 'F_ps' - The 4D matrix containing the fast-ion distribution. Order [energy,pitch,R,z] (forwards).
# 'energy' - The 1D array containing the energy grid points
# 'pitch ' - The 1D array containing the pitch grid points
# 'R' - The 1D array containing the R grid points
# 'z' - The 1D array containing the z grid points

#### Inputs (Units given when defined in script)
# Given via input file start_ps2os_template.jl, for example. 'template' should be replaced by whatever.

#### Outputs
# -

#### Saved files
# ps2os_results_[tokamak]_[TRANSP_id]_at[timepoint]s_[FI species]_[length(E_array)]x[length(pm_array)]x[length(Rm_array)].jld2
# It will have the keys
#   F_os_raw - The 1D vectorized form of the fast-ion distribution. However, 'raw' means you have to multiply with (nfast/sum(F_os_raw)) to acquire the correctly scaled vectorized distribution - Array{Int64,1}
#   nfast - The total number of fast ions in the distribution - Union{Float64, Int64}
#   numOsamples - The total number of samples of your Monte-Carlo sampling transform - Int64
#   E_array - The energy grid points of the (E,pm,Rm) grid onto which the distribution has been transformed. In keV - Array{Float64,1}
#   pm_array - The pm grid points of the (E,pm,Rm) grid onto which the distribution has been transformed - Array{Float64,1}
#   Rm_array - The Rm grid points of the (E,pm,Rm) grid onto which the distribution has been transformed. In meters - Array{Float64,1}
#   E_array_ps - The energy grid points of the original (E,p,R,z) grid of the fast-ion distribution. In keV. - Array{Float64,1}
#   p_array - The pitch grid points of the original (E,p,R,z) grid of the fast-ion distribution - Array{Float64,1}
#   R_array - The major radius grid points of the original (E,p,R,z) grid of the fast-ion distribtion. In centimeters or meters - Array{Float64,1}
#   z_array - The vertcal grid points of the original (E,p,R,z) grid of the fast-ion distribtion. In centimeters or meters - Array{Float64,1}
#   extra_kw_args - The extra keyword arguments provided by the user for the orbit integration algorithm - Dictionary
# If include1Dto3D is set to true in the input variables, there will also be a saved file called
# F_os_3D_results_[tokamak]_[TRANSP_id]_at[timepoint]s_[FI species]_[length(E_array)]x[length(pm_array)]x[length(Rm_array)].jld2
# It will have the keys
#   F_os_3D - The (E,pm,Rm) fast-ion distribution as a three-dimensional data matrix - Array{Float64,3}
#   E_array - The energy grid points of the (E,pm,Rm) grid onto which the distribution has been transformed. In keV - Array{Float64,1}
#   pm_array - The pm grid points of the (E,pm,Rm) grid onto which the distribution has been transformed - Array{Float64,1}
#   Rm_array - The Rm grid points of the (E,pm,Rm) grid onto which the distribution has been transformed. In meters - Array{Float64,1}
#   extra_kw_args - The extra keyword arguments provided by the user for the orbit integration algorithm - Dictionary
#
#### Other
# Please note that in future versions of the OWCF, the user may use a toggle input to select other transformation methods
# than simply Monte-Carlo methods. This has, however, not yet been incorporated into the OWCF.
#
# Script written by Henrik Järleblad. Last maintained 2022-08-26.
##################################################################################################################################

## --------------------------------------------------------------------------------------
# Loading packages
## --------------------------------------------------------------------------------------
verbose && println("Loading Julia packages... ")
@everywhere begin
    cd(folderpath_OWCF) # Necessary to move all the workers to the correct folder
    using EFIT # For equilibrium purposes
    using JLD2 # For saving variables to file. For back-up.
    using FileIO # -||-
    using Equilibrium # For magnetic equilibrium loading (read_geqdsk)
    using GuidingCenterOrbits # For orbit and Jacobian calculations
    using OrbitTomography # That's what this all is 'bout!
    using HDF5
    include("misc/species_func.jl") # To convert species labels to particle mass
    include("misc/availReacts.jl") # To examine fusion reaction and extract thermal and fast-ion species
    include("misc/rewriteReacts.jl") # To rewrite a fusion reaction from the A(b,c)D format to the A-b=c-D format
    include("extra/dependencies.jl")
end

## --------------------------------------------------------------------------------------
# Determining fast-ion distribution file extension (.h5 or .jld2)
## --------------------------------------------------------------------------------------
verbose && println("Determining distribution file extensions... ")
fileext_FI = (split(filepath_distr,"."))[end] # Assume last part after final '.' is the file extension
fileext_FI = lowercase(fileext_FI)

if !(fileext_FI=="h5" || fileext_FI=="jld2" || fileext_FI=="hdf5")
    println("Fast-ion distribution file format: ."*fileext_FI)
    error("Unknown fast-ion distribution file format. Please re-specify file and re-try.")
end

## --------------------------------------------------------------------------------------
# Loading tokamak equilibrium
## --------------------------------------------------------------------------------------
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

## --------------------------------------------------------------------------------------
verbose && println("Defining orbit-grid vectors... ")
if useWeightsFile
    #og = get_NES_og(filepath_W, filepath_equil, verbose=verbose)
    myfile = jldopen(filepath_W,false,false,false,IOStream)
    E_array = myfile["E_array"]
    pm_array = myfile["pm_array"]
    Rm_array = myfile["Rm_array"]
    if haskey(myfile,"reaction")
        reaction = myfile["reaction"]
    end
    if haskey(myfile,"reaction_full")
        reaction_full = myfile["reaction_full"]
    end
    close(myfile)
elseif useOrbGridFile
    verbose && println("Loading from .jld2 orbit grid file... ")
    myfile = jldopen(filepath_OG,false,false,false,IOStream)
    og = myfile["og"]
    E_array = og.energy
    pm_array = og.pitch
    Rm_array = og.r
    extra_kw_args = myfile["extra_kw_args"]
    if haskey(myfile,"FI_species")
        FI_species = myfile["FI_species"]
    end
    close(myfile)
else
    E_array = range(Emin, stop=Emax, length=nE)
    pm_array = range(pm_min, stop=pm_max, length=npm)
    if (Rm_min==nothing) || (Rm_max==nothing)
        if inclPrideRockOrbs
            # 4/5 of the distance from the HFS wall to the magnetic axis is usually enough to capture all the Pride Rock orbits
            Rm_array = range((4*M.axis[1]+minimum(wall.r))/5, stop=maximum(wall.r), length=nRm)
        else
            Rm_array = range(M.axis[1], stop=maximum(wall.r), length=nRm)
        end
    else
        Rm_array = range(Rm_min, stop=Rm_max, length=nRm)
    end
end

## ---------------------------------------------------------------------------------------------
# Determine fast-ion species from reaction
if (@isdefined reaction_full)
    thermal_species, FI_species = checkReaction(reaction_full)
    @everywhere FI_species = $FI_species # Transfer variable to all external processes
elseif (@isdefined reaction)
    FI_species = (split(reaction,"-"))[2] # Assume first species specified in reaction to be the fast-ion species. For example, in 'p-t' the 'p' will be assumed the thermal species.
    @everywhere FI_species = $FI_species # Transfer variable to all external processes
else
    FI_species = FI_species # Already defined
end

## ------
# Printing script info and inputs
println("")
println("----------------------------------------ps2os.jl----------------------------------------")
print("Tokamak specified: "*tokamak*"        ")
print("TRANSP_id specified: "*TRANSP_id*"       ")
println("Fast-ion species: "*FI_species)
println("")
println("Fast-ion distribution file specified: ")
println(filepath_distr)
println("Timepoint: "*timepoint)
println("")
println("Number of monte-carlo samples: $(numOsamples)")
println("Algorithm will save progress every $(nbatch) samples.")
println("")
println("Magnetic equilibrium file specified: ")
println("--- "*filepath_equil)
if useWeightsFile
    println("Weights file specified: ")
    println("--- "*filepath_W)
    println("Fast-ion distribution will be transformed onto same orbit-space grid as orbit weights ($(length(E_array))x$(length(pm_array))x$(length(Rm_array))).")
elseif useOrbGridFile
    println("Orbit-grid file specified: ")
    println("--- "*filepath_OG)
    println("Fast-ion distribution will be transformed onto same orbit-space grid ($(length(E_array))x$(length(pm_array))x$(length(Rm_array))) as in file.")
else
    println("Fast-ion distribution will be transformed onto a $(length(E_array))x$(length(pm_array))x$(length(Rm_array)) orbit-space grid with")
    println("--- Energy: [$(minimum(E_array)),$(maximum(E_array))] keV")
    println("--- Pitch maximum: [$(minimum(pm_array)),$(maximum(pm_array))]")
    println("--- Radius maximum: [$(minimum(Rm_array)),$(maximum(Rm_array))] meters")
end
println("")
if flip_F_ps_pitch
    println("Fast-ion distribution will be mirrored in pitch before being sampled from.")
    println("")
end
if interp
    println("Fast-ion distribution will be interpolated onto a $(nE_ps)x$(np_ps)x$(nR_ps)x$(nz_ps) (E,p,R,z) grid")
    println("before being samples from.")
    println("")
end
if distributed
    println("$(nprocs()) processors will be used for parallel computing when transforming from (E,p,R,z) to (E,pm,Rm).")
    println("")
else
    println("Single-threaded transforming from (E,p,R,z) to (E,pm,Rm).")
    println("")
end
if include1Dto3D
    println("Fast-ion (E,pm,Rm) distribution will be inflated and saved in its full 3D format.")
    println("")
end
println("If you would like to change any settings, please edit the start_ps2os_template.jl file (or equivalent)")
println("Written by Henrik Järleblad. Last maintained 2022-09-13.")
println("----------------------------------------------------------------------------------------")

## --------------------------------------------------------------------------------------
if !useOrbGridFile
    verbose && println("Calculating orbit grid... ")
    og_orbs, og = OrbitTomography.orbit_grid(M,E_array,pm_array,Rm_array; wall=wall, amu=getSpeciesAmu(FI_species), q=getSpeciesEcu(FI_species), extra_kw_args...)
    og_orbs = nothing # Minimize memory
end
## --------------------------------------------------------------------------------------
verbose && println("Loading fast-ion distribution data from file... ")
if fileext_FI=="h5" || fileext_FI=="hdf5"
    F_ps, E_array_ps, p_array, R_array, z_array = h5to4D(filepath_distr, verbose = verbose) # Load fast-ion distribution
else # Else, assume .jld2 file format
    myfile = jldopen(filepath_distr,false,false,false,IOStream)
    if haskey(myfile,"F_ps")
        F_ps = myfile["F_ps"]
    elseif haskey(myfile,"f")
        F_ps = myfile["f"]
    else
        error("Unrecognized .jld2 key for fast-ion distribution in filepath_distr. Please correct and re-try.")
    end
    E_array_ps = myfile["energy"]
    p_array = myfile["pitch"]
    R_array = myfile["R"]
    if haskey(myfile,"z")
        z_array = myfile["z"]
    else
        z_array = myfile["Z"]
    end
    close(myfile)
end

## --------------------------------------------------------------------------------------
if interp
    verbose && println("Fast-ion distribution will be interpolated onto requested (E,p,R,z) grid.")
    distr_dim = [nE_ps,np_ps,nR_ps,nz_ps]
else
    distr_dim = []
end

## --------------------------------------------------------------------------------------
verbose && println("Going into ps2os_streamlined (performance=true manually set)... ")
F_os_raw, nfast = ps2os_streamlined(F_ps, E_array_ps, p_array, R_array, z_array, filepath_equil, og; numOsamples=numOsamples, verbose=verbose, distr_dim=distr_dim, visualizeProgress=visualizeProgress, nbatch=nbatch, distributed=distributed, flip_F_ps_pitch=flip_F_ps_pitch, performance=true, GCP=getGCP(FI_species), extra_kw_args... )

## --------------------------------------------------------------------------------------
verbose && println("Saving results... ")
global filepath_tm_orig = folderpath_o*"ps2os_results_"*tokamak*"_"*TRANSP_id*"_at"*timepoint*"s_"*FI_species*"_$(length(E_array))x$(length(pm_array))x$(length(Rm_array))"
global filepath_tm = deepcopy(filepath_tm_orig)
global count = 1
while isfile(filepath_tm*".jld2") # To take care of not overwriting files. Add _(1), _(2) etc
    global filepath_tm = filepath_tm_orig*"_($(Int64(count)))"
    global count += 1 # global scope, to surpress warnings
end
global filepath_tm = filepath_tm*".jld2"
myfile = jldopen(filepath_tm,true,true,false,IOStream)
write(myfile,"F_os_raw",F_os_raw)
write(myfile,"nfast",nfast)
write(myfile,"numOsamples",numOsamples)
write(myfile,"E_array",E_array)
write(myfile,"pm_array",pm_array)
write(myfile,"Rm_array",Rm_array)
write(myfile,"E_array_ps",E_array_ps)
write(myfile,"p_array",p_array)
write(myfile,"R_array",R_array)
write(myfile,"z_array",z_array)
write(myfile,"extra_kw_args",extra_kw_args)
close(myfile)

if include1Dto3D
    F_os = (nfast/sum(F_os_raw)) .*F_os_raw
    F_os_3D = map_orbits(og, F_os)

    global filepath_output_orig = folderpath_o*"F_os_3D_"*tokamak*"_"*TRANSP_id*"_at"*timepoint*"s_"*FI_species*"_$(length(E_array))x$(length(pm_array))x$(length(Rm_array))"
    global filepath_output = deepcopy(filepath_output_orig)
    global count = 1
    while isfile(filepath_output*".jld2") # To take care of not overwriting files. Add _(1), _(2) etc
        global filepath_output = filepath_output_orig*"_($(Int64(count)))"
        global count += 1 # global scope, to surpress warnings
    end
    global filepath_output = filepath_output*".jld2"
    myfile = jldopen(filepath_output,true,true,false,IOStream)
    write(myfile, "F_os_3D", F_os_3D)
    write(myfile, "E_array", E_array)
    write(myfile, "pm_array", pm_array)
    write(myfile, "Rm_array", Rm_array)
    write(myfile, "extra_kw_args",extra_kw_args)
    close(myfile)
end

println("~~~~~~ps2os.jl is done!~~~~~~")
