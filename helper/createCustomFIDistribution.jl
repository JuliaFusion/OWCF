# This script will create a custom fast-ion (FI) distribution in (E,p,R,z) coordinates.
# It will be saved in a .jld2 file format with file keys compatible with the input requirements to e.g. calcSpec.jl.
# Currently, Gaussian and neutral beam injection (NBI) collision-physics (former NBI slowing-down) functions are supported.
# You could also choose whether the FI distribution should be created only in (E,p) as a 2D FI distribution.
# It will then be computed as constant in (R,z) for all (R,z) grid points. This is done via the constant_Rz input variable
# The NBI collision-physics function are computed using the formulas in W.G.F. Core, Nucl. Fusion 33, 829 (1993)
# and J.D. Gaffey, J. Plasma Physics 16, 149-169, (1976). In the current version of the OWCF, the collision-physics functions
# are created only in (E,p) as a 2D FI distribution, and then computed as constant in (R,z) for all (R,z) grid points.
# The filepath_equil file can be an .eqdsk file, or computed with the extra/createCustomMagneticEquilibrium.jl script.
# The filepath_thermal_profiles file can be a TRANSP .cdf file, "" (default OWCF temperature and density profiles will be used) or 
# an output file of the helper/createCustomThermalProfiles.jl script.
################################################################################

################################################################################
############################ Specify inputs ####################################
################################################################################
# PLEASE NOTE!
# You only need to specify the general inputs immidiately below and then the specific 
# inputs for the fast-ion distribution type that you need, i.e.
# if distribution_type == :gaussian, elseif distribution_type == :collisional etc.

# Start of general inputs
constant_Rz = true # Set to true, if the FI distribution should be computed as constant in (R,z) position space
distribution_type = :gaussian # Currently supported options include :gaussian, :collisional and :custom
folderpath_OWCF = "" # /path/to/the/directory/of/the/OWCF/
folderpath_out = "" # /path/to/where/you/want/your/output/file/to/be/saved.jld2
max_E = 250.0 # keV. The maximum of the energy grid points
max_p = 1.0 # -. The maximum of the pitch grid points
max_R = 3.9 # m. The maximum of the R grid points
max_z = 2.0 # m. The maximum of the z grid points
min_E = 1.0 # keV. The minimum of the energy grid points
min_p = -1.0 # -. The minimum of the pitch grid points
min_R = 1.8 # m. The minimum of the R grid points
min_z = -2.0 # m. The minimum of the z grid points
nE = 100 # Number of grid points in energy (E)
np = 50 # Number of grid points in pitch (p)
nR = 20 # Number of grid points in major radius (R). PLEASE NOTE! Even if constant_Rz is set to true, this needs to be specified!
nz = 20 # Number of grid points in vertical coordinate (z). PLEASE NOTE! Even if constant_Rz is set to true, this needs to be specified!
save_plots = true # If true, (E,p) and (R,z) plots of the computed FI distribution will be saved in .png file format. For post-analysis.
tot_N_FI = 1.0e19 # The total number of fast ions in the plasma. Obtained when integrating over all of (E,p,R,z) space
verbose = true # If true, the script will talk a lot!
# End of general inputs

# Start of specific inputs
if distribution_type == :gaussian # If you would like a gaussian FI distribution, please specify...
    floor_level = 1.0e-3 # Below floor_level*maximum(f_Gaussian), the FI distribution will be manually set to 0
    peak_E = 50.0 # keV. The energy (E) coordinate of the peak of the gaussian 
    peak_p = 0.5 # -. The pitch (E) coordinate of the peak of the gaussian 
    peak_R = 0.0 # m. The major radius (R) coordinate of the peak of the gaussian. If constant_Rz==true, this does not matter
    peak_z = 0.0 # m. The vertical coordinate (z) coordinate of the peak of the gaussian. If constant_Rz==true, this does not matter
    Sigma = [50.0, 0.05, 0.0, 0.0] # The values of Sigma for the Gaussian distribution (1/(4*pi^2)) * (det(Sigma)^(-1/2)) * exp(-0.5*transpose(X-X_peak)*(Sigma^(-1))*(X-X_peak))
elseif distribution_type == :collisional # If you would like a collision-physics FI distribution, please specify...
    assume_Ti_equalTo_Te = false # Should Ti=Te be assumed?
    assume_ni_equalTo_ne = false # Should ni=ne be assumed?
    constant_Rz = true # Currently, the collision-physics functions need to be assumed constant in Rz
    dampen = false # Set to true, if the collision-physics functions should be damped below the NBI injection energy
    FI_species = "" # The FI particle species. Valid list of species can be found in misc/species_func.jl
    filepath_equil = "" # For example "g94701_0-50.7932.eqdsk" or "solovev_equilibrium_2024-09-30.jld2"
    filepath_thermal_profiles = "" # /path/to/a/thermal/plasma/data/file. # For example "96100J01.cdf", "my_custom_thermal_profiles.jld2" or ""
    inj_E = 0.0 # keV. The injection energy (E) coordinate of the collision-physics function
    inj_p = 0.0 # -. The injection pitch (p) coordinate of the collision-physics function
    inj_R = 0.0 # m. The injection major radius (R) coordinate of the collision-physics function
    inj_z = 0.0 # m. The injection vertical (z) coordinate of the collision-physics function
    R_of_temp_n_dens = inj_R # By default, the inj_R input variable is used to determine the vertical (z) coordinate of the thermal temperature+density profiles
    z_of_temp_n_dens = inj_z # By default, the inj_z input variable is used to determine the vertical (z) coordinate of the thermal temperature+density profiles
    thermal_species = "" # The thermal ion species of interest for the temp+dens profiles. Valid list of species can be found in misc/species_func.jl
    # /\ For the thermal_species input variable, only one thermal species is currently supported.

    # ---Conditional input variables---
    # These input variables might need to be specified, depending on how the input variables above are specified.
    # Please check how you specified the input variables above, and specify the conditional input variables below accordingly.
    if (!assume_Ti_equalTo_Te || !assume_ni_equalTo_ne) && !(lowercase(filepath_thermal_profiles[end-2:end])=="cdf") # If Ti==Te is not true, or if ni=ne is not true, and the thermal profile file is not a TRANSP .cdf file, please specify...
        filepath_thermal_electron_profiles = "" # # /path/to/a/thermal/electron/data/file. Should be an output of the createCustomThermalProfiles.jl
    end
    if dampen # If dampen is set to true, please specify... 
        damp_type = :erfc # The type of damping. Can be either :erfc (complementory error-function) or :linear
        E_tail_length = nothing # The length (in keV) of the tail of the collision-physics functions below the NBI injection energy.
        # /\ If left as nothing, automatic damping below the v_L speed (W.G.F. Core, Nucl. Fusion 33, 829 (1993)) will be used
    end
    if lowercase(filepath_thermal_profiles[end-2:end])=="cdf" # If the thermal profiles data is a TRANSP file, please specify...
        timepoint = 00.000 # seconds. The timepoint of interest for the temperature+density profiles to be extracted from the TRANSP .cdf format file
    elseif filepath_thermal_profiles=="" # If the thermal profiles data is left unspecified, please specify...
        dens_on_axis = 0.0 # m^-3. The thermal density value on-axis
        temp_on_axis = 0.0 # keV. The thermal temperature value on-axis
    elseif lowercase(filepath_thermal_profiles[end-3:end])=="jld2"
        # Nothing extra to specify, if the thermal profiles data is an output file of the helper/createCustomThermalProfiles.jl script
    else
        error("Invalid thermal profile input file format! Please correct and re-try.")
    end
elseif distribution_type==:custom
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
else
    error("distribution_type input variable not correctly specified (:gaussian, :collisional, :custom). Please correct and re-try.")
end
# End of specific inputs
################################################################################
############################ Change to OWCF directory ##########################
############################ Load Julia packages ###############################
################################################################################
verbose && println("Loading Julia packages... ")
cd(folderpath_OWCF)
using Pkg
Pkg.activate(".")
using JLD2
using Statistics
using LinearAlgebra
using Interpolations
using Dates
include(folderpath_OWCF*"misc/temp_n_dens.jl")
if save_plots
    using Plots
end
if distribution_type==:collisional
    include(folderpath_OWCF*"extra/dependencies.jl")
end
################################################################################
############################ Main script #######################################
################################################################################
# Create the (equidistant) (E,p) grid
E_array = collect(range(min_E,stop=max_E,length=nE)); dE = diff(E_array)[1]
p_array = collect(range(min_p,stop=max_p,length=np)); dp = diff(p_array)[1]
verbose && println("Created (E,p) grid points. dE=$(dE) keV. dp=$(dp).")
# Create the (R,z) grid
R_array = collect(range(min_R,stop=max_R,length=nR)); dR = diff(R_array)[1]
z_array = collect(range(min_z,stop=max_z,length=nz)); dz = diff(z_array)[1]
verbose && println("Created (R,z) grid points. dR=$(dR) m. dz=$(dz) m.")

# Two different approaches, depending on if constant_Rz
if constant_Rz
    verbose && println("Assuming FI distribution to be constant in (R,z) position space.")
    if distribution_type==:gaussian
        verbose && println("Creating Gaussian distribution function... ")
        Sigmam = diagm(Sigma[1:2]) # Only the correlation lengths in (E,p) matter, because of constant_Rz. diagm(x) creates a matrix with x as diagonal 
        X_peak = [peak_E, peak_p] # The mean (peak) of the Gaussian
        F_Ep = zeros(length(E_array),length(p_array)) # The fast-ion distribution
        for (iE,E) in enumerate(E_array), (ip,p) in enumerate(p_array)
            X = [E,p]
            F_Ep[iE,ip] = inv(2*pi) * (det(Sigmam)^(-1/2)) * exp(-0.5*transpose(X .- X_peak)*inv(Sigmam)*(X .- X_peak)) # The 2D Gaussian formula
        end
        verbose && println("Setting all values below $(round(floor_level*maximum(F_Ep),sigdigits=4)) to 0 in the Gaussian FI distribution... ")
        F_Ep = map(x-> x<floor_level*maximum(F_Ep) ? 0.0 : x, F_Ep) # Set everything below floor_level*maximum(F_Ep) to 0.0, to avoid extremely small, unstable values.
    elseif distribution_type==:collisional
        verbose && println("Creating NBI collision-physics distribution function... ")
        if ((split(filepath_equil,"."))[end] == "eqdsk") || ((split(filepath_equil,"."))[end] == "geqdsk")
            M, wall = read_geqdsk(filepath_equil,clockwise_phi=false) # Assume counter-clockwise phi-direction
        else # Otherwise, assume magnetic equilibrium is a saved .jld2 file
            myfile = jldopen(filepath_equil,false,false,false,IOStream)
            M, wall = myfile["S"], myfile["wall"]
            close(myfile)
        end
        verbose && println("Successfully loaded magnetic equilibrium.")
        ψ_rz = M(R_of_temp_n_dens,z_of_temp_n_dens)
        psi_on_axis, psi_at_bdry = psi_limits(M)
        rho_of_temp_n_dens = sqrt((ψ_rz-psi_on_axis)/(psi_at_bdry-psi_on_axis)) # The formula for the normalized flux coordinate ρ_pol = sqrt((ψ-ψ_axis)/(ψ_edge-ψ_axis))
        verbose && println("Successfully computed rho_pol for the (R,z) point of temperature and density.")
        if lowercase(filepath_thermal_profiles[end-2:end])=="cdf" # If the thermal profiles data is a TRANSP file...
            T_i_vec = [(getTempProfileFromTRANSP(timepoint,filepath_thermal_profiles,thermal_species))(rho_of_temp_n_dens)] # The thermal ion temperature at the (R,z) point and timepoint of interest, extracted from a TRANSP .cdf file, and put into a single-element vector
            n_i_vec = [(getDensProfileFromTRANSP(timepoint,filepath_thermal_profiles,thermal_species))(rho_of_temp_n_dens)] # The thermal ion density at the (R,z) point and timepoint of interest, extracted from a TRANSP .cdf file, and put into a single-element vector
            n_e = n_i_vec[1] # Initially assume that n_e == n_i
            verbose && println("Successfully computed Ti, ni and ne from TRANSP .cdf file.")
        elseif filepath_thermal_profiles=="" # If the thermal profiles data is left unspecified...
            T_i_vec = [getAnalyticalTemp(temp_on_axis,rho_of_temp_n_dens)]
            n_i_vec = [getAnalyticalDens(dens_on_axis,rho_of_temp_n_dens)]
            n_e = n_i_vec[1] # Initially assume that n_e == n_i
            verbose && println("Successfully computed Ti, ni and ne from OWCF default formulas.")
        elseif lowercase(filepath_thermal_profiles[end-3:end])=="jld2" # If the thermal profiles data is an output file of the helper/createCustomThermalProfiles.jl script
            myfile = jldopen(filepath_thermal_profiles,false,false,false,IOStream)
            rho_array_from_file = myfile["rho_pol"]
            temp_array_from_file = myfile["thermal_temp"]
            dens_array_from_file = myfile["thermal_dens"]
            close(myfile)
            temp_itp = Interpolations.interpolate((rho_array_from_file,), temp_array_from_file, Gridded(Linear()))
            dens_itp = Interpolations.interpolate((rho_array_from_file,), dens_array_from_file, Gridded(Linear()))
            temp_etp = Interpolations.extrapolate(temp_itp,0)
            dens_etp = Interpolations.extrapolate(dens_itp,0)
            T_i_vec = [temp_etp(rho_of_temp_n_dens)]
            n_i_vec = [dens_etp(rho_of_temp_n_dens)]
            n_e = n_i_vec[1] # Initially assume that n_e == n_i
            verbose && println("Successfully computed Ti, ni and ne from OWCF .jld2 file.")
        else
            error("Invalid thermal profile input file format! Please correct and re-try.")
        end

        T_e = T_i_vec[1] # By default, assume Te==Ti
        if !assume_Ti_equalTo_Te # If Te==Ti should not be assumed...
            verbose && println("T_e=T_i will not be assumed.")
            if lowercase(filepath_thermal_profiles[end-2:end])=="cdf" # If the thermal profiles data is a TRANSP file...
                T_e = (getTempProfileFromTRANSP(timepoint,filepath_thermal_profiles,"e"))(rho_of_temp_n_dens)
                verbose && println("---> Successfully computed T_e from TRANSP .cdf file.")
            else
                myfile = jldopen(filepath_thermal_electron_profiles,false,false,false,IOStream)
                rho_e_array_from_file = myfile["rho_pol"]
                temp_e_array_from_file = myfile["thermal_temp"]
                close(myfile)
                temp_e_itp = Interpolations.interpolate((rho_e_array_from_file,), temp_e_array_from_file, Gridded(Linear()))
                temp_e_etp = Interpolations.extrapolate(temp_e_itp,0)
                T_e = temp_e_etp(rho_of_temp_n_dens)
                verbose && println("---> Successfully computed T_e from OWCF .jld2 file.")
            end
        end

        n_e = n_e # By default, assume that n_e==n_i
        if !assume_ni_equalTo_ne # If ne==ni should not be assumed...
            verbose && println("n_e=n_i will not be assumed.")
            if lowercase(filepath_thermal_profiles[end-2:end])=="cdf" # If the thermal profiles data is a TRANSP file...
                n_e = (getDensProfileFromTRANSP(timepoint,filepath_thermal_profiles,"e"))(rho_of_temp_n_dens)
                verbose && println("---> Successfully computed n_e from TRANSP .cdf file.")
            else
                myfile = jldopen(filepath_thermal_electron_profiles,false,false,false,IOStream)
                rho_e_array_from_file = myfile["rho_pol"]
                dens_e_array_from_file = myfile["thermal_dens"]
                close(myfile)
                dens_e_itp = Interpolations.interpolate((rho_e_array_from_file,), dens_e_array_from_file, Gridded(Linear()))
                dens_e_etp = Interpolations.extrapolate(dens_e_itp,0)
                n_e = dens_e_etp(rho_of_temp_n_dens)
                verbose && println("---> Successfully computed n_e from OWCF .jld2 file.")
            end
        end

        F_Ep = slowing_down_function(inj_E, inj_p, E_array, p_array, n_e, T_e, FI_species, [thermal_species], n_i_vec, T_i_vec; E_tail_length=E_tail_length, dampen=dampen, damp_type=damp_type)
    elseif distribution_type==:custom
        verbose && println("Creating custom distribution function... ")
        # WRITE EXTRA CODE HERE FOR THE :custom FI DISTRIBUTION FUNCTION, IF NEEDED
        # WRITE EXTRA CODE HERE FOR THE :custom FI DISTRIBUTION FUNCTION, IF NEEDED
        # WRITE EXTRA CODE HERE FOR THE :custom FI DISTRIBUTION FUNCTION, IF NEEDED
        F_Ep = zeros(length(E_array),length(p_array)) # The fast-ion distribution
        for (iE,E) in enumerate(E_array), (ip,p) in enumerate(p_array)
            F_Ep[iE,ip] = custom_FI_Ep_distribution(E,p)
        end
    else
        error("distribution_type input variable not correctly specified (:gaussian, :collisional). Please correct and re-try.")
    end
    verbose && println("Success! Normalizing FI distribution... ")
    F_Ep = tot_N_FI .*F_Ep ./ sum((dE*dp) .*F_Ep) # Re-normalize so that F_Ep integrate to tot_N_FI
    F_EpRz = repeat(F_Ep,1,1,nR,nz) # Make F_EpRz equal to F_Ep for all (R,z) grid points
else # If !constant_Rz
    if distribution_type==:gaussian
        verbose && println("Creating (E,p,R,z) Gaussian FI distribution... ")
        Sigmam = diagm(Sigma) # Only the correlation lengths in (E,p) matter, because of constant_Rz. diagm(x) creates a matrix with x as diagonal 
        X_peak = [peak_E, peak_p, peak_R, peak_z] # The mean (peak) of the Gaussian
        F_EpRz = zeros(length(E_array),length(p_array),length(R_array),length(z_array)) # The pre-allocated fast-ion distribution
        progress_length = reduce(*,[length(E_array),length(p_array),length(R_array),length(z_array)]) # The product of all (E,p,R,z) vector lengths
        global progress # Declare global scope
        progress = 1
        for (iE,E) in enumerate(E_array), (ip,p) in enumerate(p_array), (iR,R) in enumerate(R_array), (iz,z) in enumerate(z_array)
            X = [E,p,R,z]
            F_EpRz[iE,ip,iR,iz] = (1/(4*pi^2)) * (det(Sigmam)^(-1/2)) * exp(-0.5*transpose(X .- X_peak)*inv(Sigmam)*(X .- X_peak)) # The 4D Gaussian formula
            verbose && println("Creating (E,p,R,z) Gaussian FI distribution: $(round(100*progress/progress_length,sigdigits=4)) %... ")
            progress += 1
        end
        verbose && println("Setting all values below $(round(floor_level*maximum(F_EpRz),sigdigits=4)) to 0 in the (E,p,R,z) Gaussian FI distribution... ")
        F_EpRz = map(x-> x<floor_level*maximum(F_EpRz) ? 0.0 : x, F_EpRz) # Set everything below floor_level*maximum(F_EpRz) to 0.0, to avoid extremely small, unstable values.
    elseif distribution_type==:collisional
        error("!constant_Rz and distribution_type==:collisional not yet supported. Please correct and re-try.")
    elseif distribution_type==:custom
        verbose && println("Creating custom (E,p,R,z) FI distribution function... ")
        # WRITE EXTRA CODE HERE FOR THE :custom FI DISTRIBUTION FUNCTION, IF NEEDED
        # WRITE EXTRA CODE HERE FOR THE :custom FI DISTRIBUTION FUNCTION, IF NEEDED
        # WRITE EXTRA CODE HERE FOR THE :custom FI DISTRIBUTION FUNCTION, IF NEEDED
        F_EpRz = zeros(length(E_array),length(p_array),length(R_array),length(z_array)) # The pre-allocated fast-ion distribution
        progress_length = reduce(*,[length(E_array),length(p_array),length(R_array),length(z_array)]) # The product of all (E,p,R,z) vector lengths
        progress = 1
        for (iE,E) in enumerate(E_array), (ip,p) in enumerate(p_array), (iR,R) in enumerate(R_array), (iz,z) in enumerate(z_array)
            F_EpRz[iE,ip,iR,iz] = custom_FI_distribution(E,p,R,z)
        end
    else
        error("distribution_type input variable not correctly specified (:gaussian, :collisional, :custom). Please correct and re-try.")
    end
    verbose && println("Success! Normalizing FI distribution... ")
    F_EpRz = tot_N_FI .*F_EpRz ./ sum((2*pi*dE*dp*dR*dz) .*reshape(R_array,(1,1,nR,1)) .*F_EpRz) # Re-normalize so that F_EpRz integrate to tot_N_FI
end

# The date and time. For creating unique output file names
date_and_time = split("$(Dates.now())","T")[1]*"at"*split("$(Dates.now())","T")[2][1:5]

if save_plots # If plots should be saved of the computed distribution function...
    if constant_Rz
        # If the FI distribution is constant in (R,z) position space, it does not matter which (R,z) point we examine for plotting
        myplt = Plots.heatmap(E_array,p_array,F_EpRz[:,:,1,1]',xlabel="Energy [keV]",ylabel="Pitch [-]",title="f(E,p) [keV^-1]",dpi=600)
        png(myplt,folderpath_out*"createCustomFIDistribution_"*date_and_time)
    else
        F_Ep = dropdims(sum((2*pi*dR*dz) .* reshape(R_array,(1,1,length(R_array),1)) .*F_EpRz,dims=(3,4)),dims=(3,4))
        F_Rz = dropdims(sum((dE*dp) .*F_EpRz,dims=(3,4)),dims=(3,4))

        myplt_Ep = Plots.heatmap(E_array,p_array,F_Ep',xlabel="Energy [keV]",ylabel="Pitch [-]",title="f(E,p) [keV^-1]",dpi=600)
        myplt_Rz = Plots.heatmap(R_array,z_array,F_Rz',xlabel="R [m]",ylabel="z [m]",title="f(R,z) [m^-3]",aspect_ratio=:equal,dpi=600)

        png(myplt_Ep,folderpath_out*"createCustomFIDistribution_Ep"*date_and_time)
        png(myplt_Rz,folderpath_out*"createCustomFIDistribution_Rz"*date_and_time)
    end
end

# Save results
myfile = jldopen(folderpath_out*"FI_distr_$(distribution_type)_$(nE)x$(np)x$(nR)x$(nz)_"*date_and_time*".jld2",true,true,false,IOStream)
write(myfile,"F_EpRz",F_EpRz)
write(myfile,"E_array",E_array)
write(myfile,"p_array",p_array)
write(myfile,"constant_Rz",constant_Rz)
write(myfile,"R_array",R_array)
write(myfile,"z_array",z_array)
write(myfile,"distribution_type","$(distribution_type)")
write(myfile,"ntot",tot_N_FI)
if distribution_type==:gaussian
    write(myfile,"floor_level",floor_level)
    write(myfile,"Sigma",Sigma)
elseif distribution_type==:collisional
    write(myfile,"assume_Ti_equalTo_Te",assume_Ti_equalTo_Te)
    write(myfile,"assume_ni_equalTo_ne",assume_ni_equalTo_ne)
    write(myfile,"dampen",dampen)
    write(myfile,"FI_species",FI_species)
    write(myfile,"filepath_equil",filepath_equil)
    write(myfile,"filepath_thermal_profiles",filepath_thermal_profiles)
    write(myfile,"R_of_temp_n_dens",R_of_temp_n_dens)
    write(myfile,"z_of_temp_n_dens",z_of_temp_n_dens)
    write(myfile,"thermal_species",thermal_species)
    if (!assume_Ti_equalTo_Te || !assume_ni_equalTo_ne) && !(lowercase(filepath_thermal_profiles[end-2:end])=="cdf")
        write(myfile,"filepath_thermal_electron_profiles",filepath_thermal_electron_profiles)
    end
    if dampen # If dampen is set to true, please specify... 
        write(myfile,"damp_type",damp_type)
        write(myfile,"E_tail_length",E_tail_length)
    end
    if lowercase(filepath_thermal_profiles[end-2:end])=="cdf" # If the thermal profiles data is a TRANSP file, please specify...
        write(myfile,"timepoint",timepoint)
    elseif filepath_thermal_profiles=="" # If the thermal profiles data is left unspecified, please specify...
        write(myfile,"dens_on_axis",dens_on_axis)
        write(myfile,"temp_on_axis",temp_on_axis)
    else
    end
else
end
close(myfile)
println("~~~createCustomFIDistribution.jl completed successfully!~~~")