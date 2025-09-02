######################################################## createCustomFIDistr.jl #################################################################

#### Description:
# This script will create a custom fast-ion (FI) distribution in (E,p,R,z) coordinates.
# It will be saved in a .jld2 file format with file keys compatible with the input requirements to e.g. OWCF/calcSpec.jl.
# Currently, Gaussian and neutral beam injection (NBI) collision-physics (former NBI slowing-down) functions are supported.
# You could also choose whether the FI distribution should be created only in (E,p) as a 2D FI distribution.
# It will then be computed as constant in (R,z) for all (R,z) grid points. Please see the start file example at 
# OWCF/templates/start_createCustomFIDistr_template.jl for more inputs information.
#
# The NBI collision-physics function are computed using the formulas in W.G.F. Core, Nucl. Fusion 33, 829 (1993)
# and J.D. Gaffey, J. Plasma Physics 16, 149-169, (1976). In the current version of the OWCF, the collision-physics functions
# are created only in (E,p) as a 2D FI distribution, and then computed as constant in (R,z) for all (R,z) grid points.
# The filepath_equil file can be an .eqdsk file, or computed with the extra/createCustomMagneticEquilibrium.jl script.
# The filepath_thermal_profiles file can be a TRANSP .cdf file, "" (default OWCF temperature and density profiles will be used) or 
# an output file of the helper/createCustomThermalProfiles.jl script.

#### Input variables:
# The input variables are defined in e.g. the start_createCustomFIDistr_template.jl file. Please see 
# OWCF/templates/start_createCustomFIDistr_template.jl for explanation.
#
# You run the createCustomFIDistr.jl script by making a copy of the start_createCustomFIDistr_template.jl 
# file, moving it to the OWCF/ folder, specifying the inputs and executing it by executing 
# 'julia start_createCustomFIDistr_template.jl' or equivalent.

#### Outputs:
# A .jld2 file containing the created fast-ion distribution and its abscissas. If the 'filename_out' input variable has been left 
# unspecified ("") in the start file, the default output data file name format will be used: 
#   createCustomFIDistr_[distribution type]_[nE]x[np]x[nR]x[nz]_[date and time].jld2. 
# The keys of the output file are as follows:
#   - F_EpRz - The 4D array containing the fast-ion distribution data. The element with index (iE,ip,iR,iz) is the value of the 
#              fast-ion distribution at the phase-space point (E_array[iE],p_array[ip],R_array[iR],z_array[iz]). The units of the 
#              fast-ion distribution is keV^-1_m^-3. The total number of fast ions of F_EpRz is obtained by integrating 
#              over (E,p,R,z) space with the 2*pi*R Jacobian included - Array{Float64}
#   - E_array - The grid points of the energy abscissa of the fast-ion distribution - Vector{Float64}
#   - E_array_units - The units of the grid points in E_array. As of the current OWCF version, this will always be "keV" - String
#   - p_array - The grid points of the pitch (v_||/v) abscissa of the fast-ion distribution - Vector{Float64}
#   - p_array_units - The units of the grid points of p_array. As of the current OWCF version, this will always be "-" - String
#   - R_array - The grid points of the R (major radius) abscissa of the fast-ion distribution - Vector{Float64}
#   - R_array_units - The units of the grid points of R_array. As of the current OWCF version, this will always be "m" - String
#   - z_array - The grid points of the z (vertical coordinate) abscissa of the fast-ion distribution - Vector{Float64}
#   - z_array_units - The units of the grid points of z_array. As of the current OWCF version, this will always be "m" - String
#   - constant_Rz - If true, the fast-ion distribution is constant in (R,z) space - Bool
#   - distribution_type - The fast-ion distribution type. Currently supported options are :gaussian, :collisional and :custom - Symbol
#   - ntot - The total number of fast ions in the fast-ion distribution. It is F_EpRz integrated over all of (E,p,R,z) space - Float64
# Depending on the value of the 'distribution_type' input variable, different additional keys will be included in the output file.
# These additional keys are as follows:
# If distribution_type==:gaussian
#   - floor_level - If a value of the fast-ion distribution is below floor_level*maximum(f_Gaussian), the FI distribution will be 
#                   manually set to 0 for that index - Float64
#   - X_peak - The (E,p,R,z) coordinate of the peak of the Gaussian distribution - Array{Float64}
#   - Sigma - The values of the correlation Sigma for the Gaussian distribution i.e.
#             (1/(4*pi^2)) * (det(Sigma)^(-1/2)) * exp(-0.5*transpose(X-X_peak)*(Sigma^(-1))*(X-X_peak)) - Array{Float64}
# If distribution_type==:collisional
#   - assume_Ti_equalTo_Te - A bool to indicate whether it was assumed that the ion and electron temperatures were equal or not - Bool
#   - assume_ni_equalTo_ne - A bool to indicate whether it was assumed that the ion and electron densities were equal or not - Bool
#   - dampen - A bool to indicate whether the slowing-down distribution functions are dampened below the critical energy or not - Bool
#   - FI_species - The particle species of the fast ions of the fast-ion distribution. Please see OWCF/misc/species_func.jl for more - String
#   - filepath_equil - The file path to the magnetic equilibrium data used to compute the slowing-down function - String
#   - filepath_thermal_profiles - The file path to the thermal ion temperature and density profiles used to -||-  String
#   - R_of_temp_n_dens - The major radius (R) coordinate of the thermal ion temperature and density values - Float64
#   - z_of_temp_n_dens - The vertical (z) cooridnate of the -||- - Float64
#   - thermal_species - The thermal ion species of the plasma. Please see OWCF/misc/species_func.jl for more info - String
#   - filepath_thermal_electron_profiles - For some inputs, the file path to the (thermal) electron temp/dens profiles is included - String
#   - damp_type - If 'dampen' was set to true, the type of damping is included in the output file data - Symbol
#   - E_tail_length - If 'dampen' was set to true, the length of the energy tail below the critical energy is included. In keV - Float64
#   - timepoint - For some inputs, the timepoint of the tokamak discharge is included - Float64
#   - dens_on_axis - If no thermal ion profiles were specified, the thermal ion density on-axis needs to be specified. In m^-3 - Float64
#   - temp_on_axis - If no thermal ion profiles were specified, the thermal ion temperature on-axis needs to be specified. In keV - Float64

#### Other:
#

# Script written by Henrik Järleblad. Last maintained 2025-07-10.
#################################################################################################################################################

################################################################################
############################ Change to OWCF directory ##########################
############################ Load Julia packages ###############################
################################################################################
verbose && println("Loading Julia packages... ")
@everywhere begin
    using JLD2
    using Statistics
    using LinearAlgebra
    using Interpolations
    using Dates
    using ProgressMeter
    include(folderpath_OWCF*"misc/temp_n_dens.jl")
    if save_plots
        using Plots
    end
    if distribution_type==:collisional
        include(folderpath_OWCF*"extra/dependencies.jl")
    end
end
################################################################################
############################ Main script #######################################
################################################################################

if distribution_type==:collisional && !constant_Rz
    @warn "The 'distribution_type' input variable was set to $(distribution_type), but the 'constant_Rz' input variable was set to $(constant_Rz). 
    As of the current OWCF version, the 'constant_Rz' variable needs to be set to true if the 'distribution_type' input variable is set to :collisional.
    This will therefore now be enforced."
    constant_Rz = true
end

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
        if filepath_thermal_profiles==""
            filepath_thermal_profiles="nothing" # Change to "nothing" for code convenience reasons
        end
        M, wall = nothing, nothing # Initialize equilibrium variables
        try
            global M; global wall
            M, wall = read_geqdsk(filepath_equil,clockwise_phi=false) # Assume counter-clockwise phi-direction
        catch # Otherwise, assume magnetic equilibrium is a saved .jld2 file
            global M; global wall; local myfile
            myfile = jldopen(filepath_equil,false,false,false,IOStream)
            M, wall = myfile["S"], myfile["wall"]
            close(myfile)
        end
        verbose && println("Successfully loaded magnetic equilibrium.")
        ψ_rz = M(R_of_temp_n_dens,z_of_temp_n_dens)
        psi_on_axis, psi_at_bdry = psi_limits(M)
        rho_of_temp_n_dens = sqrt((ψ_rz-psi_on_axis)/(psi_at_bdry-psi_on_axis)) # The formula for the normalized flux coordinate ρ_pol = sqrt((ψ-ψ_axis)/(ψ_edge-ψ_axis))
        verbose && println("Successfully computed rho_pol for the (R,z) point of temperature and density.")
        if filepath_thermal_profiles=="nothing" # If the thermal profiles data is left unspecified...
            T_i_vec = [getAnalyticalTemp(temp_on_axis,rho_of_temp_n_dens)]
            n_i_vec = [getAnalyticalDens(dens_on_axis,rho_of_temp_n_dens)]
            n_e = n_i_vec[1] # Initially assume that n_e == n_i
            verbose && println("Successfully computed Ti, ni and ne from OWCF default formulas.")
        elseif lowercase(filepath_thermal_profiles[end-2:end])=="cdf" # If the thermal profiles data is a TRANSP file...
            if typeof(timepoint)==String
                tp = parse(Float64,replace(timepoint, "," => "."))
            end
            T_i_vec = [(getTempProfileFromTRANSP(tp,filepath_thermal_profiles,thermal_species))(rho_of_temp_n_dens)] # The thermal ion temperature at the (R,z) point and timepoint of interest, extracted from a TRANSP .cdf file, and put into a single-element vector
            n_i_vec = [(getDensProfileFromTRANSP(tp,filepath_thermal_profiles,thermal_species))(rho_of_temp_n_dens)] # The thermal ion density at the (R,z) point and timepoint of interest, extracted from a TRANSP .cdf file, and put into a single-element vector
            n_e = n_i_vec[1] # Initially assume that n_e == n_i
            verbose && println("Successfully computed Ti, ni and ne from TRANSP .cdf file.")
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
                if typeof(timepoint)==String
                    tp = parse(Float64,replace(timepoint, "," => "."))
                end
                T_e = (getTempProfileFromTRANSP(tp,filepath_thermal_profiles,"e"))(rho_of_temp_n_dens)
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
                if typeof(timepoint)==String
                    tp = parse(Float64,replace(timepoint, "," => "."))
                end
                n_e = (getDensProfileFromTRANSP(tp,filepath_thermal_profiles,"e"))(rho_of_temp_n_dens)
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
        F_max = (1/(4*pi^2)) * (det(Sigmam)^(-1/2)) # We already now the maximum value of f(E,p,R,z)
        F_min = floor_level*F_max # The numerical threshold. All f(E,p,R,z) values below this value will be set to 0
        F_EpRz = zeros(length(E_array),length(p_array),length(R_array),length(z_array)) # The pre-allocated fast-ion distribution
        progress_length = length(F_EpRz) # Total number of (E,p,R,z) grid points

        num_o_good_points = 0 # Keep track of the number of (E,p,R,z) where f(E,p,R,z)>F_min
        prog_proc = Progress(progress_length; dt=0.5, desc="Creating (E,p,R,z) Gaussian FI distribution...", color=:blue) # Progress bar
        generate_showvalues(E, p, R, z, f, gpp) = () -> [("E [keV]",Int64(round(E))), ("p [-]",round(p,digits=2)), ("R [m]",round(R,digits=1)), ("z [m]",round(z,digits=1)), ("f(E,p,R,z) [keV^-3 m^-3]",round(f,sigdigits=3)), ("Points above threshold [%]",round(gpp,digits=2))] # Progress bar function
        for (iE,E) in enumerate(E_array), (ip,p) in enumerate(p_array), (iR,R) in enumerate(R_array), (iz,z) in enumerate(z_array)
            global num_o_good_points 
            X = [E,p,R,z]
            F = F_max * exp(-0.5*transpose(X .- X_peak)*inv(Sigmam)*(X .- X_peak)) # The 4D Gaussian formula
            above_threshold = F > F_min # Is f(E,p,R,z) above the threshold?
            above_threshold && (F_EpRz[iE,ip,iR,iz] = F) # If so, replace 0.0 with f(E,p,R,z) value
            above_threshold && (num_o_good_points += 1)
            ProgressMeter.next!(prog_proc; showvalues = generate_showvalues(E, p, R, z, F, 100*num_o_good_points/progress_length))
        end
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
            global progress
            verbose && println("Fast-ion distribution creation progress: $(round(100*progress/progress_length,digits=2))... ")
            F_EpRz[iE,ip,iR,iz] = custom_FI_distribution(E,p,R,z)
            progress +=1
        end
    else
        error("distribution_type input variable not correctly specified (:gaussian, :collisional, :custom). Please correct and re-try.")
    end
    verbose && println("Success! Re-normalizing FI distribution to have $(tot_N_FI)=∫ f(E,p,R,z) 2*pi*R dEdpdRdz... ")
    F_EpRz = tot_N_FI .*F_EpRz ./ sum((2*pi*dE*dp*dR*dz) .*reshape(R_array,(1,1,nR,1)) .*F_EpRz) # Re-normalize so that F_EpRz integrate to tot_N_FI
end

# The date and time. For creating unique output file names
date_and_time = split("$(Dates.now())","T")[1]*"at"*split("$(Dates.now())","T")[2][1:5]

# Create the output data .jld2 file name, given inputs or distribution type and date+time
if !(filename_out=="") # If the 'filename_out' input variable has been specified... 
    filepath_output_orig = folderpath_out*filename_out # Use the 'filename_out' input variable to name the output data file
else # Otherwise, if the 'filename_out' input variable was left unspecified (default), use the default file name format
    filepath_output_orig = folderpath_out*"createCustomFIDistr_$(distribution_type)_$(nE)x$(np)x$(nR)x$(nz)_"*date_and_time
end
filepath_output = deepcopy(filepath_output_orig)
count = 1
while isfile(filepath_output*".jld2") # To take care of not overwriting files. Add _(1), _(2) etc
    global count; global filepath_output
    filepath_output = filepath_output_orig*"_($(Int64(count)))"
    count += 1 # global scope, to surpress warnings
end

if save_plots # If plots should be saved of the computed distribution function...
    if constant_Rz
        # If the FI distribution is constant in (R,z) position space, it does not matter which (R,z) point we examine for plotting
        myplt = Plots.heatmap(E_array,p_array,F_EpRz[:,:,1,1]',xlabel="Energy [keV]",ylabel="Pitch [-]",title="f(E,p) [keV^-1]",dpi=600)
        png(myplt,filepath_output*"_Ep_plot")
    else
        F_Ep = dropdims(sum((2*pi*dR*dz) .* reshape(R_array,(1,1,length(R_array),1)) .*F_EpRz,dims=(3,4)),dims=(3,4))
        F_Rz = dropdims(sum((dE*dp) .*F_EpRz,dims=(1,2)),dims=(1,2))

        myplt_Ep = Plots.heatmap(E_array,p_array,F_Ep',xlabel="Energy [keV]",ylabel="Pitch [-]",title="f(E,p) [keV^-1]",dpi=600, xlims=extrema(E_array), ylims=extrema(p_array))
        myplt_Rz = Plots.heatmap(R_array,z_array,F_Rz',xlabel="R [m]",ylabel="z [m]",title="f(R,z) [m^-3]", xlims=extrema(R_array), ylims=extrema(z_array), aspect_ratio=:equal, dpi=600)

        png(myplt_Ep,filepath_output*"_Ep_plot")
        png(myplt_Rz,filepath_output*"_Rz_plot")
    end
end

# Save results
myfile = jldopen(filepath_output*".jld2",true,true,false,IOStream)
write(myfile,"F_EpRz",F_EpRz)
write(myfile,"E_array",E_array)
write(myfile,"E_array_units","keV")
write(myfile,"p_array",p_array)
write(myfile,"p_array_units","-")
write(myfile,"R_array",R_array)
write(myfile,"R_array_units","m")
write(myfile,"z_array",z_array)
write(myfile,"z_array_units","m")
write(myfile,"constant_Rz",constant_Rz)
write(myfile,"distribution_type","$(distribution_type)")
write(myfile,"ntot",tot_N_FI)
if distribution_type==:gaussian
    write(myfile,"floor_level",floor_level)
    write(myfile,"X_peak",[peak_E, peak_p, peak_R, peak_z])
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
println("~~~createCustomFIDistr.jl completed successfully!~~~")