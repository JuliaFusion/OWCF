########################################################### start_solveInverseProblem_template.jl ###########################################################
# This file contains all the inputs that the script solveInverseProblem.jl needs to solve an inverse problem. That is, a problem on the form s = Wf where we 
# seek the fast-ion distribution f (Vector), given the weight functions W (Matrix) and measurements s (Vector). This file also executes the script 
# solveInverseProblem.jl after the inputs are defined. The measurement data needs to be provided together with uncertainties/errors as explained below. 
# The measurements are considered immutable. Therefore, if the weight functions are not computed for the same measurement bin centers as that of the 
# measurements s, the weight functions will be interpolated onto the measurement bin centers grid. The fast-ion distribution will be reconstructed on an 
# equidistant grid specified via the min_array, max_array and n_array input variables (please see below). The computed reconstruction(s) will be returned 
# together with useful quantities in an .jld2 file (or .hdf5, if requested)(please see below).
#
# The computed reconstruction(s) will be saved as an (N+R)-dimensional array. N is the number of independent variables upon which 
# the fast-ion distribution depend, i.e. f(x1, x2, ..., xN). R is the number of hyperparameters needed to realize the specified regularization scheme 
# (please see the 'regularization' input parameter below). The exact shape of the (N+R)-dimensional array is determined via the 'n_array' and 
# 'nr_array' input variables (please see below). 
#
### The inputs are as follows:
#
# folderpath_OWCF - The path to where the OWCF folder is saved on your computer. End with '/' (or '\' if Windows) - String
# 
# constraints - A Vector of Symbols, to specify what types of constraints that should be put on the fast-ion distribution
#               when solving the inverse problem. Currently, the following option(s) is available:
#                   - :NONNEG - If included, a non-negativity constraint is enforced on the fast-ion distribution at all 
#                               points of the grid (see the 'min_array', 'max_array' and 'n_array' input variables below) 
#                               on which the fast-ion distribution will be reconstructed. Included by default.
#               Include the above option(s) in a Vector, or leave the Vector empty - Vector{Symbol}
# coordinate_system - If an element in 'scriptSources_W' (see below) is specified as "calc2DWeights", then as of the current 
#                     OWCF version, the user needs to manually specify if the (E,p) or the (vpara,vperp) weight functions 
#                     should be used (if both are saved in the same calc2DWeights.jl output file). This is done by 
#                     specifying the 'coordinate_system' input variable to "(E,p)" or "(vpara,vperp)" - String
# excluded_measurement_intervals - A Vector of Vectors of 2-tuples, i.e. [[(X1,X2),(X3,X4),...],[(Y1,Y2),(Y3,Y4),...],...].
#                                Each element in the outermost Vector represents a set of excluded measurement intervals for the 
#                                diagnostic with the corresponding element in the 'filepaths_W' and 'filepaths_S' 
#                                input variables. Each tuple in one of the inner Vectors represents the lower (e.g. X1) and upper (X2) 
#                                bounds of a measurement interval to be excluded. As can be seen in the example above, several 
#                                measurement intervals can be excluded for the same diagnostic. If no measurement intervals are 
#                                to be excluded for a particular diagnostic, simply let the inner Vector be empty for that diagnostic.
#                                If no measurement intervals are to be excluded for any diagnostic, please specify the 
#                                'excluded_measurement_intervals' input variable as a Vector of empty Vectors. The number of (empty) 
#                                Vectors in the outermost Vector should be equal to the length of the 'filepaths_W' and 'filepaths_S'
#                                input variables - Vector{Vector{Tuple{Real,Real}}}
# excluded_measurement_units - A Vector of Strings. The units of the elements of the 'excluded_measurement_intervals' input variable.
#                              Please see OWCF/misc/convert_units.jl for lists of all accepted units of measurement in the OWCF.
#                              For each element in 'excluded_measurement_units', the tuples of the Vector with the corresponding 
#                              index in the 'excluded_measurement_intervals' input variable are assumed to have the specified units. 
#                              That is, the units of the tuples in excluded_measurement_intervals[i] are assumed to be 
#                              excluded_measurement_units[i] - Vector{String}
# exclude_zero_measurements - If true, measurements equal to 0 will be excluded from the reconstruction - Bool
# filepath_start - The path to the start file. This is set automatically, and saved in the calc2DWeights.jl output file. LEAVE UNCHANGED. - String
# filepaths_S - A Vector of Strings with filepath(s) to the signal file(s) to be used as measurements when solving the 
#               inverse problem. Has to be of equal length to the 'filepaths_W' input variable (see below). The files can be in either 
#               .jld2 or .hdf5 (.h5) format. The file(s) needs to have the following  keys (if not computed with an OWCF script):
#                   S - The measurements in the format of a Vector. The length must be equal to nEd of the weight matrix of the 
#                       corresponding file in filepaths_W. That is, e.g. length(filepaths_S[3])==size(filepaths_W[3],1) (pseudo-code) 
#                       must be true. - Vector{Float64}
#                   S_units - The units of the elements of 'S'. Please see OWCF/misc/convert_units.jl for lists of all accepted units 
#                            of measurement in the OWCF. - String
#                   err - The uncertainty of the measurements in the 'S' key. Needs to be the same length as 'S'. That is, 
#                       length(S)==length(err) must be true. - Vector{Float64}
#                   err_units - The units of the elements of 'err'. Please see OWCF/misc/convert_units.jl for lists of all accepted units 
#                              of measurement in the OWCF. - String
#                   Ed_array - The centers of the measurement bins for the measurements in 'S'. Needs to be the same length as 'S'. That is,
#                           length(S)==length(Ed_array) must be true. - Vector{Float64}
#                   Ed_array_units - The units of the elements of 'Ed_array'. Please see OWCF/misc/convert_units.jl for lists of all accepted
#                                   units of measurement in the OWCF. - String
# filepaths_W - A Vector of Strings with filepath(s) to the weight function (matrix) file(s) for solving the inverse problem. The 
#               files can be in either .jld2 or .hdf5 (.h5) format. The file(s) needs to have the following keys (if not 
#               computed with an OWCF script):
#                  W - The weight matrix in inflated format, i.e. nEd x D1 x D2 x ... x DN where nEd is the number 
#                      of measurement bins of Ed_array (see below). DN is the number of grid points in the N:th 
#                      dimension (physical quantity)(e.g. energy, pitch, pm, Rm etc) of W - Array{Float64}
#                  Ed_array - The centers of the measurement bins corresponding to the first dimension (shape) of the weight 
#                             matrix W (see above). - Vector{Float64}
#                  Ed_array_units - The units of the elements of Ed_array (see above). Please see OWCF/misc/convert_units.jl for lists of all 
#                                  accepted units of measurement in the OWCF. This is required information to be able to treat W data correctly 
#                                  - String
#                  D1_array - The grid points of the second dimension (shape) of W - Vector{Float64}
#                  D1_array_units - The units of the elements of D1_array. Please see OWCF/misc/convert_units.jl for lists of all accepted units 
#                                  of measurement in the OWCF. - String
#                  D2_array - The grid points of the third dimension (shape) of W - Vector{Float64}
#                  D2_array_units - The units of the elements of D2_array. Please see OWCF/misc/convert_units.jl for lists of all accepted units 
#                                  of measurement in the OWCF. - String
#                  ...           - .... - ...
#                  DN_array - The grid points of the N+1:th dimension (shape) of W - Vector{Float64}
#                  DN_array_units - The units of the elements of DN_array. Please see -||- - String
#               Please note that the N<=6, because the full fast-ion distribution is six-dimensional (three position-space dimensions, three 
#               velocity-space dimensions). The time dimension is currently not supported.
# folderpath_out - If specified, the output of the solveInverseProblem.jl script will be saved in this folder. - String
# gif_solutions - If set to true, a .gif file of the solutions to the inverse problems will be saved. As well as the 
#                  L-curve and S vs WF plots - Boolean
# min_array - A Vector of minimum values for the grid on which the fast-ion distribution will be reconstructed. min_array[1] must correspond
#             to the same dimension (physical quantity) as the D1_array (see 'filepaths_W' above). min_array[2] must correspond to the same 
#             dimension (physical quantity) as the D2_array. And so on, min_array[N] must correspond to the same dimension (physical quantity)
#             as the DN_array. If left unspecified, the grid of the first weight matrix in filepaths_W will be used - Vector{Float64}
# min_array_units - A Vector of Strings to specify the units of the corresponding elements in the 'min_array' input variable - Vector{String}
# max_array - A Vector of maximum values for the grid on which the fast-ion distribution will be reconstructed. max_array[1] must correspond
#             to the same dimension (physical quantity) as the D1_array (see 'filepaths_W' above). max_array[2] must correspond to the same 
#             dimension (physical quantity) as the D2_array. And so on, max_array[N] must correspond to the same dimension (physical quantity)
#             as the DN_array. If left unspecified, the grid of the first weight matrix in filepaths_W will be used - Vector{Float64}
# max_array_units - A Vector of Strings to specify the units of the corresponding elements in the 'max_array' input variable - Vector{String}
# n_array - A Vector of integer values to specify the number of grid points in each dimension (physical quantity) for the grid on which the 
#           fast-ion distribution will be reconstructed. n_array[1] must correspond to the same dimension (physical quantity) as the D1_array 
#           (see 'filepaths_W' above). n_array[2] must correspond to the same dimension (physical quantity) as the D2_array. And so on, n_array[N] 
#           must correspond to the same dimension (physical quantity) as the DN_array. If left unspecified, the grid of the first weight matrix in 
#           filepaths_W will be used - Vector{Float64}
# nr_array - A Vector of integer values to specify the number of grid points for the regularization strength (λ) for the regularization types included in the 
#           'regularization' input Vector (please see below). length(nr_array)==length(regularization) must hold. Each Integer in nr_array is assumed to specify
#           the number of λ grid points for the corresponding element in the 'regularization' input Vector. The lower and upper bounds of each λ are set 
#           automatically. Some regularization types do not need a λ (e.g. :COLLISIONS). The corresponding element in nr_array will then be treated as a 
#           dummy value. - Vector{Int64}
# noise_floor_factor - If the error (err) of the diagnostic measurements (S) has smaller values than noise_floor_factor*maximum(S), all such 
#                      values will be lifted up (changed) to noise_floor_factor*maximum(S). - Float64
# plot_solutions - If set to true, plots of the solutions to the inverse problems will be saved in .png format. As well as the 
#                  L-curve and S vs WF plots - Boolean
# regularization - A Vector of Symbols, specifying what types of regularization that should be used in the reconstruction. Currently, the 
#                  following options are available:
#                      - :ZEROTIKHONOV - If included, 0th order Tikhonov regularization will be used.
#                      - :FIRSTTIKHONOV - If included, 1st order Tikhonov regularization will be used.
#                      - :COLLISIONS - If included, the fast-ion distribution will be assumed to be spanned by a set of slowing-down functions,
#                                      incorporating collision-physics into the reconstruction. Currently, this is only available for velocity-space (2D)
#                                      reconstructions in (E,p) and (vpara,vperp).
#                      - :ICRF - If included, electromagnetic wave heating in the ion cyclotron resonance range of frequencies (ICRF) will be included 
#                                as prior information when solving the inverse problem. Currently, this is only available for velocity-space (2D) reconstructions 
#                                in (E,p) and (vpara,vperp), orbit-space reconstructions in (3D) (E,mu,Pphi), (E,Lambda,Pphi_n) and (3+1 D) (E,mu,Pphi,sigma), (E,Lambda,Pphi_n,sigma).
#                  One or several of the above options should be specified as elements in a Vector. Otherwise, empty Vector - Vector{Symbol}
# regularization_equil_filepath - If :COLLISIONS and/or :ICRF is included in the 'regularization' Vector, the filepath to a magnetic equilibrium file 
#                                 (either in .eqdsk file format or an OWCF/extra/createCustomMagneticEquilibrium.jl output file) must be specified. For :COLLISIONS, this is to 
#                                 compute necessary quantities for the creation of the slowing-down basis functions, e.g. the Spitzer slowing-down time.
#                                 FOR :ICRF, this is to compute valid orbits, if an orbit-space (3D) reconstruction is to be performed - String
# regularization_thermal_ion_temp - If :COLLISIONS is included in the 'regularization' Vector, data on the thermal ion temperatures must be specified. 
#                                   This should be in the form of a Vector of equal length to the regularization_thermal_ion_species input variable (see below).
#                                   For the elements of the Vector, the following options are currently available:
#                                       - A Float or an Int. The element corresponds to the temperature of the corresponding thermal ion 
#                                         species element in the 'regularization_thermal_ion_species' input variable (see below). The temperatures must 
#                                         be provided in keV.
#                                       - A String. The element should be a filepath to an output file of the helper/createCustomThermalProfiles.jl 
#                                         script. The element is assumed to be the temperature profile of the corresponding thermal ion species element 
#                                         in the 'regularization_thermal_ion_species' input variable (see below).
#                                       - A String. The filepath to a TRANSP output file (in .cdf file format). The thermal ion temperature profile will 
#                                         be loaded from this file. The ion species in the 'regularization_thermal_ion_species' input variable will be used 
#                                         to load the correct data from the TRANSP output file (the file is named e.g. 95679V28.cdf).
# regularization_thermal_ion_dens - If :COLLISIONS is included in the 'regularization' Vector, data on the thermal ion densities must be specified. 
#                                   This should be in the form of a Vector of equal length to the regularization_thermal_ion_species input variable (see below).
#                                   For the elements of the Vector, the following options are currently available:
#                                       - A Float or an Int. The element corresponds to the density of the corresponding thermal ion 
#                                         species element in the 'regularization_thermal_ion_species' input variable (see below). The densities must 
#                                         be provided in m^-3.
#                                       - A String. The element should be a filepath to an output file of the helper/createCustomThermalProfiles.jl 
#                                         script. The element is assumed to be the density profile of the corresponding thermal ion species element 
#                                         in the 'regularization_thermal_ion_species' input variable (see below).
#                                       - A String. The filepath to a TRANSP output file (in .cdf file format). The thermal ion density profile will 
#                                         be loaded from this file. The ion species in the 'regularization_thermal_ion_species' input variable will be used 
#                                         to load the correct data from the TRANSP output file (the file is named e.g. 95679V28.cdf).
# regularization_thermal_ion_species - If :COLLISIONS is included in the 'regularization' Vector, data on the thermal ion species must be specified. This 
#                                      is to be able to compute e.g. the Debye length needed for computing the Spitzer slowing-down time. The thermal ion 
#                                      species present in the plasma must be specified as a Vector of Strings. Please see the OWCF/misc/species_func.jl 
#                                      script for info on how the particle species should be specified as Strings.
# regularization_thermal_electron_temp - If :COLLISIONS is included in the 'regularization' Vector, data on the thermal electron temperature must be specified.
#                                        This can be done in several ways. Currently, the following options are available:
#                                           - A Float. The temperature must be provided in keV
#                                           - A String with the filepath to an output file of the helper/createCustomThermalProfiles.jl script
#                                           - A String with the filepath to a TRANSP output file (in .cdf file format, e.g. 99971L36.cdf).
# regularization_thermal_electron_dens - If :COLLISIONS is included in the 'regularization' Vector, data on the thermal electron density must be specified.
#                                        This can be done in several ways. Currently, the following options are available:
#                                           - A Float. The density must be provided in m^-3
#                                           - A String with the filepath to an output file of the helper/createCustomThermalProfiles.jl script
#                                           - A String with the filepath to a TRANSP output file (in .cdf file format, e.g. 99971L36.cdf).
# regularization_timepoint - If :COLLISIONS is included in the 'regularization' Vector, and a TRANSP output file has been specified for either the ion or 
#                            electron temperatures or densities, a timepoint of interest must be specified. Currently, the following options are available:
#                               - A Float or Int. In seconds
#                               - A String in format "XX.XXXX" or "XX,XXXX". In seconds
#                               - A String with the filepath to a TRANSP-NUBEAM output file (in .cdf file format. The file is named e.g. 95679V28_fi_1.cdf).
# regularization_wave_frequency - If :ICRF is included in the 'regularization' Vector, a wave frequency must be specified for the ICRF wave. As a Float64, in Hertz.
# regularization_cyclotron_harmonic - If :ICRF is included in the 'regularization' Vector, the cyclotron harmonic of the ICRF wave must be specified. As an Int64.
# regularization_toroidal_mode_number - If :ICRF is included in the 'regularization' Vector, and the dimensionality of the coordinate system in which the fast-ion 
#                                       distribution is to be reconstructed is greater than 2, the toroidal mode number of the ICRF wave must be specified. Int64.
# rescale_W - If set to true, the algorithm will rescale the weight functions so that the W*F_ref matrix-Vector product matches the 
#             experimental data for all diagnostics. The F_ref variable is a reference fast-ion distribution (e.g. Maxwellian distribution 
#             or a TRANSP distribution). - Bool
# rescale_W_type - If rescale_W is set to true, the type of rescaling must be specified. Currently, the available options are:
#                      - :MAXIMUM. If used, the weight functions will be rescaled so that the maximum of the WF_ref matrix-Vector product matches  
#                        the maximum of the experimental data for all diagnostics.
#                      - :MEAN. If used, the weight functions will be rescaled so that the mean of the WF_ref matrix-Vector product matches the mean 
#                        of the experimental data for all diagnostics.
# rescale_W_F_ref_source - If rescale_W is set to true, this variable determines the source of F_ref that will be used to rescale the weight functions. 
#                        The currently available options include :FILE and :GAUSSIAN - Symbol
# rescale_W_F_file_path - If rescale_W is set to true, and rescale_W_F_ref_source is set to :FILE, then rescale_W_F_file_path needs to 
#                            be specified as the filepath to a .jld2 file or .hdf5 file loadable by the JLD2to4D() or h5to4D() functions in 
#                            the OWCF/extra/dependencies.jl function library, respectively. Please check their function documentation for 
#                            file key requirements. The rescale_W_F_file_path variable can also be specified as the filepath to a TRANSP-NUBEAM 
#                            output file in .cdf format (e.g. 95679V28_fi_1.cdf), loadable by the CDFto4D() function in the 
#                            OWCF/extra/dependencies.jl function library - String
# rescale_W_F_gaus_N_FI_tot - If rescale_W is set to true and rescale_W_F_ref_source is set to :GAUSSIAN, then rescale_W_F_gaus_N_FI_tot needs to be 
#                            specified as an estimate of the total number of fast ions in the plasma - Float64
# R_of_interest - If the reconstruction is in velocity-space (2D), an (R,z) point of interest might need to be specified, to 
#                 e.g. include collision physics or ICRF as prior information, or to rescale the weight functions using a reference fast-ion distribution.
#                 Otherwise, please leave unspecified - Float64
# R_of_interest_units - The units of R_of_interest. Please see OWCF/misc/convert_units.jl for lists of all accepted units of measurement in the OWCF. - String 
# saveInHDF5format - If set to true, the output of the solveInverseProblem.jl script will be saved in .hdf5 file format, instead of .jld2 file format. - Bool 
# scriptSources_W - A Vector of Strings with the same length as the 'filepaths_W' input variable. Each element of the 
#                   'scriptSources_W' input variable informs the algorithm with which OWCF script the corresponding element 
#                   in the 'filepaths_W' input variable has been computed. For example, if scriptSources_W[2]="calcOrbWeights",
#                   then filepath_W[2] corresponds to a String with the filepath to an output file of the OWCF script 
#                   calcOrbWeights.jl. This is to enable easy I/O between OWCF scripts. If a file in filepaths_W was NOT 
#                   computed with the OWCF, just put 'none' for the corresponding element in scriptSources_W.
# tokamak - The tokamak for which the inverse problem is solved. E.g. JET, ITER etc. Will be used purelely for esthetic purposes. - String
# verbose - If set to true, the script will talk a lot! - Bool
# z_of_interest - If the reconstruction is in velocity-space (2D), an (R,z) point of interest might need to be specified, to 
#                 e.g. include collision physics or ICRF as prior information, or to rescale the weight functions using a reference fast-ion distribution.
#                 Otherwise, please leave unspecified - Float64
# z_of_interest_units - The units of z_of_interest. Please see OWCF/misc/convert_units.jl for lists of all accepted units of measurement in the OWCF. - String  
#
# btipsign - The sign of the dot product between the magnetic field and the plasma current. Can be 1 or -1. - Int64
# FI_species - If the reconstruction space is (vpara,vperp) or :COLLISIONS/:ICRF is included in regularization, the fast-ion species need to be specified. "D", "T", "alpha" etc - String
# h5file_of_nonJulia_origin - If .h5 or .hdf5 files are specified as input, and they were not created with the Julia programming language, the 
#                             h5file_of_nonJulia_origin input variable might need to be set to true. This is because some languages reverse the array 
#                             dimension order when saved as a .h5 or .hdf5 file. E.g. (E,p,R,z) sometimes become (z,R,p,E) - Bool

### Other
# Some lines below are commented out. This is because the inverse problem solving algorithm is currently single-CPU.
# In future version, multi-core processing might be supported, and those lines will then be uncommented.

# Script written by Henrik Järleblad. Last maintained 2025-04-11.
#############################################################################################################################################################

## First you have to set the system specifications
using Distributed # Needed, even though distributed might be set to false. This is to export all inputs to all workers right away, if needed.
#batch_job = false #CURRENTLY NOT AVAILABLE FOR solveInverseProblem.jl
#distributed = true #CURRENTLY NOT AVAILABLE FOR solveInverseProblem.jl
folderpath_OWCF = "" # OWCF folder path. Finish with '/'
#numOcores = 4 #CURRENTLY NOT AVAILABLE FOR solveInverseProblem.jl # When executing script via HPC cluster job, make sure you know how many cores you have requested for your batch job

## Navigate to the OWCF folder and activate the OWCF environment
cd(folderpath_OWCF)
using Pkg
Pkg.activate(".")

#CURRENTLY NOT AVAILABLE FOR solveInverseProblem.jl
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

#CURRENTLY NOT AVAILABLE FOR solveInverseProblem.jl
## If running locally and multi-threaded
#if !batch_job && distributed # Assume you are executing the script on a local laptop (/computer)
#    println("Adding processes... ")
#    addprocs(numOcores-(nprocs()-1)) # If you didn't execute this script as an HPC cluster job, then you need to add processors like this. Add all remaining available cores.
#    # The '-(nprocs()-1)' part is simply to ensure to extra processes are added, in case script needs to be restarted on a local computer
#end

## -----------------------------------------------------------------------------
@everywhere begin
    constraints = [:NONNEG]
    coordinate_system = "" # Specify as e.g. "(E,p)", "(vpara,vperp)", "(E,pm,Rm)" or "(E,Lambda,Pphi_n)"
    excluded_measurement_intervals = [[(,),...],[(,),...]]
    excluded_measurement_units = ["",""]
    exclude_zero_measurements = true
    filepath_start = Base.source_path() # For saving in output file, for posterity and reproducibility
    filepaths_S = ["",""]
    filepaths_W = ["",""]
    folderpath_out = ""
    gif_solutions = false
    min_array = []
    min_array_units = []
    max_array = []
    max_array_units = []
    n_array = []
    noise_floor_factor = 0.0
    plot_solutions = false
    nr_array = [20] # Should be equal in length to 'regularization'
    regularization = [:ZEROTIKHONOV] # Should be equal in length to 'nr_array'
    if (:COLLISIONS in regularization) || (:ICRF in regularization) # :COLLISIONS and/or :ICRF included in regularization...
        regularization_equil_filepath = "" # Please specify the file path to a magnetic equilibrium file (see template description above)
    end
    if :COLLISIONS in regularization
        regularization_thermal_ion_temp = [] # (please read start file template description above)
        regularization_thermal_ion_dens = [] # (please read start file template description above)
        regularization_thermal_ion_species = ["",""] # Or [""], ["","",""] etc 
        regularization_thermal_electron_temp = 0.0 # or "" (please read start file template description above)
        regularization_thermal_electron_dens = 0.0 # or "" (please read start file template description above)
        regularization_timepoint = 00.0000 # or "" (please read start file template description above)
    end
    if :ICRF in regularization
        regularization_wave_frequency = 0.0 # The ICRF wave frequency (Hz)
        regularization_cyclotron_harmonic = 0 # The ICRF wave cyclotron harmonic (preferably an integer)
        regularization_toroidal_mode_number = 0 # The ICRF wave toroidal mode number
    end
    rescale_W = true
    if rescale_W
        rescale_W_type = :MAXIMUM # Currently supported, :MAXIMUM and :MEAN
        rescale_W_F_ref_source = :FILE # Currently supported, :FILE and :GAUSSIAN
        if rescale_W_F_ref_source == :FILE
            rescale_W_F_file_path = ""
        end
        if rescale_W_F_ref_source == :GAUSSIAN
            rescale_W_F_gaus_N_FI_tot = 0.0e00
        end
    end
    R_of_interest = nothing
    R_of_interest_units = "m"
    saveInHDF5format = false
    scriptSources_W = ["",""]
    tokamak = ""
    verbose = false
    z_of_interest = nothing
    z_of_interest_units = "m"

    ### Extra input arguments
    btipsign = 1 # The sign of the dot product between the magnetic field and the plasma current. Most likely not needed. But please specify, if known.
    FI_species = "" # If reconstruction is in (vpara,vperp) or the 'regularization' input variable contains :COLLISIONS or :ICRF, the fast-ion species might need to be specified. Please see OWCF/misc/species_func.jl for more info on species options.
    h5file_of_nonJulia_origin = false # If rescale_W_F_file_path is an .h5 or .hdf5 file, and it was not created with the Julia programming language, this input variable should probably be set to true
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
include("solveInverseProblem.jl")

## -----------------------------------------------------------------------------
# Then you clean up after yourself
if batch_job && distributed
    for i in workers()
        rmprocs(i)
    end
end