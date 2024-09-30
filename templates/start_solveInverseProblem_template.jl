#################################### start_solveInverseProblem_template.jl #####################################
# THIS SCRIPT IS CURRENTLY UNDER CONSTRUCTION!!!
# THIS SCRIPT IS CURRENTLY UNDER CONSTRUCTION!!!
# THIS SCRIPT IS CURRENTLY UNDER CONSTRUCTION!!!
# This file contains all the inputs that the script solveInverseProblem.jl needs to solve an inverse problem.
# That is, a problem on the form s = Wf where we seek f, given W and s. This file also executes the script 
# solveInverseProblem.jl after the inputs are defined. The measurement data needs to be provided together with 
# uncertainties/errors as explained below. The computed reconstruction(s) will be returned together with useful 
# quantities in an .jld2 file (or .hdf5, if requested).
#
# The inputs are as follows:
# 
# filepaths_W - An array (of strings) with filepath(s) to the weight function file(s) for solving the inverse problem. The 
#               files can be in either .jld2 or .hdf5 (.h5) format. The file(s) needs to have the following keys (if not 
#               computed with an OWCF script):
#               W - The weight matrix in inflated format, i.e. nEd x D1 x D2 x ... x Dn where nEd is the number 
#                   of measurement bins of Ed_array (see below). Dn is the number of grid points in the n:th 
#                   dimension (physical quantity)(e.g. energy, pitch, pm, Rm etc) of W - Array{Float64}
#               Ed_array - The centers of the measurement bins corresponding to the first dimension (shape) of the weight 
#                          matrix W (see above). - Vector{Float64}
#               Ed_array_dim - The dimension (physical quantity) of the elements of Ed_array (see above). Can be, 
#                              for example, "keV", "MeV", "m", "cm", "-" (dimensionless), "a.u." (dimensionless/unknown),
#                              "keeVe" etc. Required information to be able to treat W data correctly - String
#               D1_array - The grid points of the second dimension (shape) of W - Vector{Float64}
#               D1_array_dim - The dimension (physical quantity) of the elements of D1_array. Can be, for example, "keV", 
#                              "MeV", "m", "cm" and so on (same as Ed_array_dim) - String
#               D2_array - The grid points of the third dimension (shape) of W - Vector{Float64}
#               D2_array_dim - The dimension (physical quantity) of the elements of D2_array. Can be, for example, "keV",
#                              "MeV", "m", "cm" and so on (same as Ed_array_dim) - String
#               ...          - .... - ...
#               Dn_array - The grid points of the n+1:th dimension (shape) of W - Vector{Float64}
#               Dn_array_dim - The dimension (physical quantity) of the elements of Dn_array. Can be, -||- - String
# scriptSources_W - An array (of strings) with the same length as the 'filepaths_W' input variable. Each element of the 
#                   'scriptSources_W' input variable informs the algorithm with which OWCF script the corresponding element in 
#                   the 'filepaths_W' input variable has been computed. For example, if scriptSources_W[2]="calcOrbWeights",
#                   then filepath_W[2] corresponds to a string with the filepath to an output file of the OWCF script 
#                   calcOrbWeights.jl. This is to enable easy I/O between OWCF scripts. If a file in filepaths_W was NOT 
#                   computed with the OWCF, just put 'none' for the corresponding element in scriptSources_W.
# filepaths_S - An array (of strings) with filepath(s) to the signal file(s) to be used as measurements when solving the 
#               inverse problem. The files can be in either .jld2 or .hdf5 (.h5) format. The file(s) needs to have the following keys
#               (if not computed with an OWCF script):
#               S - The measurements in the format of a vector. The length must be equal to nEd of the weight matrix of the corresponding
#               file in filepaths_W. That is, e.g. length(filepaths_S[3])==size(filepaths_W[3],1) (pseudo-code) must be true. - Vector{Float64}
#               S_dim - The dimension (physical quantity) of the elements of 'S'. Can be, for example, "keV", "MeV", "m", "cm" and so on - String
#               err - The uncertainty of the measurements in the 'S' key. Needs to be of the same length as 'S'. That is, length(S)==length(err)
#               must be true. - Vector{Float64}
#               err_dim - The dimension (physical quantity) of the elements of 'err'. Can be, for example, "keV" etc. - String
# 
#               