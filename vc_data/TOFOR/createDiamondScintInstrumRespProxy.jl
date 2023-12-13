# This script will create a proxy for the KM14 and KM15 instrumental response functions
# and save it in the OWCF/vc_data/TOFOR/ folder
folderpath_OWCF = "/home/henrikj/Codes/OWCF/"
cd(folderpath_OWCF)
using Pkg
Pkg.activate(".")
#################################################################### Load packages
using Plots
using JLD2
using Interpolations
#################################################################### Create transfer matrix
En_keV = collect(range(10000.0,stop=19000.0,length=200))
Edep_array = En_keV .- 5702.0
FWHM = 100 # keV
sigma = FWHM / (2 * sqrt(2*log(2)))

response_matrix = zeros(length(En_keV),length(Edep_array))
for i=1:length(En_keV)
    println("i: $(i)")
    x0 = En_keV[i] - 5702.0
    response_matrix[i,:] = (1/(sigma*sqrt(2*pi))) .* exp.((-1/(2*sigma^2)) .* (Edep_array .-x0).^2)
end
################################################################################ Update the response matrix proxy to be more realistic (avoid ultra-low values)
extrema(response_matrix)
input = En_keV
output = Edep_array
Plots.heatmap(input,output,response_matrix)
minimum_nonzero_value = minimum(filter(x-> x>0.0,response_matrix))
response_matrix_floord = replace(x-> (x < 0.01*maximum(response_matrix)) ? 0.0 : x, response_matrix)
minimum_nonzero_value_floord = minimum(filter(x-> x>0.0,response_matrix_floord))
extrema(response_matrix_floord)
Plots.heatmap(input,output,response_matrix_floord)
################################################################################ Update the response matrix proxy to be more realistic (avoid ultra-low values)
myfile = jldopen(folderpath_OWCF*"/vc_data/TOFOR/diamond_scintillator_response_matrix_proxy.jld2",true,true,false)
write(myfile,"input",input)
write(myfile,"output",output)
write(myfile,"response_matrix",response_matrix_floord)
close(myfile)