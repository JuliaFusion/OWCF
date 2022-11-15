############################# howToLoadFromTRANSP_FI_CDF #############################
# This short file explains, and provides examples of, how to load data from a TRANSP
# .cdf file. Both a shot file (e.g. 99965K71.cdf) and a NUBEAM file (e.g. 99965K71_fi_1.cdf).
#
# In short, to load the data from a TRANSP shot file, do the following:
# 1. Use ctrl+F to locate the identifier you seek in the OWCF/misc/transp_outputs.txt file.
#    (e.g. the total beam-target neutron data has the 'BTNTX' identifier)
#    (e.g. use ctrl+F and search for 'neutron')
# 2. Use ncread(filepath_CDF,[IDENTIFIER]) or similar to load the data you need.
#    (e.g. ncread(filepath_CDF,"BTNTX") will return the total beam-target neutron rate as a function 
#     of normalized minor plasma radius (r/a) and time t)
# 3. Use the data (e.g. read it and plot it as the examples below show)
#
# You can find more traces of information on the TRANSP identifiers in the read_nubeam() 
# function of the OWCF/extra/getEpRzFIdistrFromTRANSP.jl script.
#
# Structured by Henrik Järleblad. Last modified 2022-11-15
###################################################################################### 

## ------
folderpath_OWCF = ""
println("Loading packages... ")
cd(folderpath_OWCF)
using Pkg
Pkg.activate(".")
using HDF5
using NetCDF

## ------
filepath_cdf = ""
ncinfo(filepath_cdf)
## ------

## ------
TIME = round((ncread(filepath_cdf,"TIME"))[1],digits=4)
println("Time: $(TIME)")
## ------
filepath_FI_cdf = ""
ncinfo(filepath_FI_cdf)
## ------
ND = ncread(filepath_FI_cdf,"ND")
ND = ncread(filepath_FI_cdf,"NT")
## ------
NTOT_T_NBI = ncread(filepath_FI_cdf,"NTOT_T_NBI")
println("NTOT: $(NTOT_T_NBI[1])")
## ------
#TRANSP_RUNID = ncread(filepath_FI_cdf,"TRANSP_RUNID")
#TRANSP_RUNID = (split((NetCDF.nc_char2string(ncread(filepath_FI_cdf,"TRANSP_RUNID")))[1]))[1]
#println("TRANSP RUN-ID: "*TRANSP_RUNID)
## ------
#TOKAMAK = (split((NetCDF.nc_char2string(ncread(filepath_FI_cdf,"TOKAMAK")))[1]))[1] # This is just annoying...
#println("TOKAMAK: "*TOKAMAK)
## ------

# Investigate the ratio of beam-target neutrons to thermonuclear neutrons in a TRANSP shot
folderpath_OWCF = ""
cd(folderpath_OWCF)
using Pkg
Pkg.activate(".")
###########################################################################
using NetCDF
using Plots
###########################################################################
filepath_TRANSP = ""
###########################################################################
THERMONUCLEAR_NEUTRONS = ncread(filepath_TRANSP,"THNTX")
###########################################################################
BEAMTARGET_NEUTRONS = ncread(filepath_TRANSP,"BTNTX")
###########################################################################
TIME = ncread(filepath_TRANSP,"TIME")
println("Time: $(TIME)")
X = ncread(filepath_TRANSP,"X")
###########################################################################

for it=1:length(TIME)
    myplt = Plots.plot(X[:,it],THERMONUCLEAR_NEUTRONS[:,it],xlabel="r/a [-]",title="t= $(round(TIME[it],digits=3)) s", ylabel="Neutron yield [(cm3*s)^-1]", label="Thermonuclear")
    myplt = Plots.plot!(myplt, X[:,it],BEAMTARGET_NEUTRONS[:,it],label="BEAMTARGET")
    display(myplt)
end