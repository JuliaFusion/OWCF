## ------
folderpath_OWCF = "G:/My Drive/DTU/codes/OWCF/"
println("Loading packages... ")
cd(folderpath_OWCF)
using Pkg
Pkg.activate(".")
using HDF5
using NetCDF

## ------
#filepath_cdf = "E:/DTU/codes/OWCF/TRANSP/JET/99965/K71/99965K71_fi_2.cdf"
filepath_cdf = "G:/My Drive/DTU/codes/OWCF/TRANSP/JET/94701/V01/94701V01_fi_2.cdf"
ncinfo(filepath_cdf)
## ------
TIME = round((ncread(filepath_cdf,"TIME"))[1],digits=4)
println("Time: $(TIME)")
## ------
ND = ncread(filepath_cdf,"ND")
ND = ncread(filepath_cdf,"NT")
## ------
NTOT_T_NBI = ncread(filepath_cdf,"NTOT_T_NBI")
println("NTOT: $(NTOT_T_NBI[1])")
## ------
#TRANSP_RUNID = ncread(filepath_cdf,"TRANSP_RUNID")
#TRANSP_RUNID = (split((NetCDF.nc_char2string(ncread(filepath_cdf,"TRANSP_RUNID")))[1]))[1]
#println("TRANSP RUN-ID: "*TRANSP_RUNID)
## ------
#TOKAMAK = (split((NetCDF.nc_char2string(ncread(filepath_cdf,"TOKAMAK")))[1]))[1] # This is just annoying...
#println("TOKAMAK: "*TOKAMAK)
## ------
