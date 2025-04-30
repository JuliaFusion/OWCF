############---------------------------------------------------------------------------------------###
# A script to test a lot of the tools in the OWCF/extra/dependencies.jl script

# PLEASE NOTE! THIS SCRIPT IS UNDER ACTIVE DEVELOPMENT!!!
# PLEASE NOTE! THIS SCRIPT IS UNDER ACTIVE DEVELOPMENT!!!
# PLEASE NOTE! THIS SCRIPT IS UNDER ACTIVE DEVELOPMENT!!!
###------------------------------------------------------------------------------------------------###

include(folderpath_OWCF*"extra/dependencies.jl")
plot_test_results = plot_test_results # SET TO TRUE, VIA THE plot_test_results INPUT VARIABLE IN THE OWCF/tests/run_tests.jl SCRIPT

############---------------------------------------------------------------------------------------###
μ = [100.0, 0.6, 3.3]
σ = [25.0, 0.1, 0.1]
test_gaussian = gaussian(μ, σ; mx=[200.0, 1.0, 3.8], mn=[0.0, -1.0, 3.0], n=[30,61,64], floor_level=0.001)
if plot_test_results
    date_and_time = split("$(Dates.now())","T")[1]*"at"*split("$(Dates.now())","T")[2][1:5]
    abscissa1 = collect(range(0.0,stop=200.0,length=30)); d1 = diff(abscissa1)[1]
    abscissa2 = collect(range(-1.0,stop=1.0,length=61)); d2 = diff(abscissa2)[1]
    abscissa3 = collect(range(3.0,stop=3.8,length=64)); d3 = diff(abscissa3)[1]
    test_gaussian1 = (d1) .*dropdims(sum(test_gaussian,dims=1),dims=1)
    test_gaussian2 = (d2) .*dropdims(sum(test_gaussian,dims=2),dims=2)
    test_gaussian3 = (d3) .*dropdims(sum(test_gaussian,dims=3),dims=3)
    myplt1 = Plots.heatmap(abscissa2,abscissa3,transpose(test_gaussian1),dpi=100,title="gaussian() $(date_and_time)",titlefontsize=12)
    myplt2 = Plots.heatmap(abscissa1,abscissa3,transpose(test_gaussian2),dpi=100)
    myplt3 = Plots.heatmap(abscissa1,abscissa2,transpose(test_gaussian3),dpi=100)
    myplt = Plots.plot(myplt1,myplt2,myplt3,layout=(3,1),size=(400,1200))
    png(myplt,folderpath_OWCF*"tests/outputs/dependencies_test_gaussian")
end
###------------------------------------------------------------------------------------------------###

############---------------------------------------------------------------------------------------###
y = erf.(collect(range(-5.0,stop=5.0,length=100)); resolution=2000, sigma=1/sqrt(3))

y = erfc.(collect(range(-5.0,stop=5.0,length=100)); resolution=2000, sigma=1/sqrt(3))
###------------------------------------------------------------------------------------------------###

############---------------------------------------------------------------------------------------###
L1 = forward_difference_matrix((33,31))
###------------------------------------------------------------------------------------------------###

############---------------------------------------------------------------------------------------###
test_Rz_grid = rz_grid(2.8,4.0,33,-2.0,2.0,34; nphi=720)
###------------------------------------------------------------------------------------------------###

############---------------------------------------------------------------------------------------###
M, wall = read_geqdsk(folderpath_OWCF*"equilibrium/JET/g96100/g96100_0-53.0012.eqdsk",clockwise_phi=false)
test_FI_species = "D"; test_n_e = 1.0e20; test_T_e = 5.0; test_species_th_vec = ["D","T"]
test_n_th_vec = [0.2e20,0.8e20]; test_T_th_vec = [4.0,3.0]

test_gcp = getGCP(test_FI_species; E=100.0, p=0.6, R=3.2, z=0.3)
test_debye_length = debye_length(test_n_e, test_T_e, test_species_th_vec, test_n_th_vec, test_T_th_vec)

test_gyro_radius = gyro_radius(M,test_gcp)
test_spitzer_slowdown_time = spitzer_slowdown_time(test_n_e, test_T_e, test_FI_species, test_species_th_vec, test_n_th_vec, test_T_th_vec; plasma_model = :texas)
test_spitzer_slowdown_time = spitzer_slowdown_time(test_n_e, test_T_e, test_FI_species, test_species_th_vec, test_n_th_vec, test_T_th_vec; plasma_model = :salewski)
###------------------------------------------------------------------------------------------------###

############---------------------------------------------------------------------------------------###
test_R = 3.2
test_z = 0.4
###------------------------------------------------------------------------------------------------###

############---------------------------------------------------------------------------------------###
test_p_array = collect(range(-1.0,stop=1.0,length=10)); test_E_array = collect(range(5.0,stop=100.0,length=20))
test_icrf_freq = 32.6e6 # Hz
test_icrf_harmonic = 1

epsilon = icrf_streamlines(M, [test_E_array, test_p_array], test_FI_species, test_icrf_freq, test_icrf_harmonic; 
                           R_of_interest=test_R, z_of_interest=test_z)

dE = diff(test_E_array)[1]
dp = diff(test_p_array)[1]
E_points = [test_E_array[c[1]] for c in CartesianIndices((length(test_E_array),length(test_p_array)))][:]
p_points = [test_p_array[c[2]] for c in CartesianIndices((length(test_E_array),length(test_p_array)))][:]
dE_points = [epsilon[c][1]*dE for c in CartesianIndices((length(test_E_array),length(test_p_array)))][:]
dp_points = [epsilon[c][2]*dp for c in CartesianIndices((length(test_E_array),length(test_p_array)))][:]

if plot_test_results
    test_icrf_streamlines_plt_1 = Plots.quiver(E_points, p_points, quiver= (dE_points,dp_points),size=(900,900),arrows=Plots.arrow(:closed,:head,1,0.5))
end
###------------------------------------------------------------------------------------------------###

############---------------------------------------------------------------------------------------###
test_vpa_array = collect(range(-5.0e5,stop=5.0e5,length=40))
test_vpe_array = collect(range(0.0,stop=1.0e6,length=20))

epsilon = icrf_streamlines(M, [test_vpa_array, test_vpe_array], test_FI_species, test_icrf_freq, test_icrf_harmonic;
                           R_of_interest=test_R, z_of_interest=test_z, coord_system="(vpara,vperp)")

dvpa = diff(test_vpa_array)[1]
dvpe = diff(test_vpe_array)[2]

vpa_points = [test_vpa_array[c[1]] for c in CartesianIndices((length(test_vpa_array),length(test_vpe_array)))][:]
vpe_points = [test_vpe_array[c[2]] for c in CartesianIndices((length(test_vpa_array),length(test_vpe_array)))][:]
dvpa_points = [epsilon[c][1]*dvpa for c in CartesianIndices((length(test_vpa_array),length(test_vpe_array)))][:]
dvpe_points = [epsilon[c][2]*dvpe for c in CartesianIndices((length(test_vpa_array),length(test_vpe_array)))][:]

if plot_test_results
    test_icrf_streamlines_plt_2 = Plots.quiver(vpa_points, vpe_points, quiver= (dvpa_points,dvpe_points),size=(900,900),arrows=Plots.arrow(:closed,:head,1,0.5))
end
###------------------------------------------------------------------------------------------------###

############---------------------------------------------------------------------------------------###
test_E_array = collect(range(5.0,stop=1000.0,length=10))
test_Lambda_array = collect(range(0.0,stop=2.0,length=11))
test_Pphi_n_array = collect(range(-2.0,stop=2.0,length=12))

epsilon = icrf_streamlines(M, [test_E_array, test_Lambda_array, test_Pphi_n_array, [-1,1]], test_FI_species, test_icrf_freq, test_icrf_harmonic; 
                           coord_system="(E,Lambda,Pphi_n,sigma)")

iPphi = 6
epsilon_iPphi = epsilon[:,:,iPphi,1]

dE = diff(test_E_array)[1]
dLambda = diff(test_Lambda_array)[1]

# 2D quiver plot
E_points = [test_E_array[c[1]] for c in CartesianIndices((length(test_E_array),length(test_Lambda_array)))][:]
Lambda_points = [test_Lambda_array[c[2]] for c in CartesianIndices((length(test_E_array),length(test_Lambda_array)))][:]
dE_points = [epsilon_iPphi[c][1]*dE for c in CartesianIndices((length(test_E_array),length(test_Lambda_array)))][:]
dLambda_points = [epsilon_iPphi[c][2]*dLambda for c in CartesianIndices((length(test_E_array),length(test_Lambda_array)))][:]

if plot_test_results
    test_icrf_streamlines_plt_3 = Plots.quiver(E_points, Lambda_points, quiver = (dE_points,dLambda_points), size=(900,900),arrows=Plots.arrow(:closed,:head,1,0.5))
end
###------------------------------------------------------------------------------------------------###

############---------------------------------------------------------------------------------------###
test_FI_species_q = getSpeciesCharge(test_FI_species)
psi_axis, psi_bdry = psi_limits(M)
if psi_bdry==0
    @warn "The magnetic flux at the last closed flux surface (LCFS) is found to be 0 for the 'M' input to os2COM(). Pϕ_n=Pϕ/(q*|Ψ_w|) where Ψ_w=Ψ(mag. axis) is assumed instead of Ψ_w=Ψ(LCFS)."
    Ψ_w_norm = abs(psi_axis)
else
    Ψ_w_norm = abs(psi_bdry)
end
test_E_array = collect(range(5.0,stop=1000.0,length=10))
test_mu_array = (1000*GuidingCenterOrbits.e0*test_E_array[2]) .*collect(range(0.0,stop=2.0,length=11)) ./norm(Equilibrium.Bfield(M,magnetic_axis(M)...))
test_Pphi_array = (test_FI_species_q*Ψ_w_norm) .*collect(range(-2.0,stop=2.0,length=12))

epsilon = icrf_streamlines(M, [test_E_array, test_mu_array, test_Pphi_array], test_FI_species, test_icrf_freq, test_icrf_harmonic; 
                           coord_system="(E,mu,Pphi)")

iPphi = 6
epsilon_iPphi = epsilon[:,:,iPphi]

# 2D quiver plot
dE = diff(test_E_array)[1]
dmu = diff(test_mu_array)[1]

E_points = [test_E_array[c[1]] for c in CartesianIndices((length(test_E_array),length(test_mu_array)))][:]
mu_points = [test_mu_array[c[2]] for c in CartesianIndices((length(test_E_array),length(test_mu_array)))][:]
dE_points = [epsilon_iPphi[c][1]*dE for c in CartesianIndices((length(test_E_array),length(test_mu_array)))][:]
dmu_points = [epsilon_iPphi[c][2]*dmu for c in CartesianIndices((length(test_E_array),length(test_mu_array)))][:]

if plot_test_results
    test_icrf_streamlines_plt_4 = Plots.quiver(E_points, mu_points, quiver = (dE_points,dmu_points), size=(900,900),arrows=Plots.arrow(:closed,:head,1,0.5))
end
###------------------------------------------------------------------------------------------------###

############---------------------------------------------------------------------------------------###
if plot_test_results
    myplt = Plots.plot(test_icrf_streamlines_plt_1,
                    test_icrf_streamlines_plt_2,
                    test_icrf_streamlines_plt_3,
                    test_icrf_streamlines_plt_4,
                    layout=(2,2),size=(1800,1800))
    display(myplt)
    png(myplt,folderpath_OWCF*"tests/outputs/dependencies_test_icrf_streamlines")
end
###------------------------------------------------------------------------------------------------###

############---------------------------------------------------------------------------------------###
x = collect(range(0.0,stop=2*pi,length=100))
S = sin.(x)
S_noise, err_noise = add_noise(S, 0.05; k=0.05)
err_noise_estimate = estimate_noise(S_noise)
if plot_test_results
    test_noise_plt_1 = Plots.plot(x,S_noise,label="")
    test_noise_plt_2 = Plots.plot(x,err_noise,label="True err level")
    test_noise_plt_2 = Plots.plot!(test_noise_plt_2,x,err_noise_estimate,label="Est. err level")
end
S_noise, err_noise = add_noise(S, 0.05; k=0.15)
err_noise_estimate = estimate_noise(S_noise)
if plot_test_results
    test_noise_plt_3 = Plots.plot(x,S_noise)
    test_noise_plt_4 = Plots.plot(x,err_noise,label="True err level")
    test_noise_plt_4 = Plots.plot!(test_noise_plt_4,x,err_noise_estimate,label="Est. err level")
end
if plot_test_results
    myplt = Plots.plot(test_noise_plt_1, test_noise_plt_3, 
                       test_noise_plt_2, test_noise_plt_4,
                       layout=(2,2),size=(1600,1600))
    display(myplt)
    png(myplt,folderpath_OWCF*"tests/outputs/dependencies_test_S_noise")
end
###------------------------------------------------------------------------------------------------###

############---------------------------------------------------------------------------------------###
myfile = jldopen(folderpath_OWCF*"vc_data/TOFOR/TOFOR_instrumental_response.jld2",false,false,false,IOStream)
input_abscissa = myfile["input"]
output_abscissa = myfile["output"]
response_matrix = myfile["response_matrix"]
close(myfile)
Ed_array = collect(range(1000.0,stop=5000.0,length=100))
S_raw = exp.(((Ed_array .- 2500.0) .^2) ./(-1.0e6))
S = apply_instrumental_response(S_raw, Ed_array, input_abscissa, output_abscissa, response_matrix)

if plot_test_results
    test_instrum_resp_plt_1 = Plots.plot(Ed_array,S_raw,label="",title="Raw signal")
    test_instrum_resp_plt_2 = Plots.plot(output_abscissa,S,label="",title="Signal (w. instrum. resp.)")
    myplt = Plots.plot(test_instrum_resp_plt_1,test_instrum_resp_plt_2,layout=(1,2),size=(1600,800))
    display(myplt)
    png(myplt,folderpath_OWCF*"tests/outputs/dependencies_test_instrum_resp")
end
###------------------------------------------------------------------------------------------------###

############---------------------------------------------------------------------------------------###
if !is_vpara_vperp([test_vpa_array,test_vpe_array], ["m_s^-1","m_s^-1"])
    error("is_vpara_vperp() returned false when it should have returned true.")
end
if is_vpara_vperp([test_E_array,test_vpe_array], ["keV","m_s^-1"])
    error("is_vpara_vperp() returned true when it should have returned false.")
end
###------------------------------------------------------------------------------------------------###

############---------------------------------------------------------------------------------------###
if !is_normalized_COM([test_E_array,test_Lambda_array,test_Pphi_n_array,[-1,1]],["keV","-","-","-"])
    error("is_normalized_COM() returned false when it should have returned true.")
end
if is_normalized_COM([test_vpa_array,test_Lambda_array,test_Pphi_n_array,[-1,1]],["m_s^-1","-","-","-"])
    error("is_normalized_COM() returned true when it should have returned false.")
end
###------------------------------------------------------------------------------------------------###

############---------------------------------------------------------------------------------------###
if !is_EpRz([test_E_array,test_p_array,[2.8,3.0,3.3],[-0.5,0.0,0.5]],["keV","-","m","m"])
    error("is_EpRz() returned false when it should have returned true.")
end
if is_EpRz([test_vpa_array,test_p_array,[2.8,3.0,3.3],[-0.5,0.0,0.5]],["m_s^-1","-","m","m"])
    error("is_EpRz() returned true when it should have returned false.")
end
###------------------------------------------------------------------------------------------------###

############---------------------------------------------------------------------------------------###
L1_ICRF = icrf_regularization_matrix(M, [test_E_array, test_Pphi_n_array, test_Lambda_array, [-1,1]], test_FI_species, test_icrf_freq, test_icrf_harmonic; 
                           coord_system="(E,Pphi_n,Lambda,sigma)")
###------------------------------------------------------------------------------------------------###

############---------------------------------------------------------------------------------------###
filepath_hdf5 = jld2tohdf5(folderpath_OWCF*"tests/outputs/createCustomFIDistr_test1.jld2")
if !(h5to4D(filepath_hdf5)==JLD2to4D(folderpath_OWCF*"tests/outputs/createCustomFIDistr_test1.jld2"))
    error("I/O discrepancy between h5to4D() and JLD2to4D functions.")
end
test_F_EpRz, test_E_array, test_p_array, test_R_array, test_z_array = h5to4D(filepath_hdf5)
test_Eq_array = collect(range(extrema(test_E_array)...,length=11))
test_pq_array = collect(range(extrema(test_p_array)...,length=10))
test_Rq_array = collect(range(extrema(test_R_array)...,length=12))
test_zq_array = collect(range(extrema(test_z_array)...,length=13))
test_F_EpRz_lowres = interpFps(test_F_EpRz, test_E_array, test_p_array, test_R_array, test_z_array, test_Eq_array, test_pq_array, test_Rq_array, test_zq_array)
###------------------------------------------------------------------------------------------------###

############---------------------------------------------------------------------------------------###
dvols = get4DVols(test_Eq_array, test_pq_array, test_Rq_array, test_zq_array)
###------------------------------------------------------------------------------------------------###

############---------------------------------------------------------------------------------------###
test_E_array = collect(range(5.0,stop=1000.0,length=11))
test_pm_array = collect(range(-1.0,stop=1.0,length=19))
test_Rm_array = collect(range(magnetic_axis(M)[1]-0.1,stop=maximum(wall.r),length=20))
dE = diff(test_E_array)[1]
dpm = diff(test_pm_array)[1]
dRm = diff(test_Rm_array)[1]

extra_kw_args=Dict(:toa => true, :limit_phi => true, :max_tries => 0)
orbs, og = OrbitTomography.orbit_grid(M, test_E_array, test_pm_array, test_Rm_array; q=getSpeciesEcu(test_FI_species), amu=getSpeciesAmu(test_FI_species), wall=wall, extra_kw_args...)

test_F_os_3D = gaussian([test_E_array[3],test_pm_array[15],test_Rm_array[15]],[100.0,0.3,0.08]; mx=maximum.([test_E_array,test_pm_array,test_Rm_array]), mn=minimum.([test_E_array,test_pm_array,test_Rm_array]), n=length.([test_E_array,test_pm_array,test_Rm_array]), floor_level=0.001)
test_F_os_3D = OWCF_map_orbits(og,unmap_orbits(og,test_F_os_3D),true) # This ensures only pixels corresponding to valid orbits are kept
###------------------------------------------------------------------------------------------------###

############---------------------------------------------------------------------------------------###
test_F_COM_3D, test_E_array_os2com, test_Lambda_array_os2com, test_Pphi_n_array_os2com = os2COM(M, test_F_os_3D, test_E_array, test_pm_array, test_Rm_array, test_FI_species; needJac=true, wall=wall) # FI distribution
dE_os2com = diff(test_E_array_os2com)[1] # Should be the same as dE above
dLambda = diff(test_Lambda_array_os2com)[1]
dPphi_n = diff(test_Pphi_n_array_os2com)[1]

test_os_Ntot = sum((dE*dpm*dRm) .*test_F_os_3D)
test_COM_Ntot = sum((dE_os2com*dLambda*dPphi_n) .*(test_F_COM_3D[:,:,:,1] + test_F_COM_3D[:,:,:,2]))

Ntot_diff = abs((test_os_Ntot-test_COM_Ntot)/test_os_Ntot)
if Ntot_diff>0.03 # For this test, due to low-resolution grids, the difference between f(E,pm,Rm) and f(E,Λ,Pϕ_n;σ) in terms of total number of fast ions is approx. 2.4% (i.e. a <3% criterion will suffice to ensure correct functionality)
    error("(E,pm,Rm) → (E,Λ,Pϕ_n;σ) transform via OWCF/extra/dependencies.jl/os2COM() did not preserve the total number of fast ions, i.e. ∫f(E,pm,Rm)*dE*dpm*dRm ≠ ∫(f(E,Λ,Pϕ_n;σ=-1)+f(E,Λ,Pϕ_n;σ=+1))*dE*dΛ*dPϕ_n")
end
###------------------------------------------------------------------------------------------------###

############---------------------------------------------------------------------------------------###
if plot_test_results
    myplt1 = Plots.heatmap(test_Rm_array, test_pm_array, dropdims(sum(dE .*test_F_os_3D,dims=1),dims=1), xlabel="Rm [m]", ylabel="pm [-]")
    myplt2 = Plots.heatmap(test_Rm_array, test_E_array, dropdims(sum(dpm .*test_F_os_3D,dims=2),dims=2), xlabel="Rm [m]", ylabel="E [keV]")
    myplt3 = Plots.heatmap(test_E_array, test_pm_array, transpose(dropdims(sum(dRm .*test_F_os_3D,dims=3),dims=3)), xlabel="E [keV]", ylabel="pm [-]")
    myplt4 = Plots.heatmap(test_Pphi_n_array_os2com, test_Lambda_array_os2com, dropdims(sum(dE_os2com .*test_F_COM_3D[:,:,:,1],dims=1),dims=1), title="σ=-1", xlabel="Pϕ_n [-]", ylabel="Λ [-]")
    myplt5 = Plots.heatmap(test_Pphi_n_array_os2com, test_E_array_os2com, dropdims(sum(dLambda .*test_F_COM_3D[:,:,:,1],dims=2),dims=2), title="σ=-1", xlabel="Pϕ_n [-]", ylabel="E [keV]")
    myplt6 = Plots.heatmap(test_E_array_os2com, test_Lambda_array_os2com, transpose(dropdims(sum(dPphi_n .*test_F_COM_3D[:,:,:,1],dims=3),dims=3)), title="σ=-1", xlabel="E [keV]", ylabel="Λ [-]")
    myplt7 = Plots.heatmap(test_Pphi_n_array_os2com, test_Lambda_array_os2com, dropdims(sum(dE_os2com .*test_F_COM_3D[:,:,:,2],dims=1),dims=1), title="σ=+1", xlabel="Pϕ_n [-]", ylabel="Λ [-]")
    myplt8 = Plots.heatmap(test_Pphi_n_array_os2com, test_E_array_os2com, dropdims(sum(dLambda .*test_F_COM_3D[:,:,:,2],dims=2),dims=2), title="σ=+1", xlabel="Pϕ_n [-]", ylabel="E [keV]")
    myplt9 = Plots.heatmap(test_E_array_os2com, test_Lambda_array_os2com, transpose(dropdims(sum(dPphi_n .*test_F_COM_3D[:,:,:,2],dims=3),dims=3)), title="σ=+1", xlabel="E [keV]", ylabel="Λ [-]")

    myplt = Plots.plot(myplt1, myplt2, myplt3,
                    myplt4, myplt5, myplt6,
                    myplt7, myplt8, myplt9,
                    layout=(3,3),size=(1200,1200),dpi=200)
    display(myplt)
    png(myplt, folderpath_OWCF*"tests/outputs/dependencies_test_F_os_And_F_COM")
end
###------------------------------------------------------------------------------------------------###

############---------------------------------------------------------------------------------------###
if plot_test_results
    test_F_os_E = dropdims(sum((dpm*dRm) .*test_F_os_3D,dims=(2,3)),dims=(2,3))
    test_F_COM_E = dropdims(sum((dLambda*dPphi_n) .*(test_F_COM_3D[:,:,:,1] + test_F_COM_3D[:,:,:,2]),dims=(2,3)),dims=(2,3))
    myplt = Plots.plot(test_E_array, test_F_os_E, label="(E,pm,Rm)",linewidth=3)
    myplt = Plots.plot!(test_E_array_os2com, test_F_COM_E, label="(E,Λ,Pϕ_n;σ)",linewidth=1.5)
    myplt = Plots.plot!(xlabel="E [keV]",title="f(E)",ylabel="f [keV^-1]",dpi=100)
    display(myplt)
    png(myplt, folderpath_OWCF*"tests/outputs/dependencies_test_F_os_vs_F_COM_energy_comp")
end
###------------------------------------------------------------------------------------------------###

############---------------------------------------------------------------------------------------###
test_Lambda_array = collect(range(0.0,stop=1.5,length=29))
test_Pphi_n_array = collect(range(-1.0,stop=1.0,length=30))

test_topoMap_os = getOSTopoMap(M, test_E_array, test_pm_array, test_Rm_array; wall=wall, distinguishLost=true)
test_topoMap_COM_coGoing = getCOMTopoMap(M, test_E_array, test_Lambda_array, test_Pphi_n_array; sigma=+1)
test_topoMap_COM_ctGoing = getCOMTopoMap(M, test_E_array, test_Lambda_array, test_Pphi_n_array; sigma=-1)
test_topoMap_COM_os2COM, test_E_array_os2com, test_Lambda_array_os2com, test_Pphi_n_array_os2com = os2COM(M, test_topoMap_os, test_E_array, test_pm_array, test_Rm_array, test_FI_species; isTopoMap=true, wall=wall) # topoMap
###------------------------------------------------------------------------------------------------###

############---------------------------------------------------------------------------------------###
if plot_test_results
    pm_array_ext = vcat(test_pm_array[1]-(diff(test_pm_array))[1],test_pm_array) # Extend pm_array one row below
    Lambda_array_ext = vcat(test_Lambda_array[1]-(diff(test_Lambda_array))[1],test_Lambda_array) # Extend test_Lambda_array one row below
    Lambda_array_ext_os2com = vcat(test_Lambda_array_os2com[1]-(diff(test_Lambda_array_os2com))[1],test_Lambda_array_os2com) # Extend test_Lambda_array_os2com one row below

    iE = 1 # First energy slice
    iE_os2com = iE # First (1) is always first (1)

    topoMap_ext_1 = vcat(vcat([1,2,3,4,5,6,7,8,9],ones(length(test_Rm_array)-9))', test_topoMap_os[iE,:,:])
    topoMap_ext_2 = vcat(vcat([1,2,3,4,5,6,7,8,9],ones(length(test_Pphi_n_array)-9))', test_topoMap_COM_coGoing[iE,:,:])
    topoMap_ext_3 = vcat(vcat([1,2,3,4,5,6,7,8,9],ones(length(test_Pphi_n_array_os2com)-9))', test_topoMap_COM_os2COM[iE_os2com,:,:,2]) # coGoing
    topoMap_ext_4 = vcat(vcat([1,2,3,4,5,6,7,8,9],ones(length(test_Pphi_n_array)-9))', test_topoMap_COM_ctGoing[iE,:,:])
    topoMap_ext_5 = vcat(vcat([1,2,3,4,5,6,7,8,9],ones(length(test_Pphi_n_array_os2com)-9))', test_topoMap_COM_os2COM[iE_os2com,:,:,1]) # ctGoing

    myplt1 = Plots.heatmap(test_Rm_array, pm_array_ext, topoMap_ext_1, color=:Set1_9, legend=false, xlabel="Rm [m]", ylabel="pm", title="E: $(round(test_E_array[iE],digits=3)) keV", ylims=extrema(test_pm_array), xlims=extrema(test_Rm_array))
    myplt2 = Plots.heatmap(test_Pphi_n_array, Lambda_array_ext, topoMap_ext_2, color=:Set1_9, legend=false, xlabel="Pϕ_n [-]", ylabel="Λ [-]", title="getCOMTopoMap() σ=+1", ylims=extrema(test_Lambda_array), xlims=extrema(test_Pphi_n_array))
    myplt3 = Plots.heatmap(test_Pphi_n_array_os2com, Lambda_array_ext_os2com, topoMap_ext_3, color=:Set1_9, legend=false, xlabel="Pϕ_n [-]", ylabel="Λ [-]", title="os2COM() σ=+1", ylims=extrema(test_Lambda_array), xlims=extrema(test_Pphi_n_array))
    myplt4 = Plots.heatmap(test_Pphi_n_array, Lambda_array_ext, topoMap_ext_4, color=:Set1_9, legend=false, xlabel="Pϕ_n [-]", ylabel="Λ [-]", title="getCOMTopoMap() σ=-1", ylims=extrema(test_Lambda_array), xlims=extrema(test_Pphi_n_array))
    myplt5 = Plots.heatmap(test_Pphi_n_array_os2com, Lambda_array_ext_os2com, topoMap_ext_5, color=:Set1_9, legend=false, xlabel="Pϕ_n [-]", ylabel="Λ [-]", title="os2COM() σ=-1", ylims=extrema(test_Lambda_array), xlims=extrema(test_Pphi_n_array))

    myplt = Plots.plot(myplt1,Plots.plot(myplt2,myplt3,layout=(1,2)),Plots.plot(myplt4,myplt5,layout=(1,2)),layout=(3,1),size=(800,1200),dpi=200)
    display(myplt)
    png(myplt,folderpath_OWCF*"tests/outputs/dependencies_test_topoMap_os_1of3")
end
###------------------------------------------------------------------------------------------------###

############---------------------------------------------------------------------------------------###
if plot_test_results
    iE = Int64(round(length(test_E_array)/2)) # Middle energy slice
    iE_os2com = Int64(round(length(test_E_array_os2com)/2)) # Should be same as iE

    topoMap_ext_1 = vcat(vcat([1,2,3,4,5,6,7,8,9],ones(length(test_Rm_array)-9))', test_topoMap_os[iE,:,:])
    topoMap_ext_2 = vcat(vcat([1,2,3,4,5,6,7,8,9],ones(length(test_Pphi_n_array)-9))', test_topoMap_COM_coGoing[iE,:,:])
    topoMap_ext_3 = vcat(vcat([1,2,3,4,5,6,7,8,9],ones(length(test_Pphi_n_array_os2com)-9))', test_topoMap_COM_os2COM[iE_os2com,:,:,2]) # coGoing
    topoMap_ext_4 = vcat(vcat([1,2,3,4,5,6,7,8,9],ones(length(test_Pphi_n_array)-9))', test_topoMap_COM_ctGoing[iE,:,:])
    topoMap_ext_5 = vcat(vcat([1,2,3,4,5,6,7,8,9],ones(length(test_Pphi_n_array_os2com)-9))', test_topoMap_COM_os2COM[iE_os2com,:,:,1]) # ctGoing

    myplt1 = Plots.heatmap(test_Rm_array, pm_array_ext, topoMap_ext_1, color=:Set1_9, legend=false, xlabel="Rm [m]", ylabel="pm", title="E: $(round(test_E_array[iE],digits=3)) keV", ylims=extrema(test_pm_array), xlims=extrema(test_Rm_array))
    myplt2 = Plots.heatmap(test_Pphi_n_array, Lambda_array_ext, topoMap_ext_2, color=:Set1_9, legend=false, xlabel="Pϕ_n [-]", ylabel="Λ [-]", title="getCOMTopoMap() σ=+1", ylims=extrema(test_Lambda_array), xlims=extrema(test_Pphi_n_array))
    myplt3 = Plots.heatmap(test_Pphi_n_array_os2com, Lambda_array_ext_os2com, topoMap_ext_3, color=:Set1_9, legend=false, xlabel="Pϕ_n [-]", ylabel="Λ [-]", title="os2COM() σ=+1", ylims=extrema(test_Lambda_array), xlims=extrema(test_Pphi_n_array))
    myplt4 = Plots.heatmap(test_Pphi_n_array, Lambda_array_ext, topoMap_ext_4, color=:Set1_9, legend=false, xlabel="Pϕ_n [-]", ylabel="Λ [-]", title="getCOMTopoMap() σ=-1", ylims=extrema(test_Lambda_array), xlims=extrema(test_Pphi_n_array))
    myplt5 = Plots.heatmap(test_Pphi_n_array_os2com, Lambda_array_ext_os2com, topoMap_ext_5, color=:Set1_9, legend=false, xlabel="Pϕ_n [-]", ylabel="Λ [-]", title="os2COM() σ=-1", ylims=extrema(test_Lambda_array), xlims=extrema(test_Pphi_n_array))

    myplt = Plots.plot(myplt1,Plots.plot(myplt2,myplt3,layout=(1,2)),Plots.plot(myplt4,myplt5,layout=(1,2)),layout=(3,1),size=(800,1200),dpi=200)
    display(myplt)
    png(myplt,folderpath_OWCF*"tests/outputs/dependencies_test_topoMap_os_2of3")
end
###------------------------------------------------------------------------------------------------###

############---------------------------------------------------------------------------------------###
if plot_test_results
    #iE = end # Last energy slice
    #iE_os2com = end # Last (end) is always last (end)

    topoMap_ext_1 = vcat(vcat([1,2,3,4,5,6,7,8,9],ones(length(test_Rm_array)-9))', test_topoMap_os[end,:,:])
    topoMap_ext_2 = vcat(vcat([1,2,3,4,5,6,7,8,9],ones(length(test_Pphi_n_array)-9))', test_topoMap_COM_coGoing[end,:,:])
    topoMap_ext_3 = vcat(vcat([1,2,3,4,5,6,7,8,9],ones(length(test_Pphi_n_array_os2com)-9))', test_topoMap_COM_os2COM[end,:,:,2]) # coGoing
    topoMap_ext_4 = vcat(vcat([1,2,3,4,5,6,7,8,9],ones(length(test_Pphi_n_array)-9))', test_topoMap_COM_ctGoing[end,:,:])
    topoMap_ext_5 = vcat(vcat([1,2,3,4,5,6,7,8,9],ones(length(test_Pphi_n_array_os2com)-9))', test_topoMap_COM_os2COM[end,:,:,1]) # ctGoing

    myplt1 = Plots.heatmap(test_Rm_array, pm_array_ext, topoMap_ext_1, color=:Set1_9, legend=false, xlabel="Rm [m]", ylabel="pm", title="E: $(round(test_E_array[end],digits=3)) keV", ylims=extrema(test_pm_array), xlims=extrema(test_Rm_array))
    myplt2 = Plots.heatmap(test_Pphi_n_array, Lambda_array_ext, topoMap_ext_2, color=:Set1_9, legend=false, xlabel="Pϕ_n [-]", ylabel="Λ [-]", title="getCOMTopoMap() σ=+1", ylims=extrema(test_Lambda_array), xlims=extrema(test_Pphi_n_array))
    myplt3 = Plots.heatmap(test_Pphi_n_array_os2com, Lambda_array_ext_os2com, topoMap_ext_3, color=:Set1_9, legend=false, xlabel="Pϕ_n [-]", ylabel="Λ [-]", title="os2COM() σ=+1", ylims=extrema(test_Lambda_array), xlims=extrema(test_Pphi_n_array))
    myplt4 = Plots.heatmap(test_Pphi_n_array, Lambda_array_ext, topoMap_ext_4, color=:Set1_9, legend=false, xlabel="Pϕ_n [-]", ylabel="Λ [-]", title="getCOMTopoMap() σ=-1", ylims=extrema(test_Lambda_array), xlims=extrema(test_Pphi_n_array))
    myplt5 = Plots.heatmap(test_Pphi_n_array_os2com, Lambda_array_ext_os2com, topoMap_ext_5, color=:Set1_9, legend=false, xlabel="Pϕ_n [-]", ylabel="Λ [-]", title="os2COM() σ=-1", ylims=extrema(test_Lambda_array), xlims=extrema(test_Pphi_n_array))

    myplt = Plots.plot(myplt1,Plots.plot(myplt2,myplt3,layout=(1,2)),Plots.plot(myplt4,myplt5,layout=(1,2)),layout=(3,1),size=(800,1200),dpi=200)
    display(myplt)
    png(myplt,folderpath_OWCF*"tests/outputs/dependencies_test_topoMap_os_3of3")
end
###------------------------------------------------------------------------------------------------###

############---------------------------------------------------------------------------------------###
###------------------------------------------------------------------------------------------------###

############---------------------------------------------------------------------------------------###
###------------------------------------------------------------------------------------------------###

############---------------------------------------------------------------------------------------###
###------------------------------------------------------------------------------------------------###