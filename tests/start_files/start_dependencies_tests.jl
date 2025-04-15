# A script to test a lot of the tools in the OWCF/extra/dependencies.jl script
include(folderpath_OWCF*"extra/dependencies.jl")
plot_results = plot_test_results # SET TO TRUE, VIA THE plot_test_results INPUT VARIABLE IN THE OWCF/tests/run_tests.jl SCRIPT

############---------------------------------------------------------------------------------------###
μ = [100.0, 0.6, 3.3]
σ = [25.0, 0.1, 0.1]
myGauss = gaussian(μ, σ; mx=[200.0, 1.0, 3.8], mn=[0.0, -1.0, 3.0], n=[30,61,64], floor_level=0.001)
if plot_results
    date_and_time = split("$(Dates.now())","T")[1]*"at"*split("$(Dates.now())","T")[2][1:5]
    abscissa1 = collect(range(0.0,stop=200.0,length=30)); d1 = diff(abscissa1)[1]
    abscissa2 = collect(range(-1.0,stop=1.0,length=61)); d2 = diff(abscissa2)[1]
    abscissa3 = collect(range(3.0,stop=3.8,length=64)); d3 = diff(abscissa3)[1]
    myGauss1 = (d1) .*dropdims(sum(myGauss,dims=1),dims=1)
    myGauss2 = (d2) .*dropdims(sum(myGauss,dims=2),dims=2)
    myGauss3 = (d3) .*dropdims(sum(myGauss,dims=3),dims=3)
    myplt1 = Plots.heatmap(abscissa2,abscissa3,transpose(myGauss1),dpi=600,title="gaussian() $(date_and_time)",titlefontsize=12)
    myplt2 = Plots.heatmap(abscissa1,abscissa3,transpose(myGauss2),dpi=600)
    myplt3 = Plots.heatmap(abscissa1,abscissa2,transpose(myGauss3),dpi=600)
    myplt = Plots.plot(myplt1,myplt2,myplt3,layout=(3,1),size=(400,1200))
    png(myplt,folderpath_OWCF*"tests/outputs/dependencies_gaussian")
end
###------------------------------------------------------------------------------------------------###

############---------------------------------------------------------------------------------------###
# CONTINUE CODING HERE!!!
# CONTINUE CODING HERE!!!
# CONTINUE CODING HERE!!!
###------------------------------------------------------------------------------------------------###