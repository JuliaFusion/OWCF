######################################## run_tests.jl #############################################

# PLEASE NOTE!! THIS SCRIPTS AND RELATED SCRIPTS ARE CURRENTLY UNDER DEVELOPMENT!!!
# PLEASE NOTE!! THIS SCRIPTS AND RELATED SCRIPTS ARE CURRENTLY UNDER DEVELOPMENT!!!
# PLEASE NOTE!! THIS SCRIPTS AND RELATED SCRIPTS ARE CURRENTLY UNDER DEVELOPMENT!!!

### Description:
# This script will run all the tests of the OWCF. It will run all the start files in the folder 
# OWCF/tests/start_files/ and compare their results with the reference data in the folder 
# OWCF/tests/standards/. The results of the start files are saved in the OWCF/tests/outputs/ folder
# and will be automatically removed when all comparisons are completed. If an output result fail
# to match the standard set by the reference data, the test will fail. A test failing can mean 
# either of the following two things:
#
#   1. If you are a user of the OWCF, but NOT a developer, the installation of the OWCF might 
#      have not completely correctly. Please double-check so that you followed the installation
#      instructions correctly. If you did, there might be a problem with your Julia and/or Python
#      program. If they are ok, the failed tests might possibly, potentially, maybe be due to 
#      a mistake made by the OWCF developers. If so, please send an email to henrikj@dtu.dk or
#      anvalen@fysik.dtu.dk, or create a new issue at https://github.com/JuliaFusion/OWCF/issues.
#   2. If you are a developer of the OWCF, your changes might have broken something. Please 
#      check which test(s) failed and examine how it/they failed. If still unclear, please contact
#      henrikj@dtu.dk or anvalen@fysik.dtu.dk, to discuss why your changes caused tests to fail.
#      PLEASE NOTE! Your pull-request and merge with the OWCF main branch will NOT be approved 
#      unless you include a screenshot of your test results with your pull-request. Your screenshot
#      should clearly show that the execution of 'include("run_tests.jl")' resulted in 100% of the 
#      tests completing successfully.

### Inputs:
# folderpath_OWCF - The path to the OWCF folder on your computer - String

### Outputs:
# 

### Other:
# Please do NOT change anything (!) in the OWCF/tests/ folder. The results of the run_tests.jl 
# script should be printed to terminal only. As explained above, if you are a developer, this
# printed output should be included as a screenshot when creating a pull-request at the Github repo
# https://github.com/JuliaFusion/OWCF/pulls.

# Script written by Henrik Järleblad. Last mainted 2025-01-27.
###################################################################################################

folderpath_OWCF = "/home/henrikj/Codes/OWCF/" # The path to where the OWCF folder is saved on your computer. Remember to finish with '/'

# ------------------------------------------------------------------------------------------------#
# -------------------------------------- Start of tests ------------------------------------------#
# ----------------------- Please do not alter anything below these lines! ------------------------#
# ------------------------------------------------------------------------------------------------#

cd(folderpath_OWCF)
using Pkg
Pkg.activate(".")

using JLD2
using Dates

oldstd = stdout # To be able to surpress prints from the test scripts, but not errors

date_and_time = split("$(Dates.now())","T")[1]*"at"*split("$(Dates.now())","T")[2][1:5]
tot_number_o_tests = 20
println("---------------------------------------- The OWCF tests ----------------------------------------")
println("---------------------------------------- Current date and time: $(date_and_time)")
println("---------------------------------------- OWCF folder path: $(folderpath_OWCF)")
println("---------------------------------------- Starting the OWCF tests... ")
test_prog = 0
# ------------------------------------------------------------------------------------------------#
println("- Running 2D weight function computation test 1 (run_tests progress: $(100*test_prog/tot_number_o_tests) %)... ")
redirect_stdout(devnull); include("start_files/start_calc2DW_test_1.jl"); redirect_stdout(oldstd)
myreffile = jldopen("standards/test_1_standard.jld2",false,false,false,IOStream)
W_ref = myreffile["W"]; close(myfile)
myfile = jldopen("outputs/test_1_output.jld2",false,false,false,IOStream)
W_test = myfile["W"]; close(myfile)

diff = sum(abs.(W_test .- W_ref))
if diff<1.0e-11 # Standard set via empirical studies
    test_1_passed = true
else
    test_1_passed = false
end
test_prog += 1 # First test completed!
# ------------------------------------------------------------------------------------------------#