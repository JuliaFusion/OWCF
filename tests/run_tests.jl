######################################## run_tests.jl #############################################
#
# PLEASE NOTE!! THIS SCRIPTS AND RELATED SCRIPTS ARE CURRENTLY UNDER DEVELOPMENT!!!
# PLEASE NOTE!! THIS SCRIPTS AND RELATED SCRIPTS ARE CURRENTLY UNDER DEVELOPMENT!!!
# PLEASE NOTE!! THIS SCRIPTS AND RELATED SCRIPTS ARE CURRENTLY UNDER DEVELOPMENT!!!
#
### Description:
# This script will run all the tests of the OWCF, to ensure correct functionality. 
# It will run all the start files in the folder OWCF/tests/start_files/. If no errors are produced,
# the OWCF is assumed to function correctly. During testing, the results of the start files 
# are saved in the OWCF/tests/outputs/ folder and will be automatically removed when all tests are 
# completed. A test failing can mean either of the following two things:
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
#
### Inputs:
# folderpath_OWCF - The path to the OWCF folder on your computer - String
# plot_test_results - If true, test results will be plotted and saved in .png file format - Bool
# clear_test_outputs_folder_when_done - If true, the OWCF/tests/outputs/ folder will be cleared when 
#                                  the run_tests.jl script has finished all tests. 'false' can be useful
#                                  if test results are to be manually inspected in detail - Bool
# VERY_VERBOSE - If true, all print statements from all individual tests will be printed to terminal.
#                WARNING! The VERY_VERBOSE input variable SHOULD BE SET TO false, unless you are a 
#                developer and debugging.
#
### Outputs:
# None, after run_tests.jl has completed. However, temporary output files are saved in the OWCF/
# tests/outputs/ folder.
#
### Other:
# Please do NOT change anything (!) in the OWCF/tests/ folder. The results of the run_tests.jl 
# script should be printed to terminal only. As explained above, if you are a developer, this
# printed output should be included as a screenshot when creating a pull-request at the Github repo
# https://github.com/JuliaFusion/OWCF/pulls.

# Script written by Henrik Järleblad. Last mainted 2025-04-25.
###################################################################################################

folderpath_OWCF = "/home/henrikj/Codes/OWCF/" # The path to where the OWCF folder is saved on your computer. Remember to finish with '/'
plot_test_results = true # If set to true, test results will be plotted and saved in .png file format in the OWCF/tests/outputs/ folder
clear_test_outputs_folder_when_done = false # Should be set to true, be default
VERY_VERBOSE = false # Should be set to false, unless developer debugging

###------------------------------------------------------------------------------------------------###
###-------------------------------------- START OF TESTS ------------------------------------------###
###----------------------- PLEASE DO NOT ALTER ANYTHING BELOW THESE LINES! ------------------------###
###------------------------------------------------------------------------------------------------###
###------------------------------------------------------------------------------------------------###
###-------------------------------------- START OF TESTS ------------------------------------------###
###----------------------- PLEASE DO NOT ALTER ANYTHING BELOW THESE LINES! ------------------------###
###------------------------------------------------------------------------------------------------###
###------------------------------------------------------------------------------------------------###
###-------------------------------------- START OF TESTS ------------------------------------------###
###----------------------- PLEASE DO NOT ALTER ANYTHING BELOW THESE LINES! ------------------------###
###------------------------------------------------------------------------------------------------###
###------------------------------------------------------------------------------------------------###
###-------------------------------------- START OF TESTS ------------------------------------------###
###----------------------- PLEASE DO NOT ALTER ANYTHING BELOW THESE LINES! ------------------------###
###------------------------------------------------------------------------------------------------###
###------------------------------------------------------------------------------------------------###
###-------------------------------------- START OF TESTS ------------------------------------------###
###----------------------- PLEASE DO NOT ALTER ANYTHING BELOW THESE LINES! ------------------------###
###------------------------------------------------------------------------------------------------###
###------------------------------------------------------------------------------------------------###
###-------------------------------------- START OF TESTS ------------------------------------------###
###----------------------- PLEASE DO NOT ALTER ANYTHING BELOW THESE LINES! ------------------------###
###------------------------------------------------------------------------------------------------###
###------------------------------------------------------------------------------------------------###
###-------------------------------------- START OF TESTS ------------------------------------------###
###----------------------- PLEASE DO NOT ALTER ANYTHING BELOW THESE LINES! ------------------------###
###------------------------------------------------------------------------------------------------###
###------------------------------------------------------------------------------------------------###
###-------------------------------------- START OF TESTS ------------------------------------------###
###----------------------- PLEASE DO NOT ALTER ANYTHING BELOW THESE LINES! ------------------------###
###------------------------------------------------------------------------------------------------###
###------------------------------------------------------------------------------------------------###
###-------------------------------------- START OF TESTS ------------------------------------------###
###----------------------- PLEASE DO NOT ALTER ANYTHING BELOW THESE LINES! ------------------------###
###------------------------------------------------------------------------------------------------###
###------------------------------------------------------------------------------------------------###
###-------------------------------------- START OF TESTS ------------------------------------------###
###----------------------- PLEASE DO NOT ALTER ANYTHING BELOW THESE LINES! ------------------------###
###------------------------------------------------------------------------------------------------###
###------------------------------------------------------------------------------------------------###
###-------------------------------------- START OF TESTS ------------------------------------------###
###----------------------- PLEASE DO NOT ALTER ANYTHING BELOW THESE LINES! ------------------------###
###------------------------------------------------------------------------------------------------###
###------------------------------------------------------------------------------------------------###
###-------------------------------------- START OF TESTS ------------------------------------------###
###----------------------- PLEASE DO NOT ALTER ANYTHING BELOW THESE LINES! ------------------------###
###------------------------------------------------------------------------------------------------###
###------------------------------------------------------------------------------------------------###
###-------------------------------------- START OF TESTS ------------------------------------------###
###----------------------- PLEASE DO NOT ALTER ANYTHING BELOW THESE LINES! ------------------------###
###------------------------------------------------------------------------------------------------###
###------------------------------------------------------------------------------------------------###
###-------------------------------------- START OF TESTS ------------------------------------------###
###----------------------- PLEASE DO NOT ALTER ANYTHING BELOW THESE LINES! ------------------------###
###------------------------------------------------------------------------------------------------###
###------------------------------------------------------------------------------------------------###
###-------------------------------------- START OF TESTS ------------------------------------------###
###----------------------- PLEASE DO NOT ALTER ANYTHING BELOW THESE LINES! ------------------------###
###------------------------------------------------------------------------------------------------###
###------------------------------------------------------------------------------------------------###
###-------------------------------------- START OF TESTS ------------------------------------------###
###----------------------- PLEASE DO NOT ALTER ANYTHING BELOW THESE LINES! ------------------------###
###------------------------------------------------------------------------------------------------###
###------------------------------------------------------------------------------------------------###
###-------------------------------------- START OF TESTS ------------------------------------------###
###----------------------- PLEASE DO NOT ALTER ANYTHING BELOW THESE LINES! ------------------------###
###------------------------------------------------------------------------------------------------###
###------------------------------------------------------------------------------------------------###
###-------------------------------------- START OF TESTS ------------------------------------------###
###----------------------- PLEASE DO NOT ALTER ANYTHING BELOW THESE LINES! ------------------------###
###------------------------------------------------------------------------------------------------###
###------------------------------------------------------------------------------------------------###
###-------------------------------------- START OF TESTS ------------------------------------------###
###----------------------- PLEASE DO NOT ALTER ANYTHING BELOW THESE LINES! ------------------------###
###------------------------------------------------------------------------------------------------###
###------------------------------------------------------------------------------------------------###
###-------------------------------------- START OF TESTS ------------------------------------------###
###----------------------- PLEASE DO NOT ALTER ANYTHING BELOW THESE LINES! ------------------------###
###------------------------------------------------------------------------------------------------###
###------------------------------------------------------------------------------------------------###
###-------------------------------------- START OF TESTS ------------------------------------------###
###----------------------- PLEASE DO NOT ALTER ANYTHING BELOW THESE LINES! ------------------------###
###------------------------------------------------------------------------------------------------###
###------------------------------------------------------------------------------------------------###
###-------------------------------------- START OF TESTS ------------------------------------------###
###----------------------- PLEASE DO NOT ALTER ANYTHING BELOW THESE LINES! ------------------------###
###------------------------------------------------------------------------------------------------###
###------------------------------------------------------------------------------------------------###
###-------------------------------------- START OF TESTS ------------------------------------------###
###----------------------- PLEASE DO NOT ALTER ANYTHING BELOW THESE LINES! ------------------------###
###------------------------------------------------------------------------------------------------###
###------------------------------------------------------------------------------------------------###
###-------------------------------------- START OF TESTS ------------------------------------------###
###----------------------- PLEASE DO NOT ALTER ANYTHING BELOW THESE LINES! ------------------------###
###------------------------------------------------------------------------------------------------###
###------------------------------------------------------------------------------------------------###
###-------------------------------------- START OF TESTS ------------------------------------------###
###----------------------- PLEASE DO NOT ALTER ANYTHING BELOW THESE LINES! ------------------------###
###------------------------------------------------------------------------------------------------###
###------------------------------------------------------------------------------------------------###
###-------------------------------------- START OF TESTS ------------------------------------------###
###----------------------- PLEASE DO NOT ALTER ANYTHING BELOW THESE LINES! ------------------------###
###------------------------------------------------------------------------------------------------###

############---------------------------------------------------------------------------------------###
t_start = time() # Test script runtime start

cd(folderpath_OWCF); using Pkg; Pkg.activate(".") # Navigate to the OWCF folder, activate the OWCF environment
using Distributed # To enable test files of accessing the 'folderpath_OWCF' variable value correctly
using Dates # To enable date and time functions use
using JLD2 # To enable reading of test start files in .jld2 file format 
plot_test_results && (using Plots) # If test results are to be plotted, the Plots.jl package needs to be loaded
oldstd = stdout # To be able to suppress prints from the test scripts, but not warnings and errors
SUPPRESS_SUBTEST_PRINT = !VERY_VERBOSE

# Checking so that the OWCF/tests/outputs/ folder is empty
if !isdir(folderpath_OWCF*"tests/outputs/")
    print("The folder $(folderpath_OWCF)tests/outputs/ does not exist. Creating... ")
    mkdir(folderpath_OWCF*"tests/outputs")
    println("ok!")
end
if !isempty(readdir(folderpath_OWCF*"tests/outputs/"))
    num_o_files = length(readdir(folderpath_OWCF*"tests/outputs/")); ff = num_o_files==1 ? "file" : "files"
    println("The $(folderpath_OWCF)tests/outputs/ folder was not empty ($(num_o_files) $(ff)).")
    s = "q"
    while !(lowercase(s)=="y" || lowercase(s)=="n")
        global s
        print("Remove all files in the $(folderpath_OWCF)tests/outputs/ folder (y/n)? ")
        s = readline()
    end
    if lowercase(s)=="y"
        print("Removing all files in the $(folderpath_OWCF)tests/outputs/ folder... ")
        output_files = readdir(folderpath_OWCF*"tests/outputs/")
        for output_file in output_files
            rm(folderpath_OWCF*"tests/outputs/"*output_file)
        end
        println("Done!")
    else
        error("The run_tests.jl script cannot be executed unless the $(folderpath_OWCF)tests/outputs/ is empty. Please empty it manually, and re-start run_tests.jl.")
    end
end

date_and_time = split("$(Dates.now())","T")[1]*" at "*split("$(Dates.now())","T")[2][1:5]
test_list = readdir(folderpath_OWCF*"tests/start_files/") # All tests are in the OWCF/tests/start_files/ folder
NUMBER_OF_TESTS = length(test_list) # The total number of tests
test_prog = 0 # An integer to keep track of testing progress
err_dict = Dict() # A dictionary to keep track of which tests threw errors, and what the errors were
###------------------------------------------------------------------------------------------------###

############---------------------------------------------------------------------------------------###
println("-------------------------------------------------------------------------------- The OWCF tests --------------------------------------------------------------------------------")
println("-------------------------------------------------------------------------------- OWCF folder path: $(folderpath_OWCF)")
println("-------------------------------------------------------------------------------- Current date and time: $(date_and_time)")
println("-------------------------------------------------------------------------------- Starting the OWCF tests... ")
###------------------------------------------------------------------------------------------------###

############---------------------------------------------------------------------------------------###
for test in test_list
    global test_prog

    println("- Running $(test) (total run_tests.jl progress: $(round(100*test_prog/NUMBER_OF_TESTS,digits=3)) %)... ")
    SUPPRESS_SUBTEST_PRINT && redirect_stdout(devnull)

    try
        include(folderpath_OWCF*"tests/start_files/$(test)")
    catch e
        global err_dict
        @warn "The test $(test) unfortunately resulted in an error. Please examine error specification when test results are printed."
        err_dict["$(test)"] = "$(e)"
    end

    SUPPRESS_SUBTEST_PRINT && redirect_stdout(oldstd)
    test_prog += 1 # Test completed!
end
println("- All tests completed (total run_tests.jl progress: $(round(100*test_prog/NUMBER_OF_TESTS,digits=3)) %)... ")
###------------------------------------------------------------------------------------------------###

############---------------------------------------------------------------------------------------###
# Removing all created files from the OWCF/tests/outputs/ folder
println("--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------")
if clear_test_outputs_folder_when_done
    print("- Clearing the OWCF/tests/outputs/ folder... ")
    output_files = readdir(folderpath_OWCF*"tests/outputs/")
    for output_file in output_files
        rm(folderpath_OWCF*"tests/outputs/"*output_file)
    end
    println("Done!")
else
    println("- The $(folderpath_OWCF)tests/outputs/ folder contains all test output files. Please inspect and clear folder prior to running run_tests.jl again.")
end
###------------------------------------------------------------------------------------------------###

############---------------------------------------------------------------------------------------###
# Print all errors (if any) and total test result
println("--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------")
println("- Printing test results... ")

tot_number_o_errs = length(keys(err_dict))
println("- Total number of errors: $(tot_number_o_errs)")
if tot_number_o_errs==0
    println("- CONGRATULATIONS! The OWCF can be assumed to function correctly.")
    println("")
    println("---> IF YOU ARE A USER, please go ahead and safely use the OWCF. Enjoy!")
    println("---> IF YOU ARE A DEVELOPER, please attach a screenshot of the test results to your pull-request on Github, if you would like to merge your branch with the main OWCF branch.")
else
    println("- Printing error specification(s)... ")
    for key in keys(err_dict)
        println(" --->  $(key): $(err_dict[key])")
        println("")
    end
    println("")
    println("---> IF YOU ARE A USER, please post a new issue at https://github.com/JuliaFusion/OWCF/issues, attach the test results and describe your situation. Thank you!")
    println("---> IF YOU ARE A DEVELOPER, please investigate the error(s) (the specific OWCF scripts and lines of code (integers) are specified above), fix the errors and re-run run-tests.jl.")
end
###------------------------------------------------------------------------------------------------###

############---------------------------------------------------------------------------------------###
t_end = time() # Test script runtime end
println("-------------------------------------------------------------------------------- End of the OWCF tests")
println("-------------------------------------------------------------------------------- Total runtime: $(round(t_end-t_start)) seconds")
println("-------------------------------------------------------------------------------- OWCF folder path: $(folderpath_OWCF)")
println("-------------------------------------------------------------------------------- The OWCF tests --------------------------------------------------------------------------------")
###------------------------------------------------------------------------------------------------###