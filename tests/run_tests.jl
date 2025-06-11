######################################## run_tests.jl #############################################
#
# PLEASE NOTE!! THIS SCRIPT AND RELATED SCRIPTS ARE CURRENTLY UNDER DEVELOPMENT!!!
# PLEASE NOTE!! THIS SCRIPT AND RELATED SCRIPTS ARE CURRENTLY UNDER DEVELOPMENT!!!
# PLEASE NOTE!! THIS SCRIPT AND RELATED SCRIPTS ARE CURRENTLY UNDER DEVELOPMENT!!!

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

### Inputs:
# plot_test_results - If true, test results will be plotted and saved in .png file format - Bool
# terminate_when_first_error - If true, as soon as an error is encountered, the error and stacktrace 
#                              will be printed and the run_tests.jl script will terminate. If there 
#                              are parts of the OWCF that have not yet been tested, these will be 
#                              ignored, i.e. not tested - Bool
# clear_test_outputs_folder_when_done - If true, the OWCF/tests/outputs/ folder will be cleared and 
#                                       deleted when the run_tests.jl script has finished all tests. 
#                                       'false' can be useful if test results are to be manually 
#                                       inspected in detail - Bool
# VERY_VERBOSE - If true, all print statements from all individual tests will be printed to terminal.
#                WARNING! The VERY_VERBOSE input variable SHOULD BE SET TO false, unless you are a 
#                developer and debugging.

### Outputs:
# None, after run_tests.jl has completed. However, temporary output files are saved in the OWCF/
# tests/outputs/ folder.

### Other:
# Please do NOT change anything (!) in the OWCF/tests/ folder. The final results of the run_tests.jl 
# script should be printed to terminal only. As explained above, if you are a developer, this
# printed output should be included as a screenshot when creating a pull-request at the Github repo
# https://github.com/JuliaFusion/OWCF/pulls.
#
# Finally, please note that the run_tests.jl script will produce several data files (and .png figure 
# files, if plot_test_results==true). These will be saved in the OWCF/tests/outputs/ folder. The 
# total size of all data files and .png figure files will be in the range of tens of megabytes.

# Script written by Henrik Järleblad. Last mainted 2025-06-11.
###################################################################################################

# Inputs. To be switched freely between 'true' and 'false'
plot_test_results = false # If set to true, results of the individual tests will be plotted and saved as .png files in the OWCF/tests/outputs/ folder
terminate_when_first_error = false # If set to true, the run_tests.jl script will print the first error it encounters and skip the rest of the testing process
clear_test_outputs_folder_when_done = true # Should be set to true, be default
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

folderpath_OWCF = reduce(*,map(x-> "/"*x,split(@__DIR__,"/")[2:end-1]))*"/" # We know that the run_tests.jl file is located in the OWCF/tests/ folder. Deduce the full OWCF folder path from that information
cd(folderpath_OWCF); using Pkg; Pkg.activate(".") # Navigate to the OWCF folder, activate the OWCF environment
using Distributed # To enable test files of accessing the 'folderpath_OWCF' variable value correctly
using Dates # To enable date and time functions use
using FileIO # For writing error files

oldstdout = stdout # To be able to suppress prints from the test scripts, but not errors
oldstderr = stderr # To be able to suppress warnings from the test scripts, but not errors
SUPPRESS_SUBTEST_PRINT = !VERY_VERBOSE
###------------------------------------------------------------------------------------------------###

############---------------------------------------------------------------------------------------###
# Checking that the OWCF/tests/outputs/ folder is empty
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
        error("The run_tests.jl script cannot be executed unless the $(folderpath_OWCF)tests/outputs/ folder is empty. Please empty the folder manually, and re-start run_tests.jl.")
    end
end
println("")

# Checking that the OWCF/tests/errors/ folder is empty
if !isdir(folderpath_OWCF*"tests/errors/")
    print("The folder $(folderpath_OWCF)tests/errors/ does not exist. Creating... ")
    mkdir(folderpath_OWCF*"tests/errors")
    println("ok!")
end
if !isempty(readdir(folderpath_OWCF*"tests/errors/"))
    num_o_files = length(readdir(folderpath_OWCF*"tests/errors/")); ff = num_o_files==1 ? "file" : "files"
    println("The $(folderpath_OWCF)tests/errors/ folder was not empty ($(num_o_files) $(ff)).")
    s = "q"
    while !(lowercase(s)=="y" || lowercase(s)=="n")
        global s
        print("Remove all files in the $(folderpath_OWCF)tests/errors/ folder (y/n)? ")
        s = readline()
    end
    if lowercase(s)=="y"
        print("Removing all files in the $(folderpath_OWCF)tests/errors/ folder... ")
        error_files = readdir(folderpath_OWCF*"tests/errors/")
        for error_file in error_files
            rm(folderpath_OWCF*"tests/errors/"*error_file)
        end
        println("Done!")
    else
        error("The run_tests.jl script cannot be executed unless the $(folderpath_OWCF)tests/errors/ folder is empty. Please empty the folder manually, and re-start run_tests.jl.")
    end
end
println("")
###------------------------------------------------------------------------------------------------###

############---------------------------------------------------------------------------------------###
date_and_time = split("$(Dates.now())","T")[1]*" at "*split("$(Dates.now())","T")[2][1:5]
test_list = map(x-> "$(split(x,".")[1])", readdir(folderpath_OWCF*"tests/start_files/")) # All tests are in the OWCF/tests/start_files/ folder. Remove the ".jl" file extension
NUMBER_OF_TESTS = length(test_list) # The total number of tests
err_dict = Dict() # A dictionary to keep track of the test error file paths
test_prog = 0 # An integer to keep track of testing progress
###------------------------------------------------------------------------------------------------###

############---------------------------------------------------------------------------------------###
println("-------------------------------------------------------------------------------- The OWCF tests --------------------------------------------------------------------------------")
println("-------------------------------------------------------------------------------- OWCF folder path: $(folderpath_OWCF)")
println("-------------------------------------------------------------------------------- Current date and time: $(date_and_time)")
println("-------------------------------------------------------------------------------- Starting the OWCF tests... ")
###------------------------------------------------------------------------------------------------###

############---------------------------------------------------------------------------------------###
for test in test_list
    global test_prog # Declare test_prog as the global test_prog variable defined already before (outside of) this for-loop
    err_file = "$(folderpath_OWCF)tests/errors/$(test).err" # No global scope declaration needed, since err_file String will not be changed
    err_IO = open(err_file,"w") # No global scope declaration needed, since the err_IO channel will not be changed (only used)

    println("- Running $(test).jl (total run_tests progress: $(Int64(round(100*test_prog/NUMBER_OF_TESTS)))%)... ")

    try
        SUPPRESS_SUBTEST_PRINT && redirect_stdout(devnull) # Re-direct prints to null
        SUPPRESS_SUBTEST_PRINT && redirect_stderr(err_IO) # Re-direct warnings and errors to error output file
        # By using Base.run(), we effectively reset the namespace and unload all Julia packages every time before each test start file is run
        Base.run(`julia $(folderpath_OWCF*"tests/start_files/$(test).jl") plot_test_results $(plot_test_results)`)
    catch e
        SUPPRESS_SUBTEST_PRINT && redirect_stderr(oldstderr) # Re-direct warnings and errors back to old I/O channel (i.e. the default, visible terminal)
        global err_dict # Declare err_dict as the err_dict variable from the global scope (outside of the for-loop and the try-catch statement)
        @warn "The test $(test) unfortunately resulted in an error. Please examine the error and stack trace in $(err_file)"
        err_dict["$(test)"] = "$(err_file)" # Save the error file path to the err_dict Dictionary, with the test start file as key
        terminate_when_first_error && redirect_stdout(oldstdout) # Re-direct prints back to old I/O channel (i.e. the default, visible channel)
        terminate_when_first_error && close(err_IO) # If we break the for-loop immidiately below this line, close the error file I/O channel
        terminate_when_first_error && break # Break from the 'for test in test_list' for-loop
    end

    close(err_IO) # If terminate_when_first_error==false, we need to close the error file I/O channel here
    if !("$(test)" in keys(err_dict)) # If this test did not produce any errors...
        rm(err_file) # To avoid confusion, remove the error file
    end

    SUPPRESS_SUBTEST_PRINT && redirect_stdout(oldstdout) # Re-direct prints back to old I/O channel (i.e. the default, visible channel)
    test_prog += 1 # Test completed! Onto the next one!
end
println("- All tests completed (total run_tests.jl progress: $(Int64(round(100*test_prog/NUMBER_OF_TESTS)))%)... ")
###------------------------------------------------------------------------------------------------###

############---------------------------------------------------------------------------------------###
# Removing all created files in the OWCF/tests/outputs/ folder, if requested
println("--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------")
if clear_test_outputs_folder_when_done
    print("- Clearing the OWCF/tests/outputs/ folder... ")
    rm(folderpath_OWCF*"tests/outputs/", recursive=true)
    println("Done!")
else
    output_folder_size_in_megabytes = inv(1_000_000)*reduce(+, filesize.("$(folderpath_OWCF)tests/outputs/" .*readdir("$(folderpath_OWCF)tests/outputs/")))
    println("- The $(folderpath_OWCF)tests/outputs/ folder contains all test output files. The total size of all test output files is approx. $(round(output_folder_size_in_megabytes, digits=1)) MB. Please inspect files in folder prior to running run_tests.jl again.")
end
###------------------------------------------------------------------------------------------------###

############---------------------------------------------------------------------------------------###
# Print all errors (if any) and total test result
println("--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------")
println("- Printing test results... ")

tot_number_o_errs = length(keys(err_dict))
println("- Total number of errors: $(tot_number_o_errs)")
if tot_number_o_errs==0
    rm(folderpath_OWCF*"tests/errors/", recursive=true) # If no errors, to avoid confusion, clear and remove the OWCF/tests/errors/ folder
    println("- CONGRATULATIONS! The OWCF can be assumed to function correctly.")
    println("")
    println("---> IF YOU ARE A USER, please go ahead and safely use the OWCF. Enjoy!")
    println("---> IF YOU ARE A DEVELOPER, please attach a screenshot of the test results to your pull-request on Github, if you would like to merge your branch with the main OWCF branch.")
else
    s = tot_number_o_errs > 1 ? "s" : ""
    println("- Please examine the error$(s) and stack trace$(s) in the following file$(s) (can be opened in the same way as .txt files)... ")
    println("")
    for key in keys(err_dict)
        println(" --->  $(key): $(err_dict[key])")
        println(" ------> To re-run only this specific test, please do (with 'b' set to either true or false) 'julia $(folderpath_OWCF*"tests/start_files/$(key).jl") plot_test_results b'")
        println("")
        println("")
    end
    println("")
    println("---> IF YOU ARE A USER, please post a new issue at https://github.com/JuliaFusion/OWCF/issues, attach the test results and describe your situation. Thank you!")
    println("---> IF YOU ARE A DEVELOPER, please investigate the error$(s), fix the error$(s) and re-run run-tests.jl.")
end
###------------------------------------------------------------------------------------------------###

############---------------------------------------------------------------------------------------###
t_end = time() # Test script runtime end
println("-------------------------------------------------------------------------------- End of the OWCF tests")
println("-------------------------------------------------------------------------------- Total runtime: $(Int64(divrem(t_end-t_start,60)[1])) minutes $(Int64(round(divrem(t_end-t_start,60)[2]))) seconds")
println("-------------------------------------------------------------------------------- OWCF folder path: $(folderpath_OWCF)")
println("-------------------------------------------------------------------------------- The OWCF tests --------------------------------------------------------------------------------")
###------------------------------------------------------------------------------------------------###