######################################## run_tests.jl #############################################
# This script will install the OWCF and set all the environmental variables needed

# Script written by Henrik Järleblad. Last maintained 2025-10-07.
###################################################################################################

PYTHON_COMMAND = nothing # Please manually specify this, if you have got an error from previously running this script

t_start = time() # Installation script runtime start

if isnothing(PYTHON_COMMAND)
    s = "q"
    while !(s=="y" || s=="n")
        global s
        print("Do you run Python in your terminal with the 'python' command [y/n]? ")
        s = lowercase(readline())
    end
    if s=="y"
        PYTHON_COMMAND = "python"
    else
        s = "q"
        while !(s=="y" || s=="n")
            global s
            print("Ok. Is it the 'python3' command [y/n]? ")
            s = lowercase(readline())
        end
        if s=="y"
            PYTHON_COMMAND = "python3"
        else
            print("What is it then? (please specify): ")
            PYTHON_COMMAND = readline()
        end
    end
end

print("Do you have 'pip' installed and linked to your $(PYTHON_COMMAND) Python [y/n]? ")
s = lowercase(readline())
if !(s=="y")
    error("Please ensure that you have 'pip' installed and linked to your $(PYTHON_COMMAND) Python before installing the OWCF")
end

s = "q"
while !(s=="y" || s=="n")
    global s
    print("Are you installing the OWCF on a computational cluster [y/n]? ")
    s = lowercase(readline())
end
if s=="y"
    println("---> Restricting the CPU resources for Julia precompilation on clusters to 1 CPU... ")
    ENV["JULIA_NUM_PRECOMPILE_TASKS"]=1
end

############---------------------------------------------------------------------------------------###
folderpath_OWCF = @__DIR__ # Get the folder path to the OWCF folder
using Dates # Load the Dates Julia default package, for printing time and date
###------------------------------------------------------------------------------------------------###

############---------------------------------------------------------------------------------------###
folder_delimiter = "/" # Assume forward slash file explorer delimiter by default (Linux or Apple)
PYTHON_SUFFIX = "" # Assume no suffix by default (Linux and Apple)
if Sys.iswindows()
    println("---> Found Windows operating system. Proceeding with OWCF installation accordingly... ")
    folder_delimiter = "\\" # If Windows, we need to use \\ instead
    PYTHON_SUFFIX = ".exe" # If Windows, we need the .exe suffix for the python executable
end
###------------------------------------------------------------------------------------------------###

############---------------------------------------------------------------------------------------###
date_and_time = split("$(Dates.now())","T")[1]*" at "*split("$(Dates.now())","T")[2][1:5] # Determine the date and time
println("-------------------------------------------------------------------------------- The OWCF installation --------------------------------------------------------------------------------")
println("-------------------------------------------------------------------------------- OWCF folder path: $(folderpath_OWCF)")
println("-------------------------------------------------------------------------------- Current date and time: $(date_and_time)")
println("-------------------------------------------------------------------------------- Starting the OWCF installation... ")
###------------------------------------------------------------------------------------------------###

############---------------------------------------------------------------------------------------###
println("- Changing directory to the OWCF folder... ")
cd(folderpath_OWCF)

println("- Loading necessary Julia packages and activating the OWCF environment... ")
using Pkg
Pkg.activate(".")

println("- Determining the path to your Python executable... ")
try
    Base.run(`$(PYTHON_COMMAND) misc$(folder_delimiter)save_path_to_folder_of_python_executable.py`)
catch
    error("Could not run Python via the '$(PYTHON_COMMAND)' command. Please try manually specifying the PYTHON_COMMAND input variable in the install_OWCF.jl script, and re-run the 'julia install_OWCF.jl' command")
end
f = open("PYTHON_FOLDERPATH.txt", "r")
PYTHON_FOLDERPATH = readline(f)
close(f)
println("---> Python executable path is $(PYTHON_FOLDERPATH)$(folder_delimiter)$(PYTHON_COMMAND)$(PYTHON_SUFFIX)")
rm("$(folderpath_OWCF)$(folder_delimiter)PYTHON_FOLDERPATH.txt")
###------------------------------------------------------------------------------------------------###

############---------------------------------------------------------------------------------------###
println("- Setting the PYTHON environment variable... ")
PYTHON_EXECUTABLE = "$(PYTHON_FOLDERPATH)$(folder_delimiter)$(PYTHON_COMMAND)$(PYTHON_SUFFIX)"
ENV["PYTHON"] = PYTHON_EXECUTABLE
###------------------------------------------------------------------------------------------------###

############---------------------------------------------------------------------------------------###
println("- Instantiating and precompiling the OWCF (takes a looong time)... ")
Pkg.instantiate()
###------------------------------------------------------------------------------------------------###

############---------------------------------------------------------------------------------------###
t_end = time() # Test script runtime end
println("-------------------------------------------------------------------------------- End of the OWCF installation")
println("-------------------------------------------------------------------------------- Total runtime: $(Int64(divrem(t_end-t_start,60)[1])) minutes $(Int64(round(divrem(t_end-t_start,60)[2]))) seconds")
date_and_time = split("$(Dates.now())","T")[1]*" at "*split("$(Dates.now())","T")[2][1:5]
println("-------------------------------------------------------------------------------- Current date and time: $(date_and_time)")
println("-------------------------------------------------------------------------------- OWCF folder path: $(folderpath_OWCF)")
println("-------------------------------------------------------------------------------- The OWCF installation --------------------------------------------------------------------------------")
###------------------------------------------------------------------------------------------------###

println("")
println("(Please run 'julia $(folderpath_OWCF)$(folder_delimiter)tests$(folder_delimiter)run_tests.jl to ensure that the OWCF has been installed correctly)")