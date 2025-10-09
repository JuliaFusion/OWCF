######################################## run_tests.jl #############################################
# This script will install the OWCF and set all the environmental variables needed

# Script written by Henrik Järleblad. Last maintained 2025-10-09.
###################################################################################################

PYTHON_COMMAND = nothing # Please manually specify this, if you have got an error from previously running this script

t_start = time() # Installation script runtime start

if isnothing(PYTHON_COMMAND)
    s = "0"
    while !(s=="y" || s=="n")
        global s
        print("Do you run Python in your terminal with the 'python' command [y/n]? ")
        s = lowercase(readline())
    end
    if s=="y"
        PYTHON_COMMAND = "python"
    else
        s = "0"
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

s = "0"
while !(s=="y" || s=="n")
    print("Are you running your $PYTHON_COMMAND Python within a virtual environment [y/n]? ")
    s = lowercase(readline())
end
if s=="y"
    s == "0"
    while !(s=="y" || s=="n")
        print("Did you create your Python virtual environment with conda [y/n]? ")
        s = lowercase(readline())
    end
    if s=="y"
        @warn "The OWCF uses the PyCall.jl Julia package to connect Julia with Python. Currently, this Julia package cannot connect Julia to a Python virtual environment created with conda. After installation, correct functionality of the OWCF cannot be guaranteed."
        s = "0"
        while !(s=="y" || s=="n")
            print("Do you want to continue the installation anyway [y/n]? ")
            s = lowercase(readline())
        end
        if s=="n"
            println("The OWCF installation will now be terminated. Please try to install the OWCF by providing a Python executable not within a conda virtual environment.")
            exit()
        end
    end
    s = "0"
    while !(s=="y" || s=="n")
        print("Are you attempting to install the OWCF with your virtual environment currently activated [y/n]? ")
        s = lowercase(readline())
    end
    if s=="n"
        error("To be able to install the OWCF, Julia needs to connect to Python. However, if you want to connect Julia to a Python virtual environment, you need to activate your virtual environment before you install the OWCF. Please remember, you will also need to activate your virtual environment every time before you use the OWCF.")
    end
end

s = "0"
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

println("- Activating the OWCF environment... ")
using Pkg
Pkg.activate(".")

println("- Determining the path to your Python executable... ")
if Sys.iswindows()
    script = """
    import sys
    if hasattr(sys, "base_exec_prefix"):
        sys.stdout.write(sys.base_exec_prefix)
    else:
        sys.stdout.write(sys.exec_prefix)
    """
else
    script = """
    import sys
    if hasattr(sys, "base_exec_prefix"):
        sys.stdout.write(sys.base_prefix)
        sys.stdout.write(":")
        sys.stdout.write(sys.base_exec_prefix)
    else:
        sys.stdout.write(sys.prefix)
        sys.stdout.write(":")
        sys.stdout.write(sys.exec_prefix)
    """
end
PYTHON_FOLDERPATH = read(Base.run(`$PYTHON_COMMAND -c $script`), String) # Run the script with Python. Return the folder as a String
# Take care of Windows maybe using the folder "Program Files", which we cannot execute.
PYTHON_FOLDERPATH = replace(PYTHON_FOLDERPATH, ("Program Files" => "Program"))
# However, by default, Windows interprets commands like "C:\Program\..." as equivalent to 
# "C:\Program Files\...". So, we simply replace the "Program Files" part of the PYTHON_FOLDERPATH 
# with "Program"
PYTHON_EXECUTABLEPATH = "$(PYTHON_FOLDERPATH)$(folder_delimiter)$(PYTHON_COMMAND)$(PYTHON_SUFFIX)"
println("---> Python executable path is $PYTHON_EXECUTABLEPATH")
###------------------------------------------------------------------------------------------------###

############---------------------------------------------------------------------------------------###
println("- Setting the PYTHON environment variable... ")
ENV["PYTHON"] = PYTHON_EXECUTABLEPATH # Set the PYTHON environment variable
###------------------------------------------------------------------------------------------------###

############---------------------------------------------------------------------------------------###
println("- Instantiating and precompiling the OWCF (takes a looong time)... ")
Pkg.instantiate()
println("---> The OWCF has been successfully instantiated and all Julia packages have been precompiled.")
println("---> After the installation, please run 'julia $(folderpath_OWCF)$(folder_delimiter)tests$(folder_delimiter)run_tests.jl to ensure that the OWCF has been installed correctly.")
###------------------------------------------------------------------------------------------------###

############---------------------------------------------------------------------------------------###
t_end = time() # Test script runtime end
println("")
println("-------------------------------------------------------------------------------- End of the OWCF installation")
println("-------------------------------------------------------------------------------- Total runtime: $(Int64(divrem(t_end-t_start,60)[1])) minutes $(Int64(round(divrem(t_end-t_start,60)[2]))) seconds")
date_and_time = split("$(Dates.now())","T")[1]*" at "*split("$(Dates.now())","T")[2][1:5]
println("-------------------------------------------------------------------------------- Current date and time: $(date_and_time)")
println("-------------------------------------------------------------------------------- OWCF folder path: $(folderpath_OWCF)")
println("-------------------------------------------------------------------------------- The OWCF installation --------------------------------------------------------------------------------")
###------------------------------------------------------------------------------------------------###