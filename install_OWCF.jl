######################################## run_tests.jl #############################################
# This script will install the OWCF and set all the environmental variables needed

# Script written by Henrik Järleblad. Last maintained 2025-10-09.
###################################################################################################

PYTHON_COMMAND = nothing # Please manually specify this, if you have got an error from previously running this script

t_start = time() # Installation script runtime start

if isnothing(PYTHON_COMMAND)
    user_reply = "0"
    while !(user_reply=="y" || user_reply=="n")
        global user_reply
        print("Do you run Python in your terminal with the 'python' command [y/n]? ")
        user_reply = lowercase(readline())
    end
    if user_reply=="y"
        PYTHON_COMMAND = "python"
    else
        user_reply = "0"
        while !(user_reply=="y" || user_reply=="n")
            global user_reply
            print("Ok. Is it the 'python3' command [y/n]? ")
            user_reply = lowercase(readline())
        end
        if user_reply=="y"
            PYTHON_COMMAND = "python3"
        else
            print("What is it then? (please specify): ")
            PYTHON_COMMAND = readline()
        end
    end
end

user_reply = "0"
while !(user_reply=="y" || user_reply=="n")
    global user_reply
    print("Are you running your $PYTHON_COMMAND Python within a virtual environment [y/n]? ")
    user_reply = lowercase(readline())
end
if user_reply=="y"
    user_reply == "0"
    while !(user_reply=="y" || user_reply=="n")
        global user_reply
        print("Did you create your Python virtual environment with conda [y/n]? ")
        user_reply = lowercase(readline())
    end
    if user_reply=="y"
        @warn "The OWCF uses the PyCall.jl Julia package to connect Julia with Python. Currently, this Julia package cannot connect Julia to a Python virtual environment created with conda. After installation, correct functionality of the OWCF cannot be guaranteed."
        user_reply = "0"
        while !(user_reply=="y" || user_reply=="n")
            global user_reply
            print("Do you want to continue the installation anyway [y/n]? ")
            user_reply = lowercase(readline())
        end
        if user_reply=="n"
            println("The OWCF installation will now be terminated. Please try to install the OWCF by providing a Python executable not within a conda virtual environment.")
            exit()
        end
    end
    user_reply = "0"
    while !(user_reply=="y" || user_reply=="n")
        global user_reply
        print("Are you attempting to install the OWCF with your virtual environment currently activated [y/n]? ")
        user_reply = lowercase(readline())
    end
    if user_reply=="n"
        error("To be able to install the OWCF, Julia needs to connect to Python. However, if you want to connect Julia to a Python virtual environment, you need to activate your virtual environment before you install the OWCF. Please remember, you will also need to activate your virtual environment every time before you use the OWCF.")
    end
end

user_reply = "0"
while !(user_reply=="y" || user_reply=="n")
    global user_reply
    print("Are you installing the OWCF on a computational cluster [y/n]? ")
    user_reply = lowercase(readline())
end
if user_reply=="y"
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
folderpath_OWCF = folderpath_OWCF*folder_delimiter # Add a "/" or a "\\" to the end of the folderpath_OWCF
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

print("- The path to your Python executable is... ")
python_script = """
import sys
PYTHON_EXECUTABLE_PATH = sys.executable 
with open('PYTHON_EXECUTABLE_PATH.txt','w') as f:
    f.write(PYTHON_EXECUTABLE_PATH)
""" # The Python script needed to determine the path to the Python executable
Base.run(`$PYTHON_COMMAND -c $python_script`) # Run the script with Python
f = open("PYTHON_EXECUTABLE_PATH.txt", "r") # Read the PYTHON_EXECUTABLE_PATH.txt file
PYTHON_EXECUTABLE_PATH = readline(f) # Read the first (and only) line of the .txt file
close(f); rm("$(folderpath_OWCF)PYTHON_EXECUTABLE_PATH.txt"; force=true) # Close and remove the .txt file
PYTHON_EXECUTABLE_PATH = replace(PYTHON_EXECUTABLE_PATH, ("Program Files" => "Program")) # Take care of Windows maybe using the folder "Program Files", which we cannot execute because of whitespace
# However, by default, Windows interprets commands like "C:\Program\..." as equivalent to "C:\Program Files\...". 
# So, we simply replace the "Program Files" part of the PYTHON_EXECUTABLE_PATH with "Program"
println("$PYTHON_EXECUTABLE_PATH")
###------------------------------------------------------------------------------------------------###

############---------------------------------------------------------------------------------------###
println("- Setting the PYTHON environment variable... ")
ENV["PYTHON"] = PYTHON_EXECUTABLE_PATH # Set the PYTHON environment variable
###------------------------------------------------------------------------------------------------###

############---------------------------------------------------------------------------------------###
println("- Instantiating and precompiling the OWCF (takes a looong time)... ")
Pkg.instantiate()
println("- The OWCF has been successfully instantiated and all Julia packages have been precompiled!")
println("- Please run 'julia $(folderpath_OWCF)tests$(folder_delimiter)run_tests.jl to ensure that the OWCF has been installed correctly.")
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