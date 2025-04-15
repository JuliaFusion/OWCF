# This simple tutorial will show you how to manually create a wall for the extra/createCustomMagneticEquilibrium.jl script,
# if you don't want the algorithm to automatically decide on a wall for you.
# It is recommended that you execute this script block-by-block in VSCode.
# You do this by highlighting the code of interest, and press ctrl + enter (you can also use shift + enter, to advance line-by-line)
# Each block is separated by a line such as the one below.
# The folderpath_OWCF string you have to fill in yourself. Remember to finish with "/"
##############################################################################################################
folderpath_OWCF = "/home/henrikj/Codes/OWCF/" # First, you have to specify the path to the OWCF folder, as a string. Change this to your own OWCF path
cd(folderpath_OWCF) # Then, you make the VSCode terminal change place to the OWCF folder
using Pkg # Then, you load the Pkg.jl package which is Julia's package manager
Pkg.activate(".") # Then, you activate the virtual environment of the OWCF
##############################################################################################################
using JLD2 # Then, we load the JLD2 package, to be able to save files in the Julia format
using Plots # Then, we load the Plots.jl package, to be able to plot things
##############################################################################################################
# To create (R,z) point data, create two arrays manually and fill them with values.
# That is, replace (R1,R2,...) and (z1,z2,...) with actual (R,z) values
R = [R1,R2,...]
z = [z1,z2,...]
##############################################################################################################
# We can then save our created data via the functions provided by the JLD2.jl package
# Save as a .jld2 file
# The syntax below is the most safe way to save a .jld2 file. First, you specify the filepath (including the .jld2 file extension)
# Then, you specify if the file should be created or not (the first 'true' in the line below)
# Then, you specify if the file should be written to or not (the second 'true' in the line below)
# Then, you specify if the file should be truncated or not (the 'false' in the line below)
# Finally, you specify the output stream type. 'IOStream' is the most safe way of saving/writing. Therefore, we use it
myfile = jldopen("/path/and/filename/to/where/you/want/to/save/the/wall/data.jld2",true,true,false,IOStream) # The jldopen() function opens a Julia .jld2 file
write(myfile,"R",R) # Write the variable R to the "R" key in the .jld2 file
write(myfile,"z",z) # Write the variable z to the "z" key in the .jld2 file
close(myfile) # Don't forget to close the file 

# That concludes this very short tutorial.