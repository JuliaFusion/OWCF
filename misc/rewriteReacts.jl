########################################## rewriteReacts.jl #######################################
# This simple script contains a function that converts a fusion reaction of the format A(b,c)D to 
# the format A-b=c-D.
#
# In future versions of the OWCF, it is envisioned how several re-writing functions might be added
# to this script, and forming a small library of re-writing functions. This would ensure cross-platform
# support for the OWCF and its results files.
#
# This is simply to allow for Linux and Windows file browsers (terminal/Powershell) to be 
# able to reconize file names containing fusion reactions without producing errors.
# For example, if I would try to securely copy the file
#
# orbWeights_JET_96100J01_TOFOR_D(d,n)3He_13x11x12.jld2
#
# across computers, I would get an error. Because, as soon as the Linux or Windows file browsers
# see the '(' or the ')' symbol, they think it's an escape sequence of its own. Therefore, to be 
# able to transfer files effectively across computers, the format A-B=c-D is used instead for
# fusion reactions. At least when they need to be included in the name of a file. So, the file
#
# orbWeights_JET_96100J01_TOFOR_D-d=n-3He_13x11x12.jld2
#
# would produce no errors and would be able to be copied perfectly across computers.
#
# Written by Henrik Järleblad. Last maintained 2022-02-11
##################################################################################################

"""
    pretty2scpok(reaction_full)
    pretty2scpok(reaction_full; projVel=false)

Assume that reaction is written in the format A(b,c)D.
Rewrite to A-b=c-D format. 'scpok' stands for 'secure copy protocol (is) ok'.
Basically, that the file can be identified via a terminal/Powershell interface
without any problems. If projVel, then there is no fusion reaction and we simply 
return the 'reaction_full' string instead (which should be on the 'proj-X' format).
"""
function pretty2scpok(reaction_full::String; projVel=false)

    projVel && (return reaction_full)

    reactants = split(reaction_full,",")[1]
    products = split(reaction_full,",")[2]

    bulk_reactant = split(reactants,"(")[1]
    fast_reactant = split(reactants,"(")[2]
    fast_product = split(products,")")[1]
    bulk_product = split(products,")")[2]

    return bulk_reactant*"-"*fast_reactant*"="*fast_product*"-"*bulk_product
end