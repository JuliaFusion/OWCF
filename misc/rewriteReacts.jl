########################################## rewriteReacts.jl #######################################
# This script contains functions that re-write reactions and species, to respect conventions of 
# different codes and/or file systems.
#
# pretty2scpok()
# This function converts a fusion reaction of the format A(b,c)D to 
# the format A-b=c-D. This is simply to allow for Linux and Windows file browsers (terminal/Powershell) 
# to be able to reconize file names containing fusion reactions without producing errors.
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
#
# thermalSpecies_OWCFtoTRANSP()
# This function converts a particle species notation from OWCF to TRANSP standard.
# For example, the OWCF would write helium-3 as 3he, while TRANSP would write it as HE3.
# This conversion function is simply to be able to let the OWCF seemlessly interact with TRANSP.
#
# Written by Henrik Järleblad. Last maintained 2023-08-25
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
function pretty2scpok(reaction_full::String; projVel::Bool=false)

    projVel && (return reaction_full)

    reactants = split(reaction_full,",")[1]
    products = split(reaction_full,",")[2]

    bulk_reactant = split(reactants,"(")[1]
    fast_reactant = split(reactants,"(")[2]
    fast_product = split(products,")")[1]
    bulk_product = split(products,")")[2]

    return bulk_reactant*"-"*fast_reactant*"="*fast_product*"-"*bulk_product
end

"""
    thermalSpecies_OWCFtoTRANSP(species)

Re-write a species notation from OWCF to TRANSP standard.

Example:

julia> x = thermalSpecies_OWCFtoTRANSP("3he")
julia> println(x)
"HE3"
"""
function thermalSpecies_OWCFtoTRANSP(species::String)
    if lowercase(species)=="d"
        return "D"
    elseif lowercase(species)=="t"
        return "T"
    elseif lowercase(species)=="3he"
        return "HE3"
    elseif lowercase(species)=="h"
        return "H"
    elseif lowercase(species)=="e"
        return "E"
    elseif lowercase(species)=="4he"
        return "HE4"
    elseif lowercase(species)=="alpha"
        return "HE4"
    else
        @warn "Did not recognize particle species. Output equal to input."
        return species # If unknown, return input itself but give a warning
    end
end