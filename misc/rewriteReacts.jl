########################################## rewriteReacts.jl #######################################
# This script contains functions that re-write reactions and species, to respect conventions of 
# different codes and/or file systems.
#
# pretty2scpok()
# This function converts a fusion reaction of the format a(b,c)d to 
# the format a-b=c-d. This is simply to allow for Linux and Windows file browsers (terminal/Powershell) 
# to be able to reconize file names containing fusion reactions without producing errors.
# For example, if I would try to securely copy the file
#
# orbWeights_JET_96100J01_TOFOR_D(d,n)3He_13x11x12.jld2
#
# across computers, I would get an error. Because, as soon as the Linux or Windows file browsers
# see the '(' or the ')' symbol, they think it's an escape sequence of its own. Therefore, to be 
# able to transfer files effectively across computers, the format a-b=c-d is used instead for
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
# Written by Henrik Järleblad. Last maintained 2025-03-05
##################################################################################################

include("availReacts.jl")

"""
    pretty2scpok(reaction)
    pretty2scpok(reaction; projVel=false)

Assume that reaction is written in the format a(b,c)d, or a(b,c)d-l, where a is thermal ion, b is 
fast ion, c is fusion product particle of interest, d is fusion product particle of disinterest and 
l is the nuclear energy state of c. l can be GS, 1L or 2L, corresponding to Ground State (GS), 
1st excited energy level (1L) and 2nd excited energy level (2L). Rewrite to A-b=c-D, or A-b=c-d+l, format.
'scpok' stands for 'secure copy protocol (is) ok'. Basically, that the file can be identified via a 
terminal/Powershell interface without any problems. If projVel, then there is no fusion reaction 
and we simply return 'proj-reaction'.
"""
function pretty2scpok(reaction::String; projVel::Bool=false)

    projVel && (return "proj-"*reaction) # Shortcut, if we know it's a projected velocity form (fast-ion species only)

    thermal_reactant, fast_reactant = getFusionReactants(reaction)
    product_of_interest, product_of_disinterest = getFusionProducts(reaction)
    reaction_scpok = thermal_reactant*"-"*fast_reactant*"="*product_of_interest*"-"*product_of_disinterest

    if getReactionForm(reaction)==1
        return reaction_scpok
    elseif getReactionForm(reaction)==2
        return reaction_scpok*"+"*getEmittedParticleEnergyLevel(reaction)
    else
        return "proj-"*reaction
    end
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