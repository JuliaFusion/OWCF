################################ availReacts.jl ######################################################
# This script simply keeps track of the available fusion reactions that the OWCF can simulate.
#
# In addition, lists of reactions and particles with specific characteristics are included.
#
# PLEASE NOTE! A proton is denoted with 'h' and NOT 'p'. 

### Inputs
# fusion_reaction - Can be specified on the following forms:
#                   (1) "a(b,c)d"
#                   (2) "a(b,c)d-l"
#                   (3) "b"
#                 where a is the thermal ion, b is fast ion, c is emitted particle, d is the other 
#                 fusion product and l is the atomic nucleus energy level of the emitted particle c. 
#                 l can be GS (Ground State), 1L (1st excited energy level) and 2L (2nd excited energy level). 
#                 Form (1) should be used for all fusion reactions that do NOT result in an emitted particle 
#                 with an excited atomic nucleus level, e.g. D(T,n)4He. Form (2) should be used for all 
#                 fusion reactions that DO result in an emitted particle with an excited atomic nucleus 
#                 level, e.g. 9Be(4He,12C)n. Form (2) with l equal to GS is equivalent to form (1), i.e.
#                 "a(b,c)d-GS"=="a(b,c)d" evaluates to true. Form (3) is only used when projected velocites 
#                 are to be computed for the particle species b. PLEASE NOTE! Specify alpha particles as 
#                 '4he' or '4He' (NOT 'he4' or 'He4'). The same goes for helium-3 (specify as '3he', NOT 'he3') 
#                 etc - String

### List of particle species
# Please see the OWCF/misc/species_func.jl for a list of available particle species in the OWCF.
#
# DO NOT ALTER THIS SCRIPT. More fusion reactions might be added in the future.
#
# Script written by Henrik Järleblad. Last maintained 2025-03-03.
######################################################################################################

include("species_func.jl")

OWCF_AVAILABLE_FUSION_REACTIONS = ["D(D,n)3He-GS","T(h,g)4He-GS","h(T,g)4He-GS","D(T,n)4He-GS","T(D,n)4He-GS","D(3He,h)4He-GS","3He(D,h)4He-GS","9Be(4He,12C)n-1L","4He(9Be,12C)n-1L","9Be(4He,12C)n-2L","4He(9Be,12C)n-2L"]
OWCF_AVAILABLE_FUSION_REACTIONS_FOR_ANALYTIC_COMPUTATION = ["D(T,n)4He-GS","T(D,n)4He-GS","9Be(4He,12C)n-1L","4He(9Be,12C)n-1L","9Be(4He,12C)n-2L","4He(9Be,12C)n-2L"]

function getReactionForm(fusion_reaction::String)
    if occursin(fusion_reaction,"-")
        return 2
    end
    if occursin(fusion_reaction,"(")
        return 1
    end
    return 3
end

function getEmittedParticle(fusion_reaction::String)
    emitted_particle = fusion_reaction # Assume form (3) by default
    if occursin("(",emitted_particle) # If wrong... 
        fusion_reaction_form_1 = fusion_reaction # Assume form (1)
        if occursin("-",fusion_reaction_form_1) # If wrong...
            # Must be form (2)
            fusion_reaction_form_1, energy_level = Tuple(split(fusion_reaction,"-")) # Extract fusion reaction without energy level of the emitted particle
        end
        emitted_particle = split(split(fusion_reaction_form_1,",")[2],")")[1] # The emitted particle (")" as delimiter) obtained from the fusion products ("," as delimiter)
        return emitted_particle
    end
    return emitted_particle # If function reaches this return statement, input 'fusion_reaction' must be given on form (3) (please see file description above)
end

function getEmittedParticleEnergyLevel(fusion_reaction::String; verbose::Bool=false)
    if !occursin("-",fusion_reaction)
        verbose &&  println("No energy level specified for the $(getEmittedParticle(fusion_reaction)) particle of the $(fusion_reaction) reaction. Assuming GS (ground state).")
        return "GS"
    end
    return split(fusion_reaction,"-")[2]
end

function reactionIsAvailable(fusion_reaction::String; verbose::Bool=false)
    if getReactionForm(fusion_reaction)==3
        verbose && println("Fusion reaction is on form 3!")
        if fusion_reaction in OWCF_SPECIES
            return true
        end
    end

    if getReactionForm(fusion_reaction)==1
        verbose && println("Fusion reaction is on form 1!")
        verbose && println("Adding '-GS' to fusion reaction for check... ")
        fusion_reaction = fusion_reaction*"-GS"
    end

    if lowercase(fusion_reaction) in lowercase.(OWCF_AVAILABLE_FUSION_REACTIONS)
        return true
    else
        return false
    end
end

function reactionIsAvailableAnalytically(fusion_reaction::String; verbose::Bool=false)
    if getReactionForm(fusion_reaction)==3
        verbose && println("Expected diagnostic spectra analytically computed from projected velocity spectra is not possible.")
        return false
    end

    if getReactionForm(fusion_reaction)==1
        verbose && println("Fusion reaction a(d,c)d input has form 1. Assuming ground state for particle c and adding '-GS' for check... ")
        fusion_reaction *= "-GS"
    end

    if lowercase(fusion_reaction) in lowercase.(OWCF_AVAILABLE_FUSION_REACTIONS_FOR_ANALYTIC_COMPUTATION)
        return true
    else
        return false
    end
end

function getFusionReactants(fusion_reaction::String)
    if getReactionForm(fusion_reaction)==3
        return "proj", fusion_reaction
    end
    if getReactionForm(fusion_reaction)==2
        fusion_reaction = split(fusion_reaction,"-")[1]
    end
    reactants = split(split(fusion_reaction,",")[1],"(")
    thermal_particle = reactants[1]; fast_particle = reactants[2]
    return thermal_particle, fast_particle
end

function getFusionProducts(fusion_reaction::String)
    if getReactionForm(fusion_reaction)==3
        @warn "Only (fast) ion species specified in the fusion reaction input ($(fusion_reaction)). Returning 'proj','proj'."
        return "proj","proj"
    end
    if getReactionForm(fusion_reaction)==2
        fusion_reaction = split(fusion_reaction,"-")[1]
    end
    products = split(split(fusion_reaction,",")[1],")")
    product_of_interest = products[1]; product_of_disinterest = products[2]
    return product_of_interest, product_of_disinterest
end

function getThermalParticleSpecies(fusion_reaction::String)
    return getFusionReactants(fusion_reaction)[1]
end
getSlowParticleSpecies = getThermalParticleSpecies

function getFastParticleSpecies(fusion_reaction::String)
    return getFusionReactants(fusion_reaction)[2]
end
getEnergeticParticleSpecies = getFastParticleSpecies

function full2reactsOnly(fusion_reaction::String)
    thermal_species, FI_species = getFusionReactants(fusion_reaction)

    # Convert the full fusion reaction to reactant format (used as input to DRESS Python code)
    return thermal_species*"-"*FI_species
end
