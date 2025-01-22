################################ availReacts.jl ######################################################
# This script simply keeps track of the available fusion reactions that the OWCF can simulate via the
# interface with the DRESS code.
#
### Inputs
# reaction_full - Specified on the form a(b,c)d where a is thermal ion, b is fast ion,
#            c is emitted particle and d is the product nucleus. PLEASE NOTE! Specify alpha
#            particles as '4he' or '4He' (NOT 'he4' or 'He4'). Same goes for helium-3
#            (specify as '3he', NOT 'he3') - String
#
### List of particle species
# d (or D) - Deuterium
# t (or T) - Tritium
# p - proton
# 3he (or 3He) - Helium-3
# 4he (or 4He) - Alpha particle
# g - Gamma photon
# n - Neutron
#
# DO NOT ALTER THIS SCRIPT. More fusion reactions might be added in the future.
#
# Script written by Henrik Järleblad. Last maintained 2022-08-26.
######################################################################################################

function getEmittedParticle(reaction_full::String)
    products = split((split(reaction_full,","))[2],")")

    return products[1] # Return 'c' in the fusion reaction a(b,c)d
end

function checkReaction(reaction_full::String; verbose::Bool=false, projVelocity::Bool=false)
    if projVelocity
        thermal_species = split(reaction_full,"-")[1]
        FI_species = split(reaction_full,"-")[2]
        verbose && println("Accepted velocity projection: "*reaction_full)
        return thermal_species, FI_species
    end
    
    reactants = split((split(reaction_full,","))[1],"(")
    products = split((split(reaction_full,","))[2],")")
    if !(length(reactants)==2 && length(products)==2)
        error("Fusion reaction not specified correctly. Please use the format a(b,c)d and re-try.")
    end

    thermal_reactant = reactants[1]
    fast_reactant = reactants[2]
    emitted_particle = products[1]
    product_nucleus = products[2]

    # Account for the possibility of faulty user input. D(D,n)3He is correct. But user might input D(D,3He)n.
    if lowercase(thermal_reactant)=="d" && lowercase(fast_reactant)=="d" && ((lowercase(emitted_particle)=="n" && lowercase(product_nucleus)=="3he") || (lowercase(emitted_particle)=="3he" && lowercase(product_nucleus)=="n"))
        verbose && println("Accepted nuclear reaction: "*reaction_full)
        return "d","d","DD_neutrons"
    end

    # Account for the possibility of faulty user input. T(p,g)4He is correct. But user might input T(p,4He)g.
    if lowercase(thermal_reactant)=="t" && lowercase(fast_reactant)=="p" && ((lowercase(emitted_particle)=="g" && lowercase(product_nucleus)=="4he") || (lowercase(emitted_particle)=="4he" && lowercase(product_nucleus)=="g"))
        verbose && println("Accepted nuclear reaction: "*reaction_full)
        return "T","p","Tp_gammas"
    end

    # Account for the possibility of faulty user input. p(T,g)4He is correct. But user might input p(T,4He)g.
    if lowercase(thermal_reactant)=="p" && lowercase(fast_reactant)=="t" && ((lowercase(emitted_particle)=="g" && lowercase(product_nucleus)=="4he") || (lowercase(emitted_particle)=="4he" && lowercase(product_nucleus)=="g"))
        verbose && println("Accepted nuclear reaction: "*reaction_full)
        return "p","T","pT_gammas"
    end

    # Account for the possibility of faulty user input. D(T,n)4He is correct. But user might input D(T,4He)n.
    if lowercase(thermal_reactant)=="d" && lowercase(fast_reactant)=="t"
        if lowercase(emitted_particle)=="n" && lowercase(product_nucleus)=="4he"
            verbose && println("Accepted nuclear reaction: "*reaction_full)
            return "D","T","DT_neutrons"
        elseif lowercase(emitted_particle)=="4he" && lowercase(product_nucleus)=="n"
            verbose && println("Accepted nuclear reaction: "*reaction_full)
            return "D","T","DT_alphas"
        end
    end

    # Account for the possibility of faulty user input. T(D,n)4He is correct. But user might input T(D,4He)n.
    if lowercase(thermal_reactant)=="t" && lowercase(fast_reactant)=="d"
        if lowercase(emitted_particle)=="n" && lowercase(product_nucleus)=="4he"
            verbose && println("Accepted nuclear reaction: "*reaction_full)
            return "T","D","TD_neutrons"
        elseif lowercase(emitted_particle)=="4he" && lowercase(product_nucleus)=="n"
            verbose && println("Accepted nuclear reaction: "*reaction_full)
            return "T","D","TD_neutrons"
        end
    end

    # Account for the possibility of faulty user input. 3He(D,p)4He is correct. But user might input 3He(D,4He)p.
    if lowercase(thermal_reactant)=="3he" && lowercase(fast_reactant)=="d" && ((lowercase(emitted_particle)=="p" && lowercase(product_nucleus)=="4he") || (lowercase(emitted_particle)=="4he" && lowercase(product_nucleus)=="p"))
        verbose && println("Accepted nuclear reaction: "*reaction_full)
        return "3He","D","3HeD_protons"
    end

    # Account for the possibility of faulty user input. D(3He,p)4He is correct. But user might input D(3He,4He)p.
    if lowercase(thermal_reactant)=="D" && lowercase(fast_reactant)=="3he" && ((lowercase(emitted_particle)=="p" && lowercase(product_nucleus)=="4he") || (lowercase(emitted_particle)=="4he" && lowercase(product_nucleus)=="p"))
        verbose && println("Accepted nuclear reaction: "*reaction_full)
        return "D","3He","D3He_protons"
    end

    #!# Account for the possibility of faulty user input. 9Be(4He,12C)n is correct. 
    if lowercase(thermal_reactant)=="9be" && lowercase(fast_reactant)=="4he" && ((lowercase(emitted_particle)=="12c" && lowercase(product_nucleus)=="n") || (lowercase(emitted_particle)=="n" && lowercase(product_nucleus)=="12c"))
        verbose && println("Accepted nuclear reaction: "*reaction_full)
        return "9Be","4He","9Be4He_gammas"
    end    

    # No valid reaction detected
    error("Reaction "*reaction_full*" is not a valid fusion reaction available in the current version of the OWCF. Please re-specify and re-try.")

    return "X","X" # Will never be reached
end

function full2reactsOnly(reaction_full::String; kwargs...)
    thermal_species, FI_species, reactionName = checkReaction(reaction_full; kwargs...)

    # Convert the full fusion reaction to reactant format (used as input to Python framework)
    return thermal_species*"-"*FI_species
end

function getReactionName(reaction_full::String; kwargs...)
    thermal_species, FI_species, reaction_name = checkReaction(reaction_full; kwargs...)

    # Return the reaction type string only
    return reaction_name
end