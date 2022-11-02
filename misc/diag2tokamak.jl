# This function is purely for esthetic purposes

function diag2tokamak(diagnostic::String)
    if lowercase(diagnostic)=="ab"
        return "JET"
    elseif lowercase(diagnostic)=="tofor"
        return "JET"
    elseif lowercase(diagnostic)=="mpru"
        return "JET"
    elseif lowercase(diagnostic)=="km15"
        return "JET"
    elseif split(lowercase(diagnostic),"-")[1]=="kn3"
        return "JET"
    else
        error("Specified diagnostic not recognized with any tokamak. Please change and re-try.")
    end
end