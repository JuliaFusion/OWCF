########################################## convert_units.jl #######################################
# This script contains a record of all accepted units of measurement that can be handled by the 
# OWCF. It also contains functions that computes conversion factors between different units of 
# measurement.
#
# To make inverse units, please specify them with a "^-p" where 'p' is a specific power. For 
# example, "m^-2" is 'per meter squared'.
# 
# To make compository units, please separate them with a "_". For example, "m_s^-1" is 'meters per 
# second'.
#
# Written by Henrik Järleblad. Last maintained 2024-10-17
###################################################################################################

LENGTH_UNITS = Dict("mm"=>1.0e-3,"cm"=>1.0e-2,"dm"=>1.0e-1,"m"=>1.0e0,"km"=>1.0e3)
LENGTH_UNITS_LONG = Dict("millimeter"=>1.0e-3,"centimeter"=>1.0e-2,"decimeter"=>1.0e-1,"meter"=>1.0e0,"kilometer"=>1.0e3)

TIME_UNITS = Dict("ns"=>1.0e-9,"microsec"=>1.0e-6,"ms"=>1.0e-3,"s"=>1.0e0,"min"=>60e0,"h"=>3.6e3)
TIME_UNITS_LONG = Dict("nanosecond"=>1.0e-9,"microsecond"=>1.0e-6,"millisecond"=>1.0e-3,"second"=>1.0e0,"minute"=>60e0,"hour"=>3.6e3)

ENERGY_UNITS = Dict("ev"=>1.0e0,"kev"=>1.0e3,"mev"=>1.0e6)
ENERGY_UNITS_LONG = Dict("electronvolt"=>1.0e0,"kiloelectronvolt"=>1.0e3,"megaelectronvolt"=>1.0e6)

DIMENSIONLESS_UNITS = Dict("a.u."=>1,"au"=>1,"-"=>1,"signal"=>1,"counts"=>1)
DIMENSIONLESS_UNITS_LONG = Dict("arbitraryunits"=>1,"unknown"=>1,"signal"=>1,"counts"=>1)

OWCF_UNITS = merge(LENGTH_UNITS, TIME_UNITS, ENERGY_UNITS, DIMENSIONLESS_UNITS)
OWCF_UNITS_LONG = merge(LENGTH_UNITS_LONG, TIME_UNITS_LONG, ENERGY_UNITS_LONG, DIMENSIONLESS_UNITS_LONG)

"""
    unit_conversion_factor(unit_in, unit_out)

Take the unit_in and unit_out input variables, and compute the resulting 
unit conversion factor between the two. For example, if unit_in="m" and unit_out="km", then 
the function output would be 1.0e-3.

If keyword argument rounding is set to true, the output will be rounded to 12 significant digits.
If keyword argument verbose is set to true, the function will talk a lot!
"""
function unit_conversion_factor(unit_in_orig::String,unit_out_orig::String; rounding=true, verbose=false)

    # Convert to lowercase
    unit_in = lowercase(unit_in_orig)
    unit_out = lowercase(unit_out_orig)

    # Decompose units into components. E.g. "m_s^-1" into ["m","s^-1"]
    unit_in_components = split(unit_in,"_")
    unit_out_components = split(unit_out,"_")

    # Process all units and powers in unit_in
    unit_in_components_mag = Dict("LENGTH"=>1.0, "TIME"=>1.0, "ENERGY"=>1.0, "DIMENSIONLESS"=>1.0) # Magnitudes of unit_in components. Since conversion factor is computed via ratios, initialize 1.0's
    unit_in_components_dim = Dict("LENGTH"=>0.0, "TIME"=>0.0, "ENERGY"=>0.0, "DIMENSIONLESS"=>0.0) # Dimensions of unit_in components. Must match dimensions of unit_out components (see below)
    for unit_in_component in unit_in_components
        if occursin("^",unit_in_component) # If unit has an exponent...
            unit, power = split(unit_in_component,"^") # e.g. "s^-1" into "s","-1"
            power = parse(Float64, power)
        else
            unit, power = unit_in_component, 1.0 # Otherwise, assume exponent of 1.0
        end
        if unit in keys(LENGTH_UNITS)
            if unit_in_components_mag["LENGTH"]==0.0 # If no unit is yet registered...
                unit_in_components_mag["LENGTH"] = 1.0 # ... set the count to 1.0 so that...
            end
            unit_in_components_mag["LENGTH"] *= LENGTH_UNITS[unit]^power # ... one can factor the unit magnitudes correctly
            unit_in_components_dim["LENGTH"] += power # Add the dimension
        elseif unit in keys(TIME_UNITS)
            if unit_in_components_mag["TIME"]==0.0 # If no unit magnitude is yet registered...
                unit_in_components_mag["TIME"] = 1.0 # ... set the count to 1.0 so that...
            end
            unit_in_components_mag["TIME"] *= TIME_UNITS[unit]^power # ... one can factor the unit magnitudes correctly
            unit_in_components_dim["TIME"] += power # Add the dimension
        elseif unit in keys(ENERGY_UNITS)
            if unit_in_components_mag["ENERGY"]==0.0 # If no unit magnitude is yet registered...
                unit_in_components_mag["ENERGY"] = 1.0 # ... set the count to 1.0 so that...
            end
            unit_in_components_mag["ENERGY"] *= ENERGY_UNITS[unit]^power # ... one can factor the unit magnitudes correctly
            unit_in_components_dim["ENERGY"] += power # Add the dimension
        elseif unit in keys(DIMENSIONLESS_UNITS)
            if unit_in_components_mag["DIMENSIONLESS"]==0.0 # If no unit magnitude is yet registered...
                unit_in_components_mag["DIMENSIONLESS"] = 1.0 # ... set the count to 1.0 so that...
            end
            unit_in_components_mag["DIMENSIONLESS"] *= DIMENSIONLESS_UNITS[unit]^power # ... one can factor the unit magnitudes correctly
            unit_in_components_dim["DIMENSIONLESS"] += power # Add the dimension
        else
            verbose && println("unit_in component '"*unit*"' not found in accepted OWCF LENGTH, TIME, ENERGY and DIMENSIONLESS units. Trying LONG format... ")
            if unit in keys(LENGTH_UNITS_LONG)
                if unit_in_components_mag["LENGTH"]==0.0 # If no unit is yet registered...
                    unit_in_components_mag["LENGTH"] = 1.0 # ... set the count to 1.0 so that...
                end
                unit_in_components_mag["LENGTH"] *= LENGTH_UNITS_LONG[unit]^power # ... one can factor the unit magnitudes correctly
                unit_in_components_dim["LENGTH"] += power # Add the dimension
            elseif unit in keys(TIME_UNITS_LONG)
                if unit_in_components_mag["TIME"]==0.0 # If no unit is yet registered...
                    unit_in_components_mag["TIME"] = 1.0 # ... set the count to 1.0 so that...
                end
                unit_in_components_mag["TIME"] *= TIME_UNITS_LONG[unit]^power # ... one can factor the unit magnitudes correctly
                unit_in_components_dim["TIME"] += power # Add the dimension
            elseif unit in keys(ENERGY_UNITS_LONG)
                if unit_in_components_mag["ENERGY"]==0.0 # If no unit is yet registered...
                    unit_in_components_mag["ENERGY"] = 1.0 # ... set the count to 1.0 so that...
                end
                unit_in_components_mag["ENERGY"] *= ENERGY_UNITS_LONG[unit]^power # ... one can factor the unit magnitudes correctly
                unit_in_components_dim["ENERGY"] += power # Add the dimension
            elseif unit in keys(DIMENSIONLESS_UNITS_LONG)
                if unit_in_components_mag["DIMENSIONLESS"]==0.0 # If no unit is yet registered...
                    unit_in_components_mag["DIMENSIONLESS"] = 1.0 # ... set the count to 1.0 so that...
                end
                unit_in_components_mag["DIMENSIONLESS"] *= DIMENSIONLESS_UNITS_LONG[unit]^power # ... one can factor the unit magnitudes correctly
                unit_in_components_dim["DIMENSIONLESS"] += power # Add the dimension
            else
                error("unit_in componenent '"*unit*"' (all lowercase) not found in accepted OWCF units. Please see OWCF/misc/convert_units.jl for lists of accepted units. Please correct and re-try.")
            end
        end
    end

    # Process all units and powers in unit_out
    unit_out_components_mag = Dict("LENGTH"=>1.0, "TIME"=>1.0, "ENERGY"=>1.0, "DIMENSIONLESS"=>1.0) # Magnitudes of unit_in components. Since conversion factor is computed via ratios, initialize 1.0's
    unit_out_components_dim = Dict("LENGTH"=>0.0, "TIME"=>0.0, "ENERGY"=>0.0, "DIMENSIONLESS"=>0.0) # Dimensions of unit_in components. Must match dimensions of unit_out components (see below)
    for unit_out_component in unit_out_components
        if occursin("^",unit_out_component) # If unit has an exponent...
            unit, power = split(unit_out_component,"^") # e.g. "s^-1" into "s","-1"
            power = parse(Float64, power)
        else
            unit, power = unit_out_component, 1.0 # Otherwise, assume exponent of 1.0
        end
        if unit in keys(LENGTH_UNITS)
            if unit_out_components_mag["LENGTH"]==0.0 # If no unit is yet registered...
                unit_out_components_mag["LENGTH"] = 1.0 # ... set the count to 1.0 so that...
            end
            unit_out_components_mag["LENGTH"] *= LENGTH_UNITS[unit]^power # ... one can factor the unit magnitudes correctly
            unit_out_components_dim["LENGTH"] += power # Add the dimension
        elseif unit in keys(TIME_UNITS)
            if unit_out_components_mag["TIME"]==0.0 # If no unit magnitude is yet registered...
                unit_out_components_mag["TIME"] = 1.0 # ... set the count to 1.0 so that...
            end
            unit_out_components_mag["TIME"] *= TIME_UNITS[unit]^power # ... one can factor the unit magnitudes correctly
            unit_out_components_dim["TIME"] += power # Add the dimension
        elseif unit in keys(ENERGY_UNITS)
            if unit_out_components_mag["ENERGY"]==0.0 # If no unit magnitude is yet registered...
                unit_out_components_mag["ENERGY"] = 1.0 # ... set the count to 1.0 so that...
            end
            unit_out_components_mag["ENERGY"] *= ENERGY_UNITS[unit]^power # ... one can factor the unit magnitudes correctly
            unit_out_components_dim["ENERGY"] += power # Add the dimension
        elseif unit in keys(DIMENSIONLESS_UNITS)
            if unit_out_components_mag["DIMENSIONLESS"]==0.0 # If no unit magnitude is yet registered...
                unit_out_components_mag["DIMENSIONLESS"] = 1.0 # ... set the count to 1.0 so that...
            end
            unit_out_components_mag["DIMENSIONLESS"] *= DIMENSIONLESS_UNITS[unit]^power # ... one can factor the unit magnitudes correctly
            unit_out_components_dim["DIMENSIONLESS"] += power # Add the dimension
        else
            verbose && println("unit_out component '"*unit*"' not found in accepted OWCF LENGTH, TIME, ENERGY and DIMENSIONLESS units. Trying LONG format... ")
            if unit in keys(LENGTH_UNITS_LONG)
                if unit_out_components_mag["LENGTH"]==0.0 # If no unit is yet registered...
                    unit_out_components_mag["LENGTH"] = 1.0 # ... set the count to 1.0 so that...
                end
                unit_out_components_mag["LENGTH"] *= LENGTH_UNITS_LONG[unit]^power # ... one can factor the unit magnitudes correctly
                unit_out_components_dim["LENGTH"] += power # Add the dimension
            elseif unit in keys(TIME_UNITS_LONG)
                if unit_out_components_mag["TIME"]==0.0 # If no unit is yet registered...
                    unit_out_components_mag["TIME"] = 1.0 # ... set the count to 1.0 so that...
                end
                unit_out_components_mag["TIME"] *= TIME_UNITS_LONG[unit]^power # ... one can factor the unit magnitudes correctly
                unit_out_components_dim["TIME"] += power # Add the dimension
            elseif unit in keys(ENERGY_UNITS_LONG)
                if unit_out_components_mag["ENERGY"]==0.0 # If no unit is yet registered...
                    unit_out_components_mag["ENERGY"] = 1.0 # ... set the count to 1.0 so that...
                end
                unit_out_components_mag["ENERGY"] *= ENERGY_UNITS_LONG[unit]^power # ... one can factor the unit magnitudes correctly
                unit_out_components_dim["ENERGY"] += power # Add the dimension
            elseif unit in keys(DIMENSIONLESS_UNITS_LONG)
                if unit_out_components_mag["DIMENSIONLESS"]==0.0 # If no unit is yet registered...
                    unit_out_components_mag["DIMENSIONLESS"] = 1.0 # ... set the count to 1.0 so that...
                end
                unit_out_components_mag["DIMENSIONLESS"] *= DIMENSIONLESS_UNITS_LONG[unit]^power # ... one can factor the unit magnitudes correctly
                unit_out_components_dim["DIMENSIONLESS"] += power # Add the dimension
            else
                error("unit_out componenent '"*unit*"' (all lowercase) not found in accepted OWCF units. Please see OWCF/misc/convert_units.jl for lists of accepted units. Please correct and re-try.")
            end
        end
    end

    # CHECK THAT DIMENSIONS MATCH
    if !(unit_in_components_dim==unit_out_components_dim)
        error("unit_in ("*unit_in_orig*") and unit_out ("*unit_out_orig*") do not have the same dimensions. Please correct and re-try.")
    end

    # COMPUTE CONVERSION FACTOR
    cf = 1
    for key in keys(unit_in_components_mag)
        cf *= unit_in_components_mag[key]/unit_out_components_mag[key]
    end

    # RETURN CONVERSION FACTOR
    return rounding ? round(cf,sigdigits=12) : cf
end