########################################## convert_units.jl #######################################
# This script contains a record of all accepted units of measurement that can be handled by the 
# OWCF. It also contains functions that computes conversion factors between different units of 
# measurement. The base units are defined using the SI base units. Please see, for example,
# https://en.wikipedia.org/wiki/SI_base_unit.
#
# PLEASE NOTE! Conversion for imperial units are currently NOT supported. Only metric units.
#
# To make inverse units, please specify them with a "^-p" where 'p' is a specific power. For 
# example, "m^-2" is 'per meter squared'.
# 
# To make compository units, please separate them with a "_". For example, "m_s^-1" is 'meters per 
# second'.
#
### Other
# For angles, please use radians (rad) and NOT degrees (deg).
#
# Written by Henrik Järleblad. Last maintained 2024-10-31
###################################################################################################

TIME_UNITS = Dict("ns"=>1.0e-9,"microsec"=>1.0e-6,"ms"=>1.0e-3,"s"=>1.0e0,"min"=>60e0,"h"=>3.6e3)
TIME_UNITS_LONG = Dict("nanosecond"=>1.0e-9,"microsecond"=>1.0e-6,"millisecond"=>1.0e-3,"second"=>1.0e0,"minute"=>60e0,"hour"=>3.6e3)

LENGTH_UNITS = Dict("mm"=>1.0e-3,"cm"=>1.0e-2,"dm"=>1.0e-1,"m"=>1.0e0,"km"=>1.0e3)
LENGTH_UNITS_LONG = Dict("millimeter"=>1.0e-3,"centimeter"=>1.0e-2,"decimeter"=>1.0e-1,"meter"=>1.0e0,"kilometer"=>1.0e3)

MASS_UNITS = Dict("amu"=>1.6605402e-27,"Da"=>1.66e-27,"mg"=>1.0e-6,"g"=>1.0-3,"hg"=>1.0e-1,"kg"=>1.0e0,"t"=>1.0e3)
MASS_UNITS_LONG = Dict("atomicmassunit"=>1.6605402e-27,"dalton"=>1.66e-27,"milligram"=>1.0e-6,"gram"=>1.0-3,"hectogram"=>1.0e-1,"kilogram"=>1.0e0,"tonne"=>1.0e3)

CURRENT_UNITS = Dict("mA"=>1.0e-3,"A"=>1.0,"kA"=>1.0e3,"MA"=>1.0e6)
CURRENT_UNITS_LONG = Dict("milliampere"=>1.0e-3,"ampere"=>1.0,"kiloampere"=>1.0e3,"megaampere"=>1.0e6)

TEMPERATURE_UNITS = Dict("K"=>1.0,"dC"=>-272.15,"dF"=>-457.87)
TEMPERATURE_UNITS_LONG = Dict("kelvin"=>1.0, "degreescelcius"=>-272.15,"degreesfahrenheit"=>-457.87)

SUBSTANCE_UNITS = Dict("mol"=>1.0,"NA"=>6.02214076e23)
SUBSTANCE_UNITS_LONG = Dict("mole"=>1.0,"avogradosconstant"=>6.02214076e23)

LUMINOUSITY_UNITS = Dict("HK"=>0.92,"cp"=>0.981,"cd"=>1.0)
LUMINOUSITY_UNITS_LONG = Dict("hefnerkerze"=>0.92,"candlepower"=>0.981,"candela"=>1.0)

DIMENSIONLESS_UNITS = Dict("a.u."=>1.0,"au"=>1.0,"-"=>1.0,"signal"=>1.0,"counts"=>1.0,"count"=>1.0,"rad"=>1.0,"sr"=>1.0)
DIMENSIONLESS_UNITS_LONG = Dict("arbitraryunits"=>1.0,"unknown"=>1.0,"dimensionless"=>1.0,"signal"=>1.0,"counts"=>1.0,"count"=>1.0,"radian"=>1.0,"steradian"=>1.0)

OWCF_BASE_UNITS = merge(TIME_UNITS,LENGTH_UNITS,MASS_UNITS,CURRENT_UNITS,TEMPERATURE_UNITS,SUBSTANCE_UNITS,LUMINOUSITY_UNITS,DIMENSIONLESS_UNITS)
OWCF_BASE_UNITS_LONG = merge(TIME_UNITS_LONG,LENGTH_UNITS_LONG,MASS_UNITS_LONG,CURRENT_UNITS_LONG,TEMPERATURE_UNITS_LONG,SUBSTANCE_UNITS_LONG,LUMINOUSITY_UNITS_LONG,DIMENSIONLESS_UNITS_LONG)

# Converts to m_s^-1
SPEED_UNITS = Dict("kn"=>0.514444,"c"=>299792458.0)
SPEED_UNITS_LONG = Dict("knot"=>0.514444,"speedoflight"=>299792458.0)

# Converts to m_s^-2
ACCELERATION_UNITS = Dict("Gal"=>0.01,"g0"=>9.80665)
ACCELERATION_UNITS_LONG = Dict("galileo"=>0.01,"standardgravity"=>9.80665)

# Converts to kg_m_s^-2
FORCE_UNITS = Dict("N"=>1.0,"dyn"=>1.0e-5,"kp"=>9.80665,"lbf"=>4.448222,"pdl"=>0.138255) 
FORCE_UNITS_LONG = Dict("newton"=>1.0,"dyne"=>1.0e-5,"kilopond"=>9.80665,"poundforce"=>4.448222,"poundal"=>0.138255) 

# Converts to kg_m^-1_s^-2
PRESSURE_UNITS = Dict("Pa"=>1.0,"hPa"=>1.0e2,"mmHg"=>133.322387415,"kPa"=>1.0e3,"bar"=>1.0e5,"atm"=>101325.0)
PRESSURE_UNITS_LONG = Dict("pascal"=>1.0,"hectopascal"=>1.0e2,"millimetermercury"=>133.322387415,"kilopascal"=>1.0e3,"bar"=>1.0e5,"atmosphere"=>101325.0)

# Converts to kg_m^2_s^-2
e0_OWCF = 1.60217733e-19 # Coulombs / Joules
ENERGY_UNITS = Dict("eV"=>e0_OWCF,"keV"=>1.0e3*e0_OWCF,"MeV"=>1.0e6*e0_OWCF,"GeV"=>1.0e9*e0_OWCF,"mJ"=>1.0e-3,"J"=>1.0e0,"cal"=>4.184,"kJ"=>1.0e3,"MJ"=>1.0e6)
ENERGY_UNITS_LONG = Dict("electronvolt"=>e0_OWCF,"kiloelectronvolt"=>1.0e3*e0_OWCF,"megaelectronvolt"=>1.0e6*e0_OWCF,"gigaelectronvolt"=>1.0e9*e0_OWCF,"millijoule"=>1.0e-3,"joule"=>1.0e0,"calorie"=>4.184,"kilojoule"=>1.0e3,"megajoule"=>1.0e6)

# Converts to kg_m^2_s^-3
POWER_UNITS = Dict("W"=>1.0,"hp"=>745.7,"kW"=>1.0e3,"RT"=>3516.853,"MW"=>1.0e6,"GW"=>1.0e9)
POWER_UNITS_LONG = Dict("watt"=>1.0,"horsepower"=>745.7,"kilowatt"=>1.0e3,"refrigerationton"=>3516.853,"megawatt"=>1.0e6,"gigawatt"=>1.0e9)

# Converts to A_s
CHARGE_UNITS = Dict("e"=>1.602176634e-19,"C"=>1.0,"F"=>9.648533212e4)
CHARGE_UNITS_LONG = Dict("elementarycharge"=>1.602176634e-19,"coulomb"=>1.0,"faraday"=>9.648533212e4)

# Converts to kg_m^2_s^-3_A^-1
VOLTAGE_UNITS = Dict("V"=>1.0,"kV"=>1.0e3)
VOLTAGE_UNITS_LONG = Dict("volt"=>1.0,"kilovolt"=>1.0e3)

OWCF_UNITS = merge(OWCF_BASE_UNITS, ACCELERATION_UNITS, FORCE_UNITS, PRESSURE_UNITS, ENERGY_UNITS, POWER_UNITS, CHARGE_UNITS, VOLTAGE_UNITS)
OWCF_UNITS_LONG = merge(OWCF_BASE_UNITS_LONG, ACCELERATION_UNITS_LONG, FORCE_UNITS_LONG, PRESSURE_UNITS_LONG, ENERGY_UNITS_LONG, POWER_UNITS_LONG, CHARGE_UNITS_LONG, VOLTAGE_UNITS_LONG)


"""
units_to_base_units(units_in)

Convert units of measurement into base units of measurement, including 
a conversion factor. Please see, for example, https://en.wikipedia.org/wiki/SI_base_unit.

As an example, units_in="keV_s^-1" will result in an output of
"kg_m^2_s^-3", 1.60217733e-16
because 1.60217733e-16 is the conversion factor from eV to J=kg_m^2_s^-2.

If keyword argument verbose is set to true, the function will talk a lot!
"""
function units_to_base_units(units_in::String; verbose=false)
    units_in_components = split(units_in,"_")
    units_out = ""
    cf_out = 1.0
    for units_in_component in units_in_components
        verbose && println("Converting $(units_in_component) to base units... ")
        if occursin("^",units_in_component) # If the component has an exponent...
            unit, power = split(units_in_component,"^") # e.g. "s^-1" into "s","-1"
            power = parse(Float64, power)
        else
            unit, power = units_in_component, 1.0 # Otherwise, assume exponent of 1.0
        end
        if unit in keys(OWCF_BASE_UNITS) || unit in keys(OWCF_BASE_UNITS_LONG)# If the unit is already a base unit ...
            verbose && println("Unit $(unit) found in list of OWCF base units! Continuing... ")
            units_out *= units_in_component*"_" # Add it to the units_out...
            continue # And continue to the next units_in_component (don't perform all the checks below)
        end

        if unit in keys(SPEED_UNITS)
            units_out *= power != 1 ? "m^$(Int64(power))_s^$(Int64(-power))_" : "m_s^-1_"
            cf_out *= SPEED_UNITS[unit]^power
        elseif unit in keys(ACCELERATION_UNITS)
            units_out *= power != 1 ? "m^$(Int64(power))_s^$(Int64(-2*power))_" : "m_s^-2_"
            cf_out *= ACCELERATION_UNITS[unit]^power
        elseif unit in keys(FORCE_UNITS)
            units_out *= power != 1 ? "kg^$(Int64(power))_m^$(Int64(power))_s^$(Int64(-2*power))_" : "kg_m_s^-2_"
            cf_out *= FORCE_UNITS[unit]^power
        elseif unit in keys(PRESSURE_UNITS)
            units_out *= power != 1 ? "kg^$(Int64(power))_m^$(Int64(-power))_s^$(Int64(-2*power))_" : "kg_m^-1_s^-2_"
            cf_out *= PRESSURE_UNITS[unit]^power
        elseif unit in keys(ENERGY_UNITS)
            units_out *= power != 1 ? "kg^$(Int64(power))_m^$(Int64(2*power))_s^$(Int64(-2*power))_" : "kg_m^2_s^-2_"
            cf_out *= ENERGY_UNITS[unit]^power
        elseif unit in keys(POWER_UNITS)
            units_out *= power != 1 ? "kg^$(Int64(power))_m^$(Int64(2*power))_s^$(Int64(-3*power))_" : "kg_m^2_s^-3_"
            cf_out *= POWER_UNITS[unit]^power
        elseif unit in keys(CHARGE_UNITS)
            units_out *= power != 1 ? "A^$(Int64(power))_s^$(Int64(power))_" : "A_s_"
            cf_out *= CHARGE_UNITS[unit]^power
        elseif unit in keys(VOLTAGE_UNITS)
            units_out *= power != 1 ? "kg^$(Int64(power))_m^$(Int64(2*power))_s^$(Int64(-3*power))_A^$(Int64(-power))_" : "kg_m^2_s^-3_A^-1_"
            cf_out *= VOLTAGE_UNITS[unit]^power
        else
            verbose && println("unit component '"*unit*"' not found in accepted OWCF units. Trying LONG format... ")
            unit_lc = lowercase(unit) # For long format, lowercase is required
            if unit_lc in keys(SPEED_UNITS_LONG)
                units_out *= power != 1 ? "m^$(Int64(power))_s^$(Int64(-power))_" : "m_s^-1_"
                cf_out *= SPEED_UNITS_LONG[unit_lc]^power
            elseif unit_lc in keys(ACCELERATION_UNITS_LONG)
                units_out *= power != 1 ? "m^$(Int64(power))_s^$(Int64(-2*power))_" : "m_s^-2_"
                cf_out *= ACCELERATION_UNITS_LONG[unit_lc]^power
            elseif unit_lc in keys(FORCE_UNITS_LONG)
                units_out *= power != 1 ? "kg^$(Int64(power))_m^$(Int64(power))_s^$(Int64(-2*power))_" : "kg_m_s^-2_"
                cf_out *= FORCE_UNITS_LONG[unit_lc]^power
            elseif unit_lc in keys(PRESSURE_UNITS_LONG)
                units_out *= power != 1 ? "kg^$(Int64(power))_m^$(Int64(-power))_s^$(Int64(-2*power))_" : "kg_m^-1_s^-2_"
                cf_out *= PRESSURE_UNITS_LONG[unit_lc]^power
            elseif unit_lc in keys(ENERGY_UNITS_LONG)
                units_out *= power != 1 ? "kg^$(Int64(power))_m^$(Int64(2*power))_s^$(Int64(-2*power))_" : "kg_m^2_s^-2_"
                cf_out *= ENERGY_UNITS_LONG[unit_lc]^power
            elseif unit_lc in keys(POWER_UNITS_LONG)
                units_out *= power != 1 ? "kg^$(Int64(power))_m^$(Int64(2*power))_s^$(Int64(-3*power))_" : "kg_m^2_s^-3_"
                cf_out *= POWER_UNITS_LONG[unit_lc]^power
            elseif unit_lc in keys(CHARGE_UNITS_LONG)
                units_out *= power != 1 ? "A^$(Int64(power))_s^$(Int64(power))_" : "A_s_"
                cf_out *= CHARGE_UNITS_LONG[unit_lc]^power
            elseif unit_lc in keys(VOLTAGE_UNITS_LONG)
                units_out *= power != 1 ? "kg^$(Int64(power))_m^$(Int64(2*power))_s^$(Int64(-3*power))_A^$(Int64(-power))_" : "kg_m^2_s^-3_A^-1_"
                cf_out *= VOLTAGE_UNITS_LONG[unit_lc]^power
            else
                error("unit componenent '"*unit*"' not found in accepted OWCF units. Please see OWCF/misc/convert_units.jl for lists of accepted units. Please correct and re-try.")
            end
        end
    end
    units_out = units_out[1:end-1] # Remove "_" after final units component

    # Merge duplicate units and return String
    return units_from_dict(units_to_dict(units_out)), cf_out
end
"""
    units_conversion_factor(units_in, units_out)

Take the units_in and units_out input variables, and compute the resulting 
unit conversion factor between the two. For example, if units_in="m" and units_out="km", then 
the function output would be 1.0e-3.

If keyword argument rounding is set to true, the output will be rounded to 12 significant digits.
If keyword argument verbose is set to true, the function will talk a lot!

PLEASE NOTE! The OWCF unit notation scheme is CASE SENSITIVE. This is to avoid ambiguity for units 
such as millielectronvolt (meV) and megaelectronvolt (MeV).
"""
function units_conversion_factor(units_in::String,units_out::String; rounding=true, verbose=false)

    # Convert units into base units (https://en.wikipedia.org/wiki/SI_base_unit). Extra conversion factor (cf) might be needed (e.g. 1 keV=(1.602177e-19) J = (1.602177e-19) kg_m^2_s^-2)
    units_in, cf_in = units_to_base_units(units_in; verbose=verbose)
    units_out, cf_out = units_to_base_units(units_out; verbose=verbose)

    # Decompose units into components. E.g. "m_s^-1" into ["m","s^-1"]
    unit_in_components = split(units_in,"_")
    unit_out_components = split(units_out,"_")

    # Process all units and powers in units_in
    unit_in_components_mag = Dict("TIME"=>1.0, "LENGTH"=>1.0, "MASS"=>1.0, "CURRENT"=>1.0, "TEMPERATURE"=>1.0, "SUBSTANCE"=>1.0, "LUMINOUSITY"=>1.0, "DIMENSIONLESS"=>1.0) # Magnitudes of units_in components. Since conversion factor is computed via ratios, initialize 1.0's
    unit_in_components_dim = Dict("TIME"=>0.0, "LENGTH"=>0.0, "MASS"=>0.0, "CURRENT"=>0.0, "TEMPERATURE"=>0.0, "SUBSTANCE"=>0.0, "LUMINOUSITY"=>0.0, "DIMENSIONLESS"=>0.0) # Dimensions of units_in components. Must match dimensions of units_out components (see below)
    for unit_in_component in unit_in_components
        if occursin("^",unit_in_component) # If unit has an exponent...
            unit, power = split(unit_in_component,"^") # e.g. "s^-1" into "s","-1"
            power = parse(Float64, power)
        else
            unit, power = unit_in_component, 1.0 # Otherwise, assume exponent of 1.0
        end
        if unit in keys(TIME_UNITS)
            if unit_in_components_mag["TIME"]==0.0 # If no unit is yet registered...
                unit_in_components_mag["TIME"] = 1.0 # ... set the count to 1.0 so that...
            end
            unit_in_components_mag["TIME"] *= TIME_UNITS[unit]^power # ... one can factor the unit magnitudes correctly
            unit_in_components_dim["TIME"] += power # Add the dimension
        elseif unit in keys(LENGTH_UNITS)
            if unit_in_components_mag["LENGTH"]==0.0 # If no unit magnitude is yet registered...
                unit_in_components_mag["LENGTH"] = 1.0 # ... set the count to 1.0 so that...
            end
            unit_in_components_mag["LENGTH"] *= LENGTH_UNITS[unit]^power # ... one can factor the unit magnitudes correctly
            unit_in_components_dim["LENGTH"] += power # Add the dimension
        elseif unit in keys(MASS_UNITS)
            if unit_in_components_mag["MASS"]==0.0 # If no unit magnitude is yet registered...
                unit_in_components_mag["MASS"] = 1.0 # ... set the count to 1.0 so that...
            end
            unit_in_components_mag["MASS"] *= MASS_UNITS[unit]^power # ... one can factor the unit magnitudes correctly
            unit_in_components_dim["MASS"] += power # Add the dimension
        elseif unit in keys(CURRENT_UNITS)
            if unit_in_components_mag["CURRENT"]==0.0 # If no unit magnitude is yet registered...
                unit_in_components_mag["CURRENT"] = 1.0 # ... set the count to 1.0 so that...
            end
            unit_in_components_mag["CURRENT"] *= CURRENT_UNITS[unit]^power # ... one can factor the unit magnitudes correctly
            unit_in_components_dim["CURRENT"] += power # Add the dimension
        elseif unit in keys(TEMPERATURE_UNITS)
            if unit_in_components_mag["TEMPERATURE"]==0.0 # If no unit magnitude is yet registered...
                unit_in_components_mag["TEMPERATURE"] = 1.0 # ... set the count to 1.0 so that...
            end
            unit_in_components_mag["TEMPERATURE"] *= TEMPERATURE_UNITS[unit]^power # ... one can factor the unit magnitudes correctly
            unit_in_components_dim["TEMPERATURE"] += power # Add the dimension
        elseif unit in keys(SUBSTANCE_UNITS)
            if unit_in_components_mag["SUBSTANCE"]==0.0 # If no unit magnitude is yet registered...
                unit_in_components_mag["SUBSTANCE"] = 1.0 # ... set the count to 1.0 so that...
            end
            unit_in_components_mag["SUBSTANCE"] *= SUBSTANCE_UNITS[unit]^power # ... one can factor the unit magnitudes correctly
            unit_in_components_dim["SUBSTANCE"] += power # Add the dimension
        elseif unit in keys(LUMINOUSITY_UNITS)
            if unit_in_components_mag["LUMINOUSITY"]==0.0 # If no unit magnitude is yet registered...
                unit_in_components_mag["LUMINOUSITY"] = 1.0 # ... set the count to 1.0 so that...
            end
            unit_in_components_mag["LUMINOUSITY"] *= LUMINOUSITY_UNITS[unit]^power # ... one can factor the unit magnitudes correctly
            unit_in_components_dim["LUMINOUSITY"] += power # Add the dimension
        elseif unit in keys(DIMENSIONLESS_UNITS)
            if unit_in_components_mag["DIMENSIONLESS"]==0.0 # If no unit magnitude is yet registered...
                unit_in_components_mag["DIMENSIONLESS"] = 1.0 # ... set the count to 1.0 so that...
            end
            unit_in_components_mag["DIMENSIONLESS"] *= DIMENSIONLESS_UNITS[unit]^power # ... one can factor the unit magnitudes correctly
            unit_in_components_dim["DIMENSIONLESS"] += power # Add the dimension
        else
            verbose && println("units_in component '"*unit*"' not found in accepted OWCF base units. Trying LONG format... ")
            unit_lc = lowercase(unit) # For long format, lowercase is required
            if unit_lc in keys(TIME_UNITS_LONG)
                if unit_in_components_mag["TIME"]==0.0 # If no unit magnitude is yet registered...
                    unit_in_components_mag["TIME"] = 1.0 # ... set the count to 1.0 so that...
                end
                unit_in_components_mag["TIME"] *= TIME_UNITS_LONG[unit_lc]^power # ... one can factor the unit magnitudes correctly
                unit_in_components_dim["TIME"] += power # Add the dimension
            elseif unit_lc in keys(LENGTH_UNITS_LONG)
                if unit_in_components_mag["LENGTH"]==0.0 # If no unit magnitude is yet registered...
                    unit_in_components_mag["LENGTH"] = 1.0 # ... set the count to 1.0 so that...
                end
                unit_in_components_mag["LENGTH"] *= LENGTH_UNITS_LONG[unit_lc]^power # ... one can factor the unit magnitudes correctly
                unit_in_components_dim["LENGTH"] += power # Add the dimension
            elseif unit_lc in keys(MASS_UNITS_LONG)
                if unit_in_components_mag["MASS"]==0.0 # If no unit magnitude is yet registered...
                    unit_in_components_mag["MASS"] = 1.0 # ... set the count to 1.0 so that...
                end
                unit_in_components_mag["MASS"] *= MASS_UNITS_LONG[unit_lc]^power # ... one can factor the unit magnitudes correctly
                unit_in_components_dim["MASS"] += power # Add the dimension
            elseif unit_lc in keys(CURRENT_UNITS_LONG)
                if unit_in_components_mag["CURRENT"]==0.0 # If no unit magnitude is yet registered...
                    unit_in_components_mag["CURRENT"] = 1.0 # ... set the count to 1.0 so that...
                end
                unit_in_components_mag["CURRENT"] *= CURRENT_UNITS_LONG[unit_lc]^power # ... one can factor the unit magnitudes correctly
                unit_in_components_dim["CURRENT"] += power # Add the dimension
            elseif unit_lc in keys(TEMPERATURE_UNITS_LONG)
                if unit_in_components_mag["TEMPERATURE"]==0.0 # If no unit magnitude is yet registered...
                    unit_in_components_mag["TEMPERATURE"] = 1.0 # ... set the count to 1.0 so that...
                end
                unit_in_components_mag["TEMPERATURE"] *= TEMPERATURE_UNITS_LONG[unit_lc]^power # ... one can factor the unit magnitudes correctly
                unit_in_components_dim["TEMPERATURE"] += power # Add the dimension
            elseif unit_lc in keys(SUBSTANCE_UNITS_LONG)
                if unit_in_components_mag["SUBSTANCE"]==0.0 # If no unit magnitude is yet registered...
                    unit_in_components_mag["SUBSTANCE"] = 1.0 # ... set the count to 1.0 so that...
                end
                unit_in_components_mag["SUBSTANCE"] *= SUBSTANCE_UNITS_LONG[unit_lc]^power # ... one can factor the unit magnitudes correctly
                unit_in_components_dim["SUBSTANCE"] += power # Add the dimension
            elseif unit_lc in keys(LUMINOUSITY_UNITS_LONG)
                if unit_in_components_mag["LUMINOUSITY"]==0.0 # If no unit magnitude is yet registered...
                    unit_in_components_mag["LUMINOUSITY"] = 1.0 # ... set the count to 1.0 so that...
                end
                unit_in_components_mag["LUMINOUSITY"] *= LUMINOUSITY_UNITS_LONG[unit_lc]^power # ... one can factor the unit magnitudes correctly
                unit_in_components_dim["LUMINOUSITY"] += power # Add the dimension
            elseif unit_lc in keys(DIMENSIONLESS_UNITS_LONG)
                if unit_in_components_mag["DIMENSIONLESS"]==0.0 # If no unit magnitude is yet registered...
                    unit_in_components_mag["DIMENSIONLESS"] = 1.0 # ... set the count to 1.0 so that...
                end
                unit_in_components_mag["DIMENSIONLESS"] *= DIMENSIONLESS_UNITS_LONG[unit_lc]^power # ... one can factor the unit magnitudes correctly
                unit_in_components_dim["DIMENSIONLESS"] += power # Add the dimension
            else
                error("units_in componenent '"*unit*"' not found in accepted OWCF units. Please see OWCF/misc/convert_units.jl for lists of accepted units. Please correct and re-try.")
            end
        end
    end

    # Process all units and powers in units_out
    unit_out_components_mag = Dict("TIME"=>1.0, "LENGTH"=>1.0, "MASS"=>1.0, "CURRENT"=>1.0, "TEMPERATURE"=>1.0, "SUBSTANCE"=>1.0, "LUMINOUSITY"=>1.0, "DIMENSIONLESS"=>1.0) # Magnitudes of units_in components. Since conversion factor is computed via ratios, initialize 1.0's
    unit_out_components_dim = Dict("TIME"=>0.0, "LENGTH"=>0.0, "MASS"=>0.0, "CURRENT"=>0.0, "TEMPERATURE"=>0.0, "SUBSTANCE"=>0.0, "LUMINOUSITY"=>0.0, "DIMENSIONLESS"=>0.0) # Dimensions of units_in components. Must match dimensions of units_out components (see below)
    for unit_out_component in unit_out_components
        if occursin("^",unit_out_component) # If unit has an exponent...
            unit, power = split(unit_out_component,"^") # e.g. "s^-1" into "s","-1"
            power = parse(Float64, power)
        else
            unit, power = unit_out_component, 1.0 # Otherwise, assume exponent of 1.0
        end
        if unit in keys(TIME_UNITS)
            if unit_out_components_mag["TIME"]==0.0 # If no unit is yet registered...
                unit_out_components_mag["TIME"] = 1.0 # ... set the count to 1.0 so that...
            end
            unit_out_components_mag["TIME"] *= TIME_UNITS[unit]^power # ... one can factor the unit magnitudes correctly
            unit_out_components_dim["TIME"] += power # Add the dimension
        elseif unit in keys(LENGTH_UNITS)
            if unit_out_components_mag["LENGTH"]==0.0 # If no unit magnitude is yet registered...
                unit_out_components_mag["LENGTH"] = 1.0 # ... set the count to 1.0 so that...
            end
            unit_out_components_mag["LENGTH"] *= LENGTH_UNITS[unit]^power # ... one can factor the unit magnitudes correctly
            unit_out_components_dim["LENGTH"] += power # Add the dimension
        elseif unit in keys(MASS_UNITS)
            if unit_out_components_mag["MASS"]==0.0 # If no unit magnitude is yet registered...
                unit_out_components_mag["MASS"] = 1.0 # ... set the count to 1.0 so that...
            end
            unit_out_components_mag["MASS"] *= MASS_UNITS[unit]^power # ... one can factor the unit magnitudes correctly
            unit_out_components_dim["MASS"] += power # Add the dimension
        elseif unit in keys(CURRENT_UNITS)
            if unit_out_components_mag["CURRENT"]==0.0 # If no unit magnitude is yet registered...
                unit_out_components_mag["CURRENT"] = 1.0 # ... set the count to 1.0 so that...
            end
            unit_out_components_mag["CURRENT"] *= CURRENT_UNITS[unit]^power # ... one can factor the unit magnitudes correctly
            unit_out_components_dim["CURRENT"] += power # Add the dimension
        elseif unit in keys(TEMPERATURE_UNITS)
            if unit_out_components_mag["TEMPERATURE"]==0.0 # If no unit magnitude is yet registered...
                unit_out_components_mag["TEMPERATURE"] = 1.0 # ... set the count to 1.0 so that...
            end
            unit_out_components_mag["TEMPERATURE"] *= TEMPERATURE_UNITS[unit]^power # ... one can factor the unit magnitudes correctly
            unit_out_components_dim["TEMPERATURE"] += power # Add the dimension
        elseif unit in keys(SUBSTANCE_UNITS)
            if unit_out_components_mag["SUBSTANCE"]==0.0 # If no unit magnitude is yet registered...
                unit_out_components_mag["SUBSTANCE"] = 1.0 # ... set the count to 1.0 so that...
            end
            unit_out_components_mag["SUBSTANCE"] *= SUBSTANCE_UNITS[unit]^power # ... one can factor the unit magnitudes correctly
            unit_out_components_dim["SUBSTANCE"] += power # Add the dimension
        elseif unit in keys(LUMINOUSITY_UNITS)
            if unit_out_components_mag["LUMINOUSITY"]==0.0 # If no unit magnitude is yet registered...
                unit_out_components_mag["LUMINOUSITY"] = 1.0 # ... set the count to 1.0 so that...
            end
            unit_out_components_mag["LUMINOUSITY"] *= LUMINOUSITY_UNITS[unit]^power # ... one can factor the unit magnitudes correctly
            unit_out_components_dim["LUMINOUSITY"] += power # Add the dimension
        elseif unit in keys(DIMENSIONLESS_UNITS)
            if unit_out_components_mag["DIMENSIONLESS"]==0.0 # If no unit magnitude is yet registered...
                unit_out_components_mag["DIMENSIONLESS"] = 1.0 # ... set the count to 1.0 so that...
            end
            unit_out_components_mag["DIMENSIONLESS"] *= DIMENSIONLESS_UNITS[unit]^power # ... one can factor the unit magnitudes correctly
            unit_out_components_dim["DIMENSIONLESS"] += power # Add the dimension
        else
            verbose && println("units_in component '"*unit*"' not found in accepted OWCF base units. Trying LONG format... ")
            unit_lc = lowercase(unit) # For long format, lowercase is required
            if unit_lc in keys(TIME_UNITS_LONG)
                if unit_out_components_mag["TIME"]==0.0 # If no unit magnitude is yet registered...
                    unit_out_components_mag["TIME"] = 1.0 # ... set the count to 1.0 so that...
                end
                unit_out_components_mag["TIME"] *= TIME_UNITS_LONG[unit_lc]^power # ... one can factor the unit magnitudes correctly
                unit_out_components_dim["TIME"] += power # Add the dimension
            elseif unit_lc in keys(LENGTH_UNITS_LONG)
                if unit_out_components_mag["LENGTH"]==0.0 # If no unit magnitude is yet registered...
                    unit_out_components_mag["LENGTH"] = 1.0 # ... set the count to 1.0 so that...
                end
                unit_out_components_mag["LENGTH"] *= LENGTH_UNITS_LONG[unit_lc]^power # ... one can factor the unit magnitudes correctly
                unit_out_components_dim["LENGTH"] += power # Add the dimension
            elseif unit_lc in keys(MASS_UNITS_LONG)
                if unit_out_components_mag["MASS"]==0.0 # If no unit magnitude is yet registered...
                    unit_out_components_mag["MASS"] = 1.0 # ... set the count to 1.0 so that...
                end
                unit_out_components_mag["MASS"] *= MASS_UNITS_LONG[unit_lc]^power # ... one can factor the unit magnitudes correctly
                unit_out_components_dim["MASS"] += power # Add the dimension
            elseif unit_lc in keys(CURRENT_UNITS_LONG)
                if unit_out_components_mag["CURRENT"]==0.0 # If no unit magnitude is yet registered...
                    unit_out_components_mag["CURRENT"] = 1.0 # ... set the count to 1.0 so that...
                end
                unit_out_components_mag["CURRENT"] *= CURRENT_UNITS_LONG[unit_lc]^power # ... one can factor the unit magnitudes correctly
                unit_out_components_dim["CURRENT"] += power # Add the dimension
            elseif unit_lc in keys(TEMPERATURE_UNITS_LONG)
                if unit_out_components_mag["TEMPERATURE"]==0.0 # If no unit magnitude is yet registered...
                    unit_out_components_mag["TEMPERATURE"] = 1.0 # ... set the count to 1.0 so that...
                end
                unit_out_components_mag["TEMPERATURE"] *= TEMPERATURE_UNITS_LONG[unit_lc]^power # ... one can factor the unit magnitudes correctly
                unit_out_components_dim["TEMPERATURE"] += power # Add the dimension
            elseif unit_lc in keys(SUBSTANCE_UNITS_LONG)
                if unit_out_components_mag["SUBSTANCE"]==0.0 # If no unit magnitude is yet registered...
                    unit_out_components_mag["SUBSTANCE"] = 1.0 # ... set the count to 1.0 so that...
                end
                unit_out_components_mag["SUBSTANCE"] *= SUBSTANCE_UNITS_LONG[unit_lc]^power # ... one can factor the unit magnitudes correctly
                unit_out_components_dim["SUBSTANCE"] += power # Add the dimension
            elseif unit_lc in keys(LUMINOUSITY_UNITS_LONG)
                if unit_out_components_mag["LUMINOUSITY"]==0.0 # If no unit magnitude is yet registered...
                    unit_out_components_mag["LUMINOUSITY"] = 1.0 # ... set the count to 1.0 so that...
                end
                unit_out_components_mag["LUMINOUSITY"] *= LUMINOUSITY_UNITS_LONG[unit_lc]^power # ... one can factor the unit magnitudes correctly
                unit_out_components_dim["LUMINOUSITY"] += power # Add the dimension
            elseif unit_lc in keys(DIMENSIONLESS_UNITS_LONG)
                if unit_out_components_mag["DIMENSIONLESS"]==0.0 # If no unit magnitude is yet registered...
                    unit_out_components_mag["DIMENSIONLESS"] = 1.0 # ... set the count to 1.0 so that...
                end
                unit_out_components_mag["DIMENSIONLESS"] *= DIMENSIONLESS_UNITS_LONG[unit_lc]^power # ... one can factor the unit magnitudes correctly
                unit_out_components_dim["DIMENSIONLESS"] += power # Add the dimension
            else
                error("units_in componenent '"*unit*"' not found in accepted OWCF units. Please see OWCF/misc/convert_units.jl for lists of accepted units. Please correct and re-try.")
            end
        end
    end

    # CHECK THAT DIMENSIONS MATCH
    if !(unit_in_components_dim==unit_out_components_dim)
        error("units_in ("*units_in*") and units_out ("*units_out*") do not have the same dimensions. Please correct and re-try.")
    end

    # COMPUTE CONVERSION FACTOR
    cf = cf_in/cf_out
    for key in keys(unit_in_components_mag)
        cf *= unit_in_components_mag[key]/unit_out_components_mag[key]
    end

    # RETURN CONVERSION FACTOR
    return rounding ? round(cf,sigdigits=12) : cf
end

"""
    units_to_dict(units_in)

Return a dictionary with all the components of units_in as keys,
and the powers as values. E.g. "m_s^-1" will result in an output of 
Dict("m"=>1, "s"=>-1).
"""
function units_to_dict(units_in::String)
    unit_in_components =  split(units_in,"_")
    
    my_unit_dict = Dict()

    for unit_in_component in unit_in_components
        if occursin("^",unit_in_component) # If unit has an exponent...
            unit, power = split(unit_in_component,"^") # e.g. "s^-1" into "s","-1"
            if haskey(my_unit_dict,unit)
                my_unit_dict[unit] += parse(Float64, power)
            else
                my_unit_dict[unit] = parse(Float64, power)
            end
        else
            if haskey(my_unit_dict,unit_in_component)
                my_unit_dict[unit_in_component] += 1.0 # Otherwise, assume exponent of 1.0
            else
                my_unit_dict[unit_in_component] = 1.0 # Otherwise, assume exponent of 1.0
            end
        end
    end

    return my_unit_dict
end

"""
    units_from_dict(my_dict)

Given a dictionary produced by the units_to_dict() function (or similar), return 
the resulting units of measurement as a string. For example, if we have that 
my_dict = Dict("kg"=>1.0, "s"=>-2.0, "m"=>-3.0) then the function output will be
"kg_s^-2_m^-3"
"""
function units_from_dict(my_dict::Dict)
    units_out = ""
    for key in keys(my_dict)
        if my_dict[key]!=1.0
            units_out *= "$(key)^$(Int64(my_dict[key]))_" 
        else
            units_out *= "$(key)_"
        end
    end
    return units_out[1:end-1] # Remove the last "_"
end