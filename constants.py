# Useful physical constants, particle masses etc.
# Values are taken from NIST, unless explicitly stated differently.
# -----------------------------------------------------------------

# Speed of light in vacuum
c = 2.99792458e8          # m/s

# Elementary charge
e = 1.6021766208e-19      # C

# Vacuum permittivity
eps0 = 8.854187817e-12    # F/m

# Boltzmann constant
kb_J  = 1.38064852e-23    # J/K
kb_eV = 8.6173303e-5      # eV/K

# Atomic mass unit
u_keV = 931494.0954       # keV
u_kg  = 1.660539040e-27   # kg

# Electron mass
me = 5.48579909070e-4     # u

# Neutron mass
mn = 1.00866491588        # u

# Proton mass
mp = 1.007276466879       # u

# Deuteron mass
md = 2.013553212745       # u

# Triton mass
mt = 3.01550071632        # u

# 3-He mass
# (from NNDC Nuclear Wallet cards. Atomic binding energy neglected)
m3He = (3*u_keV + 14931.2177734375) / u_keV
m3He = m3He - 2*me        # u

# 4-He mass
m4He = 4.001506179127     # u

# 9-Be mass
# (from NNDC Nuclear Wallet cards. Atomic binding energy neglected)
m9Be = (9*u_keV + 11348.453125) / u_keV
m9Be = m9Be - 4*me        # u

# 12-C mass
# (from NNDC Nuclear Wallet cards. Atomic binding energy neglected)
m12C = 12*u_keV / u_keV
m12C = m12C - 6*me        # u

# 13-C mass
# (from NNDC Nuclear Wallet cards. Atomic binding energy neglected)
m13C = (13*u_keV + 3125.00888) / u_keV
m13C = m13C - 6*me        # u

# 10-B mass
# (from NNDC Nuclear Wallet cards. Atomic binding energy neglected)
m10B = (10*u_keV + 12050.609 ) / u_keV
m10B = m10B - 5*me        # u
