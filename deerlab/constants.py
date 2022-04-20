import numpy as np

# Mathematical constants
π = np.pi 

# Fundamental physical constants
Nav = 6.02214076e23      # Avogadro constant, mol^-1
μB = 9.2740100783e-24    # Bohr magneton, J/T (CODATA 2018 value)
μ0 = 1.25663706212e-6    # magnetic constant, N A^-2 = T^2 m^3 J^-1 (CODATA 2018)
h = 6.62607015e-34       # Planck constant, J/Hz (CODATA 2018)
ge = 2.00231930436256    # free-electron g factor (CODATA 2018 value)
hbar = h/2/π             # reduced Planck constant, J/(rad/s)

# Derived physical constants
D = (μ0/4/π)*(μB*ge)**2/hbar   # dipolar constant, rad s^-1 m^3
