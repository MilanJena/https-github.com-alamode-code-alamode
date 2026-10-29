"""Physical constants used by the ALAMODE Python tools.

CODATA 2018 recommended values (SI 2019; h, c, e and k_B are exact), mirroring
include/constants.h so that the scripts convert units exactly as alm/anphon do.
tools/tests/test_constants.py checks that the two files agree.

Note: ``ase.units`` defaults to CODATA 2014; the difference is ~1e-8 relative.
"""

import math

# Primitive constants (SI)
AMU = 1.66053906660e-27  # kg
ELECTRON_MASS = 9.1093837015e-31  # kg
C_LIGHT = 299792458.0  # m/s
H_PLANCK = 6.62607015e-34  # J s
E_CHARGE = 1.602176634e-19  # C (also J per eV)
K_BOLTZMANN = 1.380649e-23  # J/K
RYD_IN_J = 4.3597447222071e-18 / 2.0  # J (half the Hartree energy)
BOHR_IN_ANGSTROM = 0.529177210903

# Derived factors, computed as in include/constants.h
RYD_IN_EV = RYD_IN_J / E_CHARGE
AMU_RY = AMU / ELECTRON_MASS / 2.0
TIME_RY = H_PLANCK / (2.0 * math.pi) / RYD_IN_J  # s
HZ_TO_KAYSER = 1.0e-2 / (2.0 * math.pi) / C_LIGHT
RY_TO_KAYSER = HZ_TO_KAYSER / TIME_RY  # Ry (energy) -> cm^-1
KAYSER_TO_RY = TIME_RY / HZ_TO_KAYSER  # cm^-1 -> Ry
KAYSER_TO_THZ = C_LIGHT * 1.0e-10  # cm^-1 -> THz
KAYSER_TO_MEV = C_LIGHT * 1.0e2 * H_PLANCK / E_CHARGE * 1.0e3  # cm^-1 -> meV
