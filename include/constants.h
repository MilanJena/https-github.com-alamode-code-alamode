/*
 constants.h

 Copyright (c) 2014 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory 
 or http://opensource.org/licenses/mit-license.php for information.
*/

#pragma once

#include <cfloat>
#include <cmath>
#include <complex>

constexpr double pi = M_PI;
constexpr double inverse_pi = 1.0 / pi; // pi^{-1}
constexpr double tpi = 2.0 * pi;        // (2*pi)
constexpr double inv_tpi = 1.0 / tpi;   // (2*pi)^{-1}
// Physical constants: CODATA 2018 recommended values (SI 2019; h, c, e, k_B exact).
// tools/alamode_constants.py and python/alm/fcsxml.py mirror these values;
// tools/tests/test_constants.py checks that they agree.
constexpr double amu = 1.66053906660e-27;          // kg
constexpr double electron_mass = 9.1093837015e-31; // kg
constexpr double amu_ry = amu / electron_mass / 2.0;
constexpr double c_light = 299792458;             // m/s
constexpr double h_planck = 6.62607015e-34;       // J s
constexpr double e_charge = 1.602176634e-19;      // C (also J per eV)
constexpr double Ryd = 4.3597447222071e-18 / 2.0; // J (half the Hartree energy)
constexpr double time_ry = h_planck * inv_tpi / Ryd;
constexpr double Hz_to_kayser = 1.0e-2 * inv_tpi / c_light;
constexpr double Ry_to_kayser = Hz_to_kayser / time_ry; // Ry (energy) -> cm^-1
constexpr double kayser_to_Ry = time_ry / Hz_to_kayser; // cm^-1 -> Ry
constexpr double kayser_to_THz = c_light * 1.0e-10;     // cm^-1 -> THz

inline double in_kayser(const double x)
{
    return x * Ry_to_kayser;
}

constexpr double Bohr_in_Angstrom = 0.529177210903;
constexpr double Ryd_in_eV = Ryd / e_charge;
constexpr double k_Boltzmann = 1.380649e-23; // J/K
constexpr auto im = std::complex<double>(0.0, 1.0);
constexpr double eps = DBL_EPSILON;
constexpr double eps15 = 1.0e-15;
constexpr double eps12 = 1.0e-12;
constexpr double eps10 = 1.0e-10;
constexpr double eps8 = 1.0e-8;
constexpr double eps6 = 1.0e-6;
constexpr double eps4 = 1.0e-4;
