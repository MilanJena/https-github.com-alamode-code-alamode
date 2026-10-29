/*
 strain_dfc2_projection.h

 Copyright (c) 2026 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory
 or http://opensource.org/licenses/mit-license.php for information.
*/

// Projection of a real-space strain-derivative FC2 set onto index symmetry and the ASR.
// The set is dPhi2/du (Phi3 R) or d2Phi2/du2 (Phi4 R R), built by
// DerivativeIFC::compute_dV_dumn_all_real_space. Header-only (Eigen); no MPI.
//
// The derivative takes R relative to the first atom of each cluster. Exactly, the i- and
// j-centred anchorings agree by the ASR of Phi3/Phi4. With equidistant supercell images they do
// not, so the set M(i a, j b; R) is neither symmetric nor ASR-conserving. Two consequences:
// D(q) is not Hermitian, and the Hermitian solvers downstream read one triangle only. The
// projection is a chosen correction of the IFC-derived derivative, not a property of it:
//  1. M <- (M(i a, j b; R) + M(j b, i a; -R)) / 2. This averages the complete i- and j-centred
//     contributions; the j-centred one is the transposed entry of the permutation-expanded IFCs.
//  2. M <- M + delta with the smallest sum delta^2 over the support S that keeps (1) and gives
//     sum_{j,R} M(i a, j b; R) = 0 for every i, a, b (shear and antisymmetric parts included).
//     Here S is all Cartesian blocks of the atom pairs (i, j, R) of the set, closed under
//     transposition. The Lagrange solution is delta(i a, j b; R) = L_i(a, b) + L_j(b, a) on S.
//     L solves two nat x nat systems whose matrices only count the pairs of S, so the correction
//     depends on no force constant and no phonon spectrum.
// S comes from a symmetric IFC set and the norm is invariant, so the result is covariant under
// the space group. Both steps are linear and S is shared by all components, so any mixture
// sum_c u_c M_c is projected consistently. Working in real space keeps D(q) Hermitian, periodic
// and time-reversal symmetric at every q.

#pragma once

#include <Eigen/Core>
#include <Eigen/QR>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <map>
#include <vector>

namespace PHON_NS
{
struct StrainDfc2Entry
{
    int i, a, j, b;          // atoms of the primitive cell and Cartesian indices
    std::array<int, 3> R;    // lattice vector from the cell of i to that of j
    std::vector<double> val; // one value per strain component (Phi units)
};

// Returns the projected set: one entry per Cartesian block of every support pair, in a
// deterministic order. Input entries with the same key are summed. max_delta, if given,
// receives max |delta| of step 2; asr_residual the scaled ASR residual of the result
// (max |row sum| / max |M|), which is also warned about on std::cerr above 1e-8.
inline std::vector<StrainDfc2Entry> project_strain_dfc2(const std::vector<StrainDfc2Entry> &in, const int nat,
                                                        const std::size_t ncomp, double *max_delta = nullptr,
                                                        double *asr_residual = nullptr)
{
    using Key = std::array<int, 5>; // i, j, R
    std::map<Key, int> index;
    std::vector<Key> keys;
    auto add_pair = [&](const Key &k) {
        if (index.emplace(k, 0).second) keys.push_back(k);
    };
    for (const auto &e: in) {
        add_pair({e.i, e.j, e.R[0], e.R[1], e.R[2]});
        add_pair({e.j, e.i, -e.R[0], -e.R[1], -e.R[2]});
    }
    std::sort(keys.begin(), keys.end());
    for (std::size_t p = 0; p < keys.size(); ++p) index[keys[p]] = static_cast<int>(p);

    const auto npair = keys.size();
    std::vector<int> transposed(npair);
    for (std::size_t p = 0; p < npair; ++p) {
        const auto &k = keys[p];
        transposed[p] = index.at({k[1], k[0], -k[2], -k[3], -k[4]});
    }

    auto at = [ncomp](std::vector<double> &v, const std::size_t p, const int a, const int b) {
        return v.data() + (p * 9 + a * 3 + b) * ncomp;
    };

    std::vector<double> raw(npair * 9 * ncomp, 0.0);
    for (const auto &e: in) {
        const auto p = index.at({e.i, e.j, e.R[0], e.R[1], e.R[2]});
        auto *dst = at(raw, p, e.a, e.b);
        for (std::size_t c = 0; c < ncomp; ++c) dst[c] += e.val[c];
    }

    // 1. average with the transposed (j-centred) entry
    std::vector<double> sym(raw.size());
    for (std::size_t p = 0; p < npair; ++p) {
        for (int a = 0; a < 3; ++a) {
            for (int b = 0; b < 3; ++b) {
                const auto *x = at(raw, p, a, b);
                const auto *y = at(raw, transposed[p], b, a);
                auto *s = at(sym, p, a, b);
                for (std::size_t c = 0; c < ncomp; ++c) s[c] = 0.5 * (x[c] + y[c]);
            }
        }
    }

    // 2. min-norm symmetric correction restoring the row sums. The Lagrange equations for
    // L_i (3 x 3 per atom) read n_i L_i + sum_j N_ij L_j^T = -r_i (n_i: pairs of row i, N_ij: pairs
    // (i, j, R) over R). They decouple into the symmetric part S and the antisymmetric part A of L:
    //   (D + N) S = -sym(r),  (D - N) A = -anti(r),  D = diag(n_i),
    // and delta(i a, j b; R) = S_i(a, b) + S_j(a, b) + A_i(a, b) - A_j(a, b).
    // D - N is a graph Laplacian (null vector: A_i = common A) and D + N is singular on
    // bipartite components (S_i = +-S); those null vectors give delta = 0, the systems are
    // consistent, and every solution gives the same delta. The minimum-norm one is taken via
    // a complete orthogonal decomposition. Its rank threshold is Eigen's default: pivots below
    // nat * machine epsilon * the largest pivot count as zero. Pair counts are O(1)-O(1e3)
    // integers, so the gap to a true zero is wide.
    // ponytail: dense nat x nat solves, O(nat^3) time and O(nat^2) memory (about 1e4 atoms
    // before this costs ~1 GB); a sparse Laplacian solver is the upgrade path.
    Eigen::MatrixXd N = Eigen::MatrixXd::Zero(nat, nat);
    Eigen::MatrixXd rsym = Eigen::MatrixXd::Zero(nat, static_cast<Eigen::Index>(6 * ncomp));
    Eigen::MatrixXd ranti = Eigen::MatrixXd::Zero(nat, static_cast<Eigen::Index>(3 * ncomp));
    constexpr int sym_ab[6][2] = {{0, 0}, {1, 1}, {2, 2}, {1, 2}, {0, 2}, {0, 1}};
    constexpr int anti_ab[3][2] = {{1, 2}, {0, 2}, {0, 1}};
    std::vector<double> r(static_cast<std::size_t>(nat) * 9 * ncomp, 0.0); // row sums
    for (std::size_t p = 0; p < npair; ++p) {
        const auto i = keys[p][0];
        N(i, keys[p][1]) += 1.0;
        for (int ab = 0; ab < 9; ++ab) {
            const auto *s = sym.data() + (p * 9 + ab) * ncomp;
            auto *ri = r.data() + (static_cast<std::size_t>(i) * 9 + ab) * ncomp;
            for (std::size_t c = 0; c < ncomp; ++c) ri[c] += s[c];
        }
    }
    auto rv = [&](const int i, const int a, const int b, const std::size_t c) {
        return r[(static_cast<std::size_t>(i) * 9 + 3 * a + b) * ncomp + c];
    };
    for (int i = 0; i < nat; ++i) {
        for (std::size_t c = 0; c < ncomp; ++c) {
            for (int m = 0; m < 6; ++m) {
                const auto a = sym_ab[m][0], b = sym_ab[m][1];
                rsym(i, static_cast<Eigen::Index>(m * ncomp + c)) = -0.5 * (rv(i, a, b, c) + rv(i, b, a, c));
            }
            for (int m = 0; m < 3; ++m) {
                const auto a = anti_ab[m][0], b = anti_ab[m][1];
                ranti(i, static_cast<Eigen::Index>(m * ncomp + c)) = -0.5 * (rv(i, a, b, c) - rv(i, b, a, c));
            }
        }
    }
    const Eigen::MatrixXd Dg = N.rowwise().sum().asDiagonal();
    const Eigen::MatrixXd Sm = (Dg + N).completeOrthogonalDecomposition().solve(rsym);
    const Eigen::MatrixXd Am = (Dg - N).completeOrthogonalDecomposition().solve(ranti);
    // L_i(a, b) = S_i(a, b) + A_i(a, b) with S, A rebuilt from their independent components
    auto Lval = [&](const int i, const int a, const int b, const std::size_t c) {
        double v = 0.0;
        for (int m = 0; m < 6; ++m) {
            if ((sym_ab[m][0] == a && sym_ab[m][1] == b) || (sym_ab[m][0] == b && sym_ab[m][1] == a)) {
                v += Sm(i, static_cast<Eigen::Index>(m * ncomp + c));
            }
        }
        for (int m = 0; m < 3; ++m) {
            if (anti_ab[m][0] == a && anti_ab[m][1] == b) v += Am(i, static_cast<Eigen::Index>(m * ncomp + c));
            if (anti_ab[m][0] == b && anti_ab[m][1] == a) v -= Am(i, static_cast<Eigen::Index>(m * ncomp + c));
        }
        return v;
    };
    Eigen::MatrixXd L(9 * nat, static_cast<Eigen::Index>(ncomp));
    for (int i = 0; i < nat; ++i) {
        for (int a = 0; a < 3; ++a) {
            for (int b = 0; b < 3; ++b) {
                for (std::size_t c = 0; c < ncomp; ++c)
                    L(9 * i + 3 * a + b, static_cast<Eigen::Index>(c)) = Lval(i, a, b, c);
            }
        }
    }

    double dmax = 0.0;
    std::vector<StrainDfc2Entry> out;
    out.reserve(npair * 9);
    for (std::size_t p = 0; p < npair; ++p) {
        const auto &k = keys[p];
        for (int a = 0; a < 3; ++a) {
            for (int b = 0; b < 3; ++b) {
                StrainDfc2Entry e{k[0], a, k[1], b, {k[2], k[3], k[4]}, std::vector<double>(ncomp)};
                const auto *s = at(sym, p, a, b);
                for (std::size_t c = 0; c < ncomp; ++c) {
                    const auto ce = static_cast<Eigen::Index>(c);
                    const auto delta = L(9 * k[0] + 3 * a + b, ce) + L(9 * k[1] + 3 * b + a, ce);
                    dmax = std::max(dmax, std::abs(delta));
                    e.val[c] = s[c] + delta;
                }
                out.push_back(std::move(e));
            }
        }
    }
    // scaled ASR residual of the result; a rank decision gone wrong would show up here
    std::vector<double> rows(static_cast<std::size_t>(nat) * 9 * ncomp, 0.0);
    double vmax = 0.0;
    for (const auto &e: out) {
        for (std::size_t c = 0; c < ncomp; ++c) {
            rows[(static_cast<std::size_t>(e.i) * 9 + 3 * e.a + e.b) * ncomp + c] += e.val[c];
            vmax = std::max(vmax, std::abs(e.val[c]));
        }
    }
    double resid = 0.0;
    for (const auto x: rows) resid = std::max(resid, std::abs(x));
    resid = vmax > 0.0 ? resid / vmax : resid;
    if (resid > 1.0e-8) {
        std::cerr << " Warning: project_strain_dfc2: ASR residual / max|M| = " << resid
                  << " after the projection (expected round-off).\n";
    }
    if (max_delta) *max_delta = dmax;
    if (asr_residual) *asr_residual = resid;
    return out;
}
} // namespace PHON_NS
