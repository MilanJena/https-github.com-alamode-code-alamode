/*
 ifc_derivative_real_space.cpp

 Copyright (c) 2026 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory
 or http://opensource.org/licenses/mit-license.php for information.

 The static real-space kernels of DerivativeIFC (strain and displacement
 derivatives of IFC lists). They depend on the IFC lists only, not on System or
 the k-meshes, so unit tests can link this file alone.
*/

#include <algorithm>
#include <boost/sort/block_indirect_sort/block_indirect_sort.hpp>
#include <cmath>
#include <map>
#include <vector>
#include "constants.h"
#include "error.h"
#include "ifc_derivative.h"
#include "mathfunctions.h"

using namespace PHON_NS;

namespace
{

// Iterate over fcs_aligned (sorted by the first nelems indices), grouping
// consecutive entries with identical heading indices and relative vectors.
// accept(entry) filters entries; accumulate(entry) adds an accepted entry's
// contribution; flush(index_with_cell, atoms_s, relvecs, relvecs_velocity) is
// called once per closed group (and must reset the accumulator state).
template <class Accept, class Accumulate, class Flush>
void scan_fcs_groups(const std::vector<FcsArrayWithCell> &fcs_aligned, const std::size_t nelems, Accept &&accept,
                     Accumulate &&accumulate, Flush &&flush)
{
    std::vector<int> index_prev(nelems, -1), index_curr;
    std::vector<int> relvecs_int_prev(3 * (nelems - 1), 1000000), relvecs_int_curr;
    std::vector<int> index_with_cell_prev(3 * (nelems - 1) + 1, -1), index_with_cell_curr;
    std::vector<unsigned int> atoms_s_prev, atoms_s_curr;
    std::vector<Eigen::Vector3d> relvecs_prev(nelems - 1), relvecs_curr(nelems - 1);
    std::vector<Eigen::Vector3d> relvecs_vel_prev(nelems - 1), relvecs_vel_curr(nelems - 1);

    bool has_group = false;

    for (const auto &it: fcs_aligned) {

        if (!accept(it)) {
            continue;
        }

        index_curr.clear();
        relvecs_int_curr.clear();
        index_with_cell_curr.clear();
        atoms_s_curr.clear();

        index_curr.push_back(it.pairs[0].index);
        index_with_cell_curr.push_back(it.pairs[0].index);
        atoms_s_curr.emplace_back(it.atoms_s[0]);

        for (std::size_t i = 1; i < nelems; ++i) {
            index_curr.push_back(it.pairs[i].index);

            for (int j = 0; j < 3; ++j) {
                relvecs_int_curr.push_back(nint(it.relvecs[i - 1][j]));
            }

            index_with_cell_curr.push_back(it.pairs[i].index);
            index_with_cell_curr.push_back(it.pairs[i].tran);
            index_with_cell_curr.push_back(it.pairs[i].cell_s);

            atoms_s_curr.emplace_back(it.atoms_s[i]);
            relvecs_curr[i - 1] = it.relvecs[i - 1];
            relvecs_vel_curr[i - 1] = it.relvecs_velocity[i - 1];
        }

        if (!has_group || index_curr != index_prev || relvecs_int_curr != relvecs_int_prev) {
            if (has_group) {
                flush(index_with_cell_prev, atoms_s_prev, relvecs_prev, relvecs_vel_prev);
            }

            index_prev = index_curr;
            relvecs_int_prev = relvecs_int_curr;
            atoms_s_prev = atoms_s_curr;
            relvecs_prev = relvecs_curr;
            relvecs_vel_prev = relvecs_vel_curr;
            index_with_cell_prev = index_with_cell_curr;
            has_group = true;
        }

        accumulate(it);
    }

    if (has_group) {
        flush(index_with_cell_curr, atoms_s_curr, relvecs_curr, relvecs_vel_curr);
    }
}

// Rebuild the emitted pairs array from the packed index_with_cell key: the
// first atom is folded to the primitive cell (tran = cell_s = 0), the rest
// keep their translation and cell indices.
void build_pairs_vec(const std::vector<int> &index_with_cell, const std::size_t nelems,
                     std::vector<AtomCellSuper> &pairs_vec)
{
    AtomCellSuper pairs_tmp{};

    pairs_vec.clear();
    pairs_tmp.index = index_with_cell[0];
    pairs_tmp.tran = 0;
    pairs_tmp.cell_s = 0;
    pairs_vec.push_back(pairs_tmp);

    for (std::size_t i = 1; i < nelems; ++i) {
        pairs_tmp.index = index_with_cell[3 * i - 2];
        pairs_tmp.tran = index_with_cell[3 * i - 1];
        pairs_tmp.cell_s = index_with_cell[3 * i];
        pairs_vec.push_back(pairs_tmp);
    }
}

} // namespace

void DerivativeIFC::strain_fc5_parts(const std::vector<DeltaFcsStrainComponents> &groups, const double *u_tensor,
                                     const bool diag, const std::vector<double> &u0,
                                     std::vector<FcsArrayWithCell> &dfc4, std::vector<FcsArrayWithCell> &dfc3)
{
    dfc4.clear();
    dfc3.clear();
    if (groups.empty()) return;

    std::vector<std::pair<std::size_t, double>> terms;
    for (std::size_t mn = 0; mn < 9; ++mn) {
        if (diag && mn / 3 != mn % 3) continue;
        terms.emplace_back(mn, u_tensor[mn]);
    }
    extract_strain_combination(groups, terms, 1, eps15, dfc4);

    // dPhi4 : u0 on the home leg, as RealSpaceV4::q0_sweep contracts it: the strain moment
    // of each entry is measured from its home atom, which a contraction of another leg
    // would move by the separation (the two differ by separation x the FC5 ASR residual).
    // The rest, (j; l, m), is recentred on j: entries are stored with the first leg in
    // the home cell, relvecs and relvecs_velocity relative to it.
    std::map<std::vector<long>, std::size_t> slot;
    for (const auto &e: dfc4) {
        const auto w = e.fcs_val * u0[e.pairs[0].index];
        if (w == 0.0) continue;
        std::vector<AtomCellSuper> pairs(e.pairs.begin() + 1, e.pairs.end());
        std::vector<unsigned int> atoms_s(e.atoms_s.begin() + 1, e.atoms_s.end());
        const std::vector<Eigen::Vector3d> relvecs{e.relvecs[1] - e.relvecs[0], e.relvecs[2] - e.relvecs[0]};
        const std::vector<Eigen::Vector3d> relvecs_vel{e.relvecs_velocity[1] - e.relvecs_velocity[0],
                                                       e.relvecs_velocity[2] - e.relvecs_velocity[0]};
        std::vector<long> key;
        for (const auto &p: pairs) key.push_back(p.index);
        for (const auto &v: relvecs) {
            for (auto d = 0; d < 3; ++d) key.push_back(nint(v[d]));
        }
        const auto [it, added] = slot.emplace(key, dfc3.size());
        if (added) {
            dfc3.emplace_back(w, pairs, atoms_s, relvecs, relvecs_vel);
        } else {
            dfc3[it->second].fcs_val += w;
        }
    }
    dfc3.erase(std::remove_if(dfc3.begin(),
                              dfc3.end(),
                              [](const FcsArrayWithCell &e) { return std::abs(e.fcs_val) <= eps15; }),
               dfc3.end());
}

void DerivativeIFC::compute_dV_dumn_all_real_space(const std::vector<FcsArrayWithCell> &fcs_aligned,
                                                   std::vector<DeltaFcsStrainComponents> &groups, const std::size_t m,
                                                   const Eigen::Matrix3d &convmat)
{
    // Compute all 9^m strain derivatives in one pass. Each IFC contributes
    // fcs_val * prod_j vec_j[nu_j] to the 3^m components with its fixed mu indices.
    // fcs_aligned must be sorted by its first n-m indices (n is the IFC order).
    groups.clear();

    if (fcs_aligned.empty()) return;

    // touched_mu packs one bit per mu-combination, so 3^m must fit in 32 bits.
    if (m == 0 || m > 3) {
        exit("compute_dV_dumn_all_real_space", "Derivative order m must be 1, 2, or 3.");
    }

    const auto norder = fcs_aligned[0].pairs.size();

    if (m >= norder) {
        exit("compute_dV_dumn_all_real_space", "Derivative order m must be smaller than IFC order n.");
    }

    const auto nelems = norder - m;

    constexpr std::size_t pow3[4] = {1, 3, 9, 27};
    const auto ncomp = pow3[m] * pow3[m];

    std::vector<double> acc(ncomp, 0.0);
    uint32_t touched_mu = 0;
    std::vector<AtomCellSuper> pairs_vec;

    scan_fcs_groups(
        fcs_aligned,
        nelems,
        [](const FcsArrayWithCell &) { return true; },
        [&](const FcsArrayWithCell &it) {
            Eigen::Vector3d vecs[3];
            std::size_t mu[3];
            std::size_t p = 0;
            for (std::size_t j = 0; j < m; ++j) {
                mu[j] = it.coords[nelems + j];
                vecs[j] = convmat * it.relvecs_velocity[(nelems - 1) + j];
                p = p * 3 + mu[j];
            }
            touched_mu |= 1u << p;

            for (std::size_t q = 0; q < pow3[m]; ++q) {
                auto term = it.fcs_val;
                std::size_t c = 0;
                for (std::size_t j = 0; j < m; ++j) {
                    const auto nu = (q / pow3[m - 1 - j]) % 3;
                    term *= vecs[j][nu];
                    c = c * 9 + mu[j] * 3 + nu;
                }
                acc[c] += term;
            }
        },
        [&](const std::vector<int> &index_with_cell,
            const std::vector<unsigned int> &atoms_s,
            const std::vector<Eigen::Vector3d> &relvecs,
            const std::vector<Eigen::Vector3d> &relvecs_vel) {
            build_pairs_vec(index_with_cell, nelems, pairs_vec);

            groups.emplace_back();
            auto &group = groups.back();
            group.pairs = pairs_vec;
            group.atoms_s = atoms_s;
            group.relvecs = relvecs;
            group.relvecs_velocity = relvecs_vel;
            group.values = acc;
            group.touched_mu = touched_mu;

            std::fill(acc.begin(), acc.end(), 0.0);
            touched_mu = 0;
        });
}

void DerivativeIFC::extract_strain_component(const std::vector<DeltaFcsStrainComponents> &groups,
                                             const std::size_t component, const std::size_t m,
                                             const double emit_threshold, std::vector<FcsArrayWithCell> &delta_fcs)
{
    extract_strain_combination(groups, {{component, 1.0}}, m, emit_threshold, delta_fcs);
}

void DerivativeIFC::extract_strain_combination(const std::vector<DeltaFcsStrainComponents> &groups,
                                               const std::vector<std::pair<std::size_t, double>> &terms,
                                               const std::size_t m, const double emit_threshold,
                                               std::vector<FcsArrayWithCell> &delta_fcs)
{
    delta_fcs.clear();

    constexpr std::size_t pow9[4] = {1, 9, 81, 729};

    if (m == 0 || m > 3 || terms.empty()) {
        exit("extract_strain_combination", "Invalid derivative order or empty term list.");
    }

    // The mu-combination bits of the terms, matching the touched_mu convention.
    uint32_t mu_mask = 0;
    for (const auto &term: terms) {
        if (term.first >= pow9[m]) {
            exit("extract_strain_combination", "Invalid strain component index.");
        }
        std::size_t p = 0;
        for (std::size_t j = 0; j < m; ++j) {
            p = p * 3 + ((term.first / pow9[m - 1 - j]) % 9) / 3;
        }
        mu_mask |= 1u << p;
    }

    for (const auto &group: groups) {
        if (!(group.touched_mu & mu_mask)) continue;

        double val = 0.0;
        for (const auto &term: terms) {
            val += term.second * group.values[term.first];
        }
        if (emit_threshold >= 0.0 && std::abs(val) <= emit_threshold) continue;

        delta_fcs.emplace_back(val, group.pairs, group.atoms_s, group.relvecs, group.relvecs_velocity);
    }
}

void DerivativeIFC::compute_dV_dsublattice_real_space(const std::vector<FcsArrayWithCell> &fcs_aligned,
                                                      std::vector<FcsArrayWithCell> &delta_fcs,
                                                      const Eigen::VectorXd &sublattice_displacement,
                                                      const double emit_threshold)
{
    delta_fcs.clear();

    if (fcs_aligned.empty()) return;

    const auto norder = fcs_aligned[0].pairs.size();

    if (norder < 2) {
        exit("compute_dV_dsublattice_real_space", "IFC order must be at least 2.");
    }

    const auto nelems = norder - 1;

    double fcs_tmp = 0.0;
    std::vector<AtomCellSuper> pairs_vec;

    scan_fcs_groups(
        fcs_aligned,
        nelems,
        [](const FcsArrayWithCell &) { return true; },
        [&](const FcsArrayWithCell &it) {
            fcs_tmp += it.fcs_val * sublattice_displacement(it.pairs[norder - 1].index);
        },
        [&](const std::vector<int> &index_with_cell,
            const std::vector<unsigned int> &atoms_s,
            const std::vector<Eigen::Vector3d> &relvecs,
            const std::vector<Eigen::Vector3d> &relvecs_vel) {
            if (!(emit_threshold >= 0.0 && std::abs(fcs_tmp) <= emit_threshold)) {
                build_pairs_vec(index_with_cell, nelems, pairs_vec);
                delta_fcs.emplace_back(fcs_tmp, pairs_vec, atoms_s, relvecs, relvecs_vel);
            }
            fcs_tmp = 0.0;
        });
}

void DerivativeIFC::compute_dV_ddeform_real_space(const std::vector<FcsArrayWithCell> &fcs_aligned,
                                                  std::vector<FcsArrayWithCell> &delta_fcs,
                                                  const std::vector<DeformationField> &fields,
                                                  const Eigen::Matrix3d &convmat, const double emit_threshold)
{
    // Contraction of the last m legs with the displacement fields
    // d_lambda(l kappa) = sum_nu u(lambda, nu) R_nu + S(3*kappa + lambda),
    // one field per contracted leg. The relative-vector form of the affine
    // part is exact under the acoustic sum rule at every order m.
    delta_fcs.clear();

    if (fcs_aligned.empty()) return;

    const auto m = fields.size();

    if (m == 0) {
        exit("compute_dV_ddeform_real_space", "Inconsistent contraction-order input.");
    }

    const auto norder = fcs_aligned[0].pairs.size();

    if (m >= norder) {
        exit("compute_dV_ddeform_real_space", "Contraction order m must be smaller than IFC order n.");
    }

    const auto nelems = norder - m;

    auto compute_term = [&](const FcsArrayWithCell &it) {
        double term = it.fcs_val;
        for (std::size_t j = 0; j < m; ++j) {
            const auto index_tail = it.pairs[nelems + j].index;
            const auto mu = index_tail % 3;
            const Eigen::Vector3d vec = convmat * it.relvecs_velocity[(nelems - 1) + j];
            double disp = 0.0;
            for (auto nu = 0; nu < 3; ++nu) {
                disp += fields[j].displacement_gradient(mu, nu) * vec[nu];
            }
            if (fields[j].sublattice_displacement.size() != 0) {
                disp += fields[j].sublattice_displacement(index_tail);
            }
            term *= disp;
        }
        return term;
    };

    double fcs_tmp = 0.0;
    std::vector<AtomCellSuper> pairs_vec;

    scan_fcs_groups(
        fcs_aligned,
        nelems,
        [](const FcsArrayWithCell &) { return true; },
        [&](const FcsArrayWithCell &it) { fcs_tmp += compute_term(it); },
        [&](const std::vector<int> &index_with_cell,
            const std::vector<unsigned int> &atoms_s,
            const std::vector<Eigen::Vector3d> &relvecs,
            const std::vector<Eigen::Vector3d> &relvecs_vel) {
            if (!(emit_threshold >= 0.0 && std::abs(fcs_tmp) <= emit_threshold)) {
                build_pairs_vec(index_with_cell, nelems, pairs_vec);
                delta_fcs.emplace_back(fcs_tmp, pairs_vec, atoms_s, relvecs, relvecs_vel);
            }
            fcs_tmp = 0.0;
        });
}

void DerivativeIFC::compute_deformed_ifcs(const std::vector<const std::vector<FcsArrayWithCell> *> &fcs_by_order,
                                          const std::size_t target_order_index,
                                          const Eigen::Matrix3d &displacement_gradient,
                                          const Eigen::VectorXd &sublattice_displacement,
                                          const Eigen::Matrix3d &convmat, std::vector<FcsArrayWithCell> &fcs_deformed)
{
    if (target_order_index >= fcs_by_order.size()) {
        exit("compute_deformed_ifcs", "The requested IFC order is not available.");
    }

    fcs_deformed = *fcs_by_order[target_order_index];

    DeformationField field;
    field.displacement_gradient = displacement_gradient;
    field.sublattice_displacement = sublattice_displacement;

    auto inv_factorial = 1.0;
    for (std::size_t mth = 1; target_order_index + mth < fcs_by_order.size(); ++mth) {

        inv_factorial /= static_cast<double>(mth);

        const auto &fcs_higher = *fcs_by_order[target_order_index + mth];
        if (fcs_higher.empty()) continue;

        auto fcs_aligned = fcs_higher;
        sort_by_heading_indices const sorter(static_cast<unsigned int>(mth));
        boost::sort::block_indirect_sort(fcs_aligned.begin(), fcs_aligned.end(), sorter);

        std::vector<FcsArrayWithCell> delta;
        const std::vector<DeformationField> fields(mth, field);
        compute_dV_ddeform_real_space(fcs_aligned, delta, fields, convmat, eps15);

        for (auto &it: delta) {
            it.fcs_val *= inv_factorial;
        }
        fcs_deformed.insert(fcs_deformed.end(), delta.begin(), delta.end());
    }
}

void DerivativeIFC::compute_dV_dstrain_real_space(const std::vector<FcsArrayWithCell> &fcs_aligned,
                                                  std::vector<FcsArrayWithCell> &delta_fcs,
                                                  const std::vector<Eigen::Matrix3d> &strain_dirs,
                                                  const Eigen::Matrix3d &convmat, const double emit_threshold)
{
    // Compute real-space directional strain derivatives along strain_dirs[j].
    // With m = strain_dirs.size() and n the IFC order, fcs_aligned must be
    // sorted by its first n-m indices.
    if (fcs_aligned.empty()) {
        delta_fcs.clear();
        return;
    }

    const auto m = strain_dirs.size();

    if (m == 0) {
        exit("compute_dV_dstrain_real_space", "Inconsistent derivative-order input.");
    }

    const auto norder = fcs_aligned[0].pairs.size();

    if (m >= norder) {
        exit("compute_dV_dstrain_real_space", "Derivative order m must be smaller than IFC order n.");
    }

    const auto nelems = norder - m;

    delta_fcs.clear();

    auto tail_mu_matches = [&](const FcsArrayWithCell &it) {
        for (std::size_t j = 0; j < m; ++j) {
            const auto mu = it.coords[nelems + j];
            if (strain_dirs[j](mu, 0) == 0.0 && strain_dirs[j](mu, 1) == 0.0 && strain_dirs[j](mu, 2) == 0.0) {
                return false;
            }
        }
        return true;
    };

    auto compute_term = [&](const FcsArrayWithCell &it) {
        double term = it.fcs_val;
        for (std::size_t j = 0; j < m; ++j) {
            const auto mu = it.coords[nelems + j];
            const Eigen::Vector3d vec = convmat * it.relvecs_velocity[(nelems - 1) + j];
            double contracted = 0.0;
            for (int nu = 0; nu < 3; ++nu) {
                contracted += strain_dirs[j](mu, nu) * vec[nu];
            }
            term *= contracted;
        }
        return term;
    };

    double fcs_tmp = 0.0;
    std::vector<AtomCellSuper> pairs_vec;

    scan_fcs_groups(
        fcs_aligned,
        nelems,
        tail_mu_matches,
        [&](const FcsArrayWithCell &it) { fcs_tmp += compute_term(it); },
        [&](const std::vector<int> &index_with_cell,
            const std::vector<unsigned int> &atoms_s,
            const std::vector<Eigen::Vector3d> &relvecs,
            const std::vector<Eigen::Vector3d> &relvecs_vel) {
            if (!(emit_threshold >= 0.0 && std::abs(fcs_tmp) <= emit_threshold)) {
                build_pairs_vec(index_with_cell, nelems, pairs_vec);
                delta_fcs.emplace_back(fcs_tmp, pairs_vec, atoms_s, relvecs, relvecs_vel);
            }
            fcs_tmp = 0.0;
        });
}
