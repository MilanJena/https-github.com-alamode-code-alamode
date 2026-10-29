#include "rref.h"
#include <algorithm>
#include <cmath>
#include <iterator>
#include <vector>
#include "constraint.h"

#ifdef _OPENMP
#include <omp.h>
#endif

auto rref(const size_t nrows, const size_t ncols, double **mat, size_t &nrank, const double tolerance) -> void
{
    // Return the reduced row echelon form (rref) of matrix mat.
    // In addition, rank of the matrix is estimated.

    size_t jcol;
    double tmp;

    nrank = 0;

    size_t icol = 0;

    for (size_t irow = 0; irow < nrows; ++irow) {

        auto pivot = irow;

        while (std::abs(mat[pivot][icol]) < tolerance) {
            ++pivot;

            if (pivot == nrows) {
                pivot = irow;
                ++icol;

                if (icol == ncols) break;
            }
        }

        if (icol == ncols) break;

        if (std::abs(mat[pivot][icol]) > tolerance) ++nrank;

        if (pivot != irow) {
            //#pragma omp parallel for private(tmp)
            for (jcol = icol; jcol < ncols; ++jcol) {
                tmp = mat[pivot][jcol];
                mat[pivot][jcol] = mat[irow][jcol];
                mat[irow][jcol] = tmp;
            }
        }

        tmp = mat[irow][icol];
        tmp = 1.0 / tmp;
        //#pragma omp parallel for
        for (jcol = icol; jcol < ncols; ++jcol) {
            mat[irow][jcol] *= tmp;
        }

        for (auto jrow = 0; jrow < nrows; ++jrow) {
            if (jrow == irow) continue;

            tmp = mat[jrow][icol];
            //#pragma omp parallel for
            for (jcol = icol; jcol < ncols; ++jcol) {
                mat[jrow][jcol] -= tmp * mat[irow][jcol];
            }
        }
    }
}


auto rref(std::vector<std::vector<double>> &mat, const double tolerance) -> void
{
    // Return the reduced row echelon form (rref) of matrix mat.
    // In addition, rank of the matrix is estimated.

    size_t jcol;
    double tmp;

    size_t nrank = 0;
    size_t icol = 0;

    const auto nrows = mat.size();
    const auto ncols = mat[0].size();

    for (size_t irow = 0; irow < nrows; ++irow) {

        auto pivot = irow;

        while (std::abs(mat[pivot][icol]) < tolerance) {
            ++pivot;

            if (pivot == nrows) {
                pivot = irow;
                ++icol;

                if (icol == ncols) break;
            }
        }

        if (icol == ncols) break;

        if (std::abs(mat[pivot][icol]) > tolerance) ++nrank;

        if (pivot != irow) {
            for (jcol = icol; jcol < ncols; ++jcol) {
                tmp = mat[pivot][jcol];
                mat[pivot][jcol] = mat[irow][jcol];
                mat[irow][jcol] = tmp;
            }
        }

        tmp = mat[irow][icol];
        tmp = 1.0 / tmp;
        for (jcol = icol; jcol < ncols; ++jcol) {
            mat[irow][jcol] *= tmp;
        }

        for (size_t jrow = 0; jrow < nrows; ++jrow) {
            if (jrow == irow) continue;

            tmp = mat[jrow][icol];
            for (jcol = icol; jcol < ncols; ++jcol) {
                mat[jrow][jcol] -= tmp * mat[irow][jcol];
            }
        }
    }

    mat.erase(mat.begin() + nrank, mat.end());
    mat.shrink_to_fit();
}


auto rref_sparse(const size_t ncols, ConstraintSparseForm &sp_constraint, const double tolerance) -> void
{
    // This function is somewhat sensitive to the numerical accuracy.
    // The loss of numerical digits can lead to instability.
    // Column ordering may improve the stability, but I'm not sure.
    // Smaller tolerance is preferable.

    const auto nrows = sp_constraint.size();
    size_t jrow;

    // This parameter controls the stability and performance.
    // Smaller value is more stable but a little more costly.
    //    double zero_criterion = tolerance * 1.0e-3;
    constexpr auto zero_criterion = eps10;

    size_t nrank = 0;
    size_t icol = 0;
    std::vector<size_t> pivot_col;
    MapConstraintElement::iterator it_other;
    MapConstraintElement::iterator it_elem;

    for (size_t irow = 0; irow < nrows; ++irow) {

        auto pivot = irow;

        while (true) {
            it_elem = sp_constraint[pivot].find(icol);
            if (it_elem != sp_constraint[pivot].end()) {
                if (std::abs(it_elem->second) >= tolerance) {
                    break;
                }
            }

            ++pivot;
            if (pivot == nrows) {
                pivot = irow;
                ++icol;

                if (icol == ncols) break;
            }
        }

        if (icol == ncols) break;

        if (std::abs(it_elem->second) >= tolerance) ++nrank;
        pivot_col.push_back(icol);

        if (pivot != irow) {
            std::iter_swap(sp_constraint.begin() + irow, sp_constraint.begin() + pivot);
        }

        const double division_factor = 1.0 / it_elem->second;
        for (auto &[fst, snd]: sp_constraint[irow]) {
            snd *= division_factor;
        }

        for (jrow = 0; jrow < nrows; ++jrow) {
            if (jrow == irow) continue;

            it_elem = sp_constraint[jrow].find(icol);
            if (it_elem == sp_constraint[jrow].end()) continue;
            const double scaling_factor = it_elem->second;

            // Subtract irow elements from jrow
            for (const auto &[fst, snd]: sp_constraint[irow]) {
                // This part might be speeded up by using std::map::lower_bound
                // when the datatype of sp_constraint[irow] is std::map.
                if (fst < icol) {
                    continue;
                }
                it_other = sp_constraint[jrow].find(fst);
                if (it_other != sp_constraint[jrow].end()) {
                    it_other->second -= scaling_factor * snd;
                    // Delete zero elements and remove from the map.
                    // A smaller threshold is used for better stability.
                    if (std::abs(it_other->second) < zero_criterion) {
                        sp_constraint[jrow].erase(it_other);
                    }
                } else {
                    sp_constraint[jrow][fst] = -scaling_factor * snd;
                }
            }
            // Make sure to erase the icol element from the target row if it exists.
            // When the original pivot element is large, the element after subtraction can sometimes be
            // larger than the tolerance value because of the loss of significant digits.
            it_other = sp_constraint[jrow].find(icol);
            if (it_other != sp_constraint[jrow].end()) {
                sp_constraint[jrow].erase(it_other);
            }
        }
    }

    // Erase all elements smaller than the tolerance value, and the residue left of the pivot
    // (columns without a usable pivot), which the pivot normalization may have scaled above it.
    for (jrow = 0; jrow < nrows; ++jrow) {
        const auto pivot = jrow < pivot_col.size() ? pivot_col[jrow] : ncols;
        it_other = sp_constraint[jrow].begin();
        while (it_other != sp_constraint[jrow].end()) {
            if (std::abs(it_other->second) <= tolerance || it_other->first < pivot) {
                sp_constraint[jrow].erase(it_other++);
            } else {
                ++it_other;
            }
        }
    }

    // Remove emptry entries from the sp_constraint vector
    sp_constraint.erase(std::remove_if(sp_constraint.begin(),
                                       sp_constraint.end(),
                                       [](const MapConstraintElement &obj) { return obj.empty(); }),
                        sp_constraint.end());
    sp_constraint.shrink_to_fit();
}


auto rref_sparse_pivot(const size_t ncols, ConstraintSparseForm &sp_constraint, const double tolerance,
                       const bool drop_small, RrefPivotStats *stats) -> void
{
    // Gaussian elimination with maximum-magnitude row pivoting, then back-substitution to RREF.
    // Preserves the left-to-right pivot-column order used by the constraint map
    // while reducing round-off compared with first-acceptable-pivot RREF.

    const auto nrows = sp_constraint.size();
    if (nrows == 0) return;

    // Threshold for dropping numerically-zero fill, mirroring rref_sparse().
    constexpr auto zero_criterion = eps10;
    constexpr auto cancellation_criterion = 1.0e-14;

    size_t icol = 0;
    MapConstraintElement::iterator it_other;

    // Lazy column -> row occupancy index over stable row ids, so that pivot search and
    // elimination visit only rows containing icol. Lists may hold stale ids or duplicates;
    // every visit re-checks the row with find(). Rows are still swapped physically, so
    // id_at[position] / pos_of[id] track which row sits where.
    std::vector<size_t> id_at(nrows), pos_of(nrows);
    std::vector<std::vector<size_t>> col_rows(ncols);
    for (size_t r = 0; r < nrows; ++r) {
        id_at[r] = r;
        pos_of[r] = r;
        for (const auto &[fst, snd]: sp_constraint[r]) {
            if (fst < ncols) col_rows[fst].push_back(r);
        }
    }
    const auto release_column = [&](const size_t c) {
        col_rows[c].clear();
        col_rows[c].shrink_to_fit();
    };

    std::vector<size_t> pivot_col;
    std::vector<size_t> targets;
    std::vector<std::pair<size_t, double>> pivot_entries;
#ifdef _OPENMP
    std::vector<std::vector<std::pair<size_t, size_t>>> fill(omp_get_max_threads());
#else
    std::vector<std::vector<std::pair<size_t, size_t>>> fill(1);
#endif
    for (size_t irow = 0; irow < nrows; ++irow) {

        // Find the leftmost column (>= icol) that has a usable pivot among rows [irow, nrows),
        // and within that column pick the row with the maximum-magnitude entry
        // (ties -> smallest position, as in an ascending scan with strict '>').
        size_t pivot = nrows;
        double pivot_abs = 0.0;
        while (icol < ncols) {
            pivot = nrows;
            pivot_abs = 0.0;
            for (const auto id: col_rows[icol]) {
                const auto r = pos_of[id];
                if (r < irow) continue;
                const auto it = sp_constraint[r].find(icol);
                if (it != sp_constraint[r].end()) {
                    const auto a = std::abs(it->second);
                    if (a > pivot_abs || (a == pivot_abs && pivot != nrows && r < pivot)) {
                        pivot_abs = a;
                        pivot = r;
                    }
                }
            }
            if (pivot_abs >= tolerance) break; // pivot column found
            if (stats) stats->max_rejected = std::max(stats->max_rejected, pivot_abs);
            release_column(icol); // column has no pivot -> free column
            ++icol;
        }
        if (icol == ncols) break;
        if (stats) stats->min_accepted = std::min(stats->min_accepted, pivot_abs);

        // Move the chosen pivot row to position irow.
        if (pivot != irow) {
            std::iter_swap(sp_constraint.begin() + irow, sp_constraint.begin() + pivot);
            std::swap(id_at[irow], id_at[pivot]);
            pos_of[id_at[irow]] = irow;
            pos_of[id_at[pivot]] = pivot;
        }

        // Normalize the pivot row so the pivot entry becomes 1.
        const double division_factor = 1.0 / sp_constraint[irow].find(icol)->second;
        for (auto &[fst, snd]: sp_constraint[irow]) {
            snd *= division_factor;
        }

        // Eliminate column icol from the rows below containing it. Rows are independent and each
        // row sees the same operations in the same order, so the result does not depend on the
        // thread count. Rows above are reduced by back-substitution at the end: eliminating them
        // here fills them with columns that later pivots remove again.
        targets.clear();
        for (const auto id: col_rows[icol]) {
            if (pos_of[id] > irow) targets.push_back(pos_of[id]);
        }
        std::sort(targets.begin(), targets.end()); // col_rows may list a row twice
        targets.erase(std::unique(targets.begin(), targets.end()), targets.end());
        pivot_entries.clear();
        for (const auto &[fst, snd]: sp_constraint[irow]) {
            if (fst >= icol) pivot_entries.emplace_back(fst, snd);
        }
        const auto ntargets = static_cast<long>(targets.size());
#ifdef _OPENMP
#pragma omp parallel if (ntargets > 64 && ntargets * static_cast<long>(pivot_entries.size()) > 20000)
#endif
        {
#ifdef _OPENMP
            auto &new_fill = fill[omp_get_thread_num()];
#else
            auto &new_fill = fill[0];
#endif
#ifdef _OPENMP
#pragma omp for schedule(dynamic, 16)
#endif
            for (long t = 0; t < ntargets; ++t) {
                auto &row = sp_constraint[targets[t]];
                const auto it_piv = row.find(icol);
                if (it_piv == row.end()) continue;
                const double scaling_factor = it_piv->second;

                for (const auto &[fst, snd]: pivot_entries) {
                    const auto it = row.find(fst);
                    if (it != row.end()) {
                        it->second -= scaling_factor * snd;
                        if (std::abs(it->second) < zero_criterion) row.erase(it);
                    } else {
                        row[fst] = -scaling_factor * snd;
                        if (fst < ncols) new_fill.emplace_back(fst, id_at[targets[t]]);
                    }
                }
                // Make sure the pivot column is exactly cleared from this row.
                row.erase(icol);
            }
        }
        for (auto &f: fill) {
            for (const auto &[c, id]: f) col_rows[c].push_back(id);
            f.clear();
        }
        release_column(icol);
        pivot_col.push_back(icol);
    }

    // Back-substitution, last pivot first: the later pivot rows are already reduced, so one pass
    // per row removes all pivot columns right of its own pivot.
    {
        std::vector<size_t> row_of(ncols, nrows);
        for (size_t k = 0; k < pivot_col.size(); ++k) row_of[pivot_col[k]] = k;
        std::vector<std::pair<size_t, double>> hits;
        for (auto k = pivot_col.size(); k-- > 0;) {
            auto &row = sp_constraint[k];
            hits.clear();
            for (const auto &[col, val]: row) {
                if (col > pivot_col[k] && col < ncols && row_of[col] != nrows) hits.emplace_back(col, val);
            }
            for (const auto &[col, val]: hits) {
                row.erase(col);
                // As in the forward phase, skip source entries left of its pivot (residue in
                // columns without a usable pivot); only its free columns enter.
                for (const auto &[c, v]: sp_constraint[row_of[col]]) {
                    if (c > col) row[c] -= val * v;
                }
            }
            // Row k is final now: drop cancellation round-off so it stays sparse as a source for
            // the rows above. Only round-off: small but real coefficients get amplified there, so
            // the eps10 fill criterion (or dropping during the accumulation) loses accuracy.
            if (!hits.empty()) {
                for (auto it = row.begin(); it != row.end();) {
                    it = std::abs(it->second) < cancellation_criterion ? row.erase(it) : std::next(it);
                }
            }
        }
    }

    // Drop the sub-tolerance residue left of the pivot (columns without a usable pivot); rows below
    // the last pivot hold only such residue. The pivot normalization may have scaled the residue
    // above the tolerance, and a kept entry left of the pivot would be taken for the pivot by
    // Constraint::get_mapping_constraint.
    sp_constraint.resize(pivot_col.size());
    for (size_t jrow = 0; jrow < pivot_col.size(); ++jrow) {
        auto &row = sp_constraint[jrow];
        for (auto it = row.begin(); it != row.end();) {
            if (it->first < pivot_col[jrow]) {
                if (stats) stats->left_residue = true;
                it = row.erase(it);
            } else {
                ++it;
            }
        }
    }
    // Without drop_small, keep the pivot rows unrounded.
    if (!drop_small) return;

    // Erase all elements smaller than the tolerance value.
    for (size_t jrow = 0; jrow < sp_constraint.size(); ++jrow) {
        it_other = sp_constraint[jrow].begin();
        while (it_other != sp_constraint[jrow].end()) {
            if (std::abs(it_other->second) <= tolerance) {
                sp_constraint[jrow].erase(it_other++);
            } else {
                ++it_other;
            }
        }
    }

    // Remove empty entries from the sp_constraint vector.
    sp_constraint.erase(std::remove_if(sp_constraint.begin(),
                                       sp_constraint.end(),
                                       [](const MapConstraintElement &obj) { return obj.empty(); }),
                        sp_constraint.end());
    sp_constraint.shrink_to_fit();
}

namespace
{
// Replace every pivot column of `row` by its expression in the pivot rows (row_of[col] -> pivot row,
// pivot entry 1). The pivot rows hold no other pivot column, so one pass suffices.
void substitute_pivots(MapConstraintElement &row, const std::vector<const MapConstraintElement *> &row_of,
                       const double zero_criterion)
{
    std::vector<std::pair<size_t, double>> hits;
    for (const auto &[col, val]: row) {
        if (row_of[col]) hits.emplace_back(col, val);
    }
    for (const auto &[col, val]: hits) {
        row.erase(col);
        for (const auto &[c, v]: *row_of[col]) {
            if (c == col) continue;
            auto &x = row[c];
            x -= val * v;
            if (std::abs(x) < zero_criterion) row.erase(c);
        }
    }
}
} // namespace

auto rref_sparse_pivot_two_stage(const size_t ncols, ConstraintSparseForm &sp_constraint, const size_t nfirst,
                                 const double tolerance) -> bool
{
    // Same reduced row echelon form as rref_sparse_pivot (the RREF of a row space is unique),
    // but the first nfirst rows (local, e.g. symmetry) are reduced on their own first and the
    // remaining rows are reduced in the space of their free columns. A single left-to-right
    // pass mixes both blocks and its fill can exceed the final map by orders of magnitude.
    // The equivalence holds only when no threshold decision is borderline: a pivot or a rejected
    // column near the tolerance, or residue in a skipped column (which normalization can make
    // significant), may be decided differently by the two orders. Then fall back to one pass.
    // Returns false if it fell back.
    if (nfirst == 0 || nfirst >= sp_constraint.size()) {
        rref_sparse_pivot(ncols, sp_constraint, tolerance);
        return true;
    }
    constexpr auto zero_criterion = eps10;
    constexpr auto margin = 100.0;
    const auto ambiguous = [&](const RrefPivotStats &st) {
        return st.left_residue || st.min_accepted < margin * tolerance || st.max_rejected > tolerance / margin;
    };
    const auto fallback = [&](ConstraintSparseForm &&original) {
        sp_constraint = std::move(original);
        rref_sparse_pivot(ncols, sp_constraint, tolerance);
        return false;
    };
    ConstraintSparseForm original = sp_constraint;

    ConstraintSparseForm rest(std::make_move_iterator(sp_constraint.begin() + nfirst),
                              std::make_move_iterator(sp_constraint.end()));
    sp_constraint.resize(nfirst);
    RrefPivotStats stats;
    rref_sparse_pivot(ncols, sp_constraint, tolerance, false, &stats);
    if (ambiguous(stats)) return fallback(std::move(original));

    const auto pivot_index = [ncols](const ConstraintSparseForm &rows) {
        std::vector<const MapConstraintElement *> row_of(ncols, nullptr);
        for (const auto &row: rows) {
            size_t p = ncols;
            for (const auto &kv: row) p = std::min(p, kv.first);
            row_of[p] = &row;
        }
        return row_of;
    };

    auto row_of = pivot_index(sp_constraint);
    for (auto &row: rest) substitute_pivots(row, row_of, zero_criterion);
    stats = RrefPivotStats();
    rref_sparse_pivot(ncols, rest, tolerance, false, &stats);
    if (ambiguous(stats)) return fallback(std::move(original));

    row_of = pivot_index(rest);
    for (auto &row: sp_constraint) substitute_pivots(row, row_of, zero_criterion);

    sp_constraint.insert(sp_constraint.end(),
                         std::make_move_iterator(rest.begin()),
                         std::make_move_iterator(rest.end()));
    for (auto &row: sp_constraint) {
        for (auto it = row.begin(); it != row.end();) {
            it = std::abs(it->second) <= tolerance ? row.erase(it) : std::next(it);
        }
    }
    sp_constraint.erase(std::remove_if(sp_constraint.begin(),
                                       sp_constraint.end(),
                                       [](const MapConstraintElement &obj) { return obj.empty(); }),
                        sp_constraint.end());
    // Ascending pivot columns, as rref_sparse_pivot returns them.
    std::vector<std::pair<size_t, size_t>> order;
    for (size_t i = 0; i < sp_constraint.size(); ++i) {
        size_t p = ncols;
        for (const auto &kv: sp_constraint[i]) p = std::min(p, kv.first);
        order.emplace_back(p, i);
    }
    std::sort(order.begin(), order.end());
    ConstraintSparseForm sorted;
    sorted.reserve(order.size());
    for (const auto &pi: order) sorted.push_back(std::move(sp_constraint[pi.second]));
    sp_constraint = std::move(sorted);
    return true;
}

auto rref_threshold_pivot(const size_t nparam, ConstraintSparseForm &sp_constraint, const double tolerance,
                          std::vector<size_t> &pivot_col, size_t &ndropped, bool &rhs_left) -> void
{
    // Row-sequential Gauss-Jordan elimination with threshold pivoting. Each row, with the earlier pivots
    // substituted, takes its pivot among the entries >= pivot_threshold * (row maximum), leftmost first,
    // and the pivot column is eliminated from the earlier pivot rows. Columns >= nparam (constant part)
    // are never pivots. A row whose entries all fall below tolerance * (its original maximum) is
    // dependent and dropped, unless it keeps a constant part; rhs_left is set if such a row has
    // only round-off left in its parameter columns (inconsistent).
    constexpr auto pivot_threshold = 0.1;
    constexpr auto cancellation_criterion = 1.0e-14;
    constexpr auto machine_level = 1.0e-12;
    const auto nrows = sp_constraint.size();
    const auto row_max = [nparam](const MapConstraintElement &row) {
        double m = 0.0;
        for (const auto &[c, v]: row) {
            if (c < nparam) m = std::max(m, std::abs(v));
        }
        return m;
    };

    std::vector<const MapConstraintElement *> row_of(nparam + 1, nullptr);
    std::vector<std::vector<size_t>> col_rows(nparam); // pivot rows that may hold the column (lazy)
    std::vector<size_t> pivot_of_row(nrows, nparam);
    std::vector<size_t> targets;
    ndropped = 0;
    rhs_left = false;
    for (size_t irow = 0; irow < nrows; ++irow) {
        auto &row = sp_constraint[irow];
        const auto max_orig = row_max(row);
        substitute_pivots(row, row_of, cancellation_criterion * max_orig);
        const auto max_now = row_max(row);
        if (max_now <= tolerance * max_orig) {
            // A small row with a constant part is kept unless it is round-off: it is a real
            // constraint that the earlier rows do not satisfy (select_independent_rows kept it).
            const auto has_rhs = std::any_of(row.begin(), row.end(), [&](const auto &kv) {
                return kv.first >= nparam && std::abs(kv.second) > tolerance * max_orig;
            });
            if (!has_rhs || max_now <= machine_level * max_orig) {
                rhs_left = rhs_left || has_rhs;
                row.clear();
                ++ndropped;
                continue;
            }
        }
        size_t pivot = nparam;
        for (auto it = row.begin(); it != row.end();) {
            const auto a = std::abs(it->second);
            if (a < cancellation_criterion * max_now) {
                it = row.erase(it);
                continue;
            }
            if (it->first < pivot && a >= pivot_threshold * max_now) pivot = it->first;
            ++it;
        }
        const auto division_factor = 1.0 / row[pivot];
        for (auto &[c, v]: row) v *= division_factor;
        row[pivot] = 1.0;

        targets.clear();
        for (const auto r: col_rows[pivot]) targets.push_back(r);
        std::sort(targets.begin(), targets.end());
        targets.erase(std::unique(targets.begin(), targets.end()), targets.end());
        for (const auto r: targets) {
            auto &prow = sp_constraint[r];
            const auto it_piv = prow.find(pivot);
            if (it_piv == prow.end()) continue;
            const auto scaling_factor = it_piv->second;
            prow.erase(it_piv);
            for (const auto &[c, v]: row) {
                if (c == pivot) continue;
                const auto it = prow.find(c);
                if (it == prow.end()) {
                    prow[c] = -scaling_factor * v;
                    if (c < nparam) col_rows[c].push_back(r);
                } else {
                    const auto old = it->second;
                    it->second -= scaling_factor * v;
                    if (std::abs(it->second) < cancellation_criterion * std::abs(old)) prow.erase(it);
                }
            }
        }
        col_rows[pivot].clear();
        col_rows[pivot].shrink_to_fit();
        for (const auto &[c, v]: row) {
            if (c < nparam && c != pivot) col_rows[c].push_back(irow);
        }
        row_of[pivot] = &row;
        pivot_of_row[irow] = pivot;
    }

    pivot_col.clear();
    size_t nkept = 0;
    for (size_t irow = 0; irow < nrows; ++irow) {
        if (pivot_of_row[irow] == nparam) continue;
        if (nkept != irow) sp_constraint[nkept] = std::move(sp_constraint[irow]);
        pivot_col.push_back(pivot_of_row[irow]);
        ++nkept;
    }
    sp_constraint.resize(nkept);
}
