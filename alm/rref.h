#pragma once

#include <limits>
#include "fcs.h"

auto rref(const size_t nrows, const size_t ncols, double **mat, size_t &nrank,
          const double tolerance = 1.0e-12) -> void;

auto rref(std::vector<std::vector<double>> &mat, const double tolerance = 1.0e-12) -> void;

auto rref_sparse(const size_t ncols, ConstraintSparseForm &sp_constraint, const double tolerance = 1.0e-12) -> void;

// Partial-pivot RREF preserving the pivot-column order used by
// Constraint::get_mapping_constraint, with improved numerical stability.
// Threshold decisions of rref_sparse_pivot: smallest accepted pivot, largest rejected column
// maximum, and whether a pivot row kept residue left of its pivot (drop_small = false only).
struct RrefPivotStats
{
    double min_accepted = std::numeric_limits<double>::max();
    double max_rejected = 0.0;
    bool left_residue = false;
};

// drop_small = false skips the final |a| <= tolerance pruning (for chained reductions).
auto rref_sparse_pivot(const size_t ncols, ConstraintSparseForm &sp_constraint, const double tolerance = 1.0e-12,
                       const bool drop_small = true, RrefPivotStats *stats = nullptr) -> void;

// Same result as rref_sparse_pivot, but reduces rows [0, nfirst) first and the rest in their
// free columns; much less fill when the first block is local (symmetry constraints).
// Falls back to rref_sparse_pivot (returns false) if a threshold decision is borderline.
auto rref_sparse_pivot_two_stage(const size_t ncols, ConstraintSparseForm &sp_constraint, const size_t nfirst,
                                 const double tolerance = 1.0e-12) -> bool;

// Gauss-Jordan elimination of independent rows with threshold row pivoting (|a| >= 0.1 * row maximum,
// leftmost preferred); columns >= nparam hold the constant part and are never pivots. pivot_col[k] is the
// pivot (entry 1) of row k; no row holds another row's pivot. Rows below tolerance * their original
// maximum after elimination are dropped (ndropped) unless they have a constant part; rhs_left is set if
// a row with a constant part has only round-off left (inconsistent).
auto rref_threshold_pivot(const size_t nparam, ConstraintSparseForm &sp_constraint, const double tolerance,
                          std::vector<size_t> &pivot_col, size_t &ndropped, bool &rhs_left) -> void;

// void rref_sparse2(const size_t ncols,
//                   ConstraintSparseForm &sp_constraint,
//                   const double tolerance = 1.0e-12);
//
// void rref_sparse_fullpivot(const size_t ncols,
//                            ConstraintSparseForm &sp_constraint,
//                            const double tolerance = 1.0e-12);
