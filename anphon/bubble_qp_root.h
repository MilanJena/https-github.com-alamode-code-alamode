/*
 bubble_qp_root.h

 Copyright (c) 2026 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory
 or http://opensource.org/licenses/mit-license.php for information.
*/

// Quasiparticle root of F(w) = w^2 - W^2 + 2 W Re Sigma(w) = 0 (BUBBLE = 3).
// Header-only so that test/test_bubble_qp_root.cpp can drive it with a synthetic Sigma.

#pragma once

#include <algorithm>
#include <cmath>
#include <set>
#include <vector>

namespace PHON_NS
{
namespace BubbleQP
{

struct Root
{
    double omega;  // root frequency
    double self;   // Re Sigma at the root, from the same linear model: W^2 - 2 W self = omega^2
    double weight; // Z = 2 omega / F'(omega)
};

struct Result
{
    std::vector<Root> roots; // upward roots (F' > 0), increasing frequency
    int best = -1;           // index into roots of the largest Z; -1: no root
    double self = 0.0;       // Re Sigma used: roots[best].self, or the on-shell Re Sigma(W)
};

// Vertex value of the parabola through three points (x0 < x1 < x2); sign of curvature in *a_out.
inline double parabola_extremum(const double x0, const double x1, const double x2, const double y0, const double y1,
                                const double y2, double *a_out)
{
    const auto s01 = (y1 - y0) / (x1 - x0);
    const auto a = ((y2 - y1) / (x2 - x1) - s01) / (x2 - x0);
    *a_out = a;
    if (a == 0.0) return y1;
    const auto x = 0.5 * (x0 + x1 - s01 / a); // zero of s01 + a (2x - x0 - x1)
    if (x <= x0 || x >= x2) return y1;
    return y0 + (x - x0) * s01 + a * (x - x0) * (x - x1);
}

// sigma(w) returns Re Sigma at every frequency of w and must return the same values on
// every caller (it may be an MPI collective). It is called at most twice.
//  1) coarse scan of [0, wmax] with step dcoarse, plus W itself for the fallback;
//  2) the coarse intervals with a sign change of F, and the two intervals around a coarse
//     local extremum whose parabola crosses zero (a pair of roots inside one step), are
//     rescanned with a step <= dfine;
//  3) every upward crossing (F' > 0) on the fine grid is solved with Sigma linear in the
//     interval, so the root and the frequency of the corrected FC2 coincide;
//  4) the root with the largest Z = 2w / F'(w) is used; none: the on-shell Re Sigma(W).
// A root at w = 0 (F(0) = 0, F > 0 after it) is a root with Z = w1^2 / F(w1) (F is even in w).
template <class SigmaFn>
Result solve(const double omega0, const double wmax, const double dcoarse, const double dfine, SigmaFn &&sigma)
{
    const auto func = [omega0](const double w, const double s) { return w * w - omega0 * omega0 + 2.0 * omega0 * s; };

    const auto ncoarse = static_cast<int>(std::ceil(wmax / dcoarse - 1.0e-9)) + 1;
    std::vector<double> wc(ncoarse + 1);
    for (auto i = 0; i < ncoarse; ++i) wc[i] = std::min(dcoarse * static_cast<double>(i), wmax);
    wc[ncoarse] = omega0;
    const auto sc = sigma(wc);
    std::vector<double> fc(ncoarse);
    for (auto i = 0; i < ncoarse; ++i) fc[i] = func(wc[i], sc[i]);

    Result res;
    res.self = sc[ncoarse];

    if (std::abs(fc[0]) <= 1.0e-12 * omega0 * omega0 && ncoarse > 1 && fc[1] > 0.0) {
        res.roots.push_back({0.0, sc[0], wc[1] * wc[1] / fc[1]});
    }

    std::set<int> intervals;
    for (auto i = 0; i + 1 < ncoarse; ++i) {
        if ((fc[i] < 0.0) != (fc[i + 1] < 0.0)) intervals.insert(i);
    }
    for (auto i = 1; i + 1 < ncoarse; ++i) {
        const auto is_min = fc[i] <= fc[i - 1] && fc[i] <= fc[i + 1] && fc[i] >= 0.0;
        const auto is_max = fc[i] >= fc[i - 1] && fc[i] >= fc[i + 1] && fc[i] < 0.0;
        if (!is_min && !is_max) continue;
        double a;
        const auto y = parabola_extremum(wc[i - 1], wc[i], wc[i + 1], fc[i - 1], fc[i], fc[i + 1], &a);
        if ((is_min && y < 0.0) || (is_max && y >= 0.0)) {
            intervals.insert(i - 1);
            intervals.insert(i);
        }
    }

    std::vector<double> wf;
    std::vector<int> start; // first index of each sub-grid in wf
    for (const auto i: intervals) {
        const auto h = wc[i + 1] - wc[i];
        const auto nsub = std::max(1, static_cast<int>(std::ceil(h / dfine - 1.0e-9)));
        start.push_back(static_cast<int>(wf.size()));
        for (auto j = 0; j <= nsub; ++j) wf.push_back(wc[i] + h * static_cast<double>(j) / nsub);
    }
    start.push_back(static_cast<int>(wf.size()));

    if (!wf.empty()) {
        const auto sf = sigma(wf);
        for (size_t g = 0; g + 1 < start.size(); ++g) {
            for (auto i = start[g]; i + 1 < start[g + 1]; ++i) {
                const auto f0 = func(wf[i], sf[i]);
                const auto f1 = func(wf[i + 1], sf[i + 1]);
                if (!(f0 < 0.0 && f1 >= 0.0)) continue;
                // Sigma = sf[i] + k (w - wf[i]):  w^2 + 2 W k w + c = 0, convex, one root in the
                // interval (the larger one, F' = 2 sqrt(disc) > 0).
                const auto k = (sf[i + 1] - sf[i]) / (wf[i + 1] - wf[i]);
                const auto p = omega0 * k;
                const auto c = 2.0 * omega0 * (sf[i] - k * wf[i]) - omega0 * omega0;
                const auto sq = std::sqrt(std::max(p * p - c, 0.0));
                auto w = p > 0.0 ? -c / (p + sq) : -p + sq;
                w = std::min(std::max(w, wf[i]), wf[i + 1]);
                res.roots.push_back({w, sf[i] + k * (w - wf[i]), sq > 0.0 ? w / sq : 0.0});
            }
        }
    }

    std::sort(res.roots.begin(), res.roots.end(), [](const Root &x, const Root &y) { return x.omega < y.omega; });
    for (size_t r = 0; r < res.roots.size(); ++r) {
        if (res.best < 0 || res.roots[r].weight > res.roots[res.best].weight) res.best = static_cast<int>(r);
    }
    if (res.best >= 0) res.self = res.roots[res.best].self;
    return res;
}

} // namespace BubbleQP
} // namespace PHON_NS
