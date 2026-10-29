// Unit test for the BUBBLE = 3 quasiparticle root search (anphon/bubble_qp_root.h) on
// synthetic F(w) = w^2 - W^2 + 2 W Sigma(w): Sigma is built from a prescribed F.
// Build: target test_bubble_qp_root (anphon/CMakeLists.txt); run without arguments.

#include <cmath>
#include <cstdio>
#include <functional>
#include <vector>
#include "bubble_qp_root.h"

using namespace PHON_NS;

namespace
{
int failures = 0;

void check(const bool ok, const char *what)
{
    if (!ok) {
        std::printf("FAIL: %s\n", what);
        ++failures;
    }
}

constexpr double W = 100.0;    // SCPH frequency
constexpr double WMAX = 150.0; // W + max(50, W/2)
constexpr double COARSE = 2.5; // EPSILON/4 for EPSILON = 10
constexpr double FINE = 0.1;

BubbleQP::Result run(const std::function<double(double)> &f, int *ncalls)
{
    *ncalls = 0;
    const auto sigma = [&](const std::vector<double> &w) {
        ++*ncalls;
        std::vector<double> s(w.size());
        for (size_t i = 0; i < w.size(); ++i) s[i] = (f(w[i]) - w[i] * w[i] + W * W) / (2.0 * W);
        return s;
    };
    return BubbleQP::solve(W, WMAX, COARSE, FINE, sigma);
}

double used_omega(const BubbleQP::Result &r)
{
    return std::sqrt(W * W - 2.0 * W * r.self);
}
} // namespace

int main()
{
    int n;

    // single root, constant Sigma: exact
    auto r = run([](double w) { return w * w - 3600.0; }, &n);
    check(r.roots.size() == 1 && r.best == 0, "single root: count");
    check(std::abs(r.roots[0].omega - 60.0) < 1e-9 && std::abs(r.roots[0].weight - 1.0) < 1e-9,
          "single root: 60, Z = 1");
    check(n == 2, "single root: two Sigma evaluations");

    // the root of the linear-Sigma model is the frequency of the corrected FC2
    r = run([](double w) { return w * w - 3600.0 + 300.0 * std::sin(w / 3.0); }, &n);
    check(r.best >= 0, "wiggly Sigma: a root");
    if (r.best >= 0) {
        const auto &q = r.roots[r.best];
        check(std::abs(used_omega(r) - q.omega) < 1e-9, "root = sqrt(W^2 - 2 W Sigma(root))");
        check(std::abs(q.omega * q.omega - 3600.0 + 300.0 * std::sin(q.omega / 3.0)) < 1.0,
              "wiggly Sigma: F(root) ~ 0");
    }

    // several upward roots (30, 80; 50 is downward): the largest Z = 2w / F' is used
    // (Z from the secant slope of one fine interval: O(h) error)
    r = run([](double w) { return 0.01 * (w - 30.0) * (w - 50.0) * (w - 80.0); }, &n);
    check(r.roots.size() == 2, "cubic: two upward roots");
    check(r.best == 1 && std::abs(r.roots[1].omega - 80.0) < 1e-2, "cubic: root at 80 (Z 10.7 > 6)");
    check(std::abs(r.roots[0].weight - 6.0) < 0.15 && std::abs(r.roots[1].weight - 160.0 / 15.0) < 0.15,
          "cubic: Z values");

    // a down/up pair 0.6 apart inside one coarse step, all coarse values positive
    r = run([](double w) { return (w - 71.0) * (w - 71.0) - 0.09; }, &n);
    check(r.roots.size() == 1 && std::abs(r.roots[0].omega - 71.3) < 1e-2, "close pair: upward root at 71.3");

    // an up/down pair hidden below zero between coarse points, plus the main root
    r = run(
        [](double w) {
            return w < 80.0 ? 0.09 - (w - 71.0) * (w - 71.0) : -80.91 * (1.0 - (w - 80.0) * (w - 80.0) / 841.0);
        },
        &n);
    check(r.roots.size() == 2 && std::abs(r.roots[0].omega - 70.7) < 1e-2, "hidden pair below zero: root 70.7");

    // root at w = 0: F = w^2 (Z = 1)
    r = run([](double w) { return w * w; }, &n);
    check(r.best == 0 && r.roots[0].omega == 0.0 && std::abs(r.roots[0].weight - 1.0) < 1e-9, "root at 0");
    check(std::abs(r.self - 0.5 * W) < 1e-12, "root at 0: Sigma(0)");

    // no root (mode unstable at all w): on-shell fallback
    r = run([](double w) { return w * w + 100.0; }, &n);
    check(r.best < 0 && r.roots.empty(), "no root");
    check(std::abs(r.self - 100.0 / (2.0 * W) - W / 2.0) < 1e-12, "no root: on-shell Sigma(W)");
    check(n == 1, "no root: one Sigma evaluation");

    if (failures == 0) std::printf("test_bubble_qp_root: all checks passed\n");
    return failures == 0 ? 0 : 1;
}
