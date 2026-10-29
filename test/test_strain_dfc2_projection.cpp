// Unit test of the projection of a strain-derivative FC2 set (anphon/strain_dfc2_projection.h)
// on a random synthetic set. The projected set must:
//  - satisfy the index symmetry M(i a, j b; R) = M(j b, i a; -R), so D(q) is Hermitian at a generic q;
//  - satisfy the ASR per atom row and Cartesian block, per component;
//  - be left unchanged by a second projection;
//  - commute with a signed-permutation rotation plus an atom permutation (covariance);
//  - return an independently built admissible (symmetric, ASR) set unchanged;
//  - equal the orthogonal projection onto the symmetric ASR subspace from a dense reference;
//  - act per component (9 and 81 components against single-component runs);
//  - hold on a bipartite and on a disconnected support (rank-deficient systems).

#include <Eigen/Dense>
#include <cmath>
#include <complex>
#include <cstdio>
#include <map>
#include <random>
#include <vector>
#include "strain_dfc2_projection.h"

using namespace PHON_NS;

namespace
{
int nfail = 0;

void check(const bool ok, const char *what, const double value)
{
    std::printf("%-60s %10.3e %s\n", what, value, ok ? "ok" : "FAILED");
    if (!ok) ++nfail;
}

using Key = std::array<int, 7>; // i, a, j, b, R
std::map<Key, std::vector<double>> to_map(const std::vector<StrainDfc2Entry> &set)
{
    std::map<Key, std::vector<double>> m;
    for (const auto &e: set) {
        auto &v = m[{e.i, e.a, e.j, e.b, e.R[0], e.R[1], e.R[2]}];
        v.resize(e.val.size(), 0.0);
        for (std::size_t c = 0; c < e.val.size(); ++c) v[c] += e.val[c];
    }
    return m;
}

double max_difference(const std::vector<StrainDfc2Entry> &x, const std::vector<StrainDfc2Entry> &y)
{
    auto mx = to_map(x);
    auto my = to_map(y);
    double d = 0.0;
    for (const auto &[k, v]: mx) {
        const auto it = my.find(k);
        for (std::size_t c = 0; c < v.size(); ++c) {
            d = std::max(d, std::abs(v[c] - (it == my.end() ? 0.0 : it->second[c])));
        }
    }
    for (const auto &[k, v]: my) {
        if (mx.count(k)) continue;
        for (const auto x_: v) d = std::max(d, std::abs(x_));
    }
    return d;
}

// rotation S (signed permutation) of the Cartesian blocks and R, atoms permuted by perm
std::vector<StrainDfc2Entry> transform(const std::vector<StrainDfc2Entry> &set, const int S[3][3],
                                       const std::vector<int> &perm)
{
    std::vector<StrainDfc2Entry> out;
    for (const auto &e: set) {
        for (int a2 = 0; a2 < 3; ++a2) {
            for (int b2 = 0; b2 < 3; ++b2) {
                const double w = S[a2][e.a] * S[b2][e.b];
                if (w == 0.0) continue;
                StrainDfc2Entry t{perm[e.i], a2, perm[e.j], b2, {0, 0, 0}, e.val};
                for (int x = 0; x < 3; ++x) {
                    t.R[x] = 0;
                    for (int y = 0; y < 3; ++y) t.R[x] += S[x][y] * e.R[y];
                }
                for (auto &v: t.val) v *= w;
                out.push_back(t);
            }
        }
    }
    return out;
}

std::vector<StrainDfc2Entry> random_set(std::mt19937 &rng, const int nat, const std::size_t ncomp, const int n,
                                        const std::vector<std::array<int, 2>> &atom_pairs = {})
{
    std::uniform_real_distribution<double> uni(-1.0, 1.0);
    std::uniform_int_distribution<int> iat(0, nat - 1), icart(0, 2), icell(-1, 1), ipair(0, 99);
    std::vector<StrainDfc2Entry> set;
    for (int k = 0; k < n; ++k) {
        int i = iat(rng), j = iat(rng);
        if (!atom_pairs.empty()) {
            const auto &pr = atom_pairs[ipair(rng) % atom_pairs.size()];
            i = pr[0];
            j = pr[1];
        }
        StrainDfc2Entry e{i,
                          icart(rng),
                          j,
                          icart(rng),
                          {icell(rng), icell(rng), icell(rng)},
                          std::vector<double>(ncomp)};
        for (auto &v: e.val) v = uni(rng);
        set.push_back(e);
    }
    return set;
}

// max |row sum| / max |M| and the index-symmetry error of a set
std::array<double, 2> asr_and_symmetry(const std::vector<StrainDfc2Entry> &set, const int nat)
{
    const auto m = to_map(set);
    const auto ncomp = m.begin()->second.size();
    std::vector<double> rows(nat * 9 * ncomp, 0.0);
    double vmax = 0.0, sym = 0.0;
    for (const auto &[k, v]: m) {
        const auto it = m.find({k[2], k[3], k[0], k[1], -k[4], -k[5], -k[6]});
        for (std::size_t c = 0; c < ncomp; ++c) {
            rows[(k[0] * 9 + k[1] * 3 + k[3]) * ncomp + c] += v[c];
            vmax = std::max(vmax, std::abs(v[c]));
            sym = std::max(sym, std::abs(v[c] - (it == m.end() ? 0.0 : it->second[c])));
        }
    }
    double r = 0.0;
    for (const auto x: rows) r = std::max(r, std::abs(x));
    return {r / vmax, sym};
}
} // namespace

int main()
{
    constexpr int nat = 3;
    constexpr std::size_t ncomp = 4;
    std::mt19937 rng(12345);
    std::uniform_real_distribution<double> uni(-1.0, 1.0);
    std::uniform_int_distribution<int> iat(0, nat - 1), icart(0, 2), icell(-1, 1);

    std::vector<StrainDfc2Entry> raw;
    for (int n = 0; n < 200; ++n) {
        StrainDfc2Entry e{iat(rng),
                          icart(rng),
                          iat(rng),
                          icart(rng),
                          {icell(rng), icell(rng), icell(rng)},
                          std::vector<double>(ncomp)};
        for (auto &v: e.val) v = uni(rng);
        raw.push_back(e);
    }

    double dmax = 0.0;
    const auto proj = project_strain_dfc2(raw, nat, ncomp, &dmax);
    const auto m = to_map(proj);

    double sym = 0.0;
    for (const auto &[k, v]: m) {
        const auto it = m.find({k[2], k[3], k[0], k[1], -k[4], -k[5], -k[6]});
        for (std::size_t c = 0; c < ncomp; ++c) {
            sym = std::max(sym, std::abs(v[c] - (it == m.end() ? 0.0 : it->second[c])));
        }
    }
    check(sym < 1e-14, "index symmetry M(ia,jb;R) = M(jb,ia;-R)", sym);

    double asr = 0.0;
    std::vector<double> rows(nat * 9 * ncomp, 0.0);
    for (const auto &[k, v]: m) {
        for (std::size_t c = 0; c < ncomp; ++c) rows[(k[0] * 9 + k[1] * 3 + k[3]) * ncomp + c] += v[c];
    }
    for (const auto r: rows) asr = std::max(asr, std::abs(r));
    check(asr < 1e-12, "ASR sum_{j,R} M(ia,jb;R) = 0 per component", asr);

    // D(q) at a generic q, component 1: Hermitian
    const double q[3] = {0.137, -0.291, 0.418};
    std::vector<std::complex<double>> D(9 * nat * nat, 0.0);
    for (const auto &[k, v]: m) {
        const double ph = 2.0 * M_PI * (q[0] * k[4] + q[1] * k[5] + q[2] * k[6]);
        D[(3 * k[0] + k[1]) * 3 * nat + 3 * k[2] + k[3]] += v[1] * std::polar(1.0, ph);
    }
    double herm = 0.0;
    for (int x = 0; x < 3 * nat; ++x) {
        for (int y = 0; y < 3 * nat; ++y) {
            herm = std::max(herm, std::abs(D[x * 3 * nat + y] - std::conj(D[y * 3 * nat + x])));
        }
    }
    check(herm < 1e-13, "D(q) Hermitian at a generic q", herm);

    double d2 = 1.0;
    const auto again = project_strain_dfc2(proj, nat, ncomp, &d2);
    check(max_difference(proj, again) < 1e-13 && d2 < 1e-13, "idempotent (second correction vanishes)", d2);

    const int S[3][3] = {{0, -1, 0}, {1, 0, 0}, {0, 0, 1}}; // C4 about z
    const std::vector<int> perm = {2, 0, 1};
    const auto cov = max_difference(project_strain_dfc2(transform(raw, S, perm), nat, ncomp), transform(proj, S, perm));
    check(cov < 1e-12, "covariant under rotation + atom permutation", cov);


    // admissible input built independently: symmetric "springs" K on pairs (i, j, R) with the
    // compensating on-site terms; the projection must return it unchanged
    {
        std::vector<StrainDfc2Entry> ok;
        auto add = [&](const int i,
                       const int a,
                       const int j,
                       const int b,
                       const std::array<int, 3>
                           R,
                       const std::vector<double> &v) { ok.push_back({i, a, j, b, R, v}); };
        std::uniform_real_distribution<double> u2(-1.0, 1.0);
        for (int n = 0; n < 12; ++n) {
            const int i = n % nat, j = (n * 7 + 1) % nat;
            const std::array<int, 3> R = {n % 3 - 1, (n / 3) % 3 - 1, 1};
            const std::array<int, 3> mR = {-R[0], -R[1], -R[2]};
            for (int a = 0; a < 3; ++a) {
                for (int b = a; b < 3; ++b) {
                    std::vector<double> v(ncomp), mv(ncomp);
                    for (std::size_t c = 0; c < ncomp; ++c) {
                        v[c] = u2(rng);
                        mv[c] = -v[c];
                    }
                    // K(a, b) = K(b, a): blocks (a, b) and, off the diagonal, (b, a)
                    for (int t = 0; t < (a == b ? 1 : 2); ++t) {
                        const int x = t ? b : a, y = t ? a : b;
                        add(i, x, j, y, R, v);
                        add(j, y, i, x, mR, v);
                        add(i, x, i, y, {0, 0, 0}, mv);
                        add(j, y, j, x, {0, 0, 0}, mv);
                    }
                }
            }
        }
        double d = 1.0;
        const auto back = project_strain_dfc2(ok, nat, ncomp, &d);
        const auto rs = asr_and_symmetry(ok, nat);
        check(rs[0] < 1e-14 && rs[1] < 1e-14, "independent admissible input is admissible", std::max(rs[0], rs[1]));
        check(max_difference(ok, back) < 1e-13, "admissible input returned unchanged", max_difference(ok, back));
    }

    // dense reference: orthogonal projection of the raw set onto {symmetric} cap {ASR}, over the
    // same support, x = raw - pinv(C) C raw
    {
        constexpr int nat_s = 2;
        const auto small = random_set(rng, nat_s, 1, 10);
        const auto proj_s = project_strain_dfc2(small, nat_s, 1);
        std::map<Key, int> var;
        for (const auto &e: proj_s) var.emplace(Key{e.i, e.a, e.j, e.b, e.R[0], e.R[1], e.R[2]}, 0);
        int nv = 0;
        for (auto &kv: var) kv.second = nv++;
        Eigen::VectorXd x0 = Eigen::VectorXd::Zero(nv);
        for (const auto &e: small) x0(var.at({e.i, e.a, e.j, e.b, e.R[0], e.R[1], e.R[2]})) += e.val[0];
        std::vector<Eigen::VectorXd> rows_c;
        for (const auto &[k, id]: var) { // symmetry
            const auto t = var.at({k[2], k[3], k[0], k[1], -k[4], -k[5], -k[6]});
            if (t <= id) continue;
            Eigen::VectorXd c = Eigen::VectorXd::Zero(nv);
            c(id) = 1.0;
            c(t) = -1.0;
            rows_c.push_back(c);
        }
        for (int i = 0; i < nat_s; ++i) { // ASR
            for (int a = 0; a < 3; ++a) {
                for (int b = 0; b < 3; ++b) {
                    Eigen::VectorXd c = Eigen::VectorXd::Zero(nv);
                    for (const auto &[k, id]: var) {
                        if (k[0] == i && k[1] == a && k[3] == b) c(id) = 1.0;
                    }
                    rows_c.push_back(c);
                }
            }
        }
        Eigen::MatrixXd C(rows_c.size(), nv);
        for (std::size_t r = 0; r < rows_c.size(); ++r) C.row(static_cast<Eigen::Index>(r)) = rows_c[r].transpose();
        const Eigen::VectorXd ref = x0 - C.completeOrthogonalDecomposition().solve(C * x0);
        double d = 0.0;
        for (const auto &e: proj_s) {
            d = std::max(d, std::abs(e.val[0] - ref(var.at({e.i, e.a, e.j, e.b, e.R[0], e.R[1], e.R[2]}))));
        }
        check(d < 1e-12, "equals the dense min-norm (KKT) reference", d);
    }

    // per-component action: M_c = alpha_c X + beta_c Y must give alpha_c P(X) + beta_c P(Y)
    for (const std::size_t nc: {std::size_t{9}, std::size_t{81}}) {
        // X and Y share their keys: the projection depends on the support, which must be common
        const auto X = random_set(rng, nat, 1, 60);
        auto Y = X;
        std::uniform_real_distribution<double> u3(-1.0, 1.0);
        for (auto &e: Y) e.val[0] = u3(rng);
        std::vector<StrainDfc2Entry> mix;
        auto coef = [](const std::size_t c, const int w) { return std::cos(0.37 * (c + 1) * w) + 0.1 * w; };
        for (int w = 1; w <= 2; ++w) {
            for (const auto &e: w == 1 ? X : Y) {
                StrainDfc2Entry t = e;
                t.val.assign(nc, 0.0);
                for (std::size_t c = 0; c < nc; ++c) t.val[c] = coef(c, w) * e.val[0];
                mix.push_back(t);
            }
        }
        const auto pm = to_map(project_strain_dfc2(mix, nat, nc));
        const auto px = to_map(project_strain_dfc2(X, nat, 1));
        const auto py = to_map(project_strain_dfc2(Y, nat, 1));
        double d = 0.0;
        for (const auto &[k, v]: pm) {
            const auto ix = px.find(k);
            const auto iy = py.find(k);
            for (std::size_t c = 0; c < nc; ++c) {
                const double ref = coef(c, 1) * (ix == px.end() ? 0.0 : ix->second[0]) +
                                   coef(c, 2) * (iy == py.end() ? 0.0 : iy->second[0]);
                d = std::max(d, std::abs(v[c] - ref));
            }
        }
        char what[80];
        std::snprintf(what, sizeof(what), "%zu components act independently (linear mix)", nc);
        check(d < 1e-12, what, d);
    }

    // rank-deficient systems: bipartite support (pairs 0-1 only; D + N singular) and a
    // disconnected support ({0, 1} and {2, 3})
    {
        double res = 1.0;
        const auto bip = project_strain_dfc2(random_set(rng, 2, 4, 30, {{0, 1}}), 2, 4, nullptr, &res);
        const auto rs = asr_and_symmetry(bip, 2);
        check(res < 1e-12 && rs[0] < 1e-12 && rs[1] < 1e-14, "bipartite support: scaled ASR residual", rs[0]);
        const auto dis = project_strain_dfc2(random_set(rng, 4, 4, 60, {{0, 1}, {1, 1}, {2, 3}, {3, 2}, {2, 2}}),
                                             4,
                                             4,
                                             nullptr,
                                             &res);
        const auto rs2 = asr_and_symmetry(dis, 4);
        check(res < 1e-12 && rs2[0] < 1e-12 && rs2[1] < 1e-14, "disconnected support: scaled ASR residual", rs2[0]);
    }

    std::printf("max |delta| on the random set: %.3e\n", dmax);
    std::printf("%s\n", nfail ? "test_strain_dfc2_projection: FAILED" : "test_strain_dfc2_projection: pass");
    return nfail ? 1 : 0;
}
