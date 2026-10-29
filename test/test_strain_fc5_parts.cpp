/*
 test_strain_fc5_parts.cpp

 Unit test for DerivativeIFC::strain_fc5_parts, the FC5 parts that
 RELAXED_STRUCTURE = 1 adds to the IFCs of a STRAIN_FC5 = 1 state:
 dPhi4 = sum_mn u_mn dPhi4/du_mn (FC4) and dPhi4 : u0 (FC3). They must be the
 terms the SCP relaxation itself used, i.e. what RealSpaceV4::q0_sweep over the
 multi-channel dPhi4/du mesh (ScphQhaCommon::update_fc5_at_structure) adds to
 v3_renorm and q4_q0 at the same u and q0.

 A synthetic, permutation-symmetric FC5 with random values (so with a nonzero
 acoustic-sum-rule residual, which makes the strain moment depend on the atom it
 is measured from) on a two-atom orthorhombic cell; groups built by the SCPH
 builder (compute_dV_dumn_all_real_space); identity eigenvectors on a 2x2x2 mesh,
 so the q0 sweep is the mass-weighted Cartesian vertex and q0 = sqrt(m) u0.

 Checks:
 - v3 (the FC3 vertex at (Gamma, k, -k)) from dfc3 equals the sweep's v3_renorm;
 - a RealSpaceV4 of the follow-up FC4 list dfc4 reproduces the sweep's v3_renorm
   and q4_q0 (the FC4 double contraction that enters FC2);
 - both for all channels and for DIAG (u_xx, u_yy, u_zz only);
 - control: contracting the last leg instead of the home leg differs by about
   strain x u0 x the ASR residual, i.e. the test sees the difference.

 Built by the anphon CMake project as `test_strain_fc5_parts`. Exits 0 on success.
*/

#include <cmath>
#include <complex>
#include <cstdio>
#include <map>
#include <random>
#include <string>
#include <vector>
#include "ifc_derivative.h"
#include "ndarray.h"
#include "quartic_real_space.h"

using namespace PHON_NS;
using cplx = std::complex<double>;

namespace
{
int nfail = 0;

void check(const bool ok, const std::string &what)
{
    std::printf("%s: %s\n", ok ? "ok  " : "FAIL", what.c_str());
    if (!ok) ++nfail;
}

constexpr unsigned int nat = 2, ns = 3 * nat, nk = 8;
const unsigned int nki[3] = {2, 2, 2};

struct Leg
{
    unsigned int index; // 3 * atom + xyz
    int R[3];           // cell
};

// The FC5 list: every ordered 5-tuple of a symmetric tensor, the first leg in the home cell.
std::vector<FcsArrayWithCell> synthetic_fc5(const Eigen::Vector3d tau[nat], std::mt19937 &rng)
{
    std::uniform_int_distribution<int> pick_site(0, 5), pick_xyz(0, 2);
    std::uniform_real_distribution<double> val(-1.0, 1.0);
    const Leg sites[6] = {{0, {0, 0, 0}},
                          {1, {0, 0, 0}},
                          {0, {1, 0, 0}},
                          {1, {0, -1, 0}},
                          {1, {0, 0, 1}},
                          {0, {-1, 1, 0}}}; // index = atom here
    std::map<std::vector<int>, double> entries;
    for (int c = 0; c < 60; ++c) {
        std::vector<Leg> legs(5);
        for (auto &l: legs) {
            l = sites[pick_site(rng)];
            l.index = 3 * l.index + pick_xyz(rng);
        }
        const auto v = val(rng);
        std::vector<int> perm = {0, 1, 2, 3, 4};
        do {
            const auto &h = legs[perm[0]];
            std::vector<int> key;
            for (const auto p: perm) {
                key.push_back(static_cast<int>(legs[p].index));
                for (int d = 0; d < 3; ++d) key.push_back(legs[p].R[d] - h.R[d]);
            }
            entries[key] = v; // set, not add: one value per element of the symmetric tensor
        } while (std::next_permutation(perm.begin(), perm.end()));
    }
    std::vector<FcsArrayWithCell> fc5;
    for (const auto &[key, v]: entries) {
        std::vector<AtomCellSuper> pairs(5);
        std::vector<Eigen::Vector3d> relvecs, relvecs_vel;
        for (int n = 0; n < 5; ++n) {
            pairs[n].index = key[4 * n];
            pairs[n].tran = 0;
            pairs[n].cell_s = 0;
            if (n == 0) continue;
            const Eigen::Vector3d R(key[4 * n + 1], key[4 * n + 2], key[4 * n + 3]);
            relvecs.push_back(R);
            relvecs_vel.push_back(R + tau[key[4 * n] / 3] - tau[key[0] / 3]);
        }
        fc5.emplace_back(v, pairs, std::vector<unsigned int>(5, 0), relvecs, relvecs_vel);
    }
    return fc5;
}

std::string sci(const double x)
{
    char buf[32];
    std::snprintf(buf, sizeof(buf), "%.2e", x);
    return buf;
}

// max |a - b| / max |b| over two nk x ns x ns^2 arrays
double reldiff(const std::vector<cplx> &a, const std::vector<cplx> &b)
{
    double d = 0.0, s = 0.0;
    for (size_t i = 0; i < a.size(); ++i) {
        d = std::max(d, std::abs(a[i] - b[i]));
        s = std::max(s, std::abs(b[i]));
    }
    return d / std::max(s, 1.0e-300);
}

// The (Gamma, k, -k) vertex of a 3-leg list as RealSpaceV4::q0_sweep forms it (identity
// eigenvectors): v3[k][j][l, m] = 1/(4 nk) sum W(j; l, m) e^{-2 pi i k.(R_m - R_l)}
std::vector<cplx> vertex3(const std::vector<FcsArrayWithCell> &fc3, const std::vector<double> &invsqrt_mass,
                          const std::vector<Eigen::Vector3d> &xk)
{
    std::vector<cplx> v(nk * ns * ns * ns, 0.0);
    for (const auto &e: fc3) {
        const auto j = e.pairs[0].index, l = e.pairs[1].index, m = e.pairs[2].index;
        const auto w = e.fcs_val * invsqrt_mass[j / 3] * invsqrt_mass[l / 3] * invsqrt_mass[m / 3];
        const Eigen::Vector3d Rlm = e.relvecs[1] - e.relvecs[0];
        for (unsigned int k = 0; k < nk; ++k) {
            v[((k * ns + j) * ns + l) * ns + m] += w * std::exp(cplx(0.0, -2.0 * M_PI * xk[k].dot(Rlm))) / (4.0 * nk);
        }
    }
    return v;
}

struct Sweep
{
    std::vector<cplx> v3, q4;
};

Sweep sweep(quartic_rs::RealSpaceV4 &rs, const std::vector<double> &q0)
{
    NDArray<cplx, 3> v3(nk, ns, ns * ns), q4(nk, ns, ns);
    for (unsigned int k = 0; k < nk; ++k)
        for (unsigned int a = 0; a < ns; ++a) {
            for (unsigned int c = 0; c < ns * ns; ++c) v3[k][a][c] = 0.0;
            for (unsigned int b = 0; b < ns; ++b) q4[k][a][b] = 0.0;
        }
    rs.q0_sweep(q0.data(), nullptr, v3, q4, true);
    Sweep s;
    for (unsigned int k = 0; k < nk; ++k)
        for (unsigned int a = 0; a < ns; ++a) {
            for (unsigned int c = 0; c < ns * ns; ++c) s.v3.push_back(v3[k][a][c]);
            for (unsigned int b = 0; b < ns; ++b) s.q4.push_back(q4[k][a][b]);
        }
    return s;
}
} // namespace

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    std::mt19937 rng(20261012);
    const Eigen::Vector3d tau[nat] = {{0.0, 0.0, 0.0}, {0.5, 0.4, 0.3}};
    Eigen::Matrix3d lavec = Eigen::Vector3d(4.0, 4.5, 5.0).asDiagonal();
    const std::vector<double> invsqrt_mass = {1.0 / std::sqrt(1.3), 1.0 / std::sqrt(2.7)};

    auto fc5 = synthetic_fc5(tau, rng);
    // ASR residual over the last leg: sum_m Phi5(..., m mu) for fixed other legs
    std::map<std::vector<long>, double> asr;
    double vmax = 0.0;
    for (const auto &e: fc5) {
        std::vector<long> key;
        for (int n = 0; n < 4; ++n) {
            key.push_back(e.pairs[n].index);
            if (n > 0)
                for (int d = 0; d < 3; ++d) key.push_back(std::lround(e.relvecs[n - 1][d]));
        }
        key.push_back(e.pairs[4].index % 3);
        asr[key] += e.fcs_val;
        vmax = std::max(vmax, std::abs(e.fcs_val));
    }
    double res = 0.0;
    for (const auto &[k, v]: asr) res = std::max(res, std::abs(v));
    std::printf("synthetic FC5: %zu entries, ASR residual / max |Phi5| = %.3f\n", fc5.size(), res / vmax);
    check(res > 0.1 * vmax, "the synthetic FC5 violates the acoustic sum rule");

    // the SCPH builder: ScphQhaCommon::build_fc5_correction via compute_dPhi4_dumn_groups
    std::stable_sort(fc5.begin(), fc5.end(), sort_by_heading_indices(1));
    std::vector<DeltaFcsStrainComponents> groups;
    DerivativeIFC::compute_dV_dumn_all_real_space(fc5, groups, 1, lavec);

    std::vector<Eigen::Vector3d> xk;
    for (int a = 0; a < 2; ++a)
        for (int b = 0; b < 2; ++b)
            for (int c = 0; c < 2; ++c) xk.emplace_back(0.5 * a, 0.5 * b, 0.5 * c);
    std::vector<std::vector<double>> xk_rows(nk);
    std::vector<const double *> xk_ptr(nk);
    for (unsigned int k = 0; k < nk; ++k) {
        xk_rows[k] = {xk[k][0], xk[k][1], xk[k][2]};
        xk_ptr[k] = xk_rows[k].data();
    }
    std::vector<std::vector<cplx>> eye(ns, std::vector<cplx>(ns, 0.0));
    for (unsigned int s = 0; s < ns; ++s) eye[s][s] = 1.0;
    std::vector<const cplx *> eye_rows(ns);
    for (unsigned int s = 0; s < ns; ++s) eye_rows[s] = eye[s].data();
    const std::vector<const cplx *const *> evec(nk, eye_rows.data());
    const std::vector<bool> acoustic(ns, false);
    std::vector<unsigned int> kirr(nk);
    for (unsigned int k = 0; k < nk; ++k) kirr[k] = k;

    std::uniform_real_distribution<double> val(-1.0, 1.0);
    double u[9];
    for (auto &x: u) x = 0.01 * val(rng);
    std::vector<double> q0(ns), u0(ns);
    for (unsigned int a = 0; a < ns; ++a) {
        q0[a] = 0.1 * val(rng);
        u0[a] = q0[a] * invsqrt_mass[a / 3];
    }

    for (const auto diag: {false, true}) {
        const std::string tag = diag ? " (DIAG)" : " (ALL)";
        double uw[9];
        for (int mn = 0; mn < 9; ++mn) uw[mn] = diag && mn / 3 != mn % 3 ? 0.0 : u[mn];

        // the relaxation: the q0 sweep over the multi-channel dPhi4/du mesh
        quartic_rs::QuarticMesh mesh(groups, 9, invsqrt_mass, ns, nki);
        quartic_rs::RealSpaceV4
            scph(std::move(mesh), ns, nki, xk_ptr.data(), evec.data(), acoustic, 0, kirr, true, false);
        scph.set_channel_weights(uw);
        const auto ref = sweep(scph, q0);

        // the follow-up
        std::vector<FcsArrayWithCell> dfc4, dfc3;
        DerivativeIFC::strain_fc5_parts(groups, u, diag, u0, dfc4, dfc3);
        const auto d3 = reldiff(vertex3(dfc3, invsqrt_mass, xk), ref.v3);
        check(d3 < 1e-12, "FC3: dPhi4 : u0 against the relaxation's q0 sweep" + tag + ", rel. dev. " + sci(d3));
        quartic_rs::RealSpaceV4
            follow(dfc4, invsqrt_mass, ns, nki, xk_ptr.data(), evec.data(), acoustic, 0, kirr, true);
        const auto f = sweep(follow, q0);
        const auto e3 = reldiff(f.v3, ref.v3), e4 = reldiff(f.q4, ref.q4);
        check(e3 < 1e-12, "FC4: dPhi4 in the cubic q0 sweep" + tag + ", rel. dev. " + sci(e3));
        check(e4 < 1e-12, "FC4: dPhi4 : q0 q0 (enters FC2)" + tag + ", rel. dev. " + sci(e4));

        // control: the last leg contracted instead (the strain moment then measured from another atom)
        auto aligned = dfc4;
        std::stable_sort(aligned.begin(), aligned.end(), sort_by_heading_indices(1));
        std::vector<FcsArrayWithCell> last;
        DerivativeIFC::compute_dV_dsublattice_real_space(aligned,
                                                         last,
                                                         Eigen::Map<const Eigen::VectorXd>(u0.data(), ns),
                                                         -1.0);
        const auto dl = reldiff(vertex3(last, invsqrt_mass, xk), ref.v3);
        check(dl > 1e-3, "control: the last-leg contraction differs" + tag + ", rel. dev. " + sci(dl));
    }

    std::printf("STRAIN_FC5 follow-up parts --> %s\n", nfail ? "fail" : "pass");
    MPI_Finalize();
    return nfail ? 1 : 0;
}
