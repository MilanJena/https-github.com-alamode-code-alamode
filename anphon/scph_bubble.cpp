/*
 scph_bubble.cpp

 Copyright (c) 2015 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory
 or http://opensource.org/licenses/mit-license.php for information.
*/

/*
 Functions for computing bubble diagram corrections to SCPH.
 The bubble self-energy provides higher-order anharmonic corrections
 on top of the self-consistent phonon (SCPH) calculation.

 Functions included:
 - compute_free_energy_bubble_SCPH: Compute free energy from bubble diagrams
 - bubble_correction: Calculate bubble self-energy corrections to frequencies
*/

#include <algorithm>
#include <complex>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <vector>
#include "anharmonic_core.h"
#include "bubble_qp_root.h"
#include "constants.h"
#include "degeneracy_utils.h"
#include "dynamical.h"
#include "error.h"
#include "fcs_phonon.h"
#include "ifc_derivative.h"
#include "integration.h"
#include "kpoint.h"
#include "memory.h"
#include "mpi_common.h"
#include "phonon_dos.h"
#include "relaxation.h"
#include "scph.h"
#include "selfenergy.h"
#include "system.h"
#include "thermodynamics.h"

using namespace PHON_NS;

void Scph::compute_free_energy_bubble_SCPH(const unsigned int kmesh[3], std::complex<double> ****delta_dymat_scph)
{
    const auto NT = system->get_num_temperature_points();
    const auto nk_ref = dos->kmesh_dos->nk;
    const auto ns = dynamical->neval;
    NDArray<double, 3> eval;
    NDArray<std::complex<double>, 4> evec;

    if (run.my_rank == 0) {
        std::cout << '\n';
        std::cout << " -----------------------------------------------------------------\n";
        std::cout << " Calculating the vibrational free energy from the Bubble diagram \n";
        std::cout << " on top of the SCPH calculation.\n\n";
        std::cout << " This calculation requires allocation of additional memory:\n";

        size_t nsize = nk_ref * ns * ns * NT * sizeof(std::complex<double>) + nk_ref * ns * NT * sizeof(double);

        const auto nsize_dble = static_cast<double>(nsize) / 1000000000.0;
        std::cout << "  Estimated memory usage per MPI process: " << std::setw(10) << std::fixed << std::setprecision(4)
                  << nsize_dble << " GByte.\n";
        std::cout << "  To avoid possible faults associated with insufficient memory,\n"
                     "  please reduce the number of MPI processes per node and/or\n"
                     "  the number of temperature grids.\n\n";
    }

    thermodynamics->FE_bubble.resize(NT);
    eval.resize(NT, nk_ref, ns);
    evec.resize(NT, nk_ref, ns, ns); // This requires lots of RAM

    for (auto iT = 0; iT < NT; ++iT) {
        dynamical->exec_interpolation(kmesh,
                                      delta_dymat_scph[iT],
                                      nk_ref,
                                      dos->kmesh_dos->xk,
                                      dos->kmesh_dos->kvec_na,
                                      eval[iT],
                                      evec[iT],
                                      mindist_list,
                                      fcs_phonon->force_constant_with_cell[0],
                                      *dielec,
                                      *ewald);
    }

    thermodynamics->compute_FE_bubble_SCPH(eval,
                                           evec,
                                           thermodynamics->FE_bubble,
                                           *system,
                                           *dos->kmesh_dos.get(),
                                           symmetry->SymmList,
                                           *anharmonic_core,
                                           dynamical->neval,
                                           run.my_rank,
                                           run.nprocs);

    eval.clear();
    evec.clear();

    if (run.my_rank == 0) {
        std::cout << " done!\n\n";
    }
}

void Scph::bubble_correction(std::complex<double> ****delta_dymat_scph,
                             std::complex<double> ****delta_dymat_scph_plus_bubble)
{
    const auto NT = system->get_num_temperature_points();
    const auto ns = dynamical->neval;

    auto epsilon = integration->epsilon;
    const auto nk_irred_interpolate = kmesh_coarse->nk_irred;

    // The bubble self-energy is summed on KMESH_BUBBLE; the SCPH mesh and its
    // phase cache are reused when the two meshes coincide.
    MPI_Bcast(&kmesh_bubble[0], 3, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
    std::unique_ptr<KpointMeshUniform> kmesh_bubble_own;
    std::unique_ptr<PhaseFactorCache> phase_factor_own;
    std::vector<int> kmap_coarse_to_bubble_own;
    const KpointMeshUniform *kmesh_b = kmesh_dense.get();
    const PhaseFactorCache *phase_b = phase_factor.get();
    const std::vector<int> *kmap_b = &kmap_coarse_to_dense;
    if (kmesh_bubble[0] != kmesh_dense->nk_i[0] || kmesh_bubble[1] != kmesh_dense->nk_i[1] ||
        kmesh_bubble[2] != kmesh_dense->nk_i[2])
    {
        kmesh_bubble_own = std::make_unique<KpointMeshUniform>(kmesh_bubble);
        kmesh_bubble_own->setup(symmetry->SymmList,
                                system->get_reference_primcell().reciprocal_lattice_vector, // as kmesh_dense
                                symmetry->use_time_reversal && symmetry->time_reversal_sym);
        if (kpoint->get_kmap_coarse_to_dense(kmesh_coarse.get(), kmesh_bubble_own.get(), kmap_coarse_to_bubble_own) ==
            1)
        {
            exit("bubble_correction", "KMESH_BUBBLE should be an integral multiple of KMESH_INTERPOLATE");
        }
        phase_factor_own = std::make_unique<PhaseFactorCache>(kmesh_bubble_own->nk_i);
        phase_factor_own->create(true);
        kmesh_b = kmesh_bubble_own.get();
        phase_b = phase_factor_own.get();
        kmap_b = &kmap_coarse_to_bubble_own;
    }
    const auto nk_scph = kmesh_b->nk;

    NDArray<double, 2> eval;
    NDArray<double, 3> eval_bubble;
    NDArray<std::complex<double>, 3> evec;
    NDArray<double, 1> real_self;
    std::vector<std::complex<double>> omegalist;

    if (run.my_rank == 0) {
        std::cout << '\n';
        std::cout << " -----------------------------------------------------------------\n";
        std::cout << " Calculating the bubble self-energy \n";
        std::cout << " on top of the SCPH calculation.\n";
        std::cout << "  KMESH_BUBBLE: " << kmesh_b->nk_i[0] << " " << kmesh_b->nk_i[1] << " " << kmesh_b->nk_i[2]
                  << " (" << kmesh_b->nk << " points)\n\n";
    }

    eval.resize(nk_scph, ns);
    evec.resize(nk_scph, ns, ns);
    selfenergy->setup_selfenergy(dynamical->neval,
                                 integration->epsilon,
                                 thermodynamics->classical,
                                 run.my_rank,
                                 run.nprocs);

    if (run.my_rank == 0) {
        eval_bubble.resize(NT, nk_scph, ns);
        for (auto iT = 0; iT < NT; ++iT) {
            for (auto ik = 0; ik < nk_scph; ++ik) {
                for (auto is = 0; is < ns; ++is) {
                    eval_bubble[iT][ik][is] = 0.0;
                }
            }
        }
        real_self.resize(ns);
    }

    NDArray<std::vector<int>, 1> degeneracy_at_k;
    degeneracy_at_k.resize(nk_scph);

    // Relaxed run (RELAX_STR != 0). The SCP propagators already belong to the
    // relaxed structure of each temperature: the stored correction is the full
    // SCP matrix minus the reference harmonic one, so it contains the q0 and
    // strain renormalization, and the phases exp(ik.R) of fractional k do not
    // change under the homogeneous strain. Only the cubic IFCs do, to
    // Phi3 + Phi4 : d(T). Temperatures whose SCP iteration or structure did not
    // converge keep the SCPH result without the bubble correction.
    // The recorded structure is the one the optimizer accepted, one update past
    // the structure of the last SCP solve (see ScphRelaxationModel::
    // finalize_temperature); the two differ within COORD_CONV_TOL/CELL_CONV_TOL.
    const auto relaxed = relaxation->relax_str != 0;
    std::vector<double> u_tensor_all, u0_all;
    std::vector<int> use_temp(NT, 1);
    if (relaxed) {
        const auto nt = static_cast<size_t>(NT);
        const auto nsz = static_cast<size_t>(ns);
        if (nt > static_cast<size_t>(std::numeric_limits<int>::max()) / std::max(nsz, size_t{9})) {
            exit("bubble_correction", "Too many temperatures x modes to broadcast the relaxed structures.");
        }
        if (run.my_rank == 0) {
            if (relaxed_structure.u_tensor.size() != nt * 9 || relaxed_structure.u0.size() != nt * nsz ||
                relaxed_structure_recorded.size() != nt)
            {
                exit("bubble_correction", "The relaxed structures of this run are not available.");
            }
            u_tensor_all = relaxed_structure.u_tensor;
            u0_all = relaxed_structure.u0;
            for (size_t iT = 0; iT < nt; ++iT) {
                // An unrecorded temperature was never visited by the loop, so its
                // SCP matrix does not exist either: nothing to keep or correct.
                if (!relaxed_structure_recorded[iT]) {
                    exit("bubble_correction",
                         "A temperature of the grid has no relaxed structure; the structural\n"
                         " optimization loop skipped it. Check TMIN, TMAX and DT.");
                }
                const auto conv_scph = converged_scph_temp.size() != nt || converged_scph_temp[iT];
                const auto conv_str = converged_str_temp.size() != nt || converged_str_temp[iT];
                use_temp[iT] = conv_scph && conv_str ? 1 : 0;
            }
        } else {
            u_tensor_all.resize(nt * 9);
            u0_all.resize(nt * nsz);
        }
        MPI_Bcast(u_tensor_all.data(), static_cast<int>(nt * 9), MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(u0_all.data(), static_cast<int>(nt * nsz), MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(use_temp.data(), static_cast<int>(NT), MPI_INT, 0, MPI_COMM_WORLD);
        if (run.my_rank == 0) {
            std::cout << "  RELAX_STR != 0: the cubic IFCs are deformed to the relaxed structure\n"
                         "  of each temperature (Phi3 + Phi4 : d).\n\n";
        }
    }

    for (auto iT = 0; iT < NT; ++iT) {
        const auto temp = system->Tmin + system->dT * float(iT);

        if (relaxed) {
            if (!use_temp[iT]) {
                if (run.my_rank == 0) {
                    std::cout << " Temperature (K) : " << std::setw(6) << temp
                              << " : skipped (SCPH or structure not converged); the SCPH result is kept.\n\n";
                    for (auto is = 0; is < ns; ++is) {
                        for (auto js = 0; js < ns; ++js) {
                            for (auto ik = 0; ik < kmesh_coarse->nk; ++ik) {
                                delta_dymat_scph_plus_bubble[iT][is][js][ik] = delta_dymat_scph[iT][is][js][ik];
                            }
                        }
                    }
                }
                continue;
            }
            const auto offset = static_cast<size_t>(iT) * ns;
            const std::vector<double> u0_T(u0_all.begin() + offset, u0_all.begin() + offset + ns);
            std::vector<FcsArrayWithCell> fc3_deformed;
            DerivativeIFC::compute_deformed_cubic_ifcs(fcs_phonon->force_constant_with_cell,
                                                       &u_tensor_all[static_cast<size_t>(iT) * 9],
                                                       u0_T,
                                                       system->get_primcell().lattice_vector,
                                                       fc3_deformed);
            anharmonic_core->replace_cubic(fc3_deformed);
        }

        dynamical->exec_interpolation(kmesh_interpolate,
                                      delta_dymat_scph[iT],
                                      nk_scph,
                                      kmesh_b->xk,
                                      kmesh_b->kvec_na,
                                      eval,
                                      evec,
                                      mindist_list,
                                      fcs_phonon->force_constant_with_cell[0],
                                      *dielec,
                                      *ewald);

        for (unsigned int ik = 0; ik < nk_scph; ++ik) find_degenerate_groups(ns, eval[ik], degeneracy_at_k[ik]);

        if (run.my_rank == 0) std::cout << " Temperature (K) : " << std::setw(6) << temp << '\n';

        for (auto ik = 0; ik < nk_irred_interpolate; ++ik) {

            auto knum_interpolate = kmesh_coarse->kpoint_irred_all[ik][0].knum;
            auto knum = (*kmap_b)[knum_interpolate];

            if (run.my_rank == 0) {
                std::cout << "  Irred. k: " << std::setw(5) << ik + 1 << " (";
                for (auto m = 0; m < 3; ++m) std::cout << std::setw(15) << kmesh_b->xk[knum][m];
                std::cout << ")\n";
            }

            for (unsigned int snum = 0; snum < ns; ++snum) {

                if (eval[knum][snum] < eps8) {
                    if (run.my_rank == 0) real_self[snum] = 0.0;
                } else {
                    omegalist.clear();

                    if (bubble == 1) {

                        omegalist.push_back(im * epsilon);

                        auto se_bubble = selfenergy->get_bubble_selfenergy(kmesh_b,
                                                                           ns,
                                                                           eval,
                                                                           evec,
                                                                           knum,
                                                                           snum,
                                                                           temp,
                                                                           omegalist,
                                                                           phase_b,
                                                                           *anharmonic_core);

                        if (run.my_rank == 0) real_self[snum] = se_bubble[0].real();

                    } else if (bubble == 2) {

                        omegalist.push_back(eval[knum][snum] + im * epsilon);

                        auto se_bubble = selfenergy->get_bubble_selfenergy(kmesh_b,
                                                                           ns,
                                                                           eval,
                                                                           evec,
                                                                           knum,
                                                                           snum,
                                                                           temp,
                                                                           omegalist,
                                                                           phase_b,
                                                                           *anharmonic_core);

                        if (run.my_rank == 0) real_self[snum] = se_bubble[0].real();

                    } else if (bubble == 3) {
                        // Quasiparticle root of w^2 = W^2 - 2 W Re Sigma(w), W = SCPH frequency; see
                        // bubble_qp_root.h. Range [0, W + max(50 cm^-1, W/2)]: the bubble can soften a
                        // mode far below W (soft polar modes) and also stiffen it. Coarse step: at most
                        // EPSILON/4, the scale on which the Lorentzian-broadened Re Sigma can vary, and at
                        // most max(0.5 cm^-1, range/400); refinement step 0.1 cm^-1.
                        const auto omega0 = eval[knum][snum];
                        const auto cm1 = time_ry / Hz_to_kayser;
                        const auto wmax = omega0 + std::max(50.0 * cm1, 0.5 * omega0);
                        const auto dcoarse =
                            std::max(0.1 * cm1, std::min(std::max(0.5 * cm1, wmax / 400.0), 0.25 * epsilon));

                        const auto re_sigma = [&](const std::vector<double> &w) {
                            omegalist.clear();
                            for (const auto x: w) omegalist.push_back(x + im * epsilon);
                            const auto se = selfenergy->get_bubble_selfenergy(kmesh_b,
                                                                              ns,
                                                                              eval,
                                                                              evec,
                                                                              knum,
                                                                              snum,
                                                                              temp,
                                                                              omegalist,
                                                                              phase_b,
                                                                              *anharmonic_core);
                            // reduced on rank 0; every rank takes the same decisions in solve()
                            std::vector<double> re(w.size());
                            for (size_t i = 0; i < w.size(); ++i) re[i] = se[i].real();
                            MPI_Bcast(re.data(), static_cast<int>(re.size()), MPI_DOUBLE, 0, MPI_COMM_WORLD);
                            return re;
                        };
                        const auto qp = BubbleQP::solve(omega0, wmax, dcoarse, 0.1 * cm1, re_sigma);

                        if (run.my_rank == 0) {
                            real_self[snum] = qp.self;
                            std::ostringstream msg;
                            msg << std::fixed << std::setprecision(2);
                            if (qp.best < 0) {
                                msg << "BUBBLE = 3: no quasiparticle root of w^2 = W^2 - 2W Re[Self](w) in [0, "
                                    << in_kayser(wmax) << "] cm^-1 for branch " << snum + 1
                                    << " (W = " << in_kayser(omega0)
                                    << " cm^-1).\n The on-shell value Re[Self](W) is used (as BUBBLE = 2).";
                                warn("bubble_correction", msg.str().c_str());
                            } else if (qp.roots.size() > 1) {
                                msg << "BUBBLE = 3: " << qp.roots.size() << " quasiparticle roots for branch "
                                    << snum + 1 << " (W = " << in_kayser(omega0)
                                    << " cm^-1); the one with the largest Z = 2w / F'(w) is used.\n"
                                       "   root (cm^-1), Z :";
                                for (size_t r = 0; r < qp.roots.size(); ++r) {
                                    msg << "  " << in_kayser(qp.roots[r].omega) << ", " << qp.roots[r].weight
                                        << (static_cast<int>(r) == qp.best ? " (used)" : "");
                                }
                                warn("bubble_correction", msg.str().c_str());
                            }
                        }
                    }
                }
                if (run.my_rank == 0) {
                    // Explicit format: the stream state left by earlier output
                    // otherwise decides the precision (two digits in some runs).
                    const auto flags = std::cout.flags();
                    const auto prec = std::cout.precision();
                    std::cout << std::scientific << std::setprecision(6);
                    std::cout << "   branch : " << std::setw(5) << snum + 1;
                    std::cout << " omega (SC1) = " << std::setw(15) << in_kayser(eval[knum][snum]) << " (cm^-1); ";
                    std::cout << " Re[Self] = " << std::setw(15) << in_kayser(real_self[snum]) << " (cm^-1)\n";
                    std::cout.flags(flags);
                    std::cout.precision(prec);
                }
            }

            if (run.my_rank == 0) {
                // average self energy of degenerate modes
                int ishift = 0;
                double real_self_avg = 0.0;

                for (const auto &it: degeneracy_at_k[knum]) {
                    for (auto m = 0; m < it; ++m) {
                        real_self_avg += real_self[m + ishift];
                    }
                    real_self_avg /= static_cast<double>(it);

                    for (auto m = 0; m < it; ++m) {
                        real_self[m + ishift] = real_self_avg;
                    }
                    real_self_avg = 0.0;
                    ishift += it;
                }

                for (unsigned int snum = 0; snum < ns; ++snum) {
                    eval_bubble[iT][knum][snum] =
                        eval[knum][snum] * eval[knum][snum] - 2.0 * eval[knum][snum] * real_self[snum];
                    for (auto jk = 1; jk < kmesh_coarse->kpoint_irred_all[ik].size(); ++jk) {
                        auto knum2 = (*kmap_b)[kmesh_coarse->kpoint_irred_all[ik][jk].knum];
                        eval_bubble[iT][knum2][snum] = eval_bubble[iT][knum][snum];
                    }
                }

                std::cout << '\n';
            }
        }

        if (run.my_rank == 0) {
            dynamical->calc_new_dymat_with_evec(delta_dymat_scph_plus_bubble[iT],
                                                eval_bubble[iT],
                                                evec,
                                                kmesh_coarse.get(),
                                                *kmap_b,
                                                fcs_phonon->force_constant_with_cell[0],
                                                *dielec,
                                                *ewald);
        }
    }

    // Later consumers expect the reference cubic IFCs.
    if (relaxed) anharmonic_core->replace_cubic(fcs_phonon->force_constant_with_cell[1]);

    eval.clear();
    evec.clear();
    degeneracy_at_k.clear();

    eval_bubble.clear();

    if (run.my_rank == 0) {
        std::cout << " done!\n\n";
    }
}
