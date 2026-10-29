/*
 dense_symmetric_eigen.cpp

 Copyright (c) 2026 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory
 or http://opensource.org/licenses/mit-license.php for information.
*/

#include "dense_symmetric_eigen.h"
#include <mpi.h>
#include <string>
#include "error.h"

extern "C"
{
    void dsyevd_(const char *jobz, const char *uplo, int *n, double *a, int *lda, double *w, double *work, int *lwork,
                 int *iwork, int *liwork, int *info);
}

using namespace PHON_NS;

void PHON_NS::solve_dense_symmetric(int n, std::vector<double> &A, std::vector<double> &w, int num_lowest)
{
    // LAPACK backend: full spectrum regardless of num_lowest (see header).
    (void)num_lowest;

    if (n <= 0) return;
    w.resize(n);

    // Divide and conquer: much faster than dsyev for the full set of
    // eigenvectors; threaded by the LAPACK library (e.g. MKL_NUM_THREADS).
    int info = 0;
    int lwork = -1, liwork = -1;
    double work_query = 0.0;
    int iwork_query = 0;
    dsyevd_("V", "L", &n, A.data(), &n, w.data(), &work_query, &lwork, &iwork_query, &liwork, &info);
    if (info != 0) {
        exit("solve_dense_symmetric", "dsyevd workspace query failed");
    }
    lwork = static_cast<int>(work_query);
    liwork = iwork_query;
    std::vector<double> work(lwork);
    std::vector<int> iwork(liwork);
    dsyevd_("V", "L", &n, A.data(), &n, w.data(), work.data(), &lwork, iwork.data(), &liwork, &info);
    if (info != 0) {
        exit("solve_dense_symmetric", ("dsyevd failed with info = " + std::to_string(info)).c_str());
    }
}
