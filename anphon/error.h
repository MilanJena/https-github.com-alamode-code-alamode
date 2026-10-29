/*
 error.h

 Copyright (c) 2014 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory 
 or http://opensource.org/licenses/mit-license.php for information.
*/

#pragma once

#ifdef _WIN32
#include <mpi.h>
#else

#include "mpi.h"

#endif

#include <cstdlib>
#include <iostream>
#include <string>

#include "fatal_error.h"

namespace PHON_NS
{
inline void warn(const char *file, const char *message)
{
    std::cout << alamode::log_message("WARNING", file, message);
}


// Terminate the whole MPI job, not only the calling process. A failure that
// is detected on one rank only (a file that cannot be opened, an exception
// while reading it, ...) would otherwise leave the other ranks waiting in the
// next collective forever. MPI_Abort is skipped when MPI is not (or no longer)
// initialized so that the helpers stay usable in serial tools.
[[noreturn]] inline void abort_all_ranks()
{
    std::cout.flush();
    int initialized = 0, finalized = 0;
    MPI_Initialized(&initialized);
    MPI_Finalized(&finalized);
    if (initialized && !finalized) MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    std::exit(EXIT_FAILURE);
}

// Print the message, then throw alamode::Error when the module enabled the
// throw policy (alamode::throw_on_error, never set by the CLI), the run is a
// single process (MPI not running, or one rank) and not inside an OpenMP
// parallel region; otherwise abort all ranks as before. Throwing on one rank
// of a multi-rank job would leave the others in a collective.
[[noreturn]] inline void fail(const std::string &log, const std::string &text)
{
    std::cout << log;
    if (alamode::can_throw()) {
        int initialized = 0, finalized = 0, nprocs = 1;
        MPI_Initialized(&initialized);
        MPI_Finalized(&finalized);
        if (initialized && !finalized) MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
        if (nprocs == 1) alamode::terminate(text);
    }
    abort_all_ranks();
}

// The log gets the wrapped message; the exception text stays on one line.
[[noreturn]] inline void exit(const char *file, const char *message)
{
    fail(alamode::log_message("ERROR", file, message),
         alamode::format_message("ERROR in ", file, "  MESSAGE: ", message));
}

template <typename T>
[[noreturn]] void exit(const char *file, const char *message, const T info)
{
    const auto text = alamode::format_message(message, info);
    fail(alamode::log_message("ERROR", file, text), alamode::format_message("ERROR in ", file, "  MESSAGE: ", text));
}
} // namespace PHON_NS
