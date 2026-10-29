// _anphon.cpp -- nanobind binding that runs anphon on an input file.
//
// v1 surface (RELEASE_2.0_PLAN.md, P2): run_input_file(path, workdir=None), the same run as the CLI
// `anphon path` through PhononCUI. The Python layer (python/anphon/) builds on it.
//
// MPI: importing the module initializes MPI unless it is already initialized (e.g. by mpi4py, which then
// owns it); MPI is finalized at interpreter exit only if this module initialized it. Import and every run
// are collective over MPI_COMM_WORLD; a communicator argument is not offered in v1 because anphon itself
// hardcodes MPI_COMM_WORLD (phonon_cui.cpp and ~200 other sites). mpi4py, if used, must be built against
// the same MPI library.
//
// Threads: v1 supports import and run_input_file on the Python main thread only (MPI_THREAD_FUNNELED, and
// the atexit finalisation then never overlaps a run). The GIL is released during a run, so other Python
// threads keep running, but they cannot start a second run: runs are therefore also serialised, which
// anphon needs (fixed output names, process-wide chdir).
//
// Build: see python/CMakeLists.txt (nanobind_add_module + anphoncxx_static).

#include <nanobind/nanobind.h>
#include <nanobind/stl/filesystem.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/string.h>

#include <exception>
#include <filesystem>
#include <ios>
#include <iostream>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

#include "error.h"
#include "fatal_error.h"
#include "ndarray.h"
#include "phonon_cui.h"

namespace nb = nanobind;
namespace fs = std::filesystem;

namespace
{
bool mpi_initialized_here = false;

// ponytail: main thread only (see the top of the file); per-thread runs need MPI_THREAD_SERIALIZED and a
// run lock that also covers finalisation.
bool on_python_main_thread()
{
    const auto threading = nb::module_::import_("threading");
    return threading.attr("main_thread")().attr("ident").equal(threading.attr("get_ident")());
}

// Process-wide state of one run, restored on every exit path:
// - the error policy (include/fatal_error.h, include/ndarray.h): throw (-> RuntimeError) only in a
//   single-process run; otherwise, and for out-of-memory, abort all ranks;
// - the std::cout format flags: a run that fails halfway leaves its precision/width/floatfield flags
//   behind, which would garble the log of the next run in the same process.
struct RunStateScope
{
    std::ios saved_format{nullptr};
    bool saved_throw = alamode::throw_on_error;
    void (*saved_abort)() = alamode::abort_hook;
    void (*saved_oom)() = ndarray_detail::ndarray_out_of_memory_hook();

    explicit RunStateScope(const bool serial)
    {
        saved_format.copyfmt(std::cout);
        alamode::throw_on_error = serial;
        alamode::abort_hook = PHON_NS::abort_all_ranks;
        ndarray_detail::ndarray_out_of_memory_hook() = PHON_NS::abort_all_ranks;
    }

    ~RunStateScope()
    {
        std::cout.copyfmt(saved_format);
        alamode::throw_on_error = saved_throw;
        alamode::abort_hook = saved_abort;
        ndarray_detail::ndarray_out_of_memory_hook() = saved_oom;
        // anphon writes to the C++ std::cout; Python's sys.stdout has its own buffer, so text printed from
        // Python around a run may interleave out of order when stdout is not a terminal. No redirection in v1.
        std::cout.flush();
    }
};

// chdir to workdir for the run; the original cwd is restored on every exit path.
struct CwdScope
{
    fs::path saved = fs::current_path();

    explicit CwdScope(const std::optional<fs::path> &dir)
    {
        if (dir) fs::current_path(*dir);
    }

    ~CwdScope()
    {
        std::error_code ec;
        fs::current_path(saved, ec);
    }
};

[[noreturn]] void abort_parallel_run(const char *what)
{
    // as main.cpp: a failure on one rank must not leave the others waiting in a collective
    std::cout << '\n' << " ERROR: " << what << '\n';
    std::cerr << " ERROR: " << what << std::endl;
    PHON_NS::abort_all_ranks();
}

void run_input_file(const fs::path &path, const std::optional<fs::path> &workdir)
{
    if (!on_python_main_thread())
        throw std::runtime_error("anphon: run_input_file must be called from the main thread");
    int finalized = 0;
    MPI_Finalized(&finalized);
    if (finalized) throw std::runtime_error("anphon: MPI has already been finalized");
    int nprocs = 1;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    // No Python object is touched below; exceptions are C++ ones, translated by nanobind after the GIL is
    // back.
    nb::gil_scoped_release release;

    try {
        const auto input = fs::absolute(path).string();
        std::optional<fs::path> dir;
        if (workdir) {
            dir = fs::absolute(*workdir);
            if (!fs::is_directory(*dir))
                throw std::invalid_argument("anphon: workdir is not a directory: " + dir->string());
        }

        RunStateScope state(nprocs == 1);
        CwdScope cwd(dir);
        std::string arg0 = "anphon", arg1 = input;
        char *argv[] = {arg0.data(), arg1.data()};
        PHON_NS::PhononCUI phon_cui; // fresh PhononCUI (and PHON) per call
        phon_cui.run(2, argv, MPI_COMM_WORLD, /*embedded=*/true);
        // rank 0 writes the outputs: let no rank return before they are complete
        MPI_Barrier(MPI_COMM_WORLD);
    } catch (const std::exception &e) {
        if (nprocs > 1) abort_parallel_run(e.what());
        throw;
    } catch (...) {
        if (nprocs > 1) abort_parallel_run("unknown exception");
        throw std::runtime_error("anphon: unknown exception");
    }
}

// For the Python layer: a failure before a collective run on one rank must not leave the other ranks
// waiting inside it (python/anphon/_api.py, run()).
int comm_size()
{
    int finalized = 0;
    MPI_Finalized(&finalized);
    if (finalized) throw std::runtime_error("anphon: MPI has already been finalized");
    int nprocs = 1;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    return nprocs;
}

void finalize_mpi()
{
    if (!mpi_initialized_here) return;
    int finalized = 0;
    MPI_Finalized(&finalized);
    if (finalized) return;
    std::cout.flush();
    MPI_Finalize();
}
} // namespace

NB_MODULE(_anphon, m)
{
    m.doc() = "nanobind binding that runs ALAMODE anphon (experimental)";

    int initialized = 0, finalized = 0;
    MPI_Initialized(&initialized);
    MPI_Finalized(&finalized);
    if (finalized) throw nb::import_error("anphon: MPI has already been finalized");
    int provided = MPI_THREAD_SINGLE;
    if (initialized) {
        MPI_Query_thread(&provided);
    } else {
        // the thread that initialises MPI must be the one that runs and finalises it
        if (!on_python_main_thread()) throw nb::import_error("anphon: import the module from the main thread");
        // FUNNELED: MPI calls come from the main thread only (anphon calls MPI outside its OpenMP regions)
        MPI_Init_thread(nullptr, nullptr, MPI_THREAD_FUNNELED, &provided);
        mpi_initialized_here = true;
        nb::module_::import_("atexit").attr("register")(nb::cpp_function(&finalize_mpi));
    }
    if (provided < MPI_THREAD_FUNNELED) {
        throw nb::import_error("anphon: the MPI library provides less than MPI_THREAD_FUNNELED");
    }

    m.def("run_input_file",
          &run_input_file,
          nb::arg("path"),
          nb::arg("workdir") = nb::none(),
          "Run anphon on the input file `path`, as the CLI `anphon path` would, on MPI_COMM_WORLD.\n\n"
          "Relative file names inside the input (FCSFILE, ...) and all outputs are relative to `workdir`\n"
          "(an existing directory; default: the current directory). The cwd is restored afterwards.\n"
          "Collective under MPI. An error raises RuntimeError in a single-process run and aborts all\n"
          "ranks otherwise. Main thread only (RuntimeError otherwise); the GIL is released during the run.");
    m.def("comm_size", &comm_size, "Number of processes in MPI_COMM_WORLD.");
    m.def(
        "abort_all",
        [](const std::string &message) { abort_parallel_run(message.c_str()); },
        nb::arg("message"),
        "Print ' ERROR: message' and abort all MPI ranks (does not return).");
}
