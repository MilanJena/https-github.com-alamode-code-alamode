/*
 fatal_error.h

 Copyright (c) 2026 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory
 or http://opensource.org/licenses/mit-license.php for information.
*/

#pragma once

#include <cstdlib>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

#ifdef _OPENMP
#include <omp.h>
#endif

// Error policy for fatal input/runtime errors.
//
// The CLIs leave throw_on_error off: a fatal error prints its message and ends
// the process (through abort_hook when set, e.g. MPI_Abort in anphon) exactly
// as before. A Python extension module turns throw_on_error on at import, so
// the same error raises alamode::Error (RuntimeError in Python) instead of
// killing the interpreter.
//
// Contract:
// - The flag is consulted by alamode::fatal/terminate only. fatal does not
//   look at MPI: a module that runs under MPI (anphon) must set the flag only
//   when nprocs == 1 (PHON_NS::exit additionally checks the rank count).
// - Inside an OpenMP parallel region nothing throws: the error still ends the
//   process (known limit).
// - Out-of-memory paths use abort_now(), which never throws (its message must
//   not allocate, and the half-built arrays must not be unwound).
//
// The policy is per shared object: the variables and helpers below have
// hidden visibility, so each Python extension module (and each executable)
// holds its own copy and two modules in one interpreter do not share or flip
// each other's flag. `inline` alone would not give this: the static libraries
// almcxx_static / anphoncxx_static are not compiled with hidden visibility.
#if defined(__GNUC__) || defined(__clang__)
#define ALAMODE_HIDDEN __attribute__((visibility("hidden")))
#else
#define ALAMODE_HIDDEN
#endif

namespace alamode
{
ALAMODE_HIDDEN inline auto strip_blank(const std::string &text) -> std::string
{
    const auto first = text.find_first_not_of(" \n");
    if (first == std::string::npos) return "";
    return text.substr(first, text.find_last_not_of(" \n") + 1 - first);
}

// what() is the printed message without its surrounding blank lines.
struct Error: std::runtime_error
{
    explicit Error(const std::string &text) : std::runtime_error(strip_blank(text))
    {}
};

ALAMODE_HIDDEN inline bool throw_on_error = false;
ALAMODE_HIDDEN inline void (*abort_hook)() = nullptr; // run before std::exit when not throwing

// Throwing out of an OpenMP parallel region would terminate the process, so
// the policy only applies on the serial path.
ALAMODE_HIDDEN inline auto can_throw() -> bool
{
#ifdef _OPENMP
    if (omp_in_parallel()) return false;
#endif
    return throw_on_error;
}

// End the run without throwing: abort_hook (if set), then std::exit.
[[noreturn]] ALAMODE_HIDDEN inline void abort_now()
{
    std::cout.flush();
    if (abort_hook) abort_hook();
    std::exit(EXIT_FAILURE);
}

// End the run after the message has been printed: throw when the policy
// allows it, otherwise abort_now().
[[noreturn]] ALAMODE_HIDDEN inline void terminate(const std::string &text)
{
    std::cout.flush();
    if (can_throw()) throw Error(text);
    abort_now();
}

// Stream the arguments into one string with std::cout's current format flags,
// so that the text is byte-identical to streaming them into std::cout.
template <typename... Args>
ALAMODE_HIDDEN auto format_message(const Args &...args) -> std::string
{
    std::ostringstream ss;
    ss.copyfmt(std::cout);
    (ss << ... << args);
    return ss.str();
}

// Word-wrap text for the log. Every line (existing '\n' are kept) gets
// `indent` and is broken at blanks so that it fits in `width` columns;
// continuation lines also repeat the line's own leading blanks. A word longer
// than the width is not split. Columns count UTF-8 code points.
ALAMODE_HIDDEN inline auto wrap_text(const std::string &text, const std::string &indent,
                                     const std::size_t width = 90) -> std::string
{
    const auto columns = [](const std::string &s) {
        std::size_t n = 0;
        for (const unsigned char c: s) n += (c & 0xC0) != 0x80;
        return n;
    };
    std::string out;
    std::istringstream lines(text);
    std::string line;
    while (std::getline(lines, line)) {
        const auto body = line.find_first_not_of(' ');
        if (body == std::string::npos) {
            out += '\n';
            continue;
        }
        const auto lead = indent + line.substr(0, body);
        std::string cur = lead;
        std::istringstream words(line.substr(body));
        std::string word;
        while (words >> word) {
            if (cur.size() > lead.size() && columns(cur) + 1 + columns(word) > width) {
                out += cur + '\n';
                cur = lead;
            }
            if (cur.size() > lead.size()) cur += ' ';
            cur += word;
        }
        out += cur + '\n';
    }
    return out;
}

// Log form of a warning or an error: a header line, then the message wrapped
// and indented below it.
ALAMODE_HIDDEN inline auto log_message(const char *kind, const char *file, const std::string &message) -> std::string
{
    return std::string("\n ") + kind + " in " + file + '\n' + wrap_text(message, "   ");
}

// Print the streamed arguments to stdout and terminate (see above).
template <typename... Args>
[[noreturn]] ALAMODE_HIDDEN void fatal(const Args &...args)
{
    const auto text = format_message(args...);
    std::cout << text;
    terminate(text);
}
} // namespace alamode
