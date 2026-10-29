/*
 error.h

 Copyright (c) 2014-2018 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory
 or http://opensource.org/licenses/mit-license.php for information.
*/

#pragma once

#include <iostream>

#include "fatal_error.h"

namespace ALM_NS
{
inline auto warn(const char *file, const char *message) -> void
{
    std::cout << alamode::log_message("WARNING", file, message);
}

// Fatal errors: print the wrapped message, then alamode::terminate (throws
// when the Python module enabled it, std::exit otherwise; see
// include/fatal_error.h). The exception text stays on one line.
[[noreturn]] inline auto exit(const char *file, const char *message) -> void
{
    std::cout << alamode::log_message("ERROR", file, message);
    alamode::terminate(alamode::format_message("ERROR in ", file, "  MESSAGE: ", message));
}

template <typename T>
[[noreturn]] auto exit(const char *file, const char *message, const T info) -> void
{
    const auto text = alamode::format_message(message, info);
    std::cout << alamode::log_message("ERROR", file, text);
    alamode::terminate(alamode::format_message("ERROR in ", file, "  MESSAGE: ", text));
}
} // namespace ALM_NS
