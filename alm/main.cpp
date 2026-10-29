/*
 main.cpp

 Copyright (c) 2014 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory 
 or http://opensource.org/licenses/mit-license.php for information.
*/

#include <cstdlib>
#include <exception>
#include <iostream>
#include <memory>
#include "alm_cui.h"

using namespace ALM_NS;

auto main(const int argc, char **argv) -> int
{
    // Turn uncaught exceptions (HDF5 I/O, ...) into a clean error exit
    // instead of std::terminate.
    try {
        const auto alm_cui = std::make_unique<ALMCUI>();
        alm_cui->run(argc, argv);
    } catch (const std::exception &e) {
        std::cout << '\n' << " ERROR: " << e.what() << '\n';
        std::cerr << " ERROR: " << e.what() << std::endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
