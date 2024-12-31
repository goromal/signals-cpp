#pragma once

#include "signals/Signal.h"
#include "signals/Integration.h"
#include "signals/State.h"
#include "signals/Models.h"

/**
 * @mainpage signals-cpp Library Documentation
 *
 * @section intro_sec Introduction
 * Header-only templated C++ library implementing rigid-body dynamics, derivatives, integrals, and interpolation.
 *
 * @section install Installation
 * This code is meant to be built as a static library with CMake. It should be compatible with the latest versions of
 * Eigen and Boost (unit test framework only).
 * The library [manif-geom-cpp](https://github.com/goromal/manif-geom-cpp) must also be installed.
 *
 * Install with
 *
 * ```bash
 * mkdir build
 * cd build
 * cmake ..
 * make # or make install
 * ```
 *
 * By default, building will also build and run the unit tests, but this can be turned off with the CMake option
 * `BUILD_TESTS`.
 */
