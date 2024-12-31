# Signals Cpp

![example workflow](https://github.com/goromal/signals-cpp/actions/workflows/test.yml/badge.svg)

Header-only templated C++ library implementing rigid-body dynamics, derivatives, integrals, and interpolation.

View the library documentation [HERE](https://andrewtorgesen.com/signals-cpp).

## Installation

This code is meant to be built as a static library with CMake. It should be compatible with the latest versions of
Eigen and Boost (unit test framework only).
The library [manif-geom-cpp](https://github.com/goromal/manif-geom-cpp) must also be installed.

Install with

```bash
mkdir build
cd build
cmake ..
make # or make install
```

By default, building will also build and run the unit tests, but this can be turned off with the CMake option `BUILD_TESTS`.

## Docs Generation

Generate updated docs in the `docs/` directory with

```bash
doxygen Doxyfile
```

