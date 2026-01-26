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

### Signals

Manifold signal types:

```
SO3Signal<T> x;
SE3Signal<T> x;
```

### Integrators

Euler integrator:

```cpp
IntegrateEuler f;
f(SignalType &xInt, TangentSignalType x, double t, bool insertHistory = false);
f(SignalType &xInt, TangentSignalType x, double t, double dt, bool insertHistory = false);
```

Trapezoidal integrator:

```cpp
IntegrateTrapezoidal f;
f(SignalType &xInt, TangentSignalType x, double t, bool insertHistory = false);
f(SignalType &xInt, TangentSignalType x, double t, double dt, bool insertHistory = false);
```

### System Models

*Pending implementations:*

- `TranslationalDynamics1DOF<T>`
  - Point mass system confined to a straight line.
- `TranslationalDynamics2DOF<T>`
  - Point mass system confined to a plane.
- `TranslationalDynamics3DOF<T>`
  - Point mass system in a 3-dimensional space.
- `RotationalDynamics1DOF<T>`
  - Single-axis rotating mass system fixed in space.
- `RotationalDynamics3DOF<T>`
  - Rotating mass fixed in 3D space.
- `RigidBodyDynamics3DOF<T>`
  - Rigid body system confined to a plane (unicycle model).
- `RigidBodyDynamics6DOF<T>`
  - Rigid body system in 3D space.

## Docs Generation

Generate updated docs in the `docs/` directory with

```bash
doxygen Doxyfile
```

