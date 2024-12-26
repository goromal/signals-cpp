# Signals Cpp

![example workflow](https://github.com/goromal/signals-cpp/actions/workflows/test.yml/badge.svg)

Header-only templated C++ library implementing rigid-body dynamics, derivatives, integrals, and interpolation.

**Under construction**

## Implemented Types

### Signal Types

Scalar signal type:

```cpp
ScalarSignal<T> x;
```

Vector signal types:

```cpp
Vector1Signal<T> x;
Vector2Signal<T> x;
Vector3Signal<T> x;
Vector4Signal<T> x;
Vector5Signal<T> x;
Vector6Signal<T> x;
Vector7Signal<T> x;
Vector8Signal<T> x;
Vector9Signal<T> x;
Vector10Signal<T> x;
```

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

### Models

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
