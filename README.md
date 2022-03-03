# Signals Cpp

Header-only templated C++ library implementing rigid-body dynamics, derivatives, integrals, and interpolation.

## Implemented Types

### Signal Types

```cpp
SO3Signal<T> x;
SE3Signal<T> x;

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

### Dynamics Models

```cpp
SO3RigidBodyDynamics f;
f(Vector3Signal<T> &xdot, SO3Signal<T> x, Vector3Signal<T> u);
```

```cpp
SE3RigidBodyDynamics f;
f(Vector6Signal<T> &xdot, SE3Signal<T> x, Vector6Signal<T> u);
```

```cpp
SO3RigidBodySimulateEuler f;
f(SO3Signal<T> &x, Vector3Signal<T> u, double t);
```
```cpp
SE3RigidBodySimulateEuler f;
f(SE3Signal<T> &x, Vector6Signal<T> u, double t);
```
```cpp
SO3RigidBodySimulateTrapezoidal f;
f(SO3Signal<T> &x, Vector3Signal<T> u, double t);
```
```cpp
SE3RigidBodySimulateTrapezoidal f;
f(SE3Signal<T> &x, Vector6Signal<T> u, double t);
```
```cpp
SO3RigidBodySimulateRK4 f;
f(SO3Signal<T> &x, Vector3Signal<T> u, double t);
```
```cpp
SE3RigidBodySimulateRK4 f;
f(SE3Signal<T> &x, Vector6Signal<T> u, double t);
```

### Integrators

```cpp
IntegrateEuler f;
f(VectorSignal<T, n> &xInt, VectorSignal<T, n> x, double t);
f(ManifoldSignal<T, ST, n> &xInt, VectorSignal<T, n> x, double t);
```

```cpp
IntegrateTrapezoidal f;
f(VectorSignal<T, n> &xInt, VectorSignal<T, n> x, double t);
f(ManifoldSignal<T, ST, n> &xInt, VectorSignal<T, n> x, double t);
```
