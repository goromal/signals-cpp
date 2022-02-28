#pragma once
#include <Eigen/Core>
#include "signals/Signal.h"
#include "signals/SignalTypes.h"

using namespace Eigen;

template<typename T>
VectorSignal<T, SE3Signal<T>::tDim> SE3RigidBodyDynamics(const SE3Signal<T>& x, const VectorSignal<T, 6>& u)
{
    // TODO transfer settings? non-clunky computation? auto-diff?
}

// struct SE3RigidBodyDynamics
// {
//     template <typename T>
//     bool f()
// };

// struct FixedWingAttitudeEstimation
// {
//   template <typename Derived1, typename Derived2, typename Derived3>
//   bool f(Derived1 &xdot, const Derived2 &x, const Derived3 &u) const
//   {
//     typedef typename Derived1::Scalar T;
//     T phi = x(0);
//     T theta = x(1);
//     T p = u(0);
//     T q = u(1);
//     T r = u(2);
//     T Va = u(3);

//     xdot(0) = p + q * sin(phi) * tan(theta) + r * cos(phi) * tan(theta);
//     xdot(1) = q * cos(phi) - r * sin(phi);

//     return true;
//   }
// };

// TODO
