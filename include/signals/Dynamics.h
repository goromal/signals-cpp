#pragma once
#include <Eigen/Core>
#include "signals/Signal.h"
#include "signals/SignalTypes.h"

using namespace Eigen;

/**
 * Pattern for AutoDiff:
 *  - All args are Signal Types, no returns
 *  - There's a Jacobian class that's associated with a Signal-type input and output
 * 
*/

struct SO3RigidBodyDynamics
{
    typedef Vector3Signal DeltaType;
    typedef SO3Signal     StateType;
    typedef Vector3Signal InputType;

    template<typename T>
    bool f(DeltaType<T> &xdot, const StateType<T> &x, const InputType<T> &u)
    {
        // TODO
        return true;
    }
};

struct SE3RigidBodyDynamics
{
    typedef Vector6Signal DeltaType;
    typedef SE3Signal     StateType;
    typedef Vector6Signal InputType;

    template<typename T>
    bool f(DeltaType<T> &xdot, const StateType<T> &x, const InputType<T> &u)
    {
        // TODO
        return true;
    }
};
