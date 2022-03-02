#pragma once
#include "signals/Signal.h"
#include "signals/Dynamics.h"

// TODO euler and trapezoidal integration w/o dynamics, for vector types only

template<typename DynamicsType>
struct DynamicsIntegrateEuler
{
    DynamicsType f;
    template<typename T>
    bool operator()(DynamicsType::StateType<T>& x, const DynamicsType::InputType& u, const double& t)
    {
        const double dt = t - x.t();

        DynamicsType::DeltaType<T> k1, dx;

        f(k1, x, u);

        dx = k1 * dt;
        x.update(t, x() + dx(), x.dot() + dx.dot());

        return true;
    }
};

template<typename DynamicsType>
struct DynamicsIntegrateTrapezoidal
{
    DynamicsType f;
    template<typename T>
    bool operator()(DynamicsType::StateType<T>& x, const DynamicsType::InputType& u, const double& t)
    {
        const double dt = t - x.t();

        DynamicsType::DeltaType<T> k1, k2, dx;

        f(k1, x, u);
        f(k2, x + k1 * dt, u);

        dx = (k1 + k2) * dt / 2.0;
        x.update(t, x() + dx(), x.dot() + dx.dot());

        return true;
    }
};

template<typename DynamicsType>
struct DynamicsIntegrateRK4
{
    DynamicsType f;
    template<typename T>
    bool operator()(DynamicsType::StateType<T>& x, const DynamicsType::InputType& u, const double& t)
    {
        const double dt = t - x.t();

        DynamicsType::DeltaType<T> k1, k2, k3, k4, dx;

        f(k1, x, u);
        f(k2, x + k1 * dt / 2.0, u);
        f(k3, x + k2 * dt / 2.0, u);
        f(k4, x + k3 * dt, u);

        dx = (k1 + 2.0 * k2 + 2.0 * k3 + k4) * dt / 6.0;
        x.update(t, x() + dx(), x.dot() + dx.dot());

        return true;
    }
};

typedef DynamicsIntegrateEuler<SO3RigidBodyDynamics> SO3RigidBodyDynamicsIntegrateEuler;
typedef DynamicsIntegrateEuler<SE3RigidBodyDynamics> SE3RigidBodyDynamicsIntegrateEuler;

typedef DynamicsIntegrateTrapezoidal<SO3RigidBodyDynamics> SO3RigidBodyDynamicsIntegrateTrapezoidal;
typedef DynamicsIntegrateTrapezoidal<SE3RigidBodyDynamics> SE3RigidBodyDynamicsIntegrateTrapezoidal;

typedef DynamicsIntegrateRK4<SO3RigidBodyDynamics> SO3RigidBodyDynamicsIntegrateRK4;
typedef DynamicsIntegrateRK4<SE3RigidBodyDynamics> SE3RigidBodyDynamicsIntegrateRK4;
