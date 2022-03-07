#pragma once
#include <Eigen/Core>
#include "signals/Signal.h"

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

    Matrix3d J_;
    Matrix3d J_inv_;
    SO3RigidBodyDynamics(const Matrix3d& J)
    {
        J_     = J;
        J_inv_ = J.inverse();
    }

    template<typename T>
    bool
    operator()(DeltaType<T>& xdot, const StateType<T>& x, const InputType<T>& u, const bool& insertIntoHistory = false)
    {
        double t = x.t();
        xdot.update(t, x.dot(), J_inv_ * (-x.dot().cross(J_ * x.dot()) + u));

        return true;
    }
};

/**
 * Position + Velocity: World Frame
 * Rotation: Body -> World
 * Omega: Body Frame
 */
struct SE3RigidBodyDynamics
{
    typedef Vector6Signal DeltaType;
    typedef SE3Signal     StateType;
    typedef Vector6Signal InputType;

    double   m_;
    Vector3d g_;
    Matrix3d J_;
    Matrix3d J_inv_;
    SE3RigidBodyDynamics(const double& m, const Matrix3d& J, const Vector3d& g)
    {
        m_     = m;
        g_     = g;
        J_     = J;
        J_inv_ = J.inverse();
    }

    template<typename T>
    bool
    operator()(DeltaType<T>& xdot, const StateType<T>& x, const InputType<T>& u, const bool& insertIntoHistory = false)
    {
        double          t    = x.t();
        Matrix<T, 3, 1> v    = x.dot().block<3, 1>(0, 0);
        Matrix<T, 3, 1> w    = x.dot().block<3, 1>(3, 0);
        Matrix<T, 3, 1> vdot = -g + 1.0 / m_ * x().q() * u().block<3, 1>(0, 0);
        Matrix<T, 3, 1> wdot = J_inv_ * (-w.cross(J_ * w) + u().block<3, 1>(3, 0));

        xdot.update(t, (Matrix<T, 6, 1>() << v, w).finished(), (Matrix<T, 6, 1>() << vdot, wdot).finished());

        return true;
    }
};

template<typename DynamicsType>
struct SimulateEuler
{
    DynamicsType f;
    template<typename T>
    bool operator()(DynamicsType::StateType<T>&    x,
                    const DynamicsType::InputType& u,
                    const double&                  t,
                    const bool&                    insertIntoHistory = false)
    {
        const double dt = t - x.t();

        DynamicsType::DeltaType<T> k1, dx;

        f(k1, x, u);

        dx = k1 * dt;
        x.update(t, x() + dx(), x.dot() + dx.dot(), insertIntoHistory);

        return true;
    }
};

template<typename DynamicsType>
struct SimulateTrapezoidal
{
    DynamicsType f;
    template<typename T>
    bool operator()(DynamicsType::StateType<T>&    x,
                    const DynamicsType::InputType& u,
                    const double&                  t,
                    const bool&                    insertIntoHistory = false)
    {
        const double dt = t - x.t();

        DynamicsType::DeltaType<T> k1, k2, dx;

        f(k1, x, u);
        f(k2, x + k1 * dt, u);

        dx = (k1 + k2) * dt / 2.0;
        x.update(t, x() + dx(), x.dot() + dx.dot(), insertIntoHistory);

        return true;
    }
};

template<typename DynamicsType>
struct SimulateRK4
{
    DynamicsType f;
    template<typename T>
    bool operator()(DynamicsType::StateType<T>&    x,
                    const DynamicsType::InputType& u,
                    const double&                  t,
                    const bool&                    insertIntoHistory = false)
    {
        const double dt = t - x.t();

        DynamicsType::DeltaType<T> k1, k2, k3, k4, dx;

        f(k1, x, u);
        f(k2, x + k1 * dt / 2.0, u);
        f(k3, x + k2 * dt / 2.0, u);
        f(k4, x + k3 * dt, u);

        dx = (k1 + 2.0 * k2 + 2.0 * k3 + k4) * dt / 6.0;
        x.update(t, x() + dx(), x.dot() + dx.dot(), insertIntoHistory);

        return true;
    }
};

typedef SimulateEuler<SO3RigidBodyDynamics> SO3RigidBodySimulateEuler;
typedef SimulateEuler<SE3RigidBodyDynamics> SE3RigidBodySimulateEuler;

typedef SimulateTrapezoidal<SO3RigidBodyDynamics> SO3RigidBodySimulateTrapezoidal;
typedef SimulateTrapezoidal<SE3RigidBodyDynamics> SE3RigidBodySimulateTrapezoidal;

typedef SimulateRK4<SO3RigidBodyDynamics> SO3RigidBodySimulateRK4;
typedef SimulateRK4<SE3RigidBodyDynamics> SE3RigidBodySimulateRK4;
