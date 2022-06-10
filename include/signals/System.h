#pragma once
#include "signals/State.h"
#include "signals/Integration.h"

using namespace Eigen;

template<typename DynamicsType>
struct System
{
    using InputSignalType    = typename DynamicsType::InputSignalType;
    using StateDotSignalType = typename DynamicsType::StateDotSignalType;
    using StateSignalType    = typename DynamicsType::StateSignalType;

    StateSignalType    x; // e.g., sys.x(3.2).pose
    StateDotSignalType xdot;
    DynamicsType       dynamics; // Must set parameters before using

    System()
    {
        reset();
    }

    void reset()
    {
        x.reset();
        xdot.reset();
    }

    double t() const
    {
        return x.t();
    }

    template<typename IntegratorType>
    bool simulate(const InputSignalType& u,
                  const double&          tf,
                  const bool&            insertIntoHistory = false,
                  const bool&            calculateXddot    = false)
    {
        IntegratorType integrator; // TODO clunky
        if (!dynamics(xdot, x, u, x.t(), tf, insertIntoHistory, calculateXddot))
        {
            return false;
        }
        if (!integrator(x, xdot, tf, insertIntoHistory))
        {
            return false;
        }
        return true;
    }

    // TODO get statespace with auto diff
};

template<typename IST, typename SST, typename SDST, size_t d>
struct TranslationalDynamicsBase
{
    using InputSignalType    = IST;
    using StateDotSignalType = SST;
    using StateSignalType    = SDST;

    using InputType       = typename InputSignalType::BaseType;
    using StateType       = typename StateSignalType::BaseType;
    using StateDotType    = typename StateDotSignalType::BaseType;
    using StateDotDotType = typename StateDotSignalType::TangentType;

    using SpaceType = typename StateType::PoseType; // TODO replicate this pattern

    // Set these quantities to characterize the dynamics
    double    m = 1.0;
    SpaceType g;

    bool operator()(StateDotSignalType&    xdot,
                    const StateSignalType& x,
                    const InputSignalType& u,
                    const double&          t0,
                    const double&          tf,
                    const bool&            insertIntoHistory = false,
                    const bool&            calculateXddot    = false)
    {
        InputType u_k = u(t0);
        StateType x_k = x(t0);

        StateDotType xdot_k;
        xdot_k.pose  = x_k.twist;
        xdot_k.twist = -g + u_k / m;

        if (calculateXddot)
        {
            return xdot.update(tf, xdot_k, insertIntoHistory);
        }
        else
        {
            return xdot.update(tf, xdot_k, StateDotDotType::identity(), insertIntoHistory);
        }
    }
};

template<typename T>
using TranslationalDynamics1DOF =
    TranslationalDynamicsBase<ScalarSignal<T>, ScalarStateSignal<T>, ScalarStateSignal<T>, 1>;

template<typename T>
using TranslationalDynamics2DOF =
    TranslationalDynamicsBase<Vector2Signal<T>, Vector2StateSignal<T>, Vector2StateSignal<T>, 2>;

template<typename T>
using TranslationalDynamics3DOF =
    TranslationalDynamicsBase<Vector3Signal<T>, Vector3StateSignal<T>, Vector3StateSignal<T>, 3>;

template<typename T>
struct RotationalDynamics3DOF
{
    using InputSignalType    = Vector3Signal<T>;
    using StateSignalType    = SO3StateSignal<T>;
    using StateDotSignalType = Vector3StateSignal<T>;

    using InputType       = typename InputSignalType::BaseType;
    using StateType       = typename StateSignalType::BaseType;
    using StateDotType    = typename StateDotSignalType::BaseType;
    using StateDotDotType = typename StateDotSignalType::TangentType;

    using SpaceType = Vector3d;

    Matrix3d J = Matrix3d::Identity();

    bool operator()(StateDotSignalType&    xdot,
                    const StateSignalType& x,
                    const InputSignalType& u,
                    const double&          t0,
                    const double&          tf,
                    const bool&            insertIntoHistory = false,
                    const bool&            calculateXddot    = false)
    {
        InputType u_k = u(t0);
        StateType x_k = x(t0);

        StateDotDotType xdot_k;
        xdot_k.pose  = x_k.twist;
        xdot_k.twist = J.inverse() * (-x_k.twist.cross(J * x_k.twist) + u_k);

        if (calculateXddot)
        {
            return xdot.update(tf, xdot_k, insertIntoHistory);
        }
        else
        {
            return xdot.update(tf, xdot_k, StateDotDotType::identity(), insertIntoHistory);
        }
    }
};

template<typename T>
struct RigidBodyDynamics6DOF
{
    using InputSignalType    = Vector6Signal<T>;
    using StateSignalType    = SE3StateSignal<T>;
    using StateDotSignalType = Vector6StateSignal<T>;

    using InputType       = typename InputSignalType::BaseType;
    using StateType       = typename StateSignalType::BaseType;
    using StateDotType    = typename StateDotSignalType::BaseType;
    using StateDotDotType = typename StateDotSignalType::TangentType;

    using SpaceType = Vector3d;

    double    m = 1.0;
    SpaceType g = Vector3d::Zero();
    Matrix3d  J = Matrix3d::Identity();

    bool operator()(StateDotSignalType&    xdot,
                    const StateSignalType& x,
                    const InputSignalType& u,
                    const double&          t0,
                    const double&          tf,
                    const bool&            insertIntoHistory = false,
                    const bool&            calculateXddot    = false)
    {
        InputType u_k = u(t0);
        StateType x_k = x(t0);

        StateDotDotType xdot_k;
        xdot_k.pose                             = x_k.twist;
        xdot_k.twist.template block<3, 1>(0, 0) = -g + 1.0 / m * x_k.pose.q() * u_k.template block<3, 1>(0, 0);
        xdot_k.twist.template block<3, 1>(3, 0) =
            J.inverse() * (-x_k.twist.template block<3, 1>(3, 0).cross(J * x_k.twist.template block<3, 1>(3, 0)) +
                           u_k.template block<3, 1>(3, 0));

        if (calculateXddot)
        {
            return xdot.update(tf, xdot_k, insertIntoHistory);
        }
        else
        {
            return xdot.update(tf, xdot_k, StateDotDotType::identity(), insertIntoHistory);
        }
    }
};

// TODO RigidBodyDynamics6DOF [SE3]

// /**
//  * Pattern for AutoDiff:
//  *  - All args are Signal Types, no returns
//  *  - There's a Jacobian class that's associated with a Signal-type input and output
//  *
//  */

// struct SO3RigidBodyDynamics
// {
//     template<typename T>
//     using DeltaType = Vector3Signal<T>;
//     template<typename T>
//     using StateType = SO3Signal<T>;
//     template<typename T>
//     using InputType = Vector3Signal<T>;

//     Matrix3d J_;
//     Matrix3d J_inv_;
//     SO3RigidBodyDynamics(const Matrix3d& J)
//     {
//         J_     = J;
//         J_inv_ = J.inverse();
//     }

//     template<typename T>
//     bool
//     operator()(DeltaType<T>& xdot, const StateType<T>& x, const InputType<T>& u, const bool& insertIntoHistory =
//     false)
//     {
//         double t = x.t();
//         xdot.update(t, x.dot(), J_inv_ * (-x.dot().cross(J_ * x.dot()) + u));

//         return true;
//     }
// };

// /**
//  * Position + Velocity: World Frame
//  * Rotation: Body -> World
//  * Omega: Body Frame
//  */
// struct SE3RigidBodyDynamics
// {
//     template<typename T>
//     using DeltaType = Vector6Signal<T>;
//     template<typename T>
//     using StateType = SE3Signal<T>;
//     template<typename T>
//     using InputType = Vector6Signal<T>;

//     double   m_;
//     Vector3d g_;
//     Matrix3d J_;
//     Matrix3d J_inv_;
//     SE3RigidBodyDynamics(const double& m, const Matrix3d& J, const Vector3d& g)
//     {
//         m_     = m;
//         g_     = g;
//         J_     = J;
//         J_inv_ = J.inverse();
//     }

//     template<typename T>
//     bool
//     operator()(DeltaType<T>& xdot, const StateType<T>& x, const InputType<T>& u, const bool& insertIntoHistory =
//     false)
//     {
//         double          t    = x.t();
//         Matrix<T, 3, 1> v    = x.dot().template block<3, 1>(0, 0);
//         Matrix<T, 3, 1> w    = x.dot().template block<3, 1>(3, 0);
//         Matrix<T, 3, 1> vdot = -g_ + 1.0 / m_ * x().q() * u().template block<3, 1>(0, 0);
//         Matrix<T, 3, 1> wdot = J_inv_ * (-w.cross(J_ * w) + u().template block<3, 1>(3, 0));

//         xdot.update(t, (Matrix<T, 6, 1>() << v, w).finished(), (Matrix<T, 6, 1>() << vdot, wdot).finished());

//         return true;
//     }
// };

// template<typename DynamicsType>
// struct SimulateEuler
// {
//     using typename DynamicsType::DeltaType;
//     using typename DynamicsType::InputType;
//     using typename DynamicsType::StateType;

//     DynamicsType f;

//     bool operator()(StateType& x, const InputType& u, const double& t, const bool& insertIntoHistory = false)
//     {
//         const double dt = t - x.t();

//         DeltaType k1, dx;

//         f(k1, x, u);

//         dx = k1 * dt;
//         x.update(t, x() + dx(), x.dot() + dx.dot(), insertIntoHistory);

//         return true;
//     }
// };

// template<typename DynamicsType>
// struct SimulateTrapezoidal
// {
//     using typename DynamicsType::DeltaType;
//     using typename DynamicsType::InputType;
//     using typename DynamicsType::StateType;

//     DynamicsType f;

//     bool operator()(StateType& x, const InputType& u, const double& t, const bool& insertIntoHistory = false)
//     {
//         const double dt = t - x.t();

//         DeltaType k1, k2, dx;

//         f(k1, x, u);
//         f(k2, x + k1 * dt, u);

//         dx = (k1 + k2) * dt / 2.0;
//         x.update(t, x() + dx(), x.dot() + dx.dot(), insertIntoHistory);

//         return true;
//     }
// };

// template<typename DynamicsType>
// struct SimulateRK4
// {
//     using typename DynamicsType::DeltaType;
//     using typename DynamicsType::InputType;
//     using typename DynamicsType::StateType;

//     DynamicsType f;

//     bool operator()(StateType& x, const InputType& u, const double& t, const bool& insertIntoHistory = false)
//     {
//         const double dt = t - x.t();

//         DeltaType k1, k2, k3, k4, dx;

//         f(k1, x, u);
//         f(k2, x + k1 * dt / 2.0, u);
//         f(k3, x + k2 * dt / 2.0, u);
//         f(k4, x + k3 * dt, u);

//         dx = (k1 + 2.0 * k2 + 2.0 * k3 + k4) * dt / 6.0;
//         x.update(t, x() + dx(), x.dot() + dx.dot(), insertIntoHistory);

//         return true;
//     }
// };

// typedef SimulateEuler<SO3RigidBodyDynamics> SO3RigidBodySimulateEuler;
// typedef SimulateEuler<SE3RigidBodyDynamics> SE3RigidBodySimulateEuler;

// typedef SimulateTrapezoidal<SO3RigidBodyDynamics> SO3RigidBodySimulateTrapezoidal;
// typedef SimulateTrapezoidal<SE3RigidBodyDynamics> SE3RigidBodySimulateTrapezoidal;

// typedef SimulateRK4<SO3RigidBodyDynamics> SO3RigidBodySimulateRK4;
// typedef SimulateRK4<SE3RigidBodyDynamics> SE3RigidBodySimulateRK4;

template<typename T>
using Translational1DOFSystem = System<TranslationalDynamics1DOF<T>>;

typedef Translational1DOFSystem<double> Translational1DOFSystemd;