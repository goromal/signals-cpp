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

    StateSignalType    x;
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
        double t0 = x.t();
        double dt;
        if (!signal_utils::getTimeDelta(dt, t0, tf))
        {
            return false;
        }
        if (!dynamics(xdot, x, u, t0, tf, insertIntoHistory, calculateXddot))
        {
            return false;
        }
        if (!IntegratorType::integrate(x, xdot, tf, insertIntoHistory))
        {
            return false;
        }
        return true;
    }

    template<typename IntegratorType>
    bool simulate(const InputSignalType& u,
                  const double&          tf,
                  const double&          dt,
                  const bool&            insertIntoHistory = false,
                  const bool&            calculateXddot    = false)
    {
        double t_k = x.t();
        while (t_k < tf)
        {
            double dt_k;
            if (!signal_utils::getTimeDelta(dt_k, t_k, tf, dt))
            {
                return false;
            }
            double t_kp1 = t_k + dt_k;
            if (!dynamics(xdot, x, u, t_k, t_kp1, insertIntoHistory, calculateXddot))
            {
                return false;
            }
            if (!IntegratorType::integrate(x, xdot, t_kp1, insertIntoHistory))
            {
                return false;
            }
            t_k = t_kp1;
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

    using SpaceType = typename StateType::PoseType;

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

        StateDotType xdot_k;
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

        StateDotType xdot_k;
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

template<typename T>
using Translational1DOFSystem = System<TranslationalDynamics1DOF<T>>;

template<typename T>
using Translational2DOFSystem = System<TranslationalDynamics2DOF<T>>;

template<typename T>
using Translational3DOFSystem = System<TranslationalDynamics3DOF<T>>;

template<typename T>
using Rotational3DOFSystem = System<RotationalDynamics3DOF<T>>;

template<typename T>
using RigidBody6DOFSystem = System<RigidBodyDynamics6DOF<T>>;

typedef Translational1DOFSystem<double> Translational1DOFSystemd;
typedef Translational2DOFSystem<double> Translational2DOFSystemd;
typedef Translational3DOFSystem<double> Translational3DOFSystemd;
typedef Rotational3DOFSystem<double>    Rotational3DOFSystemd;
typedef RigidBody6DOFSystem<double>     RigidBody6DOFSystemd;