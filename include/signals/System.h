#pragma once
#include <optional>
#include "signals/State.h"
#include "signals/Integration.h"

using namespace Eigen;

template<typename DynamicsType>
class System
{
public:
    using InputSignalType    = typename DynamicsType::InputSignalType;
    using StateDotSignalType = typename DynamicsType::StateDotSignalType;
    using StateSignalType    = typename DynamicsType::StateSignalType;
    using ParamsType         = typename DynamicsType::ParamsType;

    StateSignalType    x;
    StateDotSignalType xdot;

    System() : params_{std::nullopt}
    {
        reset();
    }

    void setParams(const ParamsType& params)
    {
        params_ = params;
    }

    bool hasParams()
    {
        return params_.has_value();
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
        if (!hasParams())
            return false;
        double t0 = x.t();
        double dt;
        if (!signal_utils::getTimeDelta(dt, t0, tf))
        {
            return false;
        }
        if (!DynamicsType::update(xdot, x, u, t0, tf, params_.value(), insertIntoHistory, calculateXddot))
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
        if (!hasParams())
            return false;
        double t_k = x.t();
        while (t_k < tf)
        {
            double dt_k;
            if (!signal_utils::getTimeDelta(dt_k, t_k, tf, dt))
            {
                return false;
            }
            double t_kp1 = t_k + dt_k;
            if (!DynamicsType::update(xdot, x, u, t_k, t_kp1, params_.value(), insertIntoHistory, calculateXddot))
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

private:
    std::optional<ParamsType> params_;
};

struct RigidBodyParams1D
{
    RigidBodyParams1D() : m(1), g(0) {}
    double m;
    double g;
};

struct RigidBodyParams2D
{
    RigidBodyParams2D() : m(1), J(1), g(Vector2d::Zero()) {}
    double   m;
    double   J;
    Vector2d g;
};

struct RigidBodyParams3D
{
    RigidBodyParams3D() : m(1), J(Matrix3d::Identity()), g(Vector3d::Zero()) {}
    double   m;
    Matrix3d J;
    Vector3d g;
};

template<typename IST, typename SST, typename SDST, size_t d, typename PT>
struct TranslationalDynamicsBase
{
    using InputSignalType    = IST;
    using StateDotSignalType = SST;
    using StateSignalType    = SDST;

    using InputType       = typename InputSignalType::BaseType;
    using StateType       = typename StateSignalType::BaseType;
    using StateDotType    = typename StateDotSignalType::BaseType;
    using StateDotDotType = typename StateDotSignalType::TangentType;

    using ParamsType = PT;

    static bool update(StateDotSignalType&    xdot,
                       const StateSignalType& x,
                       const InputSignalType& u,
                       const double&          t0,
                       const double&          tf,
                       const ParamsType&      params,
                       const bool&            insertIntoHistory = false,
                       const bool&            calculateXddot    = false)
    {
        InputType u_k = u(t0);
        StateType x_k = x(t0);

        StateDotType xdot_k;
        xdot_k.pose  = x_k.twist;
        xdot_k.twist = -params.g + u_k / params.m;

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
    TranslationalDynamicsBase<ScalarSignal<T>, ScalarStateSignal<T>, ScalarStateSignal<T>, 1, RigidBodyParams1D>;

template<typename T>
using TranslationalDynamics2DOF =
    TranslationalDynamicsBase<Vector2Signal<T>, Vector2StateSignal<T>, Vector2StateSignal<T>, 2, RigidBodyParams2D>;

template<typename T>
using TranslationalDynamics3DOF =
    TranslationalDynamicsBase<Vector3Signal<T>, Vector3StateSignal<T>, Vector3StateSignal<T>, 3, RigidBodyParams3D>;

template<typename T>
struct RotationalDynamics1DOF
{
    using InputSignalType    = Vector1Signal<T>;
    using StateSignalType    = SO2StateSignal<T>;
    using StateDotSignalType = Vector1StateSignal<T>;

    using InputType       = typename InputSignalType::BaseType;
    using StateType       = typename StateSignalType::BaseType;
    using StateDotType    = typename StateDotSignalType::BaseType;
    using StateDotDotType = typename StateDotSignalType::TangentType;

    using ParamsType = RigidBodyParams2D;

    static bool update(StateDotSignalType&    xdot,
                       const StateSignalType& x,
                       const InputSignalType& u,
                       const double&          t0,
                       const double&          tf,
                       const ParamsType&      params,
                       const bool&            insertIntoHistory = false,
                       const bool&            calculateXddot    = false)
    {
        InputType u_k = u(t0);
        StateType x_k = x(t0);

        StateDotType xdot_k;
        xdot_k.pose  = x_k.twist;
        xdot_k.twist = u_k / params.J;

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
struct RotationalDynamics3DOF
{
    using InputSignalType    = Vector3Signal<T>;
    using StateSignalType    = SO3StateSignal<T>;
    using StateDotSignalType = Vector3StateSignal<T>;

    using InputType       = typename InputSignalType::BaseType;
    using StateType       = typename StateSignalType::BaseType;
    using StateDotType    = typename StateDotSignalType::BaseType;
    using StateDotDotType = typename StateDotSignalType::TangentType;

    using ParamsType = RigidBodyParams3D;

    static bool update(StateDotSignalType&    xdot,
                       const StateSignalType& x,
                       const InputSignalType& u,
                       const double&          t0,
                       const double&          tf,
                       const ParamsType&      params,
                       const bool&            insertIntoHistory = false,
                       const bool&            calculateXddot    = false)
    {
        InputType u_k = u(t0);
        StateType x_k = x(t0);

        StateDotType xdot_k;
        xdot_k.pose  = x_k.twist;
        xdot_k.twist = params.J.inverse() * (-x_k.twist.cross(params.J * x_k.twist) + u_k);

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
struct RigidBodyDynamics3DOF
{
    using InputSignalType    = Vector3Signal<T>;
    using StateSignalType    = SE2StateSignal<T>;
    using StateDotSignalType = Vector3StateSignal<T>;

    using InputType       = typename InputSignalType::BaseType;
    using StateType       = typename StateSignalType::BaseType;
    using StateDotType    = typename StateDotSignalType::BaseType;
    using StateDotDotType = typename StateDotSignalType::TangentType;

    using ParamsType = RigidBodyParams2D;

    static bool update(StateDotSignalType&    xdot,
                       const StateSignalType& x,
                       const InputSignalType& u,
                       const double&          t0,
                       const double&          tf,
                       const ParamsType&      params,
                       const bool&            insertIntoHistory = false,
                       const bool&            calculateXddot    = false)
    {
        InputType u_k = u(t0);
        StateType x_k = x(t0);

        StateDotType xdot_k;
        xdot_k.pose                             = x_k.twist;
        xdot_k.twist.template block<2, 1>(0, 0) = -(x_k.pose.q().inverse() * params.g) -
                                                  x_k.twist(2) * Matrix<T, 2, 1>(-x_k.twist(1), x_k.twist(0)) +
                                                  1.0 / params.m * u_k.template block<2, 1>(0, 0);
        xdot_k.twist(2) = u_k(2) / params.J;

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

    using ParamsType = RigidBodyParams3D;

    static bool update(StateDotSignalType&    xdot,
                       const StateSignalType& x,
                       const InputSignalType& u,
                       const double&          t0,
                       const double&          tf,
                       const ParamsType&      params,
                       const bool&            insertIntoHistory = false,
                       const bool&            calculateXddot    = false)
    {
        InputType u_k = u(t0);
        StateType x_k = x(t0);

        // The entire xdot vector is in the body-frame, since we're using it as a local perturbation
        // vector for the oplus and ominus operations from SE(3). i.e., we increment the translation
        // vector and attitude jointly, unlike with many filter/controller implementations.
        StateDotType xdot_k;
        xdot_k.pose = x_k.twist;
        xdot_k.twist.template block<3, 1>(0, 0) =
            -(x_k.pose.q().inverse() * params.g) -
            x_k.twist.template block<3, 1>(3, 0).cross(x_k.twist.template block<3, 1>(0, 0)) +
            1.0 / params.m * u_k.template block<3, 1>(0, 0);
        xdot_k.twist.template block<3, 1>(3, 0) =
            params.J.inverse() *
            (-x_k.twist.template block<3, 1>(3, 0).cross(params.J * x_k.twist.template block<3, 1>(3, 0)) +
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
using Rotational1DOFSystem = System<RotationalDynamics1DOF<T>>;

template<typename T>
using Rotational3DOFSystem = System<RotationalDynamics3DOF<T>>;

template<typename T>
using RigidBody3DOFSystem = System<RigidBodyDynamics3DOF<T>>;

template<typename T>
using RigidBody6DOFSystem = System<RigidBodyDynamics6DOF<T>>;

typedef Translational1DOFSystem<double> Translational1DOFSystemd;
typedef Translational2DOFSystem<double> Translational2DOFSystemd;
typedef Translational3DOFSystem<double> Translational3DOFSystemd;
typedef Rotational1DOFSystem<double>    Rotational1DOFSystemd;
typedef Rotational3DOFSystem<double>    Rotational3DOFSystemd;
typedef RigidBody3DOFSystem<double>     RigidBody3DOFSystemd;
typedef RigidBody6DOFSystem<double>     RigidBody6DOFSystemd;
