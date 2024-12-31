#pragma once
#include <optional>
#include "signals/State.h"
#include "signals/Integration.h"

using namespace Eigen;

/**
 * @brief Base type for all model simulators.
 *
 * Provides methods for initialization, reset, and simulation of dynamic models.
 *
 * **Derived types**:
 *
 * - `Translational1DOFModel(d)`
 * - `Translational2DOFModel(d)`
 * - `Translational3DOFModel(d)`
 * - `Rotational1DOFModel(d)`
 * - `Rotational3DOFModel(d)`
 * - `RigidBody3DOFModel(d)`
 * - `RigidBody6DOFModel(d)`
 */
template<typename DynamicsType>
class Model
{
public:
    using InputSignalType    = typename DynamicsType::InputSignalType;
    using StateDotSignalType = typename DynamicsType::StateDotSignalType;
    using StateSignalType    = typename DynamicsType::StateSignalType;
    using ParamsType         = typename DynamicsType::ParamsType;

    /**
     * @brief Model state signal.
     */
    StateSignalType x;
    /**
     * @brief Time derivative of the model state signal.
     */
    StateDotSignalType xdot;

    Model() : params_{std::nullopt}
    {
        reset();
    }

    /**
     * @brief Initialize the model with any required parameters.
     *
     * The required parameters are determined by the model / dynamics specialization.
     */
    void setParams(const ParamsType& params)
    {
        params_ = params;
    }

    /**
     * @brief Verify that the model has parameters explicity set by `setParams()`.
     */
    bool hasParams()
    {
        return params_.has_value();
    }

    /**
     * @brief Zero out the model state and derivative variables and reset simulation time to zero.
     */
    void reset()
    {
        x.reset();
        xdot.reset();
    }

    /**
     * @brief Get the current simulation time.
     */
    double t() const
    {
        return x.t();
    }

    /**
     * @brief Simulate the system response to an input over a specified time interval.
     * @param u The input signal, which should be defined up to tf.
     * @param tf The time to simulate to. Ideally the delta from the current time is small.
     * @param insertIntoHistory Whether to store the result in state memory.
     * @param calculateXddot Whether to use finite differencing to calculate the second time derivative of the state.
     * @returns Whether the simulation was successful.
     */
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

    /**
     * @brief Simulate the system response to an input over a specified time interval, chunked up into smaller
     * integration increments.
     * @param u The input signal, which should be defined up to tf.
     * @param tf The time to simulate to.
     * @param dt Time delta length by which to chunk up the integrations. Ideally this is small.
     * @param insertIntoHistory Whether to store the result in state memory.
     * @param calculateXddot Whether to use finite differencing to calculate the second time derivative of the state.
     * @returns Whether the simulation was successful.
     */
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

/**
 * @brief Parameters for a 1D rigid body model.
 *
 * Picture a point mass existing on a line, whether horizontal or vertical.
 */
struct RigidBodyParams1D
{
    RigidBodyParams1D() : m(1), g(0) {}
    /**
     * @brief Model mass.
     */
    double m;
    /**
     * @brief Gravitational constant.
     *
     * Essentially defines which way is "down." Set to zero if e.g., the 1D dimension is horizontal.
     */
    double g;
};

/**
 * @brief Parameters for a 2D rigid body model.
 *
 * Picture a mass (not necessarily a point mass) confined to a 2D plane.
 */
struct RigidBodyParams2D
{
    RigidBodyParams2D() : m(1), J(1), g(Vector2d::Zero()) {}
    /**
     * @brief Model mass.
     */
    double m;
    /**
     * @brief Moment of inertia.
     *
     * The moment of inertia is about the axis coming out of the 2D plane.
     */
    double J;
    /**
     * @brief Gravitational vector.
     *
     * Essentially defines which way is "down." Set to all zeroes if e.g., the 2D plane represents flat ground.
     */
    Vector2d g;
};

/**
 * @brief Parameters for a 3D rigid body model.
 *
 * Picture a mass (not necessarily a point mass) free to move around 3D space.
 */
struct RigidBodyParams3D
{
    RigidBodyParams3D() : m(1), J(Matrix3d::Identity()), g(Vector3d::Zero()) {}
    /**
     * @brief Model mass.
     */
    double m;
    /**
     * @brief Moment of inertia.
     *
     * Moments of inertia about all three principal axes, represented as \f$\boldsymbol{J}\in\mathbb{R}^{3\times 3}\f$.
     */
    Matrix3d J;
    /**
     * @brief Gravitational vector.
     *
     * Essentially defines which way is "down." Set to all zeroes if there's no gravity.
     */
    Vector3d g;
};

/**
 * @brief Base class for defining the dynamics for 1D, 2D, and 3D *point masses*.
 */
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

    /**
     * @brief Update a provided state time derivative given an input and time interval.
     * @param xdot The state time derivative signal to update.
     * @param x The state signal to reference for the dynamics.
     * @param u The input signal to reference for the dynamics.
     * @param t0 The time at which to sample the state and input.
     * @param tf The time at which to modify the state time derivative.
     * @param params The rigid body model parameters.
     * @param insertIntoHistory Whether to insert the answer into state time derivative signal history.
     * @param calculateXddot Whether to use finite differencing to calculate the second time derivative of the state.
     */
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

/**
 * @brief Definition of the dynamics for planar rotation-only motion of a mass.
 */
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

    /**
     * @brief Update a provided state time derivative given an input and time interval.
     * @param xdot The state time derivative signal to update.
     * @param x The state signal to reference for the dynamics.
     * @param u The input signal to reference for the dynamics.
     * @param t0 The time at which to sample the state and input.
     * @param tf The time at which to modify the state time derivative.
     * @param params The 2D rigid body model parameters.
     * @param insertIntoHistory Whether to insert the answer into state time derivative signal history.
     * @param calculateXddot Whether to use finite differencing to calculate the second time derivative of the state.
     */
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

/**
 * @brief Definition of the dynamics for 3D rotation-only motion of a mass.
 */
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

    /**
     * @brief Update a provided state time derivative given an input and time interval.
     * @param xdot The state time derivative signal to update.
     * @param x The state signal to reference for the dynamics.
     * @param u The input signal to reference for the dynamics.
     * @param t0 The time at which to sample the state and input.
     * @param tf The time at which to modify the state time derivative.
     * @param params The 3D rigid body model parameters.
     * @param insertIntoHistory Whether to insert the answer into state time derivative signal history.
     * @param calculateXddot Whether to use finite differencing to calculate the second time derivative of the state.
     */
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

/**
 * @brief Definition of the dynamics for planar motion of a mass that's allowed to rotate.
 */
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

    /**
     * @brief Update a provided state time derivative given an input and time interval.
     * @param xdot The state time derivative signal to update.
     * @param x The state signal to reference for the dynamics.
     * @param u The input signal to reference for the dynamics.
     * @param t0 The time at which to sample the state and input.
     * @param tf The time at which to modify the state time derivative.
     * @param params The 2D rigid body model parameters.
     * @param insertIntoHistory Whether to insert the answer into state time derivative signal history.
     * @param calculateXddot Whether to use finite differencing to calculate the second time derivative of the state.
     */
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

/**
 * @brief Definition of the dynamics for a 3D rigid body that can rotate about any axis.
 */
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

    /**
     * @brief Update a provided state time derivative given an input and time interval.
     * @param xdot The state time derivative signal to update.
     * @param x The state signal to reference for the dynamics.
     * @param u The input signal to reference for the dynamics.
     * @param t0 The time at which to sample the state and input.
     * @param tf The time at which to modify the state time derivative.
     * @param params The 3D rigid body model parameters.
     * @param insertIntoHistory Whether to insert the answer into state time derivative signal history.
     * @param calculateXddot Whether to use finite differencing to calculate the second time derivative of the state.
     */
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

#define MAKE_MODEL(ModelBaseName, ModelDOF)                                                                            \
    template<typename T>                                                                                               \
    using ModelBaseName##ModelDOF##Model = Model<ModelBaseName##Dynamics##ModelDOF<T>>;                                \
    typedef ModelBaseName##ModelDOF##Model<double> ModelBaseName##ModelDOF##Modeld;

MAKE_MODEL(Translational, 1DOF)
MAKE_MODEL(Translational, 2DOF)
MAKE_MODEL(Translational, 3DOF)
MAKE_MODEL(Rotational, 1DOF)
MAKE_MODEL(Rotational, 3DOF)
MAKE_MODEL(RigidBody, 3DOF)
MAKE_MODEL(RigidBody, 6DOF)
