#pragma once
#include "signals/Signal.h"

/**
 * @brief Base type for all integrators.
 *
 * Provides methods for incremental (e.g., for control) or holistic (e.g., for open-loop simulation) integrations of
 * arbitrary signal types.
 *
 * **Derived types**:
 *
 * - `EulerIntegrator`
 * - `TrapezoidalIntegrator`
 * - `SimpsonIntegrator`
 */
template<typename IntegratorType>
struct Integrator
{
    /**
     * @brief Integrate a signal from the current time to the specified end time.
     * @param xInt The output signal representing the integration.
     * @param x The input signal to be integrated.
     * @param tf The time to integrate to. Ideally the delta from the current time is small.
     * @param insertIntoHistory Whether to store the result in xInt's memory.
     * @returns Whether the integration was successful.
     */
    template<typename BaseSignalSpec, typename TangentSignalSpec>
    static bool integrate(Signal<BaseSignalSpec, TangentSignalSpec>&          xInt,
                          const Signal<TangentSignalSpec, TangentSignalSpec>& x,
                          const double&                                       tf,
                          const bool&                                         insertIntoHistory = false)
    {
        double t0 = xInt.t();
        double dt;
        if (!signal_utils::getTimeDelta(dt, t0, tf))
        {
            return false;
        }
        return IntegratorType::integrate(xInt, x, t0, tf, insertIntoHistory);
    }

    /**
     * @brief Integrate a signal from the current time to the specified end time, chunked up into smaller integration
     * increments.
     * @param xInt The output signal representing the integration.
     * @param x The input signal to be integrated.
     * @param tf The time to integrate to.
     * @param dt Time delta length by which to chunk up the integrations. Ideally this is small.
     * @param insertIntoHistory Whether to store the result in xInt's memory.
     * @returns Whether the integration was successful.
     */
    template<typename BaseSignalSpec, typename TangentSignalSpec>
    static bool integrate(Signal<BaseSignalSpec, TangentSignalSpec>&          xInt,
                          const Signal<TangentSignalSpec, TangentSignalSpec>& x,
                          const double&                                       tf,
                          const double&                                       dt,
                          const bool&                                         insertIntoHistory = false)
    {
        double t_k     = xInt.t();
        bool   success = true;
        while (t_k < tf && success)
        {
            double dt_k;
            if (signal_utils::getTimeDelta(dt_k, t_k, tf, dt))
            {
                double t_kp1 = t_k + dt_k;
                success &= IntegratorType::integrate(xInt, x, t_k, t_kp1, insertIntoHistory);
                t_k = t_kp1;
            }
            else
            {
                success = false;
            }
        }
        return success;
    }
};

/**
 * @brief Specification for numerically integrating a black box function using Euler's method.
 */
struct EulerIntegratorSpec
{
    /**
     * @brief Euler integration implementation.
     * @param xInt The output signal representing the integration.
     * @param x The input signal to be integrated.
     * @param t0 The start time for the integration.
     * @param tf The time to integrate to. Ideally the delta from the start time is small.
     * @param insertIntoHistory Whether to store the result in xInt's memory.
     * @returns Whether the integration was successful.
     *
     * Euler integration increments the current integral by
     *
     * \f$x(t_f)\Delta t\f$
     */
    template<typename BaseSignalSpec, typename TangentSignalSpec>
    static bool integrate(Signal<BaseSignalSpec, TangentSignalSpec>&          xInt,
                          const Signal<TangentSignalSpec, TangentSignalSpec>& x,
                          const double&                                       t0,
                          const double&                                       tf,
                          const bool&                                         insertIntoHistory)
    {
        double dt = tf - t0;
        return xInt.update(tf, xInt() + x(tf) * dt, x(tf), insertIntoHistory);
    }
};

/**
 * @brief Specification for numerically integrating a black box function using the Trapezoidal method.
 */
struct TrapezoidalIntegratorSpec
{
    /**
     * @brief Trapezoidal integration implementation.
     * @param xInt The output signal representing the integration.
     * @param x The input signal to be integrated.
     * @param t0 The start time for the integration.
     * @param tf The time to integrate to. Ideally the delta from the start time is small.
     * @param insertIntoHistory Whether to store the result in xInt's memory.
     * @returns Whether the integration was successful.
     *
     * Trapezoidal integration increments the current integral by
     *
     * \f$\frac{x(t_0)+x(t_f)}{2}\Delta t\f$
     */
    template<typename BaseSignalSpec, typename TangentSignalSpec>
    static bool integrate(Signal<BaseSignalSpec, TangentSignalSpec>&          xInt,
                          const Signal<TangentSignalSpec, TangentSignalSpec>& x,
                          const double&                                       t0,
                          const double&                                       tf,
                          const bool&                                         insertIntoHistory)
    {
        double dt = tf - t0;
        return xInt.update(tf, xInt() + (x(t0) + x(tf)) * dt / 2.0, x(tf), insertIntoHistory);
    }
};

/**
 * @brief Specification for numerically integrating a black box function using Simpson's method.
 */
struct SimpsonIntegratorSpec
{
    /**
     * @brief Simpson integration implementation.
     * @param xInt The output signal representing the integration.
     * @param x The input signal to be integrated.
     * @param t0 The start time for the integration.
     * @param tf The time to integrate to. Ideally the delta from the start time is small.
     * @param insertIntoHistory Whether to store the result in xInt's memory.
     * @returns Whether the integration was successful.
     *
     * Simpson's method for integration increments the current integral by
     *
     * \f$\frac{x(t_0)+4x((t_0+t_f)/2)+x(t_f)}{6}\Delta t\f$
     */
    template<typename BaseSignalSpec, typename TangentSignalSpec>
    static bool integrate(Signal<BaseSignalSpec, TangentSignalSpec>&          xInt,
                          const Signal<TangentSignalSpec, TangentSignalSpec>& x,
                          const double&                                       t0,
                          const double&                                       tf,
                          const bool&                                         insertIntoHistory)
    {
        double dt = tf - t0;
        return xInt.update(tf,
                           xInt() + (x(t0) + 4.0 * x((t0 + tf) / 2.0) + x(tf)) * dt / 6.0,
                           x(tf),
                           insertIntoHistory);
    }
};

#define MAKE_INTEGRATOR(IntegratorName) typedef Integrator<IntegratorName##IntegratorSpec> IntegratorName##Integrator;

MAKE_INTEGRATOR(Euler)
MAKE_INTEGRATOR(Trapezoidal)
MAKE_INTEGRATOR(Simpson)
