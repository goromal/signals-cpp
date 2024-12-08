#pragma once
#include "signals/Signal.h"

template<typename IntegratorType>
struct Integrator
{
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

struct EulerIntegratorSpec
{
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

struct TrapezoidalIntegratorSpec
{
    template<typename BaseSignalSpec, typename TangentSignalSpec>
    static bool integrate(Signal<BaseSignalSpec, TangentSignalSpec>&          xInt,
                          const Signal<TangentSignalSpec, TangentSignalSpec>& x,
                          const double&                                       t0,
                          const double&                                       tf,
                          const bool&                                         insertIntoHistory)
    {
        double dt = tf - t0;
        return xInt.update(tf, xInt() + (xInt.dot() + x(tf)) * dt / 2.0, x(tf), insertIntoHistory);
    }
};

// TODO SimpsonIntegratorSpec

typedef Integrator<EulerIntegratorSpec>       EulerIntegrator;
typedef Integrator<TrapezoidalIntegratorSpec> TrapezoidalIntegrator;
