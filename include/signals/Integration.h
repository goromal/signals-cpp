#pragma once
#include "signals/Signal.h"

template<typename IntegratorType>
class Integrator
{
private:
    bool
    _get_dt(double& dt, const double& t0, const double& tf, const double& dt_max = std::numeric_limits<double>::max())
    {
        if (t0 >= tf || t0 < 0)
        {
            return false;
        }
        dt = std::min(tf - t0, dt_max);
        return true;
    }

public:
    Integrator() {}

    template<typename BaseSignalSpec, typename TangentSignalSpec>
    bool operator()(Signal<BaseSignalSpec, TangentSignalSpec>&          xInt,
                    const Signal<TangentSignalSpec, TangentSignalSpec>& x,
                    const double&                                       t,
                    const bool&                                         insertIntoHistory = false)
    {
        double dt;
        if (!_get_dt(dt, xInt.t(), t))
        {
            return false;
        }
        return IntegratorType::Integrate(xInt, x, xInt.t(), t, insertIntoHistory);
    }

    template<typename BaseSignalSpec, typename TangentSignalSpec>
    bool operator()(Signal<BaseSignalSpec, TangentSignalSpec>&          xInt,
                    const Signal<TangentSignalSpec, TangentSignalSpec>& x,
                    const double&                                       t,
                    const double&                                       dt,
                    const bool&                                         insertIntoHistory = false)
    {
        double t_k     = xInt.t();
        bool   success = true;
        while (t_k < t && success)
        {
            double dt_k;
            if (_get_dt(dt_k, t_k, t, dt))
            {
                double t_kp1 = t_k + dt_k;
                success &= IntegratorType::Integrate(xInt, x, t_k, t_kp1, insertIntoHistory);
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

struct EulerIntegrator
{
    template<typename BaseSignalSpec, typename TangentSignalSpec>
    static bool Integrate(Signal<BaseSignalSpec, TangentSignalSpec>&          xInt,
                          const Signal<TangentSignalSpec, TangentSignalSpec>& x,
                          const double&                                       t0,
                          const double&                                       tf,
                          const bool&                                         insertIntoHistory)
    {
        double dt = tf - t0;
        return xInt.update(tf, xInt() + x(tf) * dt, x(tf), insertIntoHistory);
    }
};

struct TrapezoidalIntegrator
{
    template<typename BaseSignalSpec, typename TangentSignalSpec>
    static bool Integrate(Signal<BaseSignalSpec, TangentSignalSpec>&          xInt,
                          const Signal<TangentSignalSpec, TangentSignalSpec>& x,
                          const double&                                       t0,
                          const double&                                       tf,
                          const bool&                                         insertIntoHistory)
    {
        double dt = tf - t0;
        return xInt.update(tf, xInt() + (xInt.dot() + x(tf)) * dt / 2.0, x(tf), insertIntoHistory);
    }
};

typedef Integrator<EulerIntegrator>       IntegrateEuler;
typedef Integrator<TrapezoidalIntegrator> IntegrateTrapezoidal;
