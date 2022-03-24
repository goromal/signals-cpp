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
    template<typename BaseSignalSpec, typename TangentSignalSpec>
    bool operator()(Signal<BaseSignalSpec, TangentSignalSpec>&          xInt,
                    const Signal<TangentSignalSpec, TangentSignalSpec>& x,
                    const double&                                       t,
                    const bool&                                         insertIntoHistory = false)
    {
        double dt;
        if (!_get_dt(dt, xInt.t(), t))
            return false;
        return IntegratorType::Integrate(xInt, x, t, dt, insertIntoHistory);
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
                t_k += dt_k;
                success &= IntegratorType::Integrate(xInt, x, t_k, dt_k, insertIntoHistory);
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
                          const double&                                       t,
                          const double&                                       dt,
                          const bool&                                         insertIntoHistory)
    {
        xInt.update(t, xInt() + x(t) * dt, x(t), insertIntoHistory);
        return true;
    }
};

struct TrapezoidalIntegrator
{
    template<typename BaseSignalSpec, typename TangentSignalSpec>
    static bool Integrate(Signal<BaseSignalSpec, TangentSignalSpec>&          xInt,
                          const Signal<TangentSignalSpec, TangentSignalSpec>& x,
                          const double&                                       t,
                          const double&                                       dt,
                          const bool&                                         insertIntoHistory)
    {
        xInt.update(t, xInt() + (xInt.dot() + x(t)) * dt / 2.0, x(t), insertIntoHistory);
        return true;
    }
};

typedef Integrator<EulerIntegrator>       IntegrateEuler;
typedef Integrator<TrapezoidalIntegrator> IntegrateTrapezoidal;
