#pragma once
#include "signals/Signal.h"

template<typename IntegratorType>
class IntegratorBase
{
private:
    IntegratorType integrator_;

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
    template<typename T, typename ST, size_t n>
    bool operator()(ManifoldSignal<T, ST, n>& xInt,
                    const VectorSignal<T, n>& x,
                    const double&             t,
                    const bool&               insertIntoHistory = false)
    {
        double dt;
        if (!_get_dt(dt, xInt.t(), t))
            return false;
        return integrator_(xInt, x, dt, insertIntoHistory);
    }

    template<typename T, size_t n>
    bool operator()(VectorSignal<T, n>&       xInt,
                    const VectorSignal<T, n>& x,
                    const double&             t,
                    const bool&               insertIntoHistory = false)
    {
        double dt;
        if (!_get_dt(dt, xInt.t(), t))
            return false;
        return integrator_(xInt, x, dt, insertIntoHistory);
    }

    template<typename T, typename ST, size_t n>
    bool operator()(ManifoldSignal<T, ST, n>& xInt,
                    const VectorSignal<T, n>& x,
                    const double&             t,
                    const double&             dt,
                    const bool&               insertIntoHistory = false)
    {
        double t_k     = xInt.t();
        bool   success = true;
        while (t_k < t && success)
        {
            double dt_k;
            if (_get_dt(dt_k, t_k, t, dt))
            {
                success &= integrator_(xInt, x, t_k, dt_k, insertIntoHistory);
                t_k += dt_k;
            }
            else
            {
                success = false;
            }
        }
        return success;
    }

    template<typename T, size_t n>
    bool operator()(VectorSignal<T, n>&       xInt,
                    const VectorSignal<T, n>& x,
                    const double&             t,
                    const double&             dt,
                    const bool&               insertIntoHistory = false)
    {
        double t_k     = xInt.t();
        bool   success = true;
        while (t_k < t && success)
        {
            double dt_k;
            if (_get_dt(dt_k, t_k, t, dt))
            {
                success &= integrator_(xInt, x, t_k, dt_k, insertIntoHistory);
                t_k += dt_k;
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
    template<typename T, typename ST, size_t n>
    bool operator()(ManifoldSignal<T, ST, n>& xInt,
                    const VectorSignal<T, n>& x,
                    const double&             t,
                    const double&             dt,
                    const bool&               insertIntoHistory)
    {
        xInt.update(t, xInt(t) + x(t) * dt, x(t), insertIntoHistory);
        return true;
    }

    template<typename T, size_t n>
    bool operator()(VectorSignal<T, n>&       xInt,
                    const VectorSignal<T, n>& x,
                    const double&             t,
                    const double&             dt,
                    const bool&               insertIntoHistory)
    {
        xInt.update(t, xInt(t) + x(t) * dt, x(t), insertIntoHistory);
        return true;
    }
};

struct TrapezoidalIntegrator
{
    template<typename T, typename ST, size_t n>
    bool operator()(ManifoldSignal<T, ST, n>& xInt,
                    const VectorSignal<T, n>& x,
                    const double&             t,
                    const double&             dt,
                    const bool&               insertIntoHistory)
    {
        xInt.update(t, xInt(t) + (xInt.dot(t) + x(t)) * dt / 2.0, x(t), insertIntoHistory);
        return true;
    }

    template<typename T, size_t n>
    bool operator()(VectorSignal<T, n>&       xInt,
                    const VectorSignal<T, n>& x,
                    const double&             t,
                    const double&             dt,
                    const bool&               insertIntoHistory)
    {
        xInt.update(t, xInt(t) + (xInt.dot(t) + x(t)) * dt / 2.0, x(t), insertIntoHistory);
        return true;
    }
};

typedef IntegratorBase<EulerIntegrator>       IntegrateEuler;
typedef IntegratorBase<TrapezoidalIntegrator> IntegrateTrapezoidal;
