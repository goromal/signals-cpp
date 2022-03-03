#pragma once
#include "signals/Signal.h"

struct IntegrateEuler
{
    template<typename T, typename ST, size_t n>
    bool operator()(ManifoldSignal<T, ST, n>& xInt,
                    const VectorSignal<T, n>& x,
                    const double&             t,
                    const bool&               insertIntoHistory = false)
    {
        const double dt = t - xInt.t();
        xInt.update(t, xInt() + x() * dt, x(), insertIntoHistory);

        return true;
    }

    template<typename T, size_t n>
    bool operator()(VectorSignal<T, n>&       xInt,
                    const VectorSignal<T, n>& x,
                    const double&             t,
                    const bool&               insertIntoHistory = false)
    {
        const double dt = t - xInt.t();
        xInt.update(t, xInt() + x() * dt, x(), insertIntoHistory);

        return true;
    }
}

struct IntegrateTrapezoidal
{
    template<typename T, typename ST, size_t n>
    bool operator()(ManifoldSignal<T, ST, n>& xInt,
                    const VectorSignal<T, n>& x,
                    const double&             t,
                    const bool&               insertIntoHistory = false)
    {
        const double dt = t - xInt.t();
        xInt.update(t, xInt() + (xInt.dot() + x()) * dt / 2.0, x(), insertIntoHistory);

        return true;
    }

    template<typename T, size_t n>
    bool operator()(VectorSignal<T, n>&       xInt,
                    const VectorSignal<T, n>& x,
                    const double&             t,
                    const bool&               insertIntoHistory = false)
    {
        const double dt = t - xInt.t();
        xInt.update(t, xInt() + (xInt.dot() + x()) * dt / 2.0, x(), insertIntoHistory);

        return true;
    }
}
