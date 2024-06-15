#pragma once
#include <algorithm>

namespace signal_utils
{
inline bool __getTimeDelta(double&       dt,
                           const double& t0,
                           const double& tf,
                           const double& dt_max = std::numeric_limits<double>::max())
{
    if (t0 >= tf || t0 < 0)
    {
        return false;
    }
    dt = std::min(tf - t0, dt_max);
    return true;
}

} // end namespace signal_utils
