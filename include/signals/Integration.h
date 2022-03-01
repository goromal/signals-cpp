#pragma once

// TODO lower-order integrators

template <typename DynamicsType>
struct DynamicsIntegrateRK4
{
    template<typename T>
    bool operator()(DynamicsType::StateType<T> &x, const DynamicsType::InputType &u, const double &dt)
    {
        // TODO += the time already in the signal
        // do we need to update the time logic in Signal?
        return true;
    }
};
