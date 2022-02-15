#pragma once

#include "signals/SignalTypes.h"

template<typename SignalType>
class Derivative
{
private:
    double sigma_;
    bool initialized_;
    SignalType::TangentType derivCurr_;
    SignalType::TangentType derivPrev_;
    SignalType valPrev_;

public:
    Derivative(void) : sigma_{0.05}
    {
        reset();
    }
    Derivative(const double sigma) : sigma_{sigma}
    {
        reset();
    }
    // TODO
};

