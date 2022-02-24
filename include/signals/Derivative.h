#pragma once

#include "signals/SignalTypes.h"

template<typename SignalType>
class Derivative
{
public:
    using T = SignalType::ScalarType;
    using TangentType = SignalType::TangentType;
    Derivative(void) : sigma_{0.05}
    {
        reset();
    }
    Derivative(const double sigma) : sigma_{sigma}
    {
        reset();
    }
    void reset()
    {
        derivCurr_.setZero();
        derivPrev_.setZero();
        valPrev_.setZero();
        initialized_ = false;
    }
    void update(const SignalType& val, const double& dt)
    {
        if (initialized_)
        {
            derivCurr_ = (2. * sigma_ - dt) / (2. * sigma_ + dt) * derivPrev_ +
                         2. / (2. * sigma_ + dt) * (val - valPrev_);
            derivPrev_ = derivCurr_;
            valPrev_ = val;
        }
        else
        {
            derivCurr_.setZero();
            derivPrev_ = derivCurr_;
            valPrev_ = val;
            initialized_ = true;
        }
    }
    inline TangentType get() { return derivCurr_; }

private:
    double sigma_;
    bool initialized_;
    TangentType derivCurr_;
    TangentType derivPrev_;
    SignalType valPrev_;
};

typedef Derivative<Vector1dSignal> Vector1dSignalDerivative;
typedef Derivative<Vector2dSignal> Vector2dSignalDerivative;
typedef Derivative<Vector3dSignal> Vector3dSignalDerivative;
typedef Derivative<Vector4dSignal> Vector4dSignalDerivative;
typedef Derivative<Vector5dSignal> Vector5dSignalDerivative;
typedef Derivative<Vector6dSignal> Vector6dSignalDerivative;
typedef Derivative<Vector7dSignal> Vector7dSignalDerivative;

typedef Derivative<SO3dSignal> SO3dSignalDerivative;
typedef Derivative<SE3dSignal> SE3dSignalDerivative;
