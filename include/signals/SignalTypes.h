#pragma once

#include <Eigen/Core>
#include <SO3.h>
#include <SE3.h>

#include "Signal.h"

using namespace Eigen;

template<typename T, size_t n>
class VectorSignal : public Signal<T, Matrix<T, n, 1>, n>
{
    SignalType getZeroSignal() override
    {
        return SignalType::Zero();
    }
    SignalType getNansSignal() override
    {
        return SignalType::Zero(); // TODO fix
    }
    TangentType getZeroTangent() override
    {
        return TangentType::Zero();
    }
    TangentType getNansTangent() override
    {
        return TangentType::Zero(); // TODO fix
    }
};

template<typename T, typename ManifType, size_t n>
class ManifoldSignal : public Signal<T, ManifType, n>
{
    SignalType getZeroSignal() override
    {
        return SignalType::identity();
    }
    SignalType getNansSignal() override
    {
        return SignalType::identity(); // TODO fix
    }
    TangentType getZeroTangent() override
    {
        return TangentType::Zero();
    }
    TangentType getNansTangent() override
    {
        return TangentType::Zero(); // TODO fix
    }
};

template<typename T>
using SO3Signal = ManifoldSignal<T, SO3<T>, 3>;

template<typename T>
using SE3Signal = ManifoldSignal<T, SE3<T>, 6>;

typedef VectorSignal<double, 1>  Vector1dSignal;
typedef VectorSignal<double, 2>  Vector2dSignal;
typedef VectorSignal<double, 3>  Vector3dSignal;
typedef VectorSignal<double, 4>  Vector4dSignal;
typedef VectorSignal<double, 5>  Vector5dSignal;
typedef VectorSignal<double, 6>  Vector6dSignal;
typedef VectorSignal<double, 7>  Vector7dSignal;
typedef VectorSignal<double, 8>  Vector8dSignal;
typedef VectorSignal<double, 9>  Vector9dSignal;
typedef VectorSignal<double, 10> Vector10dSignal;

typedef SO3Signal<double> SO3dSignal;
typedef SE3Signal<double> SE3dSignal;
