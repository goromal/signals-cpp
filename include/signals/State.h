#pragma once
#include "signals/Signal.h"

using namespace Eigen;

template<typename T, typename PoseTypeSpec, size_t PoseDim, typename TwistTypeSpec, size_t TwistDim>
struct State
{
    using PoseType  = typename PoseTypeSpec::Type;
    using TwistType = typename TwistTypeSpec::Type;

    PoseType  pose;
    TwistType twist;

    State() {}

    State(T* arr) : pose(arr), twist(arr + PoseDim) {}

    State(const State& other)
    {
        this->pose  = other.pose;
        this->twist = other.twist;
    }

    static State identity()
    {
        State x;
        x.pose  = PoseTypeSpec::ZeroType();
        x.twist = TwistTypeSpec::ZeroType();
        return x;
    }
};

template<typename T, typename PTS, size_t PD, typename TTS, size_t TD>
State<T, PTS, PD, TTS, TD> operator+(const State<T, PTS, PD, TTS, TD>& l, const State<T, TTS, TD, TSS, TD>& r)
{
    State<T, PTS, PD, TTS, TD> lpr = l;
    lpr.pose += r.pose;
    lpr.twist += r.twist;
    return lpr;
}

template<typename T, typename PTS, size_t PD, typename TTS, size_t TD>
State<T, TTS, TD, TSS, TD> operator-(const State<T, PTS, PD, TTS, TD>& l, const State<T, PTS, PD, TTS, TD>& r)
{
    State<T, TTS, TD, TTS, TD> lmr;
    lmr.pose  = l.pose - r.pose;
    lmr.twist = l.twist - r.twist;
    return lmr;
}

template<typename T, typename PTS, size_t PD, typename TTS, size_t TD>
State<T, PTS, PD, TTS, TD> operator*(const double& l, const State<T, PTS, PD, TTS, TD>& r)
{
    State<T, PTS, PD, TTS, TD> lr = r;
    lr.pose *= l;
    lr.twist *= l;
    return lr;
}

template<typename T, typename PTS, size_t PD, typename TTS, size_t TD>
State<T, PTS, PD, TTS, TD> operator*(const State<T, PTS, PD, TTS, TD>& l, const double& r)
{
    State<T, PTS, PD, TTS, TD> lr = l;
    l.pose *= r;
    l.twist *= r;
    return lr;
}

template<typename T>
using ScalarStateType = State<T, ScalarSignalSpec<T>, 1, ScalarSignalSpec<T>, 1>;

template<typename T, size_t d>
using VectorStateType = State<T, VectorSignalSpec<T, d>, d, VectorSignalSpec<T, d>, d>;

template<typename T, typename ManifoldType, size_t PD, size_t TD>
using ManifoldStateType = State<T, ManifoldSignalSpec<ManifoldType>, PD, VectorSignalSpec<T, TD>, TD>;

template<typename T>
struct ScalarStateSignalSpec
{
    using Type = ScalarStateType<T>;
    static Type ZeroType()
    {
        return Type::identity();
    }
    static Type NansType()
    {
        return Type::identity(); // TODO fix
    }
};

template<typename T, size_t d>
struct VectorStateSignalSpec
{
    using Type = VectorStateType<T, d>;
    static Type ZeroType()
    {
        return Type::identity();
    }
    static Type NansType()
    {
        return Type::identity(); // TODO fix
    }
};

template<typename T, typename ManifoldType, size_t PD, size_t TD>
struct ManifoldStateSignalSpec
{
    using Type = ManifoldStateType<T, ManifoldType, PD, TD>;
    static Type ZeroType()
    {
        return Type::identity();
    }
    static Type NansType()
    {
        return Type::identity(); // TODO fix
    }
};

template<typename T, size_t d>
using VectorStateSignal = Signal<VectorStateSignalSpec<T, d>, VectorStateSignalSpec<T, d>>;

template<typename T, typename ManifoldType, size_t PD, size_t TD>
using ManifoldStateSignal = Signal<ManifoldStateSignalSpec<T, ManifoldType, PD, TD>, VectorStateSignalSpec<T, TD>>;

template<typename T>
using ScalarStateSignal = Signal<ScalarStateSignalSpec<T>, ScalarStateSignalSpec<T>>;
template<typename T>
using Vector1StateSignal = VectorStateSignal<T, 1>;
template<typename T>
using Vector2StateSignal = VectorStateSignal<T, 2>;
template<typename T>
using Vector3StateSignal = VectorStateSignal<T, 3>;
template<typename T>
using Vector4StateSignal = VectorStateSignal<T, 4>;
template<typename T>
using Vector5StateSignal = VectorStateSignal<T, 5>;
template<typename T>
using Vector6StateSignal = VectorStateSignal<T, 6>;
template<typename T>
using Vector7StateSignal = VectorStateSignal<T, 7>;
template<typename T>
using Vector8StateSignal = VectorStateSignal<T, 8>;
template<typename T>
using Vector9StateSignal = VectorStateSignal<T, 9>;
template<typename T>
using Vector10StateSignal = VectorStateSignal<T, 10>;
template<typename T>
using SO3StateSignal = ManifoldStateSignal<T, SO3<T>, 4, 3>;
template<typename T>
using SE3StateSignal = ManifoldStateSignal<T, SE3<T>, 7, 6>;

typedef ScalarStateSignal<double>   ScalardStateSignal;
typedef Vector1StateSignal<double>  Vector1dStateSignal;
typedef Vector2StateSignal<double>  Vector2dStateSignal;
typedef Vector3StateSignal<double>  Vector3dStateSignal;
typedef Vector4StateSignal<double>  Vector4dStateSignal;
typedef Vector5StateSignal<double>  Vector5dStateSignal;
typedef Vector6StateSignal<double>  Vector6dStateSignal;
typedef Vector7StateSignal<double>  Vector7dStateSignal;
typedef Vector8StateSignal<double>  Vector8dStateSignal;
typedef Vector9StateSignal<double>  Vector9dStateSignal;
typedef Vector10StateSignal<double> Vector10dStateSignal;
typedef SO3StateSignal<double>      SO3dStateSignal;
typedef SE3StateSignal<double>      SE3dStateSignal;
