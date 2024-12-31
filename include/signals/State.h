#pragma once
#include "signals/Signal.h"

using namespace Eigen;

/**
 * @brief Base type for all Model state representations.
 *
 * Provides a convenient union of the "pose" and "twist" (or time derivative of pose) components of a state vector and
 * defines arithmetic operations for that union.
 *
 * A state is *not* a Signal type in itself, which is why separate `*Signal` types are derived below.
 *
 * **Derived types**:
 *
 * - `Scalar(d)State` and `Scalar(d)StateSignal`
 * - `Vector1(d)State` and `Vector1(d)StateSignal`
 * - `Vector2(d)State` and `Vector2(d)StateSignal`
 * - `Vector3(d)State` and `Vector3(d)StateSignal`
 * - `Vector4(d)State` and `Vector4(d)StateSignal`
 * - `Vector5(d)State` and `Vector5(d)StateSignal`
 * - `Vector6(d)State` and `Vector6(d)StateSignal`
 * - `Vector7(d)State` and `Vector7(d)StateSignal`
 * - `Vector8(d)State` and `Vector8(d)StateSignal`
 * - `Vector9(d)State` and `Vector9(d)StateSignal`
 * - `Vector10(d)State` and `Vector10(d)StateSignal`
 * - `SO2(d)State` and `SO2(d)StateSignal`
 * - `SO3(d)State` and `SO3(d)StateSignal`
 * - `SE2(d)State` and `SE2(d)StateSignal`
 * - `SE3(d)State` and `SE3(d)StateSignal`
 */
template<typename T, typename PoseTypeSpec, size_t PoseDim, typename TwistTypeSpec, size_t TwistDim>
struct State
{
    using PoseType  = typename PoseTypeSpec::Type;
    using TwistType = typename TwistTypeSpec::Type;

    /**
     * @brief TODO
     */
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

    static State nans()
    {
        State x;
        x.pose  = PoseTypeSpec::NansType();
        x.twist = TwistTypeSpec::NansType();
        return x;
    }

    State& operator*=(const double& s)
    {
        pose *= s;
        twist *= s;
        return *this;
    }

    template<typename T2>
    State& operator+=(const State<T2, TwistTypeSpec, TwistDim, TwistTypeSpec, TwistDim>& r)
    {
        pose += r.pose;
        twist += r.twist;
        return *this;
    }
};

template<typename T, typename PTS, size_t PD, typename TTS, size_t TD>
State<T, PTS, PD, TTS, TD> operator+(const State<T, PTS, PD, TTS, TD>& l, const State<T, TTS, TD, TTS, TD>& r)
{
    State<T, PTS, PD, TTS, TD> lpr = l;
    lpr.pose += r.pose;
    lpr.twist += r.twist;
    return lpr;
}

template<typename T, typename PTS, size_t PD, typename TTS, size_t TD>
State<T, TTS, TD, TTS, TD> operator-(const State<T, PTS, PD, TTS, TD>& l, const State<T, PTS, PD, TTS, TD>& r)
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
    lr.pose *= r;
    lr.twist *= r;
    return lr;
}

template<typename T>
using ScalarStateType = State<T, ScalarSignalSpec<T>, 1, ScalarSignalSpec<T>, 1>;

template<typename T>
inline std::ostream& operator<<(std::ostream& os, const ScalarStateType<T>& x)
{
    os << "ScalarStateType: pose=" << x.pose << "; twist=" << x.twist;
    return os;
}

template<typename T>
ScalarStateType<T> operator/(const ScalarStateType<T>& l, const double& r)
{
    ScalarStateType<T> lr = l;
    lr.pose /= r;
    lr.twist /= r;
    return lr;
}

template<typename T, size_t d>
using VectorStateType = State<T, VectorSignalSpec<T, d>, d, VectorSignalSpec<T, d>, d>;

template<typename T, size_t d>
inline std::ostream& operator<<(std::ostream& os, const VectorStateType<T, d>& x)
{
    os << "VectorStateType: pose=" << x.pose.transpose() << "; twist=" << x.twist.transpose();
    return os;
}

template<typename T, size_t d>
VectorStateType<T, d> operator/(const VectorStateType<T, d>& l, const double& r)
{
    VectorStateType<T, d> lr = l;
    lr.pose /= r;
    lr.twist /= r;
    return lr;
}

template<typename T, typename ManifoldType, size_t PD, size_t TD>
using ManifoldStateType = State<T, ManifoldSignalSpec<ManifoldType>, PD, VectorSignalSpec<T, TD>, TD>;

template<typename T, typename ManifoldType, size_t PD, size_t TD>
inline std::ostream& operator<<(std::ostream& os, const ManifoldStateType<T, ManifoldType, PD, TD>& x)
{
    os << "VectorStateType: pose=" << x.pose << "; twist=" << x.twist;
    return os;
}

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
        return Type::nans();
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
        return Type::nans();
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
        return Type::nans();
    }
};

template<typename T>
using ScalarStateSignal = Signal<ScalarStateSignalSpec<T>, ScalarStateSignalSpec<T>>;

template<typename T>
inline std::ostream& operator<<(std::ostream& os, const ScalarStateSignal<T>& x)
{
    os << "ScalarStateSignal at t=" << x.t() << ": pose=" << x().pose << "; twist=" << x().twist;
    return os;
}

template<typename T, size_t d>
using VectorStateSignal = Signal<VectorStateSignalSpec<T, d>, VectorStateSignalSpec<T, d>>;

template<typename T, size_t d>
inline std::ostream& operator<<(std::ostream& os, const VectorStateSignal<T, d>& x)
{
    os << "VectorStateSignal at t=" << x.t() << ": pose=" << x().pose.transpose()
       << "; twist=" << x().twist.transpose();
    return os;
}

template<typename T, typename ManifoldType, size_t PD, size_t TD>
using ManifoldStateSignal = Signal<ManifoldStateSignalSpec<T, ManifoldType, PD, TD>, VectorStateSignalSpec<T, TD>>;

template<typename T, typename ManifoldType, size_t PD, size_t TD>
inline std::ostream& operator<<(std::ostream& os, const ManifoldStateSignal<T, ManifoldType, PD, TD>& x)
{
    os << "ManifoldStateSignal at t=" << x.t() << ": pose=" << x().pose << "; twist=" << x().twist;
    return os;
}

#define MAKE_VECTOR_STATES(Dimension)                                                                                  \
    template<typename T>                                                                                               \
    using Vector##Dimension##State = VectorStateType<T, Dimension>;                                                    \
    template<typename T>                                                                                               \
    using Vector##Dimension##StateSignal = VectorStateSignal<T, Dimension>;                                            \
    typedef Vector##Dimension##State<double>       Vector##Dimension##dState;                                          \
    typedef Vector##Dimension##StateSignal<double> Vector##Dimension##dStateSignal;

#define MAKE_MANIF_STATES(Manif, Dimension, TangentDimension)                                                          \
    template<typename T>                                                                                               \
    using Manif##State = ManifoldStateType<T, Manif<T>, Dimension, TangentDimension>;                                  \
    template<typename T>                                                                                               \
    using Manif##StateSignal = ManifoldStateSignal<T, Manif<T>, Dimension, TangentDimension>;                          \
    typedef Manif##State<double>       Manif##dState;                                                                  \
    typedef Manif##StateSignal<double> Manif##dStateSignal;

template<typename T>
using ScalarState = ScalarStateType<T>;
typedef ScalarState<double>       ScalardState;
typedef ScalarStateSignal<double> ScalardStateSignal;
MAKE_VECTOR_STATES(1)
MAKE_VECTOR_STATES(2)
MAKE_VECTOR_STATES(3)
MAKE_VECTOR_STATES(4)
MAKE_VECTOR_STATES(5)
MAKE_VECTOR_STATES(6)
MAKE_VECTOR_STATES(7)
MAKE_VECTOR_STATES(8)
MAKE_VECTOR_STATES(9)
MAKE_VECTOR_STATES(10)
MAKE_MANIF_STATES(SO2, 2, 1)
MAKE_MANIF_STATES(SO3, 4, 3)
MAKE_MANIF_STATES(SE2, 4, 3)
MAKE_MANIF_STATES(SE3, 7, 6)
