#pragma once
#include <algorithm>
#include <limits>
#include <Eigen/Core>
#include <SO2.h>
#include <SO3.h>
#include <SE2.h>
#include <SE3.h>
#include "Utils.h"

using namespace Eigen;

/**
 * @brief Methods for interpolating signal values between discrete time points.
 */
enum InterpolationMethod
{
    ZERO_ORDER_HOLD, /**< Hold the previous value constant until the next data point. */
    LINEAR,          /**< Linearly interpolate between adjacent data points. */
    CUBIC_SPLINE     /**< Use cubic spline interpolation for smooth transitions. */
};

/**
 * @brief Methods for extrapolating signal values outside the defined time range.
 */
enum ExtrapolationMethod
{
    NANS,    /**< Return NaN (Not a Number) values outside the defined range. */
    ZEROS,   /**< Return zero values outside the defined range. */
    CLOSEST  /**< Return the closest boundary value (first or last data point). */
};

/**
 * @brief Methods for computing numerical derivatives of signals.
 */
enum DerivativeMethod
{
    DIRTY,       /**< Use a discrete-time dirty derivative filter for smoothing. */
    FINITE_DIFF  /**< Use simple finite difference approximation. */
};

/**
 * @brief Template class for time-series signals with interpolation, extrapolation, and derivative capabilities.
 *
 * The Signal class represents a time-varying signal that stores both the value and its time derivative.
 * It supports various interpolation methods for querying values at arbitrary time points, configurable
 * extrapolation behavior outside the defined time range, and automatic derivative computation.
 *
 * @tparam BaseSignalSpec Type specification for the signal values (can be scalar, vector, or manifold types).
 * @tparam TangentSignalSpec Type specification for the signal derivatives (typically vector types).
 *
 * **Example usage:**
 * ```cpp
 * ScalardSignal mySignal;
 * mySignal.setInterpolationMethod(InterpolationMethod::LINEAR);
 * mySignal.update(1.0, 5.0, true);  // time=1.0, value=5.0, store in history
 * mySignal.update(2.0, 7.0, true);  // time=2.0, value=7.0, store in history
 * double interpValue = mySignal(1.5);  // Query interpolated value at t=1.5
 * ```
 */
template<typename T, typename BaseSignalSpec, typename TangentSignalSpec>
class Signal
{
public:
    using BaseType    = typename BaseSignalSpec::Type;
    using TangentType = typename TangentSignalSpec::Type;

    /**
     * @brief Data point structure storing time, value, and derivative.
     */
    struct SignalDP
    {
        double      t;     /**< Time stamp. */
        BaseType    x;     /**< Signal value at time t. */
        TangentType xdot;  /**< Time derivative of signal at time t. */
    };

    /**
     * @brief Comparator for sorting signal data points by time.
     */
    struct
    {
        bool operator()(SignalDP a, SignalDP b) const
        {
            return a.t < b.t;
        }
    } SignalDPComparator;

    InterpolationMethod interpolationMethod;  /**< Current interpolation method. */
    ExtrapolationMethod extrapolationMethod;  /**< Current extrapolation method. */
    DerivativeMethod    derivativeMethod;     /**< Current derivative computation method. */

    /**
     * @brief Default constructor initializing with LINEAR interpolation, ZEROS extrapolation, and DIRTY derivative.
     */
    Signal()
    {
        interpolationMethod = InterpolationMethod::LINEAR;
        extrapolationMethod = ExtrapolationMethod::ZEROS;
        derivativeMethod    = DerivativeMethod::DIRTY;
        reset();
    }

    /**
     * @brief Copy constructor.
     * @param other The signal to copy from.
     */
    Signal(const Signal& other)
    {
        this->interpolationMethod = other.interpolationMethod;
        this->extrapolationMethod = other.extrapolationMethod;
        this->derivativeMethod    = other.derivativeMethod;
        this->t_                  = other.t_;
        this->x_                  = other.x_;
        this->xdot_               = other.xdot_;
        this->signalHistory_      = other.signalHistory_;
        this->needsSort_          = other.needsSort_;
    }

    /**
     * @brief Create a new signal representing the time derivative of this signal.
     * @return A signal containing the derivative values from this signal's history.
     */
    Signal<T, TangentSignalSpec, TangentSignalSpec> dotSignal()
    {
        Signal<T, TangentSignalSpec, TangentSignalSpec> signalDot;
        for (auto signalDP : signalHistory_)
        {
            signalDot.update(signalDP.t, signalDP.xdot, true);
        }
        signalDot.update(t(), dot());
        return signalDot;
    }

    template<typename S, typename BSS, typename TSS>
    friend Signal<S, BSS, TSS> operator+(const Signal<S, BSS, TSS>& l, const Signal<S, TSS, TSS>& r);

    template<typename S, typename BSS, typename TSS>
    friend Signal<S, TSS, TSS> operator-(const Signal<S, BSS, TSS>& l, const Signal<S, BSS, TSS>& r);

    template<typename S, typename BSS, typename TSS>
    friend Signal<S, BSS, TSS> operator*(const double& l, const Signal<S, BSS, TSS>& r);

    template<typename S, typename BSS, typename TSS>
    friend Signal<S, BSS, TSS> operator*(const Signal<S, BSS, TSS>& l, const double& r);

    /**
     * @brief Get the current time of the signal.
     * @return The current time value.
     */
    double t() const
    {
        return t_;
    }

    /**
     * @brief Get the current signal value.
     * @return The current value.
     */
    BaseType operator()() const
    {
        return x_;
    }

    /**
     * @brief Get the current time derivative of the signal.
     * @return The current derivative value.
     */
    TangentType dot() const
    {
        return xdot_;
    }

    /**
     * @brief Get the signal value at a specific time using interpolation/extrapolation.
     * @param t The time to query.
     * @return The interpolated/extrapolated signal value at time t.
     */
    BaseType operator()(const double& t) const
    {
        return xAt(t);
    }

    /**
     * @brief Get the signal derivative at a specific time using interpolation/extrapolation.
     * @param t The time to query.
     * @return The interpolated/extrapolated derivative value at time t.
     */
    TangentType dot(const double& t) const
    {
        return xDotAt(t);
    }

    /**
     * @brief Get signal values at multiple time points.
     * @param t Vector of time points to query.
     * @return Vector of interpolated/extrapolated signal values.
     */
    std::vector<BaseType> operator()(const std::vector<double>& t) const
    {
        std::vector<BaseType> xInterp;
        for (size_t i = 0; i < t.size(); i++)
        {
            xInterp.push_back(xAt(t[i]));
        }
        return xInterp;
    }

    /**
     * @brief Get signal derivatives at multiple time points.
     * @param t Vector of time points to query.
     * @return Vector of interpolated/extrapolated derivative values.
     */
    std::vector<TangentType> dot(const std::vector<double>& t) const
    {
        std::vector<TangentType> xDotInterp;
        for (size_t i = 0; i < t.size(); i++)
        {
            xDotInterp.push_back(xDotAt(t[i]));
        }
        return xDotInterp;
    }

    /**
     * @brief Set the interpolation method for this signal.
     * @param method The interpolation method to use.
     */
    void setInterpolationMethod(InterpolationMethod method)
    {
        interpolationMethod = method;
    }

    /**
     * @brief Set the extrapolation method for this signal.
     * @param method The extrapolation method to use.
     */
    void setExtrapolationMethod(ExtrapolationMethod method)
    {
        extrapolationMethod = method;
    }

    /**
     * @brief Set the derivative computation method for this signal.
     * @param method The derivative method to use.
     */
    void setDerivativeMethod(DerivativeMethod method)
    {
        derivativeMethod = method;
    }

    /**
     * @brief Reset the signal to initial state, clearing all history.
     */
    void reset()
    {
        t_    = -1.0;
        x_    = BaseSignalSpec::ZeroType();
        xdot_ = TangentSignalSpec::ZeroType();
        signalHistory_.clear();
        needsSort_ = false;
    }

    /**
     * @brief Update the signal with a new value at a given time, computing derivative automatically.
     * @param _t The time of the new data point.
     * @param _x The value at time _t.
     * @param insertHistory Whether to store this data point in history for interpolation.
     * @return true if update was successful, false otherwise.
     */
    bool update(const double& _t, const BaseType& _x, bool insertHistory = false)
    {
        TangentType _xdot;
        if (!calculateDerivative(_t, _x, _xdot))
        {
            return false;
        }
        return update(_t, _x, _xdot, insertHistory);
    }

    /**
     * @brief Update the signal with a new value and derivative at a given time.
     * @param _t The time of the new data point.
     * @param _x The value at time _t.
     * @param _xdot The derivative at time _t.
     * @param insertHistory Whether to store this data point in history for interpolation.
     * @return true if update was successful, false otherwise.
     */
    bool update(const double& _t, const BaseType& _x, const TangentType& _xdot, bool insertHistory = false)
    {
        t_    = _t;
        x_    = _x;
        xdot_ = _xdot;
        if (insertHistory)
        {
            if (insertIntoHistory(_t, _x, _xdot))
            {
                if (needsSort_)
                {
                    std::sort(signalHistory_.begin(), signalHistory_.end(), SignalDPComparator);
                    needsSort_ = false;
                }
                return true;
            }
            else
            {
                return false;
            }
        }
        return true;
    }

    /**
     * @brief Update the signal with a history of values, computing derivatives automatically.
     * @param _tHistory Vector of time points.
     * @param _xHistory Vector of values corresponding to each time point.
     * @return true if update was successful, false otherwise.
     */
    bool update(const std::vector<double>& _tHistory, const std::vector<BaseType>& _xHistory)
    {
        size_t nTH = _tHistory.size();
        size_t nXH = _xHistory.size();
        if (nTH != nXH)
        {
            return false;
        }
        for (size_t i = 0; i < nXH; i++)
        {
            TangentType _xdot;
            if (!calculateDerivative(_tHistory[i], _xHistory[i], _xdot))
            {
                return false;
            }
            if (!insertIntoHistory(_tHistory[i], _xHistory[i], _xdot))
            {
                return false;
            }
        }
        if (needsSort_)
        {
            std::sort(signalHistory_.begin(), signalHistory_.end(), SignalDPComparator);
            needsSort_ = false;
        }
        return true;
    }

    /**
     * @brief Update the signal with a history of values and derivatives.
     * @param _tHistory Vector of time points.
     * @param _xHistory Vector of values corresponding to each time point.
     * @param _xdotHistory Vector of derivatives corresponding to each time point.
     * @return true if update was successful, false if vector sizes don't match.
     */
    bool update(const std::vector<double>&      _tHistory,
                const std::vector<BaseType>&    _xHistory,
                const std::vector<TangentType>& _xdotHistory)
    {
        size_t nTH  = _tHistory.size();
        size_t nXH  = _xHistory.size();
        size_t nXDH = _xdotHistory.size();
        if (nXH != nXDH || nXDH != nTH)
        {
            return false;
        }
        for (size_t i = 0; i < nXH; i++)
        {
            if (!insertIntoHistory(_tHistory[i], _xHistory[i], _xdotHistory[i]))
            {
                return false;
            }
        }
        if (needsSort_)
        {
            std::sort(signalHistory_.begin(), signalHistory_.end(), SignalDPComparator);
            needsSort_ = false;
        }
        return true;
    }

    static inline BaseType baseZero()
    {
        return BaseSignalSpec::ZeroType();
    }

    static inline TangentType tangentZero()
    {
        return TangentSignalSpec::ZeroType();
    }

    static inline T baseNorm(const BaseType& x)
    {
        return BaseSignalSpec::Norm(x);
    }

    static inline T tangentNorm(const TangentType& x)
    {
        return TangentSignalSpec::Norm(x);
    }

private:
    double                t_;
    BaseType              x_;
    TangentType           xdot_;
    std::vector<SignalDP> signalHistory_;
    bool                  needsSort_;

    bool insertIntoHistory(const double& _t, const BaseType& _x, const TangentType& _xdot)
    {
        if (_t > t_)
        {
            t_    = _t;
            x_    = _x;
            xdot_ = _xdot;
        }
        if (signalHistory_.size() > 0)
        {
            double mostRecentTime = signalHistory_[signalHistory_.size() - 1].t;
            if (_t == mostRecentTime)
            {
                return false;
            }
            else if (_t < mostRecentTime)
            {
                needsSort_ = true;
            }
        }
        signalHistory_.push_back({_t, _x, _xdot});
        return true;
    }

    BaseType xAt(const double& t) const
    {
        if (t == t_)
        {
            return x_;
        }
        else if (signalHistory_.size() > 0)
        {
            int idx = getInterpIndex(t);
            if (idx < 0 || idx >= static_cast<int>(signalHistory_.size()))
            {
                switch (extrapolationMethod)
                {
                case ExtrapolationMethod::NANS:
                    return BaseSignalSpec::NansType();
                    break;
                case ExtrapolationMethod::CLOSEST:
                    if (idx < 0)
                    {
                        return signalHistory_[0].x;
                    }
                    else
                    {
                        return signalHistory_[signalHistory_.size() - 1].x;
                    }
                    break;
                case ExtrapolationMethod::ZEROS:
                default:
                    return BaseSignalSpec::ZeroType();
                    break;
                }
            }
            else
            {
                BaseType    y = xAtIdx(idx);
                TangentType dy;
                switch (interpolationMethod)
                {
                case InterpolationMethod::ZERO_ORDER_HOLD: {
                    dy = TangentSignalSpec::ZeroType();
                    break;
                }
                case InterpolationMethod::LINEAR: {
                    double   t1 = tAtIdx(idx);
                    double   t2 = tAtIdx(idx + 1);
                    BaseType y1 = xAtIdx(idx);
                    BaseType y2 = xAtIdx(idx + 1);
                    dy          = (t - t1) / (t2 - t1) * (y2 - y1);
                    break;
                }
                case InterpolationMethod::CUBIC_SPLINE: {
                    double   t0 = tAtIdx(idx - 1);
                    double   t1 = tAtIdx(idx);
                    double   t2 = tAtIdx(idx + 1);
                    double   t3 = tAtIdx(idx + 2);
                    BaseType y0 = xAtIdx(idx - 1);
                    BaseType y1 = xAtIdx(idx);
                    BaseType y2 = xAtIdx(idx + 1);
                    BaseType y3 = xAtIdx(idx + 2);
                    dy =
                        (t - t1) / (t2 - t1) *
                        ((y2 - y1) + (t2 - t) / (2. * (t2 - t1) * (t2 - t1)) *
                                         (((t2 - t) * (t2 * (y1 - y0) + t0 * (y2 - y1) - t1 * (y2 - y0))) / (t1 - t0) +
                                          ((t - t1) * (t3 * (y2 - y1) + t2 * (y3 - y1) - t1 * (y3 - y2))) / (t3 - t2)));
                    break;
                }
                }
                return y + dy;
            }
        }
        else
        {
            return BaseSignalSpec::ZeroType();
        }
    }

    TangentType xDotAt(const double& t) const
    {
        if (t == t_)
        {
            return xdot_;
        }
        else if (signalHistory_.size() > 0)
        {
            int idx = getInterpIndex(t);
            if (idx < 0 || idx >= static_cast<int>(signalHistory_.size()))
            {
                switch (extrapolationMethod)
                {
                case ExtrapolationMethod::NANS:
                    return TangentSignalSpec::NansType();
                    break;
                case ExtrapolationMethod::CLOSEST:
                    if (idx < 0)
                    {
                        return signalHistory_[0].xdot;
                    }
                    else
                    {
                        return signalHistory_[signalHistory_.size() - 1].xdot;
                    }
                    break;
                case ExtrapolationMethod::ZEROS:
                default:
                    return TangentSignalSpec::ZeroType();
                    break;
                }
            }
            else
            {
                TangentType y = xDotAtIdx(idx);
                TangentType dy;
                switch (interpolationMethod)
                {
                case InterpolationMethod::ZERO_ORDER_HOLD: {
                    dy = TangentSignalSpec::ZeroType();
                    break;
                }
                case InterpolationMethod::LINEAR: {
                    double      t1 = tAtIdx(idx);
                    double      t2 = tAtIdx(idx + 1);
                    TangentType y1 = xDotAtIdx(idx);
                    TangentType y2 = xDotAtIdx(idx + 1);
                    dy             = (t - t1) / (t2 - t1) * (y2 - y1);
                    break;
                }
                case InterpolationMethod::CUBIC_SPLINE: {
                    double      t0 = tAtIdx(idx - 1);
                    double      t1 = tAtIdx(idx);
                    double      t2 = tAtIdx(idx + 1);
                    double      t3 = tAtIdx(idx + 2);
                    TangentType y0 = xDotAtIdx(idx - 1);
                    TangentType y1 = xDotAtIdx(idx);
                    TangentType y2 = xDotAtIdx(idx + 1);
                    TangentType y3 = xDotAtIdx(idx + 2);
                    dy =
                        (t - t1) / (t2 - t1) *
                        ((y2 - y1) + (t2 - t) / (2. * (t2 - t1) * (t2 - t1)) *
                                         (((t2 - t) * (t2 * (y1 - y0) + t0 * (y2 - y1) - t1 * (y2 - y0))) / (t1 - t0) +
                                          ((t - t1) * (t3 * (y2 - y1) + t2 * (y3 - y1) - t1 * (y3 - y2))) / (t3 - t2)));
                    break;
                }
                }
                return y + dy;
            }
        }
        else
        {
            return TangentSignalSpec::ZeroType();
        }
    }

    // Implementation note: signal history size must > 0
    int getInterpIndex(const double& t) const
    {
        if (t < signalHistory_[0].t)
        {
            return -1;
        }
        else if (t > signalHistory_[signalHistory_.size() - 1].t)
        {
            return static_cast<int>(signalHistory_.size());
        }
        else
        {
            for (size_t i = 0; i < signalHistory_.size(); i++)
            {
                double t_i   = signalHistory_[i].t;
                double t_ip1 = tAtIdx(i + 1);
                if (t_i <= t && t_ip1 > t)
                {
                    return static_cast<int>(i);
                }
            }
            return -1;
        }
    }

    double tAtIdx(const int& idx) const
    {
        if (idx < 0)
        {
            return signalHistory_[0].t + static_cast<double>(idx);
        }
        else if (idx >= signalHistory_.size())
        {
            return signalHistory_[signalHistory_.size() - 1].t +
                   static_cast<double>(idx - static_cast<int>(signalHistory_.size() - 1));
        }
        else
        {
            return signalHistory_[idx].t;
        }
    }

    BaseType xAtIdx(const int& idx) const
    {
        if (idx < 0)
        {
            return signalHistory_[0].x;
        }
        else if (idx >= signalHistory_.size())
        {
            return signalHistory_[signalHistory_.size() - 1].x;
        }
        else
        {
            return signalHistory_[idx].x;
        }
    }

    TangentType xDotAtIdx(const int& idx) const
    {
        if (idx < 0)
        {
            return signalHistory_[0].xdot;
        }
        else if (idx >= signalHistory_.size())
        {
            return signalHistory_[signalHistory_.size() - 1].xdot;
        }
        else
        {
            return signalHistory_[idx].xdot;
        }
    }

    // Implementation note: should always be followed by a call to update()
    bool calculateDerivative(const double& _t, const BaseType& _x, TangentType& _xdot)
    {
        const static double sigma = 0.05;
        if (_t <= t_)
        {
            return false;
        }
        if (t_ >= 0)
        {
            const double      dt = _t - t_;
            const TangentType dx = _x - x_;
            switch (derivativeMethod)
            {
            case DerivativeMethod::DIRTY:
                _xdot = (2. * sigma - dt) / (2. * sigma + dt) * xdot_ + 2. / (2. * sigma + dt) * dx;
                break;
            case DerivativeMethod::FINITE_DIFF:
                _xdot = dx / dt;
                break;
            }
        }
        else
        {
            _xdot = TangentSignalSpec::ZeroType();
        }
        return true;
    }
};

/**
 * @brief Add a tangent signal to a base signal.
 * @param l The left-hand side signal (base type).
 * @param r The right-hand side signal (tangent type) to add.
 * @return A new signal with r added to l's values and derivatives.
 */
template<typename T, typename BaseSignalSpec, typename TangentSignalSpec>
Signal<T, BaseSignalSpec, TangentSignalSpec> operator+(const Signal<T, BaseSignalSpec, TangentSignalSpec>&    l,
                                                       const Signal<T, TangentSignalSpec, TangentSignalSpec>& r)
{
    Signal<T, BaseSignalSpec, TangentSignalSpec> lpr = l;
    lpr.x_ += r(l.t());
    lpr.xdot_ += r.dot(l.t());
    for (auto& signalDP : lpr.signalHistory_)
    {
        signalDP.x += r(signalDP.t);
        signalDP.xdot += r.dot(signalDP.t);
    }
    return lpr;
}

/**
 * @brief Subtract one signal from another, producing a tangent signal.
 * @param l The left-hand side signal.
 * @param r The right-hand side signal to subtract.
 * @return A new tangent signal representing the difference l - r.
 */
template<typename T, typename BaseSignalSpec, typename TangentSignalSpec>
Signal<T, TangentSignalSpec, TangentSignalSpec> operator-(const Signal<T, BaseSignalSpec, TangentSignalSpec>& l,
                                                          const Signal<T, BaseSignalSpec, TangentSignalSpec>& r)
{
    Signal<T, TangentSignalSpec, TangentSignalSpec> lmr;
    lmr.interpolationMethod = l.interpolationMethod;
    lmr.extrapolationMethod = l.extrapolationMethod;
    lmr.derivativeMethod    = l.derivativeMethod;
    lmr.needsSort_          = l.needsSort_;
    lmr.t_                  = l.t();
    lmr.x_                  = l() - r(l.t());
    lmr.xdot_               = l.dot() - r.dot(l.t());
    std::vector<double>                           tHistory;
    std::vector<typename TangentSignalSpec::Type> xHistory;
    std::vector<typename TangentSignalSpec::Type> xdotHistory;
    for (auto& signalDP : l.signalHistory_)
    {
        tHistory.push_back(signalDP.t);
        xHistory.push_back(signalDP.x - r(signalDP.t));
        xdotHistory.push_back(signalDP.xdot - r.dot(signalDP.t));
    }
    lmr.update(tHistory, xHistory, xdotHistory);
    return lmr;
}

/**
 * @brief Multiply a signal by a scalar from the left.
 * @param l The scalar multiplier.
 * @param r The signal to multiply.
 * @return A new signal with all values and derivatives scaled by l.
 */
template<typename T, typename BaseSignalSpec, typename TangentSignalSpec>
Signal<T, BaseSignalSpec, TangentSignalSpec> operator*(const double&                                       l,
                                                       const Signal<T, BaseSignalSpec, TangentSignalSpec>& r)
{
    Signal<T, BaseSignalSpec, TangentSignalSpec> lr = r;
    lr.x_ *= l;
    lr.xdot_ *= l;
    for (auto& signalDP : lr.signalHistory_)
    {
        signalDP.x *= l;
        signalDP.xdot *= l;
    }
    return lr;
}

/**
 * @brief Multiply a signal by a scalar from the right.
 * @param l The signal to multiply.
 * @param r The scalar multiplier.
 * @return A new signal with all values and derivatives scaled by r.
 */
template<typename T, typename BaseSignalSpec, typename TangentSignalSpec>
Signal<T, BaseSignalSpec, TangentSignalSpec> operator*(const Signal<T, BaseSignalSpec, TangentSignalSpec>& l,
                                                       const double&                                       r)
{
    Signal<T, BaseSignalSpec, TangentSignalSpec> lr = l;
    lr.x_ *= r;
    lr.xdot_ *= r;
    for (auto& signalDP : lr.signalHistory_)
    {
        signalDP.x *= r;
        signalDP.xdot *= r;
    }
    return lr;
}

/**
 * @brief Type specification for scalar-valued signals.
 * @tparam T Scalar type (e.g., double, float).
 */
template<typename T>
struct ScalarSignalSpec
{
    using Type = T;
    /**
     * @brief Returns zero value for the scalar type.
     */
    static Type ZeroType()
    {
        return (T)0.0;
    }
    /**
     * @brief Returns NaN value for the scalar type.
     */
    static Type NansType()
    {
        return (T)1. / 0.;
    }
    static T Norm(const Type& a)
    {
        return a;
    }
};

/**
 * @brief Type specification for vector-valued signals.
 * @tparam T Scalar element type (e.g., double, float).
 * @tparam d Dimension of the vector.
 */
template<typename T, size_t d>
struct VectorSignalSpec
{
    using Type = Matrix<T, d, 1>;
    /**
     * @brief Returns zero vector.
     */
    static Type ZeroType()
    {
        return Type::Zero();
    }
    /**
     * @brief Returns vector filled with NaN values.
     */
    static Type NansType()
    {
        return Type::Constant(std::numeric_limits<T>::quiet_NaN());
    }
    static T Norm(const Type& a)
    {
        return a.norm();
    }
};

/**
 * @brief Type specification for manifold-valued signals (e.g., SO2, SO3, SE2, SE3).
 * @tparam ManifoldType The manifold type (e.g., SO2<double>, SE3<float>).
 */
template<typename T, typename ManifoldType>
struct ManifoldSignalSpec
{
    using Type = ManifoldType;
    /**
     * @brief Returns identity element of the manifold.
     */
    static Type ZeroType()
    {
        return Type::identity();
    }
    /**
     * @brief Returns manifold element with NaN values.
     */
    static Type NansType()
    {
        return Type::nans();
    }
    static T Norm(const Type& a)
    {
        return ManifoldType::Log(a).norm();
    }
};

template<typename T>
using ScalarSignal = Signal<T, ScalarSignalSpec<T>, ScalarSignalSpec<T>>;

template<typename T>
inline std::ostream& operator<<(std::ostream& os, const ScalarSignal<T>& x)
{
    os << "ScalarSignal at t=" << x.t() << ": " << x();
    return os;
}

template<typename T, size_t d>
using VectorSignal = Signal<T, VectorSignalSpec<T, d>, VectorSignalSpec<T, d>>;

template<typename T, size_t d>
inline std::ostream& operator<<(std::ostream& os, const VectorSignal<T, d>& x)
{
    os << "VectorSignal at t=" << x.t() << ": " << x().transpose();
    return os;
}

template<typename T, typename ManifoldType, size_t d>
using ManifoldSignal = Signal<T, ManifoldSignalSpec<T, ManifoldType>, VectorSignalSpec<T, d>>;

template<typename T, typename ManifoldType, size_t d>
inline std::ostream& operator<<(std::ostream& os, const ManifoldSignal<T, ManifoldType, d>& x)
{
    os << "ManifoldSignal at t=" << x.t() << ": " << x();
    return os;
}

#define MAKE_VECTOR_SIGNAL(Dimension)                                                                                  \
    template<typename T>                                                                                               \
    using Vector##Dimension##Signal = VectorSignal<T, Dimension>;                                                      \
    typedef Vector##Dimension##Signal<double> Vector##Dimension##dSignal;

#define MAKE_MANIF_SIGNAL(Manif, Dimension)                                                                            \
    template<typename T>                                                                                               \
    using Manif##Signal = ManifoldSignal<T, Manif<T>, Dimension>;                                                      \
    typedef Manif##Signal<double> Manif##dSignal;

typedef ScalarSignal<double> ScalardSignal;
MAKE_VECTOR_SIGNAL(1)
MAKE_VECTOR_SIGNAL(2)
MAKE_VECTOR_SIGNAL(3)
MAKE_VECTOR_SIGNAL(4)
MAKE_VECTOR_SIGNAL(5)
MAKE_VECTOR_SIGNAL(6)
MAKE_VECTOR_SIGNAL(7)
MAKE_VECTOR_SIGNAL(8)
MAKE_VECTOR_SIGNAL(9)
MAKE_VECTOR_SIGNAL(10)
MAKE_MANIF_SIGNAL(SO2, 1)
MAKE_MANIF_SIGNAL(SO3, 3)
MAKE_MANIF_SIGNAL(SE2, 3)
MAKE_MANIF_SIGNAL(SE3, 6)
