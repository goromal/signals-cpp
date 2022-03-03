#pragma once
#include <algorithm>
#include <Eigen/Core>
#include <SO3.h>
#include <SE3.h>

using namespace Eigen;

template<typename T, typename ST, size_t tDim>
class Signal
{
public:
    typedef T                  ScalarType;
    typedef ST                 SignalType;
    typedef Matrix<T, tDim, 1> TangentType;
    const static size_t        tangentDim = tDim;

    struct SignalDP
    {
        double      t;
        SignalType  x;
        TangentType xdot;
    };
    struct
    {
        bool operator()(SignalDP a, SignalDP b) const
        {
            return a.t < b.t;
        }
    } SignalDPComparator;

    enum Interpolation
    {
        ZERO_ORDER_HOLD,
        LINEAR,
        CUBIC_SPLINE
    } interpolationMethod;
    enum Extrapolation
    {
        NANS,
        ZEROS,
        CLOSEST
    } extrapolationMethod;
    enum Derivative
    {
        DIRTY,
        FINITE_DIFF
    } derivativeMethod;

    Signal()
    {
        interpolationMethod = Interpolation::LINEAR;
        extrapolationMethod = Extrapolation::ZEROS;
        derivativeMethod    = Derivative::DIRTY;
        reset();
    }

    double t()
    {
        return t_;
    }
    SignalType operator()()
    {
        return x_;
    }
    TangentType dot()
    {
        return xdot_;
    }
    SignalType operator()(const double& t)
    {
        xAt(t);
    }
    TangentType dot(const double& t)
    {
        return xDotAt(t);
    }
    std::vector<SignalType> operator()(const std::vector<double>& t)
    {
        std::vector<SignalType> xInterp;
        for (size_t i = 0; i < t.size(); i++)
        {
            xInterp.push_back(xAt(t[i]));
        }
        return xInterp;
    }
    std::vector<TangentType> dot(const std::vector<double>& t)
    {
        std::vector<TangentType> xDotInterp;
        for (size_t i = 0; i < t.size(); i++)
        {
            xDotInterp.push_back(xDotAt(t[i]));
        }
        return xDotInterp;
    }

    void setInterpolationMethod(Interpolation method)
    {
        interpolationMethod = method;
    }
    void setExtrapolationMethod(Extrapolation method)
    {
        extrapolationMethod = method;
    }
    void setDerivativeMethod(Derivative method)
    {
        derivativeMethod = method;
    }

    void reset()
    {
        t_ = -1.0;
        x_ = SignalType();
        xdot_.setZero();
        signalHistory_.clear();
        needsSort_ = false;
    }

    bool update(const double& _t, const SignalType& _x, bool insertHistory = false)
    {
        TangentType _xdot;
        if (!calculateDerivative(_t, _x, _xdot))
        {
            return false;
        }
        return update(_t, _x, _xdot, insertHistory);
    }

    bool update(const double& _t, const SignalType& _x, const TangentType& _xdot, bool insertHistory = false)
    {
        t_    = _t;
        x_    = _x;
        xdot_ = _xdot;
        if (insertHistory)
        {
            return insertIntoHistory(_t, _x, _xdot);
        }
        return true;
    }

    bool update(const std::vector<double>& _tHistory, const std::vector<SignalType>& _xHistory)
    {
        size_t nTH = _tHistory.size();
        size_t nXH = _xHistory.size();
        if (nTH != nXH)
        {
            return false;
        }
        std::vector<TangentType> _xdotHistory;
        for (size_t i = 0; i < nTH; i++)
        {
            TangentType _xdot;
            if (!calculateDerivative(_tHistory[i], _xHistory[i], _xdot))
            {
                return false;
            }
            _xdotHistory.push_back(_xdot);
        }
        return update(_tHistory, _xHistory, _xdotHistory);
    }

    bool update(const std::vector<double>&      _tHistory,
                const std::vector<SignalType>&  _xHistory,
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
        return true;
    }

    Signal& operator=(const Signal& other)
    {
        interpolationMethod = other.interpolationMethod;
        extrapolationMethod = other.extrapolationMethod;
        derivativeMethod    = other.derivativeMethod;
        t_                  = other.t_;
        x_                  = other.x_;
        xdot_               = other.xdot_;
        signalHistory_      = other.signalHistory_;
        needsSort_          = other.needsSort_;
        return *this;
    }

    template<typename T2, typename ST2, size_t tDim2>
    friend Signal<T2, ST2, tDim2> operator+(const Signal<T2, ST2, tDim2>&                                 l,
                                            const Signal<T2, Signal<T2, ST2, tDim2>::TangentType, tDim2>& r)
    {
        Signal<T2, ST2, tDim2> lpr = l;
        lpr.x += r.x(l.t());
        lpr.xdot += r.xdot(l.t());
        for (auto signalDP : lpr.signalHistory_)
        {
            signalDP.x += r.x(signalDP.t);
            signalDP.xdot += r.xdot(signalDP.t);
        }
        return lpr;
    }

    template<typename T2, typename ST2, size_t tDim2>
    friend Signal<T2, Signal<T2, ST2, tDim2>::TangentType, tDim2> operator-(const Signal<T2, ST2, tDim2>& l,
                                                                            const Signal<T2, ST2, tDim2>& r)
    {
        Signal<T2, Signal<T2, ST2, tDim2>::TangentType, tDim2> lmr;
        lmr.interpolationMethod = l.interpolationMethod;
        lmr.extrapolationMethod = l.extrapolationMethod;
        lmr.derivativeMethod    = l.derivativeMethod;
        lmr.needsSort_          = l.needsSort_;
        lmr.t_                  = l.t();
        lmr.x_                  = l.x() - r.x(l.t());
        lmr.xdot_               = l.dot() - r.dot(l.t());
        std::vector<double>                                                              tHistory;
        std::vector<Signal<T2, Signal<T2, ST2, tDim2>::TangentType, tDim2>::SignalType>  xHistory;
        std::vector<Signal<T2, Signal<T2, ST2, tDim2>::TangentType, tDim2>::TangentType> xdotHistory;
        for (auto signalDP : l.signalHistory_)
        {
            tHistory.push_back(signalDP.t);
            xHistory.push_back(signalDP.x - r.x(signalDP.t));
            xdotHistory.push_back(signalDP.xdot - r.xdot(signalDP.t));
        }
        lmr.update(tHistory, xHistory, xdotHistory);
        return lmr;
    }

    template<typename T2, typename ST2, size_t tDim2>
    friend Signal<T2, ST2, tDim2> operator*(const double& l, const Signal<T2, ST2, tDim2>& r)
    {
        Signal<T2, ST2, tDim2> lr = r;
        lr.x_ *= l;
        lr.xdot_ *= l;
        for (auto signalDP : lr.signalHistory_)
        {
            signalDP.x *= l;
            signalDP.xdot *= l;
        }
        return lr;
    }

    template<typename T2, typename ST2, size_t tDim2>
    friend Signal<T2, ST2, tDim2> operator*(const Signal<T2, ST2, tDim2>& l, const double& r)
    {
        Signal<T2, ST2, tDim2> lr = l;
        lr.x_ *= r;
        lr.xdot_ *= r;
        for (auto signalDP : lr.signalHistory_)
        {
            signalDP.x *= r;
            signalDP.xdot *= r;
        }
        return lr;
    }

protected:
    double                t_;
    SignalType            x_;
    TangentType           xdot_;
    std::vector<SignalDP> signalHistory_;
    bool                  needsSort_;

    bool insertIntoHistory(const double& _t, const SignalType& _x, const TangentType& _xdot)
    {
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

    SignalType xAt(const double& t)
    {
        if (signalHistory_.size() > 0)
        {
            int idx = getInterpIndex(t);
            if (idx < 0 || idx > static_cast<int>(signalHistory_.size()))
            {
                switch (extrapolationMethod)
                {
                case Extrapolation::NANS:
                    return getNansSignal();
                    break;
                case Extrapolation::CLOSEST:
                    if (idx < 0)
                    {
                        return signalHistory_[0].x;
                    }
                    else
                    {
                        return signalHistory_[signalHistory_.size() - 1].x;
                    }
                    break;
                case Extrapolation::ZEROS:
                default:
                    return getZeroSignal();
                    break;
                }
            }
            else
            {
                SignalType  y = xAtIdx(idx);
                TangentType dy;
                switch (interpolationMethod)
                {
                case Interpolation::ZERO_ORDER_HOLD:
                {
                    dy = getZeroTangent();
                    break;
                }
                case Interpolation::LINEAR:
                {
                    double     t1 = tAtIdx(idx);
                    double     t2 = tAtIdx(idx + 1);
                    SignalType y1 = xAtIdx(idx);
                    SignalType y2 = xAtIdx(idx + 1);
                    dy            = (t - t1) / (t2 - t1) * (y2 - y1);
                    break;
                }
                case Interpolation::CUBIC_SPLINE:
                {
                    double     t0 = tAtIdx(idx - 1);
                    double     t1 = tAtIdx(idx);
                    double     t2 = tAtIdx(idx + 1);
                    double     t3 = tAtIdx(idx + 2);
                    SignalType y0 = xAtIdx(idx - 1);
                    SignalType y1 = xAtIdx(idx);
                    SignalType y2 = xAtIdx(idx + 1);
                    SignalType y3 = xAtIdx(idx + 2);
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
            return SignalType();
        }
    }

    TangentType xDotAt(const double& t)
    {
        if (signalHistory_.size() > 0)
        {
            int idx = getInterpIndex(t);
            if (idx < 0 || idx > static_cast<int>(signalHistory_.size()))
            {
                switch (extrapolationMethod)
                {
                case Extrapolation::NANS:
                    return getNansTangent();
                    break;
                case Extrapolation::CLOSEST:
                    if (idx < 0)
                    {
                        return signalHistory_[0].xdot;
                    }
                    else
                    {
                        return signalHistory_[signalHistory_.size() - 1].xdot;
                    }
                    break;
                case Extrapolation::ZEROS:
                default:
                    return getZeroTangent();
                    break;
                }
            }
            else
            {
                TangentType y = xdotAtIdx(idx);
                TangentType dy;
                switch (interpolationMethod)
                {
                case Interpolation::ZERO_ORDER_HOLD:
                {
                    dy = getZeroTangent();
                    break;
                }
                case Interpolation::LINEAR:
                {
                    double      t1 = tAtIdx(idx);
                    double      t2 = tAtIdx(idx + 1);
                    TangentType y1 = xdotAtIdx(idx);
                    TangentType y2 = xdotAtIdx(idx + 1);
                    dy             = (t - t1) / (t2 - t1) * (y2 - y1);
                    break;
                }
                case Interpolation::CUBIC_SPLINE:
                {
                    double      t0 = tAtIdx(idx - 1);
                    double      t1 = tAtIdx(idx);
                    double      t2 = tAtIdx(idx + 1);
                    double      t3 = tAtIdx(idx + 2);
                    TangentType y0 = xdotAtIdx(idx - 1);
                    TangentType y1 = xdotAtIdx(idx);
                    TangentType y2 = xdotAtIdx(idx + 1);
                    TangentType y3 = xdotAtIdx(idx + 2);
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
            return getZeroTangent();
        }
    }

    // Implementation note: signal history size must > 0
    int getInterpIndex(const double& t)
    {
        if (needsSort_)
        {
            std::sort(signalHistory_.begin(), signalHistory_.end(), SignalDPComparator);
            needsSort_ = false;
        }
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
        }
    }

    double tAtIdx(const int& idx)
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

    SignalType xAtIdx(const int& idx)
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

    TangentType xDotAtIdx(const int& idx)
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

    virtual SignalType  getZeroSignal()  = 0;
    virtual SignalType  getNansSignal()  = 0;
    virtual TangentType getZeroTangent() = 0;
    virtual TangentType getNansTangent() = 0;

    // Implementation note: should always be followed by a call to update()
    bool calculateDerivative(const double& _t, const SignalType& _x, const TangentType& _xdot)
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
            case Derivative::DIRTY:
                _xdot = (2. * sigma - dt) / (2. * sigma + dt) * xdot_ + 2. / (2. * sigma + dt) * dx;
                break;
            case Derivative::FINITE_DIFF:
                _xdot = dx / dt;
                break;
            }
        }
        else
        {
            _xdot = getZeroTangent();
        }
        return true;
    }
};

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
template<typename T>
using Vector1Signal = VectorSignal<T, 1>;
template<typename T>
using Vector2Signal = VectorSignal<T, 2>;
template<typename T>
using Vector3Signal = VectorSignal<T, 3>;
template<typename T>
using Vector4Signal = VectorSignal<T, 4>;
template<typename T>
using Vector5Signal = VectorSignal<T, 5>;
template<typename T>
using Vector6Signal = VectorSignal<T, 6>;
template<typename T>
using Vector7Signal = VectorSignal<T, 7>;
template<typename T>
using Vector8Signal = VectorSignal<T, 8>;
template<typename T>
using Vector9Signal = VectorSignal<T, 9>;
template<typename T>
using Vector10Signal = VectorSignal<T, 10>;

typedef SO3Signal<double>      SO3dSignal;
typedef SE3Signal<double>      SE3dSignal;
typedef Vector1Signal<double>  Vector1dSignal;
typedef Vector2Signal<double>  Vector2dSignal;
typedef Vector3Signal<double>  Vector3dSignal;
typedef Vector4Signal<double>  Vector4dSignal;
typedef Vector5Signal<double>  Vector5dSignal;
typedef Vector6Signal<double>  Vector6dSignal;
typedef Vector7Signal<double>  Vector7dSignal;
typedef Vector8Signal<double>  Vector8dSignal;
typedef Vector9Signal<double>  Vector9dSignal;
typedef Vector10Signal<double> Vector10dSignal;
