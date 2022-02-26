#pragma once
#include <Eigen/Core>
// TODO implement vector sorting, or maybe use different data structure?
// https://en.cppreference.com/w/cpp/algorithm/sort

using namespace Eigen;

// TODO INTEGRATION SUPPORT (in separate file, for vector signals only)

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

    enum Interpolation
    {
        ZERO_ORDER_HOLD,
        LINEAR,
        CUBIC_SPLINE
    } interpolationMethod;
    enum Extrapolation
    {
        NANS,
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
        extrapolationMethod = Extrapolation::NANS;
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

    //   template<typename T2, typename ST2, size_t tDim2>
    //   friend Signal<T2, ST2, tDim2> operator+(const Signal<T2, ST2, tDim2>&                                 l,
    //                                           const Signal<T2, Signal<T2, ST2, tDim2>::TangentType, tDim2>& r)
    //   {
    //       // TODO don't forget histories! use l's time (what if l needs sort? This should be friend
    //   }

    //  template<typename T2, typename ST2, size_t tDim2>
    //  friend Signal<T2, Signal<T2, ST2, tDim2>::TangentType, tDim2> operator-(const Signal<T2, ST2, tDim2>& l,
    //                                                                          const Signal<T2, ST2, tDim2>& r)
    //  {
    //      // TODO
    //  }

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

private:
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
                // TODO switch for extrapolation method
            }
            else
            {
                // TODO
                // switch for interpolation method
                // get needed ts and ys for base and interp
                // call one of the overloaded (TangentType and SignalType) interp methods
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
                // TODO switch for extrapolation method
            }
            else
            {
                // TODO
                // switch for interpolation method
                // get needed ts and ys for base and interp
                // call one of the overloaded (TangentType and SignalType) interp methods
            }
        }
        else
        {
            return TangentType::Zero();
        }
    }

    // Implementation note: signal history size must > 0
    int getInterpIndex(const double& t)
    {
        if (needsSort_)
        {
            // TODO sort vector by time
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
        // TODO
    }

    SignalType xAtIdx(const int& idx)
    {
        // TODO
    }

    TangentType xDotAtIdx(const int& idx)
    {
        // TODO
    }

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
            _xdot = TangentType::Zero();
        }
        return true;
    }
};
