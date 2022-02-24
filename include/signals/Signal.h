#pragma once
#include <Eigen/Core>
#include <SO3.h>
#include <SE3.h>
// TODO implement vector sorting, or maybe use different data structure?
// https://en.cppreference.com/w/cpp/algorithm/sort

using namespace Eigen;

// TODO INTEGRATION SUPPORT (in separate file, for vector signals only)

template<typename T, typename ST, size_t tDim>
class Signal
{
public:
   typedef T                ScalarType;
   typedef ST               SignalType;
   typedef Matrix<T,tDim,1> TangentType; 
   const static size_t tangentDim = tDim;
   
   struct SignalDP
   {
      double t;
      SignalType x;
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
      ZEROS,
      CLOSEST
   } extrapolationMethod;
   enum Derivative
   {
      DIRTY,
      FORWARD_DIFF,
      CENTRAL_DIFF
   } derivativeMethod;
   
   Signal() 
   {
      interpolationMethod = Interpolation::LINEAR;
      extrapolationMethod = Extrapolation::NANS;
      derivativeMethod = Derivative::DIRTY;
      reset();
   }
   
   double t() { return t; }
   SignalType operator()() { return x; }
   TangentType dot() { return xdot; }
   SignalType operator()(const double &t) { xAt(t); }
   TangentType dot(const double &t) { return xDotAt(t); }
   std::vector<SignalType> operator()(const std::vector<double> &t)
   {
      std::vector<SignalType> xInterp;
      for (size_t i = 0; i < t.size(); i++)
      {
         xInterp.push_back(xAt(t[i]));
      }
      return xInterp;
   }
   std::vector<TangentType> dot(const std::vector<double> &t)
   {
      std::vector<TangentType> xDotInterp;
      for (size_t i = 0; i < t.size(); i++)
      {
         xDotInterp.push_back(xDotAt(t[i]));
      }
      return xDotInterp;
   }
   
   void setInterpolationMethod(Interpolation method) { interpolationMethod = method; }
   void setExtrapolationMethod(Extrapolation method) { extrapolationMethod = method; }
   void setDerivativeMethod(Derivative method) { derivativeMethod = method; }
   
   void reset()
   {
      t = -1.0;
      x = SignalType();
      xdot.setZero();
      tHistory.clear();
      xHistory.clear();
      xdotHistory.clear();
      needsSort = false;
   }
   
   bool update(const double &_t, const SignalType &_x, bool insertHistory=false)
   {
      TangentType _xdot;
      if (!calculateDerivative(_t, _x, _xdot))
      {
         return false;
      } 
      return update(_t, _x, _xdot, insertHistory);
   }
   
   bool update(const double &_t, const SignalType &_x, const TangentType &_xdot, bool insertHistory=false)
   {
      t = _t;
      x = _x;
      xdot = _xdot;
      if (insertHistory)
      {
         return insertIntoHistory(_t, _x, _xdot);
      }
      return true; 
   }
   
   bool update(const std::vector<double> &_tHistory, const std::vector<SignalType> &_xHistory)
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
   
   bool update(const std::vector<double> &_tHistory, const std::vector<SignalType> &_xHistory,
               const std::vector<TangentType> &_xdotHistory)
   {
      size_t nTH = _tHistory.size();
      size_t nXH = _xHistory.size();
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
   
   Signal& operator= (const Signal &other)
   {
      x = other.x;
      xdot = other.xdot;
      // TODO history and other private types
      return *this;
   }

private:
   double t;
   SignalType x;
   TangentType xdot;
   std::vector<SignalDP> signalHistory;
   bool needsSort;
   
   bool insertIntoHistory(const double &_t, const SignalType &_x, const TangentType &_xdot)
   {
      if (signalHistory.size() > 0)
      {
         double mostRecentTime = signalHistory[signalHistory.size()-1].t;
         if (_t == mostRecentTime)
         {
            return false;
         }
         else if (_t < mostRecentTime)
         {
            needsSort = true;
         }
      }
      signalHistory.push_back({_t, _x, _xdot});
      return true;
   }
   
   SignalType xAt(const double &t)
   {
      // TODO
   }
   
   TangentType xDotAt(const double &t)
   {
      // TODO
   }
   
   bool calculateDerivative(const double &_t, const SignalType &_x, const TangentType &_xdot)
   {
      // TODO
   } 
};

template<typename T, typename ST, size_t tDim>
Signal<T,ST,tDim> operator+(const Signal<T,ST,tDim> &l, const Signal<T,Signal<T,ST,tDim>::TangentType,tDim> &r)
{
   // TODO don't forget histories! use l's time (what if l needs sort? This should be friend
}

template<typename T, typename ST, size_t tDim>
Signal<T,Signal<T,ST,tDim>::TangentType,tDim> operator-(const Signal<T,ST,tDim> &l, const Signal<T,ST,tDim> &r)
{
   // TODO
}

template<typename T, typename ST, size_t tDim>
Signal<T,ST,tDim> operator*(const double &l, const Signal<T,ST,tDim> &r)
{
   // TODO
}

template<typename T, typename ST, size_t tDim>
Signal<T,ST,tDim> operator*(const Signal<T,ST,tDim> &l, const double &r)
{
   // TODO
}
