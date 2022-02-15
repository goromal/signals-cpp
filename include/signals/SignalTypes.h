#pragma once

#include <Eigen/Core>
#include <SO3.h>
#include <SE3.h>

using namespace Eigen;

template<typename T>
class Signal
{
public:
    using ScalarType = T;

    size_t dimensionality() { return dim_; }
    size_t tangentDimensionality() { return tan_dim_; }
    virtual void setZero() = 0;

protected:
    size_t dim_;
    size_t tan_dim_;
};

template<typename T, int dim>
class VectorSignal final : public Signal<T>
{
public:
    using TangentType = VectorSignal;

    VectorSignal() : v_{Vector<T, dim>::Zero()}, dim_{dim}, tan_dim_{dim} {}
    VectorSignal(const Vector<T, dim> &v) : v_{v}, dim_{dim}, tan_dim_{dim} {}
    void setZero() override { v_.setZero(); }
    inline const Vector<T, dim> get() const { return v_; }
    inline Vector<T, dim> get() { return v_; }
    VectorSignal& operator= (const VectorSignal &other)
    {
        v_ = other.get();
        return *this;
    }
    VectorSignal operator+ (const TangentType &other)
    {
        return VectorSignal(get() + other.get());
    }
    VectorSignal& operator+= (const TangentType &other)
    {
        v_ += other.get();
        return *this;
    }
    TangentType operator- (const VectorSignal &other)
    {
        return VectorSignal(get() - other.get());
    }
    VectorSignal& operator-= (const VectorSignal &other)
    {
        v_ -= other.get();
        return *this;
    }

private:
    Vector<T, dim> v_;
};

template<typename T, typename Base, int tangentDim>
class ManifoldSignal : public Signal<T>
{
public:
    using TangentType = VectorSignal<T,tangentDim>;

    ManifoldSignal() : X_{Base<T>::Identity()}, tan_dim_{tangentDim} {}
    ManifoldSignal(const Base<T> &X) : X_{X}, tan_dim_{tangentDim} {}
    void setZero() override { X_ = Base<T>::Identity(); }
    inline const Base<T> get() const { return X_; }
    inline Base<T> get() { return X_; } 
    ManifoldSignal& operator= (const ManifoldSignal &other) 
    { 
        X_ = other.get();
        return *this; 
    }
    ManifoldSignal operator* (const ManifoldSignal &other) const
    {
        return ManifoldSignal(get() * other.get());
    }
    ManifoldSignal& operator*= (const ManifoldSignal &other)
    {
        X_ *= other.get();
        return *this;
    }
    ManifoldSignal operator+ (const TangentType &other) const
    {
        return ManifoldSignal(get() + other.get());
    }
    ManifoldSignal& operator+= (const TangentType &other)
    {
        X_ += other.get();
        return *this;
    }
    TangentType operator- (const ManifoldSignal& other) const
    {
        return TangentType(get() - other.get());
    }

protected:
    Base<T> X_;
};

template<typename T>
class SO3Signal final : public ManifoldSignal<T,SO3,3>
{
public:
    SO3Signal() : dim_{4}, ManifoldSignal<T,SO3,3>() {}
    SO3Signal(const SO3<T> &X) : dim_{4}, ManifoldSignal<T,SO3,3>(X) {}
    VectorSignal<T,3> operator* (const VectorSignal<T,3> &v)
    {
        return VectorSignal<T,3>(get() * v.get());
    }
};

template<typename T>
class SE3Signal final : public ManifoldSignal<T,SE3,6>
{
public:
    SE3Signal() : dim_{7}, ManifoldSignal<T,SE3,6>() {}
    SE3Signal(const SE3<T> &X) : dim_{7}, ManifoldSignal<T,SE3,6>(X) {}
    VectorSignal<T,3> operator* (const VectorSignal<T,3> &v)
    {
        return VectorSignal<T,3>(get() * v.get());
    }
};

