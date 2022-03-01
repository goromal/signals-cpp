#pragma once

#include <Eigen/Core>
#include <SO3.h>
#include <SE3.h>

#include "Signal.h"

using namespace Eigen;

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
