#include <boost/test/unit_test.hpp>
#include <Eigen/Core>
#include "signals/Signals.h"

using namespace Eigen;

using Vector1d = Matrix<double, 1, 1>;

BOOST_AUTO_TEST_SUITE(TestSignal)

BOOST_AUTO_TEST_CASE(TestLinearInterpolation)
{
    ScalardSignal v;
    v.setInterpolationMethod(InterpolationMethod::LINEAR);
    v.setExtrapolationMethod(ExtrapolationMethod::ZEROS);

    BOOST_CHECK(v.update(0., 2., 1., true));
    BOOST_CHECK(v.update(1., 3., 1., true));
    BOOST_CHECK(v.update(2., 4., 1., true));

    BOOST_CHECK_CLOSE(v(1.5), 3.5, 1e-8);
    BOOST_CHECK_CLOSE(v(0.), 2., 1e-8);
    BOOST_CHECK_CLOSE(v(2.01), 0., 1e-8);
    BOOST_CHECK_CLOSE(v.dot(1.5), 1., 1e-8);
    BOOST_CHECK_CLOSE(v.dot(2.01), 0., 1e-8);

    v.setExtrapolationMethod(ExtrapolationMethod::CLOSEST);

    BOOST_CHECK_CLOSE(v(2.01), 4., 1e-8);
    BOOST_CHECK_CLOSE(v.dot(2.01), 1., 1e-8);

    SO3dSignal q;
    q.setInterpolationMethod(InterpolationMethod::LINEAR);
    q.setExtrapolationMethod(ExtrapolationMethod::ZEROS);

    BOOST_CHECK(q.update(0., SO3d::fromEuler(0., 0., 0.), true));
    BOOST_CHECK(q.update(1., SO3d::fromEuler(0., 0., M_PI), true));

    SO3d q_half = SO3d::fromEuler(0., 0., M_PI / 2.);

    BOOST_CHECK_CLOSE(q(0.5).w(), q_half.w(), 1e-8);
    BOOST_CHECK_CLOSE(q(0.5).x(), q_half.x(), 1e-8);
    BOOST_CHECK_CLOSE(q(0.5).y(), q_half.y(), 1e-8);
    BOOST_CHECK_CLOSE(q(0.5).z(), q_half.z(), 1e-8);
    BOOST_CHECK_CLOSE(q(-0.5).w(), 1., 1e-8);
    BOOST_CHECK_CLOSE(q(-0.5).x(), 0., 1e-8);
    BOOST_CHECK_CLOSE(q(-0.5).y(), 0., 1e-8);
    BOOST_CHECK_CLOSE(q(-0.5).z(), 0., 1e-8);
}

BOOST_AUTO_TEST_CASE(TestDirtyDerivative)
{
    Vector2dSignal v;
    double         a = 1.75;
    double         b = -0.5;
    for (size_t i = 0; i < 100; i++)
    {
        v.update(static_cast<double>(i), Vector2d(i * a, i * b), true);
    }
    BOOST_CHECK_CLOSE(v.dot(33.3)(0), a, 1.);
    BOOST_CHECK_CLOSE(v.dot(66.6)(1), b, 1.);

    Vector1dSignal u;
    double         dt = 0.05;
    for (size_t i = 0; i < 1000; i++)
    {
        double t = i * dt;
        u.update(t, Vector1d(sin(t)), true);
    }
    BOOST_CHECK_CLOSE(u.dot()(0), cos(u.t()), 5.);
}

BOOST_AUTO_TEST_CASE(TestSetEquality)
{
    Vector1dSignal v1;
    v1.update(0., Vector1d(4.), true);
    v1.update(10., Vector1d(-10.), true);

    Vector1dSignal v2 = v1;

    BOOST_CHECK_CLOSE(v1(5.)(0), v2(5.)(0), 1e-8);
    BOOST_CHECK_CLOSE(v1.dot(5.)(0), v2.dot(5.)(0), 1e-8);
    BOOST_CHECK_CLOSE(v1(20.)(0), v2(20.)(0), 1e-8);
    BOOST_CHECK_CLOSE(v1.dot(20.)(0), v2.dot(20.)(0), 1e-8);
}

BOOST_AUTO_TEST_CASE(TestScaledVectorPlusMinus)
{
    Vector1dSignal u, v1, v2, v3;
    double         dt = 0.01;
    for (size_t i = 0; i < 1000; i++)
    {
        double t = i * dt;
        v1.update(t, Vector1d(sin(t)), true);
        v2.update(t, Vector1d(cos(t)), true);
        v3.update(t, Vector1d(sin(2. * t)), true);
    }
    u = v1 + 2. * v2 - 3. * v3;

    double t_test = 5.0004;
    BOOST_CHECK_CLOSE(u(t_test)(0), sin(t_test) + 2. * cos(t_test) - 3. * sin(2. * t_test), 1e-3);
}

BOOST_AUTO_TEST_CASE(TestScaledManifoldPlusMinus)
{
    SO3dSignal u, v1, v2, v3;
    double     dt = 0.01;
    for (size_t i = 0; i < 1000; i++)
    {
        double t = i * dt;
        v1.update(t, SO3d::fromEuler(0.1, -0.1, 0), true);
        v2.update(t, SO3d::fromEuler(0.2, -0.3, 0), true);
        v3.update(t, SO3d::fromEuler(-0.1, 0.24, 0.1), true);
    }
    u = v1 + (2. * v2 - v3);

    double t_test = 5.0004;
    SO3d   q_est  = u(t_test);
    SO3d   q_tru =
        SO3d::fromEuler(0.1, -0.1, 0) + (2. * SO3d::fromEuler(0.2, -0.3, 0) - SO3d::fromEuler(-0.1, 0.24, 0.1));

    BOOST_CHECK_CLOSE(q_est.w(), q_tru.w(), 1e-3);
    BOOST_CHECK_CLOSE(q_est.x(), q_tru.x(), 1e-3);
    BOOST_CHECK_CLOSE(q_est.y(), q_tru.y(), 1e-3);
    BOOST_CHECK_CLOSE(q_est.z(), q_tru.z(), 1e-3);
}

BOOST_AUTO_TEST_SUITE_END()
