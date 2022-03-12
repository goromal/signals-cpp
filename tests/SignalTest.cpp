#include <boost/test/unit_test.hpp>
#include <Eigen/Core>
#include "signals/Signals.h"

using namespace Eigen;

BOOST_AUTO_TEST_SUITE(TestSignal)

BOOST_AUTO_TEST_CASE(TestLinearInterpolation)
{
    Vector1dSignal v;
    v.setInterpolationMethod(Vector1dSignal::Interpolation::LINEAR);
    v.setExtrapolationMethod(Vector1dSignal::Extrapolation::ZEROS);

    BOOST_CHECK(v.update(0., Vector1d(2.), Vector1d(1.), true));
    BOOST_CHECK(v.update(1., Vector1d(3.), Vector1d(1.), true));
    BOOST_CHECK(v.update(2., Vector1d(4.), Vector1d(1.), true));

    BOOST_CHECK_CLOSE(v(1.5), 3.5, 1e-8);
    BOOST_CHECK_CLOSE(v(0.), 2., 1e-8);
    BOOST_CHECK_CLOSE(v(2.01), 0., 1e-8);
    BOOST_CHECK_CLOSE(v.dot(1.5), 1., 1e-8);
    BOOST_CHECK_CLOSE(v.dot(2.01), 0., 1e-8);

    SO3dSignal q;
    q.setInterpolationMethod(SO3dSignal::Interpolation::LINEAR);
    q.setExtrapolationMethod(SO3dSignal::Extrapolation::ZEROS);

    BOOST_CHECK(q.update(0., SO3d::fromEuler(0., 0., 0.), true));
    BOOST_CHECK(q.update(1., SO3d::fromEuler(0., 0., M_PI), true));

    SO3d q_half = SO3d::fromEuler(0., 0., M_PI/2.);

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
    // TODO
}

BOOST_AUTO_TEST_CASE(TestSetEquality)
{
    // TODO
}

BOOST_AUTO_TEST_CASE(TestScaledPlus)
{
    // TODO
}

BOOST_AUTO_TEST_CASE(TestScaledMinus)
{
    // TODO
}

BOOST_AUTO_TEST_CASE(TestNothing)
{
    Vector3dSignal v;
    std::cout << v.update(0, Vector3d(3., 2., 1.)) << std::endl;
    std::cout << v().transpose() << std::endl;

    Vector3dSignal u = 3.0 * v;
    std::cout << u().transpose() << std::endl;
}

BOOST_AUTO_TEST_SUITE_END()
