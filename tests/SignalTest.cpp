#include <boost/test/unit_test.hpp>
#include <Eigen/Core>
#include "signals/Signals.h"

using namespace Eigen;

BOOST_AUTO_TEST_SUITE(TestSignal)

BOOST_AUTO_TEST_CASE(TestLinearInterpolation)
{
    // TODO
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
