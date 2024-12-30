#include <boost/test/unit_test.hpp>
#include <Eigen/Core>
#include "signals/Function.h"

using namespace Eigen;

BOOST_AUTO_TEST_SUITE(TestFunction) // ^^^^ TODO

BOOST_AUTO_TEST_CASE(TestEulerIntegrator)
{
    const double dt = 0.001;
    double       t  = 0.;
    double       tf = 5. * M_PI / 4.;

    ScalardSignal v_ref, v_int;

    while (t <= tf)
    {
        BOOST_CHECK(v_ref.update(t, sin(t), cos(t), true));
        t += dt;
    }

    ScalardSignal v_ref_dot = v_ref.dotSignal();

    BOOST_CHECK(v_int.update(0., 0.));

    BOOST_CHECK(EulerIntegrator::integrate(v_int, v_ref_dot, t, dt, true));

    BOOST_CHECK_CLOSE(v_int(0.), v_ref(0.), 1.);
    BOOST_CHECK_CLOSE(v_int(M_PI / 4.), v_ref(M_PI / 4.), 1.);
    BOOST_CHECK_CLOSE(v_int(M_PI / 2.), v_ref(M_PI / 2.), 1.);
    BOOST_CHECK_CLOSE(v_int(3. * M_PI / 4.), v_ref(3. * M_PI / 4.), 1.);
}

BOOST_AUTO_TEST_CASE(TestTrapezoidalIntegrator)
{
    const double dt = 0.001;
    double       t  = 0.;
    double       tf = 5. * M_PI / 4.;

    ScalardSignal v_ref, v_int;

    while (t <= tf)
    {
        BOOST_CHECK(v_ref.update(t, sin(t), cos(t), true));
        t += dt;
    }

    ScalardSignal v_ref_dot = v_ref.dotSignal();

    BOOST_CHECK(v_int.update(0., 0.));

    BOOST_CHECK(TrapezoidalIntegrator::integrate(v_int, v_ref_dot, t, dt, true));

    BOOST_CHECK_CLOSE(v_int(0.), v_ref(0.), 1.);
    BOOST_CHECK_CLOSE(v_int(M_PI / 4.), v_ref(M_PI / 4.), 1.);
    BOOST_CHECK_CLOSE(v_int(M_PI / 2.), v_ref(M_PI / 2.), 1.);
    BOOST_CHECK_CLOSE(v_int(3. * M_PI / 4.), v_ref(3. * M_PI / 4.), 1.);
}

BOOST_AUTO_TEST_SUITE_END()
