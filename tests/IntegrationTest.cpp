#include <boost/test/unit_test.hpp>
#include <Eigen/Core>
#include "signals/Signals.h"

using namespace Eigen;

BOOST_AUTO_TEST_SUITE(TestIntegration)

BOOST_AUTO_TEST_CASE(TestEulerIntegrator)
{
    IntegrateEuler integrateEuler;

    const double dt = 0.001;
    double       t  = 0.;
    double       tf = 5. * M_PI / 4.;

    ScalardSignal v_ref, v_int;

    while (t <= tf)
    {
        v_ref.update(t, sin(t), cos(t), true);
        t += dt;
    }

    ScalardSignal v_ref_dot = v_ref.dotSignal();

    v_int.update(0., 0.);

    integrateEuler(v_int, v_ref_dot, t, dt, true);

    BOOST_CHECK_CLOSE(v_int(0.), v_ref(0.), 1.);
    BOOST_CHECK_CLOSE(v_int(M_PI / 4.), v_ref(M_PI / 4.), 1.);
    BOOST_CHECK_CLOSE(v_int(M_PI / 2.), v_ref(M_PI / 2.), 1.);
    BOOST_CHECK_CLOSE(v_int(3. * M_PI / 4.), v_ref(3. * M_PI / 4.), 1.);
}

BOOST_AUTO_TEST_CASE(TestTrapezoidalIntegrator)
{
    IntegrateTrapezoidal integrateTrapezoidal;

    const double dt = 0.001;
    double       t  = 0.;
    double       tf = 5. * M_PI / 4.;

    ScalardSignal v_ref, v_int;

    while (t <= tf)
    {
        v_ref.update(t, sin(t), cos(t), true);
        t += dt;
    }

    ScalardSignal v_ref_dot = v_ref.dotSignal();

    v_int.update(0., 0.);

    integrateTrapezoidal(v_int, v_ref_dot, t, dt, true);

    BOOST_CHECK_CLOSE(v_int(0.), v_ref(0.), 1.);
    BOOST_CHECK_CLOSE(v_int(M_PI / 4.), v_ref(M_PI / 4.), 1.);
    BOOST_CHECK_CLOSE(v_int(M_PI / 2.), v_ref(M_PI / 2.), 1.);
    BOOST_CHECK_CLOSE(v_int(3. * M_PI / 4.), v_ref(3. * M_PI / 4.), 1.);
}

BOOST_AUTO_TEST_SUITE_END()
