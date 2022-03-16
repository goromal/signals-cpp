#include <boost/test/unit_test.hpp>
#include <Eigen/Core>
#include "signals/Signals.h"

using namespace Eigen;

using Vector1d = Matrix<double, 1, 1>;

BOOST_AUTO_TEST_SUITE(TestIntegration)

BOOST_AUTO_TEST_CASE(TestEulerIntegrator)
{
    IntegrateEuler integrateEuler;

    const double dt = 0.01;
    double       t  = 0.;

    Vector1dSignal v_ref, v_int;
    SO3dSignal     q_ref, q_int;

    while (t <= M_PI)
    {
        v_ref.update(t, Vector1d(sin(t)), true);
        q_ref.update(t, SO3d::fromAxisAngle(Vector3d(0., 1., 0.), t), true);
        t += dt;
    }

    Vector1dSignal v_ref_dot = v_ref.dotSignal();
    Vector3dSignal q_ref_dot = q_ref.dotSignal();

    v_int.update(0., Vector1d(0.));
    q_int.update(0., SO3d::identity());

    integrateEuler(v_int, v_ref_dot, t, dt);
    integrateEuler(q_int, q_ref_dot, t, dt);

    // TODO
}

BOOST_AUTO_TEST_CASE(TestTrapezoidalIntegrator)
{
    IntegrateTrapezoidal integrateTrapezoidal;
    const double         dt = 0.01;
    // TODO
}

BOOST_AUTO_TEST_SUITE_END()
