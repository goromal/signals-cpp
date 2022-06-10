#include <boost/test/unit_test.hpp>
#include <Eigen/Core>
#include "signals/Signals.h"

#include <iostream> // ----

using namespace Eigen;

BOOST_AUTO_TEST_SUITE(TestDynamics)

BOOST_AUTO_TEST_CASE(TestTransSim)
{
    Translational1DOFSystemd sys;
    ScalardSignal            u;

    double t   = 0;
    double x_d = 1;
    double kp  = 1;
    double kd  = 0.5;

    BOOST_CHECK(sys.x.update(t, ScalardState::identity()));

    const double       dt        = 0.01;
    const unsigned int num_iters = 1e4;

    for (unsigned int i = 0; i < num_iters; i++)
    {
        BOOST_CHECK(u.update(t, kp * (x_d - sys.x().pose) - kd * sys.x().twist));
        t += dt;
        BOOST_CHECK(sys.simulate<IntegrateEuler>(u, t));
    }

    BOOST_CHECK_CLOSE(sys.x().pose, 1., 1.);
    BOOST_CHECK_LT(sys.x().twist, 1e-8);

    // TODO lump integration with u and sys2
}

BOOST_AUTO_TEST_CASE(TestSO3Sim)
{
    // TODO
}

BOOST_AUTO_TEST_CASE(TestSE3Sim)
{
    // SE3RigidBodySimulateEuler simulator;
    // TODO
}

BOOST_AUTO_TEST_SUITE_END()
