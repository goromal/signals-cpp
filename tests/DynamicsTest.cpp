#include <boost/test/unit_test.hpp>
#include <Eigen/Core>
#include "signals/Signals.h"

using namespace Eigen;

BOOST_AUTO_TEST_SUITE(TestDynamics)

BOOST_AUTO_TEST_CASE(TestTransSim)
{
    Translational1DOFSystemd sys, sys2;
    ScalardSignal            u;

    RigidBodyParams1D p;
    p.m = 1.0;
    p.g = 0.0;
    sys.setParams(p);
    sys2.setParams(p);

    double t   = 0;
    double x_d = 1;
    double kp  = 1;
    double kd  = 0.5;

    BOOST_CHECK(sys.x.update(t, ScalardState::identity()));
    BOOST_CHECK(sys2.x.update(t, ScalardState::identity()));

    const double       dt        = 0.01;
    const unsigned int num_iters = 1e4;

    for (unsigned int i = 0; i < num_iters; i++)
    {
        BOOST_CHECK(u.update(t, kp * (x_d - sys.x().pose) - kd * sys.x().twist, true));
        t += dt;
        BOOST_CHECK(sys.simulate<EulerIntegrator>(u, t));
    }

    BOOST_CHECK_CLOSE(sys.x().pose, 1., 1.);
    BOOST_CHECK_LT(sys.x().twist, 1e-8);

    BOOST_CHECK(sys2.simulate<EulerIntegrator>(u, t, dt));

    BOOST_CHECK_CLOSE(sys2.x().pose, 1., 1.);
    BOOST_CHECK_LT(sys2.x().twist, 1e-8);
}

BOOST_AUTO_TEST_CASE(TestSO2Sim)
{
    Rotational1DOFSystemd sys;
    Vector1dSignal        u;
    
    RigidBodyParams2D p;
    p.m = 1.0;
    p.J = 1.0;
    p.g = Vector2d::Zero();
    sys.setParams(p);

    double t       = 0;
    double angle_d = 1.0;
    SO2d   x_d     = SO2d::fromAngle(angle_d);
    double kp      = 1;
    double kd      = 0.5;

    BOOST_CHECK(sys.x.update(t, SO2dState::identity()));

    const double       dt        = 0.01;
    const unsigned int num_iters = 1e4;

    for (unsigned int i = 0; i < num_iters; i++)
    {
        BOOST_CHECK(u.update(t, kp * (x_d - sys.x().pose) - kd * sys.x().twist));
        t += dt;
        BOOST_CHECK(sys.simulate<EulerIntegrator>(u, t));
    }

    BOOST_CHECK_CLOSE(sys.x().pose.w(), x_d.w(), 1.);
    BOOST_CHECK_CLOSE(sys.x().pose.x(), x_d.x(), 1.);
    BOOST_CHECK_LT(sys.x().twist.x(), 1e-8);
}

BOOST_AUTO_TEST_CASE(TestSO3Sim)
{
    Rotational3DOFSystemd sys;
    Vector3dSignal        u;
    
    RigidBodyParams3D p;
    p.m = 1.0;
    p.J = Matrix3d::Identity();
    p.g = Vector3d::Zero();
    sys.setParams(p);

    double t       = 0;
    double roll_d  = 1.0;
    double pitch_d = -2.0;
    double yaw_d   = 0.5;
    SO3d   x_d     = SO3d::fromEuler(roll_d, pitch_d, yaw_d);
    double kp      = 1;
    double kd      = 0.5;

    BOOST_CHECK(sys.x.update(t, SO3dState::identity()));

    const double       dt        = 0.01;
    const unsigned int num_iters = 1e4;

    for (unsigned int i = 0; i < num_iters; i++)
    {
        BOOST_CHECK(u.update(t, kp * (x_d - sys.x().pose) - kd * sys.x().twist));
        t += dt;
        BOOST_CHECK(sys.simulate<EulerIntegrator>(u, t));
    }

    BOOST_CHECK_CLOSE(sys.x().pose.w(), x_d.w(), 1.);
    BOOST_CHECK_CLOSE(sys.x().pose.x(), x_d.x(), 1.);
    BOOST_CHECK_CLOSE(sys.x().pose.y(), x_d.y(), 1.);
    BOOST_CHECK_CLOSE(sys.x().pose.z(), x_d.z(), 1.);
    BOOST_CHECK_LT(sys.x().twist.x(), 1e-8);
    BOOST_CHECK_LT(sys.x().twist.y(), 1e-8);
    BOOST_CHECK_LT(sys.x().twist.z(), 1e-8);
}

BOOST_AUTO_TEST_CASE(TestSE3Sim)
{
    RigidBody6DOFSystemd sys;
    Vector6dSignal       u;
    
    RigidBodyParams3D p;
    p.m = 1.0;
    p.J = Matrix3d::Identity();
    p.g = Vector3d::Zero();
    sys.setParams(p);

    double t       = 0;
    double tx_d    = 0.5;
    double ty_d    = -3.0;
    double tz_d    = 2.0;
    double roll_d  = 1.0;
    double pitch_d = -2.0;
    double yaw_d   = 0.5;
    SE3d   x_d     = SE3d::fromVecAndQuat(Vector3d(tx_d, ty_d, tz_d),
                                          SO3d::fromEuler(roll_d, pitch_d, yaw_d));
    double kp      = 1;
    double kd      = 0.5;

    BOOST_CHECK(sys.x.update(t, SE3dState::identity()));

    const double       dt        = 0.01;
    const unsigned int num_iters = 1e4;

    for (unsigned int i = 0; i < num_iters; i++)
    {
        BOOST_CHECK(u.update(t, kp * (x_d - sys.x().pose) - kd * sys.x().twist));
        t += dt;
        BOOST_CHECK(sys.simulate<EulerIntegrator>(u, t));
    }

    BOOST_CHECK_CLOSE(sys.x().pose.t().x(), x_d.t().x(), 1.);
    BOOST_CHECK_CLOSE(sys.x().pose.t().y(), x_d.t().y(), 1.);
    BOOST_CHECK_CLOSE(sys.x().pose.t().z(), x_d.t().z(), 1.);
    BOOST_CHECK_CLOSE(sys.x().pose.q().w(), x_d.q().w(), 1.);
    BOOST_CHECK_CLOSE(sys.x().pose.q().x(), x_d.q().x(), 1.);
    BOOST_CHECK_CLOSE(sys.x().pose.q().y(), x_d.q().y(), 1.);
    BOOST_CHECK_CLOSE(sys.x().pose.q().z(), x_d.q().z(), 1.);
    BOOST_CHECK_LT(sys.x().twist(0), 1e-8);
    BOOST_CHECK_LT(sys.x().twist(1), 1e-8);
    BOOST_CHECK_LT(sys.x().twist(2), 1e-8);
    BOOST_CHECK_LT(sys.x().twist(3), 1e-8);
    BOOST_CHECK_LT(sys.x().twist(4), 1e-8);
    BOOST_CHECK_LT(sys.x().twist(5), 1e-8);
}

BOOST_AUTO_TEST_SUITE_END()
