#include <Eigen/Dense>

using namespace Eigen;

Vector3d shp2pck(double ETsecs, double angle0, double angle1, double angle2, double spin2, double spindot2) {

  // return Vector3d = RA, Dev, W0

//function [RA, Dec, W0] = shape2pck(ETsecs, angle0, angle1, angle2, spin2, spindot2)
    
    // Equatorial vs. Ecliptic. Shape uses this value (IAU2006)
#define EPS  84381.406/3600*M_PI/180

    Matrix3d rot2eq, r, rt, rot_eq2eqx;
  double epoch, phi, theta, psi, cphi, ctheta, cpsi, sphi, stheta, spsi;
  double spin_rate, spin_accel;
  double alpha_z, delta_z, W, W0;

  Vector3d x_ecr, x_ec, x_eq, z_ec, z_eq, x_eqx_eq, y_eqx_eq, z_eqx_eq, x_bennu_eqx;
  Vector3d xhat(1,0,0), zhat(0,0,1);


    rot2eq << 1, 0 ,0,
          0, cos(EPS), sin(EPS),
          0, -sin(EPS), cos(EPS);

    // Time and rotation rate
    epoch = ETsecs; // Sec past J2000. Should this routine do the ET-UT J2000?
    phi = angle0* M_PI / 180.;
    theta = angle1* M_PI / 180.;
    psi = angle2* M_PI / 180.;
    ctheta = cos(theta);
    stheta = sin(theta);
    cphi = cos(phi);
    sphi = sin(phi);
    cpsi = cos(psi);
    spsi = sin(psi);
    spin_rate = spin2 / 86400; // in deg/sec
    spin_accel = spindot2 / 86400 / 86400;

    // ZXZ euler angles to rotation matrix (paraphrasing shape code)
    r << cpsi*cphi-ctheta*sphi*spsi, -spsi*cphi-ctheta*sphi*cpsi,  stheta*sphi,
         cpsi*sphi+ctheta*cphi*spsi, -spsi*sphi+ctheta*cphi*cpsi, -stheta*cphi,
         spsi*stheta,                 cpsi*stheta,                 ctheta;
    rt = r.transpose();
    
    
// X-axis of Bennu in equatorial frame
x_ecr = xhat.transpose() * rt;

x_ec = x_ecr/x_ecr.norm();

x_eq = x_ec.transpose() * rot2eq.transpose();

// Z-axis of Bennu in equatorial frame
z_ec = zhat.transpose()* rt;

z_eq = z_ec.transpose()* rot2eq.transpose();

alpha_z = atan2(z_eq(1),z_eq(0)) * 180 / M_PI;
delta_z = asin(z_eq(2)) * 180 / M_PI;



// Bennu equinox frame: X = equinox, Z = spin pole, Y = Z x X
x_eqx_eq = zhat.cross(z_eq); x_eqx_eq = x_eqx_eq/x_eqx_eq.norm();
y_eqx_eq = z_eq.cross(x_eqx_eq);
z_eqx_eq = z_eq;
rot_eq2eqx << x_eqx_eq,y_eqx_eq,z_eqx_eq;

// Bennu x-axis in equinox frame
x_bennu_eqx = rot_eq2eqx * x_eq;

// Angle from equinox to x-axis
W = atan2(x_bennu_eqx(1), x_bennu_eqx(0)) * 180/M_PI;

W0 = fmod(W - spin_rate * epoch - 0.5 * spin_accel * epoch * epoch, 360);
return(Vector3d(alpha_z, delta_z, W0));

}

