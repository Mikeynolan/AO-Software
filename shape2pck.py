from math import pi, sin, cos, atan2, asin
import numpy as np

# arrays are transposed wrt matlab

def shape2pck(ETsecs, angle0, angle1, angle2, spin2, spindot2):
    # Equatorial vs. Ecliptic. Shape uses this value (IAU2006)
    eps = 84381.406/3600*pi/180
    rot2eq = np.array([
      [1,  0,        0       ],
      [0,  cos(eps), sin(eps)],
      [0, -sin(eps), cos(eps)]
    ])
    #rot2eq = rot2eq.transpose()
    # Time and rotation rate
    epoch = ETsecs # Sec past J2000. Should this routine do the ET-UT J2000?
    phi = angle0 * pi / 180
    theta = angle1 * pi / 180
    psi = angle2 * pi / 180
    ctheta = cos(theta)
    stheta = sin(theta)
    cphi = cos(phi)
    sphi = sin(phi)
    cpsi = cos(psi)
    spsi = sin(psi)
    spin_rate = spin2 / 86400; # in deg/sec
    spin_accel = spindot2 / 86400 / 86400
    xhat = np.array([1,0,0])
    zhat = np.array([0,0,1])
    # ZXZ euler angles to rotation matrix (paraphrasing shape code)
    r = np.array([
      [cpsi*cphi-ctheta*sphi*spsi, -spsi*cphi-ctheta*sphi*cpsi,  stheta*sphi],
      [cpsi*sphi+ctheta*cphi*spsi, -spsi*sphi+ctheta*cphi*cpsi, -stheta*cphi],
      [spsi*stheta,                 cpsi*stheta,                 ctheta     ]
    ])
    rt = r.transpose()
    # X-axis of Bennu in equatorial frame
    x_ecr = np.matmul(xhat, rt)
    
    x_ec = x_ecr/np.linalg.norm(x_ecr)
    x_eq = np.matmul(x_ec, rot2eq)
    
    # Z-axis of Bennu in equatorial frame
    z_ec = np.matmul(zhat, rt)
    z_eq = np.matmul(z_ec, rot2eq)
    
    alpha_z = atan2(z_eq[1],z_eq[0]) * 180 / pi
    delta_z = asin(z_eq[2]) * 180 / pi
    
    # Bennu equinox frame: X = equinox, Z = spin pole, Y = Z x X
    x_eqx_eq = np.cross(zhat, z_eq)
    x_eqx_eq = x_eqx_eq/np.linalg.norm(x_eqx_eq)
    y_eqx_eq = np.cross(z_eq, x_eqx_eq)
    z_eqx_eq = z_eq
    rot_eq2eqx = np.array([x_eqx_eq,y_eqx_eq,z_eqx_eq])
    
    # Bennu x-axis in equinox frame
    x_bennu_eqx = np.matmul(rot_eq2eqx, x_eq)
    
    # Angle from equinox to x-axis
    W = atan2(x_bennu_eqx[1], x_bennu_eqx[0]) * 180/pi
    
    W0 = (W - spin_rate * epoch - 0.5 * spin_accel * epoch * epoch) % 360
    return(alpha_z, delta_z, W0)
