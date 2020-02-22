from math import pi, sin, cos, atan2, asin, acos
import numpy as np

# some arrays are transposed wrt matlab

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

def pck2shape(RA0, DEC0, W0):
    #
    # This will give values for an epoch of J2000. It could do a more
    # useful epoch, or that could happen elsewhere. Note that J2000 isn't a
    # legal date in shape, which requires integer seconds. Since it's J2000,
    # don't care about W1 and W2
    #I am as literally as possible undoing shape2pck, which comes from
    #Steve Chesley's compute_W0
    eps = 84381.406/3600*pi/180
    rot2eq = np.array([
      [1,  0,        0       ],
      [0,  cos(eps), sin(eps)],
      [0, -sin(eps), cos(eps)]
    ])
    # We're going the other way, so transpose
    rotfeq = rot2eq.transpose()
    a = RA0 * pi / 180.
    d = DEC0 * pi / 180
    w = W0 * pi / 180
    xhat = np.array([1,0,0])
    zhat = np.array([0,0,1])

    z_eq =  np.array([cos(a)*cos(d), sin(a)*cos(d), sin(d)]);
    # Bennu equinox frame: X = equinox, Z = spin pole, Y = Z x X
    x_eqx_eq = np.cross(zhat, z_eq)
    x_eqx_eq = x_eqx_eq/np.linalg.norm(x_eqx_eq)
    y_eqx_eq = np.cross(z_eq, x_eqx_eq)
    z_eqx_eq = z_eq
    rot_eq2eqx = np.array([x_eqx_eq,y_eqx_eq,z_eqx_eq])

    #Bennu x-axis: Xz is 0 in both coordinate systems because it is by
    #definition the cross-product of the Zs

    x_bennu_eqx = np.array([cos(w), sin(w), 0])
    x_eq = np.matmul(rot_eq2eqx.transpose(), x_bennu_eqx)
    y_eq = np.cross(z_eq, x_eq)
    
    z_ec = np.matmul(z_eq, rotfeq)
    x_ec = np.matmul(x_eq, rotfeq)
    y_ec = np.matmul(y_eq, rotfeq)

    theta = acos(z_ec[2])
    psi = atan2(x_ec[2], y_ec[2])
    phi = atan2(z_ec[0], -z_ec[1])

    angle0 = phi * 180 / pi
    angle1 = theta * 180 / pi
    angle2 = psi * 180 / pi

    return(angle0, angle1, angle2)
