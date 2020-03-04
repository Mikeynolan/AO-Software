#!/usr/bin/env python3
"""Convert spin block between shape and pck.

shape2pck converts the spin parameters back and forth from shape to pck
(surprise) formats. It takes an input file either on the command line
or as stdin and outputs the spin state portion of the other format. It
undertands spin accelerations (YORP), but no other changing spin state
parameters (precession in pck, npa, libration, spin impulses in mod files)
It's probably not robust against nonstandard formatting. In particular,
there can be no blank lines in the spin block.'
"""
from math import pi, sin, cos, atan2, asin, acos
import numpy as np
from astropy.time import Time
import re
import datetime
import sys
import argparse

# search for either pck or mod marker.


# \begindata
#    BODY2101955_POLE_RA    = (  85.4567       0.              0. )
#    BODY2101955_POLE_DEC   = ( -60.3574       0.              0. )
#    BODY2101955_PM         = ( 135.8156  2011.14645095 1.815e-06 )
#    BODY2101955_LONG_AXIS  = (   0.                              )
# \begintext
def main():
    """Process pck or mod file.

    Determnine which file type input, parse it, and call the appropriate
    translator. Output to stdout.

    Returns
    -------
    None.

    """
    parser = argparse.ArgumentParser(
        description="convert spin blocks between shape and spice")
    parser.add_argument('infile', nargs='?', type=argparse.FileType('r'),
                        help="File to convert, if none will try stdin",
                        default=sys.stdin)

    parser.add_argument("-e", "--epoch", action="store",
                        help="t0 to use when generating a shape spin block."
                        "Default is 2001-01-00T00:00:00. "
                        "Ignored when coverting to PCK.",
                        default="2000-01-01T00:00:00")
    args = parser.parse_args()

    RE_PCK = re.compile(r'\\begindata')
    RE_MOD = re.compile(r'{SPIN STATE}')
    ftype = 0
    it = args.infile
    for line in it:
        mo = re.search(RE_PCK, line)
        if mo:
            ftype = 1
            break
        mo = re.search(RE_MOD, line)
        if mo:
            ftype = 2
            break

    # THis is a perl-y way to do it.
    if (0 == ftype):
        sys.exit('No spin state found')

    # We just found one or the other of the markers.
    if (1 == ftype):
        # PCK. Parse and send to spck2shape
        # Could have multiple \beginttext and \begindata
        # Could also get spicypy, but that seems like overkill to read
        # three lines
        indata = True  # Did a seek[0] becsue I cant move it
        RA0 = 999
        DEC0 = 999
        W0 = 999
        for line in it:
            vals = line.split()
            if indata and vals:
                if vals[0].endswith("_POLE_RA"):
                    RA0 = float(vals[3])
                elif vals[0].endswith("_POLE_DEC"):
                    DEC0 = float(vals[3])
                elif vals[0].endswith("_PM"):
                    W0 = float(vals[3])
                    W1 = float(vals[4])
                    W2 = float(vals[5])
                elif vals[0].startswith("\\begintext"):
                    indata = False
            elif vals and vals[0].startswith("\\begindata"):
                indata = True
        # Done parsing pck
        it.close()
        if RA0 < 999 and DEC0 < 999 and W0 < 999:
            angles = pck2shape(RA0, DEC0, W0)
            writemod(args.epoch, W1, W2, angles)
        else:
            print("Didn't find spin parameters in PCK file")
            exit(1)
    elif 2 == ftype:
        # mod file.
        # {SPIN STATE}
        # 2001  3 18  0  0  0 {yyyy mo dd hh mm ss of t0}
        #  c    60.0000000000 {angle 0 (deg) lambda=330.000000}
        #  c   179.0000000000 {angle 1 (deg) beta=-89.000000}
        #  f   320.4246664927 {angle 2 (deg)}
        #  c     0.0000000000 {spin 0 (deg/day)}
        #  c     0.0000000000 {spin 1 (deg/day)}
        #  c   712.1444532272 {spin 2 (deg/day) P=12.132370}
        #  c     0.5511373633 {moment of inertia 0}
        #  c     1.2996277328 {moment of inertia 1}
        #  c     1.3736408397 {moment of inertia 2}
        #  c     0.0000000000 {spin 0 dot (deg/day/day)}
        #  c     0.0000000000 {spin 1 dot (deg/day/day)}
        #  c     0.0000000000 {spin 2 dot (deg/day/day)}
        #  c     0.0000000000 {Libration Amplitude (degrees)}
        #  c     0.0000000000 {Libration Frequency (degrees/day)}
        #  c     0.0000000000 {Libration Phase (degrees)}
        #                   0 {number of spin impulses}
        line = next(it)
        vals = line.split()
        epoch = [int(i) for i in vals[0:6]]
        epoch = Time(datetime.datetime(*epoch))
        # PCK needs J2000 TDB.
        daysJ2000 = (epoch - Time(datetime.datetime(2000, 1, 1, 12, 0, 0),
                                  scale='tdb')).jd
        line = next(it)
        vals = line.split()
        angle0 = float(vals[1])
        line = next(it)
        vals = line.split()
        angle1 = float(vals[1])
        line = next(it)
        vals = line.split()
        angle2 = float(vals[1])
        line = next(it)
        vals = line.split()
        spin0 = float(vals[1])
        line = next(it)
        vals = line.split()
        spin1 = float(vals[1])
        line = next(it)
        vals = line.split()
        spin2 = float(vals[1])
        line = next(it)  # Skip moments of inertia
        line = next(it)
        line = next(it)
        line = next(it)
        # ancient files don't have YORP parameters
        if line.index("spin 0 dot") > 0:
            vals = line.split()
            spin0dot = float(vals[1])
            line = next(it)
            vals = line.split()
            spin1dot = float(vals[1])
            line = next(it)
            vals = line.split()
            spin2dot = float(vals[1])
        else:
            spin0dot = 0
            spin1dot = 0
            spin2dot = 0

        # Check for invalid ones
        if 0 != spin0 or 0 != spin1 or 0 != spin0dot or 0 != spin1dot:
            sys.exit("NPA rotation can't be converted to a text pck")
        it.close()
        stuff = shape2pck(daysJ2000, angle0, angle1, angle2, spin2, spin2dot)
        writepck(stuff, spin2, spin2dot)

# some arrays are transposed wrt matlab


def shape2pck(ETdays, angle0, angle1, angle2, spin2, spindot2):
    """Convert shape parameters to pck parameters.

    Take in shape euler angles etc., convert to J2000 coordinates and epoch,
    and produce output RA, Dec, and phase parameters. This was conv erted from
    Steve Chesley's compute_W0 pretty literally.
    """
    # Equatorial vs. Ecliptic. Shape uses this value (IAU2006)
    eps = 84381.406/3600*pi/180
    rot2eq = np.array([
      [1,         0,        0],
      [0,  cos(eps), sin(eps)],
      [0, -sin(eps), cos(eps)]
    ])
    # rot2eq = rot2eq.transpose()
    # Time and rotation rate
    epoch = ETdays  # Days past J2000 TDB.
    phi = angle0 * pi / 180
    theta = angle1 * pi / 180
    psi = angle2 * pi / 180
    ctheta = cos(theta)
    stheta = sin(theta)
    cphi = cos(phi)
    sphi = sin(phi)
    cpsi = cos(psi)
    spsi = sin(psi)
    spin_rate = spin2  # in deg/day
    spin_accel = spindot2
    xhat = np.array([1, 0, 0])
    zhat = np.array([0, 0, 1])
    # ZXZ euler angles to rotation matrix (paraphrasing shape code)
    r = np.array([
      [cpsi*cphi-ctheta*sphi*spsi, -spsi*cphi-ctheta*sphi*cpsi,  stheta*sphi],
      [cpsi*sphi+ctheta*cphi*spsi, -spsi*sphi+ctheta*cphi*cpsi, -stheta*cphi],
      [spsi*stheta,                 cpsi*stheta,                 ctheta]
    ])
    rt = r.transpose()
    # X-axis of Bennu in equatorial frame
    x_ecr = np.matmul(xhat, rt)

    x_ec = x_ecr/np.linalg.norm(x_ecr)
    x_eq = np.matmul(x_ec, rot2eq)

    # Z-axis of Bennu in equatorial frame
    z_ec = np.matmul(zhat, rt)
    z_eq = np.matmul(z_ec, rot2eq)

    alpha_z = atan2(z_eq[1], z_eq[0]) * 180 / pi
    delta_z = asin(z_eq[2]) * 180 / pi

    # Bennu equinox frame: X = equinox, Z = spin pole, Y = Z x X
    x_eqx_eq = np.cross(zhat, z_eq)
    x_eqx_eq = x_eqx_eq/np.linalg.norm(x_eqx_eq)
    y_eqx_eq = np.cross(z_eq, x_eqx_eq)
    z_eqx_eq = z_eq
    rot_eq2eqx = np.array([x_eqx_eq, y_eqx_eq, z_eqx_eq])

    # Bennu x-axis in equinox frame
    x_bennu_eqx = np.matmul(rot_eq2eqx, x_eq)

    # Angle from equinox to x-axis
    W = atan2(x_bennu_eqx[1], x_bennu_eqx[0]) * 180/pi
    W0 = (W - spin_rate * epoch - 0.5 * spin_accel * epoch * epoch) % 360

    return(alpha_z, delta_z, W0)


def pck2shape(RA0, DEC0, W0):
    """Convert pck parameter to shape angles.

    This will give values for an epoch of J2000. That's not legal inb a mod
    file, so need to fix it up later. Since it's J2000,
    don't care about W1 and W2
    I am as literally as possible undoing shape2pck, which comes from
    Steve Chesley's compute_W0
    """
    eps = 84381.406/3600*pi/180
    rot2eq = np.array([
      [1,         0,        0],
      [0,  cos(eps), sin(eps)],
      [0, -sin(eps), cos(eps)]
    ])
    # We're going the other way, so transpose
    rotfeq = rot2eq.transpose()
    a = RA0 * pi / 180.
    d = DEC0 * pi / 180
    w = W0 * pi / 180
    zhat = np.array([0, 0, 1])

    z_eq = np.array([cos(a)*cos(d), sin(a)*cos(d), sin(d)])
    # Bennu equinox frame: X = equinox, Z = spin pole, Y = Z x X
    x_eqx_eq = np.cross(zhat, z_eq)
    x_eqx_eq = x_eqx_eq/np.linalg.norm(x_eqx_eq)
    y_eqx_eq = np.cross(z_eq, x_eqx_eq)
    z_eqx_eq = z_eq
    rot_eq2eqx = np.array([x_eqx_eq, y_eqx_eq, z_eqx_eq])

    # Bennu x-axis: Xz is 0 in both coordinate systems because it is by
    # definition the cross-product of the Zs

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


def writepck(stuff, spin2, spin2dot):
    """Write out the pck block.

    Write out the pck-formatted values. Use lots of digits, because the
    user can always trim but can't add more.
    """
    RA0 = stuff[0]
    DEC0 = stuff[1]
    W0 = stuff[2]
    W2 = spin2dot/2  # squared term of polynomial, not accel

    print("\\begindata")
    print(f'BODY2XXXXXX_POLE_RA   = ( {RA0:14.10f}    0.0    0.0 )')
    print(f'BODY2XXXXXX_POLE_DEC  = ( {DEC0:14.10f}    0.0    0.0 )')
    print(f'BODY2XXXXXX_PM        = ( {W0:14.10f} {spin2:.10f} {W2:7e} )')
    print(f'BODY2XXXXXX_LONG_AXIS = (   0.                         )')
    print("\\begintext")


def writemod(ep, W1, W2, angles):
    """
    Write out the mod file splibn block.

    Parameters
    ----------
    ep : string ISO date
        Date for use in setting the t0 in the spin block. THe angle is rotated
        to this date from J2000 as provided.
    W1 : float
        Spin rate in degrees/day.
    W2 : float
        Spin rate change in deg/day/day. This is the squared term, not the
        acceleration.
    angles : list of floats
        contains the euler angles for J2000.

    Returns
    -------
    None.
    """
    J2000 = Time(datetime.datetime(2000, 1, 1, 12, 0, 0), scale='tdb')
    moddate = Time(ep, format='isot', scale='utc')
    # truncate seconds if provided
    moddate = Time(int(moddate.to_value(format='unix')), format='unix',
                   scale='utc')
    diff = (moddate - J2000).jd
    d = moddate.to_value('datetime')
    w0 = angles[2]
    w = (w0 + W1*diff + W2*diff*diff) % 360
    elon = (angles[0] - 90) % 360
    elat = (90 - angles[1])
    per = 24 / (W1 / 360)
    accel = W2 * 2  # accel, not squared term
    print(f"{{SPIN STATE}}")
    print(f"    {d.year:4d} {d.month:2d} {d.day:2d} "
          f"{d.hour:2d} {d.minute:2d} {d.second:2d} "
          f"{{yyyy mo dd hh mm ss of t0}}")
    print(f" c   {angles[0]:14.10f} {{angle 0 (deg) lambda={elon:10.6f}}}")
    print(f" c   {angles[1]:14.10f} {{angle 1 (deg) beta={elat:10.6f}}}")
    print(f" c   {w:14.10f} {{angle 2 (deg)}}")
    print(f" c     0.0000000000 {{spin 0 (deg/day)}}")
    print(f" c     0.0000000000 {{spin 1 (deg/day)}}")
    print(f" c   {W1:.10f} {{spin 2 (deg/day) P={per:.6f}}}")
    print(f" c     1.0000000000 {{moment of inertia 0}}")
    print(f" c     1.1000000000 {{moment of inertia 1}}")
    print(f" c     1.2000000000 {{moment of inertia 1}}")
    print(f" c     0.0000000000 {{spin 0 dot (deg/day/day)}}")
    print(f" c     0.0000000000 {{spin 1 dot (deg/day/day)}}")
    print(f" c   {accel:14.10f} {{spin 2 dot (deg/day/day)}}")
    print(f" c     0.0000000000 {{Libration Amplitude (degrees)}}")
    print(f" c     0.0000000000 {{Libration Frequency (degrees/day)}}")
    print(f" c     0.0000000000 {{Libration Phase (degrees)}}")
    print(f"                  0 {{number of spin impulses}}")


if __name__ == '__main__':
    main()
