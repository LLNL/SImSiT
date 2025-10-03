"""
Utility functions for telescope observation simulation and image generation.
"""

import contextlib
import sys

import astropy.units as u
import galsim
import numpy as np
import ssapy
from astropy.coordinates import GCRS, AltAz
from astropy.table import Table
from astropy.time import Time

from .tracker import transform_wcs


def draw_stars(stars, *, t0, exptime, wcs0, tracker, psf, image, nsplit=2):
    """
    Parameters
    ----------
    stars : astropy.table.Table
        Columns for
            ra, dec : float
                degrees
            nphot : float
    exptime : float
        Seconds
    wcs0 : galsim.GSFitsWCS
    tracker : xfiles.Tracker
    t0 : astropy.time.Time
    psf : galsim.GSObject
    image : galsim.Image
    nsplit : int

    Returns
    -------
    stars : astropy.table.Table
        Columns for
            ra, dec : float
                degrees
            nphot : float
            x, y : float
                Image coordinates
    """
    # Figure out coordinates first so we can use vectorization.  The goal is to
    # determine the radec at t=t0 that yields the same image position as
    # actual radec at later times.
    dt = exptime/nsplit * u.s
    ras = []
    decs = []
    xs = []
    ys = []
    for isplit in range(nsplit+1):
        t = t0 + dt*isplit
        boresight = tracker.get_boresight(t)
        rot_sky_pos = tracker.get_rot_sky_pos(t)
        wcs = transform_wcs(wcs0, boresight, rot_sky_pos)
        x, y = wcs.radecToxy(stars['ra'], stars['dec'], units='degrees')
        ra, dec = wcs0.xyToradec(x, y, units='radians')
        ras.append(ra)
        decs.append(dec)
        xs.append(x)
        ys.append(y)

    # Add midpoint xy to star table
    tmid = t0 + 0.5*exptime*u.s
    boresight = tracker.get_boresight(tmid)
    rot_sky_pos = tracker.get_rot_sky_pos(tmid)
    wcs = transform_wcs(wcs0, boresight, rot_sky_pos)
    stars['x'], stars['y'] = wcs.radecToxy(
        stars['ra'], stars['dec'], units='degrees'
    )

    for istar, star in enumerate(stars):
        for isplit in range(1, nsplit+1):
            obj = galsim.Convolve(
                getLineGSObject(
                    galsim.CelestialCoord(
                        ras[isplit-1][istar]*galsim.radians,
                        decs[isplit-1][istar]*galsim.radians
                    ),
                    galsim.CelestialCoord(
                        ras[isplit][istar]*galsim.radians,
                        decs[isplit][istar]*galsim.radians
                    )
                ),
                psf
            )
            xy = galsim.PositionD(
                0.5*(xs[isplit-1][istar] + xs[isplit][istar]),
                0.5*(ys[isplit-1][istar] + ys[isplit][istar])
            )
            local_wcs = wcs0.local(xy)
            stamp = (obj*star['nphot']/nsplit).drawImage(
                wcs=local_wcs,
                center=xy,
                method='phot'
            )
            bounds = stamp.bounds & image.bounds
            if bounds.area() > 0:
                image[bounds] += stamp[bounds]

    # Compute xy positions in FITS coords which start with (1.0, 1.0) in the
    # center of the lower left pixel.
    stars['x_FITS'] = stars['x'] - image.xmin + 1
    stars['y_FITS'] = stars['y'] - image.ymin + 1

    return stars


def draw_sat(
    orbit, *,
    t0, exptime,
    wcs0, tracker,
    psf, image,
    observer, nphot, nsplit=10,
    propagator=None
):
    """
    Parameters
    ----------
    orbit : ssapy.Orbit
    t0 : astropy.time.Time
    exptime : float
        Seconds
    wcs0 : galsim.GSFitsWCS
    tracker : xfiles.Tracker
    psf : galsim.GSObject
    image : galsim.Image
    observer : ssapy.EarthObserver
    nphot : float
    nsplit : int
    propagator : ssa.Propagator instance

    Returns
    -------
    sats : astropy.table.Table
        Columns for
            ra0, dec0, ra1, dec1 : float
                degrees
            x0, y0, x1, y1 : float
                Image coordinates
    """
    table = Table()
    ra0 = []
    dec0 = []
    ra1 = []
    dec1 = []
    xs0 = []
    ys0 = []
    xs1 = []
    ys1 = []
    # only difference compared to draw_stars is that sat position needs to be
    # recomputed too
    dt = exptime/nsplit
    ts = t0 + np.linspace(0, exptime, nsplit+1)*u.s

    for idx, orb in enumerate(orbit):
        ra, dec, _ = ssapy.radec(orb, ts, observer=observer,
                                propagator=propagator)

        xy0 = wcs0.posToImage(
            galsim.CelestialCoord(
                ra[0]*galsim.radians,
                dec[0]*galsim.radians
            )
        )
        ra0.append(ra[0])
        dec0.append(dec[0])
        ra1.append(ra[-1])
        dec1.append(dec[-1])
        xs0.append(xy0.x)
        ys0.append(xy0.y)

        for i in range(1, nsplit+1):
            t = t0 + i*dt*u.s
            boresight = tracker.get_boresight(t)
            rot_sky_pos = tracker.get_rot_sky_pos(t)
            wcst = transform_wcs(wcs0, boresight, rot_sky_pos)
            xy1 = wcst.posToImage(
                galsim.CelestialCoord(
                    ra[i]*galsim.radians,
                    dec[i]*galsim.radians
                )
            )
            obj = galsim.Convolve(
                getLineGSObject(
                    wcs0.posToWorld(xy0),
                    wcs0.posToWorld(xy1)
                ),
                psf
            )
            xy = (xy0+xy1)/2
            local_wcs = wcs0.local(xy)
            stamp = (obj*nphot[idx]/nsplit).drawImage(
                wcs=local_wcs,
                center=xy,
                method='phot'
            )
            bounds = stamp.bounds & image.bounds
            if bounds.area() > 0:
                image[bounds] += stamp[bounds]
            xy0 = xy1
        xs1.append(xy1.x)
        ys1.append(xy1.y)

    table['x0'] = np.array([x for x in xs0])
    table['y0'] = np.array([y for y in ys0])
    table['ra0'] = np.array([np.rad2deg(ra) for ra in ra0])
    table['dec0'] = np.array([np.rad2deg(dec) for dec in dec0])
    table['ra1'] = np.array([np.rad2deg(ra) for ra in ra1])
    table['dec1'] = np.array([np.rad2deg(dec) for dec in dec1])
    table['x1'] = np.array([x for x in xs1])
    table['y1'] = np.array([y for y in ys1])
    table['nphot'] = np.array([n for n in nphot])

    # Add in FITS xy
    table['x0_FITS'] = np.array([x0 - image.xmin + 1 for x0 in table['x0']])
    table['x1_FITS'] = np.array([x1 - image.xmin + 1 for x1 in table['x1']])
    table['y0_FITS'] = np.array([y0 - image.ymin + 1 for y0 in table['y0']])
    table['y1_FITS'] = np.array([y1 - image.ymin + 1 for y1 in table['y1']])

    return table, wcst


def random_dark_time(*, t_day, site, rng):
    """Find a random time during the next night.

    Parameters
    ----------
    t_day : astropy.time.Time
        A point during the day before target sunset.
    site : astroplan.site
    rng : np.random.Generator

    Returns
    -------
    t_dark : astropy.time.Time
        Random time during the night following t_previous_day.
    """
    sunset = site.sun_set_time(t_day, which='next')
    sunrise = site.sun_rise_time(sunset, which='next')
    return Time(rng.uniform(sunset.gps, sunrise.gps), format='gps')


def random_boresight(*, observer, t0, horizon, rng):
    """Uniformly pick boresight direction from spherical cap above desired
    horizon.

    Parameters
    ----------
    observer : ssapy.EarthObserver
    t0 : astropy.time.Time
    horizon : galsim.Angle
    rng : np.random.Generator

    Returns
    -------
    boresight : galsim.CelestialCoord
    """
    # Pick a boresight uniformly in spherical cap 20 degrees above
    # horizon
    zmin = np.sin(horizon)
    z = rng.uniform(zmin, 1)
    alt = np.arctan(z/np.sqrt(1-z**2))
    az = rng.uniform(0, 2*np.pi)

    # Transform to ra/dec
    sc = AltAz(
        az*u.rad, alt*u.rad, obstime=t0, location=observer._location
    )
    radec = sc.transform_to(GCRS())
    boresight = galsim.CelestialCoord(
        radec.ra.rad*galsim.radians,
        radec.dec.rad*galsim.radians
    )
    return boresight


def random_disk(center, radius, rng):
    """Return point roughly distributed uniformly in spherical cap.

    Notes
    -----
    Approximation is better for smaller radii

    Parameters
    ----------
    center : galsim.CelestialCoord
    radius : galsim.Angle
    rng : np.random.Generator

    Returns
    -------
    point : galsim.CelestialCoord
    """
    arclength = np.sqrt(rng.uniform(0, radius.rad**2))*galsim.radians
    return center.greatCirclePoint(
        galsim.CelestialCoord.from_xyz(*rng.normal(size=3)),
        arclength
    )


def generate_orbit(height, coord, heading, vperp, vpar, observer, t0):
    """Get an orbit from observational parameters and some velocity assumptions.

    Parameters
    ----------
    height : float
        height of sat above Earth's surface in meters
    coord : galsim.CelestialCoord
        Line of sight from observer to satellite
    heading : galsim.Angle
        Positional (measured from North through East) of satellite's velocity
    vperp : float
        Sat velocity magnitude perpendicular to r_sat in m/s
    vpar : float
        Sat velocity parallel to r_sat in m/s
    observer : ssapy.EarthObserver
    t0 : astropy.time.Time

    Returns
    -------
    out : ssapy.Orbit
    """

    r_obs, _ = observer.getRV(t0)
    los = np.array(coord.get_xyz())
    b = 2*np.dot(los, r_obs)
    c = (
        np.sum(r_obs**2) -
        (ssapy.constants.WGS72_EARTH_RADIUS+height)**2
    )
    disc = b**2 - 4*c
    rho = (-b + np.sqrt(disc))/2
    r_sat = r_obs + rho*los

    # Get North in plane perp to los
    alpha = coord.ra
    delta = coord.dec
    north_los = np.array([
        -np.cos(delta) * np.cos(alpha),
        -np.cos(delta) * np.sin(alpha),
        np.sin(delta),
    ])
    axis = np.array(coord.get_xyz())
    # Rotate to heading.
    n_head = north_los*np.cos(heading)
    n_head += np.cross(axis, north_los) * np.sin(heading)
    n_head += axis*np.dot(axis, north_los)*(1 - np.cos(heading))

    # project into plane perp to r_sat
    n_sat = ssapy.utils.normed(r_sat)
    vpar3 = n_sat*vpar
    vperp3 = n_head - n_sat*np.dot(n_head, n_sat)
    vperp3 *= vperp/ssapy.utils.norm(vperp3)
    v_sat = vperp3 + vpar3
    return ssapy.Orbit(r_sat, v_sat, t0)


def getLineGSObject(p0, p1):
    """Make a line GSObject connecting two celestial coordinates.

    Parameters
    ----------
    p0, p1 : galsim.CelestialCoord
        Coordinates to connect.

    Returns
    -------
    obj : GSObject
    """
    ncp = galsim.CelestialCoord(0*galsim.degrees, 90*galsim.degrees)
    q = p0.angleBetween(ncp, p1)
    length = p0.distanceTo(p1)/galsim.arcsec
    return galsim.Box(
        1e-12,
        max(length, 1e-12),
    ).rotate(q)


# Put a tqdm progressbar on last line, but still be able to have a scrolling print above that.
# https://stackoverflow.com/questions/36986929/redirect-print-command-in-python-script-through-tqdm-write/37243211#37243211
class DummyFile(object):
    file = None
    def __init__(self, file):
        self.file = file

    def write(self, x):
        from tqdm import tqdm

        # only print every other output for some reason???
        tqdm.write(x, end="", file=self.file)

    def __eq__(self, other):
        return other is self.file


@contextlib.contextmanager
def nostdout():
    save_stdout = sys.stdout
    sys.stdout = DummyFile(sys.stdout)
    yield
    sys.stdout = save_stdout
