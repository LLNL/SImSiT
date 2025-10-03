"""
World Coordinate System (WCS) generation with radial optical distortion.
"""

import galsim
import numpy as np


def radialWCS(
    th, dthdr,
    world_origin,
    rot_sky_pos=0*galsim.degrees,
    n=10, order=3, verbose=False
):
    """Make a WCS from a radial distortion polynomial

    Parameters
    ----------
    th : array
        Field angles in degrees
    dthdr: array
        Radial plate scale in arcsec per pixel
    world_origin:  galsim.CelestialCoord
        Origin of ra, dec
    rot_sky_pos: galsim.Angle
        Rotation angle in radians
    n: int
        Number of control points to use
    order: int
        Order of SIP part of fitted WCS

    Returns
    -------
        wcs : galsim.GSFitsWCS
    """
    import warnings

    from scipy.integrate import IntegrationWarning, quad
    from scipy.interpolate import interp1d

    u = np.deg2rad(np.linspace(-th[-1], th[-1], n))
    u, v = np.meshgrid(u, u)
    rho = np.hypot(u, v)
    w = rho <= np.deg2rad(th[-1])
    u = u[w]
    v = v[w]
    rho = rho[w]

    interp = interp1d(th, dthdr, kind='cubic')  # deg -> arcsec/pix
    integrand = lambda arcsec: 1./interp(arcsec/3600)  # arcsec -> pix/arcsec

    x = np.empty_like(u)
    y = np.empty_like(u)

    for idx in np.ndindex(u.shape):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=IntegrationWarning)
            r = quad(integrand, 0, np.rad2deg(rho[idx])*3600)[0]
        x[idx] = r*u[idx]/rho[idx]
        y[idx] = r*v[idx]/rho[idx]

    sth, cth = np.sin(rot_sky_pos), np.cos(rot_sky_pos)
    R = np.array([[cth, -sth], [sth, cth]])
    x, y = R @ np.array([x, y])

    ra, dec = world_origin._deproject(u, v, projection='postel')
    wcs = galsim.FittedSIPWCS(x, y, ra, dec, order=order)

    if verbose:
        # Report residual
        x1, y1 = wcs.radecToxy(ra, dec, units="radians")
        print(
            "WCS inversion x residuals:"
            f"  mean = {np.mean(x-x1)}"
            f"  std = {np.std(x-x1)}"
        )
        print(
            "WCS inversion y residuals:"
            f"  mean = {np.mean(y-y1)}"
            f"  std = {np.std(y-y1)}"
        )

    return wcs
