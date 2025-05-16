from abc import ABC, abstractmethod

import galsim
import numpy as np
import ssapy


class Tracker(ABC):
    @abstractmethod
    def get_boresight(self, time):
        pass

    @abstractmethod
    def get_rot_sky_pos(self, time):
        pass


class OrbitTracker(Tracker):
    """ Track a ssapy.Orbit

    Parameters
    ----------
    orbit : ssapy.Orbit
    observer : ssapy.EarthObserver
    t0 : astropy.time.Time
    rot_sky_pos0 : galsim.Angle
        Initial rotation between camera and sky
    mount : {'EQ'}
        Mount type for rotator tracking.
    propagator : ssa.Propagator
        Propagator to use for satellite motion
    """
    def __init__(self, *, orbit, observer, t0, rot_sky_pos0, mount='EQ',
                 propagator=None):
        if mount.upper() not in ['EQ']:
            raise ValueError(f"Unknown mount type: {mount}")
        self.orbit = orbit
        self.observer = observer
        self.t0 = t0
        self.rot_sky_pos0 = rot_sky_pos0
        self.mount = mount
        self.propagator = propagator

    def get_boresight(self, time):
        ra, dec, _ = ssapy.radec(self.orbit, time, observer=self.observer,
                                 propagator=self.propagator)
        return galsim.CelestialCoord(ra*galsim.radians, dec*galsim.radians)

    def get_rot_sky_pos(self, time):
        # EQ mounted telescope naturally preserves rotSkyPos
        # TODO: other mount types
        return self.rot_sky_pos0


class InertialTracker(Tracker):
    """ Track by slewing at a constant rate w.r.t. the sky.

    Note: sidereal tracking is obtained by setting rot_rate = 0

    Parameters
    ----------
    t0 : astropy.time.Time
    boresight0 : galsim.CelestialCoord
        Initial boresight at time t0.
    rot_sky_pos0 : galsim.Angle
        Initial rotation between camera and sky
    rot_axis : galsim.CelestialCoord
    rot_rate : galsim.Angle
        Rotation angle over 1 second.
    mount : {'EQ'}
        Mount type for rotator tracking.
    """
    def __init__(
        self, *, t0, boresight0, rot_sky_pos0, rot_axis, rot_rate, mount='EQ'
    ):
        if mount.upper() not in ['EQ']:
            raise ValueError(f"Unknown mount type: {mount}")
        self.t0 = t0
        self.boresight0 = boresight0
        self.rot_sky_pos0 = rot_sky_pos0
        self.rot_axis = rot_axis
        self.rot_rate = rot_rate
        self.mount = mount

        if self.rot_rate != 0:
            self._xyz0 = np.array(self.boresight0.get_xyz())
            axis_xyz = np.array(self.rot_axis.get_xyz())
            self._cross = np.cross(axis_xyz, self._xyz0)
            self._dot = np.dot(axis_xyz, self._xyz0)*axis_xyz

    def get_boresight(self, time):
        if self.rot_rate == 0.0:
            return self.boresight0
        theta = ((time-self.t0).sec)*self.rot_rate
        # Rodriguez rotation formula:
        # https://en.wikipedia.org/wiki/Axis%E2%80%93angle_representation
        xyz = (
            np.cos(theta)*self._xyz0 +
            np.sin(theta)*self._cross +
            (1-np.cos(theta))*self._dot
        )
        return galsim.CelestialCoord.from_xyz(*xyz)

    def get_rot_sky_pos(self, time):
        # EQ mounted telescope naturally preserves rotSkyPos
        # TODO: other mount types
        return self.rot_sky_pos0


def SiderealTracker(boresight0, rot_sky_pos0):
    return InertialTracker(
        t0=None,
        boresight0=boresight0,
        rot_sky_pos0=rot_sky_pos0,
        rot_axis=None,
        rot_rate=0.0,
        mount='EQ'
    )

def transform_wcs(
    wcs0, boresight1, rot_sky_pos1
):
    """ Transform WCS by recentering and reorienting it.

    Note: assumes that initial boresight is at wcs0.center and that initial
    rot_sky_pos is obtainable from the wcs0.cd matrix.

    Parameters
    ----------
    wcs0 : galsim.GSFitsWCS
        Initial WCS
    boresight1 : galsim.CelestialCoord
        New boresight
    rot_sky_pos1 : galsim.Angle
        New rot_sky_pos

    Returns
    -------
    out : galsim.GSFitsWCS
    """
    rot_sky_pos0 = np.arctan2(-wcs0.cd[0,1], -wcs0.cd[0,0])

    dth = rot_sky_pos1.rad - rot_sky_pos0

    sth, cth = np.sin(dth), np.cos(dth)

    wcs = wcs0.copy()
    R = np.array([[cth, -sth], [sth, cth]])
    wcs.cd = R @ wcs.cd
    if hasattr(wcs, 'cdinv'):
        del wcs.cdinv
    wcs.center = boresight1
    return wcs
