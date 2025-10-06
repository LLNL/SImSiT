"""
Telescope tracking system implementations for observation simulations.
"""

from abc import ABC, abstractmethod

import galsim
import numpy as np
import ssapy


class Tracker(ABC):
    """Abstract base class for telescope tracking systems.
    
    Defines the interface for tracking systems that determine telescope
    pointing direction and camera orientation as a function of time.
    """
    
    @abstractmethod
    def get_boresight(self, time):
        """Get the telescope boresight direction at a given time.
        
        Parameters
        ----------
        time : astropy.time.Time
            Time at which to compute boresight
            
        Returns
        -------
        galsim.CelestialCoord
            Boresight pointing direction
        """
        pass

    @abstractmethod
    def get_rot_sky_pos(self, time):
        """Get the camera rotation angle relative to sky at a given time.
        
        Parameters
        ----------
        time : astropy.time.Time
            Time at which to compute rotation
            
        Returns
        -------
        galsim.Angle
            Position angle of camera "up" direction with respect to North
        """
        pass


class OrbitTracker(Tracker):
    """Track a satellite orbit with the telescope.
    
    Implements tracking of a moving satellite by continuously updating the
    telescope boresight to follow the satellite's predicted position based
    on its orbit.

    Parameters
    ----------
    orbit : ssapy.Orbit
        Satellite orbit to track
    observer : ssapy.EarthObserver
        Ground-based observer location
    t0 : astropy.time.Time
        Initial time for tracking
    rot_sky_pos0 : galsim.Angle
        Initial rotation angle between camera and sky coordinate systems
    mount : {'EQ'}, optional
        Mount type for rotator tracking (default: 'EQ' for equatorial mount)
    propagator : ssapy.Propagator, optional
        Propagator to use for satellite motion prediction
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
        """Get boresight direction by computing satellite position at given time.
        
        Parameters
        ----------
        time : astropy.time.Time
            Time at which to compute satellite position
            
        Returns
        -------
        galsim.CelestialCoord
            Right ascension and declination of satellite
        """
        ra, dec, _ = ssapy.radec(self.orbit, time, observer=self.observer,
                                 propagator=self.propagator)
        return galsim.CelestialCoord(ra*galsim.radians, dec*galsim.radians)

    def get_rot_sky_pos(self, time):
        """Get camera rotation angle.
        
        For equatorial mounts, the rotation angle is naturally preserved
        during tracking.
        
        Parameters
        ----------
        time : astropy.time.Time
            Time at which to compute rotation (unused for EQ mounts)
            
        Returns
        -------
        galsim.Angle
            Camera rotation angle (constant for EQ mounts)
        """
        # EQ mounted telescope naturally preserves rotSkyPos
        # TODO: other mount types
        return self.rot_sky_pos0


class InertialTracker(Tracker):
    """Track by slewing at a constant rate with respect to the inertial sky.
    
    Implements tracking where the telescope rotates around a fixed axis at a
    constant angular rate. Sidereal tracking is a special case with zero
    rotation rate.

    Parameters
    ----------
    t0 : astropy.time.Time
        Initial time for tracking
    boresight0 : galsim.CelestialCoord
        Initial boresight direction at time t0
    rot_sky_pos0 : galsim.Angle
        Initial rotation angle between camera and sky coordinate systems
    rot_axis : galsim.CelestialCoord
        Axis of rotation in sky coordinates
    rot_rate : galsim.Angle
        Rotation angle per second
    mount : {'EQ'}, optional
        Mount type for rotator tracking (default: 'EQ' for equatorial mount)
        
    Notes
    -----
    Uses the Rodriguez rotation formula to compute boresight at arbitrary times.
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
        """Get boresight direction at given time using Rodriguez rotation formula.
        
        Parameters
        ----------
        time : astropy.time.Time
            Time at which to compute boresight
            
        Returns
        -------
        galsim.CelestialCoord
            Boresight direction rotated by (time - t0) * rot_rate around rot_axis
        """
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
        """Get camera rotation angle.
        
        For equatorial mounts, the rotation angle is naturally preserved
        during tracking.
        
        Parameters
        ----------
        time : astropy.time.Time
            Time at which to compute rotation (unused for EQ mounts)
            
        Returns
        -------
        galsim.Angle
            Camera rotation angle (constant for EQ mounts)
        """
        # EQ mounted telescope naturally preserves rotSkyPos
        # TODO: other mount types
        return self.rot_sky_pos0


def SiderealTracker(boresight0, rot_sky_pos0):
    """Create a sidereal tracker (stationary with respect to stars).
    
    Convenience function that creates an InertialTracker with zero rotation
    rate, resulting in standard sidereal tracking where the telescope remains
    fixed relative to the celestial sphere.
    
    Parameters
    ----------
    boresight0 : galsim.CelestialCoord
        Fixed boresight direction
    rot_sky_pos0 : galsim.Angle
        Fixed rotation angle between camera and sky coordinate systems
        
    Returns
    -------
    InertialTracker
        Tracker configured for sidereal tracking (zero rotation rate)
    """
    return InertialTracker(
        t0=None,
        boresight0=boresight0,
        rot_sky_pos0=rot_sky_pos0,
        rot_axis=None,
        rot_rate=0.0,
        mount='EQ'
    )

def transform_wcs(wcs0, boresight1, rot_sky_pos1):
    """Transform WCS by recentering and reorienting it.
    
    Creates a new WCS by updating the center position and rotating the
    coordinate system. Useful for propagating WCS through time as the
    telescope tracks.

    Parameters
    ----------
    wcs0 : galsim.GSFitsWCS
        Initial WCS to transform
    boresight1 : galsim.CelestialCoord
        New boresight (field center)
    rot_sky_pos1 : galsim.Angle
        New camera rotation angle with respect to North

    Returns
    -------
    galsim.GSFitsWCS
        Transformed WCS with updated center and orientation
        
    Notes
    -----
    Assumes that:
        - Initial boresight is at wcs0.center
        - Initial rot_sky_pos is encoded in the wcs0.cd matrix
    The transformation applies a rotation matrix to the CD matrix and
    updates the reference position.
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