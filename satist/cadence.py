"""
Observation cadence strategies for satellite tracking.
"""

import numpy as np
from astropy.time import Time, TimeDelta


class SimpleCadence:
    """ Uniformly take exposures during satellite pass.

    Parameters
    ----------
    exptime : astropy.time.TimeDelta
        Time between shutter open and shutter close.
    delay : astropy.time.TimeDelta
        Time between successive shutter openings in seconds.
    """
    def __init__(self, exptime, delay):
        self.exptime = TimeDelta(exptime)
        self.delay = TimeDelta(delay)

    def get_times(self, pass_):
        """ Get exposure start and end times.

        Parameters
        ----------
        pass_ : dict
            From ssapy.compute.refine_pass

        Returns
        -------
        tStart, tEnd : astropy.time.Time
            Exposure start/end times
        """
        if not pass_['illumAtStart'] and not pass_['illumAtEnd']:
            return Time([], format='mjd'), Time([], format='mjd')

        if pass_['illumAtStart']:
            tStart = pass_['tStart']
        else:
            tStart = pass_['tTerminator']

        if pass_['illumAtEnd']:
            tEnd = pass_['tEnd']
        else:
            tEnd = pass_['tTerminator']

        ts = Time(
            np.arange(
                tStart.gps,
                tEnd.gps,
                self.delay.sec
            ),
            format='gps'
        )
        return ts, ts + self.exptime
