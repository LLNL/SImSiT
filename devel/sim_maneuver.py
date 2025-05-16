#!/usr/bin/env python
# encoding: utf-8
# Design
# - Pick sats from Celestrak DB
# - Pick randomish values for area/mass/drag-coef/solar-rad-coef
# - Filter to LEO: period < 2 hours
# - Observation sites are GEODSS + national labs
# - Propagate with ~realistic propagator for 24 hours
#   - 20,20 gravity
#   - sun/moon point sources
#   - drag from Harris-Priester atmospheric density model
#   - (constant) solar radiation pressure, accel is directly away from Sun
# - For ~half of sats, pick a random time between [10-20] hours to make an
#   impulse maneuver.
# - For each site, pick between 5 and 10 instances where alt is above 20 deg
#   - If fewer than 5 are available, use all.
# - For each instance, record paired observations (ra/dec/time tuples) with
#   \Delta t = 3 seconds
#
# - Output tables
# - Sat table:
#   - satID
#   - TLE
#   - area/mass/CD/CR
#   - time of impulse
#   - \Delta v of impulse
# - Observation table (per satellite):
#   - ra/dec
#   - time
#   - site, site_r, site_v tuples
# - Ephemeris table (per satellite):
#   - r
#   - v
#   - t
# - Additional questions/comments:
# - How should we simulate endpoint uncertainty?
# - No aberration
#   - BTW, think that only diurnal aberration affects sats, annual+diurnal
#     affects stars
# - No speed-of-light time delay
# - No atm refraction
# - ra/dec topocentric parallel to ICRS, so no precession/nutation/TEME.

import os
import time
import glob
from textwrap import dedent

import numpy as np
from astropy.time import Time
from astropy.table import Table
from tqdm import tqdm

import ssa


def ten(sixty):
    import re
    ds = 1
    d, m, s = [float(i) for i in re.split(':', sixty)]
    if str(d)[0] == '-':
        ds = -1
        d = abs(d)
    return ds*(d+m/60+s/3600)


# Use to read 3 lines at a time from TLE files.
def group(iterator, count):
    itr = iter(iterator)
    while True:
        try:
            yield tuple([next(itr) for i in range(count)])
        except StopIteration:
            return


def readTLEs(fn):
    from sgp4.io import twoline2rv
    from sgp4.propagation import sgp4
    from sgp4.earth_gravity import wgs84
    tles = []
    satids = []
    names = []
    rs = []
    vs = []
    ts = []
    with open(fn) as f:
        lines = [l.rstrip() for l in f.readlines()]
    for line0, line1, line2 in group(lines, 3):
        sat = twoline2rv(line1, line2, wgs84)
        r, v = sgp4(sat, 0)
        t = Time(sat.jdsatepoch, format='jd')
        tles.append((line1, line2))
        satids.append(int(line1[2:7]))
        names.append(line0)
        rs.append(r)
        vs.append(v)
        ts.append(t)

    return tles, satids, names, rs, vs, ts


def main(args):
    # National labs + GEODSS
    sitetxt = dedent(
    """ LLNL    37:41:13.07 -121:42:21.15  185.0
        SNLCA   37:40:44.25 -121:42:20.29  195.0
        SNLNM   35:03:14.10 -106:31:49.44 1673.0
        LBNL    37:52:33.25 -122:15:00.20  261.0
        SLAC    37:25:11.47 -122:12:09.74   83.0
        ORNL    35:55:51.96  -84:18:35.92  261.0
        ANL     41:43:05.81  -87:58:43.93  220.0
        Ames    42:01:47.89  -93:38:53.95  304.0
        BNL     40:51:51.35  -72:52:30.61   20.0
        PPPL    40:21:00.14  -74:36:10.92   31.0
        PNNL    46:20:42.70 -119:16:45.23  120.0
        FNAL    41:50:26.43  -88:16:45.82  233.0
        TJNAF   37:05:50.37  -76:29:13.86   19.0
        LANL    35:50:38.61 -106:17:13.78 2163.0
        NREL    39:44:26.61 -105:10:07.01 1760.0
        SRNL    33:20:37.72  -81:44:05.92  112.0
        NETL    39:40:13.54  -79:58:30.44  296.0
        INL     43:32:28.20 -112:30:14.50 1593.0
        Maui    20:42:29.89 -156:15:27.17 3028.0
        Socorro 33:49:01.92 -106:39:35.64 1525.0
        DG      -7:24:42.21   72:27:10.64    4.0"""
    )
    sites = []
    for line in sitetxt.split('\n'):
        name, lat, lon, elev = line.split()
        site = ssa.EarthObserver(ten(lon), ten(lat), float(elev))
        site._location
        site.name = name
        sites.append(site)

    if args.nsite is not None:
        sites = sites[:args.nsite]

    satTable = Table(
        names=[
            'satID', 'j', 'r', 'v', 't', 'TLE',
            'area', 'mass', 'CD', 'CR',
            't_impulse', 'dv_impulse'
        ],
        dtype=["i4","i4","3f8","3f8","f8","2U69","f8","f8","f8","f8","f8","3f8"]
    )

    # Read in TLEs
    tles = []
    satIDs = []
    names = []
    rs = []
    vs = []
    ts = []
    print("Reading TLEs")
    for fn in glob.glob("data/*.txt"):
        tle, sid, n, r, v, t = readTLEs(fn)
        tles.extend(tle)
        satIDs.extend(sid)
        names.extend(n)
        rs.extend(r)
        vs.extend(v)
        ts.extend(t)
    rs = np.array(rs)*1e3
    vs = np.array(vs)*1e3
    ts = np.array(ts)
    tles = np.array(tles)
    satIDs = np.array(satIDs)
    names = np.array(names)
    orbits = ssa.Orbit(rs, vs, ts)

    # Filter to LEO
    w = orbits.period < 7200
    tles = tles[w]
    satIDs = satIDs[w]
    names = names[w]
    orbits = orbits[w]

    np.random.seed(57721)

    # Use a ~realistic propagator.
    accel = (
        ssa.AccelHarmonic(20, 20) + ssa.AccelMoon() + ssa.AccelSun()
        + ssa.AccelDrag() + ssa.AccelSolRad()
    )

    for j in range(args.nsat):
        # Pick a random orbit to simulate
        i = np.random.choice(len(orbits))
        orbit = orbits[i]
        name = names[i]
        satID = satIDs[i]
        tle = tles[i]

        # Pick some randomish values for mass, area, drag-coefficient,
        # reflection coefficient.
        mass = np.exp(np.random.uniform(np.log(1.5e0), np.log(1.5e4)))
        area = (mass)**(2/3) * np.random.uniform(0.01, 0.1)
        CD = np.random.uniform(2.0, 2.3)
        CR = np.random.uniform(1.2, 1.7)

        # Randomly decide if a maneuver happens or not, as well as
        # magnitude.  For now, assume direction is directly prograde.
        do_maneuver = np.random.choice([True, False], p=[0.5, 0.5])
        if do_maneuver:
            dt_maneuver = 10*3600 + np.random.uniform(10*3600.)
            # between 1 cm/s and 1 m/s
            maneuver_mag = np.random.uniform(0.01, 1.0)
        else:
            dt_maneuver = 86401.0
            maneuver_mag = 0.0

        propagator = ssa.RK78Propagator(
            accel, h=10.0, mass=mass, area=area, CD=CD, CR=CR
        )

        t0 = orbit.t
        eph_t_catalog = t0 + np.arange(0, 86400, 10)
        eph_r_catalog = []
        eph_v_catalog = []
        eph_alt_catalog = {site.name:[] for site in sites}
        eph_az_catalog = {site.name:[] for site in sites}

        obs_t_catalog = []
        obs_site_catalog = []
        obs_ra_catalog = []
        obs_dec_catalog = []
        obs_range_catalog = []
        obs_alt_catalog = []
        obs_az_catalog = []

        # Assemble ephemeris, being careful to apply impulse at appropriate
        # t_maneuver
        print("Computing ephemeris")
        t_prop = eph_t_catalog[eph_t_catalog < t0+dt_maneuver]
        t_prop = np.concatenate([t_prop, [t0+dt_maneuver]])

        r, v = ssa.rv(
            orbit,
            t_prop,
            propagator=propagator
        )
        eph_r_catalog.extend(r[:-1])  # don't include t_maneuver
        eph_v_catalog.extend(v[:-1])
        for site in tqdm(sites):
            alt, az = ssa.altaz(
                orbit,
                t_prop,
                observer=site,
                propagator=propagator
            )
            eph_alt_catalog[site.name].extend(alt[:-1])
            eph_az_catalog[site.name].extend(az[:-1])

        if do_maneuver:
            dv = maneuver_mag * ssa.utils.normed(eph_v_catalog[-1])
            perturbed_orbit = ssa.Orbit(
                eph_r_catalog[-1],
                eph_v_catalog[-1] + dv,
                t0+dt_maneuver
            )

            t_prop = eph_t_catalog[eph_t_catalog > t0+dt_maneuver]
            t_prop = np.concatenate([[t0+dt_maneuver], t_prop])
            r, v = ssa.rv(
                perturbed_orbit,
                t_prop,
                propagator=propagator
            )
            eph_r_catalog.extend(r[1:])  # don't include t_maneuver
            eph_v_catalog.extend(v[1:])
            print("Computing perturbed emphemeris")
            for site in tqdm(sites):
                alt, az = ssa.altaz(
                    perturbed_orbit,
                    t_prop,
                    observer=site,
                    propagator=propagator
                )
                eph_alt_catalog[site.name].extend(alt[1:])
                eph_az_catalog[site.name].extend(az[1:])
        else:
            dv = (0,0,0)

        # Have our ephemeris, now can go back and determine when visible and
        # make some "observations".

        print("Observing")
        for site in tqdm(sites):
            n_obs = np.random.randint(5, 10)

            w = eph_alt_catalog[site.name] > np.deg2rad(20.0)
            w = np.random.permutation(w.nonzero()[0])[:n_obs]
            t_obs = eph_t_catalog[w]
            t_obs += np.random.uniform(-5, 5, len(t_obs))
            t_obs = np.hstack([t_obs, t_obs+3])
            t_obs = np.sort(t_obs)

            wpre = t_obs < t0+dt_maneuver
            alt = np.zeros(len(t_obs))
            az = np.zeros_like(alt)
            ra = np.zeros_like(alt)
            dec = np.zeros_like(alt)
            range_ = np.zeros_like(alt)

            if np.any(wpre):
                alt[wpre], az[wpre] = ssa.altaz(
                    orbit,
                    t_obs[wpre],
                    observer=site,
                    propagator=propagator
                )
                ra[wpre], dec[wpre], range_[wpre] = ssa.radec(
                    orbit,
                    t_obs[wpre],
                    observer=site,
                    propagator=propagator
                )
            if np.any(~wpre):
                alt[~wpre], az[~wpre] = ssa.altaz(
                    perturbed_orbit,
                    t_obs[~wpre],
                    observer=site,
                    propagator=propagator
                )
                ra[~wpre], dec[~wpre], range_[~wpre] = ssa.radec(
                    perturbed_orbit,
                    t_obs[~wpre],
                    observer=site,
                    propagator=propagator
                )

            obs_t_catalog.extend(t_obs)
            obs_site_catalog.extend([site.name]*len(ra))
            obs_ra_catalog.extend(ra)
            obs_dec_catalog.extend(dec)
            obs_range_catalog.extend(range_)
            obs_alt_catalog.extend(alt)
            obs_az_catalog.extend(az)

        # Can add this now...
        satTable.add_row((
            satID, j, orbit.r, orbit.v, orbit.t, tle,
            area, mass, CD, CR, t0+dt_maneuver, dv
        ))

        # Write out observation and ephemeris tables here.
        obsTable = Table()
        obsTable['obs_site'] = obs_site_catalog
        obsTable['obs_t'] = obs_t_catalog
        obsTable['obs_ra'] = obs_ra_catalog
        obsTable['obs_dec'] = obs_dec_catalog
        obsTable['obs_range'] = obs_range_catalog
        obsTable['obs_az'] = obs_az_catalog
        obsTable['obs_alt'] = obs_alt_catalog

        obsFile = os.path.join(args.outdir, f"{j:03d}_observation.fits")
        obsTable.write(obsFile, format='fits')

        ephTable = Table()
        ephTable['eph_t'] = eph_t_catalog
        ephTable['eph_r'] = eph_r_catalog
        ephTable['eph_v'] = eph_v_catalog

        ephFile = os.path.join(args.outdir, f"{j:03d}_emphemeris.fits")
        ephTable.write(ephFile, format='fits')

    # Write out sat table here.
    satFile = os.path.join(args.outdir, f"sat.fits")
    satTable.write(satFile, format='fits')


if __name__ == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("--nsat", default=20, type=int)
    parser.add_argument("--nsite", default=None, type=int)
    parser.add_argument("--outdir", default='output', type=str)
    args = parser.parse_args()

    main(args)
