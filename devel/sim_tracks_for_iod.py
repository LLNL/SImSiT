#!/usr/bin/env python
# encoding: utf-8
#
# Simulated tracks for IOD validation tests
# [modified from the sim_maneuver.py script]
#
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

import numpy as np
from astropy.time import Time
from astropy.table import Table
from tqdm import tqdm

import ssa
import xfiles.utils as xutils


def read_TLEs():
    tles = []
    satIDs = []
    names = []
    rs = []
    vs = []
    ts = []
    print("Reading TLEs")
    for fn in glob.glob("tles/*.txt"):
        tle, sid, n, r, v, t = xutils.readTLEs(fn)
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

    return tles, satIDs, names, orbits


def main(args):
    sites = xutils.get_sites()
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
    tles, satIDs, names, orbits = read_TLEs()

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

    track_id_ndx = 1

    for j in range(args.nsat):
        # Pick a random orbit to simulate
        i = np.random.choice(len(orbits))
        orbit = orbits[i]
        name = names[i]
        satID = satIDs[i]
        tle = tles[i]

        print("===== Simulating satellite {}, {:d} / {:d} =====".format(names[i], 
            j+1, args.nsat))

        # Pick some randomish values for mass, area, drag-coefficient,
        # reflection coefficient.
        mass = np.exp(np.random.uniform(np.log(1.5e0), np.log(1.5e4)))
        area = (mass)**(2/3) * np.random.uniform(0.01, 0.1)
        CD = np.random.uniform(2.0, 2.3)
        CR = np.random.uniform(1.2, 1.7)

        # Use Kepler propagator for faster development and testing
        propagator = ssa.KeplerianPropagator()
        #        
        # propagator = ssa.RK78Propagator(
        #     accel, h=10.0, mass=mass, area=area, CD=CD, CR=CR
        # )

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
        obs_track_catalog = []

        # Assemble ephemeris
        print("Computing ephemeris")
        t_prop = eph_t_catalog

        r, v = ssa.rv(
            orbit,
            t_prop,
            propagator=propagator
        )
        eph_r_catalog.extend(r)
        eph_v_catalog.extend(v)
        for site in tqdm(sites):
            alt, az = ssa.altaz(
                orbit,
                t_prop,
                observer=site,
                propagator=propagator
            )
            eph_alt_catalog[site.name].extend(alt)
            eph_az_catalog[site.name].extend(az)
            dv = (0,0,0)

        # Have our ephemeris, now can go back and determine when visible and
        # make some "observations".

        print("Observing")
        for site in tqdm(sites):
            # Pick between 5 and 10 instances where observations are possible.
            # Use all available instances if fewer than 5 are available.
            n_obs = np.random.randint(5, 10)

            # Require alitutude greater than 20 degrees for good observing
            w = eph_alt_catalog[site.name] > np.deg2rad(20.0)
            w = np.random.permutation(w.nonzero()[0])[:n_obs]
            t_obs = eph_t_catalog[w]
            n = len(t_obs)
            t_obs += np.random.uniform(-5, 5, len(t_obs))
            t_obs = np.hstack([t_obs + i * args.exposure_time_sec
                               for i in range(args.num_obs_per_track)])
            t_obs = np.sort(t_obs)

            track_id = np.hstack([np.ones(args.num_obs_per_track, dtype=int)* (track_id_ndx + i)
                                  for i in range(n)])
            track_id_ndx += n

            alt = np.zeros(len(t_obs))
            az = np.zeros_like(alt)
            ra = np.zeros_like(alt)
            dec = np.zeros_like(alt)
            range_ = np.zeros_like(alt)

            # FIXME: This error happens here
            #    ValueError: zero-size array to reduction operation maximum which has no identity
            #    Possibly because t_obs has length zero?
            alt, az = ssa.altaz(
                orbit,
                t_obs,
                observer=site,
                propagator=propagator
            )
            ra, dec, range_ = ssa.radec(
                orbit,
                t_obs,
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
            obs_track_catalog.extend(track_id)

        # Can add this now...
        satTable.add_row((
            satID, j, orbit.r, orbit.v, orbit.t, tle,
            area, mass, CD, CR, t0, dv
        ))

        # Write out observation and ephemeris tables here.
        obsTable = Table()
        obsTable['obs_site'] = obs_site_catalog
        obsTable['obs_track_id'] = obs_track_catalog
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
    parser.add_argument("--nsat", default=2, type=int)
    parser.add_argument("--nsite", default=None, type=int)
    parser.add_argument("--exposure_time_sec", default=30, type=float)
    parser.add_argument("--num_obs_per_track", default=10, type=int)
    parser.add_argument("--outdir", default='data/sim_for_iod', type=str)
    args = parser.parse_args()

    main(args)
