#!/usr/bin/env python
# encoding: utf-8
#
# Simulated tracks for IOD validation tests
# [modified from the sim_maneuver.py script]
#
# Design
# - Pick sats from GMM description of elements
# - Pick randomish values for area/mass/drag-coef/solar-rad-coef
# - Filter to LEO: perigee < 8060 km.
# - Observation sites are GEODSS + national labs
# - Propagate with ~realistic propagator for 24 hours
#   - 20,20 gravity
#   - sun/moon point sources
#   - drag from Harris-Priester atmospheric density model
#   - (constant) solar radiation pressure, accel is directly away from Sun
# - For each site, pick between 5 and 10 instances where alt is above 20 deg
#   - Observe every 10 s while satellite is up.
# - For each instance, record paired observations (ra/dec/time tuples) with
#   \Delta t = 3 seconds
#
# - Output tables
# - Sat table:
#   - satID
#   - TLE
#   - area/mass/CD/CR
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
import pickle

import numpy as np
from astropy import constants as const
from astropy.time import Time, TimeDelta
from astropy.table import Table
from astropy.io import fits
from astropy import units as u
from tqdm import tqdm
from textwrap import dedent

import ssapy
from ssapy.accel import AccelConstNTW

def ten(sixty):
    import re
    ds = 1
    d, m, s = [float(i) for i in re.split(':', sixty)]
    if str(d)[0] == '-':
        ds = -1
        d = abs(d)
    return ds*(d+m/60+s/3600)


def orbits_from_gmm(gmm, n, t0, leo=True):
    extrafac = 10
    log10a, log10perigee, ii = gmm.sample(int(n*10))[0].T
    aa = 10.**log10a # [m]
    perigee = 10.**log10perigee
    ee = 1-perigee/aa
    m = ((ee > 0) & (aa > ssapy.constants.WGS84_EARTH_RADIUS + 175) &
         (ii >= 0) & (ii <= np.pi))
    if leo:
        period = 2*np.pi * np.sqrt((aa**3)/(const.G.value*const.M_earth.value)) # [s]
        m = m & (period < 7200)
        # m = m & (perigee < 8060000)
    aa = aa[m][:n]
    ee = ee[m][:n]
    ii = ii[m][:n]
    raan, omega, phase = np.random.uniform(0, 2*np.pi, (3, n))
    return ssapy.orbit.Orbit.fromKeplerianElements(
        aa, ee, ii, omega, raan, phase, t0)


def orbits_from_file(filename):
    orbcat = Table.read(filename)
    return ssapy.orbit.Orbit.fromKeplerianElements(
        orbcat['a'], orbcat['e'], orbcat['i'], orbcat['omega'],
        orbcat['raan'], orbcat['phase'], orbcat['epoch'])



# lon, lat, height
sitelonlat = np.array(
    [[ 1.38081303e+02,  3.60786683e+01,  1.00024174e+02],
     [ 1.70465796e+02, -4.39866313e+01,  1.02898778e+03],
     [ 1.16334579e+02, -3.15140027e+01,  2.64973390e+02],
     [-7.07641333e+01, -3.04704471e+01,  1.57402603e+03],
     [-1.56268884e+02,  2.09090663e+01,  2.30996859e+02],
     [-1.19411818e+02,  3.70706047e+01,  2.30978019e+02],
     [ 2.99144082e+01, -2.30878942e+01,  9.22993878e+02],
     [ 2.59734056e+01,  3.50703962e+01,  3.90006008e+02],
     [ 1.51736426e+02, -2.67330801e+01,  5.43983008e+02],
     [-1.49188680e+02,  6.43002551e+01,  1.92013251e+02],
     [-6.82983235e+01,  7.65701971e+01,  4.22885872e+02],
     [-2.32599846e+00,  3.81657617e+01,  1.62199038e+03],
     [ 1.49637885e+01,  3.74655042e+01,  6.10014910e+01],
     [-1.25890725e+02,  4.91438380e+01,  9.98110168e+00],
     [ 1.17141082e+01,  5.09816675e+01,  3.31000164e+02],
     [-5.46080512e+01,  4.89568654e+01,  1.27970881e+02],
     [-1.43783037e+00,  5.11444232e+01,  9.19923723e+01],
     [-7.85224609e+01,  3.80332540e+01,  2.63971420e+02]])


def main(args):
    sites = []
    for i, (lon, lat, elev) in enumerate(sitelonlat):
        site = ssapy.EarthObserver(lon, lat, elev, fast=True)
        site._location
        site.name = '%d' % i
        site.id = i
        sites.append(site)

    if args.nsite is not None:
        sites = sites[:args.nsite]

    satTable = Table(
        names=[
            'satID', 'r', 'v', 't',
            'area', 'mass', 'CD', 'CR',
        ],
        dtype=["i4", "3f8","3f8","f8","f8","f8","f8","f8"]
    )

    if args.tstart is None:
        t0 = Time('2021-01-01T00:00:00').gps
    else:
        t0 = Time(tstart).gps
    np.random.seed(57722)
    if args.orbfile is not None:
        orbits = orbits_from_file(args.orbfile)
        if args.nsat > len(orbits):
            args.nsat = len(orbits)
            print(f'Only simulating {args.nsat} orbits')
    elif args.gmm is not None:
        gmm = pickle.load(open(args.gmm, 'rb'))
        orbits = orbits_from_gmm(gmm, args.nsat, t0)
    else:
        raise ValueError('Must set either gmm or orbfile.')

    # Use a ~realistic propagator.
    accel = (
        ssapy.AccelHarmonic(20, 20) + ssapy.AccelMoon() + ssapy.AccelSun()
#        + ssapy.AccelDrag()
# drag computation is hanging sometimes; commenting out.
        + ssapy.AccelSolRad()
    )

    track_id_ndx = 1

    for i in range(args.nsat):
        print("===== Simulating satellite {:d} / {:d} =====".format(
            i+1, args.nsat))

        # Pick some randomish values for mass, area, drag-coefficient,
        # reflection coefficient.
        mass = np.exp(np.random.uniform(np.log(1.5e0), np.log(1.5e4)))
        area = (mass)**(2/3) * np.random.uniform(0.01, 0.1)
        CD = np.random.uniform(2.0, 2.3)
        CR = np.random.uniform(1.2, 1.7)

        orbit = orbits[i]

        if args.tstart is None:
            tstart = t0 + np.random.uniform(0, 86400)
        else:
            tstart = t0

        if args.dvamp > 0:
            accelwmaneuver = accel
            expectednkick = args.ndvperorbit*args.nday*86400/orbit.period
            nkick = np.random.poisson(expectednkick)
            kicktimes = np.random.uniform(tstart, tstart+args.nday*86400, nkick)
            kickamps = args.dvamp * np.random.randn(nkick)
            burnlength = 10  # 10 s burns
            for j in range(nkick):
                direction = ssapy.utils.normed(np.random.randn(3))
                accelwmaneuver += AccelConstNTW(
                    direction*kickamps[j]/burnlength,
                    [kicktimes[j]-burnlength/2, kicktimes[j]+burnlength/2])
            # this set up does not try to avoid having satellites simultaneously
            # accelerating in two different directions.  That won't happen much
            # since 10 second burns are very short relative to a period,
            # but it may happen.
        else:
            accelwmaneuver = accel

        # Use Kepler propagator for faster development and testing
        # propagator = ssapy.KeplerianPropagator()
        ode_kwargs = dict(method='DOP853',
                          rtol=1e-9, atol=(1e-1, 1e-1, 1e-1, 1e-4, 1e-4, 1e-4))
        propagator = ssapy.propagator.SciPyPropagator(
            accelwmaneuver, ode_kwargs=ode_kwargs)
        propkw = dict(mass=mass, area=area, CD=CD, CR=CR)

        orbit.propkw = propkw

        passes = ssapy.compute.find_passes(orbit, sites, tstart, 86400*args.nday,
                                           dt=30, propagator=propagator)
        # defaults to Keplerian propagator, accepts propagator argument.
        obswindows = []
        for observer in passes:
            for time in passes[observer]:
                passdat = ssapy.compute.refine_pass(orbit, observer, time,
                                                    propagator=propagator)
                # defaults to Keplerian propagator, accepts propagator argument.
                passdat['observer'] = observer
                obswindows.append(passdat)

        obsdtype = [
            ('satID', 'i4'), ('t', 'f8'), ('r', '3f8'), ('v', '3f8'),
            ('ra', 'f8'), ('dec', 'f8'), ('slant', 'f8'),
            ('pmra', 'f8'), ('pmdec', 'f8'), ('slantrate', 'f8'),
            ('rStation', '3f8'), ('vStation', '3f8'), ('sensor', 'i4'),
            ('illum', 'bool'), ('sunalt', 'f8'), ('lon', 'f8'), ('lat', 'f8'),
            ('elev', 'f8')]
        allobs = []

        # observe satellite every 5 s in each window.
        minwindow = 0.2*u.min
        nsec = 5
        for window in obswindows:
            if window['duration'] < minwindow:
                continue
            tstart = (window['tStart'] + TimeDelta(
                np.random.uniform(0, 1)*minwindow/2)).gps
            tend = (window['tEnd'] - TimeDelta(
                np.random.uniform(0, 1)*minwindow/2)).gps
            tobs = np.arange(tstart, tend, nsec)
            rr, vv = ssapy.compute.rv(orbit, tobs, propagator=propagator)
            ra, dec, slant, pmra, pmdec, slantrate = (
                ssapy.compute.radec(orbit, tobs, observer=window['observer'],
                                  propagator=propagator, rate=True))
            rStation, vStation = window['observer'].getRV(tobs)
            obs = np.zeros(len(tobs), dtype=obsdtype)
            obs['t'] = tobs
            obs['r'] = rr
            obs['v'] = vv
            obs['ra'] = np.degrees(ra)
            obs['dec'] = np.degrees(dec)
            obs['slant'] = slant
            obs['pmra'] = np.degrees(pmra)
            obs['pmdec'] = np.degrees(pmdec)
            obs['slantrate'] = slantrate
            obs['satID'] = i+1
            obs['rStation'] = rStation
            obs['vStation'] = vStation
            obs['sensor'] = window['observer'].id
            obs['lon'] = window['observer'].lon
            obs['lat'] = window['observer'].lat
            obs['elev'] = window['observer'].elevation
            obs['sunalt'] = [np.degrees(window['observer'].sunAlt(t))
                             for t in tobs]
            if window['illumAtStart'] == window['illumAtEnd']:
                obs['illum'] = window['illumAtStart']
            elif window['illumAtStart']:
                obs['illum'] = tobs < window['tTerminator'].gps
            else:
                obs['illum'] = tobs > window['tTerminator'].gps
            allobs.append(obs)
        allobs = np.concatenate(allobs)
        satfilename = f'sat-obs-{i+1:04d}.fits'
        hdul = fits.HDUList()
        hduobs = fits.BinTableHDU(allobs)
        hduobs.header['mass'] = mass
        hduobs.header['area'] = area
        hduobs.header['CR'] = CR
        hduobs.header['CD'] = CD
        hduobs.name = 'OBS'

        tobs = np.arange(tstart, tstart+86400*args.nday, 5)
        rv = ssapy.compute.rv(orbit, tobs, propagator=propagator)
        rvtab = np.zeros(len(tobs), dtype=[('r', '3f8'), ('v', '3f8'),
                                           ('t', 'f8')])
        rvtab['r'] = rv[0]
        rvtab['v'] = rv[1]
        rvtab['t'] = tobs
        hdurvt = fits.BinTableHDU(rvtab)
        hdurvt.name = 'RV'
        hdul.append(hduobs)
        hdul.append(hdurvt)
        
        if not os.path.exists(args.datapath):
            os.makedirs(args.datapath)
        filepath = os.path.join(args.datapath, satfilename)
        hdul.writeto(filepath, overwrite=args.overwrite)


if __name__ == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("--nsat", default=2, type=int)
    parser.add_argument("--nsite", default=None, type=int)
    parser.add_argument("--gmm", default='branches/track_images/gmm100.pkl', type=str,
                        help='filename of GMM describing satellite distribution')
    parser.add_argument("--datapath", default='branches/track_images/track_obs/', type=str,
                        help='Path all track .fits files will be stored')               
    parser.add_argument("--orbfile", default=None, type=str)
    parser.add_argument("--overwrite", default=False, action='store_true')
    parser.add_argument("--dvamp", default=0, type=float,
                        help='Amplitude of velocity kicks in m/s')
    parser.add_argument("--ndvperorbit", default=0, type=int,
                        help='Number of velocity kicks per orbit')
    parser.add_argument("--nday", default=1, help='Number of days to simulate.')
    parser.add_argument('--tstart', default=None,
                        help='Time to start simulation.  If none, every sat '
                        'will start at a different time.')
    args = parser.parse_args()

    main(args)
