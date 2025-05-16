import yaml
import galsim
import numpy as np
from astropy.io import fits
from astropy import units as u
from astropy.time import Time
import ssapy
import xfiles
import xfiles.simulate
import os
import re


def generate_correlated_series(times, correlationtime,
                               mean, amplitude, noise, rng):
    if len(times) > 1000:
        raise ValueError('bailing out; this sounds expensive')
    dt = times[:, None] - times[None, :]
    covar = np.exp(-(dt/correlationtime)**2)
    covar *= amplitude**2
    covar += np.eye(covar.shape[0])*noise**2
    chol = np.linalg.cholesky(covar)
    xx = rng.normal(scale=1, size=len(times))
    yy = chol.dot(xx)
    yy += mean
    return yy


def mag_lambsph(a, r_sat, range_sat = 35786000, phi=0):
    """
    Apparent magnitude of a spherical satellite.

    @param a        albedo of satellite
    @param r_sat    radius of the sphaerical satellite [units: meters]
    @param range_sat    distance of the satellite from the observer [units: meters]
    @param phi      solar phase angle. opening angle between
                    sun-satellite-observer, i.e 0 when linear geometric order
                    is sun, earth, then satellite, and pi when order is sun,
                    satellite, then earth. [units: radians]

    @return Apparent visual magnitude of a satellite approximated as a Lambertian sphere
    """
    max_phi = np.pi
    if np.max(phi) >= max_phi:
        raise ValueError("Error: currently phi must be less than pi/2. Exiting.")
    flux = ( (2./3. * a * r_sat**2) / (np.pi * range_sat**2) *
            (np.sin(phi) + (np.pi - phi) * np.cos(phi))
            )
    # apparent visual magnitude of the Lambertian sphere
    # The value -26.74 is the apparent visual magnitude of the sun
    mag = -26.74 - 2.5 * np.log10(flux)
    return mag


def simulate_one_track_images(filename, configfn, truthallfile, nobs, npass=1, iobs=0, stride=1):
    config = yaml.safe_load(open(configfn, 'r'))
    rng = np.random.default_rng(config['seed'] + iobs)
    obs = fits.getdata(filename, 'OBS')
    rv = fits.getdata(filename, 'RV')
    # night time observations only, when the satellite is not in shadow.
    mgood = (obs['illum'] > 0) & (obs['sunalt'] < -10)
    obs = obs[mgood]

    if npass != 1:
        raise ValueError('only npass = 1 supported at present.')
    # choose an observation with at least 12 hr of simulation afterward
    m = np.max(rv['t']) - obs['t'] > 3600*12
    idx = rng.choice(np.flatnonzero(m))
    mkeep = ((obs['sensor'] == obs['sensor'][idx])
             & (np.abs(obs['t'] - obs['t'][idx]) < 20*60))
    # keep observations of this pass; defined as within 20 min of this
    # observation
    obs = obs[mkeep]

    # cut down on volume if stride is > 1
    obs = obs[::stride]

    exptime = config['cadence']['exptime'][-1]

    u_to_sd = 0.2887

    instrument = xfiles.Instrument.fromConfig(config['instrument'])
    zpnom = instrument.compute_LSST_scaled_zp()

    # okay, let's generate some conditions
    conddict = dict()
    crange = config['conditions']['psf_fwhm_range']
    conddict['fwhm'] = [
        (crange[0] + crange[1])/2,
        (crange[1] - crange[0])*u_to_sd,
        0.1]
    crange = config['conditions']['sky_range']
    conddict['sky_mag'] = [
        (crange[0] + crange[1])/2,
        (crange[1] - crange[0])*u_to_sd,
        0.03]
    crange = config['conditions']['zp_range']
    conddict['zp'] = [
        zpnom + (crange[0] + crange[1])/2,
        (crange[1] - crange[0])*u_to_sd,
        0.02]
    conddict['dx'] = [
        0, config['tracker']['jitter'], 1./60/60]
    conddict['dy'] = [
        0, config['tracker']['jitter'], 1./60/60]
    # by having rot ~ 0 always, are we essentially saying that
    # we have no rotator and are using an equatorial mount telescope?
    conddict['rot'] = [
        0, 1., 1./60./60.]
    for xstr in 'xyz':
        conddict[f'v_perturb_{xstr}'] = [
            0, config['tracker']['error']['v_perturb'],
            config['tracker']['error']['v_perturb'] / 100]

    conditions = np.zeros(len(obs), dtype=[
            (name, 'f4') for name in conddict])
    for sensor0 in np.unique(obs['sensor']):
        m = obs['sensor'] == sensor0
        for name in conddict:
            xx = generate_correlated_series(
                obs['t'][m], 6*60*60,
                conddict[name][0], conddict[name][1], conddict[name][2], rng)
            conditions[name][m] = xx
    conditions['fwhm'] = np.hypot(conditions['fwhm'], 0.5)
    # don't let seeing ever get much better than 0.5".

    size = np.random.beta(1, 3)*(config['sat']['size_range'][1] - config['sat']['size_range'][0]) + config['sat']['size_range'][0]
    albedo = rng.uniform(*config['sat']['albedo_range'])

    if nobs is not None:
        num_obs = int(nobs)
    else:
        num_obs = len(obs)

    for i in range(num_obs):
        print(f'obs {i+1} of {num_obs} for sat {iobs}')
        zp = conditions['zp'][i]
        sky = conditions['sky_mag'][i]
        observer = ssapy.EarthObserver(
                obs['lon'][i], obs['lat'][i], obs['elev'][i])
        t0 = Time(obs['t'][i], format='gps')
        orbit = ssapy.Orbit(obs['r'][i], obs['v'][i], t=t0)
        # this won't have the right propagator information!
        # but it will be the correct osculating orbit, and
        # we only need it to be right for the integration time
        # but we'll do the wrong thing if the maneuver is during
        # the exposure.
        robs, vobs = observer.getRV(t0)

        distance = np.sqrt(np.sum((obs['r'][i] - robs)**2))
        sunpos = ssapy.utils.sunPos(t0)
        solphaseangle = ssapy.compute.unitAngle3(
            ssapy.utils.normed(obs['r'][i] - robs),
            ssapy.utils.normed(obs['r'][i] - sunpos))
        sat_mag = mag_lambsph(albedo, size, distance, solphaseangle)
        print(f'mag: {sat_mag:.2f}, dist (km): {distance/1000:.0f}, '
              f'size (m): {size:.1f}, albedo: {albedo:.2f}')

        boresight = ssapy.utils.lb2unit(np.radians(obs['ra'][i]),
                                        np.radians(obs['dec'][i]))
        zhat = np.array([0, 0, 1])
        rAlpha = ssapy.utils.normed(np.cross(zhat, boresight))
        rDec = ssapy.utils.normed(np.cross(boresight, rAlpha))
        newboresight = ssapy.utils.unit2lb(
            boresight + conditions['dx'][i]*rAlpha + conditions['dy'][i]*rDec)
        newboresight = galsim.CelestialCoord(
            newboresight[0]*galsim.radians, newboresight[1]*galsim.radians)

        tracker_cfg = config['tracker']
        tracking_error = tracker_cfg['error']
        rot_sky_pos0 = conditions['rot'][i] * galsim.degrees
        if tracker_cfg['type'] == 'orbit':
            t_prev = t0 - tracking_error['t_rewind']*u.min
            r_prev, v_prev = ssapy.rv(orbit, t_prev)
            v_prev += np.array([
                conditions['v_perturb_x'][i],
                conditions['v_perturb_y'][i],
                conditions['v_perturb_z'][i]])
            tracking_orbit = ssapy.Orbit(r_prev, v_prev, t_prev)
            tracker = xfiles.OrbitTracker(
                orbit=tracking_orbit, observer=observer, t0=t0,
                rot_sky_pos0=rot_sky_pos0,
                propagator=ssapy.KeplerianPropagator()
                )
        elif tracker_cfg['type'] == 'sidereal':
            tracker = xfiles.SiderealTracker(
                newboresight,
                rot_sky_pos0
                )

        parameters = dict(
            orbit=orbit, observer=observer, t0=t0, instrument=instrument,
            exptime=exptime, tracker=tracker, sat_mag=sat_mag,
            sky_sb=sky, psf_fwhm=conditions['fwhm'][i], zp=zp, i_obs=i,
            propagator=ssapy.KeplerianPropagator())
        hdu, hduTrueWCS, hduTrueWCST, sample_doc, true_hdulist = xfiles.simulate.make_image(
            config, rng, parameters)

        hdu.header['SENSOR'] = obs['sensor'][i]
        data_path = os.path.join(config['outdir'], 'public')
        if not os.path.exists(data_path):
            os.makedirs(data_path)
        image_path = os.path.join(data_path, f'{iobs:04d}_{i+1:03d}.fits')
        hdu.writeto(image_path, overwrite=True)
        hduTrueWCS.writeto(os.path.join(config['outdir'],"private",f'{iobs:04d}_{i+1:03d}.wcs.fits'),overwrite=True)
        hduTrueWCST.writeto(os.path.join(config['outdir'],"private",f'{iobs:04d}_{i+1:03d}.wcst.fits'),overwrite=True)
        true_hdulist[0].name = f"SAT_{iobs:04d}_{i+1:03d}"
        true_hdulist[1].name = f"STAR_{iobs:04d}_{i+1:03d}"
        truthallfile.append(true_hdulist[0])
        truthallfile.append(true_hdulist[1])

    obs = np.array(obs)[i:i+1]
    import numpy.lib.recfunctions as rfn
    # rename some columns for compatibility with IOD scoring script
    obs = rfn.append_fields(obs, ['rx', 'ry', 'rz', 'vx', 'vy', 'vz'],
                            [obs['r'][:, 0], obs['r'][:, 1], obs['r'][:, 2],
                             obs['v'][:, 0], obs['v'][:, 1], obs['v'][:, 2]])
    return obs, truthallfile  # contains state at end of sim, needed for truth



def make_sample_submission(endstate, rvt):
    noise = np.random.randn(6)
    # just the state for now
    return dict(rx=float(endstate['r'][0, 0] + noise[0]),
                ry=float(endstate['r'][0, 1] + noise[1]),
                rz=float(endstate['r'][0, 2] + noise[2]),
                vx=float(endstate['v'][0, 0] + noise[3]),
                vy=float(endstate['v'][0, 1] + noise[4]),
                vz=float(endstate['v'][0, 2] + noise[5]),
                t=float(endstate['t'][0]))


def simulate_many_track_images(ephemfiles, configfile, npass, stride, nobs):
    config = yaml.safe_load(open(configfile, 'r'))
    truthhdul = fits.HDUList()
    truthallhdul = fits.HDUList()
    sample_docs = [dict(branch=config['meta']['branch'], 
                        competitor_name=config['meta']['competitor_name'],
                        display_true_name=config['meta']['display_true_name'])]
    truthpath = os.path.join(config['outdir'], 'private')
    if not os.path.exists(truthpath):
        os.makedirs(truthpath)
    truthfile = os.path.join(truthpath, 'truth.fits')
    truthallfile = os.path.join(truthpath, 'truth_all.fits')
    for ephemfile in ephemfiles:
        iobs = int(os.path.basename(ephemfile)[-9:-5].lstrip('0'))
        endstate, truthallhdul = simulate_one_track_images(ephemfile, configfile, truthallhdul, nobs, npass,
                                             iobs=iobs, stride=stride)
        rvt = fits.getdata(ephemfile, 'RV')
        statehdu = fits.BinTableHDU(endstate)
        statehdu.name = f'SAT_{iobs:04d}'
        truthhdul.append(statehdu)
        rvthdu = fits.BinTableHDU(rvt)
        rvthdu.name = f'SAT_RVT_{iobs:04d}'
        truthhdul.append(rvthdu)
        sample_sat_dicts = [make_sample_submission(endstate, rvt)]
        idx = int(re.split('-|\.', os.path.basename(ephemfile))[-2])
        # write in this format for compatibility with IOD scoring file
        sample_docs.append(dict(file=f'SAT_{idx:04d}',
                                IOD=sample_sat_dicts))
    truthallhdul.writeto(truthallfile, overwrite=True)
    truthhdul.writeto(truthfile, overwrite=True)
    ndemo = config['n_demo']
    truthsamplepath = os.path.join(config['outdir'], 'public',
                                   f'truth_{ndemo}.fits')
    truthhdul[:ndemo*2+1].writeto(truthsamplepath, overwrite=True)
    sampallpath = os.path.join(config['outdir'], 'private',
                               'sample_submission.yaml')
    yaml.safe_dump_all(sample_docs, open(sampallpath, 'w'))
    sampsomepath = os.path.join(config['outdir'], 'public',
                                f'sample_submission_{ndemo}.yaml')
    yaml.safe_dump_all(sample_docs[:ndemo+1], open(sampsomepath, 'w'))


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(
        description='Simulate images following satellite.')
    parser.add_argument('ephemfiles', type=str, nargs='+',
                        help='file names with sat ephemerides')
    parser.add_argument('--config', type=str,
                        help='config file name')
    parser.add_argument('--npass', type=int, default=1,
                        help='number of overhead passes')
    parser.add_argument('--stride', type=int, default=1,
                        help='stride to use when going through observations.')
    parser.add_argument('--nobs', type=str, default=None)
    args = parser.parse_args()
    simulate_many_track_images(args.ephemfiles, configfile=args.config,
                               npass=args.npass, stride=args.stride, nobs=args.nobs)
    xfiles.simulate.make_sky_flat(args.config)
