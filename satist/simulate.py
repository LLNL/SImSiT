"""
Generate synthetic telescope observations with stars and satellites.
"""

import os

import astroplan
import astropy.io.fits as fits
import astropy.units as u
import galsim
import glob
import numpy as np
import ssapy
import yaml
from astropy.time import Time
from tqdm import tqdm
import random

import satist as xfiles


def make_image(config, rng, parameters):
    """Generate a single simulated telescope image with stars and satellites.
    
    Creates a synthetic image including star field, satellite streaks, realistic
    noise, and vignetting. Generates both the challenge image with perturbed WCS
    and truth data including accurate positions and magnitudes.
    
    Parameters
    ----------
    config : dict
        Configuration dictionary containing:
            - n_sat : int
                Number of satellites per image
            - catalog : dict
                Catalog configuration with fields:
                    - develop : bool (optional)
                        Use mock catalog if True
                    - gaia_dir : str
                        Directory of GAIA data
                    - min_snr : float
                        Minimum SNR for stars to include
            - sat : dict
                Satellite configuration
            - tracker : dict
                Tracker configuration
    rng : numpy.random.Generator
        Random number generator
    parameters : dict
        Observation parameters containing:
            - observer : ssapy.EarthObserver
            - t0 : astropy.time.Time
            - instrument : Instrument
            - exptime : float
            - orbit : list of ssapy.Orbit
            - tracker : Tracker
            - zp : float
            - sat_mag : list of float
            - sky_sb : float
            - psf_fwhm : float
            - i_obs : int
            - propagator : ssapy.Propagator
    
    Returns
    -------
    hdu : astropy.io.fits.PrimaryHDU
        FITS HDU containing challenge image with perturbed WCS
    hduTrueWCS : astropy.io.fits.PrimaryHDU
        FITS HDU with true WCS at start of exposure
    hduTrueWCST : astropy.io.fits.PrimaryHDU
        FITS HDU with true WCS at end of exposure
    sample_docs : dict
        Sample submission format with perturbed positions and fluxes
    truth_hdulist : list
        List of FITS HDUs containing true satellite and star catalogs
    """
    observer = parameters['observer']
    t0 = parameters['t0']
    instrument = parameters['instrument']
    exptime = parameters['exptime']
    orbits = parameters['orbit']
    tracker = parameters['tracker']
    zp = parameters['zp']
    sat_mags = parameters['sat_mag']
    sky_sb = parameters['sky_sb']
    sky_phot = 10**(-0.4*(sky_sb - zp))
    psf_fwhm = parameters['psf_fwhm']
    psf = galsim.Kolmogorov(fwhm=psf_fwhm, 
                            gsparams=galsim.GSParams(folding_threshold=1e-3))
    i_obs = parameters['i_obs']
    propagator = parameters['propagator']
    nsat = config['n_sat']


    ####################################################################
    # Setup tracking
    ####################################################################

    # Because we're using an imperfect tracker, our boresight at t=t0 is
    # no longer at boresight0.  Let's fix this.
    boresight0 = tracker.get_boresight(t0)
    # rot_sky_pos0 = config['tracker']['rot_sky_pos0'] * galsim.degrees
    rot_sky_pos0 = tracker.get_rot_sky_pos(t0)
    wcs0 = instrument.get_wcs(
        boresight=boresight0, rot_sky_pos=rot_sky_pos0
    )

    print(f"{t0.iso = }")
    print(f"{boresight0.ra = }")
    print(f"{boresight0.dec = }")
    print(f"{rot_sky_pos0 = }")
    for orbit in orbits:
        print(f"{orbit.period = } s")
        print(f"{orbit.e = }")
        print(f"{orbit.i = } rad")

    ####################################################################
    # Get star catalog and trim faint sources
    ####################################################################

    if config['catalog'].get("develop", False):
        catalog = xfiles.MockStarCatalog()
    else:
        catalog = xfiles.GaiaStarCatalog(config['catalog']['gaia_dir'])
    cushion = 0.05 * galsim.degrees
    radius = instrument.field_radius + cushion
    p0 = tracker.get_boresight(t0)
    p1 = tracker.get_boresight(t0+exptime*u.s)
    stars = catalog.get_stars(p0, p1, radius, t0)
    stars['nphot'] = 10**(-0.4*(stars['i_mag'] - zp)) * exptime
    length = p0.distanceTo(p1) / galsim.arcsec / instrument.pix_size
    stars['SNR'] = instrument.streak_snr(
        nphot=stars['nphot'],
        length=length,
        psf_fwhm=psf_fwhm,
        sky_phot=sky_phot
    )
    stars = stars[stars['SNR'] > config['catalog']['min_snr']]
    print(f"{len(stars) = }")

    ####################################################################
    # Generate a mag/nphot for satellite
    ####################################################################

    sat_nphots = []
    for sat_mag in sat_mags:
        if sat_mag is not None:
            sat_nphots.append(10**(-0.4*(sat_mag - zp)) * exptime)
        else:
            sat_nphots.append(0)
    print(f"{sat_mags = }")
    print(f"{sat_nphots = }")

    ####################################################################
    # Run the simulation!
    ####################################################################

    image = instrument.init_image(sky_phot=sky_phot, exptime=exptime)
    stars = xfiles.tools.draw_stars(
        stars,
        t0=t0, exptime=exptime,
        wcs0=wcs0, tracker=tracker,
        psf=psf, image=image
    )
    sats, wcst = xfiles.tools.draw_sat(
        orbits,
        t0=t0, exptime=exptime,
        wcs0=wcs0, tracker=tracker,
        psf=psf, image=image,
        observer=observer,
        nphot=sat_nphots,
        propagator=propagator
    )

    # TODO: change here
    sats['i_mag'] = sat_mags
    instrument.apply_vignetting(image)
    instrument.apply_noise(image, sky_phot, rng)

    ####################################################################
    # Format output sample submission files
    ####################################################################

    sample_star_dicts = []
    star_indices = np.argsort(stars['i_mag'])
    sigma_pix = 5.0  # lie by this amount
    sigma_sky = np.deg2rad(5.0/3600) # 5 arcsec -> rad
    for i in star_indices[:10]:
        star = stars[i]
        sample_star_dicts.append({
            'x':float(star['x_FITS'] + rng.uniform(0, sigma_pix)),
            'y':float(star['y_FITS'] + rng.uniform(0, sigma_pix)),
            'flux':float(star['nphot'] * rng.uniform(1, 1.1)),
            'ra': float(star['ra'] + rng.uniform(0, sigma_sky) % 360),
            'dec':float(star['dec'] + rng.uniform(0, sigma_sky)),
            'mag':float(star['i_mag'] + rng.uniform(0, 0.1))
        })

    # Add satellites
    sample_sat_dicts = []
    for sat in sats:
        sample_sat_dicts.append({
            'x0':float(sat['x0_FITS'] + rng.uniform(0, sigma_pix)),
            'y0':float(sat['y0_FITS'] + rng.uniform(0, sigma_pix)),
            'x1':float(sat['x1_FITS'] + rng.uniform(0, sigma_pix)),
            'y1':float(sat['y1_FITS'] + rng.uniform(0, sigma_pix)),
            'flux':float(sat['nphot']) * rng.uniform(1, 1.1),
            'ra0': float(sat['ra0'] + rng.uniform(0, sigma_pix) % 360),
            'dec0':float(sat['dec0'] + rng.uniform(0, sigma_pix)),
            'ra1': float(sat['ra1'] + rng.uniform(0, sigma_pix) % 360),
            'dec1':float(sat['dec1'] + rng.uniform(0, sigma_pix)),
            'mag':float(sat['i_mag']),
        })

    sample_docs = {
        'file':f"{i_obs+1:04d}.fits",
        'sats':sample_sat_dicts,
        'stars':sample_star_dicts
    }

    ####################################################################
    # Store truth catalogs and write out true WCS
    ####################################################################

    star_hdu = fits.table_to_hdu(stars)
    star_hdu.name = f"STAR_{i_obs+1:04d}"
    sat_hdu = fits.table_to_hdu(sats)
    sat_hdu.name = f"SAT_{i_obs+1:04d}"

    truth_hdulist = [sat_hdu, star_hdu]

    # Write True WCS to private section
    hduTrueWCS = fits.PrimaryHDU()
    hduTrueWCST = fits.PrimaryHDU()
    hduTrueWCS.header['TIMESYS']='TAI'
    hduTrueWCS.header['DATE-BEG']=t0.tai.isot
    hduTrueWCS.header['DATE-END']=(t0+exptime*u.s).tai.isot
    hduTrueWCS.header['TELAPSE']=exptime
    hduTrueWCS.header['TIMEUNIT']='s'
    obsx = observer._location.x.to(u.m).value
    obsy = observer._location.y.to(u.m).value
    obsz = observer._location.z.to(u.m).value
    hduTrueWCS.header['OBSGEO-X']=(obsx, 'ITRS (m)')
    hduTrueWCS.header['OBSGEO-Y']=(obsy, 'ITRS (m)')
    hduTrueWCS.header['OBSGEO-Z']=(obsz, 'ITRS (m)')
    wcs0.writeToFitsHeader(hduTrueWCS.header, image.bounds)
    wcst.writeToFitsHeader(hduTrueWCST.header, image.bounds)

    ####################################################################
    # Write out challenge image with mangled WCS
    ####################################################################

    # Delete SIP and mangle the non-SIP part of true WCS to get
    # something suitable for the challenge.
    hdu = fits.PrimaryHDU(image.array)
    hdu.header['TIMESYS']='TAI'
    hdu.header['DATE-BEG']=t0.tai.isot
    hdu.header['DATE-END']=(t0+exptime*u.s).tai.isot
    hdu.header['TELAPSE']=exptime
    hdu.header['TIMEUNIT']='s'
    hdu.header['OBSGEO-X']=(obsx, 'ITRS (m)')
    hdu.header['OBSGEO-Y']=(obsy, 'ITRS (m)')
    hdu.header['OBSGEO-Z']=(obsz, 'ITRS (m)')
    hdu.header['MAGZP'] = zp + rng.uniform(-0.1, 0.1)

    htmp = hduTrueWCS.header
    crpix = htmp['CRPIX1'], htmp['CRPIX2']
    cd = np.array([
        [htmp['CD1_1'], htmp['CD1_2']],
        [htmp['CD2_1'], htmp['CD2_2']]
    ])
    crval1, crval2 = htmp['CRVAL1'], htmp['CRVAL2']
    shift_scale = 1/60.  # 1 arcmin shift
    crval1 += rng.uniform(-0.5*0.5)*shift_scale
    crval2 += rng.uniform(-0.5*0.5)*shift_scale
    rot_shift = np.deg2rad(rng.uniform(-0.5, 0.5))  # 1 degree rotation
    cr, sr = np.cos(rot_shift), np.sin(rot_shift)
    cd = cd @ np.array([[cr, sr], [-sr, cr]])

    htmpt = hduTrueWCST.header
    crpixt = htmpt['CRPIX1'], htmpt['CRPIX2']
    cdt = np.array([
        [htmpt['CD1_1'], htmpt['CD1_2']],
        [htmpt['CD2_1'], htmpt['CD2_2']]
    ])
    crval1t, crval2t = htmpt['CRVAL1'], htmpt['CRVAL2']
    shift_scale = 1/60.  # 1 arcmin shift
    crval1t += rng.uniform(-0.5*0.5)*shift_scale
    crval2t += rng.uniform(-0.5*0.5)*shift_scale
    rot_shiftt = np.deg2rad(rng.uniform(-0.5, 0.5))  # 1 degree rotation
    crt, srt = np.cos(rot_shift), np.sin(rot_shift)
    cdt = cdt @ np.array([[crt, srt], [-srt, crt]])

    hdu.header['CTYPE1'] = 'RA---TAN'
    hdu.header['CTYPE2'] = 'DEC--TAN'
    hdu.header['CRPIX1'], hdu.header['CRPIX2'] = crpix
    hdu.header['CD1_1'] = cd[0,0]
    hdu.header['CD1_2'] = cd[0,1]
    hdu.header['CD2_1'] = cd[1,0]
    hdu.header['CD2_2'] = cd[1,1]
    hdu.header['CUNIT1'] = 'deg'
    hdu.header['CUNIT2'] = 'deg'
    hdu.header['CRVAL1'] = crval1
    hdu.header['CRVAL2'] = crval2

    return hdu, hduTrueWCS, hduTrueWCST, sample_docs, truth_hdulist


def simulate(config):
    """Run full simulation to generate synthetic telescope observations.
    
    Generates a dataset of simulated telescope images with stars and satellites,
    including both public challenge data (with perturbed WCS) and private truth
    data. Each observation varies in observing conditions, telescope pointing,
    and satellite configurations.
    
    Parameters
    ----------
    config : dict
        Configuration dictionary containing:
            - outdir : str
                Output directory for generated files
            - instrument : dict
                Instrument configuration
            - n_obs : int
                Number of observations to generate
            - n_demo : int
                Number of demo observations to include in public data
            - seed : int
                Random seed for reproducibility
            - cadence : dict
                Observation cadence with 'exptime' list
            - sites : list
                List of observing site names
            - conditions : dict
                Observing conditions ranges:
                    - psf_fwhm_range : tuple
                    - sky_range : tuple
                    - zp_range : tuple
            - sat : dict
                Satellite configuration with 'mag_range'
            - n_sat : int
                Number of satellites per image (1 or 2)
            - tracker : dict
                Tracker configuration
            - catalog : dict
                Star catalog configuration
            - meta : dict
                Metadata to include in output
    
    Notes
    -----
    Creates the following directory structure:
        outdir/
            public/
                ####.fits (challenge images)
                sample_submission_N.yaml (demo sample)
                truth_N.fits (demo truth catalog)
                sky_flat.fits (optional, via make_sky_flat)
            private/
                ####.wcs.fits (true WCS at t0)
                ####.wcst.fits (true WCS at t1)
                sample_submission.yaml (complete truth)
                truth.fits (complete truth catalog)
    """
    for d in ["public", "private"]:
        dir_ = os.path.join(config['outdir'], d)
        if not os.path.exists(dir_):
            os.makedirs(dir_)

    instrument = xfiles.Instrument.fromConfig(config['instrument'])

    sample_docs = [config['meta']]
    truth_demo_hdul = fits.HDUList()
    truth_hdul = fits.HDUList()
    propagator = ssapy.KeplerianPropagator()

    for i_obs in tqdm(range(config['n_obs'])):
        with xfiles.tools.nostdout():
            # Use seed + i_obs so we could theoretically start in the middle.
            rng = np.random.default_rng(config['seed']+i_obs)

            ####################################################################
            # Setup the observing configuration
            ####################################################################

            exptime = rng.choice(config['cadence']['exptime'])
            site = astroplan.Observer.at_site(rng.choice(config['sites']))
            observer = ssapy.EarthObserver(
                lon=site.location.lon.to(u.deg).value,
                lat=site.location.lat.to(u.deg).value,
                elevation=site.location.height.to(u.m).value
            )
            cond_cfg = config['conditions']
            psf_fwhm = rng.uniform(*cond_cfg['psf_fwhm_range'])
            sky_sb = rng.uniform(*cond_cfg['sky_range'])
            zp_offset = rng.uniform(*cond_cfg['zp_range'])
            zp = instrument.compute_LSST_scaled_zp() + zp_offset
            print()
            print(f"{exptime = :.2f} s")
            print(f"{site.name = }")
            print(f"{psf_fwhm = :.2f} arcsec")
            print(f"{sky_sb = :.2f} mag / arcsec^2")
            print(f"{zp_offset = :.2f} mag")
            print(f"{zp = :.2f} mag")
            print(f"{sky_sb = :.2f} mag / arcsec^2")

            ####################################################################
            # Pick a time when it's night out, pick a direction to point,
            # and generate a LEO satellite nearby.
            ####################################################################

            t_day = Time("2010-01-01") + rng.uniform(0, 365)*u.d
            t0 = xfiles.tools.random_dark_time(t_day=t_day, site=site, rng=rng)

            boresight0 = xfiles.tools.random_boresight(
                observer=observer, t0=t0, horizon=np.deg2rad(20), rng=rng
            )
            sat_height = rng.uniform(400e3, 800e3)
            sat_coord = xfiles.tools.random_disk(
                boresight0, 0.2*instrument.field_radius, rng
            )
            heading = rng.uniform(0, 2*np.pi)
            vperp = rng.normal(7800, 10)  # nearly circular LEO
            vpar = rng.normal(10)

            tmid = t0 + 0.5*exptime*u.s
            orbit = xfiles.tools.generate_orbit(
                sat_height, sat_coord, heading, vperp, vpar, observer, tmid
            )

            if config['sat']['mag_range'] is not None:
                sat_mag = rng.uniform(*config['sat']['mag_range'])
            else:
                sat_mag = None

            if config['n_sat'] == 1:
                orbits = [orbit]
                sat_mags = [sat_mag]
            elif config['n_sat'] == 2:
                # the current values here are for CSOs (change back to *6 for CSOs)
                # np.random.beta(9.7*10**-6,2.4*10**-5)*1.45*10**-4)
                #ra_offset = random.choice((-1, 1)) * np.random.beta(2,5)*6 * 4.84814*10**-6
                #dec_offset = random.choice((-1, 1)) * np.random.beta(2,5)*6 * 4.84814*10**-6
                ra_offset = random.choice((-1, 1)) * psf_fwhm * random.uniform(.2, 1.5) # (.2, 1.5) 
                dec_offset = random.choice((-1, 1)) * psf_fwhm * random.uniform(.2, 1.5) # (.2, 1.5)
                print(f"{ra_offset = } arcsec")
                print(f"{dec_offset = } arcsec")
                orbit2 = xfiles.tools.generate_orbit(
                sat_height, galsim.CelestialCoord(
                    ra=(sat_coord.ra.rad+(ra_offset*4.84814*10**-6))*galsim.radians, 
                    dec=(sat_coord.dec.rad+(dec_offset*4.84814*10**-6))*galsim.radians), 
                    heading, vperp, vpar, observer, tmid
                )
                # magnitude of satellite 2 should be within 1.5 mags of satellite 1, not going over bounds
                if sat_mag - 1.5 < config['sat']['mag_range'][0]:
                    mag_bound_1 = -(sat_mag - config['sat']['mag_range'][0])
                else: 
                    mag_bound_1 = -1.5
                if sat_mag + 1.5 > config['sat']['mag_range'][1]:
                    mag_bound_2 = config['sat']['mag_range'][1] - sat_mag
                else: 
                    mag_bound_2 = 1.5
                
                mag_offset = random.uniform(mag_bound_1,mag_bound_2)
                sat_mag2 = sat_mag + mag_offset
                print(f"{mag_offset = } mag")
                orbits = [orbit, orbit2]
                sat_mags = [sat_mag, sat_mag2]
            else:
                print("Not a currently supported number of satellites per image (n_sat)")

            tracker_cfg = config['tracker']
            tracking_error = tracker_cfg['error']
            rot_sky_pos0 = config['tracker']['rot_sky_pos0'] * galsim.degrees
            if tracker_cfg['type'] == 'orbit':
                t_prev = t0 - tracking_error['t_rewind']*u.min
                r_prev, v_prev = ssapy.rv(orbit, t_prev, propagator=propagator)
                v_prev += rng.normal(
                    scale=tracking_error['v_perturb'],
                    size=3
                )
                tracking_orbit = ssapy.Orbit(r_prev, v_prev, t_prev)
                tracker = xfiles.OrbitTracker(
                    orbit=tracking_orbit, observer=observer, t0=t0,
                    rot_sky_pos0=rot_sky_pos0,
                    propagator=propagator
                )
            elif tracker_cfg['type'] == 'sidereal':
                tracking_boresight = xfiles.tools.random_disk(
                    boresight0,
                    instrument.field_radius*tracking_error['boresight'],
                    rng
                )
                tracker = xfiles.SiderealTracker(
                    tracking_boresight,
                    rot_sky_pos0
                )

            parameters = dict(
                orbit=orbits, observer=observer, t0=t0, instrument=instrument,
                exptime=exptime, tracker=tracker, sat_mag=sat_mags,
                sky_sb=sky_sb, psf_fwhm=psf_fwhm, zp=zp, i_obs=i_obs,
                propagator=propagator)

            hdu, hduTrueWCS, hduTrueWCST, sample_doc, true_hdulist = make_image(
                config, rng, parameters)
            truth_hdul.append(true_hdulist[0])
            truth_hdul.append(true_hdulist[1])
            if i_obs <= (config['n_demo']-1):
                truth_demo_hdul.append(true_hdulist[0])
                truth_demo_hdul.append(true_hdulist[1])

            sample_docs.append(sample_doc)

            hdu.writeto(
                os.path.join(
                    config['outdir'],
                    "public",
                    f"{i_obs+1:04d}.fits"
                ),
                overwrite=True
            )
            
            hduTrueWCS.writeto(
                os.path.join(
                    config['outdir'],
                    "private",
                    f"{i_obs+1:04d}.wcs.fits"
                ),
                overwrite=True
            )

            hduTrueWCST.writeto(
                os.path.join(
                    config['outdir'],
                    "private",
                    f"{i_obs+1:04d}.wcst.fits"
                ),
                overwrite=True
            )

    # Write out complete set of sample submissions privately
    with open(
        os.path.join(
            config['outdir'],
            "private",
            "sample_submission.yaml"
        ),
        "w"
    ) as f:
        yaml.safe_dump_all(sample_docs, f)

    # Write out first n_demo sample submissions publicly
    publicsamplefn = os.path.join(
        config['outdir'],
        "public",
        f"sample_submission_{config['n_demo']}.yaml"
    )
    with open(publicsamplefn, "w") as f:
        yaml.safe_dump_all(sample_docs[:config['n_demo']+1], f)

    truth_hdul.writeto(
        os.path.join(
            config['outdir'],
            "private",
            "truth.fits"
        ),
        overwrite=True
    )
    truth_demo_hdul.writeto(
        os.path.join(
            config['outdir'],
            "public",
            f"truth_{config['n_demo']}.fits"
        ),
        overwrite=True
    )


def make_sky_flat(configfile):
    """Generate a median sky flat field calibration image.
    
    Creates a normalized flat field by computing the median of multiple images
    from the dataset. Uses up to 1000 randomly sampled images to create a robust
    flat field correction.
    
    Parameters
    ----------
    configfile : str
        Path to YAML configuration file containing:
            - outdir : str
                Output directory containing public data
            - meta : dict
                Metadata with 'branch' field
            - instrument : dict
                Instrument config with 'image_shape'
    
    Notes
    -----
    The output flat field is normalized such that:
        1. Each pixel is divided by the median of its flattened image
        2. The median across all images is computed
        3. Final normalization to mean = 1.0
    
    Output is written to: {outdir}/public/sky_flat.fits
    """
    config = yaml.safe_load(open(configfile, 'r'))
    data_path = os.path.join(config['outdir'], 'public')

    if (config['meta']['branch'] == 'sidereal_track') or (config['meta']['branch'] == 'target_track'):
        file_format = "????_???.fits"
    else:
        file_format = "????.fits"

    all_files = glob.glob(os.path.join(data_path, file_format))
    if len(all_files) >= 1000:
        files = random.sample(all_files, 1000)
    else:
        files = all_files
    all_images = []

    for file in tqdm(files, desc="Creating sky flat", leave=True):
        image = fits.getdata(file)
        image_flattened = image.flatten()
        normalized_image = image_flattened / np.median(image_flattened)
        all_images.append(normalized_image)

    median_image = np.median(np.vstack(all_images), axis=0)
    final_image = median_image / np.mean(median_image)
    sky_flat = np.reshape(final_image, (config['instrument']['image_shape'][1], config['instrument']['image_shape'][0]))

    hdu = fits.PrimaryHDU(sky_flat)
    hdu.writeto(os.path.join(config['outdir'], "public", "sky_flat.fits"), overwrite=True)