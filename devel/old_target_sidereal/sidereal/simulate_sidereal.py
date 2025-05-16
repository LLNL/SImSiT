import os
import numpy as np
import yaml
import galsim
from scipy.interpolate import interp1d
import astropy.units as u
import tqdm
import astroplan
from astropy.time import Time
import ssa
import astropy.io.fits as fits
import glob


def sky_brightness(config, rng):
    """
    Parameters
    ----------
    config : dict
        configuration from yaml file
    rng : galsim.BaseDeviate
        Random number generator

    Returns
    -------
    sky_phot : float
        Photons per second per sq arcsec
    sky_mag : float
        Mag per sq arcsec
    zp : float
        Zeropoint mag
    """

    ud = galsim.UniformDeviate(rng)

    # Scaling from simple estimates of Rubin Observatory sky brightness
    # Numbers from https://github.com/jmeyers314/LSST_ETC
    etendue = 319.  # m^2 deg^2
    fov = 9.6 # deg^2
    area = etendue / fov  # m^2

    # Add some variability to the background.  Always assume ssa site is worse
    # than Cerro Pachon.
    B = config['dark_sky_i_mag'] - ud()*config['sky_range']

    # zeropoint in photons per second arriving from a zeroth magnitude AB source
    s0 = area*1.249 # For LSST zeropoint
    # i band sky brightnesses in AB mag / arcsec^2.
    # from http://www.lsst.org/files/docs/gee_137.28.pdf
    # Rescale LSST area to ~1m scope area (with 0.1 fractional obscuration)
    area_ratio = (1**2-0.1**2)/(8.36**2*(1-0.61**2))
    zp = 24 - ud()*config['zp_range'] + 2.5*np.log10(s0*area_ratio)

    # Sky brightness in photons per arcsec^2 per second
    sbar = 10**(-(B - zp) / 2.5)

    print()
    print("Sky brightness")
    print(f"{sbar:.2f} sky photons per square arcsecond per second")
    print(f"{B:.2f} magnitudes per square arcsecond")
    print(f"{zp:.2f} zero point")
    return sbar, B, zp


def polyWCS(th, dthdr, world_origin, theta=0, n=10, order=3, verbose=False):
    """
    Make a WCS from distortion polynomial

    Parameters:
        th:  Field angles in degrees
        dthdr:  Radial plate scale in arcsec per pixel
        world_origin:  Origin of ra, dec
        theta: Rotation angle in radians
        n: number of control points to use
        order: order of SIP part of fitted WCS

    Returns:
        wcs
    """
    from scipy.integrate import quad, IntegrationWarning
    import warnings

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

    sth, cth = np.sin(theta), np.cos(theta)
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


def generate_time_and_observer(config, rng):
    """Generate an observation time and observer.

    Parameters
    ----------
    config : dict
        configuration from yaml file
    rng : galsim.BaseDeviate
        Random number generator

    Returns
    -------
    time : astropy.time.Time
    observer : ssa.EarthObserver
    """
    # Site and time.
    ud = galsim.UniformDeviate(rng)
    isite = int(np.floor(ud()*len(config['sites'])))
    site_name = config['sites'][isite]
    site = astroplan.Observer.at_site(site_name)
    time = Time("2010-01-01") + ud()*365*u.d
    sunset = site.sun_set_time(time, which='next')
    sunrise = site.sun_rise_time(sunset, which='next')
    time = ud()*(sunrise-sunset)+sunset
    observer = ssa.EarthObserver(
        lon=site.location.lon.to(u.deg).value,
        lat=site.location.lat.to(u.deg).value,
        elevation=site.location.height.to(u.m).value
    )
    print()
    print(f"site: {site_name}")
    print(f"sunset:  {sunset.iso}")
    print(f"sunrise: {sunrise.iso}")
    print(f"time:    {time.iso}")
    return time, observer


def generate_wcs(config, rng, time, observer):
    """Generate a WCS

    Parameters
    ----------
    config : dict
        configuration from yaml file
    rng : galsim.BaseDeviate
        Random number generator
    time : astropy.time.Time
    observer : ssa.EarthObserver

    Returns
    -------
    wcs : galsim.WCS
    """
    from astropy.coordinates import AltAz, GCRS
    ud = galsim.UniformDeviate(rng)

    if config['develop']:
        # Restrict ra/dec range so we don't need to download the whole Gaia
        # catalog during development.
        ra = np.deg2rad(ud())
        dec = np.deg2rad(ud())
    else:
        # grab random point above 20 degrees altitude for the boresight.
        # Cylinder is an equal area projection from a sphere, so can uniformly
        # generate on cylinder and then project back to the sphere.
        zmin = np.sin(np.deg2rad(20))
        z = (1-zmin)*ud() + zmin
        az = 2*np.pi*ud()
        alt = np.arctan(z/np.sqrt(1-z**2))

        sc = AltAz(
            az*u.rad,
            alt*u.rad,
            obstime=time,
            location=observer._location
        )
        radec = sc.transform_to(GCRS)
        ra = radec.ra.rad
        dec = radec.dec.rad

    boresight = galsim.CelestialCoord(ra*galsim.radians, dec*galsim.radians)
    print()
    print("Generating WCS")
    print(f"Boresight = {boresight}")
    # uniform random rotator angle
    rot = 2*np.pi*ud()
    wcs = polyWCS(
        config['distortion']['th'],
        np.array(config['distortion']['dthdr'])*config['pixel_scale'],
        boresight,
        rot,
        n=10, order=3
    )
    return wcs


def init_image(config, wcs, exp_time, sky_phot):
    """Initialize an image.

    Parameters
    ----------
    config : dict
        configuration from yaml file
    wcs : galsim.WCS
    exp_time : image exposure time
    sky_phot : float
        Background flux in photons per sq arcsec

    Returns
    -------
    image : galsim.Image
        Initialized image with background sky value filled in.
    """
    nx, ny = config['fov']
    bounds = galsim.BoundsI(-nx//2, nx//2-1, -ny//2, ny//2-1)
    image = galsim.Image(bounds, wcs=wcs)
    # Sky background
    wcs.makeSkyImage(image, sky_phot*exp_time)
    print()
    print("Initializing image")
    return image


def apply_vignetting(config, image):
    """Apply vignetting to an image.

    Parameters
    ----------
    config : dict
        configuration from yaml file
    image : galsim.Image
        Image to be modified in place.
    """
    vigfun = interp1d(
        config['vignetting']['th'],
        config['vignetting']['unvig'],
        kind='cubic'
    )
    # Vignette the background
    # Should really use wcs here to use angular distances, but for now we'll
    # cheat and just use pixel distances.
    xx, yy = image.get_pixel_centers()
    rr = np.hypot(xx, yy)  # dist from center in pixels
    rr *= np.sqrt(image.wcs.pixelArea(galsim.PositionD(0, 0)))  # -> arcsec
    image.array[:] *= vigfun(rr/3600)
    print()
    print("Applying vignetting to image")


def get_stars(config, image, time, observer):
    """Get star catalog from GAIA and convert from ICRF to apparent
    (accounting for proper motion, parallax, and aberration).

    Parameters
    ----------
    config : dict
        configuration from yaml file
    image : galsim.Image
        Image defining bounds and wcs used to search Gaia catalog.
    time : astropy.time.Time
    observer : ssa.EarthObserver

    Returns
    -------
    stars : astropy.table.Table
        Table with columns
            - ra
            - dec
            - i_flux
            - pm_ra
            - pm_dec
            - parallax
    """
    import esutil
    from astropy.table import Table, vstack
    from ssa.utils import catalog_to_apparent
    center = image.wcs.posToWorld(image.true_center)
    htm = esutil.htm.HTM(7)
    corners = [
        (image.bounds.xmin, image.bounds.ymin),
        (image.bounds.xmin, image.bounds.ymax),
        (image.bounds.xmax, image.bounds.ymin),
        (image.bounds.xmax, image.bounds.ymax)
    ]
    max_dist = 0*galsim.degrees
    for x, y in corners:
        sky = image.wcs.posToWorld(galsim.PositionD(x, y))
        dist = center.distanceTo(sky)
        if dist > max_dist:
            max_dist = dist
    radius = max_dist.deg + 0.1
    shards = htm.intersect(center.ra.deg, center.dec.deg, radius)
    gaia_dir = config['gaia_dir']

    table_list = []
    for shard in shards:
        file = f"{gaia_dir}/{shard}.fits"
        sdata = Table.read(file)
        data = Table()
        data['ra'], data['dec'] = catalog_to_apparent(
            sdata['coord_ra'].to(u.rad).value,
            sdata['coord_dec'].to(u.rad).value,
            time,
            observer=None,  # we don't want to include diurnal aberration
            pmra=sdata['pm_ra'].to(u.mas).value,
            pmdec=sdata['pm_dec'].to(u.mas).value,
            parallax=sdata['parallax'].to(u.arcsec).value,
        )
        try:
            data['x'], data['y'] = image.wcs.radecToxy(
                data['ra'], data['dec'], units='radians'
            )
        except:
            import ipdb; ipdb.set_trace()
        data['i_mag'] = sdata['phot_g_mean_flux'].to(u.ABmag).value

        mask = data['x'] >= image.bounds.xmin
        mask &= data['x'] <= image.bounds.xmax
        mask &= data['y'] >= image.bounds.ymin
        mask &= data['y'] <= image.bounds.ymax
        mask &= sdata['phot_g_mean_flux'] > 0.0
        table_list.append(data[mask])
    print()
    print("Fetching stars from GAIA catalog")
    return vstack(table_list)


def generate_psf(config, rng):
    """Generate PSF

    Parameters
    ----------
    config : dict
        configuration from yaml file
    rng : galsim.BaseDeviate
        Random number generator

    Returns
    -------
    psf : galsim.GSObject
    """
    ud = galsim.UniformDeviate(rng)
    fwhm = ud()
    fwhm *= (config['psf_fwhm'][1]-config['psf_fwhm'][0])
    fwhm += config['psf_fwhm'][0]
    psf = galsim.Kolmogorov(
        fwhm=fwhm,
        gsparams=galsim.GSParams(folding_threshold=1e-3)
    )
    print()
    print("Generating PSF")
    print(f"FWHM = {fwhm:.2f} arcsec")
    return psf


def draw_stars(image, stars, psf):
    """Draw stars onto image

    Parameters
    ----------
    image : galsim.Image
        Image defining bounds and wcs used to search Gaia catalog.
    stars : astropy.table.Table
        Output from get_stars()
    psf : galsim.GSObject
        PSF to use.
    """
    print()
    print("Drawing stars")
    for star in tqdm.tqdm(stars):
        x, y = star['x'], star['y']
        local_wcs = image.wcs.local(galsim.PositionD(x, y))
        stamp = (psf*star['nphot']).drawImage(
            wcs=local_wcs,
            center=(x, y),
        )
        bounds = stamp.bounds & image.bounds
        image[bounds] += stamp[bounds]


def get_streaks(config, exp_time, rng, image):
    """Get streak parameters

    Parameters
    ----------
    config : dict
        configuration from yaml file
    exp_time : image exposure time
    rng : galsim.BaseDeviate
        Random number generator
    image : galsim.Image
        Image defining bounds and wcs used to search Gaia catalog.

    Returns
    -------
    streaks : astropy.table.Table
        Table with columns:
            - ra1 : rad
            - dec1 : rad
            - ra2 : rad
            - dec2 : rad
            - mag : float
    """
    print()
    print("Generating streak")
    # For now, always generate exactly one streak
    from astropy.table import Table
    ud = galsim.UniformDeviate(rng)

    streak_mag = ud() * (config['streak_mag'][1] - config['streak_mag'][0])
    streak_mag += config['streak_mag'][0]

    streak_rate = ud() * (config['streak_rate'][1] - config['streak_rate'][0])
    streak_rate += config['streak_rate'][0]

    scale = np.sqrt(image.wcs.pixelArea(galsim.PositionD(0, 0)))  # arcsec/pix
    streak_length = streak_rate * exp_time / scale

    streak_angle = ud()*2*np.pi

    print(f"streak magnitude = {streak_mag:.2f}")
    print(f"streak length = {streak_length:.2f} pix")

    xc = (ud()-0.5)*config['position_box']
    yc = (ud()-0.5)*config['position_box']

    x0 = streak_length/2*np.cos(streak_angle) + xc
    y0 = streak_length/2*np.sin(streak_angle) + yc
    x1 = -streak_length/2*np.cos(streak_angle) + xc
    y1 = -streak_length/2*np.sin(streak_angle) + yc

    ra0, dec0 = image.wcs.xyToradec(x0, y0, units='radians')
    ra1, dec1 = image.wcs.xyToradec(x1, y1, units='radians')

    out = Table()
    out['x0'] = [x0]
    out['y0'] = [y0]
    out['x1'] = [x1]
    out['y1'] = [y1]
    out['ra0'] = [ra0]
    out['dec0'] = [dec0]
    out['ra1'] = [ra1]
    out['dec1'] = [dec1]
    out['mag'] = [streak_mag]
    return out


def draw_streaks(image, streaks, psf, n_split=10):
    """Draw satellite streak

    Parameters
    ----------
    image : galsim.Image
        Image defining bounds and wcs used to search Gaia catalog.
    streaks : astropy.table.Table
        Output from get_streaks
    psf : galsim.GSObject
        PSF to use.
    n_split : int
        Number of piecewise segments into which to split streak
    """
    print()
    print("Drawing streaks")
    ncp = galsim.CelestialCoord(0*galsim.degrees, 90*galsim.degrees)

    for istreak, streak in enumerate(streaks):
        point0 = galsim.CelestialCoord(
            streak['ra0']*galsim.radians,
            streak['dec0']*galsim.radians
        )
        point1 = galsim.CelestialCoord(
            streak['ra1']*galsim.radians,
            streak['dec1']*galsim.radians
        )
        nphot = streak['nphot']

        print(
            f"Drawing streak {istreak} of {len(streaks)}"
            f" with {nphot:.0f} photons"
        )

        dist = point0.distanceTo(point1)
        p0 = point0
        for i in tqdm.trange(1, n_split+1):
            p1 = point0.greatCirclePoint(point1, dist*i/n_split)
            q = p0.angleBetween(ncp, p1)
            pos0 = image.wcs.toImage(p0)
            pos1 = image.wcs.toImage(p1)
            x0, y0 = pos0.x, pos0.y
            x1, y1 = pos1.x, pos1.y
            xmid = (x0+x1)/2
            ymid = (y0+y1)/2
            local_wcs = image.wcs.local(galsim.PositionD(xmid, ymid))
            length = dist/n_split/galsim.arcsec
            box = galsim.Box(
                1e-12,
                length,
                flux=nphot/n_split
            ).rotate(q)
            obj = galsim.Convolve(box, psf)
            stamp = obj.drawImage(
                wcs=local_wcs,
                center=(xmid, ymid),
            )
            bounds = stamp.bounds & image.bounds
            if bounds.area():
                image[bounds] += stamp[bounds]
            p0 = p1


def make_image(config, exp_time, rng):
    """Make image

    Parameters
    ----------
    config : dict
        configuration from yaml file
    exp_time : image exposure time
    rng : galsim.BaseDeviate
        Random number generator

    Returns
    -------
    image : galsim.Image
    stars : astropy.table.Table
    streaks : astropy.table.Table
    time : astropy.time.Time
    observer : ssa.EarthObserver  
    zp : float
    """
    # Fiducial sky background for Cerro Pachon
    sky_phot, sky_mag, zp = sky_brightness(config, rng)

    time, observer = generate_time_and_observer(config, rng)
    wcs = generate_wcs(config, rng, time, observer)
    image = init_image(config, wcs, exp_time, sky_phot)
    pix_std = np.sqrt(np.mean(image.array))  # Assuming Poisson stats
    print(f"mean pix count = {np.mean(image.array):.2f}")

    psf = generate_psf(config, rng)

    streaks = get_streaks(config, exp_time, rng, image)
    streaks['nphot'] = 10**(-(streaks['mag'] - zp) / 2.5)
    streaks['nphot'] *= exp_time

    draw_streaks(image, streaks, psf)

    stars = get_stars(config, image, time, observer)
    stars['nphot'] = 10**(-(stars['i_mag'] - zp) / 2.5)
    stars['nphot'] *= exp_time

    # Remove stars with fewer photons than pixel standard deviation
    phot_thresh = pix_std * config['faint_star_skip']
    wkeep = stars['nphot'] > phot_thresh
    print(
        f"Removing {len(wkeep)-sum(wkeep)} of {len(stars)} stars "
        f"with fewer than {phot_thresh:.0f} photons"
    )
    stars = stars[wkeep]

    draw_stars(image, stars, psf)

    apply_vignetting(config, image)

    noise = galsim.CCDNoise(
        rng,
        sky_level=0,
        gain=config['gain'],
        read_noise=config['read_noise']
    )
    image.addNoise(noise)

    # x/y in stars and streaks currently set with origin at image center.
    # but for the truth catalog, we'll instead set the origin to the lower left.
    stars['x'] -= image.xmin-1
    stars['y'] -= image.ymin-1
    streaks['x0'] -= image.xmin-1
    streaks['y0'] -= image.ymin-1
    streaks['x1'] -= image.xmin-1
    streaks['y1'] -= image.ymin-1
    stars['ra'] = np.degrees(stars['ra'])
    stars['dec'] = np.degrees(stars['dec'])
    streaks['ra0'] = np.degrees(streaks['ra0'])
    streaks['ra1'] = np.degrees(streaks['ra1'])
    streaks['dec0'] = np.degrees(streaks['dec0'])
    streaks['dec1'] = np.degrees(streaks['dec1'])

    return image, stars, streaks, time, observer, zp

def make_sky_flat(config):
    """
    Generates a median sky flat .fits file

    Parameters
    ----------
    config : dict
        configuration from yaml file
    """

    files = glob.glob(os.path.join('public/', "???.fits"))
    all_images = []

    for file in tqdm.tqdm(files, desc="Creating sky flat", leave=True):
        image = fits.getdata(file)
        image_flattened = image.flatten()
        normalized_image = image_flattened / np.median(image_flattened)
        all_images.append(normalized_image)

    median_image = np.median(np.vstack(all_images), axis=0)
    final_image = median_image / np.mean(median_image)
    sky_flat = np.reshape(final_image, (config["fov"][1], config["fov"][0]))

    hdu = fits.PrimaryHDU(sky_flat)
    hdu.writeto(f"public/sky_flat.fits", overwrite=True)

def simulate_sidereal(args):
    for d in ["public", "private"]:
        if not os.path.exists(d):
            os.makedirs(d)

    config = yaml.safe_load(open(args.config, 'r'))
    rng = galsim.BaseDeviate(config['seed'])
    
    ud = galsim.UniformDeviate(rng)

    sample_docs = [{'branch' : config['branch'], 'competitor_name' : config['competitor_name'], 'display_true_name' : config['display_true_name']}]
    truth_demo_hdul = fits.HDUList()
    truth_hdul = fits.HDUList()

    for i_obs in tqdm.tqdm(range(config['n_obs'])):
    # for i_obs in range(config['n_obs']):
        print()
        print()
        print()
        print("------------------------------------------------")
        print(f"Creating image {i_obs+1} of {config['n_obs']}")
        print()
        exp_time = config['exp_time'][int(np.floor(ud()*len(config['exp_time'])))]
        img, stars, streaks, time, observer, zp = make_image(config, exp_time, rng)
        stars.meta = {}
        obsx = observer._location.x.to(u.m).value
        obsy = observer._location.y.to(u.m).value
        obsz = observer._location.z.to(u.m).value
        ud = galsim.UniformDeviate(rng)

        # Add 10 brightest stars to sample submission file
        sample_star_dicts = []

        star_indices = np.argsort(stars['i_mag'])
        sigma_pix = 5.0
        sigma_sky = np.deg2rad(5.0/3600) # 5 arcsec -> rad
        for i in star_indices[:10]:
            star = stars[i]
            sample_star_dicts.append({
                'x':float(star['x'] + ud()*sigma_pix),
                'y':float(star['y'] + ud()*sigma_pix),
                'flux':float(star['nphot'] * (1+0.1*ud())),
                'ra': float(np.rad2deg(star['ra']  + ud()*sigma_sky) % 360),
                'dec':float(np.rad2deg(star['dec'] + ud()*sigma_sky)),
                'mag':float(star['i_mag'] + 0.1*ud())
            })

        # Add satellites
        sample_sat_dicts = []
        for streak in streaks:
            sample_sat_dicts.append({
                'x0':float(streak['x0'] + ud()*sigma_pix),
                'y0':float(streak['y0'] + ud()*sigma_pix),
                'x1':float(streak['x1'] + ud()*sigma_pix),
                'y1':float(streak['y1'] + ud()*sigma_pix),
                'flux':float(streak['nphot']),
                'ra0': float(np.rad2deg(streak['ra0']  + ud()*sigma_sky) % 360),
                'dec0':float(np.rad2deg(streak['dec0'] + ud()*sigma_sky)),
                'ra1': float(np.rad2deg(streak['ra1']  + ud()*sigma_sky) % 360),
                'dec1':float(np.rad2deg(streak['dec1'] + ud()*sigma_sky)),
                'mag':float(streak['mag']),
            })

        sample_docs.append({
            'file':f"{i_obs:03d}.fits",
            'sats':sample_sat_dicts,
            'stars':sample_star_dicts
        })

        # Write out first n_demo sample submissions publicly
        if i_obs == (config['n_demo']-1):
            fn = "public/sample_submission_"+str(config["n_demo"])+"_sidereal.yaml"
            with open(fn, "w") as f:
                yaml.safe_dump_all(sample_docs, f)

        # Now for the truth catalogs
        star_hdu = fits.table_to_hdu(stars)
        star_hdu.name = f"STAR_{i_obs:03d}"
        sat_hdu = fits.table_to_hdu(streaks)
        sat_hdu.name = f"SAT_{i_obs:03d}"

        if i_obs <= (config['n_demo']-1):
            truth_demo_hdul.append(star_hdu)
            truth_demo_hdul.append(sat_hdu)
        truth_hdul.append(sat_hdu)
        truth_hdul.append(star_hdu)

        # Write True WCS to private section
        # hduTrueWCS = fits.PrimaryHDU(img.array)
        hduTrueWCS = fits.PrimaryHDU()
        hduTrueWCS.header['TIMESYS']='TAI'
        hduTrueWCS.header['DATE-BEG']=time.tai.isot
        hduTrueWCS.header['DATE-END']=(time+exp_time*u.s).tai.isot
        hduTrueWCS.header['TELAPSE']=exp_time
        hduTrueWCS.header['TIMEUNIT']='s'
        hduTrueWCS.header['OBSGEO-X']=(obsx, 'ITRS (m)')
        hduTrueWCS.header['OBSGEO-Y']=(obsy, 'ITRS (m)')
        hduTrueWCS.header['OBSGEO-Z']=(obsz, 'ITRS (m)')
        img.wcs.writeToFitsHeader(hduTrueWCS.header, img.bounds)
        hduTrueWCS.writeto(f"private/{i_obs:03d}.wcs.fits", overwrite=True)

        # Delete SIP and mangle a non-sip part of true WCS to get something
        # suitable for the challenge.
        hdu = fits.PrimaryHDU(img.array)
        hdu.header['TIMESYS']='TAI'
        hdu.header['DATE-BEG']=time.tai.isot
        hdu.header['DATE-END']=(time+exp_time*u.s).tai.isot
        hdu.header['TELAPSE']=exp_time
        hdu.header['TIMEUNIT']='s'
        hdu.header['OBSGEO-X']=(obsx, 'ITRS (m)')
        hdu.header['OBSGEO-Y']=(obsy, 'ITRS (m)')
        hdu.header['OBSGEO-Z']=(obsz, 'ITRS (m)')
        hdu.header['MAGZP'] = zp + (1-2*ud())*config['magzperr']

        htmp = hduTrueWCS.header
        crpix = htmp['CRPIX1'], htmp['CRPIX2']
        cd = np.array([
            [htmp['CD1_1'], htmp['CD1_2']],
            [htmp['CD2_1'], htmp['CD2_2']]
        ])
        crval = htmp['CRVAL1'], htmp['CRVAL2']
        shift_scale = 1/60.  # 1 arcmin shift
        crval += np.array([(ud()-0.5)*shift_scale, (ud()-0.5)*shift_scale])
        rot_shift = np.deg2rad(1*(ud()-0.5))  # 1 degree rotation
        cr, sr = np.cos(rot_shift), np.sin(rot_shift)
        cd = cd @ np.array([[cr, sr], [-sr, cr]])

        hdu.header['CTYPE1'] = 'RA---TAN'
        hdu.header['CTYPE2'] = 'DEC--TAN'
        hdu.header['CRPIX1'], hdu.header['CRPIX2'] = crpix
        hdu.header['CD1_1'] = cd[0,0]
        hdu.header['CD1_2'] = cd[0,1]
        hdu.header['CD2_1'] = cd[1,0]
        hdu.header['CD2_2'] = cd[1,1]
        hdu.header['CUNIT1'] = 'deg'
        hdu.header['CUNIT2'] = 'deg'
        hdu.header['CRVAL1'], hdu.header['CRVAL2'] = crval

        hdu.writeto(f"public/{i_obs:03d}.fits", overwrite=True)

    # Write out complete set of sample submissions privately
    fn = "private/sample_submission.yaml"
    with open(fn, "w") as f:
        yaml.safe_dump_all(sample_docs, f)

    truth_hdul.writeto("private/truth.fits", overwrite=True)
    truth_demo_hdul.writeto("public/truth_"+str(config["n_demo"])+"_sidereal.fits", overwrite=True)

    make_sky_flat(config)


if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument(
        "--config",
        type=str,
        default="simulate_sidereal_config.yaml"
    )
    args = parser.parse_args()

    simulate_sidereal(args)
