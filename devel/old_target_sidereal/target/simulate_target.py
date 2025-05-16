import yaml
import galsim
import numpy as np
from astropy.table import Table, vstack
import astropy.units as u
import astropy.io.fits as fits
import os
import tqdm


class CubicWCS(galsim.wcs.CelestialWCS):
    """WCS that includes a cubic radial term along with a 'gnomiic' projection to sky coordinates.

    This wcs goes from x,y -> u,v -> ra, dec

    The x,y to u,v transformation is determined by
        u = scale*x*(1 + r3*(x^2 + y^2)),
        v = scale*y*(1 + r3*(x^2 + y^2))

    Where scale will be the the pixel scale at the center of the image and r3 determines
    the amount of radial distortion.  For a value of r3=3e-8 the pixel scale changes
    by ~10% at a radius of ~ 1400.

    The inverse from u,v -> x,y is done via Cardano's method for solving cubic equations.

    Note if an image with this wcs is written to a file, galsim cannot currently read the
    file because all wcs classes are expected to be in the galsim namespace

    Parameters:
        scale:            The nominal pixel scale
        r3:               The amount of radial distortion
        origin:           Origin position for the image coordinate system.
        world_origin:     Origin position in ra,dec

    """

    def __init__(self, scale, r3, origin, world_origin):
        self._scale = scale
        self._r3 = r3
        self._color = None
        self._set_origin(origin)
        self._world_origin = world_origin

        self._q = 1./(self._r3*self._scale)
        self._p = 1./self._r3

        self._torad = galsim.arcsec / galsim.radians

    def _ufunc(self, x, y):
        return self._scale*x*(1. + self._r3 * (x**2 + y**2))

    def _vfunc(self, x, y):
        return self._scale*y*(1. + self._r3 * (x**2 + y**2))

    def _xfunc(self, u, v):
        wsq = u*u + v*v
        if wsq == 0.:
            return 0.
        else:
            w = np.sqrt(wsq)
            temp = (np.sqrt(self._q**2*wsq/4 + self._p**3/27))
            r = (temp + self._q*w/2)**(1./3) - (temp - self._q*w/2)**(1./3)
            return u*r/w

    def _yfunc(self, u, v):
        wsq = u*u + v*v
        if wsq == 0.:
            return 0.
        else:
            w = np.sqrt(wsq)
            temp = (np.sqrt(self._q**2*wsq/4 + self._p**3/27))
            r = (temp + self._q*w/2)**(1./3) - (temp - self._q*w/2)**(1./3)
            return v*r/w

    def _radec_func(self, x, y):
        return self._world_origin.deproject_rad(self._ufunc(x, y)*self._torad,
                                                self._vfunc(x, y)*self._torad)

    def _xy_func(self, ra, dec):
        u, v = self._world_origin.project_rad(ra, dec)

        u /= self._torad
        v /= self._torad
        x, y = (self._xfunc(u, v), self._yfunc(u, v))
        return (x, y)

    @property
    def origin(self):
        """The input radec_func
        """
        return self._origin

    @property
    def world_origin(self):
        """The input radec_func
        """
        return self._world_origin

    @property
    def radec_func(self):
        """The input radec_func
        """
        return self._radec_func

    @property
    def xy_func(self):
        """The input xy_func
        """
        return self._xy_func

    def _radec(self, x, y, color=None):
        try:
            return self._radec_func(x, y)
        except Exception as e:
            try:
                world = [self._radec(x1, y1) for (x1, y1) in zip(x, y)]
            except Exception:  # pragma: no cover
                # Raise the original one if this fails, since it's probably more relevant.
                raise e
            ra = np.array([w[0] for w in world])
            dec = np.array([w[1] for w in world])
            return ra, dec

    def _xy(self, ra, dec, color=None):
        try:
            return self._xy_func(ra, dec)
        except Exception as e:
            try:
                image = [self._xy(ra1, dec1) for (ra1, dec1) in zip(ra, dec)]
            except Exception:  # pragma: no cover
                # Raise the original one if this fails, since it's probably more relevant.
                raise e
            x = np.array([w[0] for w in image])
            y = np.array([w[1] for w in image])
            return x, y

    def _newOrigin(self, origin, world_origin):
        return CubicWCS(self._scale, self._r3, origin, world_origin)

    def _withOrigin(self, origin, world_origin, color):
        return self._newOrigin(origin, world_origin)

    def _writeHeader(self, header, bounds):
        header["GS_WCS"] = ("CubicWCS", "GalSim WCS name")
        header["GS_X0"] = (self.origin.x, "GalSim image origin x")
        header["GS_Y0"] = (self.origin.y, "GalSim image origin y")
        header["GS_RA0"] = (self.world_origin.ra.deg, "GalSim world origin x")
        header["GS_DEC0"] = (self.world_origin.dec.deg,
                             "GalSim world origin x")
        header["GS_SCALE"] = (self._scale, "Nominal pixel scale")
        header["GS_R3"] = (self._r3, "Cubic coefficient")

        return self.affine(bounds.true_center)._writeLinearWCS(header, bounds)

    @staticmethod
    def _readHeader(header):
        x0 = header["GS_X0"]
        y0 = header["GS_Y0"]
        ra0 = header["GS_RA0"]
        dec0 = header["GS_DEC0"]
        scale = header['GS_SCALE']
        r3 = header['GS_R3']

        return CubicWCS(scale, r3, galsim.PositionD(x0, y0),
                        galsim.CelestialCoord(ra0*galsim.degrees, dec0*galsim.degrees))
        return None

    def copy(self):
        ""
        return CubicWCS(self._scale, self._r3, self.origin, self.world_origin)

    def __eq__(self, other):
        ""
        return (self is other or
                (isinstance(other, CubicWCS) and
                 self._scale == other._scale and
                 self._xy_r3 == other._r3 and
                 self.origin == other.origin and
                 self.world_origin == other.world_origin
                 ))

    def __repr__(self):
        return "galsim.CubicWCS(%r, %r, %r %r)" % (self._scale, self._r3, self.origin, self.world_origin)

    def __hash__(self): return hash(repr(self))

    def __getstate__(self):
        d = self.__dict__.copy()
        return d

    def __setstate__(self, d):
        self.__dict__ = d


def to_mag(flux, zp):
    """Helper function to from flux -> magnitude"""
    return -2.5*np.log10(flux) + zp


def to_flux(mag, zp):
    """Helper function to from magnitude -> flux"""
    return 10**(-0.4*(mag-zp))


def get_gaia_stars(image, gaia_dir, depth=7):
    """Get gaia stars that overlap an image

    Given an image with a valid wcs this function will find all overlapping htm indexes and
    extract the stars that overlap the region.  It expects the stars have already been divided
    into separate files for each htm index.  This is the current format for Rubin references
    catalogs.

    """
    try:
        import esutil
    except:
        raise Exception('The package esutil is required for the htm computations')
    htm = esutil.htm.HTM(depth)
    wcs = image.wcs
    center = wcs.posToWorld(image.true_center)
    max_dist = 0*galsim.radians

    # Find the corner that is furthest away from the center
    corners = [
        (image.bounds.xmin, image.bounds.ymin),
        (image.bounds.xmax, image.bounds.ymin),
        (image.bounds.xmax, image.bounds.ymax),
        (image.bounds.xmin, image.bounds.ymax)
    ]

    for x, y in corners:
        sky = wcs.posToWorld(galsim.PositionD(x, y))
        dist = center.distanceTo(sky)
        if dist > max_dist:
            max_dist = dist
    radius = max_dist.deg + 0.1
    shards = htm.intersect(center.ra.deg, center.dec.deg, radius)

    table_list = []
    for shard in shards:
        file = f"{gaia_dir}/{shard}.fits"
        if os.path.exists(file) is False:
            continue
        print('Reading gaia file', file)
        data = Table.read(file)
        mask = np.zeros(len(data), dtype=bool)
        for i, d in enumerate(data):
            world = galsim.CelestialCoord(
                d['coord_ra']*galsim.radians, d['coord_dec']*galsim.radians)
            pos = wcs.posToImage(world)
            if image.bounds.includes(pos):
                mask[i] = True
        table_list.append(data[mask])

    table = vstack(table_list)

    return table


def simulate_image(sat_mag_min, sat_mag_max, snr, psf_type, position_box,
                   site, rng, center, use_local, gaia_dir=None, 
                   n_stars=0, star_len=50, star_mag_min=16, star_mag_max=22):
    """
    Simulate an image with one satellite and mutliple stars in "target" mode where
    the satellite is a PSF.

    The satellite will be simulated with the given magnitude and the noise will
    be adjusted such that the object has the specified S/N ratio.  The object
    will be randolmly positioned around the central area of the image determined
    by the position_box parameter.

    There are three options for a wcs:
      - CubicWCS - if site['r3'] > 0 and center is not None.  This is a
                   radial distortion + a gnomic projection.
      - TanWCS - if site['r3'] <= 0 and center is not None.  This does not
                 include off-diagonal terms in the cd matrix.
      - PixelScale - if site['r3'] <=0 and center=None.  This is a constant
                     pixel scale over the whole image.

    If the use_local option is selected and the stars are streaks, the star images
    will use a locally linear approximation to the wcs.  This can increase the
    speed of the simulations.


    For star objects if the gaia_dir is specified and the wcs is either a
    CubicWCS or TanWCS, then objects in the gaia catalogs will be used and the
    star_mag parameters will be ignored.  Otherwise stars at random positions will
    be generated.  Streaks for Gaia stars will be simulated with the center of
    the streak corresponding to the Gaia position


    Parameters:
        sat_mag_min:      Minimum magnitude sampled from for satellite magnitude.
        sat_mag_max:      Maximum magnitude sampled from for satellite magnitude.
        snr:              Flux S/N ratio.  The noise is adjusted to get this S/N for the satellite.
        psf_type:         PSF type, currently only Gaussian is allowed.
        position_box:     Randomly place satellite in this box around the center of the image.
        site:             Parameters of the telescope.  Must currently include:
                          fov, scale, r3, ZP, psf_size
        rng:              galsim.BaseDeviate for random number generation.
        use_local:        Use local wcs when simulating star streaks.
        gaia_dir:         Location of gaia files.
        n_stars:          Number of stars to simulate.  Ignored if gaia_dir specified.
        star_len:         Length of star streaks to simulate.
        star_mag_min:     Minimum magnitude to sample from for star magnitude.
                          Ignored if gaia_dir specified.
        star_mag_min:     Maximum magnitude to sample from for star magnitude.
                          Ignored if gaia_dir specified.

    Returns:
        an image of the scene and a catalog with information on the simulated sources.
    """

    for d in ["public", "private"]:
        if not os.path.exists(d):
            os.makedirs(d)

    sample_star_dicts = []
    sample_sat_dicts = []

    table = Table(
        names=[
            'object_type', 'x', 'y', 'x0', 'y0', 'x1', 'y1', 'flux', 'snr'
        ],
        dtype=["10str", "f8", "f8", "f8", "f8", "f8", "f8", "f8", "f8"]
    )

    x_size_arcmin, y_size_arcmin = site['fov']
    x_size = x_size_arcmin*60 / site['scale']
    y_size = y_size_arcmin*60 / site['scale']

    image = galsim.ImageF(x_size, y_size)

    if site['r3'] > 0 and center:
        wcs = CubicWCS(site['scale'], site['r3'], image.true_center, center)
    elif center:
        wcs = galsim.TanWCS(galsim.AffineTransform(site['scale'], 0, 0, site['scale'],
                                                   origin=image.true_center),
                            world_origin=center)
    else:
        wcs = galsim.PixelScale(site['scale'])

    image.wcs = wcs

    psf = None
    psf_size_min, psf_size_max = site['psf_size']
    psf_size = np.random.uniform(psf_size_min, psf_size_max)
    if psf_type == 'Gaussian':
        psf = galsim.Gaussian(fwhm=psf_size)
    else:
        raise Exception(f'Not a valid psf_type: {psf_type}')

    angle = np.random.uniform(0, 2*np.pi)
    # Draw primary object with some offset from the center
    cx = image.true_center.x
    cy = image.true_center.y

    cbox = galsim.BoundsI(int(cx), int(cx), int(cy), int(cy))
    cbox = cbox.withBorder(int((position_box*60/site['scale'])/2))
    x = np.random.uniform(cbox.xmin, cbox.xmax)
    y = np.random.uniform(cbox.ymin, cbox.ymax)


    sat_mag = np.random.uniform(sat_mag_min, sat_mag_max)
    sat_flux = to_flux(sat_mag, site['ZP'])

    sat = psf.withFlux(sat_flux)
    #print(f"Adding satellite at: {x}, {y} with mag: {mag}, snr: {snr}, PSF size: {psf_size}")

    # draw image with Gaussian noise to determine the variance needed to acheive the desired S/N
    test_noise = galsim.GaussianNoise(rng, sigma=1)
    local_wcs = local_wcs = wcs.local(galsim.PositionD(cx, cy))
    test_image = sat.drawImage(wcs=local_wcs)
    var = test_image.addNoiseSNR(test_noise, snr, preserve_flux=True)

    # account for noise in variable pixel size
    noise = None
    if False:
        var_image = galsim.ImageF(xsize, ysize)
        var_image = wcs.makeSkyImage(var_image, 1)
        var_image /= np.min(var.image)

        var_image *= var
        noise = galsim.VariableGaussianNoise(rng, var_image)
    else:
        noise = test_noise.withVariance(var)

    image.addNoise(noise)
    sat.drawImage(image, center=(x, y), add_to_image=True)

    table.add_row(('sat', x, y, x, y, x, y, sat_flux, snr))

    sample_sat_dicts.append({
        'x0' : float(x + rng_dicts.normal(scale=1.0)),
        'y0' : float(y + rng_dicts.normal(scale=1.0)),
        'x1' : float(x + rng_dicts.normal(scale=1.0)),
        'y1' : float(y + rng_dicts.normal(scale=1.0)),
        'flux' : float(sat_flux + rng_dicts.normal(scale=sat_flux/20.)),
    })

    if gaia_dir is None or center is None:
        star_mag = np.random.uniform(star_mag_min, star_mag_max, size=n_stars)
        star_fluxes = to_flux(star_mag, site['ZP'])
        sx = np.random.uniform(image.bounds.xmin, image.bounds.xmax, size=n_stars)
        sy = np.random.uniform(image.bounds.ymin, image.bounds.ymax, size=n_stars)
        sx0 = sx - star_length/2 * np.cos(angle)
        sx1 = sx + star_length/2 * np.cos(angle)
        sy0 = sy - star_length/2 * np.sin(angle)
        sy1 = sy + star_length/2 * np.sin(angle)
    else:
        cat = get_gaia_stars(image, gaia_dir)

        star_mag = cat['phot_g_mean_flux'].to(u.ABmag).value
        ra = cat['coord_ra']
        dec = cat['coord_dec']
        mask = (star_mag > star_mag_min) & (star_mag < star_mag_max)
        star_mag = star_mag[mask]
        ra = ra[mask]
        dec = dec[mask]
        sx = []
        sy = []
        sx0 = []
        sy0 = []
        sx1 = []
        sy1 = []
        star_fluxes = to_flux(star_mag, site['ZP'])

        for rra, ddec in zip(ra, dec):
            sky = galsim.CelestialCoord(rra*galsim.radians, ddec*galsim.radians)
            pos = wcs.posToImage(sky)
            sx.append(pos.x)
            sy.append(pos.y)
            sx0.append(pos.x - star_length/2 * np.cos(angle))
            sx1.append(pos.x + star_length/2 * np.cos(angle))
            sy0.append(pos.y - star_length/2 * np.sin(angle))
            sy1.append(pos.y + star_length/2 * np.sin(angle))

    for i in range(len(star_mag)):

        star_flux = to_flux(star_mag[i], site['ZP'])

        #print(f"Adding star at: {sx[i]}, {sy[i]} with mag: {star_mag[i]}, length: {star_len}")
        box = galsim.Box(height=0.01, width=star_len*site['scale'],
                         flux=star_flux).rotate(angle*galsim.radians)
        star = galsim.Convolve([box, psf])

        if use_local:
            local_wcs = wcs.local(galsim.PositionD(sx[i], sy[i]))
            stamp = star.drawImage(wcs=local_wcs, center=(sx[i], sy[i]))
            bounds = stamp.bounds & image.bounds
            image[bounds] += stamp[bounds]
        else:
            star.drawImage(image, center=(sx[i], sy[i]), add_to_image=True)

        table.add_row(('star', sx[i], sy[i], sx0[i], sy0[i], sx1[i], sy1[i], star_fluxes[i], np.nan))

        sample_star_dicts.append({
            'x': float(sx[i] + rng_dicts.normal(scale=1.0)),
            'y': float(sy[i] + rng_dicts.normal(scale=1.0)),
            'flux': float(star_fluxes[i] +
                          rng_dicts.normal(scale=star_fluxes[i]/20.))
            })
            #Add to sample_star_dicts for endpoints
            #'x0': float(sx0[i] + rng_dicts.normal(scale=1.0)),
            #'x1': float(sx1[i] + rng_dicts.normal(scale=1.0)),
            #'y0': float(sy0[i] + rng_dicts.normal(scale=1.0)),
            #'y1': float(sy1[i] + rng_dicts.normal(scale=1.0)),

    return image, table, sample_star_dicts, sample_sat_dicts

if __name__ == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument('config', default='simulate_target_config.yaml', help='config file')
    #parser.add_argument("--outdir", default='output/', type=str)
    args = parser.parse_args()

    #if os.path.exists(args.outdir) is False:
    #    os.makedirs(args.outdir)
    params = yaml.safe_load(open(args.config, 'r'))
    site = params['telescopes'][params['sims']['site']]

    np.random.seed(params['sims']['seed'])
    rng = galsim.BaseDeviate(params['sims']['seed'])
    rng_dicts = np.random.default_rng(params['sims']['seed'])
    snr_min, snr_max = params['sims']['snr']
    psf_type = params['sims']['psf_type']
    position_box = params['sims']['position_box']
    use_gaia = params['sims']['use_gaia']
    gaia_dir = None
    if use_gaia:
        gaia_dir = params['sims']['gaia_dir']

    if params['sims']['center']:
        center = galsim.CelestialCoord(params['sims']['center'][0]*galsim.degrees,
                                       params['sims']['center'][1]*galsim.degrees)
    else:
        center = None

    star_mag_min, star_mag_max = params['sims']['star_mag']
    sat_mag_min, sat_mag_max = params['sims']['mag']
    star_length_min, star_length_max = params['sims']['star_length']
    n_star = params['sims']['n_star']
    use_local = params['sims']['use_local']
    ims = []

    truth5_hdul = fits.HDUList()
    truth_hdul = fits.HDUList()
    sample_docs = [{'branch' : params['branch'], 'competitor_name' : params['competitor_name'], 'display_true_name' : params['display_true_name']}]

    for i in tqdm.tqdm(range(params['sims']['n_obs']+5)):
        snr = np.random.uniform(snr_min, snr_max)
        star_length = np.random.uniform(star_length_min, star_length_max)
        im, truth_table, sample_star_dicts, sample_sat_dicts = simulate_image(
            sat_mag_min, sat_mag_max, snr, psf_type, position_box, site,
            rng, center, use_local, gaia_dir, n_star,
            star_length, star_mag_min, star_mag_max)
        ims.append(im)
        im.write(f'public/{i:03d}.fits')

        sats = truth_table[truth_table['object_type']=='sat']
        satsx0 = sats['x0']
        satsx1 = sats['x1']
        satsy0 = sats['y0']
        satsy1 = sats['y1']
        sat_fluxes = sats['flux']

        stars = truth_table[truth_table['object_type']=='star']
        sx = stars['x']
        sy = stars['y']
        sx0 = stars['x0']
        sx1 = stars['x1']
        sy0 = stars['y0']
        sy1 = stars['y1']
        star_fluxes = stars['flux']

        star_cols = fits.ColDefs([
            fits.Column(name='x', format='E', array=sx),
            fits.Column(name='y', format='E', array=sy),
            fits.Column(name='flux', format='E', array=star_fluxes)
        ])
            #Add to star_cols for endpoints
            #fits.Column(name='x0', format='E', array=sx0),
            #fits.Column(name='y0', format='E', array=sy0),
            #fits.Column(name='x1', format='E', array=sx1),
            #fits.Column(name='y1', format='E', array=sy1),

        sat_cols = fits.ColDefs([
            fits.Column(name='x0', format='E', array=satsx0),
            fits.Column(name='y0', format='E', array=satsy0),
            fits.Column(name='x1', format='E', array=satsx1),
            fits.Column(name='y1', format='E', array=satsy1),
            fits.Column(name='flux', format='E', array=sat_fluxes),
        ])

        sat_hdu = fits.BinTableHDU.from_columns(sat_cols, name=f"SAT_{i:04d}")
        star_hdu = fits.BinTableHDU.from_columns(star_cols, name=f"STAR_{i:04d}")

        if i < 5:
            truth5_hdul.append(star_hdu)
            truth5_hdul.append(sat_hdu)
        truth_hdul.append(star_hdu)
        truth_hdul.append(sat_hdu)

        rng_dicts.shuffle(sample_star_dicts)
        sample_docs.append({'file':f"{i:03d}.fits", 'sats':sample_sat_dicts, 'stars':sample_star_dicts})

        if i == 4:
            with open("public/sample_submission_5_target.yaml", "w") as f:
                yaml.safe_dump_all(sample_docs, f)

    with open("private/sample_submission.yaml", "w") as f:
        yaml.safe_dump_all(sample_docs, f)

    truth_hdul.writeto("private/truth.fits", overwrite=True)
    truth5_hdul.writeto("public/truth_5_target.fits", overwrite=True)

