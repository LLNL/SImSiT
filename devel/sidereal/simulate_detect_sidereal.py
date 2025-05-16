import numpy as np
import astropy.io.fits as fits
import galsim
import yaml
import matplotlib.pyplot as plt
import os
import tqdm


def simulate_detect_sidereal(args):
    for d in ["public", "private"]:
        if not os.path.exists(d):
            os.makedirs(d)

    config = yaml.safe_load(open(args.config, 'r'))
    nx, ny = config['fov']
    area = nx * ny

    rng = np.random.default_rng(config['seed'])

    sample_docs = [{'branch' : config['branch'], 'competitor_name' : config['competitor_name'], 'display_true_name' : config['display_true_name']}]
    truth5_hdul = fits.HDUList()
    truth_hdul = fits.HDUList()

    for i_obs in tqdm.tqdm(range(config['n_obs']+config['n_demo'])):
        star_density = np.exp(rng.uniform(*np.log(config['star_density'])))
        n_star = int(star_density * area)
        psf_fwhm = rng.uniform(*config['psf_fwhm'])

        image = galsim.ImageD(nx, ny, scale=1)
        image.array[:] = rng.normal(scale=config['noise'], size=(ny, nx))

        sample_star_dicts = []
        star_xs = rng.uniform(0, nx, n_star)
        star_ys = rng.uniform(0, ny, n_star)
        star_fluxes = np.exp(rng.uniform(*np.log(config['star_flux']), n_star))
        star_fluxes = np.sort(star_fluxes)[::-1]  # brightest first
        for i_star, (x, y, flux) in enumerate(zip(star_xs, star_ys, star_fluxes)):
            star = galsim.Gaussian(fwhm=psf_fwhm)*flux
            star.drawImage(image=image, method='phot', add_to_image=True, center=(x, y))
            if i_star < 50:
                sample_star_dicts.append({
                    'x':float(x+rng.normal(scale=1)),
                    'y':float(y+rng.normal(scale=1)),
                    'flux':float(flux+rng.normal(scale=flux/20)),
                })

        # Add streak
        streak_x = image.center.x + rng.uniform(-0.5, 0.5)*config['position_box']
        streak_y = image.center.y + rng.uniform(-0.5, 0.5)*config['position_box']
        streak_angle = rng.uniform(0, 2*np.pi)
        streak_length = rng.uniform(*config['streak_length'])
        streak_x0 = streak_x - 0.5*streak_length*np.cos(streak_angle)
        streak_x1 = streak_x + 0.5*streak_length*np.cos(streak_angle)
        streak_y0 = streak_y - 0.5*streak_length*np.sin(streak_angle)
        streak_y1 = streak_y + 0.5*streak_length*np.sin(streak_angle)
        streak_flux = np.exp(rng.uniform(*np.log(config['streak_flux'])))
        streak = galsim.Convolve(galsim.Box(streak_length, 1e-8), galsim.Gaussian(fwhm=psf_fwhm))
        streak *= streak_flux
        streak = streak.rotate(streak_angle*galsim.radians)
        streak.drawImage(image=image, method='phot', add_to_image=True, center=(streak_x, streak_y))

        # Write the image publicly
        hdu = fits.PrimaryHDU(image.array)
        hdu.header['XANGLE'] = np.rad2deg(streak_angle) + rng.normal(scale=10.0)
        hdu.header['XLENGTH'] = streak_length + rng.normal(scale=streak_length/10)
        hdul = fits.HDUList([hdu])
        hdul.writeto(f"public/{i_obs:03d}.fits", overwrite=True)

        star_cols = fits.ColDefs([
            fits.Column(name='x', format='E', array=star_xs),
            fits.Column(name='y', format='E', array=star_ys),
            fits.Column(name='flux', format='E', array=star_fluxes)
        ])
        star_hdu = fits.BinTableHDU.from_columns(star_cols, name=f"STAR_{i_obs:04d}")
        # First 5 observations have public truth tables
        if i_obs < 5:
            truth5_hdul.append(star_hdu)
        truth_hdul.append(star_hdu)

        sat_cols = fits.ColDefs([
            fits.Column(name='x0', format='E', array=np.array([streak_x0])),
            fits.Column(name='y0', format='E', array=np.array([streak_y0])),
            fits.Column(name='x1', format='E', array=np.array([streak_x1])),
            fits.Column(name='y1', format='E', array=np.array([streak_y1])),
            fits.Column(name='flux', format='E', array=np.array([streak_flux])),
        ])
        sat_hdu = fits.BinTableHDU.from_columns(sat_cols, name=f"SAT_{i_obs:04d}")
        if i_obs < 5:
            truth5_hdul.append(sat_hdu)
        truth_hdul.append(sat_hdu)

        sample_sat_dicts = []
        sample_sat_dicts.append({
            'x0': float(streak_x0 + rng.normal(scale=1.0)),
            'y0': float(streak_y0 + rng.normal(scale=1.0)),
            'x1': float(streak_x1 + rng.normal(scale=1.0)),
            'y1': float(streak_y1 + rng.normal(scale=1.0)),
            'flux': float(streak_flux + rng.normal(scale=streak_flux/20.))
        })
        # randomize the order of the sample stars, since we don't know what order
        # the competitors will use.
        rng.shuffle(sample_star_dicts)
        sample_docs.append({'file':f"{i_obs:03d}.fits", 'sats':sample_sat_dicts, 'stars':sample_star_dicts})
        # Write out first 5 sample submissions publicly
        if i_obs == (config['n_demo']-1):
            with open("public/sample_submission_5_detect_sidereal.yaml", "w") as f:
                yaml.safe_dump_all(sample_docs, f)

    # Write out complete catalog sample submission privately (so we can test
    # the scoring)
    with open("private/sample_submission.yaml", "w") as f:
        yaml.safe_dump_all(sample_docs, f)

    truth_hdul.writeto("private/truth.fits", overwrite=True)
    truth5_hdul.writeto("public/truth_5_detect_sidereal.fits", overwrite=True)


if __name__ == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument('config', default='simulate_detect_sidereal_config.yaml', help='config file')
    args = parser.parse_args()
    simulate_detect_sidereal(args)
