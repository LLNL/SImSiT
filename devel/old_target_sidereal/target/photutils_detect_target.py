import os
import glob
import numpy as np
import photutils
import astropy
from astropy.io import fits


def fit_im(im, isig=None, dqim=None):
    if dqim is None:
        dqim = np.zeros_like(im, dtype='i4')
    if isig is None:
        pass
    im = im.copy()
    mask = photutils.make_source_mask(im, nsigma=2, npixels=5, 
                                      dilate_size=11)
    mean, median, std = astropy.stats.sigma_clipped_stats(
        im, sigma=3.0, mask=mask)
    if isig is None:
        isig = np.zeros_like(im, dtype='f4') + 1/std
    daofind = photutils.DAOStarFinder(fwhm=6.0, threshold=5*std)
    sources = daofind(im - median)
    if sources is None:
        return []
    sources['xcentroid'] += 1  # convert to galsim standard?
    sources['ycentroid'] += 1  # convert to galsim standard?
    sources = sources[np.abs(sources['roundness2']) < 0.4]
    sources = sources[np.abs(sources['roundness1']) < 0.5]
    bd = 3
    onedge = ((sources['xcentroid'] < bd+1) |
              (sources['ycentroid'] < bd+1) |
              (sources['xcentroid'] > im.shape[1]-bd) |
              (sources['ycentroid'] > im.shape[0]-bd))
    sources = sources[~onedge]
    return sources


def main():
    from argparse import ArgumentParser
    from tqdm import tqdm
    import yaml

    parser = ArgumentParser()
    parser.add_argument("--public_directory", type=str,
                        default="public/")
    parser.add_argument("--outfile", type=str,
                        default="photutils_detect_target.yaml")
    args = parser.parse_args()

    files = sorted(glob.glob(os.path.join(args.public_directory,
                                          "???.fits")))
    docs = [{'branch': target, 'competitor_name' : 'LLNL_Dev', 'display_true_name' : True}]

    for filename in tqdm(files):
        im = fits.getdata(filename)
        sources = fit_im(im)
        if len(sources) > 0:
            out = [dict(x0=float(sources['xcentroid'][i]), 
                        y0=float(sources['ycentroid'][i]),
                        x1=float(sources['xcentroid'][i]), 
                        y1=float(sources['ycentroid'][i]),
                        flux=float(sources['flux'][i])) 
                   for i in range(len(sources))]
        else:
            out = []
        docs.append({'file': os.path.basename(filename), 'sats':out})

    with open(args.outfile, "w") as f:
        yaml.safe_dump_all(docs, f)


if __name__ == '__main__':
    main()

