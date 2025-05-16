import os
import glob
import numpy as np
from astropy.io import fits

# what should we do here?
# subtract a sky background
# do a convolution with an X arcsecond Gaussian
# find peaks
# throw out peaks that aren't peaky enough
# compute centroids & fluxes


def dummy_detect(filename, sigma_val=5):
    im = fits.getdata(filename)

    # Calculate the maximum pixel value and standard deviation of image.
    idx = np.unravel_index(np.argmax(im), im.shape)
    max_value = im[idx[0]][idx[1]]

    std_dev = np.std(im)
    # Determine if the image was a detection, given the set threshold.
    if max_value/std_dev >= sigma_val:
        return idx[1]+1, idx[0]+1, max_value
    else:
        return None, None, None


def main():
    from argparse import ArgumentParser
    from tqdm import tqdm
    import yaml

    parser = ArgumentParser()
    parser.add_argument("--public_directory", type=str,
                        default="public/")
    parser.add_argument("--outfile", type=str,
                        default="dummy_detect_target.yaml")
    args = parser.parse_args()

    files = sorted(glob.glob(os.path.join(args.public_directory,
                                          "???.fits")))
    docs = [{'branch' : 'Detect Target', 'competitor_name' : 'LLNL_Dev', 'display_true_name' : True}]

    for filename in tqdm(files):
        x, y, flux = dummy_detect(filename)
        if x is not None:
            out = [dict(x0=float(x), y0=float(y),
                        x1=float(x), y1=float(y),
                        flux=float(flux))]
        else:
            out = []
        docs.append({'file': os.path.basename(filename), 'sats':out})

    with open(args.outfile, "w") as f:
        yaml.safe_dump_all(docs, f)


if __name__ == '__main__':
    main()
