import numpy as np
from skimage import exposure
import astropy.io.fits as fits
import galsim


# Here throught FRT() taken from https://github.com/guynir42/pyradon
def empty(array):
    if array is not None and np.asarray(array).size>0:
        return False
    else:
        return True

# -padding (True): adds zero padding to the active axis to fill up powers of 2.
def padMatrix(M):
    N = M.shape[0]
    dN = int(2**np.ceil(np.log2(N))-N)
#    print "must add "+str(dN)+" lines..."
    M = np.vstack((M, np.zeros((dN, M.shape[1]), dtype=M.dtype)))
    return M

# -expand (False): adds zero padding to the sides of the passive axis to allow
#                  for corner-crossing streaks
def expandMatrix(M):
    Z = np.zeros((M.shape[0],M.shape[0]), dtype=M.dtype)
    M = np.hstack((Z,M,Z))
    return M

# optimized to do array math in place.
def shift_add_in_place(M1, M2, gap, output):
    output[:] = M1
    if gap > 0:
        output[gap:] += M2[:-gap]
    elif gap < 0:
        output[:gap] += M2[-gap:]
    else:
        output += M2

def getPartialDims(M, log_level):
    x = M.shape[1]
    y = M.shape[0]

    y = int(y/2**log_level)
    z = int(2**(log_level+1)-1)

    return (z,y,x)

def getNumLogFoldings(M):
    return int(np.ceil(np.log2(M.shape[0])))

def FRT(
    M_in,
    transpose=False,
    expand=False,
    padding=True,
    partial=False,
    output=None
):
    """ Fast Radon Transform (FRT) of the input matrix M_in (must be 2D numpy
        array)
    Additional arguments:
     -transpose (False): transpose M_in (replace x with y) to check all the
                         other angles.
     -expand (False): adds zero padding to the sides of the passive axis to
                      allow for corner-crossing streaks
     -padding (True): adds zero padding to the active axis to fill up powers of
                      2.
     -partial (False): use this to save second output, a list of Radon partial
                       images (useful for calculating variance at different
                       length scales)
     -output (None): give the a pointer to an array with the right size, for FRT
                     to put the return value into it.  Note that if partial=True
                     then output must be a list of arrays with the right
                     dimensions.
    """
    ############### CHECK INPUTS AND DEFAULTS #################################

    if empty(M_in):
        return

    if M_in.ndim>2:
        raise Exception("FRT cannot handle more dimensions than 2D")

    ############## PREPARE THE MATRIX #########################################

    M = np.array(M_in) # keep a copy of M_in to give to finalizeFRT

    np.nan_to_num(M, copy=False)  # get rid of NaNs (replace with zeros)

    if transpose:
        M = M.T

    if padding:
        M = padMatrix(M)

    if expand:
        M = expandMatrix(M)

    ############## PREPARE THE MATRIX #########################################

    Nfolds = getNumLogFoldings(M)
    (Nrows, Ncols) = M.shape

    M_out = []

    if not empty(output): # will return the entire partial transform list
        if partial:
            for m in range(2,Nfolds+1):
                if output[m-1].shape!=getPartialDims(M,m):
                    raise RuntimeError(
                        "Wrong dimensions of output array["+str(m-1)+"]: "
                        +str(output[m-1].shape)
                        +", should be "
                        +str(getPartialDims(M,m))
                    )

            M_out = output

        else:
            if output.shape[0]!=2*M.shape[0]-1:
                raise RuntimeError(
                    "Y dimension of output ("
                    +str(output.shape[0])
                    +") is inconsistent with (padded and doubled) input ("
                    +str(M.shape[0]*2-1)
                    +")"
                )
            if output.shape[1]!=M.shape[1]:
                raise RuntimeError(
                    "X dimension of output ("
                    +str(output.shape[1])
                    +") is inconsistent with (expanded?) input ("
                    +str(M.shape[1])
                    +")"
                )

    dx = np.array([0], dtype='int64')

    M = M[np.newaxis,:,:]

    for m in range(1, Nfolds+1): # loop over logarithmic steps

        M_prev = M
        dx_offset = int(len(dx)/2)

        Nrows = M_prev.shape[1]

        max_dx = 2**(m)-1
        dx = np.arange(-max_dx, max_dx+1)
        if partial and not empty(output):
            # we already have memory allocated for this result
            M = M_out[m-1]
        else:
            # make a new array each time
            M = np.zeros((len(dx), Nrows//2, Ncols), dtype=M.dtype)

        counter = 0

        # shifts are independent of i, j, so precompute before loop
        prev_dxs = (dx.astype(float)/2).astype(int)
        prev_js = prev_dxs + dx_offset
        gap_xs = dx - prev_dxs

        # loop over pairs of rows (number of rows in new M)
        for i in range(Nrows//2):
            for j, (dxj, j_in_prev, gap_x) in enumerate(
                zip(dx, prev_js, gap_xs)
            ):
                M1 = M_prev[j_in_prev, counter, :]
                M2 = M_prev[j_in_prev, counter+1,:]
                shift_add_in_place(M1, M2, -gap_x, M[j,i,:])
            counter+=2

        # only append to the list if it hasn't been given from the start using
        # "output"
        if partial and empty(output):
            M_out.append(M)

#     end of loop on m

    if not partial:
        # we don't care about partial transforms, we were not given an array to
        # fill
        M_out = M[:,0,:] # lose the empty dim

        if not empty(output): # do us a favor and also copy it into the array
            np.copyto(output, M_out)
            # this can be made more efficient if we use the "output" array as
            # target for assignment at the last iteration on m.
            # this will save an allocation and a copy of the array.
            # however, this is probably not very expensive and not worth
            # the added complexity of the code

    return M_out

def refine(image, x0, y0, x1, y1, show=False):
    """Refine a detected streak by perturbing one endpoints along the streak
    direction.
    """
    dx = x1 - x0
    dy = y1 - y0
    if np.abs(dy) < np.abs(dx):
        y1s, x1s = refine(image.T, y0, x0, y1, x1, show=show)
        return x1s, y1s
    ys = np.arange(y0, y0+2*dy, np.sign(dy))
    xs = ((ys-y0)*dx/dy+x0).astype(int)

    w = (ys >= 0)
    w &= (ys < image.shape[0])
    w &= (xs >= 0)
    w &= (xs < image.shape[1])

    signal = np.zeros(len(xs), dtype=float)
    signal[w] = image[ys[w], xs[w]]
    signal = np.cumsum(signal)/np.sqrt(np.arange(2*np.abs(dy))+1)

    # ignore first quarter..., as we insist streak doesn't shrink by more than
    # 0.5x
    start = np.abs(dy)//2
    best_dy = np.sign(dy)*(np.argmax(signal[start:])+start)
    if show:
        import matplotlib.pyplot as plt
        plt.plot(signal)
    return int(x0+(best_dy*dx/dy)), y0+best_dy

def findStreak(image, histEq=True, show=False):
    """Find a single streak using a partial Radon transform.
    """
    original_image = image
    if histEq:
        image = image - np.min(image)
        image /= np.max(image)
        image = exposure.equalize_hist(image)
        image -= np.mean(image)
        var = np.var(image.ravel())
    else:
        var = np.var(image[image<np.percentile(image, 95)])

    # First no transpose
    prt = FRT(image, partial=True, expand=False)
    snrs = np.array([np.nanmax(rt) for rt in prt])
    snrs /= np.sqrt(var*2**np.arange(2, len(prt)+2))

    best_snr_idx = np.nanargmax(snrs)
    best_snr = np.nanmax(snrs)
    subrow, subimage, col = np.unravel_index(
        np.nanargmax(prt[best_snr_idx]),
        prt[best_snr_idx].shape
    )
    subrows, subimages, cols = prt[best_snr_idx].shape
    x0 = col
    dx = subrow - subrows//2
    x1 = col + dx
    y0 = subimage*(2**(best_snr_idx+1))
    y1 = (subimage + 1)*2**(best_snr_idx+1) - 1

    # Now check transpose
    prt = FRT(image.T, partial=True, expand=False)
    snrs = np.array([np.nanmax(rt) for rt in prt])
    snrs /= np.sqrt(var*2**np.arange(2, len(prt)+2))
    if np.nanmax(snrs) > best_snr:
        best_snr_idx = np.nanargmax(snrs)
        best_snr = np.nanmax(snrs)
        subrow, subimage, col = np.unravel_index(
            np.nanargmax(prt[best_snr_idx]),
            prt[best_snr_idx].shape
        )
        subrows, subimages, cols = prt[best_snr_idx].shape
        x0 = col
        dx = subrow - subrows//2
        x1 = col + dx
        y0 = subimage*(2**(best_snr_idx+1))
        y1 = (subimage + 1)*2**(best_snr_idx+1) - 1
        x0, y0 = y0, x0
        x1, y1 = y1, x1

    # Refine estimate
    x0_original, y0_original = x0, y0
    x1_original, y1_original = x1, y1
    x1, y1 = refine(image, x0, y0, x1, y1)
    x0, y0 = refine(image, x1, y1, x0, y0)

    if show:
        vmin, vmax = np.percentile(original_image, [1, 99])
        plt.figure(figsize=(22, 22))
        plt.imshow(original_image, vmin=vmin, vmax=vmax)

        plt.scatter(
            [x0_original, x1_original],
            [y0_original, y1_original],
            c='r'
        )
        plt.scatter(
            [x0, x1],
            [y0, y1],
            facecolors='none',
            edgecolors='b',
            s=100
        )

        plt.show()
    return x0, y0, x1, y1


def estimateFlux(image, x0, y0, x1, y1):
    dx = x1 - x0
    dy = y1 - y0
    if np.abs(dy) < np.abs(dx):
        flux = estimateFlux(image.T, y0, x0, y1, x1)
        return flux
    xss = []
    yss = []
    ys = np.arange(y0, y0+2*dy, np.sign(dy))
    xs = ((ys-y0)*dx/dy+x0).astype(int)
    for ddx in range(-6, 7):
        yss.append(ys)
        xss.append(xs+ddx)
    xss = np.array(xss).ravel()
    yss = np.array(yss).ravel()

    w = (yss >= 0)
    w &= (yss < image.shape[0])
    w &= (xss >= 0)
    w &= (xss < image.shape[1])

    return np.sum(image[yss[w], xss[w]])


def main():
    from argparse import ArgumentParser
    from tqdm import tqdm
    import glob
    import os
    import yaml

    parser = ArgumentParser()
    parser.add_argument("--public_directory", type=str, default="public/")
    parser.add_argument("--outfile", type=str, default="radon_detect.yaml")
    args = parser.parse_args()

    files = glob.glob(os.path.join(args.public_directory, "???.fits"))

    docs = [{'branch' : 'Detect Sidereal', 'competitor_name' : 'Competitor A', 'display_true_name' : True}]

    for file in tqdm(files):
        image = fits.getdata(file)
        x0, y0, x1, y1 = findStreak(image)
        flux = estimateFlux(image, x0, y0, x1, y1)
        out = {
            'x0':float(x0),
            'y0':float(y0),
            'x1':float(x1),
            'y1':float(y1),
            'flux':float(flux)
        }
        docs.append({'file':os.path.basename(file), 'sats':[out]})
    with open(args.outfile, "w") as f:

        yaml.safe_dump_all(docs, f)


if __name__ == '__main__':
    main()
