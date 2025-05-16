import os
import yaml
import numpy as np
import glob
from datetime import datetime
from astropy.io import fits
from astropy.wcs import WCS, utils
from astropy import wcs
from astropy.time import Time
import matplotlib.pyplot as plt
from astropy.visualization import (ImageNormalize, ZScaleInterval)
from matplotlib.patches import Circle
from astropy.coordinates import SkyCoord
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.patches as patches
import warnings
warnings.filterwarnings('ignore', category=wcs.FITSFixedWarning)

def check_endpoints(path, branch, file_format, num_img, total_obs, verbose=False):
    all_files = glob.glob(os.path.join(path, 'public', file_format))
    files = np.random.choice(glob.glob(os.path.join(path, 'public', file_format)), num_img)
    pdf = PdfPages('tests/check_endpoints_'+str(branch)+'.pdf')

    if (branch == 'sidereal_track') or (branch == 'target_track'):
        all_nums = []
        for number in range(1, total_obs+1):
            all_nums.append(pad_number(number))

        last_obs = []
        for num in all_nums:
            all_obs = sorted([f for f in all_files if f.split("/")[4].startswith(num)])
            last_obs.append(sorted(all_obs)[-1])
    
    for idx, file in enumerate(files):
        print("*========== "+str(file)+" ==========*")
        if (branch == 'sidereal_track') or (branch == 'target_track'):
            num = file.split("/")[-1][:8]
        else:
            num = file.split("/")[-1][:4]
        img_hdul = fits.open(file)
        img = img_hdul[0].data
        
        if verbose:
            print(file)
            print(img_hdul[0].header)
 
        truth_wcs_hdul = fits.open(os.path.join(path, 'private', str(num)+'.wcs.fits'))
        truth_wcst_hdul = fits.open(os.path.join(path, 'private', str(num)+'.wcst.fits'))
        hdr = truth_wcs_hdul[0].header
        hdrt = truth_wcst_hdul[0].header
        w = WCS(hdr)
        wt = WCS(hdrt)

        fig, ax = plt.subplots(1)
        ax.imshow(img, cmap='gray', norm=ImageNormalize(img, interval=ZScaleInterval()))

        if (branch == 'sidereal_track') or (branch == 'target_track'):
            truth_all_hdul = fits.open(os.path.join(path, 'private', 'truth_all.fits'))
            if file in last_obs:
                truth_hdul = fits.open(os.path.join(path, 'private', 'truth.fits'))
                # Get ra, dec value from truth file (this is the ra, dec of satellite in the LAST observation)
                ra = truth_hdul[str('SAT_'+str(num[:4]))].data.ra[0]
                dec = truth_hdul[str('SAT_'+str(num[:4]))].data.dec[0]
                ra_dec_pix = w.wcs_world2pix([[ra, dec]], 1)
                print("truth.fits - ra0, dec0:                 "+str(ra_dec_pix[0][0])+", "+str(ra_dec_pix[0][1]))
                ax.add_patch(patches.Circle((ra_dec_pix[0][0], ra_dec_pix[0][1]), 60, color='cyan', fill=False, linewidth=2, label="truth.fits: ra/dec"))
        else:
            truth_all_hdul = fits.open(os.path.join(path, 'private', 'truth.fits'))

        # Get ra, dec from truth file (this is the ra and dec of the endpoints in each image)
        ra0 = truth_all_hdul[str('SAT_'+str(num))].data.ra0[0]
        ra1 = truth_all_hdul[str('SAT_'+str(num))].data.ra1[0]
        dec0 = truth_all_hdul[str('SAT_'+str(num))].data.dec0[0]
        dec1 = truth_all_hdul[str('SAT_'+str(num))].data.dec1[0]
        ra_dec_0_pix = w.wcs_world2pix([[ra0, dec0]], 1)
        ra_dec_1_pix = w.wcs_world2pix([[ra1, dec1]], 1)
        ra_dec_1_pix_t = wt.wcs_world2pix([[ra1, dec1]], 1)

        # get x, y data from truth file (this is the x and y of the endpoints in each image)
        x0 = truth_all_hdul[str('SAT_'+str(num))].data.x0_FITS[0]
        x1 = truth_all_hdul[str('SAT_'+str(num))].data.x1_FITS[0]
        y0 = truth_all_hdul[str('SAT_'+str(num))].data.y0_FITS[0]
        y1 = truth_all_hdul[str('SAT_'+str(num))].data.y1_FITS[0]
            
        #ax.add_patch(patches.Circle((x0,y0), 20, color='blue', fill=False, linewidth=2, label="truth_all.fits: x/y"))
        #ax.add_patch(patches.Circle((x1,y1), 20, color='blue', fill=False, linewidth=2))
        #ax.add_patch(patches.Circle((ra_dec_0_pix[0][0], ra_dec_0_pix[0][1]), 35, color='pink', fill=False, linewidth=2, label="truth_all.fits: ra0/dec0"))
        #ax.add_patch(patches.Circle((ra_dec_1_pix[0][0], ra_dec_1_pix[0][1]), 35, color='orange', fill=False, linewidth=2, label="truth_all.fits: ra1/dec1"))
        #ax.add_patch(patches.Circle((ra_dec_1_pix_t[0][0], ra_dec_1_pix_t[0][1]), 50, color='red', fill=False, linewidth=2, label="truth_all.fits tranformed: ra1/dec1"))
        
        print("truth_all.fits transformed - ra1, dec1: "+str(ra_dec_1_pix_t[0][0])+", "+str(ra_dec_1_pix_t[0][1]))
        print("truth_all.fits - x1, y1:                "+str(x1)+", "+str(y1))
        print("*******************************************************************************")
        print("truth_all.fits - ra0, dec0: "+str(ra_dec_0_pix[0][0])+", "+str(ra_dec_0_pix[0][1]))
        print("truth_all.fits - x0, y0:    "+str(x0)+", "+str(y0))
        
        ax.invert_yaxis()
        #plt.legend()
        plt.title(file)
        plt.show()
        pdf.savefig(fig)
    pdf.close()

def pad_number(number):
    if number < 10:
        num = '000'+str(number)
    elif number > 999:
        num = number
    elif number < 100:
        num = '00'+str(number)
    else:
        num = '0'+str(number)
    return num


def check_times(path, branch, total_obs):
    files = glob.glob(os.path.join(path, 'public', "????_???.fits"))
    truth_hdul = fits.open(os.path.join(path, 'private', 'truth.fits'))
    
    all_nums = []
    for number in range(1, total_obs+1):
        all_nums.append(pad_number(number))

    last_obs = []
    for num in all_nums:
        all_obs = [f for f in files if f.split("/")[4].startswith(num)]
        if all_obs:
            last_obs.append(sorted(all_obs)[-1])
    
    for file in last_obs:
        sat = file.split("/")[4][:8]
        img_hdul = fits.open(file)

        truth_wcs_hdul = fits.open(os.path.join(path, 'private', str(sat)+'.wcs.fits'))
        wcs_start_time = truth_wcs_hdul[0].header["DATE-BEG"]
        wcs_end_time = truth_wcs_hdul[0].header["DATE-END"]

        truth_time = Time(truth_hdul[str('SAT_'+str(sat[:4]))].data.t[0], format='gps').to_value('isot')
        t_start = Time(img_hdul[0].header["DATE-BEG"], format='isot', scale='utc')
        t_end = Time(img_hdul[0].header["DATE-END"], format='isot', scale='utc')

        assert wcs_start_time == truth_time, "WCS exposure start time and time in truth file do not match" 
        assert truth_time == t_start, "Time in truth file and image exposure start time do not match"
        assert wcs_end_time == t_end, "WCS exposure end time and image exposure end time do not match"
    

def check_num_obs(path, branch, total_obs):
    files = glob.glob(os.path.join(path, 'public', "????_???.fits"))
    num_obs = []
    for number in range(total_obs):
        num = pad_number(number)
        amount = 0
        for file in files:
            file_name = file.split("/")[4]
            if file_name.startswith(num):
                amount += 1
        num_obs.append(amount)

    print("Total satellites: "+str(len(num_obs)))
    print("Minimum n_obs = "+str(min(num_obs)))
    print("Average n_obs = "+str(np.average(num_obs)))
    print("Median n_obs = "+str(np.median(num_obs)))
    print("Maximum n_obs = "+str(max(num_obs)))

    plt.figure()
    plt.hist(num_obs)
    plt.title(branch)
    plt.show()
    plt.savefig('tests/num_obs_'+str(branch))


if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument(
        "config",
        type=str
    )
    parser.add_argument(
        "num_obs",
        type=int
    )
    args = parser.parse_args()
    nobs = args.num_obs
    config = yaml.safe_load(open(args.config, 'r'))
    branch = config['meta']['branch']
    path = config['outdir']

    if (branch == 'sidereal_track') or (branch == 'target_track'):
        check_endpoints(path, branch, "????_???.fits", 3, nobs)
        check_num_obs(path, branch, nobs)
        check_times(path, branch, nobs)
    elif (branch == 'sidereal') or (branch == 'target'):
        check_endpoints(path, branch, "????.fits", 3, nobs)
    else:
        print("We do not support that branch yet")
