import yaml
import numpy as np
from datetime import datetime
import astropy.io.fits as fits
import math
import pandas as pd
from os import path
import sys
from astropy.coordinates import SkyCoord

def check_submission(args):

    # Let yaml raise exception here
    with open(args.submission, 'r') as f:
        docs = yaml.safe_load_all(f)
        docs = [doc for doc in docs]
    
    # Checking metadata
    assert isinstance(docs[0]['branch'], str)
    assert isinstance(docs[0]['competitor_name'], str)
    assert isinstance(docs[0]['display_true_name'], bool)

    # Initialize parameters that control printing
    sat_missing_param_calibrate = False
    star_missing_param_detect = False
    star_missing_param_calibrate = False

    # Initialize parameters that control scoring
    calibrate_sat = False
    detect_star = False
    calibrate_star = False

    # Checking minimum information is supplied, and verifying types.
    for idx, doc in enumerate(docs[1:]):
        # Let python raise error if key not found.
        assert isinstance(doc['file'], str)

        for idy, sat in enumerate(doc['sats']):

            sat_detect_var = ['flux', 'x0', 'x1', 'y0', 'y1']
            if all(k in sat for k in sat_detect_var):
                for var in sat_detect_var:
                    assert isinstance(sat[var], float)

            else:
                # Exit the script if any sat fails to have the minimum information supplied
                print("ERROR: Sat detection parameters not found, or some parameters missing. The minimum requirements to be scored on are not met.")
                sys.exit()

            sat_calibrate_var = ['dec0', 'dec1', 'mag', 'ra0', 'ra1']
            if all(k in sat for k in sat_calibrate_var):
                for var in sat_calibrate_var:
                    assert isinstance(sat[var], float)
                if idx == len(doc)-2 and idy == len(doc['sats'])-1: 
                    if not sat_missing_param_calibrate:
                        calibrate_sat = True
            else:
                if not sat_missing_param_calibrate:
                    sat_missing_param_calibrate = True

        if 'stars' in doc:
            for idz, star in enumerate(doc['stars']):

                star_detect_var = ['flux', 'x', 'y']
                if all(k in star for k in star_detect_var):
                    if idx == len(doc)-2 and idz == len(doc['stars'])-1:
                        if not star_missing_param_detect:
                            detect_star = True
                    for var in star_detect_var:
                        assert isinstance(star[var], float)
                else:
                    if not star_missing_param_detect:
                        star_missing_param_detect = True

                star_calibrate_var = ['dec', 'mag', 'ra']
                if all(k in star for k in star_calibrate_var):
                    if idx == len(doc)-2 and idz == len(doc['stars'])-1:
                        if not star_missing_param_calibrate:
                            calibrate_star = True
                    for var in star_calibrate_var:
                        assert isinstance(star[var], float)
                else:
                    if not star_missing_param_calibrate:
                        star_missing_param_calibrate = True

    if detect_star or calibrate_star:
        print("Optional star parameters found.")
    if calibrate_sat:
        print("Scoring on detection and calibration.")
    else:
        print("Scoring on detection only.")

    return calibrate_sat, detect_star, calibrate_star

def score_submission(args, csat, dstar, cstar):
    submission = [s for s in yaml.safe_load_all(open(args.submission, 'r'))]
    truth_hdul = fits.open(args.truth)

    # Start with RMSE of satellite endpoints
    pos_sq_error_thresh = 400.0  # 400 pix^2, i.e., 20 pix distance
    n_true = 0
    n_sub = 0
    n_true_found_detect = 0
    n_true_found_calibrate = 0
    n_sub_assoc_detect = 0
    n_sub_assoc_calibrate = 0
    xy_se = 0
    radec_se = 0
    flux_se = 0
    mag_se = 0

    competitor_name = submission[0]['competitor_name']
    competitor_name_split = competitor_name.split()
    if len(competitor_name_split) > 1:
        competitor_name = "_".join(competitor_name_split)

    if int((len(truth_hdul)-1)/2) == 10:
        print("Scoring 10 truth submissions.")
        subs_to_score = [s for s in submission[1:] if s['file'] in ['000.fits', '001.fits', '002.fits', '003.fits', '004.fits', '005.fits', '006.fits', '007.fits', '008.fits', '009.fits']]
    else: 
        print(f"Scoring {len(submission)-11} submissions.")
        subs_to_score = [s for s in submission[1:] if s['file'] not in ['000.fits', '001.fits', '002.fits', '003.fits', '004.fits', '005.fits', '006.fits', '007.fits', '008.fits', '009.fits']]

    # Loop over images skipping metadata and first 10 submissions
    for sidx, sdoc in enumerate(subs_to_score):   
        idx = int(sdoc['file'].replace(".fits", ""))
        truth_sats = truth_hdul[f"SAT_{idx:03d}"].data
        n_sub += len(sdoc['sats'])
        n_true += len(truth_sats)
        sub_associated_xy = [False]*n_sub
        sub_associated_radec = [False]*n_sub
        # Find best submitted satellite for each true satellite
        # Keep track of whether or not each true sat is found, and also
        # if each submitted sat is associated with a true sat.
        # Note, it *is* possible for a single submitted satellite to get
        # associated with multiple true satellites.

        # Start by loading all submitted sats into numpy arrays
        sdict = {}

        for sat in sdoc['sats']:
            for k in ['dec0', 'dec1', 'flux', 'mag', 'ra0', 'ra1', 'x0', 'x1', 'y0', 'y1']:
                if k in sat.keys():
                    sdict[k] = []
                    sdict[k].append(sat[k])
                    sdict[k] = np.array(sdict[k])

        for x0, y0, x1, y1, ra0, dec0, ra1, dec1, mag, nphot in truth_sats:
            if csat:

                c0_obs = SkyCoord(sdict['ra0'], sdict['dec0'], unit="deg")
                c1_obs = SkyCoord(sdict['ra1'], sdict['dec1'], unit="deg")
                c0_true = SkyCoord(ra0, dec0, unit="rad")
                c1_true = SkyCoord(ra1, dec1, unit="rad")
                sep0 = c0_obs.separation(c0_true)
                sep1 = c1_obs.separation(c1_true)
                sq_error= sep0.degree**2 + sep1.deg**2
                idx= np.argmin(sq_error)
                min_sq_error = sq_error[idx]
                # We need to check for 0 <-> 1 label ambiguity too.
                sep0T = c0_obs.separation(c1_true)
                sep1T = c1_obs.separation(c0_true)
                sq_errorT= sep0T.degree**2 + sep1T.deg**2
                if np.min(sq_errorT) < min_sq_error:
                    idx = np.argmin(sq_errorT)
                    min_sq_error = sq_errorT[idx]
                if min_sq_error < pos_sq_error_thresh:
                    radec_se += min_sq_error
                    mag_se += np.square(sdict['mag'][idx] - mag)
                    n_true_found_calibrate += 1
                    sub_associated_radec[idx] = True
                n_sub_assoc_calibrate += np.sum(sub_associated_radec)
            
            sq_error = np.sum([
                np.square(sdict['x0']-x0),
                np.square(sdict['y0']-y0),
                np.square(sdict['x1']-x1),
                np.square(sdict['y1']-y1)
            ], axis=0)
            idx = np.argmin(sq_error)
            min_sq_error = sq_error[idx]
            # We need to check for 0 <-> 1 label ambiguity too.
            sq_errorT = np.sum([
                np.square(sdict['x0']-x1),
                np.square(sdict['y0']-y1),
                np.square(sdict['x1']-x0),
                np.square(sdict['y1']-y0)
            ], axis=0)
            if np.min(sq_errorT) < min_sq_error:
                idx = np.argmin(sq_errorT)
                min_sq_error = sq_errorT[idx]
            if min_sq_error < pos_sq_error_thresh:
                xy_se += min_sq_error
                flux_se += np.square((sdict['flux'][idx] - nphot)/nphot*(2.5/np.log(10)))
                n_true_found_detect += 1
                sub_associated_xy[idx] = True
            n_sub_assoc_detect += np.sum(sub_associated_xy)

    print()
    print("Scores:")
    print(f"calibrate completeness: {n_true_found_calibrate/n_true}")
    print(f"detect completeness: {n_true_found_detect/n_true}")
    print(f"calibrate fpr: {1-(n_sub_assoc_calibrate/n_sub)}")
    print(f"detect fpr: {1-(n_sub_assoc_detect/n_sub)}")
    print(f"radec_rmse: {np.sqrt(radec_se/n_true_found_calibrate)}")
    print(f"xy_rmse: {np.sqrt(xy_se/n_true_found_detect)}")
    print(f"mag_rmse: {np.sqrt(mag_se/n_true_found_calibrate)}")
    print(f"flux_rmse: {np.sqrt(flux_se/n_true_found_detect)}")

    scores = pd.DataFrame({'comp_name': competitor_name, 'sub_filename': args.submission, 
                            'date': datetime.utcnow().isoformat(), 'detect_completeness': n_true_found_detect/n_true, 
                            'calibrate_completeness': n_true_found_calibrate/n_true,
                            'detect_fpr': 1-(n_sub_assoc_detect/n_sub), 'calibrate_fpr': 1-(n_sub_assoc_calibrate/n_sub),
                            'xy_RMSE': np.sqrt(xy_se/n_true_found_detect), 'radec_RMSE': np.sqrt(radec_se/n_true_found_calibrate),
                            'flux_RMSE': np.sqrt(flux_se/n_true_found_detect),
                            'mag_RMSE': np.sqrt(mag_se/n_true_found_calibrate)}, index=[0])

    # Save scores to file for the dashboard
    if path.exists('score_history_calibrate_sidereal.txt'):
        scores.to_csv('score_history_calibrate_sidereal.txt', mode='a', index=False, header=False)
    else:
        scores.to_csv("score_history_calibrate_sidereal.txt", index=False)


if __name__ == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("submission")
    parser.add_argument("truth")
    args = parser.parse_args()

    csat, dstar, cstar = check_submission(args)
    score_submission(args, csat, dstar, cstar)
