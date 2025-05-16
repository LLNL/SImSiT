import yaml
import numpy as np
from datetime import datetime
import astropy.io.fits as fits
import pandas as pd
from os import path

def check_submission(args):
    # Let yaml raise exception here
    print(f"Opening file: {args.submission}")
    with open(args.submission, 'r') as f:
        docs = yaml.safe_load_all(f)
        docs = [doc for doc in docs]

    # Checking metadata
    assert isinstance(docs[0]['branch'], str)
    assert isinstance(docs[0]['competitor_name'], str)
    assert isinstance(docs[0]['display_true_name'], bool)
    # Do the format check first, since it's fast. Starting after metadata.
    for doc in docs[1:]:
        # Let python raise error if key not found.
        assert isinstance(doc['file'], str)
        for sat in doc['sats']:
            for var in ['x0', 'y0', 'x1', 'y1', 'flux']:
                assert isinstance(sat[var], float)
        if 'stars' in doc:
            for star in doc['stars']:
                for var in ['x', 'y', 'flux']:
                    assert isinstance(star[var], float)

def score_submission(args):
    submission = [s for s in yaml.safe_load_all(open(args.submission, 'r'))]
    truth_hdul = fits.open(args.truth)

    # Skip branch, competitor_name, and display_name_true, as well as first 5 submissions
    print(f"scoring {len(submission)-6} submissions")

    # Start with RMSE of satellite endpoints
    pos_sq_error_thresh = 400.0  # 400 pix^2, i.e., 20 pix distance
    n_true = 0
    n_sub = 0
    n_true_found = 0
    n_sub_assoc = 0
    pos_se = 0
    mag_se = 0

    competitor_name = submission[0]['competitor_name']
    competitor_name_split = competitor_name.split()
    if len(competitor_name_split) > 1:
        competitor_name = "_".join(competitor_name_split)

    if int((len(truth_hdul)-1)/2) == 5:
        print("Scoring 5 truth submissions.")
        subs_to_score = [s for s in submission[1:] if s['file'] in ['000.fits', '001.fits', '002.fits', '003.fits', '004.fits']]
    else: 
        print(f"Scoring {len(submission)-6} submissions.")
        subs_to_score = [s for s in submission[1:] if s['file'] not in ['000.fits', '001.fits', '002.fits', '003.fits', '004.fits']]

    # Loop over images skipping metadata and first 5 submissions
    for sdoc in subs_to_score: 
        idx = int(sdoc['file'].replace(".fits", ""))
        truth_sats = truth_hdul[f"SAT_{idx:04d}"].data
        n_sub += len(sdoc['sats'])
        n_true += len(truth_sats)
        sub_associated = [False]*n_sub
        # Find best submitted satellite for each true satellite
        # Keep track of whether or not each true sat is found, and also
        # if each submitted sat is associated with a true sat.
        # Note, it *is* possible for a single submitted satellite to get
        # associated with multiple true satellites.

        # Exit early if no satellites are present.
        if len(sdoc['sats']) == 0:
            continue

        # Start by loading all submitted sats into numpy arrays
        sdict = {}
        for k in ['x0', 'x1', 'y0', 'y1', 'flux']:
            sdict[k] = []
        for sat in sdoc['sats']:
            for k in ['x0', 'x1', 'y0', 'y1', 'flux']:
                sdict[k].append(sat[k])
        for k in ['x0', 'x1', 'y0', 'y1', 'flux']:
            sdict[k] = np.array(sdict[k])

        for x0, y0, x1, y1, flux in truth_sats:
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
                pos_se += min_sq_error
                mag_se += np.square((sdict['flux'][idx] - flux)/flux*(2.5/np.log(10)))
                n_true_found += 1
                sub_associated[idx] = True
        n_sub_assoc += np.sum(sub_associated)

    print(f"completeness: {n_true_found/n_true}")
    print(f"fpr: {1-(n_sub_assoc/n_sub)}")
    print(f"pos_rmse: {np.sqrt(pos_se/n_true_found)}")
    print(f"mag_rmse: {np.sqrt(mag_se/n_true_found)}")

    scores = pd.DataFrame({'comp_name': competitor_name, 'sub_filename': args.submission, 
                            'date': datetime.utcnow().isoformat(), 'completeness': n_true_found/n_true, 
                            'fpr': 1-(n_sub_assoc/n_sub), 'pos_RMSE': np.sqrt(pos_se/n_true_found), 
                            'mag_RMSE': np.sqrt(mag_se/n_true_found)}, index=[0])

    # Save scores to file for the dashboard
    if path.exists('score_history_target.txt'):
        scores.to_csv('score_history_target.txt', mode='a', index=False, header=False)
    else:
        scores.to_csv("score_history_target.txt", index=False)


if __name__ == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("submission")
    parser.add_argument("truth")
    args = parser.parse_args()

    check_submission(args)
    score_submission(args)
