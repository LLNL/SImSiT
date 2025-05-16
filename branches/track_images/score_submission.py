import yaml
import numpy as np
from datetime import datetime
import astropy.io.fits as fits
from astropy.stats import mad_std
import math
import progressbar
from xfiles import uncertainty_metrics
from ssapy import utils
import pdb
import ssapy
import pandas as pd
from os import path
import glob
from astropy.table import Table
import os


prop_times = [5, 15, 30, 60]

def check_submission(args):
    # Let yaml raise exception here
    print(f"Opening file: {args.submission}")
    with open(args.submission, 'r') as f:
        docs = yaml.safe_load_all(f)
        docs = [doc for doc in docs]

    print(f"Found {len(docs)-1} submissions.")
    # Checking metadata
    assert isinstance(docs[0]['branch'], str)
    assert isinstance(docs[0]['competitor_name'], str)
    assert isinstance(docs[0]['display_true_name'], bool)

    branch = docs[0]['branch']
    if (branch != "target_track") and (branch != "sidereal_track"):
        print("ERROR: Invalid branch name")
        sys.exit()

    # Do the format check first, since it's fast. Starting after metadata.
    # We're requiring here that the order of the entries
    # in the file matches the order we saved in truth_times.
    for doc in docs[1:]:
        # Let python raise error if key not found.
        assert isinstance(doc['file'], str)
        if len(doc['IOD']) > 1:
            raise AssertionError('multiple ODs in one IOD submission?')
        for sat in doc['IOD']:
            for var in ['rx', 'ry', 'rz', 't', 'vx', 'vy', 'vz']:
                assert isinstance(sat[var], float)
    return branch
    
def score_submission_rmse(args, submission):
    truth_hdul = fits.open(args.truth)

    print(f"scoring {len(submission)} submissions on rmse")

    # Start with RMSE of x, y, z, vx, vy, vz
    names = ['rx', 'ry', 'rz', 'vx', 'vy', 'vz']
    allnames = [n+'sub' for n in names] + [n+'true' for n in names]
    data = np.zeros(len(submission), dtype=[(n, 'f8') for n in allnames])

    for i, sdoc in enumerate(submission):
        idx = int(sdoc['file'].replace('SAT_', ''))
        truthrv = truth_hdul[f"SAT_{idx:04d}"].data
        for n in names:
            data[n+'true'][i] = truthrv[n]
            data[n+'sub'][i] = sdoc['IOD'][0][n]
    rmsedict = dict()
    madstddict = dict()
    for n in names:
        rmse = np.sqrt(np.sum((data[n+'sub']-data[n+'true'])**2)/len(data))
        madstd = mad_std(data[n+'sub']-data[n+'true'])
        rmsedict[n] = rmse
        madstddict[n] = madstd
        print(f'{n}_rmse: {rmse}')
        print(f'{n}_madstd: {madstd}')

    return rmsedict


def score_submission_loc(args, submission):
    truth_hdul = fits.open(args.truth)

    print(f"scoring {len(submission)} submissions on location, using submitted predicted locations")

    totals = {}
    for t in prop_times:
        totals[t] = {}
        totals[t]["within_distance"] = 0
        totals[t]["total_sats"] = 0
    chi2dict = {t: list() for t in prop_times}


    bar = progressbar.ProgressBar(maxval=len(submission)).start()

    for sub_id, sdoc in enumerate(submission):
        bar.update(sub_id)
    
        idx = int(sdoc['file'].replace(".csv", "").replace('SAT_', ''))
        truth_rvt = truth_hdul[f"SAT_{idx:04d}"].data

        # Get truth values
        truth_v = (truth_rvt['vx'][0], truth_rvt['vy'][0], truth_rvt['vz'][0])
        truth_r = (truth_rvt['rx'][0], truth_rvt['ry'][0], truth_rvt['rz'][0])
        truth_t = truth_rvt['t'][0]

        #Get submission values
        sub_v = (sdoc['IOD'][0]['vx'], sdoc['IOD'][0]['vy'], sdoc['IOD'][0]['vz'])
        sub_r = (sdoc['IOD'][0]['rx'], sdoc['IOD'][0]['ry'], sdoc['IOD'][0]['rz'])
        sub_t = sdoc['IOD'][0]['t']

        #Determine orbits
        truth_orbit = ssapy.orbit.Orbit(truth_r,truth_v,truth_t)
        sub_orbit = ssapy.orbit.Orbit(sub_r,sub_v,sub_t)

        r_earth = 6.371e6  # meters

        for t in prop_times:
            #Compute position of propagated IOD
            r_truth = ssapy.compute.rv(truth_orbit,truth_t+(t*60),propagator=ssapy.propagator.SGP4Propagator())[0]
            try:
                r_sub = ssapy.compute.rv(sub_orbit,sub_t+(t*60),propagator=ssapy.propagator.SGP4Propagator())[0]
            except ValueError:
                print('Hit ValueError on ssapy.compute.rv. Removing satellite and continuing.')
                continue

            totals[t]["total_sats"]+=1
            r_ground = ssapy.utils.normed(r_truth)*r_earth
            #Check if predicted location is within .5 deg of truth position (for telescope FOV)
            if (ssapy.utils.unitAngle3(ssapy.utils.normed(r_truth), ssapy.utils.normed(r_sub - r_ground))*180/np.pi)<=.5:
                totals[t]["within_distance"]+=1

    scores = {}
    for t in prop_times:
        scores[t] = totals[t]["within_distance"]/totals[t]["total_sats"]
    
    return scores

if __name__ == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("submission")
    parser.add_argument("truth")
    args = parser.parse_args()

    branch = check_submission(args)

    submission = [s for s in yaml.safe_load_all(open(args.submission, 'r'))]

    competitor_name = submission[0]['competitor_name']
    competitor_name_split = competitor_name.split()
    if len(competitor_name_split) > 1:
        competitor_name = "_".join(competitor_name_split)

    truth_hdul = fits.open(args.truth)
    if int((len(truth_hdul)-1)/2) < 5:
        print("Scoring < 5 truth submissions.")
        subs_to_score = [s for s in submission[1:]]
    elif int((len(truth_hdul)-1)/2) == 5:
        print("Scoring 5 truth submissions.")
        subs_to_score = [s for s in submission[1:] if s['file'] in ['SAT_0000.', 'SAT_0001', 'SAT_0002', 'SAT_0003.', 'SAT_0004']]
    else: 
        print(f"Scoring {len(submission)-6} submissions.")
        subs_to_score = [s for s in submission[1:] if s['file'] not in ['SAT_0000.', 'SAT_0001', 'SAT_0002', 'SAT_0003', 'SAT_0004']]
    # do not score metadata entry or first five example entries
    rmsedict = score_submission_rmse(args, subs_to_score)
    scores_loc = score_submission_loc(args, subs_to_score)

    names = ['rx', 'ry', 'rz', 'vx', 'vy', 'vz']

    scores_header = pd.DataFrame({'comp_name': competitor_name, 'branch': branch, 'sub_filename': args.submission, 'date': datetime.utcnow().isoformat()}, index=[0])
    
    dict_name = []
    dict_data = []

    for n in names:
        # Append RMSE scores
        dict_name.append(str(n)+"_RMSE")
        dict_data.append(rmsedict[n])
    for t in prop_times:
        # Append scores based on their submitted locations
        dict_name.append('our_prop_'+str(t))
        dict_data.append(scores_loc[t])
    
    scores = pd.DataFrame([dict_data], columns=dict_name, index=[0])
    final_scores = pd.concat([scores_header, scores], axis=1)

    # Save scores to file for the dashboard
    if path.exists('score_history.txt'):
        final_scores.to_csv('score_history.txt', mode='a', index=False, header=False)
    else:
        final_scores.to_csv("score_history.txt", index=False)
