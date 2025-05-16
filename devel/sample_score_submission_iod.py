import yaml
import numpy as np
from datetime import datetime
import astropy.io.fits as fits
from astropy.stats import mad_std
import ssa
import math
import progressbar
import uncertainty_metrics
from ssa import utils
import pdb
import pandas as pd
from os import path


prop_times = [5, 15, 30, 60]


def score_submission_rmse(args, submission):
    truth_hdul = fits.open(args.truth)

    print(f"scoring {len(submission)} submissions on rmse")

    # Start with RMSE of x, y, z, vx, vy, vz
    names = ['rx', 'ry', 'rz', 'vx', 'vy', 'vz']
    allnames = [n+'sub' for n in names] + [n+'true' for n in names]
    data = np.zeros(len(submission), dtype=[(n, 'f8') for n in allnames])

    for i, sdoc in enumerate(submission):
        idx = int(sdoc['file'].replace(".csv", "").replace('SAT_', ''))
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


def score_submission_our_loc(args, submission):
    truth_hdul = fits.open(args.truth)

    print(f"scoring {len(submission)} submissions on location using our propagator")
    
    totals = {}
    for t in prop_times:
        totals[t] = {}
        totals[t]["within_distance"]=0
        totals[t]["total_sats"]=0

    bar = progressbar.ProgressBar(maxval=len(submission)).start()

    for sub_id, sdoc in enumerate(submission[1:]):
        bar.update(sub_id)
        
        idx = int(sdoc['file'].replace(".csv", "").replace('SAT_', ''))
        truthrv = truth_hdul[f"SAT_{idx:04d}"].data

        #Get truth values
        truth_v = (truthrv['vx'][0], truthrv['vy'][0], truthrv['vz'][0])
        truth_r = (truthrv['rx'][0], truthrv['ry'][0], truthrv['rz'][0])
        truth_t = truthrv['t'][0]

        #Get submission values
        sub_v = (sdoc['IOD'][0]['vx'], sdoc['IOD'][0]['vy'], sdoc['IOD'][0]['vz'])
        sub_r = (sdoc['IOD'][0]['rx'], sdoc['IOD'][0]['ry'], sdoc['IOD'][0]['rz'])
        sub_t = sdoc['IOD'][0]['t']

        #Determine orbits
        truth_orbit = ssa.orbit.Orbit(truth_r,truth_v,truth_t)
        sub_orbit = ssa.orbit.Orbit(sub_r,sub_v,sub_t)

        r_earth = 6.371e6 #meters

        for t in prop_times:
            #Compute position of propagated IOD
            r_truth = ssa.compute.rv(truth_orbit,truth_t+(t*60),propagator=ssa.propagator.SGP4Propagator())[0]
            try:
                r_sub = ssa.compute.rv(sub_orbit,sub_t+(t*60),propagator=ssa.propagator.SGP4Propagator())[0]
            except ValueError:
                print('Hit ValueError on ssa.compute.rv. Removing satellite and continuing.')
                continue

            totals[t]["total_sats"]+=1
            r_ground = ssa.utils.normed(r_truth)*r_earth
            #Check if predicted location is within .5 deg of truth position (for telescope FOV)
            if (ssa.utils.unitAngle3(ssa.utils.normed(r_truth), ssa.utils.normed(r_sub - r_ground))*180/np.pi)<=.5:
                totals[t]["within_distance"]+=1

    scores = {}
    for t in prop_times:
        scores[t] = totals[t]["within_distance"]/totals[t]["total_sats"]

    return scores

def score_submission_their_loc(args, submission):
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
        truthrv = truth_hdul[f"SAT_{idx:04d}"].data

        # Get truth values
        truth_v = (truthrv['vx'][0], truthrv['vy'][0], truthrv['vz'][0])
        truth_r = (truthrv['rx'][0], truthrv['ry'][0], truthrv['rz'][0])
        truth_t = truthrv['t'][0]

        # Determine orbits
        truth_orbit = ssa.orbit.Orbit(truth_r,truth_v,truth_t)

        r_earth = 6.371e6  # meters
        propagator = ssa.propagator.SGP4Propagator()

        for t in prop_times:
            # Compute position of propagated IOD
            r_truth = ssa.compute.rv(truth_orbit, truth_t+(t*60),
                                     propagator=propagator)[0]
            predloc = sdoc['predicted_downrange_location'][0][t]
            preduv = ssa.utils.lb2unit(np.radians(predloc['ra']),
                                       np.radians(predloc['dec']))
            totals[t]["total_sats"]+=1

            # Check if predicted location is within .5 deg of truth position
            # (for telescope FOV)
            r_ground = preduv*r_earth
            from ssa.utils import normed, unitAngle3
            dang = unitAngle3(normed(r_truth - r_ground), preduv)*180/np.pi
            if dang <= 0.5:
                totals[t]["within_distance"] += 1
            rtrue,  dtrue = ssa.utils.unit2lb(r_truth - r_ground)
            cd = np.cos(np.radians(predloc['dec']))
            raresid = ((((predloc['ra'] - np.degrees(rtrue))+180)%360)-180)*cd
            decresid = predloc['dec'] - np.degrees(dtrue)
            resid = np.array([raresid, decresid])
            chi2 = resid.dot(np.linalg.inv(sdoc['predicted_location_covariance'][0][t])).dot(resid)
            chi2dict[t].append(chi2)

    loc_covar_scores = dict()
    for t in prop_times:
        ndof = 2
        cvmp = uncertainty_metrics.cvm_chi2_test(
            np.array(chi2dict[t]), ndof, alpha=True)
        nquant = 5
        zquant = uncertainty_metrics.pearsons_chi(
            np.array(chi2dict[t]), ndof, nquant)
        loc_covar_scores[t] = cvmp  # + (z for z in zquant)

    scores = {}
    for t in prop_times:
        scores[t] = totals[t]["within_distance"]/totals[t]["total_sats"]
        scores[str(t)+'_cov'] = loc_covar_scores[t]
    
    return scores

def score_submission_covar(args, submission):
    
    #############################################################################################
    # Please refer to the appendix of README_iod.md for more information on scoring covariances #
    #############################################################################################

    print(f"scoring {len(submission)} submissions on covariance")
    truth_hdul = fits.open(args.truth)

    names = ['rx', 'ry', 'rz', 'vx', 'vy', 'vz']
    allnames = [n+'sub' for n in names] + [n+'true' for n in names]
    dtype = [(n, 'f8') for n in allnames]
    dtype += [('state_vector_covariance', 'f4', (6, 6))]
    data = np.zeros(len(submission), dtype=dtype)

    for i, sdoc in enumerate(submission):    # Loop over IODs
        idx = int(sdoc['file'].replace(".csv", "").replace('SAT_', ''))
        truthrv = truth_hdul[f"SAT_{idx:04d}"].data
        data['state_vector_covariance'][i] = sdoc['state_vector_covariance'][0]
        for n in names:
            data[n+'true'][i] = truthrv[n]
            # not possible for there to be more than one OD in a file?
            data[n+'sub'][i] = sdoc['IOD'][0][n]
    resid = np.array([data[n+'sub']-data[n+'true'] for n in names]).T
    chi2 = uncertainty_metrics.chi2(resid, data['state_vector_covariance'])
    cvmp = uncertainty_metrics.cvm_chi2_test(chi2, 6, alpha=True)
    nquant = 5
    zquant = uncertainty_metrics.pearsons_chi(chi2, 6, nquant)

    scores = dict()
    scores['cvmp'] = cvmp
    for i in range(nquant):
        scores[f'zq{i}'] = zquant[i]

    return scores

if __name__ == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("submission")
    parser.add_argument("truth")
    args = parser.parse_args()
    submission = [s for s in yaml.safe_load_all(open(args.submission, 'r'))]

    competitor_name = submission[0]['competitor_name']
    competitor_name_split = competitor_name.split()
    if len(competitor_name_split) > 1:
        competitor_name = "_".join(competitor_name_split)

    subs_to_score = [s for s in submission[1:] if s['file'] in ['SAT_0000.csv', 'SAT_0001.csv', 'SAT_0002.csv', 'SAT_0003.csv', 'SAT_0004.csv']]
    # do not score metadata entry or first five example entries
    rmsedict = score_submission_rmse(args, subs_to_score)
    scores_our_loc = score_submission_our_loc(args, subs_to_score)

    #initialize scores_their_loc
    scores_their_loc = {}
    state_covar = {}
    for k in range(6):
        state_covar[k] = float('inf')
    for t in prop_times:
        scores_their_loc[str(t)+"_cov"] = float('inf')
        scores_their_loc[t] = float('inf')
    if subs_to_score[0].get('predicted_downrange_location') is not None:
        if subs_to_score[0].get('predicted_location_covariance') is not None:
            scores_their_loc = score_submission_their_loc(args, subs_to_score)
        else:
            print('Predicted locations provided, but no covariance given. Score = 0.0')
    else:
        print('No predicted downrange locations and covariance given. Score = 0.0')
    if subs_to_score[0].get('state_vector_covariance') is not None:
        state_covar = score_submission_covar(args, subs_to_score)
    else:
        print('No state vector covariance given. Score = 0.0')

    names = ['rx', 'ry', 'rz', 'vx', 'vy', 'vz']

    scores_header = pd.DataFrame({'comp_name': competitor_name, 'sub_filename': args.submission, 'date': datetime.utcnow().isoformat()}, index=[0])
    
    dict_name = []
    dict_data = []

    for n in names:
        # Append RMSE scores
        dict_name.append('RMSE_'+str(n))
        dict_data.append(rmsedict[n])
    for t in prop_times:
        # Append scores based on our location
        dict_name.append('our_prop_'+str(t))
        dict_data.append(scores_our_loc[t])
    for t in prop_times:
        # Append scores based on their submitted locations
        dict_name.append('prop_'+str(t))
        dict_data.append(scores_their_loc[t])
    for t in prop_times:
        # Append covariance scores based on their submitted locations
        dict_name.append('cvm_'+str(t))
        dict_data.append(scores_their_loc[str(t)+"_cov"])
    for k in state_covar:
        # Append state_covar scores
        dict_name.append(str(k))
        dict_data.append(state_covar[k])
    
    scores = pd.DataFrame([dict_data], columns=dict_name, index=[0])
    final_scores = pd.concat([scores_header, scores], axis=1)

    # Save scores to file for the dashboard
    if path.exists('score_history_iod.txt'):
        final_scores.to_csv('score_history_iod.txt', mode='a', index=False, header=False)
    else:
        final_scores.to_csv("score_history_iod.txt", index=False)

