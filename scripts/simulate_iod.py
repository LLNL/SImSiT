import os
import random
import pandas as pd
import numpy as np
from astropy.time import Time
import astropy.units as u
import astropy.io.fits as fits
import yaml
import progressbar
import ssapy


def parse_float(x):
    try:
        x = float(x)
    except ValueError:
        x = float('nan')
    return x


def parse_odjob(filename):
    odjob_dtype = np.dtype([
        ('time', np.float64),
        ('tag', np.int32),
        ('sensor', np.int32),
        ('ra', np.float64),
        ('dec', np.float64),
        ('px', np.float64),
        ('py', np.float64),
        ('pz', np.float64),
        ('prov', 'S20'),
        ('r', np.float64),
        ('rr', np.float64),
        ('rcs', np.float64),
        ('snr', np.float64),
        ('mv', np.float64),
    ])
    with open(filename, 'r') as fd:
        lines = fd.readlines()[2:]  # Skip first 2 lines
    arr = np.empty(len(lines), dtype=odjob_dtype)
    for i, line in enumerate(lines):
        line = line.replace('[', '')  # Get rid of '[' and ']' to make normal CSV
        line = line.replace(']', '')
        time, tag, sensor, ra, dec, px, py, pz, prov, r, rr, rcs, snr, mv = line.split(',')
        timeunix = float(time)
        tag = int(tag)
        sensor = int(sensor)
        ra = float(ra)
        dec = float(dec)
        px = float(px)
        py = float(py)
        pz = float(pz)
        r = parse_float(r)
        rr = parse_float(rr)
        rcs = parse_float(rcs)
        snr = parse_float(snr)
        mv = parse_float(mv)
        arr[i] = timeunix, tag, sensor, ra, dec, px, py, pz, prov, r, rr, rcs, snr, mv
    arr['time'] = Time(arr['time'], format='unix').gps
    return arr


def read_tle(sat_name, tle_filename):
    """
    Get the TLE data from the file for the satellite with the given name

    :param sat_name: NORAD name of the satellite
    :type sat_name: str

    :param tle_filename: Path and name of file where TLE is
    :type tle_filename: str

    :return: Both lines of the TLE for the satellite
    :rtype: (str, str)

    :raises IOError: when TLE file is invalid
    :raises KeyError: when satellite is not found in TLE file
    """
    with open(tle_filename) as tle_f:
        lines = [l.rstrip() for l in tle_f.readlines()]

    index_found = False
    for idx, line in enumerate(lines):
        if line.split()[1] == str(sat_name)+'U':
            sat_line_ind = idx
            index_found = True
            break
    if index_found is False:
        raise KeyError(
            "No satellite '{}' in file '{}'".format(sat_name, tle_filename))
    try:
        line1 = lines[sat_line_ind]
        line2 = lines[sat_line_ind + 1]
        return line1, line2
    except IndexError:
        raise IOError("Incorrectly formatted TLE file")


def main():
    # Make the public and private directories to store the final .csv files
    for d in ["branches/iod/public/iod_datasets", "branches/iod/private"]:
        if not os.path.exists(d):
            os.makedirs(d)

    # Get all of the catalog names from ODJobs and remove the readme from the
    # list
    catalogs = os.listdir(path='branches/iod/ODJobs_Simulated_Data_20001_20005_v2/')
    catalogs.remove('_README.md')
    print('Total catalogs in ODJobs: ' + str(len(catalogs)))

    # Get only the files that start with "20001_" because that is the truth
    # file we have
    # Also use optical only catalogs
    catalogs_new = []

    for file in catalogs:
        cat_mask = file.startswith('20001_') and '_opt_' in file
        if cat_mask == True:
            catalogs_new.append(file)
    print('Number of 20001 catalogs in ODJobs: ' + str(len(catalogs_new)))

    # Randomly shuffle catalogs, with seed 28
    random.Random(28).shuffle(catalogs_new)

    print('Final catalogs in dataset: ' + str(len(catalogs_new)))

    # Define an empty array to save satellite ID's for looking up in truth
    # table
    sat_id = []
    max_time = []

    print('Saving truth sat IDs')

    # Get the relevant information from each satellite and save as a .csv
    # under private
    for idx, file in enumerate(catalogs_new):
        data = parse_odjob('branches/iod/ODJobs_Simulated_Data_20001_20005_v2/'+str(file))
        sat_id.append(data['tag'][0])
        max_time.append(np.max(data['time']))

        time = data['time']
        ra = data['ra']
        dec = data['dec']
        px = data['px']*1000
        py = data['py']*1000
        pz = data['pz']*1000

        df = pd.DataFrame({'time': time, 'ra':ra, 'dec':dec, 'px':px, 'py':py, 'pz':pz})
        df.to_csv('branches/iod/public/iod_datasets/SAT_'+ "%04d"% (idx) + '.csv', index_label=False, index=False)

    # Save the satellite ID's privately
    np.savetxt("branches/iod/private/sat_id_list.csv", np.array(sat_id), fmt='%d',
               delimiter=",")
    np.savetxt("branches/iod/private/odjob_list.csv", np.array(catalogs_new), fmt='%s',
               delimiter=',')

    truth5_hdul = fits.HDUList()
    truth_hdul = fits.HDUList()
    rng = np.random.default_rng(28)
    sample_docs = [{'branch' : 'IOD', 'competitor_name' : 'Competitor A', 'display_true_name' : True}]

    print('Calculating IODs from TLEs')

    bar = progressbar.ProgressBar(maxval=len(catalogs_new)).start()

    for idx, sat in enumerate(sat_id):
        bar.update(idx)
        # read the TLE for each satellite
        tle_tuple = read_tle(sat, 'branches/iod/20001_tim_tle.txt')
        # Calculate r, v, t for each satellite
        orbit = ssapy.Orbit.fromTLETuple(tle_tuple)
        r_new, v_new = ssapy.rv(orbit, max_time[idx],
                              propagator=ssapy.SGP4Propagator())

        # Save the truth data in .fits files
        cols = fits.ColDefs([
            fits.Column(name='rx', format='E', array=np.array([r_new[0]])),
            fits.Column(name='ry', format='E', array=np.array([r_new[1]])),
            fits.Column(name='rz', format='E', array=np.array([r_new[2]])),
            fits.Column(name='vx', format='E', array=np.array([v_new[0]])),
            fits.Column(name='vy', format='E', array=np.array([v_new[1]])),
            fits.Column(name='vz', format='E', array=np.array([v_new[2]])),
            fits.Column(name='t', format='D',
                        array=np.array([max_time[idx]])),
            ])

        hdu = fits.BinTableHDU.from_columns(cols, name=f"SAT_{idx:04d}")

        if idx < 5:
            truth5_hdul.append(hdu)
        truth_hdul.append(hdu)

        sample_sat_dicts = []
        sample_sat_dicts.append({
            'rx': float(r_new[0] + rng.normal(scale=1.0)),
            'ry': float(r_new[1] + rng.normal(scale=1.0)),
            'rz': float(r_new[2] + rng.normal(scale=1.0)),
            'vx': float(v_new[0] + rng.normal(scale=1.0)),
            'vy': float(v_new[1] + rng.normal(scale=1.0)),
            'vz': float(v_new[2] + rng.normal(scale=1.0)),
            't': float(max_time[idx])
            })

        sample_docs.append({'file':f"SAT_{idx:04d}", 'IOD':sample_sat_dicts})

    with open("branches/iod/private/sample_submission_truth.yaml", "w") as f:
        yaml.safe_dump_all(sample_docs, f)

    truth_hdul.writeto("branches/iod/private/truth.fits", overwrite=True)
    truth5_hdul.writeto("branches/iod/public/truth_5.fits", overwrite=True)

if __name__ == "__main__":
    main()
