import yaml
import os
import glob
import numpy as np
from astropy.table import Table

def check_submission(args):
    # Let yaml raise exception here
    print(f"Opening file: {args.infile}")
    with open(args.infile, 'r') as f:
        docs = yaml.safe_load_all(f)
        docs = [doc for doc in docs]

    files = glob.glob('public/iod_datasets/SAT_*.csv')
    truth_times = {}
    for fn in files:
        tab = Table.read(fn)
        truth_times[os.path.basename(fn)] = np.max(tab['time'])

    print(f"Found {len(docs)-1} submissions.")
    # Checking metadata
    assert isinstance(docs[0]['branch'], str)
    assert isinstance(docs[0]['competitor_name'], str)
    assert isinstance(docs[0]['display_true_name'], bool)

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
                assert (sat['t'] -
                        truth_times[doc['file']] <= .001)
        if doc.get('predicted_downrange_location') is not None:
            for loc in doc['predicted_downrange_location']:
                for time in [5, 15, 30, 60]:
                    for pos in ['ra', 'dec']:
                        assert isinstance(loc[time][pos], float)
        else:
            print('No predicted downrange locations supplied')
        if doc.get('predicted_location_covariance') is not None:
            for cov in doc['predicted_location_covariance']:
                for time in [5, 15, 30, 60]:
                    assert isinstance(cov[time], list)
                    for cov_1 in cov[time]:
                        assert isinstance(cov_1, list)
                        for cov_2 in cov_1:
                            assert isinstance(cov_2, float)
                
        else:
            print('No predicted location covariances supplied')
        if doc.get('state_vector_covariance') is not None:
            for cov in doc['state_vector_covariance']:
                assert isinstance(cov, list)
                for val in cov:
                    assert isinstance(val, list)
                    for val_2 in val:
                        assert isinstance(val_2, float)
        else:
            print('No state vector covariances supplied')


if __name__ == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("infile")
    parser.add_argument("--plot", action='store_true')
    parser.add_argument("--public_directory", default="public")
    args = parser.parse_args()
    check_submission(args)
