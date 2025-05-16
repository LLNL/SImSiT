import yaml
import sys

def check_submission(args):

    # Let yaml raise exception here
    print(f"Opening file: {args.infile}")
    with open(args.infile, 'r') as f:
        docs = yaml.safe_load_all(f)
        docs = [doc for doc in docs]

    # First section of yaml file is metadata. Don't include in submission count.
    print(f"Found {len(docs)-1} submissions.")
    
    # Checking metadata
    assert isinstance(docs[0]['branch'], str)
    assert isinstance(docs[0]['competitor_name'], str)
    assert isinstance(docs[0]['display_true_name'], bool)

    sat_missing_param = False
    star_missing_param_detect = False
    star_missing_param_calibrate = False

    # Checking minimum information is supplied, and verifying types.
    for idx, doc in enumerate(docs[1:]):
        # Let python raise error if key not found.
        assert isinstance(doc['file'], str)

        for idy, sat in enumerate(doc['sats']):

            sat_detect_var = ['flux', 'x0', 'x1', 'y0', 'y1']
            if all(k in sat for k in sat_detect_var):
                for var in sat_detect_var:
                    assert isinstance(sat[var], float)
                if idx == len(doc)-2 and idy == len(doc['sats'])-1: 
                    print("Sat detection parameters found.")
            else:
                # Exit the script if any sat fails to have the minimum information supplied
                print("ERROR: Sat detection parameters not found, or some parameters missing. The minimum requirements to be scored on are not met.")
                sys.exit()

            sat_calibrate_var = ['dec0', 'dec1', 'mag', 'ra0', 'ra1']
            if all(k in sat for k in sat_calibrate_var):
                for var in sat_calibrate_var:
                    assert isinstance(sat[var], float)
                if idx == len(doc)-2 and idy == len(doc['sats'])-1: 
                    if not sat_missing_param:
                        print("Sat calibration parameters found.")
            else:
                if not sat_missing_param:
                    print("No sat calibration parameters found, or some parameters missing. Scoring on detection only.")
                    sat_missing_param = True

        if 'stars' in doc:
            for idz, star in enumerate(doc['stars']):
                star_detect_var = ['flux', 'x', 'y']
                if all(k in star for k in star_detect_var):
                    if idx == len(doc)-2 and idz == len(doc['stars'])-1:
                        if not star_missing_param_detect:
                            print("Star detection parameters found.")
                    for var in star_detect_var:
                        assert isinstance(star[var], float)
                else:
                    if not star_missing_param_detect:
                        star_missing_param_detect = True
                star_calibrate_var = ['dec', 'mag', 'ra']
                if all(k in star for k in star_calibrate_var):
                    if idx == len(doc)-2 and idz == len(doc['stars'])-1:
                        if not star_missing_param_calibrate:
                            print("Star calibration parameters found.")
                    for var in star_calibrate_var:
                        assert isinstance(star[var], float)
                else:
                    if not star_missing_param_calibrate:
                        star_missing_param_calibrate = True



if __name__ == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("infile")
    args = parser.parse_args()
    check_submission(args)
