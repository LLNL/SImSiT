# **Branch Overview**
This branch is derived from the ODJobs data, and consists of 445 tracks of optical, angle-only observations of varying lengths with noise properties consistent with typical target-tracking images. This branch is designed to test the accuracy with which algorithms can make initial orbit determinations (IOD), inferring the position and velocity of a satellite from optical-only observations over a fraction of an orbital period. Truth information is provided for the first 5 satellites, therefore competitor scores are calculated on the remaining 440 satellite IODs.

# **Data and Formatting**

This directory contains the following:
- ***iod_datasets***:
    - A folder containing 445 catalogs in comma seperated variable (.csv) format, one catalog for each satellite (named 0000 through 0444), with the following columns:
    - **t**: The time of that observation in GPS seconds. (max(t), the last observation, is the time that you want to calculate the IOD at). 
    - **ra**: The right ascension of the satellite during the given observation in degrees (topocentric)
    - **dec**: The declination of the satellite during the given observation in degrees (topocentric)
    - **px**: The x position (in meters, GCRF) of the observer at that observation
    - **py**: The y position (in meters, GCRF) of the observer at that observation 
    - **pz**: The z position (in meters, GCRF) of the observer at that observation
- ***truth_5.fits***: 
    - This is a FITS format file with 5 extensions inside (each titled "SAT_[XXXX]"), corresponding to the truth data from the first 5 .csv catalogs.
    - Each extension contains the satellites' IOD at the time of the last observation in the data file for each satellite:
        - **rx**: x component of the satellite location in meters (GCRF)
        - **ry**: y component of the satellite location in meters (GCRF)
        - **rz**: z component of the satellite location in meters (GCRF)
        - **vx**: x component of the satellite velocity in meters/second (GCRF)
        - **vy**: y component of the satellite velocity in meters/second (GCRF)
        - **vz**: z component of the satellits velocity in meters/second (GCRF)
        - **t**: Time of the epoch in GPS seconds. (This should match the time of the last observation in the data file for that satellite.)
- ***sample_submission_5.yaml*** 
    - This is a .yaml file that shows what a sample submission to be given back to us should look like (including optional fields).
- ***score_submission.py***:
    - This script will test your submission file format, and will score the first 5 submissions. You can run this with:

        `python score_submission.py [your_submission_file] truth_5.fits`

# **Submission Format**

The format for IOD submissions should include:
- **branch**: The branch that you intend to submit the submission to for scoring. ('iod')
- **competitor_name**: The name you would like us to use to identify your submission file. This will be the name displayed on the dashboard.
- **display_true_name**: A boolean deciding whether you would like your `competitor_name` publically displayed on the dashboard. If false, your score will be assigned an anonymous `competitor_name`.
- **IOD**: 
  - the satellite state (r, v) composed each of 3 component tuples (x, y, and z components) in GCRF coordinates and MKS units.
  - epoch for each satellite in units of GPS seconds. This should match the time of the last observation in the corresponding satellite's catalog.
- **file**: The current catalog filename
- **predicted_downrange_location (optional)**: 
  - Predicted location (ra and dec) of the satellite downrange 5, 15, 30 and 60 minutes after the IOD epoch.
- **predicted_location_covariance (optional)**:
  - covariances describing the uncertainties associated with future downrange angular positions. Covariance format is described in the **appendix**.
- **state_vector_covariance (optional)**:
  - covariances describing the uncertainties associated with the predicted satellite state. Covariance format is described in the **appendix**.

The following is an example of what the IOD output should look like, before being submitted for scoring:
   
    branch: IOD
    competitor_name: Competitor A
    display_true_name: true
    ---
    IOD:
    - rx: -1537532.386353053
      ry: -6096250.110596744
      rz: -2932138.755877168
      t: 1261905338.001
      vx: -311.47745042782947
      vy: 3343.468207729812
      vz: -6799.62176813366
    file: SAT_0000.csv
    predicted_downrange_location:
    - 5:
        dec: -43.53189885830184
        ra: -107.91569546166455
      15:
        dec: -78.58191049565286
        ra: -142.54281274857303
      30:
        dec: -41.828333297984834
        ra: 86.42787859450476
      60:
        dec: 68.6751345185049
        ra: 59.241353858585796
    predicted_location_covariance:
    - 5:
      - - 0.0015969151712406703
        - 0.001675564215813318
      - - 0.001675564215813318
        - 0.0018524397935565096
      15:
      - - 0.001254350083664949
        - 9.23377673931266e-05
      - - 9.23377673931266e-05
        - 0.0003252086528691818
      30:
      - - 0.006958373179426201
        - -0.038822584981181865
      - - -0.038822584981181865
        - 0.22495388818961692
      60:
      - - 0.9427673378963979
        - -2.405701377607157
      - - -2.405701377607157
        - 6.150014585774827
    state_vector_covariance:
    - - - 16436394.221447604
        - 2736349.250020679
        - 5185419.779898066
        - 9556.914609126696
        - -31237.786818765395
        - 80851.45467810427
      - - 2736349.250020679
        - 457110.06632237736
        - 862903.6326622862
        - 1620.25695961848
        - -5171.651501185036
        - 13450.424836641443
      - - 5185419.779898067
        - 862903.6326622862
        - 1638434.6306020964
        - 2935.4533322657076
        - -9889.659954989256
        - 25550.432948030717
      - - 9556.914609126696
        - 1620.25695961848
        - 2935.4533322657076
        - 11.287088246693255
        - -15.659582604787422
        - 45.05632880534354
      - - -31237.786818765395
        - -5171.651501185036
        - -9889.659954989254
        - -15.65958260478742
        - 60.717231465657946
        - -154.52044794196547
      - - 80851.45467810426
        - 13450.424836641443
        - 25550.432948030713
        - 45.05632880534353
        - -154.52044794196547
        - 398.6250618811772
    ---
    IOD:
      ...
    ...

Where the number after "SAT" corresponds to the same number of the catalog that was used to determine the IOD.

# ***Scores:***
The summary score that will be reflected on the main page of the dashboard for this branch is `our_prop_30` (This is the percent of predicted satellites that are within the telescope FOV after 30 minutes using our propagator).

**The following are the score entires in `score_history.txt`, by index:**
- **[comp_name]**: Competitor name
- **[branch]**: Branch you wish to submit to
- **[sub_filename]**: Submission filename
- **[date]**: Date of submission
- **[rx_RMSE]**: Root Mean Square Error of rx
- **[ry_RMSE]**: Root Mean Square Error of ry
- **[rz_RMSE]**: Root Mean Square Error of rz
- **[vx_RMSE]**: Root Mean Square Error of vx
- **[vy_RMSE]**: Root Mean Square Error of vy
- **[vz_RMSE]**: Root Mean Square Error of vz
- **[our_prop_5]**: Percent of predicted satellite locations within telescope FOV 5 minutes downrange using our propagator
- **[our_prop_15]**: Percent of predicted satellite locations within telescope FOV 15 minutes downrange using our propagator
- **[our_prop_30]**: Percent of predicted satellite locations within telescope FOV 30 minutes downrange using our propagator
- **[our_prop_60]**: Percent of predicted satellite locations within telescope FOV 60 minutes downrange using our propagator

***Note***: **The following scores may be "inf" if you have not supplied us the required information for scoring**
- **[prop_5]**: Percent of predicted satellite locations within telescope FOV 5 minutes downrange using your supplied downrange ra and dec
- **[prop_15]**: Percent of predicted satellite locations within telescope FOV 15 minutes downrange using your supplied downrange ra and dec
- **[prop_30]**: Percent of predicted satellite locations within telescope FOV 30 minutes downrange using your supplied downrange ra and dec
- **[prop_60]**: Percent of predicted satellite locations within telescope FOV 60 minutes downrange using your supplied downrange ra and dec

**The following 4 columns are Cramer von-Mises probabilities (CvM)**
- **[cvm_5]**: CvM for predicted location 5 minutes downrange using your supplied downrange ra and dec
- **[cvm_15]**: CvM for predicted location 15 minutes downrange using your supplied downrange ra and dec
- **[cvm_30]**: CvM for predicted location 30 minutes downrange using your supplied downrange ra and dec
- **[cvm_60]**: CvM for predicted location 60 minutes downrange using your supplied downrange ra and dec
- **[cvmp]**: CvM for the accuracy of the state estimate (rx, ry, rz, vx, vy, vz), giving the likelihood that the true observations were drawn from Gaussian distributions centered at submitted locations with submitted covariances

**The last 5 columns are the pulls which would be used to compute a Pearson chi^2 statistic for the submission. Each column is a z-score for the number of estimates that have a chi^2 in the bin indicated:**
- **[zq0]**: zq for the 0-20% percentile of the distribution
- **[zq1]**: zq for the 20-40% percentile of the distribution
- **[zq2]**: zq for the 40-60% percentile of the distribution
- **[zq3]**: zq for the 60-80% percentile of the distribution
- **[zq4]**: zq for the 80-100% percentile of the distribution


# **Appendix**

## ***Submission Files***:
- **Predicted_location_covariance**: For IOD downrange angle predictions, the matrix should be a two by two matrix with entries in the following order: [ra, dec] (zero indexed) with units in degrees. Repeat for each downrange time interval.
    - The covariance in the ra direction should be scaled like cos(dec), so that an uncertainty ellipse of the same angular size should correspond to a convariance matrix with the same determinant anywhere on the sky.
- **state_vector_covariance**: For IOD state predictions, the matrix should be a six by six matrix with entries in the following order: [x, y, z, vx, vy, vz] (zero indexed) with units compatible with the state parameters.
    - Ex: entry [2, 4] in the matrix should be the covariance of the z estimate with the vy estimate, with units of m * (m/s).

## ***Scoring***:
- Covariance matrices are scored based on the Cramer-von Mises (CvM) statistic and Pearson chi squared statistics, following Horwood et al. (2014). 
    - Chi^2 is computed for each reported position given the ground truth and the estimated covariance matrix. The resulting chi^2 distribution is compared with a chi^2 distribution with two degrees of freedom for the downrange angular predictions.
      - The Cramer-von Mises statistic is converted to a likelihood and used to assess how compatible these distributions are:
          - Scores range from .001 (poor match), to 1 (perfect match).
    - For state predictions we report z-scores for the fraction of chi^2 landing in bins corresponding to chi^2 in the following percentiles: [0, 0.2], [0.2, 0.4], [0.4, 0.6], [0.6, 0.8], [0.8, 1.0]. 
        - Scores of 0 are good
        - Positive and negative values indicate that more or fewer values that expected are falling in those bins.