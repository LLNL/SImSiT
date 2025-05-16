# **Branch Overview**
This branch provides simulated "stacks" of images corresponding to ground-based observations for 105 different LEO satellites. A .5s exposure is taken once every 30 seconds throughout the time the satellite is overhead. 
This branch is designed to test the accuracy with which algorithms can make initial orbit determinations (IOD), inferring the position and velocity of a satellite by accurately measuring the location of endpoints for the satellite being tracked in each stack of images. Truth information is provided for the first 5 satellites, therefore competitor scores are calculated on the remaining 100 satellite IODs.

- Images include optical distortion, a variable surface brightness background, and vignetting.
- Star positions and fluxes are derived from Gaia DR2. 
  - Simulated star positions account for proper motion, parallax, and annual aberration. 
  - Simulated fluxes are derived from Gaia G-band measurements. (For simulation purposes, we use flux_i == flux_G)
- Pixel coordinate convention follows FITS:  the center of the lower left pixel is (1.0, 1.0)
- Important header keywords:
  - DATE-BEG, DATE-END - time of shutter open/close.  Shutter travel assumed to be instantaneous
  - OBSGEO-X, OBSGEO-Y, OBSGEO-Z - ITRS position of observatory
  - Approximate TAN wcs keywords: (CTYPE, CRPIX, CD, CUNIT, CRVAL)

# **Data and Formatting**

This directory contains the following:
- ***images***:
    - A folder containing simulated FITs images 105 satellite tracks. Each satellite has images for 1 overhead pass over a ground sensor, with observations taking on the format "[XXXX]_[YYY].fits" where [XXXX] is the satellite number and [YYY] is the observation number. (ex: 0003_007.fits is the 7th observation of the 3rd satellite.)
- ***sky_flat.fits***:
  - A flat field image computed by normalizing each individual challenge image, taking the median over the course of the night, and dividing out by the image mean.
- ***truth_5.fits***:
    - This is a FITS format file with 5 extensions inside, corresponding to the truth data from the first 5 satellites.
    - There are 5 extensions titled "SAT_[XXXX]" which correspond to the satellite truth data for the LAST OBSERVATION in each of the first 5 FITS image stacks, each containing the following:
        - **satID**: Satellite ID
        - **t**: The time of that observation in GPS seconds. (The time that you want to calculate the IOD at). 
        - **r**: Satellite location (rx, ry, rz) in meters (GCRF)
        - **v**: Satellite velocity (vx, vy, vz) in meters/second (GCRF)
        - **ra**: The right ascension of the satellite at time t, in degrees
        - **dec**: The declination of the satellite at time t, in degrees
        - **slant**: distance from the observer to the target, in meters
        - **pmra**: Proper motion in right ascension in radians/second
        - **pmdec**: Proper motion in declination in radians/second
        - **slantrate**: derivative of distance from observer in meters/second
        - **rStation**: Ground station location (rx, ry, rz) in meters (GCRF)
        - **vStation**: Ground station velocity (vx, vy, vz) in meters (GCRF)
        - **sensor**: Sensor ID
        - **illum**: Boolean value. True = satellite is illuminated. False = satellite in Earth's shadow.
        - **sunalt**: Altitude of the sun in radians
        - **lon**: Geodetic latitude of ground observation location in degrees
        - **lat**: Geodetic latitude in ground observation location in degrees
        - **elev**: Elevation of ground observation location in meters
        - **rx**: x component of the satellite location in meters (GCRF)
        - **ry**: y component of the satellite location in meters (GCRF)
        - **rz**: z component of the satellite location in meters (GCRF)
        - **vx**: x component of the satellite velocity in meters/second (GCRF)
        - **vy**: y component of the satellite velocity in meters/second (GCRF)
        - **vz**: z component of the satellits velocity in meters/second (GCRF)
    - There are 5 extensions titled "SAT_RVT_[XXXX]" corresponding to states where the first 5 satellites are propagated throughout time, each containing the following columns:
        - **r**: Satellite location (rx, ry, rz) in meters (GCRF)
        - **v**: Satellite velocity (vx, vy, vz) in meters/second (GCRF)
        - **t**: The time of the observation in GPS seconds.
- ***sample_submission_5.yaml***:
    - This is a .yaml file that shows what a sample submission to be given back to us should look like. 
    - Note: Only the brightest 10 stars are included
- ***score_submission.py***:
    - This script will test your submission file format, and will score the first 5 submissions. You can run this with:

        `python score_submission.py [your_submission_file] truth_5.fits`

# **Submission Format**

The format for target tracking submissions should include:
- **branch**: The branch that you intend to submit the submission to for scoring. ('target_track' or 'sidereal_track')
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

The following is an example of what the target tracking output should look like for the `sidereal_track` branch, before being submitted for scoring:
   

    branch: sidereal_track
    competitor_name: LLNL_dev
    display_true_name: true
    ---
    IOD:
    - rx: -1078188.235811465
      ry: -416352.80151879176
      rz: 7204171.363008768
      t: 1293621426.1067815
      vx: -3698.323018743952
      vy: 6369.638573203506
      vz: -193.1157223371356
    file: SAT_0001
    ---
    IOD:
    - rx: -2289445.197358595
      ry: -887395.1346812476
      rz: 6772431.399440748
      t: 1293637708.0508995
      vx: 4552.377873945053
      vy: 5429.810393862256
      vz: 2248.09634539989
    file: SAT_0002
    ---
    IOD:
      ...
    ...

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