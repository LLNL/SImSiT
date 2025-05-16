# **Branch Overview**
This branch provides simulated images of ground-based sidereal tracking observations for 505 different satellites. This branch is designed to test the accuracy with which algorithms can measure the location of endpoints of the satellite being tracked in the image relative to withheld truth information. Truth information is provided for the first 5 images, therefore competitor scores are calculated on the remaining 500 images.

- Images include optical distortion, a constant surface brightness background, and vignetting.
- Star positions and fluxes are derived from Gaia DR2. 
  - Simulated star positions account for proper motion, parallax, and annual aberration. 
  - Simulated fluxes are derived from Gaia G-band measurements. (For simulation purposes, we use flux_i == flux_G)
- Pixel coordinate convention follows FITS:  the center of the lower left pixel is (1.0, 1.0)
- Important header keywords:
  - DATE-BEG, DATE-END - time of shutter open/close.  shutter travel assumed to be instantaneous
  - OBSGEO-X, OBSGEO-Y, OBSGEO-Z - ITRS position of observatory
  - Approximate TAN wcs keywords: (CTYPE, CRPIX, CD, CUNIT, CRVAL)

# **Data and Formatting**

This directory contains the following:
- ***calibrate_sidereal_images***:
    - A folder containing 505 simulated sidereal tracking FITS images of satellite streaks with varying lengths, orientations, and signal-to-noise ratios.
- ***truth_5_calibrate_sidereal.fits***:
    - This is a FITS format file with 10 extensions inside, corresponding to the truth data from the first 5 FITS images.
    - There are 5 extensions titled "SAT_[XXX]" which correspond to the satellite truth data for the first 5 images, each containing the following 5 columns:
        - **x0**: The x position of one satellite streak endpoint (in pixel coordinates)
        - **y0**: The y position of one satellite streak endpoint (in pixel coordinates)
        - **x1**: The x position of the other satellite streak endpoint (in pixel coordinates)
        - **y1**: The y position of the other satellite streak endpoint (in pixel coordinates)
        - **ra0**: The right ascension of one satellite streak endpoint (in radians)
        - **dec0**: The declination of one satellite streak endpoint (in radians)
        - **ra1**: The right ascension of the other satellite streak endpoint (in radians)
        - **dec1**: The declination of the other satellite streak endpoint (in radians)
        - **mag**: The magnitude of the satellite in the i-band (Equivelant to the Gaia G-band for these simulations)
        - **nphot**: The flux of the satellite (in units of election counts)
    - There are 5 extensions titled "STAR_[XXX]" corresponding to the star truth data in each of the first 5 images, each containing the following 3 columns:
        - **ra**: The right ascension of the center of the star (in radians)
        - **dec**: The declination of the center of the star (in radians)
        - **x**: The x position of the center of the star (in pixel coordinates)
        - **y**: The y position of the center of the star (in pixel coordinates)
        - **i_mag**: The magnitude of the star in the i-band (Equivelant to the Gaia G-band for these simulations)
        - **nphot**: The flux of the star (in units of election counts)
- ***sample_submission_5_calibrate_sidereal.yaml***:
    - This is a .yaml file that shows what a sample submission to be given back to us should look like. 
- ***check_submission_calibrate_sidereal.py***: 
    - A python script that you can run using 

        `python check_submission_calibrate_sidereal.py [submission_file_to_check]` 
        
        that will let you know if your submission is formatted properly for scoring. 
- ***sample_score_submission_calibrate_sidereal.py***:
    - A python script that you can run using

        `python sample_score_submission_calibrate_sidereal.py [submission_file] [truth_file]`

    - You can test the scoring with the files we provided you:
    ***sample_submission_5_calibrate_sidereal.yaml*** and ***truth_5_calibrate_sidereal.fits***
    - Then you can try scoring the first 5 of your own submission file:
     ***[your_submission_file].yaml*** and ***truth_5_calibrate_sidereal.fits***

# **Submission Format**

A .yaml file should be submitted with the following information at the top of the file:
- **branch**: The branch that you intend to submit the submission to for scoring.
- **competitor_name**: The name you would like us to use to identify your submission file. This will be the name displayed on the dashboard.
- **display_true_name**: A boolean deciding whether you would like your `competitor_name` publically displayed on the dashboard. If false, your score will be assigned an anonymous `competitor_name`.

The following information should be included afterwards, for each image:
- **file**: The current image filename
- **sats**: The information for the satellite in the image:
  - **dec0**: The declination of one satellite streak endpoint (in degrees)
  - **dec1**: The declination of the other satellite streak endpoint (in degrees)
  - **flux**: The satellite flux (in analog to digital units (ADU))
  - **mag**: The magnitude of the satellite
  - **ra0**: The right ascension of one satellite streak endpoint (in degrees)
  - **ra1**: The right ascension of the other satellite streak endpoint (in degrees)
  - **x0**: The x position of one satellite streak endpoint (in pixel coordinates)
  - **x1**: The x position of the other satellite streak endpoint (in pixel coordinates)
  - **y0**: The y position of one satellite streak endpoint (in pixel coordinates)
  - **y1**: The y position of the other satellite streak endpoint (in pixel coordinates)
- **stars (optional)**: The information for the star streaks in the image:
  - **dec**: The declination of the center of the star (in degrees)
  - **flux**: The flux of the star (in analog to digital units (ADU))
  - **mag**: The magnitude of the star
  - **ra**: The right ascension of the center of the star (in degrees)
  - **x**: The x position of the center of the star (in pixel coordinates)
  - **y**: The y position of the center of the star (in pixel coordinates)

You can use `check_submission_calibrate_sidereal.py` to format your calibrate branch sidereal tracking submissions properly.

The submission file should have the following format (stars optional):

    branch: Calibrate Sidereal
    competitor_name: Competitor A
    display_true_name: true
    file: 000.fits
    sats:
    - dec0: 0.9824085674177758
      dec1: 0.9827694532849132
      flux: 7457.090027770151
      mag: 13.494536657817662
      ra0: 0.7990369541568129
      ra1: 0.7961051351350783
      x0: 915.2029741109612
      x1: 1089.8194690755304
      y0: 724.7517276023891
      y1: 582.6143955011476
    stars:
    - dec: 0.984135803379026
      flux: 1189804.0268676954
      mag: 8.050434365599326
      ra: 0.7917760447120902
      x: 1436.2030497375654
      y: 441.1471290131511
    - dec: 0.9828443340978492
      flux: 777127.4472634476
      mag: 8.567482602510315
      ra: 0.787298190854206
      x: 1507.5133555757925
      y: 67.04054447333851
    - dec: 0.9814851323107805
      flux: 245450.66978989646
      mag: 9.755139911280292
      ra: 0.789036863363843
      x: 1281.0602737235993
      y: 54.51005867202332
      ...
    ---
    file: 001.fits
      ...
    ...
    

# **Scores**
There will be a summary score reflected on the main page of the dashboard for this branch. The score that we use to generate that summary score is `pos_RMSE` (This is the Root Mean Square Error of the submitted position).

**The following are the score entires in `score_history_calibrate_sidereal.txt`, by index:**
- **[comp_name]**: Competitor name
- **[sub_filename]**: Submission filename
- **[date]**: Date of submission
- **[completeness]**: The number of true satellites found in the submission divided by the number of true satellites
- **[fpr]**: False positive rate: The number of false satellites found in the submission (that could not be associated with a true satellite) divided by the number of total satellites found in the submission
- **[pos_RMSE]**: Root Mean Square Error (RMSE) of the submitted position
- **[mag_RMSE]**: Root Mean Square Error (RMSE) of the submitted magnitude