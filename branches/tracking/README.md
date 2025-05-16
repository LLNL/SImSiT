# **Branch Overview**
This branch provides simulated images of ground-based observations for 1,010 different satellites. This branch is designed to test the accuracy with which algorithms can measure the location of endpoints of the satellite being tracked in the image relative to withheld truth information. Truth information is provided for the first 10 images, therefore competitor scores are calculated on the remaining 1,000 images.

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
    - A folder containing 1,010 simulated FITS images of satellites with varying signal-to-noise ratios. Files are named [XXXX].fits, where the [XXXX] is the satellite number. (ex: 0235.fits is satellite number 235)
- ***sky_flat.fits***:
  - A flat field image computed by normalizing each individual challenge image, taking the median over the course of the night, and dividing out by the image mean.
- ***truth_10.fits***:
    - This is a FITS format file with 10 extensions inside, corresponding to the truth data from the first 10 FITS images.
    - There are 10 extensions titled "SAT_[XXXX]" which correspond to the satellite truth data for the first 10 images, each containing the following:
        - **x0**: The x position of one satellite streak endpoint (in pixel coordinates)
        - **y0**: The y position of one satellite streak endpoint (in pixel coordinates)
        - **x1**: The x position of the other satellite streak endpoint (in pixel coordinates)
        - **y1**: The y position of the other satellite streak endpoint (in pixel coordinates)
        - **ra0**: The right ascension of one satellite streak endpoint (in degrees)
        - **dec0**: The declination of one satellite streak endpoint (in degrees)
        - **ra1**: The right ascension of the other satellite streak endpoint (in degrees)
        - **dec1**: The declination of the other satellite streak endpoint (in degrees)
        - **mag**: The magnitude of the satellite in the i-band (Equivelant to the Gaia G-band for these simulations)
        - **nphot**: The flux of the satellite (in units of election counts)
    - There are 10 extensions titled "STAR_[XXXX]" corresponding to the star truth data in each of the first 10 images, each containing the following columns:
        - **ra**: The right ascension of the center of the star (in degrees)
        - **dec**: The declination of the center of the star (in degrees)
        - **x**: The x position of the center of the star (in pixel coordinates)
        - **y**: The y position of the center of the star (in pixel coordinates)
        - **i_mag**: The magnitude of the star in the i-band (Equivelant to the Gaia G-band for these simulations)
        - **nphot**: The flux of the star (in units of election counts)
- ***sample_submission_10.yaml***:
    - This is a .yaml file that shows what a sample submission to be given back to us should look like. 
    - Note: Only the brightest 10 stars are included
- ***score_submission.py***:
    - This script will test your submission file format, and will score the first 10 submissions. You can run this with:

        `python score_submission.py [your_submission_file] truth_10.fits`

# **Submission Format**

A .yaml file should be submitted with the following information at the top of the file:
- **branch**: The branch that you intend to submit the submission to for scoring. ('target' or 'sidereal')
- **competitor_name**: The name you would like us to use to identify your submission file. This will be the name displayed on the dashboard.
- **display_true_name**: A boolean deciding whether you would like your `competitor_name` publically displayed on the dashboard. If false, your score will be assigned an anonymous `competitor_name`.

The following information should be included afterwards, for each image:
- **file**: The current image filename
- **sats**: The information for the satellite in the image. "Detection" scores are given based on x, y, and flux submissions, while "detection" scores are given based on ra, dec, and mag quantities. At MINIMUM one set of submissions must be supplied.
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
- **stars (OPTIONAL)**: The information for the star streaks in the image:
  - **dec**: The declination of the center of the star (in degrees)
  - **flux**: The flux of the star (in analog to digital units (ADU))
  - **mag**: The magnitude of the star
  - **ra**: The right ascension of the center of the star (in degrees)
  - **x**: The x position of the center of the star (in pixel coordinates)
  - **y**: The y position of the center of the star (in pixel coordinates)

See ***sample_submission_10.yaml*** for an example submission file. 

# **Scores**
There will be a summary score reflected on the main page of the dashboard for this branch. The score that we use to generate that summary score is `pos_RMSE` (This is the Root Mean Square Error of the submitted position).

**The following are the entires in `score_history.txt`, by index:**
- **[comp_name]**: Competitor name
- **[branch]**: Branch you wish to submit to
- **[sub_filename]**: Submission filename
- **[date]**: Date of submission
- **[detect_completeness]**: The number of true satellites found in the submission divided by the number of true satellites. (Using x0, y0, x1, y0, and flux)
- **[calibrate_completeness]**: The number of true satellites found in the submission divided by the number of true satellites. (Using ra0, dec0, ra1, dec0, and mag)
- **[detect_fpr]**: False positive rate: The number of false satellites found in the submission (that could not be associated with a true satellite) divided by the number of total satellites found in the submission. (Using x0, y0, x1, y0, and flux)
- **[calibrate_fpr]**: False positive rate: The number of false satellites found in the submission (that could not be associated with a true satellite) divided by the number of total satellites found in the submission. (Using ra0, dec0, ra1, dec0, and mag)
- **[xy_RMSE]**: Root Mean Square Error (RMSE) of the submitted xy positions
- **[radec_RMSE]**: Root Mean Square Error (RMSE) of the submitted ra and dec positions
- **[flux_RMSE]**: Root Mean Square Error (RMSE) of the submitted flux
- **[mag_RMSE]**: Root Mean Square Error (RMSE) of the submitted magnitude