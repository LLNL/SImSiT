# **Branch Overview**
This branch provides simulated images of ground-based sidereal tracking observations for 505 different satellites. This branch is designed to test the accuracy with which algorithms can measure the location of endpoints of the satellite being tracked in the image relative to withheld truth information. Truth information is provided for the first 5 images, therefore competitor scores are calculated on the remaining 500 images.

# **Data and Formatting**

This directory contains the following:
- ***detect_sidereal_images***:
    - A folder containing 505 simulated sidereal tracking FITS images of satellite streaks with varying lengths, orientations, and signal-to-noise ratios.
- ***truth_5_detect_sidereal.fits***:
    - This is a FITS format file with 10 extensions inside, corresponding to the truth data from the first 5 FITS images.
    - There are 5 extensions titled "SAT_[XXXX]" which correspond to the satellite truth data for the first 5 images, each containing the following 5 columns:
        - **x0**: The x position of one satellite streak endpoint (in pixel coordinates)
        - **y0**: The y position of one satellite streak endpoint (in pixel coordinates)
        - **x1**: The x position of the other satellite streak endpoint (in pixel coordinates)
        - **y1**: The y position of the other satellite streak endpoint (in pixel coordinates)
        - **flux**: The flux of the satellite (in analog to digital units (ADU))
    - There are 5 extensions titled "STAR_[XXXX]" corresponding to the star truth data in each of the first 5 images, each containing the following 3 columns:
        - **x**: The x position of the center of the star (in pixel coordinates)
        - **y**: The y position of the center of the star (in pixel coordinates)
        - **flux**: The flux of the star (in analog to digital units (ADU))
- ***sample_submission_5_detect_sidereal.yaml***:
    - This is a .yaml file that shows what a sample submission to be given back to us should look like. 
- ***check_submission_detect_sidereal.py***: 
    - A python script that you can run using 

        `python check_submission_detect_sidereal.py [submission_file_to_check]` 
        
        that will let you know if your submission is formatted properly for scoring. 
- ***sample_score_submission_detect_sidereal.py***:
    - A python script that you can run using

        `python sample_score_submission_detect_sidereal.py [submission_file] [truth_file]`

    - You can test the scoring with the files we provided you:
    ***sample_submission_5_detect_sidereal.yaml*** and ***truth_5_detect_sidereal.fits***
    - Then you can try scoring the first 5 of your own submission file:
     ***[your_submission_file].yaml*** and ***truth_5_detect_sidereal.fits***

# **Submission Format**

A .yaml file should be submitted with the following information at the top of the file:
- **branch**: The branch that you intend to submit the submission to for scoring.
- **competitor_name**: The name you would like us to use to identify your submission file. This will be the name displayed on the dashboard.
- **display_true_name**: A boolean deciding whether you would like your `competitor_name` publically displayed on the dashboard. If false, your score will be assigned an anonymous `competitor_name`.

The following information should be included afterwards, for each image:
- **file**: The current image filename
- **sats**: The information for the satellite in the image:
  - **flux**: The satellite flux (in analog to digital units (ADU))
  - **x0**: The x position of one satellite streak endpoint (in pixel coordinates)
  - **x1**: The x position of the other satellite streak endpoint (in pixel coordinates)
  - **y0**: The y position of one satellite streak endpoint (in pixel coordinates)
  - **y1**: The y position of the other satellite streak endpoint (in pixel coordinates)
- **stars (optional)**: The information for the star streaks in the image:
  - **flux**: The flux of the star (in analog to digital units (ADU))
  - **x**: The x position of the center of the star (in pixel coordinates)
  - **y**: The y position of the center of the star (in pixel coordinates)

You can use `check_submission_detect_sidereal.py` to format your detection branch sidereal tracking submissions properly.

The submission file should have the following format (stars optional):

    branch: Detect Sidereal
    competitor_name: Competitor A
    display_true_name: true
    file: 000.fits
    sats:
    - flux: 58625.247235005365
      x0: 745.0992620337679
      x1: 651.7329351747386
      y0: 673.184121283623
      y1: 680.6109279569921
    stars:
    - flux: 97251.24373052332
      x: 900.1646259076072
      y: 957.4608780491294
    - flux: 107305.30596431502
      x: 439.8455792717815
      y: 340.2567502148446
    - flux: 100140.76596160399
      x: 968.8341348194286
      y: 101.77436929642184
      ...
    ---
    file: 001.fits
      ...
    ...
    

# **Scores**
There will be a summary score reflected on the main page of the dashboard for this branch. The score that we use to generate that summary score is `pos_RMSE` (This is the Root Mean Square Error of the submitted position).

**The following are the entires in `score_history_detect_sidereal.txt`, by index:**
- **[comp_name]**: Competitor name
- **[sub_filename]**: Submission filename
- **[date]**: Date of submission
- **[completeness]**: The number of true satellites found in the submission divided by the number of true satellites
- **[fpr]**: False positive rate: The number of false satellites found in the submission (that could not be associated with a true satellite) divided by the number of total satellites found in the submission
- **[pos_RMSE]**: Root Mean Square Error (RMSE) of the submitted position
- **[mag_RMSE]**: Root Mean Square Error (RMSE) of the magnitude, determined from the submitted flux