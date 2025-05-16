
# **Data and Formatting**

This page describes the data, and the format of the data that we have prepackaged for download. All of our datasets can be downloaded here: *[Insert link to top level of directory (xfiles)]*. 

## **Detection Pipeline**

The detection pipeline consists of two branches: target tracking and sidereal tracking. Each branch has 500 unique .fits images, each containing a satellite and stars. Detection algorithms will be tested on this data, and will then be scored on whether the detection was accurately predicted.

### **Target Tracking**

The target tracking detection dataset and files are located here: *[Insert link to target tracking folder (xfiles/detection/target/)]*.

After downloading the above zipped folder, you will have the following :

- 500 simulated target tracking .fits images of satellites with varying signal-to-noise ratios.
- ***truth_5_detect_target.fits***:
    - This is a .fits format file with 10 extensions inside, corresponding to the truth data from the first 5 .fits images.
    - There are 5 extensions titled "SAT_[XXXX]" which correspond to the satellite truth data for the first 5 images, each containing the following:
        - **x**: The x value of the satellite position (in pixel coordinates)
        - **y**: The y value of the satellite position (in pixel coordinates)
        - **flux**: The flux of the satellite
    - There are 5 extensions titled "STAR_[XXXX]" corresponding to the star truth data in each of the first 5 images, each containing the following 3 columns:
        - **x0**: The x position of the first endpoint of the star streak (in pixel coordinates)
        - **x1**: The x position of the second endpoint of the star streak (in pixel coordinates)
        - **y0**: The y position of the first endpoint of the star streak (in pixel coordinates)
        - **y1**: The y position of the second endpoint of the star streak (in pixel coordinates)
        - **flux**: The flux of the star
- ***sample_submission_5_detect_target.yaml***:
    - This is a .yaml file that shows what a sample submission back to us should look like. 
- ***check_submission_detect_target.py***: 
    - A python script that you can run using:

        `python check_submission_detect_target.py [submission_file_to_check]` 
     
        that will let you know if your submission is formatted properly for scoring.

### **Sidereal Tracking**

The sidereal tracking detection dataset and files are located here: *[Insert link to target tracking folder (xfiles/detection/sidereal/)]*.

After downloading the above zipped folder, you will have data in the following format:

- 500 simulated sidereal tracking .fits images of satellites streaks with varying lengths, orientations, and signal-to-noise ratios.
- ***truth_5_detect_sidereal.fits***:
    - This is a .fits format file with 10 extensions inside, corresponding to the truth data from the first 5 .fits images.
    - There are 5 extensions titled "SAT_[XXXX]" which correspond to the satellite truth data for the first 5 images, each containing the folling 5 columns:
        - **x0**: The x position of the first satellite streak endpoint (in pixel coordinates)
        - **x1**: The x position of the second satellite streak endpoint (in pixel coordinates
        - **y0**: The y position of the first satellite streak endpoint (in pixel coordinates)
        - **y1**: The x position of the second satellite streak endpoint (in pixel coordinates)
        - **flux**: The flux of the satellite
    - There are 5 extensions titled "STAR_[XXXX]" corresponding to the star truth data in each of the first 5 images, each containing the following 3 columns:
        - **x**: The x position of the center of the star streak (in pixel coordinates)
        - **y**: The y position of the center of the star streak (in pixel coordinates)
        - **flux**: The flux of the star
- ***sample_submission_5_detect_sidereal.yaml***:
    - This is a .yaml file that shows what a sample submission back to us should look like. 
- ***check_submission_detect_sidereal.py***: 
    - A python script that you can run using 
    
        `python check_submission_detect_sidereal.py [submission_file_to_check]` 
        
        that will let you know if your submission is formatted properly for scoring. 

## **IOD Pipeline**

*[Explain more about this pipeline]* This section specifies the data and data format for IOD determination. This data is specifically for angles-only EO observations at this time. 

### **Target Tracking**

The target tracking IOD dataset is located here: *[Insert link to target tracking folder (xfiles/iod/target/)]* 

This data consists of tracks of angles-only observations of varying lengths with noise properties consistent with target-tracking image metric measurements.

After downloading the above dataset you will have the following:

- 500 catalogs in comma seperated variable (.csv) format, one catalog for each satellite (named 0000 through 0499), with the following columns:
    - **t**: The time of that observation in GPS seconds. (The last observation time is the time that you want to calculate the IOD at). 
    - **ra**: The right ascension of the satellite during the given observation in GCRS coordinates in degrees
    - **dec**: The declination of the satellite during the given observation in GCRS coordinates in degrees
    - **px**: The x position of the observer at that observation *[km]*
    - **py**: The y position of the observer at that observation *[km]*
    - **pz**: The z position of the observer at that observation *[km]*
- ***truth_5_iod_target.fits***: 
    - Truth tables for the 5 satellites, each containing the satellites' IOD at the time of the last observation in the data file for each satellite:
        - **rx**: x component of the satellite location in *[units]*
        - **ry**: y component of the satellite location in *[units]*
        - **rz**: z component of the satellite location in *[units]*
        - **vx**: x component of the satellite velocity in *[units]*
        - **vy**: y component of the satellite velocity in *[units]*
        - **vz**: z component of the satellits velocity in *[units]*
        - **t**: Time of the epoch in GPS seconds. (This should match the time of the last observation in the data file for that satellite.)
- ***sample_submission_5_iod_target.fits*** 
    - Sample submission for the first 5 satellite IOD's
- ***check_submission_iod_target.py***: 
    - A python script that you can run using 
    
        `python check_submission_iod_target.py [submission_file_to_check]` 
        
        that will let you know if your submission is formatted properly for scoring. 

