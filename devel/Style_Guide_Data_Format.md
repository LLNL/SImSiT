
# **Data and Formatting**

This page describes the data, and the format of the data that we have prepackaged for download. All of our datasets can be downloaded here: *[Insert link]*. 

## **Detection Pipeline**

The detection pipeline consists of four branches: target tracking and sidereal tracking techniques, each with a "Detect" and a "Detect & Calibrate" branch. Each branch has 500 unique .fits images, each containing a satellite and stars. Detection algorithms will be tested on this data, and will then be scored on whether the detection was accurately predicted.

### **Target Tracking**

The target tracking detection dataset and files are located here: *[Insert link to zipped folder]* and includes the files for the "Detect" and the "Detect & Calibrate" branches.

After downloading the above zipped folder, you will have the following :

- 500 simulated target tracking .fits images (*combination of simulated and real still?*) of satellites with varying signal-to-noise ratios.
- ***truth5.fits***:
    - This is a .fits format file with 10 extensions inside, corresponding to the truth data from the first 5 .fits images.
    - There are 5 extensions titled "SAT_[XXXX]" which correspond to the satellite truth data for the first 5 images, each containing the folling 4 columns:
        - **x**: The x value of the satellite position (in pixel coordinates)
        - **y**: The y value of the satellite position (in pixel coordinates)
        - **flux**: The flux of the satellite
        - **snr**: The image signal-to-noise ratio
    - There are 5 extensions titled "STAR_[XXXX]" corresponding to the star truth data in each of the first 5 images, each containing the following 3 columns:
        - **x**: The x position of the center of the star streak (in pixel coordinates)
        - **y**: The y position of the center of the star streak (in pixel coordinates)
        - **flux**: The flux of the star
- ***sample_submission5.yaml***:
    - This is a .yaml file that shows what a sample submission back to us should look like. 
- ***check_submission.py***: 
    - A python script that you can run using `python check_submission.py [submission_file_to_check]` that will let you know if your submission is formatted properly to send back to us. 

#############Edit this section when the time comes##################
- ***[Insert Detect & Calbirate Branch Filename]***: 
    - The following files are for the target tracking "Detect & Calibrate" Branch:
        - ***[INSERT FILENAME]***
            - A single comma seperated variable (.csv) file that contains the truth information for the first 5 target tracking images, with 4 columns:
                - The image filename
                - satellite position (ra)
                - satellite position (dec)
                - Satellite magnitude 
            - ***[Insert Filename]*** 
                - A .yaml file that is a sample submission for the first 5 target tracking detect & calibrate branch images
            - ***[Insert Filename]***
                - A python script that verifies the submission contains the correct data and is the right format for the target tracking detect & calibrate branch
            - ***[Insert Filename]***
                - A python script that performs scoring on the first 5 target tracking detect & calibrate images
            - ***[Insert Filename]***
                - A python script showing how to assemble submission files for the target tracking detect & calibrate branch

### **Sidereal Tracking**

The sidereal tracking detection dataset and files are located here: *[Insert link to zipped folder]* and includes the files for the "Detect" and the "Detect & Calibrate" branches.

After downloading the above zipped folder, you will have data in the following format:

- 500 simulated sidereal tracking .fits images (*combination of simulated and real still?*) of satellites streaks with varying lengths, orientations, and signal-to-noise ratios.
- ***truth5.fits***:
    - This is a .fits format file with 10 extensions inside, corresponding to the truth data from the first 5 .fits images.
    - There are 5 extensions titled "SAT_[XXXX]" which correspond to the satellite truth data for the first 5 images, each containing the folling 5 columns:
        - **x0**: The x value of the first satellite streak endpoint (in pixel coordinates)
        - **y0**: The y value of the first satellite streak endpoint (in pixel coordinates)
        - **x1**: The x value of the second satellite streak endpoint (in pixel coordinates)
        - **y1**: The x value of the second satellite streak endpoint (in pixel coordinates)
        - **flux**: The flux of the satellite
    - There are 5 extensions titled "STAR_[XXXX]" corresponding to the star truth data in each of the first 5 images, each containing the following 3 columns:
        - **x**: The x position of the center of the star streak (in pixel coordinates)
        - **y**: The y position of the center of the star streak (in pixel coordinates)
        - **flux**: The flux of the star
- ***sample_submission5.yaml***:
    - This is a .yaml file that shows what a sample submission back to us should look like. 
- ***check_submission.py***: 
    - A python script that you can run using `python check_submission.py [submission_file_to_check]` that will let you know if your submission is formatted properly to send back to us. 

#############Edit this section when the time comes##################
- ***[Insert Detect & Calbirate Branch Filename]***: 
    - The following files are for the sidereal tracking "Detect & Calibrate" Branch:
        - ***[INSERT FILENAME]***
            - A single comma seperated variable (.csv) file that contains the truth information for the first 5 sidereal tracking images, with 4 columns:
                - The image filename
                - satellite position (ra)
                - satellite position (dec)
                - Satellite magnitude 
            - ***[Insert Filename]*** 
                - A .yaml file that is a sample submission for the first 5 sidereal tracking detect & calibrate branch images
            - ***[Insert Filename]***
                - A python script that verifies the submission contains the correct data and is the right format for the sidereal tracking detect & calibrate branch
            - ***[Insert Filename]***
                - A python script that performs scoring on the first 5 sidereal tracking detect & calibrate images
            - ***[Insert Filename]***
                - A python script showing how to assemble submission files for the sidereal tracking detect & calibrate branch

## **IOD Pipeline**

*[Explain more about this pipeline]* This section specifies the data and data format for IOD determination. This data is specifically for EO observations at this time. 

### **Target Tracking**

The target tracking IOD dataset is located here: *[Insert link]* 

This data consists of tracks of angles-only observations of varying lengths with noise properties consistent with target-tracking image metric measurements.

After downloading the above dataset you will have the following:

- ***[INSERT FOLDER TITLE]*** : 500 catalogs in the comma separated variable (.csv) format, one catalog for each satellite, with 
    - The catalogs have columns as follows:
        - ra
        - dec 
        - time *[units]*
        - sensor location *[units]*
        - *[Other potential columns]*
- ***[Insert Filename]***: 
    - Truth tables for the 5 satellites:
        - *[INFO]*
- ***[Insert Filename]*** 
    - Sample submission for the first 5 satellite IOD's

### **Sidereal Tracking**

The sidereal tracking IOD dataset is located here: *[Insert link]*

This data consists of tracks of angles-only observations of varying lengths with noise properties consistent with sidereal-tracking image metric measurements (streak images).

This dataset contains the same information as listed in the Target Tracking section above, using the following data:

Catalog : *[Insert foldername]*

Truth table : *[Insert Filename]*

Sample submission : *[Insert Filename]*


## **Maneuver Pipeline**
