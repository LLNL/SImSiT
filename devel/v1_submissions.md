
# **Submission Format**

This page describes the format that your files should be in before being returned to us for scoring.  

## **Detection Pipeline**

### **Target Tracking**

A .yaml file should be submitted with the following information for each image:
- **file**: The current image filename
- **sats**: The information for the satellite in the image
  - **flux**: The satellite flux
  - **x**: The x position of the satellite (in pixel coordinates)
  - **y**: The y position of the satellite (in pixel coordinates)
- **(Optional) stars**: The information for the star streaks in the image
  - **flux**: The flux of the star
  - **x0**: The x position of the first endpoint of the star streak (in pixel coordinates)
  - **x1**: The x position of the second endpoint of the star streak (in pixel coordinates)
  - **y0**: The y position of the first endpoint of the star streak (in pixel coordinates)
  - **y1**: The y position of the second endpoint of the star streak (in pixel coordinates)

You can use the *[insert link (check_submission_detect_target.py)]* file to format your detection branch target tracking submissions properly. 

The submission file should have the following format (stars optional):

    file: 000.fits
    sats:
    - flux: 64786.66258436015
      x: 508.9425769817726
      y: 489.91057109837544
    stars: 
    - flux: 97733.41294322626
      x0: 561.6798168573041
      y0: 834.6402846654115
      x1: 439.8455792717815
      y1: 340.2567502148446
    - flux: 103997.20692109835
      x0: 836.5198689474443
      y0: 450.69324341254344
      x1: 968.8341348194286
      y1: 101.77436929642184
    ---
    file: 001.fits
      ...
    ...

### **Sidereal Tracking**

A .yaml file should be submitted with the following information for each image:
- **file**: The current image filename
- **sats**: The information for the satellite in the image:
  - **flux**: The satellite flux
  - **x0**: The x position of the first streak endpoint of the satellite (in pixel coordinates)
  - **y0**: The y position of the first streak endpoint of the satellite (in pixel coordinates)
  - **x1**: The x position of the second streak endpoint of the satellite (in pixel coordinates)
  - **y1**: The y position of the second streak endpoint of the satellite (in pixel coordinates)
- **(Optional) stars**: The information for the star streaks in the image:
  - **flux**: The flux of the star
  - **x**: The x position of the star (in pixel coordinates)
  - **y**: The x position of the star (in pixel coordinates)

You can use the *[insert link (check_submission_detect_sidereal.py)]* file to format your detection branch sidereal tracking submissions properly.

The submission file should have the following format (stars optional):

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
    ---
    file: 001.fits
      ...
    ...

## **IOD Pipeline**

### **Target Tracking**

The format for target tracking IOD submissions should include:
- the satellite state (r, v) composed each of 3 component tuples (x, y, and z components) in Geocentric Celestial Reference System (GCRS) [*cgs units*].
- epoch for each satellite in units of GPS seconds. This should match the time of the last observation in the corresponding satellite's catalog.

You can use the *[insert link (check_submission_iod_target.py)]* file to format your target tracking IOD submissions properly.

The following is an example of what the IOD output should look like, before being submitted for scoring:

    IOD:
    - rx: 
      ry: 
      rz: 
      t:
      vx: 
      vy: 
      vz: 
    sat: SAT_0000
    ---
    IOD:
      ...
    ...

Where the number after "SAT" corresponds to the same number of the catalog that was used to determine the IOD.

