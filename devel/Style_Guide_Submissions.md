
# **Submission Format**

This page describes the format that your files should be in before being uploaded to us for scoring.  

## **Detection Pipeline**

### **1. Target Tracking - Detect Branch**

A .yaml file should be submitted with the following information for each image:
- **file**: The current image filename
- **sats**: The information for the satellite in the image:
  - **flux**: The satellite flux
  - **x**: The x position of the satellite (in pixel coordinates)
  - **y**: The y position of the satellite (in pixel coordinates)
- **stars**: The information for the star streaks in the image:
  - **flux**: The flux of the star
  - **x**: The x position of the center of the star streak (in pixel coordinates)
  - **y**: The x position of the center of the star streak (in pixel coordinates)

You can use the *[insert file name and link]* file to format your detection branch target tracking submissions properly. 

#############Edit this, now we are only using the x, y position of the center of the star streak?

The submission file should have the following format:

    file : 000.fits
    sats :
    - flux : 64786.66258436015
      x : 508.9425769817726
      y : 489.91057109837544
    stars : 
    - flux : 97733.41294322626
      x0 : 561.6798168573041
      x1 : 175.93168378833437
      y0 : 834.6402846654115
      y1 : 565.6853272128955
    - flux : 103997.20692109835
      x0 : 836.5198689474443
      x1 : 135.949891056725
      y0 : 450.69324341254344
      y1 : 908.4001418248758

    file : 001.fits
      ...
    ...

### **2. Sidereal Tracking - Detect Branch**

A .yaml file should be submitted with the following information for each image:
- **file**: The current image filename
- **sats**: The information for the satellite in the image:
  - **flux**: The satellite flux
  - **x0**: The x position of the first streak endpoint of the satellite (in pixel coordinates)
  - **y0**: The y position of the first streak endpoint of the satellite (in pixel coordinates)
  - **x1**: The x position of the second streak endpoint of the satellite (in pixel coordinates)
  - **y1**: The y position of the second streak endpoint of the satellite (in pixel coordinates)
- **stars**: The information for the star streaks in the image:
  - **flux**: The flux of the star
  - **x**: The x position of the star (in pixel coordinates)
  - **y**: The x position of the star (in pixel coordinates)

You can use the *[insert file name and link]* file to format your detection branch sidereal tracking submissions properly.

The submission file should have the following format:

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

    file : 001.fits
      ...
    ...

### **3. Target Tracking - Detect & Calibrate Branch**

A .yaml file should be submitted with the following information for each image:
- Satellite location (position in ra and dec)
- Satellite magnitude
- Location of each star in the image (The ra and dec position of each streak endpoint) 
- magnitude of each star in the image

You can use the *[insert file name and link]* file to format your detect & calibrate branch target tracking submissions properly. 

The submission file should have the following format:

    file : 000.fits
    sats :
    - mag : 
      x : 
      y : 
    stars : 
    - mag : 
      x0 : 
      x1 : 
      t0 : 
      y1 : 
    - mag : 
      x0 : 
      x1 : 
      y0 : 
      y1 : 

    file : 001.fits
      ...
    ...

### **4. Sidereal Tracking - Detect & Calibrate Branch**

A .yaml file should be submitted with the following information for each image:
- Satellite location (The ra and dec position of each streak endpoint)
- Satellite magnitude
- Location of each star in the image (ra and dec position) 
- Magnitude value for each star in the image 

You can use the *[insert file name and link]* file to format your detect & calibrate branch sidereal tracking submissions properly.

The submission file should have the following format:

    file: 000.fits
    sats:
    - mag: 
      x0: 
      x1: 
      y0: 
      y1: 
    stars:
    - mag: 
      x: 
      y: 
    - mag: 
      x: 
      y: 
    - mag: 
      x: 
      y: 

    file : 001.fits
      ...
    ...


## **IOD Pipeline**

The format for target tracking IOD submissions should include:
- the satellite state (r, v) composed each of 3 component tuples (x, y, and z components) in Geocentric Celestial Reference System (GCRS) units of m and m/s respectively.
- epoch for each satellite in units of GPS seconds

NOTE: The only difference between the target tracking, sidereal tracking, and fractional tracking IOD methods is the number of observations of the input data.

*[Get more information from Nate. Make sure the above statement is still correct]*

You can use the *[insert file name and link]* file to format your target tracking IOD submissions properly, the *[insert]* file to format your sidereal tracking IOD submissions properly, and the *[insert]* file to format your fractional tracking IOD submissions properly. *[They may all be the same and can use one file?]*

The following is an example of what the IOD output should look like, before being submitted for scoring:

    Sat001 :
    - rx : 
      ry : 
      rz : 
      vx : 
      vy : 
      vz : 
      t : 
    
    Sat002:
      ...
    ...

Where the number after "Sat" corresponds to the same number of the catalog that was used to determine the IOD.

## **Maneuver Pipeline**

