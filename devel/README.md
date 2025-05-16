># **Overview**
This package contains a variety of challenge data sets designed to test the robustness and accuracy of algorithms relevant to satellite detection, calibration, and characterization. Truth information is provided for the first 5-10 images/catalogs in each branch. The branches currently available are the following:

1. **Target Tracking**
    - Includes simulated images of ground-based target tracking observations for 1,010 different targets
    - This branch tests on "detection" tasks using x,y locations and flux.
    - The README is located at `target/README_target.md`

2. **Sidereal Tracking**
    - Simulated images of ground-based sidereal tracking observations for 1,010 different targets
    - Includes "detection" and optional "calibrate" information.
    - The README is located at `sidereal/README_sidereal.md`

3. **IOD**
    - Tracks of optical, angle-only observations of varying lengths for 445 targets
    - This branch tests on 
    - The README is located at `iod/README_iod.md`

# **Information**

- A dashboard with competitor scores for each branch is located at https://xfiles.llnl.gov.
    - To request an account please email xfiles@llnl.gov.
- For a detailed list of what changed between data versions, please consult `change_log.md`.

