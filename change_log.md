### **Note:** Please use data, and submission formatting, from the current version at all times. Previous versions will become obsolete, and no longer supported, after some time.

<br>

# **Change Log**
###
### **(Current version) v3.0 updates**
- Added the sidereal_track branch
- Added the target_track branch
- Generate ephemerides for satellites sampled randomly from public catalogs
- Update satellite size (beta distribution), albedo (uniform distribution), and brightness (Lambertian sphere) distributions
- Update image seeing, transparency, sky brightness, and pointing offsets to change smoothly as a function of time
### **v2.3 updates**
- Target branch now uses Gaia catalog
- Re-simulated both the Target and Sidereal branches
### **v2.2 updates**
- Added variable integration times to sidereal branch
- Added production of sky_flat.fits on the sidereal branch
- Fixed ra/dec scoring on sidereal branch
- Combined checking and scoring scripts in every branch
- Added variable sky brightness to sidereal branch
- Added flux calibration factor to sidereal branch image headers
- Re-simulated sidereal branch with all changes, increasing from 505 images to 1,010
- Add 10 truth objects to each applicable file (previously 5), and updated file names to reflect
### **v2.1 updates**
- Updated scoring scripts
- Replaced 4 scores in the calibrate/sidereal, detect/sidereal, and detect/target branches with "completeness" and "false positive rate"
- Column labels added to all score_history.txt files
- Added 5 more simulated images to calibrate branch
- Added a section to each branch's README file that explains what each score is
### **v2.0 updates**
- Added a "Calibrate - Sidereal Tracking" branch.
- Added optional "predicted downrange location" values to IOD branch submission files, giving the option to supply ra and dec values of the satellite at four times downrange (5, 15, 30, and 60 minutes).
- Added the option to supply covariances on the IOD branch, for both predicted downrange location, as well as predicted state vector.