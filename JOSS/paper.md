---
title: 'SatIST: The Satellite Image Simulation Toolkit'
tags:
  - Python
  - astronomy
  - satellites
  - space domain awareness
authors:
  - name: LLNL SatIST team 
    equal-contrib: true
    affiliation: 1

affiliations:
 - name: Space Science Institute, Lawrence Livermore National Laboratory, 7000 East Ave., Livermore, CA 94550, USA
   index: 1

date: 17 October 2024
bibliography: paper.bib

aas-doi:
aas-journal:
---

# Summary

The Satellite Image Simulation Toolkit (`SatIST`) is a Python software package designed 
to generate diverse and realistic satellite imaging scenarios. `SatIST` is built on top
of `GalSim` [@Rowe2015] and `SSAPy` [@SSAPy2023, @Yeager20203] and serves as a toolkit 
for simulating imaging data that supports the development and testing of algorithms used in 
satellite detection, calibration, and characterization. SatIST provides a suite of 
simulation tools that allow users to replicate various satellite observation conditions, 
including sidereal (where the sensor tracks fixed stars, causing the satellite to appear 
as a streak in the image) target tracking (where the sensor follows the satellite 
during exposure). By enabling the creation of scenarios that mimic real-world satellite 
operations, SatIST facilitates advancements in satellite image data processing and the study 
of satellite behavior under different observational parameters. `SatIST` has been used to 
benchmark and test machine learning alogrithms that aim to identify closely
separated satilites in astronomical images [@Pruett2023].


![Example of simulated image output from `SatIST`. Left: An example image taken with sidereal (star) tracking, 
where the fast-moving satellite appears as a streak across the image. Right: An example image taken with 
satellite tracking, where the stars appear as streaks, and the satellite appears as a point source.](example_satist_out.png)

# Statement of need


# Method


`SatIST` combines the functionality of `GalSim` and `SSAPy`. `GalSim`
enables users to customize a wide range of observing parameters, 
such as detector settings, optical system properties, and realistic atmospheric conditions. 
Specifcally users can manually configure sensor parameters such as vignetting and optical distortion 
based on real, characterized sensors or simulateand test future sensor designs.

`SatIST` generates realistic astronmical background scenes using star positions, fluxes, and colors
from the Gaia Data Release 2 catalog [@GaiaCollaboration2018]. Gaia optical magnitudes can be transformed
into other passbands such as Short Wave Infrared (SWIR) to generate releastics astronomical scenes in different
wavelengths.  

Through its integration, the suite facilitates the generation of simulated images for any pairing of ground-based sensors and satellite systems.




# Acknowledgements

`SatIST` depends on numpy [@Harris2020], scipy [@Virtanen2020], asdf [@Greenfield2015], matplotlib [@Hunter2007], and scikit-learn [@sklearn_api].
This work was performed under the auspices of the U.S.
Department of Energy by Lawrence Livermore National
Laboratory (LLNL) under Contract DE-AC52-07NA27344. The document number is LLNL-JRNL-XXXX and the code number is LLNL-CODE-XXXX.



# References
