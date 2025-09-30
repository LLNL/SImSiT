==========
Background
==========

Summary
-------

The Satellite Image Simulation Toolkit (SImSiT) is a Python software package designed 
to generate diverse and realistic satellite imaging scenarios. SImSiT is built on top 
of GalSim (:cite:t:`Rowe2015`) and SSAPy (:cite:t:`SSAPy2023` :cite:t:`@Yeager2023`) and serves as a toolkit for 
simulating imaging data that supports the development and testing of algorithms for 
satellite detection, calibration, and characterization. SImSiT provides a suite of 
simulation tools that allow users to replicate various satellite observation conditions, 
including: Sidereal tracking (where the sensor tracks fixed stars, causing the satellite 
to appear as a streak in the image), and Target tracking (where the sensor follows the 
satellite during exposure), as shown in Figure 1. By enabling the creation of scenarios 
that mimic real-world satellite operations, SImSiT facilitates advancements in 
satellite image data processing and the study of satellite behavior under different
observational parameters. SImSiT has been used to benchmark and test machine 
learning algorithms that aim to identify closely separated satellites in
ground-based astronomical images (:cite:t:`Pruett2023`). 

Statement of need
-----------------

The astronomy community has built many specialized image simulation tools. These include tools specific 
to a given mission or telescope, such as `Roman-I-Sim` <https://romanisim.readthedocs.io/en/latest/>_`, 
which generates high-fidelity images for the
upcoming Nancy Grace Roman Space Telescope :cite:t:`Spergel2015`, as well as end-to-end image simulations 
of the ground-based Vera C. Rubin Observatory's LSSTCam (e.g.,:cite:t:`Peng2013`).
The community has also built more general software tools that can simulate detailed telescope designs 
using ray-tracing methods—PhoSim (:cite:t:`Peterson2015`) —and generate general astronomical scenes, for example, 
GalSim (:cite:t:`@Rowe2015`) and SkyMaker (:cite:t:`Bertin2009`). Moreover, software packages that model detailed observational 
characteristics of satellites have also been developed (:cite:t:`Fankhauser2023`). However, no software packages that 
combine general image simulations of satellites with realistic orbits currently exist—SImSiT fills this gap.

With the advent of large astronomical imaging surveys such as LSST :cite:t:`Ivezi2019` and The Zwicky Transient 
Facility (:cite:t:`Bellm2019`) in the coming years, and the proliferation of satellites visible in those surveys (:cite:t:`Mroz2022`), 
simulating images with satellites has become increasingly important for the broad astronomy and space communities. 

Method
------

SImSiT combines the functionality of GalSim and SSAPy. GalSim enables users to customize a wide 
range of observing parameters, such as detector settings, optical system properties, and realistic 
atmospheric conditions. Specifically, users can manually configure sensor parameters, such as vignetting
and optical distortion, based on real and characterized sensors or simulate and test future sensor designs.

The functionality provided by SSAPy enables users to simulate the trajectories of satellites in a wide 
variety of realistic orbits, such as low-Earth or geostationary orbits. Users can also use this 
functionality to simulate arbitrary satellite maneuvers between sensor exposures.

SImSiT generates realistic astronomical background scenes using star positions, 
fluxes, and colors from the Gaia Data Release 2 catalog (:cite:t:`GaiaCollaboration2018`). 
SImSiT also provides tools to convert Gaia optical magnitudes into other passbands, 
such as Short Wave Infrared (SWIR), to generate realistic astronomical scenes 
in different wavelengths. SImSiT models satellites as a Lambertian sphere with 
spherically symmetric emissivity. Users can specify the brightness of a satellite in a 
given photometric filter, its size (radius), and albedo. 
