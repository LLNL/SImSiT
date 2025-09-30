# SImSiT

[![Documentation Status](https://readthedocs.org/projects/popclass/badge/?version=latest)](https://simsit.readthedocs.io/en/latest/?badge=latest)

The Satellite Image Simulation Toolkit (SImSiT) is a Python software package designed 
to generate diverse and realistic satellite imaging scenarios. SImSiT is built on top 
of [GalSim](https://github.com/GalSim-developers/GalSim) and [SSAPy](https://github.com/LLNL/SSAPy) and serves as a toolkit for 
simulating imaging data that supports the development and testing of algorithms for 
satellite detection, calibration, and characterization.

For more details on the project including the installation, contributing, and the getting started guide see the [documentation](https://simsit.readthedocs.io).

![image info](./docs/source/example_satist_out.png)

Example of simulated image output from SImSiT of satellites in low-Earth orbit observed in the optical i-band filter using a ground-based sensor:
Left: An example image taken with sidereal (star) tracking, where the fast-moving satellite appears as a streak across the image.
Right: An example image taken with satellite tracking, where the stars appear as streaks, and the satellite appears as a point source.

## License

SImSiT is distributed under the terms of the MIT license. All new contributions must be made under the MIT license.

See Link to [license](https://github.com/LLNL/simsit/blob/main/LICENSE) and [NOTICE](https://github.com/LLNL/simsit/blob/main/NOTICE) for details.

SPDX-License-Identifier: MIT

LLNL-CODE-2009738
