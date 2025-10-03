===============
Getting Started
===============

This page shows a user how to genreate simluated images with satillies using SImSiT.

Simulation
==========

To simulate the sidereal branch:
---------------------------------

- Edit ``config/sidereal.yaml``  
  - (For full simulation: ``n_obs = 1010``, ``n_demo = 10``, ``develop=False``)
- Run: ``python scripts/simulate.py config/sidereal.yaml``
- This will create a folder with all simulated data under ``branches/tracking/sidereal``

To simulate the target branch:
------------------------------

- Edit ``config/target.yaml``  
  - (For full simulation: ``n_obs = 1010``, ``n_demo = 10``, ``develop=False``)
- Run: ``python scripts/simulate.py config/target.yaml``
- This will create a folder with all simulated data under ``branches/tracking/target``
- To simulate target tracking images with no satellite, run:  
  ``python scripts/simulate.py config/empty_target.yaml``

To simulate the IOD branch:
---------------------------

- **NOTE:** To simulate this branch you must download  
  ``ODJobs_Simulated_Data_20001_20005_v2/`` from the eSTM Teams page  
  (found under General/Files/) and place it in the ``branches/iod/`` directory.
- Run: ``python scripts/simulate_iod.py``
- This will create a folder with all simulated data under ``branches/iod``
- (See ``branches/iod/README_private.md`` for more info)

To simulate the track images branches:
--------------------------------------

- Run: ``python scripts/simulate_tracks.py``  
  - (For full simulation pass ``--nsat 105``)

- **Sidereal tracking:**  
  - Run: ``python scripts/simulate_track_images.py --config config/sidereal_track.yaml branches/track_images/track_obs/sat-obs-*.fits``  
    - (For full simulation pass ``--stride 12``)  
    - To run as the non-tracking branch, run with ``--nobs=1`` (this simulates just one image per satellite)
  - This will create a folder with all simulated data under ``branches/track_images/sidereal_track``

- **Target tracking:**  
  - Run: ``python scripts/simulate_track_images.py --config config/target_track.yaml branches/track_images/track_obs/sat-obs-*.fits``  
    - (For full simulation pass ``--stride 12``)
  - This will create a folder with all simulated data under ``branches/track_images/target_track``