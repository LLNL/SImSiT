># **Overview**
This package contains simulation scripts to create a variety of challenge data sets designed to test the robustness and accuracy of algorithms relevant to satellite detection, calibration, and characterization.

# **Installation**

- Install ssapy
- Install lsst_sphgeom
- run: `pip install -r requirements.txt`
- run: `python setup.py install`

# **Simulation**

**To simulate the sidereal branch:**
    
- edit `config/sidereal.yaml`
    - (For full simulation: n_obs = 1010, n_demo = 10, develop=False)
- run: `python scripts/simulate.py config/sidereal.yaml`
- This will create a folder with all simulated data under `branches/tracking/sidereal`

**To simulate the target branch:**
    
- edit `config/target.yaml`
    - (For full simulation: n_obs = 1010, n_demo = 10, develop=False)
- run: `python scripts/simulate.py config/target.yaml`
- This will create a folder with all simulated data under `branches/tracking/target`
- To simulate target tracking images with no satellite, run `python scripts/simulate.py config/empty_target.yaml`

**To simulate the IOD branch:**
    
- NOTE: To simulate this branch you must download `ODJobs_Simulated_Data_20001_20005_v2/` from the eSTM Teams page (found under General/Files/) and place it in the `branches/iod/` directory
- run: `python scripts/simulate_iod.py`
- This will create a folder with all simulated data under `branches/iod`
- (See `branches/iod/README_private.md` for more info)

**To simulate the track images branches:**

- run: `python scripts/simulate_tracks.py`
    - (For full simulation pass `--nsat 105`)
    
- Sidereal tracking:
    - run: `python scripts/simulate_track_images.py --config config/sidereal_track.yaml branches/track_images/track_obs/sat-obs-*.fits`
        - (For full simulation pass `--stride 12`)
        - To run as the non-tracking branch, run with `--nobs=1` (this simulates just one image per satellite)
    - This will create a folder with all simulated data under `branches/track_images/sidereal_track`
- Target tracking:
    - run: `python scripts/simulate_track_images.py --config config/target_track.yaml branches/track_images/track_obs/sat-obs-*.fits`
        - (For full simulation pass `--stride 12`)
    - This will create a folder with all simulated data under `branches/track_images/target_track`

**After all simulations are complete:**

- Check data: `python tests/check_data.py config/[config_file] [num_sats]`
- Package data: `./package_data.sh`

# **Information**

- A dashboard with competitor scores for each branch is located at https://xfiles.llnl.gov.
- Current slides reside on the eSTM Teams page (found under General/Files/Presentations/xfiles_v[X.X].pptx)
- To request an account or provide feedback, please email xfiles@llnl.gov
- For a detailed list of what changed between data versions, please consult `change_log.md`

#### **Errors:**

- If there is a gcc error when installing `lsst_sphgeom`:
    - Install using `CC=/usr/tce/packages/gcc/gcc-8.1.0/bin/gcc python -m pip install git+https://github.com/lsst/sphgeom`
    - Everytime you run `simulate.py` include `LD_PRELOAD=/usr/tce/packages/gcc/gcc-8.1.0/lib64/libstdc++.so.6` before 