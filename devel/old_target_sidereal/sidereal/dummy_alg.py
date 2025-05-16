# Start with dummy_detect_target_alg.py
    # Get x, y and flux and pass to Bob's code

# Run Bob's code
    # get all star information
    # Add flux calibration into header, so that we can use it later and/or calculate mag from flux on the fly
    # convert each x, y into ra, dec

# Turn all information into correct yaml format to generate submission
    # Include header with branch, comp_name, etc
    # star ra, dec = coord_ra, coord_dec
    # star mag = base_PsfFlux_mag
    # star flux = base_PsfFlux_instFlux
    # sat x,  y from dummy script
    # sat flux from dummy script
    # sat ra, dec = xy after wcs
    # sat mag = flux after calibration factor applied
    # Make sure all units are converted into correct format for scoring