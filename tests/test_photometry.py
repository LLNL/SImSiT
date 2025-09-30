import pytest
import numpy as np
from astropy.table import Table
from satist.photometry import (
    convert_gaia_magnitude,
    get_gaia_sources,
    get_gaia_magnitude_histogram,
    vega_to_ab_offset
)


# Test for convert_gaia_magnitude
def test_convert_gaia_magnitude():
    gaia_g = 15.0
    gaia_bp = 16.0
    gaia_rp = 14.0

    # Test for 2MASS_Ks filter
    target_filter = '2MASS_Ks'
    result = convert_gaia_magnitude(gaia_g, gaia_bp, gaia_rp, target_filter)
    expected = 0.1885 + (-2.092 * 2.0) + (0.1345 * (2.0 ** 2))
    assert np.isclose(result, expected), f"Expected {expected}, got {result}"

    # Test for SDSS12_i filter
    target_filter = 'SDSS12_i'
    result = convert_gaia_magnitude(gaia_g, gaia_bp, gaia_rp, target_filter)
    expected = 0.29676 + (-0.64728 * 2.0) + (0.10141 * (2.0 ** 2))
    assert np.isclose(result, expected), f"Expected {expected}, got {result}"

def test_get_gaia_sources():
    ra_deg = 10.0
    dec_deg = 10.0
    radius_deg = 1.0
    star_lim = 10

    results = get_gaia_sources(ra_deg, dec_deg, radius_deg, gaia_data_release=3, star_lim=star_lim)
    assert isinstance(results, Table), "Expected results to be an Astropy Table"
    assert len(results) <= star_lim, f"Expected at most {star_lim} results, got {len(results)}"


def test_get_gaia_magnitude_histogram():
    nbins = 50
    results = get_gaia_magnitude_histogram(gaia_data_release=3, star_lim=100, nbins=nbins)

    assert isinstance(results, dict), "Expected results to be a dictionary"
    for filter in ['G', 'BP', 'RP']:
        assert filter in results, f"Expected filter {filter} in results"
        assert 'bins' in results[filter], f"Expected 'bins' key in {filter} results"
        assert 'values' in results[filter], f"Expected 'values' key in {filter} results"
        assert len(results[filter]['bins']) == nbins, f"Expected {nbins} bins for {filter}, got {len(results[filter]['bins'])}"


# Test for vega_to_ab_offset
def test_vega_to_ab_offset():
    offsets = {
        '2MASS_Ks': 1.85,
        '2MASS_H': 1.39,
        '2MASS_J': 0.91,
        'G': 0.105,
        'BP': 0.0292,
        'RP': 0.3542
    }

    for filter_name, expected_offset in offsets.items():
        result = vega_to_ab_offset(filter_name)
        assert np.isclose(result, expected_offset), f"Expected {expected_offset} for {filter_name}, got {result}"
