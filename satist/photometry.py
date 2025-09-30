"""
Functions to transform fluxes and mangituedes between different photometric systems.

See: https://gea.esac.esa.int/archive/documentation/GDR2/Data_processing/chap_cu5pho/sec_cu5pho_calibr/ssec_cu5pho_PhotTransf.html#Ch5.F11
"""

#import spextre
from astroquery.gaia import Gaia
import astropy.units as u
from astropy.coordinates import SkyCoord
import numpy as np 

def convert_gaia_magnitude(gaia_g, gaia_bp, gaia_rp, target_filter='2MASS_Ks'):
    """
    Photometric transforms from the Gaia Data Release 2 data G, G_BP, G_RP into a target
    different photometric system.

    Args:
        gaia_g (float): Magnitude in the Gaia G passband.
        gaia_bp (float): Magnitude in the Gaia RP passband.
        gaia_rp (float): Magnitude in the Gaia BP passpand.
        target_fileter (str): Target photometric system to convert 
        Gaia magnitudes.

    Returns:
        Magniude in the passband specified by the target_filter argument.
    
    """
    coeff_dict = {
        '2MASS_Ks': [0.1885, -2.092, 0.1345],
        '2MASS_H' : [0.1621, -1.968, 0.1328],
        '2MASS_J' : [0.01883,-1.394, 0.07893],
        'SDSS12_i': [0.29676,-0.64728, 0.10141]
    }

    bp_minus_rp = gaia_bp - gaia_rp
    coeffs = coeff_dict[target_filter]
    return coeffs[0] + coeffs[1] * bp_minus_rp + coeffs[2] * (bp_minus_rp**2)


def spectrum_to_magnitude():
    return True


def get_gaia_sources(ra_deg, dec_deg, radius_deg, gaia_data_release=3, star_lim=-1):
    """
    Get all gaia sources  within a circle centered on ra_deg, dec_deg with
    radius radaius deg.

    Args:
        ra_deg (float): Right ascension in degress on the ICRF.
        ded_deg (float): Declination in degrees on the ICRF.
        radius_deg (float): Radius in degrees of the region to search.
        gaia_data_release (int): Gaia data release to be used (2 or 3),
        star_lim (int): Limit number of stars to be returned default -1 means no
        limit.
    Returns:
        results (astropy.Table): results of query that contain the gaia colum names.
    """
    Gaia.ROW_LIMIT = star_lim
    Gaia.MAIN_GAIA_TABLE = f"gaiadr{gaia_data_release}.gaia_source"

    coord = SkyCoord(ra=ra_deg, dec=dec_deg, unit=(u.degree, u.degree), frame='icrs')
    query = Gaia.cone_search_async(coord, radius=u.Quantity(radius_deg, u.deg))
    results = query.get_results()
    return results


def get_gaia_magnitude_histogram(gaia_data_release=3, star_lim=-1, nbins=100):
    """
    Get Histrogram of Gaia Mangnitudes in G, BP and RP.

    Args:
        gaia_data_release (int): Gaia data release to be used (2 or 3),
        star_lim (int): Limit number of stars to be returned default -1 means no
        limit.
        filter: (str): Gaia photometric filter either G, BP, or RP
        nbins (int): Number of bins in the magnitude histogram.
    Returns:

    """
    catalog = get_gaia_sources(10, 10, 1, gaia_data_release=3, star_lim=-1)
    gaia_filters = ['g', 'bp', 'rp']
    results = {}

    for filter in gaia_filters: 
        magnitude_dist = catalog[f'phot_{filter}_mean_mag'].data
        magnitude_dist = magnitude_dist[~np.isnan(magnitude_dist)]
        min_mag, max_mag = np.min(magnitude_dist), np.max(magnitude_dist)
        bins = np.linspace(min_mag, max_mag, num=nbins)
        values, _ = np.histogram(magnitude_dist, bins=bins)
        results[filter.upper()] = {}
        results[filter.upper()]['bins'] = bins
        results[filter.upper()]['values'] = values 

    return results 


def vega_to_ab_offset(filter_name):
    """
    Vega to AB magnitude system offset terms.
    For 2MASS See: https://www.astronomy.ohio-state.edu/martini.10/usefuldata.html
    For Gaia See: https://gea.esac.esa.int/archive/documentation/GDR2/Data_processing/chap_cu5pho/sec_cu5pho_calibr/ssec_cu5pho_calibr_extern.html
    Table 5.2 and 5.3.

    Args:
        filter_name (str): name of filter to get magnitude offset for.
    Returns:
        offset (float): magnitude offset (add this to VEGA to get AB)
    """

    offsets = {'2MASS_Ks': 1.85,
                '2MASS_H': 1.39,
                '2MASS_J': 0.91,
                'G': 0.105, 
                'BP': 0.0292,  
                'RP': 0.3542}
    
    return offsets[filter_name]








