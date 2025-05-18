"""
Functions to transform fluxes and mangituedes between Gaia and different photometric systems.

See: https://gea.esac.esa.int/archive/documentation/GDR2/Data_processing/chap_cu5pho/sec_cu5pho_calibr/ssec_cu5pho_PhotTransf.html#Ch5.F11
"""

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
