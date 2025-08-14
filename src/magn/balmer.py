import numpy as np


def whitford(wave=None):
    """Return Whitford extinction

    Parameters
    ----------

    wave : np.float32, or ndarray of np.float32
        wavelength(s) in Ang

    Returns
    -------

    ex : np.float32, or ndarray of np.float32
        extinction in mag
"""
    if(type(wave) == np.ndarray):
        ex = np.zeros(len(wave), dtype=np.float32)
        iwave = np.where(1. / wave * 10000. <= 2.29)[0]
        ex[iwave] = 0.74 / wave[iwave] * 10000. - 0.34
        iwave = np.where(1. / wave * 10000. > 2.29)[0]
        ex[iwave] = 0.43 / wave[iwave] * 10000. + 0.37
    else:
        if(1./ wave * 10000. <= 2.29):
            ex = 0.74 / wave * 10000. - 0.34
        else:
            ex = 0.43 / wave * 10000. + 0.37
    return(ex)
    
def balmer_to_av(hahb=None, hahb_nominal=2.85, minav=0.):
    """Return AV estimated given Balmer decrement

    Parameters
    ----------

    hahb : np.float32
        Ha / Hb ratio

    hahb_nominal : np.float32
        assumed Ha / Hb ratio in absence of dust

    minav : np.float32
        minimum A_V 

    Returns
    -------

    av : np.float32

    Notes
    -----

    Assuming Whitford (1958) extinction curve, as parametrized by Miller & Mathews (1972),
    return the AV associated with a given Ha/Hb
"""
    ex_ha = whitford(wave=6563.)
    ex_hb = whitford(wave=4861.)

    # reddening of Halpha / Hbeta observed
    hahb_mag = 2.5 * np.log10(hahb / hahb_nominal)  # dust makes this number > 0

    # extinction curve normalization to match
    exnorm = hahb_mag / (ex_hb - ex_ha)

    a_v = whitford(wave=5500.) * exnorm

    if(type(a_v) == np.ndarray):
        ilz = np.where(a_v < minav)[0]
        a_v[ilz] = minav
    else:
        if(a_v < minav):
            a_v = minav

    return(a_v)
    
def dust_correction(a_v=None, hahb=None, wave=None, hahb_nominal=2.85, minav=0.):
    """Return dust corrections estimated given Balmer decrement

    Parameters
    ----------

    a_v : np.float32
        A_v (or None if hahb is given)

    hahb : np.float32
        Ha / Hb ratio (or None if hahb is given)

    wave : np.float32, or ndarray of np.float32
        wavelength(s) in Ang

    hahb_nominal : np.float32
        assumed Ha / Hb ratio in absence of dust

    minav : np.float32
        minimum A_V 

    Returns
    -------

    factor : np.float32
        factor to multiply fluxes at each wavelength by

    Notes
    -----

    If A_V is supplied, uses that, otherwise requires hahb input.

    Assuming Whitford (1958) extinction curve, as parametrized by Miller & Mathews (1972),
    return the factor to multiply fluxes at each wavelength
"""
    if(a_v is None):
        a_v = balmer_to_av(hahb=hahb, hahb_nominal=hahb_nominal, minav=minav)

    ex_v = whitford(wave=5500.)

    exnorm = a_v / ex_v
    
    ex = whitford(wave=wave)
    a = ex * exnorm

    factor = 10.**(0.4 * a)

    return(factor)
