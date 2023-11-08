import numpy as np


def vdisp_to_logmbh(vdisp=None):
    """Apply M-sigma to MaNGA velocity dispersion

    Parameters
    ----------

    vdisp : np.float32, or ndarray of np.float32
        velocity dispersion

    Returns
    -------

    logmbh : np.float32, or ndarray of np.float32
        log10 of the black hole mass (solar masses)

    Notes
    -----

    Uses Gultekin et al (2009):
        log10 M = 8.12 + 4.24 log10(vdisp)
"""
    logmbh = 8.12 + 4.24 * np.log10(vdisp / 200.)
    return(logmbh)


def logmbh_to_logledd(logmbh=None):
    """Calculate Eddington ratio for a black hole mass

    Parameters
    ----------

    logmbh : np.float32, or ndarray of np.float32
        log10 of the black hole mass (solar masses)

    Returns
    -------

    logledd : np.float32, or ndarray of np.float32
        log10 of the Eddington luminosity in erg/s
"""
    logledd = logmbh + 38 + np.log10(1.26)
    return(logledd)


def logo3_to_logbolo(logo3=None):
    """Conversion from OIII luminosity to bolometric

    Parameters
    ----------

    logo3 : np.float32
        log10 OIII in erg/s

    Returns
    -------

    logbolo : np.float32
        log10 bolometric in erg/s

    Notes
    -----
    
    Uses Lamastra et al (2009)
"""
    logbolo = 40. + np.log10(112.) + 1.2 * (logo3 - 40.)
    return(logbolo)


def logw2_to_logbolo(logw2=None):
    """Conversion from W2 luminosity to bolometric

    Parameters
    ----------

    logw2 : np.float32
        nu Lnu in W2 in erg/s

    Returns
    -------

    logbolo : np.float32
        log10 bolometric in erg/s

    Notes
    -----
    
    Uses Stern (2015), then multiplies by 20 per Comerford
"""
    loglx = 40.981 + 1.024 * (logw2 - 41.) - 0.047 * (logw2 - 41.)**2
    logbolo = loglx + np.log10(20.)
    return(logbolo)

