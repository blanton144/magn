import numpy as np


def vdisp_to_logmbh(vdisp=None, version='gultekin09'):
    """Apply M-sigma to MaNGA velocity dispersion

    Parameters
    ----------

    vdisp : np.float32, or ndarray of np.float32
        velocity dispersion (km/s)

    version : str
        version of relationship

    Returns
    -------

    logmbh : np.float32, or ndarray of np.float32
        log10 of the black hole mass (solar masses)

    Notes
    -----

    Uses relation of form: 
        log10 M = alpha + beta * log10(vdisp / 200)

    Versions available are:
      'gultekin09' : Gultekin et al (2009): alpha=8.12, beta=4.24
      'kormendy13' : Kormendy & Ho (2013): alpha=8.5, beta=4.4
      'graham13' : Graham & Scott (2013): alpha=8.14, beta=5.2
      'greene06' : Greene & Ho (2006): alpha=7.86, beta=3.65
"""
    if(version == 'gultekin09'):
        alpha = 8.12
        beta = 4.24
    elif(version == 'kormendy13'):
        alpha = 8.5
        beta = 4.4
    elif(version == 'graham13'):
        alpha = 8.14
        beta = 5.2
    elif(version == 'greene06'):
        alpha = 7.86
        beta = 3.65
    else:
        raise ValueError("No such version {v}".format(v=version))

    logmbh = alpha + beta * np.log10(vdisp / 200.)
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


def loglum_to_logbolo(loglum=None, version=''):
    """Conversion from OIII or Hb luminosity to bolometric

    Parameters
    ----------

    loglum : np.float32
        log10 OIII or Hb in erg/s

    Returns
    -------

    logbolo : np.float32
        log10 bolometric in erg/s

    Notes
    -----
    
    Available versions:
      'pennell17_nodust_o3' - based on OIII uncorrected for dust
      'heckman04_nodust_o3' - based on OIII uncorrected for dust
      'heckman04_dust_o3' - based on OIII corrected for dust
      'netzer19_dust_hb' - based on Hb corrected for dust
"""
    if(version == 'pennell17_nodust_o3'):
       logbolo = 46.058 + 0.5617 * (loglum - 42.5)
    elif(version == 'heckman04_nodust_o3'):
       logbolo = 3.544 + loglum
    elif(version == 'heckman04_dust_o3'):
       logbolo = 2.778 + loglum
    elif(version == 'netzer19_dust_hb'):
       logbolo = 45.661 + 1.18 * (loglum - 42.)
    else:
        raise ValueError("No such version {v}".format(v=version))
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

