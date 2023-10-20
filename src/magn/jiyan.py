import os
import glob
import numpy as np
import astropy.io.fits
import matplotlib.pyplot as plt
import scipy.interpolate
import fitsio
import mnsa.mnsautils


def read_central_flux(version=None, dr17=None):
    """Read central flux file
    
    Parameters
    ----------

    version : str
        catalog version

    dr17 : bool
        if set, use DR17 data

    Returns
    -------

    cf : dict of ndarrays of nd.float32
        central fluxes for channels

    cf_ivar : dict of ndarrays of np.float32
        inverse variance of cf

    cf_nbad : dict of ndarrays of np.int32
        number of bad pixels in cf determination

    good : ndarray of bool
        whether observation should be considered good

    Comments
    --------

    "good" requires an inverse variance positive and fewer than 4 bad pixels.
    
""" 
    drpall = fitsio.read(os.path.join(os.getenv('MNSA_DATA'),
                                      '{v}', 'manga', 'redux', '{v}',
                                      'drpall-{v}.fits').format(v=version))

    cffile = os.path.join(os.getenv('MNSA_DATA'), version + '.analysis',
                          'central_flux',
                          'central-flux-{version}.fits'.format(version=version))
    if(dr17):
        cffile = cffile.replace('central-flux-', 'central-flux-dr17-')
    central_flux = fitsio.read(cffile, ext='FLUXES')
    channels = fitsio.read(cffile, ext='CHANNELS_EMLINE_GFLUX')

    ichannel = dict()
    for channel in ['Hb-4862', 'OIII-5008', 'Ha-6564', 'NII-6585',
                    'SII-6718', 'SII-6732']:
        ichannel[channel] = np.where(channels['channel'] == channel)[0][0]

    cf = dict()
    cf_ivar = dict()
    cf_nbad = dict()
    for channel in ichannel:
        cf[channel] = central_flux['central_emline_gflux'][:, ichannel[channel]]
        cf_ivar[channel] = central_flux['central_emline_gflux_ivar'][:, ichannel[channel]]
        cf_nbad[channel] = central_flux['central_emline_gflux_nbad'][:, ichannel[channel]]
    cf['SII-6718+6732'] = cf['SII-6718'] + cf['SII-6732']
    cf_ivar['SII-6718+6732'] = np.zeros(len(cf['SII-6718']), dtype=np.float32)
    igd = np.where((cf_ivar['SII-6718'] > 0) & (cf_ivar['SII-6732'] > 0))[0]
    cf_ivar['SII-6718+6732'][igd] = 1. / (1. / cf_ivar['SII-6718'][igd] +
                                          1. / cf_ivar['SII-6732'][igd])

    good = np.ones(len(central_flux), dtype=bool)
    good = good & (drpall['z'] > 0)
    for channel in ichannel:
        good = (good & (cf_ivar[channel] > 0.) & (cf_nbad[channel] < 4))

    return(cf, cf_ivar, cf_nbad, good)


def ratios_to_pspace(n2ha=None, s2ha=None, o3hb=None,
                     n2ha_err=None, s2ha_err=None, o3hb_err=None):
    """Converts log line ratios to Ji & Yan P1/P2/P3 space

    Parameters
    ----------

    n2ha, n2ha_err : np.float32 or ndarray of np.float32
        log [NII] 6583 / H-alpha flux ratio

    s2ha, s2ha_err : np.float32 or ndarray of np.float32
        log ([SII] 6716 + 6730) / H-alpha flux ratio

    o3hb, o3hb_err : np.float32 or ndarray of np.float32
        log ([OIII] 5008) / H-beta flux ratio

    Returns
    -------

    p1 : np.float32 or ndarray of np.float32
        P1 component (SF or AGN-ness)

    p2 : np.float32 or ndarray of np.float32
        P2 component (metallicity-ish)

    p3 : np.float32 or ndarray of np.float32
        P3 component (ionization-ish?)

    p1_err : np.float32 or ndarray of np.float32
        error in P1 component (SF or AGN-ness)

    p2_err : np.float32 or ndarray of np.float32
        error in P2 component (metallicity-ish)

    p3_err : np.float32 or ndarray of np.float32
        error in P3 component (ionization-ish?)

    Notes
    -----

    p1_err, p2_err, p3_err are only returned if all input errors are supplied.
"""
    p1 = 0.63 * n2ha + 0.51 * s2ha + 0.59 * o3hb
    p2 = - 0.63 * n2ha + 0.78 * s2ha
    p3 = - 0.46 * n2ha - 0.37 * s2ha + 0.81 * o3hb
    if((n2ha_err is None) |
       (s2ha_err is None) |
       (o3hb_err is None)):
        return(p1, p2, p3)
    p1_err = np.sqrt(0.63**2 * n2ha_err**2 +
                     0.51**2 * s2ha_err**2 +
                     0.59**2 * o3hb_err**2)
    p2_err = np.sqrt(0.63**2 * n2ha_err**2 +
                     0.78**2 * s2ha_err**2)
    p3_err = np.sqrt(0.46**2 * n2ha_err**2 +
                     0.37**2 * s2ha_err**2 +
                     0.81**2 * o3hb_err**2)
    return(p1, p2, p3, p1_err, p2_err, p3_err)


def would_be_agn(redshift=None, cf=None, cf_ivar=None, agn_ratios=None,
                 log_luminosity_o3=None):
    """Return whether a given luminosity will pass the limit to be classified as AGN

    Parameters
    ----------

    redshift : np.float32
        redshift of galaxy

    cf : dict of np.float32
        fluxes for each line

    cf_ivar : dict of np.float32
        inverse variance of fluxes for each line

    agn_ratios : dict of np.float32
        ratios of lines relative to [OIII] for typical AGN

    log_luminosity_o3 : np.float32
        limiting [OIII] luminosity
    
    Returns
    -------

    detected : ndarray of bool
        are the necessary lines detected

    isagn : ndarray of bool
        is this an identified AGN

    lines : ndarray
        lines structure returned by select_agn()

    Notes
    -----

    cf and agn_ratios must have the keys 'NII-6585', 'SII-6718+6732',
    'SII-67I8', 'SII-6732', 'OIII-5008', 'Ha-6564', and 'Hb-4862'

    cf_ivar must have the keys 'NII-6585', 'SII-67I8', 'SII-6732',
    'OIII-5008', 'Ha-6564', and 'Hb-4862'
"""
    logterm = mnsa.mnsautils.log_flux_to_luminosity(redshift=redshift)
    log_flux_o3 = log_luminosity_o3 - logterm

    ncf = dict()
    for channel in cf:
        ncf[channel] = cf[channel] + (10.**log_flux_o3) * agn_ratios[channel]
    detected, isagn, lines = select_agn(cf=ncf, cf_ivar=cf_ivar,
                                        good=np.ones(1, dtype=bool))
    return(detected, isagn, lines)


def select_agn(cf=None, cf_ivar=None, good=None, p1_limit=-0.3, p3_limit=0.5):
    """Select AGN using P1 and P3

    Parameters
    ----------

    cf : dict of ndarrays of np.float32
        fluxes for each line

    cf_ivar : dict of ndarrays of np.float32
        inverse variance of fluxes for each line

    good : ndarray of bool
        whether to check each galaxy

    p1_limit : np.float32
        lower limit of P1 to be an AGN

    p3_limit : np.float32
        lower limit of P3 to be an AGN

    Returns
    -------

    detected : ndarray of bool
        are the necessary lines detected

    isagn : ndarray of bool
        is this an identified AGN

    lines : ndarray
       contains line measurements and errors

    Notes
    -----

    cf must have the keys 'NII-6585', 'SII-6718+6732',
    'SII-67I8', 'SII-6732', 'OIII-5008', 'Ha-6564', and 'Hb-4862'

    cf_ivar must have the keys 'NII-6585', 'SII-67I8', 'SII-6732',
    'OIII-5008', 'Ha-6564', and 'Hb-4862'

    lines has: 'n2ha', 'n2ha_err', 's2ha', 's2ha_err', 'o3hb', 'o3hb_err',
      'p1', 'p1_err', 'p2', 'p2_err', 'p3', 'p3_err'
"""
    detected = good.copy()
    for channel in cf_ivar:
        detected[good] = (detected[good] & (cf[channel][good] *
                                            np.sqrt(cf_ivar[channel][good]) > 2.))

    n2ha = np.log10(cf['NII-6585'][detected] / cf['Ha-6564'][detected])
    n2ha_err = (1. / np.log(10.)) * np.sqrt(1. / (cf_ivar['NII-6585'][detected] * cf['NII-6585'][detected]**2) +
                                            1. / (cf_ivar['Ha-6564'][detected] * cf['Ha-6564'][detected]**2))
                                                  
    s2ha = np.log10(cf['SII-6718+6732'][detected] / cf['Ha-6564'][detected])
    s2ha_err = (1. / np.log(10.)) * np.sqrt(1. / (cf_ivar['SII-6718+6732'][detected] * cf['SII-6718+6732'][detected]**2) +
                                            1. / (cf_ivar['Ha-6564'][detected] * cf['Ha-6564'][detected]**2))

    o3hb = np.log10(cf['OIII-5008'][detected] / cf['Hb-4862'][detected])
    o3hb_err = (1. / np.log(10.)) * np.sqrt(1. / (cf_ivar['OIII-5008'][detected] * cf['OIII-5008'][detected]**2) +
                                            1. / (cf_ivar['Hb-4862'][detected] * cf['Hb-4862'][detected]**2))

    p1, p2, p3, p1_err, p2_err, p3_err = ratios_to_pspace(n2ha=n2ha,
                                                          s2ha=s2ha,
                                                          o3hb=o3hb,
                                                          n2ha_err=n2ha_err,
                                                          s2ha_err=s2ha_err,
                                                          o3hb_err=o3hb_err)
    
    lines_dtype = np.dtype([('n2ha', np.float32),
                            ('n2ha_err', np.float32),
                            ('s2ha', np.float32),
                            ('s2ha_err', np.float32),
                            ('o3hb', np.float32),
                            ('o3hb_err', np.float32),
                            ('p1', np.float32),
                            ('p1_err', np.float32),
                            ('p2', np.float32),
                            ('p2_err', np.float32),
                            ('p3', np.float32),
                            ('p3_err', np.float32)])

    lines = np.zeros(len(good), dtype=lines_dtype)
    lines['n2ha'] = - 9999
    lines['n2ha_err'] = - 9999
    lines['s2ha'] = - 9999
    lines['s2ha_err'] = - 9999
    lines['o3hb'] = - 9999
    lines['o3hb_err'] = - 9999
    lines['p1'] = - 9999
    lines['p1_err'] = - 9999
    lines['p2'] = - 9999
    lines['p2_err'] = - 9999
    lines['p3'] = - 9999
    lines['p3_err'] = - 9999

    lines['n2ha'][detected] = n2ha
    lines['n2ha_err'][detected] = n2ha_err
    lines['s2ha'][detected] = s2ha
    lines['s2ha_err'][detected] = s2ha_err
    lines['o3hb'][detected] = o3hb
    lines['o3hb_err'][detected] = o3hb_err
    lines['p1'][detected] = p1
    lines['p1_err'][detected] = p1_err
    lines['p2'][detected] = p2
    lines['p2_err'][detected] = p2_err
    lines['p3'][detected] = p3
    lines['p3_err'][detected] = p3_err

    isagn = ((detected) & (lines['p1'] > p1_limit) &
             (lines['p3'] > p3_limit))

    return(detected, isagn, lines)


def find_o3_limit(redshift=None, cf=None, cf_ivar=None, agn_ratios=None):
    """Find [OIII] luminosity limit

    Parameters
    ----------

    redshift : np.float32
        redshift of galaxy

    cf : dict of ndarray of np.float32
        fluxes for each line (1-element arrays)

    cf_ivar : dict of ndarray of np.float32
        inverse variance of fluxes for each line (1-element arrays)

    agn_ratios : dict of ndarray of np.float32
        ratios of lines relative to [OIII] for typical AGN

    Returns
    -------

    log_luminosity_o3_limit : np.float32
        limiting [OIII] luminosity

    Notes
    -----

    cf and agn_ratios must have the keys 'NII-6585', 'SII-6718+6732',
    'SII-67I8', 'SII-6732', 'OIII-5008', 'Ha-6564', and 'Hb-4862'

    cf_ivar must have the keys 'NII-6585', 'SII-67I8', 'SII-6732',
    'OIII-5008', 'Ha-6564', and 'Hb-4862'
"""
    bounds = np.array([15., 50.], dtype=np.float32)
    detected, isagn0, lines = would_be_agn(redshift=redshift, cf=cf, cf_ivar=cf_ivar, agn_ratios=agn_ratios,
                                           log_luminosity_o3=bounds[0])
    if(isagn0 == True):
        return(- 9999.)
    detected, isagn1, lines = would_be_agn(redshift=redshift, cf=cf, cf_ivar=cf_ivar, agn_ratios=agn_ratios,
                                           log_luminosity_o3=bounds[1])
    if(isagn1 == False):
        return(- 9999.)
    while((bounds[1] - bounds[0]) > 1.e-4):
        bounds_middle = 0.5 * (bounds[1] + bounds[0])
        detected, isagn, lines = would_be_agn(redshift=redshift, cf=cf, cf_ivar=cf_ivar, agn_ratios=agn_ratios,
                                              log_luminosity_o3=bounds_middle)
        if(isagn):
            bounds[1] = bounds_middle
        else:
            bounds[0] = bounds_middle
    log_luminosity_o3_limit = 0.5 * (bounds[1] + bounds[0])
    return(log_luminosity_o3_limit)


class JiYan(object):
    """Class for Ji & Yan models

    Parameters
    ----------

    Attributes
    ----------

    Notes
    -----

"""
    def __init__(self):
        self.name = None

        self.lines = dict()
        self.lines['o3'] = ['O__3__5006']
        self.lines['hb'] = ['H__1__4861']
        self.lines['n2'] = ['N__2__6583']
        self.lines['ha'] = ['H__1__6562']
        self.lines['s2'] = ['S__2__6716', 'S__2__6730']
        self.lines['o2'] = ['BLND__3727']
        self.lines['o1'] = ['O__1__6300']

        self.raw_grids = None
        self.raw_metallicity_o = None
        self.raw_log_ionization = None

        self.grids = None
        self.metallicity_o = None
        self.log_ionization = None
        return

    def set_pspace(self):
        self.n2ha = np.log10(self.grids['n2'] / self.grids['ha'])
        self.o3hb = np.log10(self.grids['o3'] / self.grids['hb'])
        self.s2ha = np.log10(self.grids['s2'] / self.grids['ha'])

        self.p1, self.p2, self.p3 = ratios_to_pspace(n2ha=self.n2ha, s2ha=self.s2ha, o3hb=self.o3hb)
        return

    def read_model(self, modelname=None):
        """Read in a particular model"""
        # Load photoionization model grids
        self.name = modelname
        modeldir = os.path.join(os.getenv('MNSA_DIR'), 'data', 'jiyan')
        nameList = []
        n_model = len(glob.glob(os.path.join(modeldir, modelname) + '*'))
        for i in np.arange(n_model, dtype=np.int32):
            nameList += [os.path.join(modeldir, modelname) + str(i) + '_line.fits']

        # Hack for metallicity grid
        if('bpl' in modelname):
            self.raw_metallicity_o = np.array([-0.75, -0.5, -0.25, 0.,
                                               0.25, 0.5, 0.75],
                                              dtype=np.float32)
        else:
            self.raw_metallicity_o = np.array([-1.3, -0.7, -0.4, 0., 0.3, 0.5],
                                              dtype=np.float32)

        for iname, name in enumerate(nameList):
            hdu = astropy.io.fits.open(name)

            # extract emission line fluxes
            tmp = dict()
            for line in self.lines:
                tmp[line] = 0.
                for d in self.lines[line]:
                    tmp[line] = tmp[line] + hdu[1].data[d]

            if(self.raw_grids is None):
                self.raw_grids = dict()
                for line in self.lines:
                    nper = len(tmp[line])
                    self.raw_grids[line] = np.zeros((n_model, nper),
                                                    dtype=np.float32)

            for line in self.lines:
                self.raw_grids[line][iname, :] = tmp[line]

            self.raw_log_ionization = hdu[1].data['IONIZATION']

        self.grids = dict()
        refine = 5
        for line in self.lines:
            n0_raw = self.raw_grids[line].shape[0]
            n1_raw = self.raw_grids[line].shape[1]
            x0_raw = np.arange(n0_raw, dtype=np.float32)
            x1_raw = np.arange(n1_raw, dtype=np.float32)
            ginterp = scipy.interpolate.RegularGridInterpolator((x0_raw, x1_raw),
                                                                np.log10(self.raw_grids[line]),
                                                                method='cubic')
            ointerp = scipy.interpolate.CubicSpline(x0_raw,
                                                    self.raw_metallicity_o)
            iinterp = scipy.interpolate.CubicSpline(x1_raw,
                                                    self.raw_log_ionization)

            n0 = (n0_raw - 1) * refine + 1
            x0 = np.arange(n0, dtype=np.float32)
            x0 = x0 / np.float32(n0 - 1) * (n0_raw - 1)
            n1 = (n1_raw - 1) * refine + 1
            x1 = np.arange(n1, dtype=np.float32)
            x1 = x1 / np.float32(n1 - 1) * (n1_raw - 1)
            x1g, x0g = np.meshgrid(x1, x0)
            self.grids[line] = 10.**(ginterp((x0g, x1g)))
            self.metallicity_o = ointerp(x0)
            self.log_ionization = iinterp(x1)

        return

    def plot_mesh(self, x=None, y=None, mask=None, **kwargs):
        """Plot a 2D mesh

        Parameters
        ----------

        x : 2D ndarray of np.float32
           X-axis values

        y : 2D ndarray of np.float32
           Y-axis values

        mask : 2D ndarray of bool
           plot only these points
"""
        n0 = x.shape[0]
        n1 = x.shape[1]
        if(mask is not None):
            outx = np.ma.masked_where(mask == False, x)
            outy = np.ma.masked_where(mask == False, y)
        else:
            outx = x
            outy = y

        for i0 in np.arange(n0, dtype=np.int32):
            plt.plot(outx[i0, :], outy[i0, :], **kwargs)
        for i1 in np.arange(n1, dtype=np.int32):
            plt.plot(outx[:, i1], outy[:, i1], **kwargs)

        return
