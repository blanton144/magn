# MaNGA Narrow-Line AGN

## Introduction

This directory contains the files used in the analysis of Blanton et al 
(2026) of narrow-line AGN in the MaNGA data set, using the techniques of 
Ji & Yan (2020). These include the catalogs of AGN luminosities, thresholds,
and other derived quantites, and the inferred distributions of luminosities
and Eddington ratios (ERs).

Please read the paper before using the data set. In particular, the 
"detected" luminosity and ER distributions are strongly affected by 
selection effects so do not reflect the true distributions of those 
quantities.

Most of the analysis here is based on an alternative set of MaNGA cubes
using the methods of Liu et al. (2020). The differences between these
results and using the standard DR17 cubes is fairly minor, but include
the analysis with the latter for comparison. 

Where used in the paper, the star-formation rates and stellar masses 
are based on the DR17 Pipe3D analysis and use `log_SFR_SF` and `log_Mass`
respectively.

Our analysis results are in the following subdirectories:

* **standard** contains our standard results, using Liu et al (2020) cubes. 

* **dr17** contains the results using the DR17 version of the MaNGA cubes.

* **nsigma3** contains the results using a 3-sigma criterion for the detection of lines

* **p1p3** contains the results using an alternative selection criterion in P1-P3 space

Within each subdirectory there is a single AGN selection, and multiple results for 
Fagn as a function of luminosity or Eddington ratio, based on different methods to 
calculate these quantities.

The file naming convention and file format is the same within each subdirectory, and
each has the following files:

| File name                                 | Contents                         |
|-------------------------------------------|----------------------------------|
| `manga-nlr-agn-0.3.2.fits`                | AGN luminosities and thresholds  |
| `manga-nlr-agn-params-0.3.2.fits`         | Eddington ratios and thresholds  |
| `manga-nlr-lum-[lum]-fagn-0.3.2.fits`     | Detected luminosity distribution | 
| `manga-nlr-lum-[lum]-Fagn-0.3.2.fits`     | Modeled luminosity distribution  |
| `manga-nlr-edr-[ms]-[lb]-fagn-0.3.2.fits` | Detected ER distribution         |
| `manga-nlr-edr-[ms]-[lb]-Fagn-0.3.2.fits` | Modeled ER distribution          |

Various choices can be made for:

* Line luminosity [lum]: `hb_corr` or `o3`

* M-sigma relation [ms]: `kormendy13`, `gultekin09`, `graham13`, or `greene06`

* Bolometric luminosity [lb]: `netzer19_dust_hb`, `heckman04_nodust_o3`, or `heckman04_dust_o3`

## AGN luminosities and thresholds

The `manga-nlr-agn-0.3.2.fits` files have header parameters and a single table in 
`HDU2`, whose rows refer to the same MaNGA Plate-IFUs in the same order as the 
DRPAll file in the DR17 MaNGA release.

The header parameters are:

| Name                 |  Description                                      |
|----------------------|---------------------------------------------------|
| `MNSA_CAT`           | Catalog version                                   |
| `P1_THRESHOLD`       | Threshold used for P1                             |
| `P3_THRESHOLD`       | Threshold used for P3                             |
| `NSIGMA_THRESHOLD`   | Detection threshold for lines                     |
| `OIII-5008`          | Relative flux of [OIII] 5008 for canonical AGN    |
| `NII-6585`           | Relative flux of [NII] 6585 for canonical AGN     |
| `SII-6718`           | Relative flux of [SII] 6718 for canonical AGN     |
| `SII-6732`           | Relative flux of [SII] 6732 for canonical AGN     |
| `SII-6718+6732`      | Relative flux of [SII] 6718,6732for canonical AGN |
| `HA-6564`            | Relative flux of Halpha for canonical AGN         |
| `HB-6564`            | Relative flux of Hbeta for canonical AGN          |
| `N2HA`               | Log NII / Halpha for canonical AGN                |
| `S2HA`               | Log SII / Halpha for canonical AGN Catalog        |
| `O3HB`               | Log OIII / Hbeta for canonical AGN Catalog        |
| `O3HA`               | Log OIII / Halpha for canonical AGN Catalog       |
| `O3_DUST_CORRECTION` | Median OIII dust correction for AGN               |
| `HA_DUST_CORRECTION` | Median Halpha dust correction for AGN             |
| `HB_DUST_CORRECTION` | Median Halpha dust correction for AGN version     |
| `P1`                 | P1 for canonical AGN                              |
| `P2`                 | P2 for canonical AGN                              |
| `P3`                 | P3 for canonical AGN                              |

The table contains:

| Name                               | Data Type | Unit                | Description                                          |
|------------------------------------|-----------|---------------------|------------------------------------------------------|
| `plateifu`                         | `str`     |                     | MaNGA Plate-IFU                                      |
| `z`                                | `float32` |                     | Redshift from DRPAll                                 |
| `isagn`                            | `bool`    |                     | Is an AGN?                                           |
| `good`                             | `bool`    |                     | Are measurements good?                               |
| `detected`                         | `bool`    |                     | Are all lines detected?                              |
| `p1`                               | `float32` | dex                 | P1 (Ji and Yan 2020)                                 |
| `p1_err`                           | `float32` | dex                 | Error in P1                                          |
| `p2`                               | `float32` | dex                 | P2 (Ji and Yan 2020)                                 |
| `p2_err`                           | `float32` | dex                 | Error in P2                                          |
| `p3`                               | `float32` | dex                 | P3 (Ji and Yan 2020)                                 |
| `p3_err`                           | `float32` | dex                 | Error in P3                                          |
| `o3`                               | `float32` | 10^{-17} erg/cm^2/s | Flux in [OIII] 5008                                  |
| `o3_err`                           | `float32` | 10^{-17} erg/cm^2/s | Error in `o3`                                        |
| `n2ha`                             | `float32` | dex                 | Log of [NII] 6585 / Halpha ratio                  |
| `n2ha_err`                         | `float32` | dex                 | Error in `n2ha`                                      |
| `s2ha`                             | `float32` | dex                 | Log of [SII] 6718,6732 / Halpha ratio             |
| `s2ha_err`                         | `float32` | dex                 | Error in `s2ha`                                      |
| `o3hb`                             | `float32` | dex                 | Log of [OIII] 5008 / Hbeta ratio                  |
| `o3hb_err`                         | `float32` | dex                 | Error in `o3hb`                                      |
| `hahb`                             | `float32` | dex                 | Log of Halpha / Hbeta ratio                    |
| `hahb_err`                         | `float32` | dex                 | Error in `hahb`                                      |
| `hahb_av`                          | `float32` | mag                 | $A_V$ inferred from `hahb`                           |
| `o3_dust_correction`               | `float32` |                     | Dust correction for [OIII]                           |
| `ha_dust_correction`               | `float32` |                     | Dust correction for Halpha                        |
| `hb_dust_correction`               | `float32` |                     | Dust correction for Hbeta                         |
| `log_luminosity_o3`                | `float32` | log10(erg/s)        | Log of [OIII] luminosity                             |
| `log_luminosity_o3_threshold`      | `float32` | log10(erg/s)        | Log of [OIII] luminosity threshold                   |
| `log_luminosity_o3_corr`           | `float32` | log10(erg/s)        | Log of [OIII] dust-corrected luminosity              |
| `log_luminosity_o3_corr_threshold` | `float32` | log10(erg/s)        | Log of [OIII] dust-corrected luminosity threshold    |
| `log_luminosity_ha`                | `float32` | log10(erg/s)        | Log of Halpha luminosity                          |
| `log_luminosity_ha_threshold`      | `float32` | log10(erg/s)        | Log of Halpha luminosity threshold                |
| `log_luminosity_ha_corr`           | `float32` | log10(erg/s)        | Log of Halpha dust-corrected luminosity           |
| `log_luminosity_ha_corr_threshold` | `float32` | log10(erg/s)        | Log of Halpha dust-corrected luminosity threshold |
| `log_luminosity_hb`                | `float32` | log10(erg/s)        | Log of Hbeta luminosity                           |
| `log_luminosity_hb_threshold`      | `float32` | log10(erg/s)        | Log of Hbeta luminosity threshold                 |
| `log_luminosity_hb_corr`           | `float32` | log10(erg/s)        | Log of Hbeta dust-corrected luminosity            |
| `log_luminosity_hb_corr_threshold` | `float32` | log10(erg/s)        | Log of Hbeta dust-corrected luminosity threshold  |

## AGN bolometric luminosies and black hole masses

The `manga-nlr-agn-params-0.3.2.fits` files have a three tables.

`HDU2` has inferred parameters regarding the AGN black hole masses and bolometric
luminosities. Each row refers to the same MaNGA Plate-IFU as the corresponding
row of the `manga-nlr-agn-params-0.3.2.fits` file

| Name                | Data Type    | Unit             | Description                                                     |
|---------------------|--------------|------------------|-----------------------------------------------------------------|
| `plateifu`          | `str`        |                  | MaNGA Plate-IFU                                                 |
| `vdisp`             | `float32`    | km/s             | Velocity dispersion within 1Re                                  |
| `logmbh`            | `float32[4]` | log solar masses | Black hole masses from M-sigma relations                        |
| `logledd`           | `float32[4]` | log erg/s        | Eddington luminosities based on black hole masses               |
| `logbolo`           | `float32[4]` | log erg/s        | Bolometric luminosities from different indicators               |
| `logbolo_threshold` | `float32[4]` | log erg/s        | Bolometric luminosity thresholds from different indicators      |

`HDU3` has a table with the names of the four M-sigma relations. The quantites in the
`logmbh` and `logledd` arrays in `HDU2` are in the same order as they are listed in `HDU3`.


| Name                | Data Type    | Unit             | Description                  |
|---------------------|--------------|------------------|------------------------------|
| `mbh_version`       | `str`        |                  | Name of bolometric indicator |

The M-sigma relations are:

* `gultekin09` from Gultekin et al (2009)
* `kormendy13`: from Kormendy et al. (2013)
* `graham13`: from Graham et al. (2013)
* `greene06`: from Greene et al. (2006)

`HDU4` has a table with the names of the four bolometric luminosity estimates
as described in the paper. The quantities in the `logbolo` and `logbolo_threshold`
arrays in `HDU2` are in the same order as they are listed in `HDU4`.


| Name                | Data Type    | Unit             | Description                  |
|---------------------|--------------|------------------|------------------------------|
| `bolo_version`      | `str`        |                  | Name of bolometric indicator |

The bolometric luminosity estimates are:

* `pennell17_nodust_o3` (not used in the paper)
* `heckman04_nodust_o3`: from Heckman et al. (2004)
* `heckman04_dust_o3`: from Kauffmann et al. (2009)
* `netzer19_dust_hb`: from Netzer (2019)

## Luminosity and Eddington Ratio Distributions

We report the cumulative luminosity and Eddington ratio distributions,
both for the raw detected set of AGN and for the selection-corrected
model inference based on a Schechter function, as a function of 
stellar mass and star formation rate. We recommend using a single 
luminosity or Eddington ratio value for each stellar mass from 
these files, not using the function as a whole (since the errors are
highly correlated between luminosity bins).

The `manga-nlr-lum-*` files refer to luminosity distributions, and the
`manga-nlr-edr-*` files refer to Eddington ratio distributions.

The indicator `fagn` refers to a raw detected distribution, whereas
`Fagn` refers to a selection-corrected distribution.

The HDUs are most easily referenced by name, and are

* `FAGN_LAMBDA_CS` or `FAGN_LUM_CS`: the grid of Eddington ratios or
  luminosities, respectively, on which the function is defined. These 
  values are expressed in log-base-10, and the units of the luminosities
  are erg/s.

* `MASS_BIN_CENTERS`: the stellar mass bin centers (in
  log-base-10 of stellar masses).

* `SSFR_BIN_CENTERS`: the specific SFR bin centers (in
  log-base-10 of per-year units).

* `MASS_BIN_EDGES`: the stellar mass bin edges (in
  log-base-10 of stellar masses).

* `SSFR_BIN_EDGES`: the specific SFR bin edges (in
  log-base-10 of per-year units).

* `LOG_FAGN`: the estimate of the cumulative distribution in
  luminosity or Eddington ratio, for each stellar mass bin, sSFR bin,
  and luminosity or Eddington ratio value (i.e. it is an
  `NMASSxNSSFRxNL` array). That is, it is the log-base-10 of the
  fraction of galaxies in each stellar mass and sSFR bin with an AGN
  above each Eddington ratio or luminosity value.

* `LOG_FAGN_STD`: for the inferred model distribution based on the
  Schechter function fits only, the 1-sigma error estimate in dex of
  the `LOG_FAGN` values, based on the standard deviation of values in
  the MCMC chain for the .

* `LOG_FAGN_MEDIAN`: for the inferred model distribution based on
  the Schechter function fits only, the median of the MCMC distribution.

* `LOG_FAGN_LOW`: for the raw detected distributions only, the 1-sigma
  lower limit based on the Clopper-Pearson method for binomial statistics.

* `LOG_FAGN_HIGH`: for the raw detected distributions only, the 1-sigma
  upper limit based on the Clopper-Pearson method for binomial statistics.

