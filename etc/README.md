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

Subdirectories are as follows:

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
| `manga-nlr-edr-[ms]-[lb]-Fagn-0.3.2.fits` | Detected ER distribution         |

Various choices can be made for:

* Line luminosity [lum]: `hb_corr` or `o3`

* M-sigma relation [ms]: `kormendy13`, `gultekin09`, `graham13`, or `greene06`

* Bolometric luminosity [lb]: `netzer19_dust_hb`, `heckman04_nodust_o3`, or `heckman04_dust_o3`

## AGN luminosities and thresholds

The `manga-nlr-agn-0.3.2.fits` files have a single table in `HDU2`:

| Name                               | Data Type | Unit                    | Description                                          |
|------------------------------------|-----------|-------------------------|------------------------------------------------------|
| `plateifu`                         | `str`     |                         | MaNGA Plate-IFU                                      |
| `z`                                | `float32` |                         | Redshift from DRPAll                                 |
| `isagn`                            | `bool`    |                         | Is an AGN?                                           |
| `good`                             | `bool`    |                         | Are measurements good?                               |
| `detected`                         | `bool`    |                         | Are all lines detected?                              |
| `p1`                               | `float32` | dex                     | P1 (Ji and Yan 2020)                                 |
| `p1_err`                           | `float32` | dex                     | Error in P1                                          |
| `p2`                               | `float32` | dex                     | P2 (Ji and Yan 2020)                                 |
| `p2_err`                           | `float32` | dex                     | Error in P2                                          |
| `p3`                               | `float32` | dex                     | P3 (Ji and Yan 2020)                                 |
| `p3_err`                           | `float32` | dex                     | Error in P3                                          |
| `o3`                               | `float32` | $10^{-17}$ erg/cm$ ^2$/s | Flux in [OIII] 5008                                  |
| `o3_err`                           | `float32` | $10^{-17}$ erg/cm^2/s | Error in `o3`                                        |
| `n2ha`                             | `float32` | dex                     | Log of [NII] 6585 / H$\alpha$ ratio                  |
| `n2ha_err`                         | `float32` | dex                     | Error in `n2ha`                                      |
| `s2ha`                             | `float32` | dex                     | Log of [SII] 6718,6732 / H$\alpha$ ratio             |
| `s2ha_err`                         | `float32` | dex                     | Error in `s2ha`                                      |
| `o3hb`                             | `float32` | dex                     | Log of [OIII] 5008 / H$\beta$ ratio                  |
| `o3hb_err`                         | `float32` | dex                     | Error in `o3hb`                                      |
| `hahb`                             | `float32` | dex                     | Log of H$\alpha$ / H$\beta$ ratio                    |
| `hahb_err`                         | `float32` | dex                     | Error in `hahb`                                      |
| `hahb_av`                          | `float32` | mag                     | $A_V$ inferred from `hahb`                           |
| `o3_dust_correction`               | `float32` |                         | Dust correction for [OIII]                           |
| `ha_dust_correction`               | `float32` |                         | Dust correction for H$ \alpha$                        |
| `hb_dust_correction`               | `float32` |                         | Dust correction for H$\beta$                         |
| `log_luminosity_o3`                | `float32` | $log_{10}$(erg/s)       | Log of [OIII] luminosity                             |
| `log_luminosity_o3_threshold`      | `float32` | $log_{10}$(erg/s)       | Log of [OIII] luminosity threshold                   |
| `log_luminosity_o3_corr`           | `float32` | $log_{10}$(erg/s)       | Log of [OIII] dust-corrected luminosity              |
| `log_luminosity_o3_corr_threshold` | `float32` | $log_{10}$(erg/s)       | Log of [OIII] dust-corrected luminosity threshold    |
| `log_luminosity_ha`                | `float32` | $log_{10}$(erg/s)       | Log of H$\alpha$ luminosity                          |
| `log_luminosity_ha_threshold`      | `float32` | $log_{10}$(erg/s)       | Log of H$\alpha$ luminosity threshold                |
| `log_luminosity_ha_corr`           | `float32` | $log_{10}$(erg/s)       | Log of H$\alpha$ dust-corrected luminosity           |
| `log_luminosity_ha_corr_threshold` | `float32` | $log_{10}$(erg/s)       | Log of H$\alpha$ dust-corrected luminosity threshold |
| `log_luminosity_hb`                | `float32` | $log_{10}$(erg/s)       | Log of H$\beta$ luminosity                           |
| `log_luminosity_hb_threshold`      | `float32` | $log_{10}$(erg/s)       | Log of H$\beta$ luminosity threshold                 |
| `log_luminosity_hb_corr`           | `float32` | $log_{10}$(erg/s)       | Log of H$\beta$ dust-corrected luminosity            |
| `log_luminosity_hb_corr_threshold` | `float32` | $log_{10}$(erg/s)       | Log of H$\beta$ dust-corrected luminosity threshold  |

