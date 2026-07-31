# Upstream Data Acquisition and Preparation

This directory contains documentation and optional scripts for obtaining and preparing the original data sources used in the study.

The default reproducibility workflow does **not** begin from all original downloaded archives. Instead, it begins from the prepared input files distributed under `Data/`.

The main workflow is executed from the repository root using:

```bash
REPRODUCTION_MODE=full \
Rscript --vanilla run_all.R
```

## 1. REACHES Documentary Temperature Data

The REACHES dataset is reconstructed from historical documentary records and contains information including:

- documentary temperature level;
- calendar year;
- longitude; and
- latitude.

Source:

[NOAA Paleoclimatology Study 37720](https://www.ncei.noaa.gov/access/paleo-search/study/37720)

The prepared REACHES workbook used by the analysis is:

```text
Data/temperature index value.v1.xlsx
```

This file is used by the scripts for:

- descriptive analysis;
- interval-censored covariance calibration;
- kriging;
- functional clustering;
- documentary-coverage analysis; and
- held-out internal validation.

The current reproducibility workflow begins from the prepared Excel workbook above. Any cleaning, variable selection, or formatting applied to the originally downloaded REACHES files should be documented here if the prepared workbook is regenerated.

## 2. Last Millennium Ensemble Data

The CESM Last Millennium Ensemble contains simulated climate fields indexed by:

- longitude;
- latitude;
- ensemble member;
- calendar month; and
- climate variable.

Source:

[CESM Last Millennium Ensemble](https://www.cesm.ucar.edu/community-projects/lme)

The analysis uses the near-surface air-temperature variable:

```text
TREFHT
```

The original temperature values are stored in Kelvin.

### Original NetCDF preparation

The optional script:

```text
Data/Get_data/GetLME.R
```

documents the conversion from the original LME NetCDF files to the prepared member-level files used by the public workflow.

The original NetCDF files are not included in this repository because of their size. Consequently, `GetLME.R` is not called by the default `run_all.R` workflow.

To reproduce the NetCDF-to-CSV preparation, users must:

1. obtain the required LME NetCDF files from the CESM source;
2. place them in the input location expected by `GetLME.R`;
3. verify that the filenames and time ranges match those specified in the script;
4. run `GetLME.R`; and
5. save the resulting prepared member files under:

```text
Data/LME data/
```

The distributed reproducibility workflow begins with:

```text
Data/LME data/a1.csv.gz
Data/LME data/a2.csv.gz
Data/LME data/a3.csv.gz
Data/LME data/a4.csv.gz
Data/LME data/a5.csv.gz
Data/LME data/a6.csv.gz
Data/LME data/a7.csv.gz
Data/LME data/a8.csv.gz
Data/LME data/a9.csv.gz
Data/LME data/a10.csv.gz
Data/LME data/a11.csv.gz
Data/LME data/a12.csv.gz
Data/LME data/a13.csv.gz
```

These files are processed by:

```text
Code/DataPreparation/prepare_lme_annual.R
```

That script:

- reads the 13 prepared ensemble-member files;
- converts monthly values to annual temperature fields;
- constructs the 1368–1911 analysis period;
- converts temperature from Kelvin to degrees Celsius;
- calculates the ensemble mean; and
- spatially interpolates the LME fields to the nominal GHCN station coordinates used in the city-specific analyses.

For full upstream reproducibility, this documentation should also record:

- the exact NetCDF filenames;
- the LME dataset version or access date;
- ensemble-member identifiers;
- temporal coverage of each file;
- the original longitude and latitude conventions; and
- file checksums where available.

## 3. GHCN Instrumental Temperature Data

The instrumental validation data are obtained from the NOAA Global Historical Climatology Network Monthly dataset.

Source directory:

[NOAA GHCN-Monthly Version 4](https://www.ncei.noaa.gov/pub/data/ghcn/v4/)

The quality-controlled average-temperature archive used as the upstream source is:

```text
ghcnm.tavg.latest.qcf.tar.gz
```

### Archive extraction

The archive can be extracted using standard decompression tools. For example:

```bash
tar -xzf ghcnm.tavg.latest.qcf.tar.gz
```

The extracted GHCN files are not necessarily ordinary comma-separated CSV files. They follow the GHCN-Monthly distribution format and require parsing, station selection, missing-value handling, temperature scaling, and reshaping before they can be used by the analysis.

The prepared file used by the current workflow is:

```text
Data/GHCNv4.xlsx
```

The GHCN preparation procedure should document:

- the exact GHCN archive version or download date;
- the selected station identifiers;
- station names;
- station longitude and latitude;
- the temperature variable used;
- the original temperature scaling;
- missing-value codes;
- quality-control flags;
- monthly-to-annual aggregation rules;
- minimum number of available months required for an annual mean;
- the treatment of duplicate or replacement station records; and
- the exact procedure used to produce `GHCNv4.xlsx`.

The city-specific analyses use the following nominal GHCN station coordinates:

| City | Longitude | Latitude |
|---|---:|---:|
| Beijing | 116.283 | 39.933 |
| Shanghai | 121.433 | 31.167 |
| Hong Kong | 114.167 | 22.333 |

The REACHES and LME city-specific series are evaluated at these same nominal coordinates.

## 4. Historical Gridded Population Data

The historical population data are obtained from the **Gridded Population Dataset in the Traditional Cultivated Region of China from 1776 to 1953**.

Source:

[CASEarth Data Sharing and Service Portal](https://data.casearth.cn/dataset/6538b2b4819aec0f262199b4)

The CASEarth data portal is primarily presented in Chinese. Users may need to register for an account or log in before downloading the data.

The dataset was developed by Zhang et al. (2022) using prefecture-level historical population estimates examined and corrected by historians. A random forest population-distribution model was used to allocate the prefecture-level population totals to a 10 km by 10 km spatial grid based on terrain, climate, river, and city-related environmental factors.

The dataset contains population-density rasters for six historical time slices:

```text
1776_pd.tif
1820_pd.tif
1851_pd.tif
1880_pd.tif
1910_pd.tif
1953_pd.tif
```

The raster values represent:

`population density in persons per square kilometre`

The geographic coverage corresponds to the traditional cultivated region of China, consisting of 18 historical provinces and excluding Taiwan Prefecture.


## 4. Historical Gridded Population Data

The historical population data are obtained from the **Gridded Population Dataset in the Traditional Cultivated Region of China from 1776 to 1953**.

Source:

[CASEarth Data Sharing and Service Portal](https://data.casearth.cn/dataset/6538b2b4819aec0f262199b4)

The CASEarth dataset page is presented in Chinese. The original population raster files can be downloaded directly from the dataset page by following the portal's download instructions.

The dataset was developed by Zhang et al. (2022) using prefecture-level historical population estimates examined and corrected by historians. A random forest population-distribution model was used to allocate prefecture-level population totals to a 10 km by 10 km spatial grid based on terrain, climate, river, and city-related environmental factors.

The dataset covers the traditional cultivated region of China and provides gridded population-density reconstructions for six historical time slices:


```text
1776
1820
1851
1880
1910
1953
```

The downloaded dataset contains the following GeoTIFF files:

```text
1776_pd.tif
1820_pd.tif
1851_pd.tif
1880_pd.tif
1910_pd.tif
1953_pd.tif
```

The raster values represent population density in persons per square kilometre.

The geographic coverage corresponds to the traditional cultivated region of China, consisting of 18 historical provinces and excluding Taiwan Prefecture.



## Reproducibility Boundary

The default public workflow begins from the prepared files distributed under `Data/`, including:

```text
Data/temperature index value.v1.xlsx
Data/GHCNv4.xlsx
Data/LME data/a1.csv.gz
...
Data/LME data/a13.csv.gz
```

The scripts and documentation under `Data/Get_data/` describe optional upstream acquisition and preparation steps.

Because some original archives are not redistributed, users reproducing the entire process from the original data sources must first obtain those files from the cited providers and follow the documented preprocessing steps.
