# Code

This directory contains the data-preparation, statistical-analysis, validation, and figure-generation scripts used in the paper:

> **Reconstructing East Asian Temperatures from 1368 to 1911 Using Historical Documents, Climate Models, and Data Assimilation**

The preferred method for reproducing the analysis is to use the master workflow from the repository root:

```bash
REPRODUCTION_MODE=full \
Rscript --vanilla run_all.R
```

Detailed descriptions of the input data are available in:

```text
Data/README.md
```

Descriptions of the optional upstream data-acquisition and preparation procedures are available in:

```text
Data/Get_data/README.md
```

## Directory Organization

```text
Code/
├── Analysis/
│   ├── coverage_population_GAM.R
│   └── interval_compatibility.R
├── DataPreparation/
│   ├── prepare_lme_annual.R
│   └── prepare_worldclim_bias_correction.R
├── Supplementary/
│   ├── FigureS1.R
│   ├── FigureS2.R
│   └── FigureS3.R
├── Assimilation_grid.R
├── Figure10.R
├── Figure2.R
├── Figure3.R
├── Figure4.R
├── Figure5.R
├── Figure6.R
├── Figure7-8.R
├── Figure7e.R
├── Figure9abc.R
├── Figure9d.R
├── Get_tempe_all_data.R
├── Prior.R
├── README.md
└── prepare_calibration.R
```

All scripts should be executed from the repository root. The scripts use project-relative paths through `here::here()` and should not require manually changing the working directory.

## Master Workflow

The master workflow is defined in:

```text
run_all.R
```

It executes each script in a separate R session, stops when a step fails, verifies key outputs, and records the execution status and runtime of each step.

A full manuscript reproduction is run with:

```bash
REPRODUCTION_MODE=full \
Rscript --vanilla run_all.R
```

A faster smoke test is available with:

```bash
REPRODUCTION_MODE=smoke \
Rscript --vanilla run_all.R
```

Smoke mode uses a reduced penalty-search grid in `Prior.R` and is intended only to test whether the software workflow executes successfully.

**Smoke-mode results should not be used as manuscript results.**

Individual workflow steps can be executed using the `--only` argument. For example:

```bash
Rscript --vanilla run_all.R --only=figure2
```

```bash
Rscript --vanilla run_all.R --only=kriging
```

```bash
Rscript --vanilla run_all.R --only=prior
```

```bash
Rscript --vanilla run_all.R --only=figure9d
```

The execution record is written to:

```text
Output/Logs/run_manifest.csv
```

Individual script logs are written under:

```text
Output/Logs/
```

## Data-Preparation and Core Analysis Scripts

### `DataPreparation/prepare_lme_annual.R`

This script reads the 13 prepared LME ensemble-member files:

```text
Data/LME data/a1.csv.gz
...
Data/LME data/a13.csv.gz
```

It:

- constructs annual temperature fields for 1368–1911;
- converts temperature from Kelvin to degrees Celsius;
- calculates the LME ensemble mean;
- interpolates the LME fields to the nominal GHCN station coordinates for Beijing, Shanghai, and Hong Kong; and
- saves the annual and city-specific LME products.

Important outputs include:

```text
Output/Intermediate/LME/lme_annual_1368_1911.rds
Output/Intermediate/LME/lme_city3_annual_1368_1911.csv
Output/Intermediate/LME/lme_ensemble_mean_1368_1911.csv
```

The optional upstream NetCDF-to-CSV conversion is documented separately in:

```text
Data/Get_data/GetLME.R
```

The original NetCDF files are not required by the default workflow.


### `DataPreparation/prepare_worldclim_bias_correction.R`

This script estimates location-specific climatological bias corrections for the CESM Last Millennium Ensemble (LME) using WorldClim v2.1 as an external reference.

It compares WorldClim v2.1 monthly mean temperatures (1970–2000, 30-arcsecond resolution) with the corresponding climatology from 13 CESM-LME ensemble members.

For each location \(s\), the additive correction is defined as

\[
d(s) = \overline{W}(s) - \overline{L}_{1970:2000}(s),
\]

where \(\overline{W}(s)\) is the WorldClim annual mean temperature and \(\overline{L}_{1970:2000}(s)\) is the corresponding LME climatological mean.

The script:

- extracts modern (1970–2000) temperatures from the original LME NetCDF files;
- applies spatial interpolation consistent with the existing LME preparation procedure;
- extracts WorldClim temperatures at the 266 LME grid locations and three validation stations (Beijing, Shanghai, and Hong Kong);
- computes location-specific additive temperature corrections;
- identifies locations without valid WorldClim values, primarily ocean grid cells;
- provisionally estimates missing corrections using inverse-distance weighting (IDW) from four nearby valid land locations; and
- saves the climatologies, correction offsets, and diagnostic tables.

Among the 89 LME grid cells used in Figure 6 and Supplementary Figure S1, 78 have directly available WorldClim values, while 11 require IDW-based corrections.

For locations without valid WorldClim values, correction offsets are estimated using inverse-distance weighting (IDW) from four nearby land-grid locations. These interpolated corrections should be interpreted with caution, particularly for offshore locations.

This script computes correction offsets only; it does not yet apply them to the historical LME series, quantile mapping, or data assimilation.

#### Required external data

The script requires:

1. The original CESM-LME NetCDF files containing TREFHT for 1970–2000, for all 13 ensemble members.
2. The 12 WorldClim v2.1 monthly temperature GeoTIFF files (`wc2.1_30s_tavg_01.tif` through `wc2.1_30s_tavg_12.tif`).

These large external datasets are not stored in this repository.

The input directories can be configured using the environment variables `LME_NC_DIR` and `WORLDCLIM_DIR`.

#### Execution

Run from the repository root:

```bash
Rscript --vanilla Code/DataPreparation/prepare_worldclim_bias_correction.R
```

#### Outputs

The generated files are saved under `Output/LME_Bias_Correction/`, including:

```text
locations_all_269.csv
worldclim_monthly_all_269.csv
worldclim_masked_locations.csv
lme_annual_1970_2000_all_269.csv
worldclim_offsets_all_269.csv
worldclim_offsets_grid266.csv
worldclim_offsets_city3.csv
```

These outputs are generated locally and are excluded from Git version control.

The script runs independently of the default `run_all.R` workflow. Its correction offsets can be used to adjust the historical LME temperature fields before quantile mapping and data assimilation.


### `prepare_calibration.R`

This script estimates the calibration quantities used by the interval-censored covariance and kriging procedures.

Input:

```text
Data/temperature index value.v1.xlsx
```

Primary output:

```text
Output/Intermediate/calibration_parameters.rds
```

This file contains the calibrated covariance parameters and the plotting objects used by `Figure4.R`.

### `Get_tempe_all_data.R`

This script performs interval-censored kriging of the REACHES documentary temperature index.

Inputs include:

```text
Output/Intermediate/calibration_parameters.rds
```

The script generates:

- predictions on the regular spatial analysis grid;
- predictions at the native LME grid locations;
- direct predictions at the nominal GHCN station coordinates for Beijing, Shanghai, and Hong Kong;
- prediction uncertainties;
- numerical diagnostics; and
- metadata describing the covariance parameters and spatial targets.

Important outputs include:

```text
Output/Intermediate/REACHES/reaches_kriging_lme_grid_mean.csv
Output/Intermediate/REACHES/reaches_kriging_lme_grid_variance.csv
Output/Intermediate/REACHES/reaches_kriging_city3_mean.csv
Output/Intermediate/REACHES/reaches_kriging_city3_sd.csv
Output/Intermediate/REACHES/reaches_kriging_city3_metadata.csv
Output/Intermediate/REACHES/reaches_kriging_metadata.csv
```

The city-specific REACHES predictions are kriged directly at the nominal GHCN station coordinates rather than extracted from the nearest regular-grid cell.

### `Prior.R`

This script estimates the ridge-penalized nonstationary AR(1) parameters used as the dynamic prior in the city-specific assimilation analyses.

Input:

```text
Output/Intermediate/LME/lme_city3_annual_1368_1911.csv
```

It estimates the annual parameter sequences for Beijing, Shanghai, and Hong Kong and writes:

```text
Output/Intermediate/Prior/mtB.csv
Output/Intermediate/Prior/muB.csv
Output/Intermediate/Prior/rtB.csv

Output/Intermediate/Prior/mtS.csv
Output/Intermediate/Prior/muS.csv
Output/Intermediate/Prior/rtS.csv

Output/Intermediate/Prior/mtH.csv
Output/Intermediate/Prior/muH.csv
Output/Intermediate/Prior/rtH.csv
```

The final manuscript workflow uses the complete penalty-search grid under:

```text
REPRODUCTION_MODE=full
```

The reduced grid under `REPRODUCTION_MODE=smoke` is used only for software testing.

### `Assimilation_grid.R`

This script performs the full-grid spatial assimilation used for the supplementary regional summaries.

Inputs include:

```text
Output/Figure6/Figure6_cluster_assignments.csv
Output/Intermediate/LME/lme_annual_1368_1911.rds
Output/Intermediate/REACHES/reaches_kriging_lme_grid_mean.csv
Output/Intermediate/REACHES/reaches_kriging_lme_grid_variance.csv
Output/Intermediate/REACHES/reaches_kriging_metadata.csv
```

Important outputs include:

```text
Output/Intermediate/Assimilation/assimilated_posterior_lme_grid_mean.csv
Output/Intermediate/Assimilation/assimilated_posterior_lme_grid_variance.csv
Output/Intermediate/Assimilation/assimilated_posterior_lme_grid_diagnostics.csv
Output/Intermediate/Assimilation/assimilated_posterior_lme_grid_metadata.csv
```

The full-grid assimilation uses the fixed spatial penalty values documented in the script.

## Main-Manuscript Figure Scripts

### `Figure2.R`

Input:

```text
Data/temperature index value.v1.xlsx
```

Outputs:

```text
Output/Figure2/Figure2(a).png
Output/Figure2/Figure2(b).png
Output/Figure2/Figure2(c).png
```

The panels show:

- annual counts of REACHES documentary records;
- empirical frequencies of documentary temperature levels; and
- the spatial distribution of REACHES records.

### `Figure3.R`

Inputs:

```text
Output/Intermediate/LME/lme_city3_annual_1368_1911.csv
```

Outputs:

```text
Output/Figure3/Figure3(a).png
Output/Figure3/Figure3(b).png
Output/Figure3/Figure3(c).png
```

These panels provide the city-specific REACHES and LME comparisons for Beijing, Shanghai, and Hong Kong.

### `Figure4.R`

Input:

```text
Output/Intermediate/calibration_parameters.rds
```

Outputs:

```text
Output/Figure4/Figure4(a).png
Output/Figure4/Figure4(b).png
Output/Figure4/Figure4(c).png
```

`prepare_calibration.R` must be completed before running `Figure4.R`.

### `Figure5.R`

Inputs include:

```text
Output/Intermediate/calibration_parameters.rds
```

Outputs:

```text
Output/Figure5/Figure5(a1).png
Output/Figure5/Figure5(a2).png
Output/Figure5/Figure5(a3).png
Output/Figure5/Figure5(b1).png
Output/Figure5/Figure5(b2).png
Output/Figure5/Figure5(b3).png
```

`prepare_calibration.R` must be completed before running `Figure5.R`.

### `Figure6.R`

Input:

```text
Output/Intermediate/REACHES/reaches_kriging_lme_grid_mean.csv
```

Outputs:

```text
Output/Figure6/Figure6_cluster_map.jpg
Output/Figure6/Figure6_cluster_1.jpg
Output/Figure6/Figure6_cluster_2.jpg
Output/Figure6/Figure6_cluster_3.jpg
Output/Figure6/Figure6_cluster_4.jpg
Output/Figure6/Figure6_cluster_assignments.csv
```

The current analysis selects four clusters using the stated model-selection procedure.

The cluster assignments saved by this script are also used by `Supplementary/FigureS1.R`.

### `Figure7-8.R`

Inputs:

```text
Output/Intermediate/LME/lme_city3_annual_1368_1911.csv
Output/Intermediate/REACHES/reaches_kriging_city3_mean.csv
Output/Intermediate/REACHES/reaches_kriging_city3_sd.csv
```

Outputs:

```text
Output/Figure7-8/Figure7a.png
Output/Figure7-8/Figure7b.png
Output/Figure7-8/Figure7c.png
Output/Figure7-8/Figure7d.png

Output/Figure7-8/Figure8a.png
Output/Figure7-8/Figure8b.png
Output/Figure7-8/Figure8c.png
Output/Figure7-8/Figure8d.png
```

The script performs the city-specific uncertainty-aware quantile-mapping calculations and produces the corresponding Figure 7 and Figure 8 panels.

### `Figure7e.R`

Inputs:

```text
Output/Intermediate/REACHES/reaches_kriging_lme_grid_mean.csv
Output/Intermediate/REACHES/reaches_kriging_lme_grid_variance.csv
Output/Intermediate/LME/lme_annual_1368_1911.rds
```

Output:

```text
Output/Figure7-8/Figure7e.png
```

This script produces the full-grid quantile-mapping result shown in Figure 7(e).

### `Figure9abc.R`

Inputs:

```text
Output/Intermediate/Prior/mtB.csv
Output/Intermediate/Prior/muB.csv
Output/Intermediate/Prior/rtB.csv

Output/Intermediate/Prior/mtS.csv
Output/Intermediate/Prior/muS.csv
Output/Intermediate/Prior/rtS.csv

Output/Intermediate/Prior/mtH.csv
Output/Intermediate/Prior/muH.csv
Output/Intermediate/Prior/rtH.csv
```

Main-manuscript outputs:

```text
Output/Figure9/Figure9a.png
Output/Figure9/Figure9b.png
Output/Figure9/Figure9c.png
```

Supplementary outputs:

```text
Output/Supplementary/FigureS4a.png
Output/Supplementary/FigureS4b.png
Output/Supplementary/FigureS4c.png

Output/Supplementary/FigureS5a.png
Output/Supplementary/FigureS5b.png
Output/Supplementary/FigureS5c.png
```

Additional plot data are saved in:

```text
Output/Figure9/Figure9abc_parameter_plot_data.csv
```

### `Figure9d.R`

Inputs include:

```text
Output/Intermediate/LME/lme_city3_annual_1368_1911.csv
Output/Intermediate/REACHES/reaches_kriging_city3_mean.csv
Output/Intermediate/REACHES/reaches_kriging_city3_sd.csv
Output/Intermediate/REACHES/reaches_kriging_metadata.csv
Output/Intermediate/Prior/
```

Main-manuscript output:

```text
Output/Figure9/Figure9d.png
```

Supplementary outputs:

```text
Output/Supplementary/FigureS4d.png
Output/Supplementary/FigureS5d.png
```

The script runs the annual Kalman filter and Rauch–Tung–Striebel smoother for all three cities.

Years without documentary information are retained in the annual state sequence and are handled through prediction-only Kalman steps.

The script also generates:

```text
Data/Valid/assimilation_Beijing.csv
Data/Valid/assimilation_Shanghai.csv
Data/Valid/assimilation_HongKong.csv
```

and diagnostic outputs under:

```text
Output/Intermediate/Assimilation/
```

The diagnostics include:

- calibrated REACHES temperature;
- effective proxy temperature;
- dynamic prior prediction;
- filtered temperature;
- smoothed temperature;
- LME ensemble-mean temperature;
- measurement-equation parameters;
- Kalman observation weights;
- filter-identity checks;
- bias and RMSE summaries; and
- anomaly-based comparisons.

### `Figure10.R`

Inputs:

```text
Data/Valid/assimilation_Beijing.csv
Data/Valid/assimilation_Shanghai.csv
Data/Valid/assimilation_HongKong.csv
Data/GHCNv4.xlsx
```

Outputs:

```text
Output/Figure10/Figure10a.png
Output/Figure10/Figure10b.png
Output/Figure10/Figure10c.png
Output/Figure10/Figure10_metrics.csv
```

The script compares the city-specific reconstructed temperatures with the GHCN instrumental records and reports:

- correlation;
- mean bias;
- root mean-square error;
- anomaly correlation; and
- anomaly root mean-square error.

## Supplementary Figure Scripts

### `Supplementary/FigureS1.R`

Inputs include:

```text
Output/Figure6/Figure6_cluster_assignments.csv
Output/Intermediate/REACHES/reaches_kriging_lme_grid_mean.csv
Output/Intermediate/LME/lme_ensemble_mean_1368_1911.csv
Output/Intermediate/Assimilation/assimilated_posterior_lme_grid_mean.csv
```

The script produces twelve functional-boxplot panels:

```text
Output/Supplementary/FigureS1_panels/FigureS1_reaches_cluster_1.png
Output/Supplementary/FigureS1_panels/FigureS1_reaches_cluster_2.png
Output/Supplementary/FigureS1_panels/FigureS1_reaches_cluster_3.png
Output/Supplementary/FigureS1_panels/FigureS1_reaches_cluster_4.png

Output/Supplementary/FigureS1_panels/FigureS1_lme_cluster_1.png
Output/Supplementary/FigureS1_panels/FigureS1_lme_cluster_2.png
Output/Supplementary/FigureS1_panels/FigureS1_lme_cluster_3.png
Output/Supplementary/FigureS1_panels/FigureS1_lme_cluster_4.png

Output/Supplementary/FigureS1_panels/FigureS1_posterior_cluster_1.png
Output/Supplementary/FigureS1_panels/FigureS1_posterior_cluster_2.png
Output/Supplementary/FigureS1_panels/FigureS1_posterior_cluster_3.png
Output/Supplementary/FigureS1_panels/FigureS1_posterior_cluster_4.png
```

Additional outputs include:

```text
Output/Intermediate/FigureS1/FigureS1_cluster_assignments.csv
Output/Intermediate/FigureS1/FigureS1_diagnostics.csv
Output/Intermediate/FigureS1/FigureS1_panel_files.csv
```

Figure S1 uses the same four-cluster assignments generated by `Figure6.R`.

### `Supplementary/FigureS2.R`

Input:

```text
Data/temperature index value.v1.xlsx
```

Outputs:

```text
Output/Supplementary/FigureS2a.png
Output/Supplementary/FigureS2b.png
Output/Supplementary/FigureS2c.png
Output/Supplementary/FigureS2d.png
Output/Supplementary/FigureS2e.png
```

### `Supplementary/FigureS3.R`

Inputs:

```text
Data/Valid/assimilation_Beijing.csv
Data/Valid/assimilation_Shanghai.csv
Data/Valid/assimilation_HongKong.csv
```

Outputs:

```text
Output/Supplementary/FigureS3a.png
Output/Supplementary/FigureS3b.png
Output/Supplementary/FigureS3c.png
Output/Supplementary/FigureS3_plot_data.csv
```



## Additional Analysis and Table Scripts

### `Analysis/coverage_population_GAM.R`

This script performs the documentary-coverage and historical-population analysis described in Section 5.1.

Inputs include:

```text
Data/temperature index value.v1.xlsx
Output/Intermediate/LME/lme_annual_1368_1911.rds
Data/Population/
```

Important outputs include:

```text
Output/Tables/Table2_coverage_GAM.csv
Output/Tables/Table2_Hgu_coefficient.csv
Output/Tables/TableS2_M4_parametric_coefficients.csv
Output/Tables/TableS2_M4_parametric_coefficients_formatted.csv
```

The script fits the baseline, climate-adjusted, population-adjusted, and full population-adjusted climate models.

### `Analysis/interval_compatibility.R`

This script performs the held-out interval-compatibility validation described in Section 5.2.

Inputs include:

```text
Data/temperature index value.v1.xlsx
```

Outputs:

```text
Output/Validation/interval_compatibility_summary.csv
Output/Validation/interval_compatibility_summary_formatted.csv
```

The validation procedure:

- holds out complete site-year observations;
- prevents held-out records from entering model fitting;
- uses great-circle spatial distances;
- evaluates nominal compatibility levels;
- calculates ranked probability scores;
- calculates classification accuracy; and
- calculates negative log scores.

## Recommended Execution Order

The recommended order is encoded in `run_all.R`.

The WorldClim-based LME bias-correction procedure can be executed separately using `DataPreparation/prepare_worldclim_bias_correction.R`. It requires the original CESM-LME NetCDF files and WorldClim GeoTIFF files and is not currently part of the default `run_all.R` workflow.


The main dependency sequence is:

```text
prepare_lme_annual.R
        ↓
prepare_calibration.R
        ↓
Get_tempe_all_data.R
        ↓
Figure and validation scripts
        ↓
Prior.R
        ↓
Figure9abc.R
        ↓
Figure9d.R
        ↓
Figure10.R and FigureS3.R
```

The spatial supplementary workflow follows:

```text
Figure6.R
        ↓
Assimilation_grid.R
        ↓
Supplementary/FigureS1.R
```

The complete execution order, logs, and output checks are managed automatically by:

```text
run_all.R
```

## Generated Files

Generated files are written under:

```text
Output/
```

The `Output/` directory is excluded from version control and is created automatically when the workflow is run.

Reference copies of the final manuscript and supplementary figures are provided separately under:

```text
Figure/
```

The generated outputs can be compared with the corresponding reference copies when checking the reproduced results.

The relationship between manuscript items, scripts, inputs, and outputs is documented in:

```text
output_manifest.csv
```

at the repository root.
