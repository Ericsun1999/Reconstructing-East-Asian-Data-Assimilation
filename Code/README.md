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
│   └── prepare_lme_annual.R
├── Supplementary/
│   ├── FigureS1.R
│   ├── FigureS2.R
│   └── FigureS5.R
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

The cluster assignments saved by this script are also used by `Supplementary/FigureS5.R`.

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
Output/Supplementary/FigureS3a.png
Output/Supplementary/FigureS3b.png
Output/Supplementary/FigureS3c.png

Output/Supplementary/FigureS4a.png
Output/Supplementary/FigureS4b.png
Output/Supplementary/FigureS4c.png
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
Output/Supplementary/FigureS3d.png
Output/Supplementary/FigureS4d.png
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

Input:

```text
Data/temperature index value.v1.xlsx
```

Outputs:

```text
Output/Supplementary/FigureS1a.png
Output/Supplementary/FigureS1b.png
Output/Supplementary/FigureS1c.png
Output/Supplementary/FigureS1d.png
Output/Supplementary/FigureS1e.png
```

### `Supplementary/FigureS2.R`

Inputs:

```text
Data/Valid/assimilation_Beijing.csv
Data/Valid/assimilation_Shanghai.csv
Data/Valid/assimilation_HongKong.csv
```

Outputs:

```text
Output/Supplementary/FigureS2a.png
Output/Supplementary/FigureS2b.png
Output/Supplementary/FigureS2c.png
Output/Supplementary/FigureS2_plot_data.csv
```

### `Supplementary/FigureS5.R`

Inputs include:

```text
Output/Figure6/Figure6_cluster_assignments.csv
Output/Intermediate/REACHES/reaches_kriging_lme_grid_mean.csv
Output/Intermediate/LME/lme_ensemble_mean_1368_1911.csv
Output/Intermediate/Assimilation/assimilated_posterior_lme_grid_mean.csv
```

The script produces twelve functional-boxplot panels:

```text
Output/Supplementary/FigureS5_panels/FigureS5_reaches_cluster_1.png
Output/Supplementary/FigureS5_panels/FigureS5_reaches_cluster_2.png
Output/Supplementary/FigureS5_panels/FigureS5_reaches_cluster_3.png
Output/Supplementary/FigureS5_panels/FigureS5_reaches_cluster_4.png

Output/Supplementary/FigureS5_panels/FigureS5_lme_cluster_1.png
Output/Supplementary/FigureS5_panels/FigureS5_lme_cluster_2.png
Output/Supplementary/FigureS5_panels/FigureS5_lme_cluster_3.png
Output/Supplementary/FigureS5_panels/FigureS5_lme_cluster_4.png

Output/Supplementary/FigureS5_panels/FigureS5_posterior_cluster_1.png
Output/Supplementary/FigureS5_panels/FigureS5_posterior_cluster_2.png
Output/Supplementary/FigureS5_panels/FigureS5_posterior_cluster_3.png
Output/Supplementary/FigureS5_panels/FigureS5_posterior_cluster_4.png
```

Additional outputs include:

```text
Output/Intermediate/FigureS5/FigureS5_cluster_assignments.csv
Output/Intermediate/FigureS5/FigureS5_diagnostics.csv
Output/Intermediate/FigureS5/FigureS5_panel_files.csv
```

Figure S5 uses the same four-cluster assignments generated by `Figure6.R`.

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
Figure10.R and FigureS2.R
```

The spatial supplementary workflow follows:

```text
Figure6.R
        ↓
Assimilation_grid.R
        ↓
Supplementary/FigureS5.R
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
