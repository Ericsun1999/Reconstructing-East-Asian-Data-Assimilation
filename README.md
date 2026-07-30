# Reconstructing-East-Asian-Data-Assimilation 

This repository contains the code and prepared data used to reproduce the computational results in the paper:

> **Reconstructing East Asian Temperatures from 1368 to 1911 Using Historical Documents, Climate Models, and Data Assimilation**  
> Eric Sun, Kuan-hui Elaine Lin, Wan-Ling Tseng, Pao K. Wang, and Hsin-Cheng Huang  
> [arXiv preprint](http://arxiv.org/abs/2410.21790)

The repository includes implementations of:

- interval-censored kriging for the REACHES documentary temperature index;
- uncertainty-aware quantile mapping;
- ridge-penalized nonstationary AR(1) estimation;
- Kalman filtering and Rauch–Tung–Striebel smoothing;
- functional clustering and spatial regionalization;
- documentary-coverage and historical-population analysis; and
- held-out interval-compatibility validation.

The workflow generates Figures 2–10 in the main manuscript, Figures S1–S5 in the Supplementary Material, and the associated computational tables and diagnostic outputs.

## Repository organization

```text
.
├── Code/                   Analysis and figure-generation scripts
├── Data/                   Input data and documented prepared datasets
├── Figure/                 Reference copies of Figures 2–10 and Figures S1–S5
├── Output/                 Generated outputs; created by the workflow
├── renv/                   Project environment infrastructure
├── .Rprofile               Activates the project-specific renv environment
├── renv.lock               Locked R package environment
├── run_all.R               Master reproducibility workflow
├── output_manifest.csv     Mapping from manuscript items to code and outputs
└── README.md
```

Detailed descriptions of the analysis scripts and their usage are provided in:

```text
Code/README.md
```

Descriptions, provenance, formats, and usage instructions for the input data are provided in:

```text
Data/README.md
```

The `Figure/` directory contains reference copies of the manuscript and supplementary figures so that the reported figures can be viewed without first running the computational workflow.

The `Output/` directory is created by the workflow and does not need to exist before the analysis is run. After executing `run_all.R`, newly generated figures, tables, intermediate files, diagnostics, and execution logs are saved under `Output/`.


## Software Environment

The analysis is implemented in R.

The package versions used by the project are recorded in:

```text
renv.lock
```

From the repository root, install `renv` when necessary:

```bash
Rscript -e 'install.packages("renv", repos = "https://cloud.r-project.org")'
```

Restore the project-specific R environment with:

```bash
Rscript -e 'renv::restore(prompt = FALSE)'
```

The project `.Rprofile` activates the local `renv` environment automatically when the repository is opened as an R project.

The R version, operating system, and loaded package versions for a completed run can be inspected using:

```r
sessionInfo()
```

## Reproducing the Complete Analysis

All commands should be executed from the repository root, which contains:

```text
run_all.R
Code/
Data/
renv.lock
```

Run the complete manuscript workflow with:

```bash
REPRODUCTION_MODE=full \
Rscript --vanilla run_all.R
```

The full mode uses the complete penalty-parameter search grid used for the manuscript analysis.

The workflow:

1. checks the required input files;
2. executes each analysis script in a separate R session;
3. stops immediately when a step fails;
4. checks that key output files were generated;
5. records the execution time and status of each step; and
6. stores individual log files for all workflow stages.

Execution logs are saved under:

```text
Output/Logs/
```

The main workflow execution record is:

```text
Output/Logs/run_manifest.csv
```

This file reports:

- workflow step;
- generating script;
- start and end times;
- elapsed execution time;
- exit status;
- success status;
- missing expected outputs; and
- corresponding log file.

## Optional Smoke Test

A faster smoke-test mode is available for checking whether the complete software workflow executes successfully:

```bash
REPRODUCTION_MODE=smoke \
Rscript --vanilla run_all.R
```

Smoke mode uses a reduced penalty-parameter search grid in `Prior.R` to shorten execution time.

Smoke mode is intended only for software and workflow testing.

**Smoke-mode results should not be used as the reproduced manuscript results. All manuscript results must be generated using `REPRODUCTION_MODE=full`.**

## Running Individual Workflow Steps

Individual stages can be tested separately using the `--only` argument.

For example:

```bash
Rscript --vanilla run_all.R --only=figure2
```

```bash
Rscript --vanilla run_all.R --only=calibration
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

Available workflow identifiers are defined in the `workflow` object in:

```text
run_all.R
```

## LME Input Boundary

The distributed reproducibility workflow begins with the 13 prepared, gzip-compressed Last Millennium Ensemble member files:


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

These files are the inputs to:

```text
Code/DataPreparation/prepare_lme_annual.R
```

The optional script:

```text
Data/Get_data/GetLME.R
```

documents the upstream conversion from the original NetCDF files to the prepared member-level CSV files.

The original NetCDF files are not redistributed in this repository because of their size. Consequently, `Data/Get_data/GetLME.R` is not called by the default `run_all.R` workflow.

The data documentation describes:

- the original LME data source;
- expected NetCDF filenames;
- ensemble-member structure;
- temperature variable;
- temporal coverage;
- spatial coordinates;
- monthly-to-annual conversion; and
- conversion from Kelvin to degrees Celsius.

## Generated Outputs

The workflow creates outputs under:

```text
Output/
```

Major output directories include:


```text
Output/Figure2/
Output/Figure3/
Output/Figure4/
Output/Figure5/
Output/Figure6/
Output/Figure7-8/
Output/Figure9/
Output/Figure10/
Output/Supplementary/
Output/Tables/
Output/Validation/
Output/Intermediate/
Output/Logs/
```

The file:

```text
output_manifest.csv
```

maps each manuscript or supplementary figure, table, and essential intermediate product to:

- its generating script;
- its primary inputs;
- its output files;
- its workflow step;
- its output type;
- whether it is required for the manuscript;
- whether it was successfully generated; and
- whether it was verified against the revised manuscript.

## Main and Supplementary Figures

The provided scripts generate:

- Figures 2–10 in the main manuscript;
- Figures S1–S5 in the Supplementary Material;
- city-specific assimilation diagnostics;
- validation summaries;
- documentary-coverage and population-model tables; and
- intermediate products required to reproduce the final analyses.

The revised functional clustering analysis uses four clusters selected by the stated model-selection procedure.

Figure S5 reuses the same cluster assignments generated for Figure 6.

## City-Specific Spatial Locations

The city-specific REACHES and LME analyses use the same nominal GHCN station coordinates:

| City | Longitude | Latitude |
|---|---:|---:|
| Beijing | 116.283 | 39.933 |
| Shanghai | 121.433 | 31.167 |
| Hong Kong | 114.167 | 22.333 |

REACHES predictions are kriged directly at these coordinates.

The LME temperature fields are spatially interpolated to the same nominal coordinates.

The use of common nominal coordinates ensures that the city-specific REACHES, LME, assimilation, and GHCN validation series are evaluated at consistent locations.

## Reproducibility Records

After a complete run, the following files provide the main computational record:

```text
Output/Logs/run_manifest.csv
output_manifest.csv
renv.lock
```


The individual execution logs are stored under:

```text
Output/Logs/
```

Generated figures and numerical results should be checked against the corresponding items in the final manuscript before creating a versioned release.

## Full and Smoke Modes

The workflow supports two execution modes.

### Full mode

```bash
REPRODUCTION_MODE=full \
Rscript --vanilla run_all.R
```

Full mode:

- uses the complete penalty search grid;
- reproduces the final manuscript analysis;
- should be used for all reported figures and numerical results; and
- may require substantially more computation time.

### Smoke mode


```bash
REPRODUCTION_MODE=smoke \
Rscript --vanilla run_all.R
```

Smoke mode:

- uses a reduced penalty search grid;
- checks whether the end-to-end workflow executes;
- is useful for identifying missing files, packages, and hidden dependencies; and
- does not reproduce the final manuscript model-selection procedure.

## Reproducibility Status

The complete workflow has been tested locally from the repository root under the project-specific `renv` environment.

The full workflow should be considered successful when every row of:

```text
Output/Logs/run_manifest.csv
```

satisfies:

```text
exit_status = 0
success = TRUE
missing_outputs = empty
```

A successful smoke-mode run verifies software integration.

A successful full-mode run verifies the complete local scientific workflow using the manuscript penalty-search settings.

## Citation

A version-specific software citation and archived DOI will be added when the final reproducibility release is created.




