# Code

This directory contains the data-preparation, statistical-analysis, validation, and figure-generation scripts used in the paper:

> **Reconstructing East Asian Temperatures from 1368 to 1911 Using Historical Documents, Climate Models, and Data Assimilation**

For the current two-version manuscript reconstruction, the
recommended entry point is the versioned controller from the
repository root:

```bash
REPRODUCTION_MODE=full \
Rscript --vanilla run_all.R --all
```

The independent LME and dynamic prior are estimated once and reused,
while the inferred-zero primary figures and Original sensitivity
figures remain in different output workspaces. Legacy single-step
execution is available through `run_all.R --only=<step>`.

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
├── Analysis/                      Manuscript tables, validation and comparisons
│   ├── coverage_population_GAM.R
│   ├── interval_compatibility.R
│   ├── snapshot_augmentation_reconstruction.R
│   └── compare_augmentation_reconstructions.R
├── DataPreparation/               All source-data preparation
│   ├── prepare_lme_annual.R
│   ├── prepare_worldclim_bias_correction.R
│   └── prepare_precipitation_augmented_reaches.R
├── SupportingAnalyses/            Optional diagnostics and exploration
│   ├── README.md                  Explanations and optional commands
│   └── (five supporting .R scripts)
├── Workflow/
│   └── run_single_version.R       Internal worker called by run_all.R
├── Supplementary/                 Supplementary figures S1–S3 scripts
├── Figure2.R ... Figure10.R        Main figure scripts
├── Figure9d.R                     Measurement variance helper definitions
├── cv_observation_variance.R      Required Original-input blocked-CV estimator
├── Get_tempe_all_data.R
├── Prior.R
└── prepare_calibration.R
```

All commands below are run from the repository root. The five
SupportingAnalyses scripts are optional diagnostics or exploratory
methods and are not required by the full reproduction; see
[their documentation](SupportingAnalyses/README.md). The Original-input
observation-variance CV estimator is a required upstream step and is
called automatically from the root `run_all.R`.

## Versioned manuscript reconstruction (primary and sensitivity)

Use the repository-root dual-version `run_all.R` for manuscript
reproduction, rather than two sequential unversioned runs. The **primary**
reconstruction includes the documented site-years with inferred zero
temperature indices; **Original** uses only the recorded temperature
indices and serves as the sensitivity reference. Inferred zeros are
distinguished from observed temperature categories in the input
provenance and descriptive maps.

From the repository root:

~~~bash
git pull
REPRODUCTION_MODE=full Rscript --vanilla run_all.R --all
~~~

This executes the shared 13-member LME annual preparation, WorldClim
bias correction and `Prior.R` **once**. Both documentary-input
workspaces then reuse the exact same prepared LME and prior files
while separately recalibrating REACHES, recomputing complete spatial
kriging, and regenerating Figures 2–10, Supplementary Figures S1–S5,
associated tables, diagnostics and Kalman posteriors.

Outputs are **physically separate, including intermediate data**:

~~~text
Output/Reconstructions/
  shared/Output/Intermediate/LME/              one common LME
  shared/Output/Intermediate/Prior/            one common prior
  inferred-zero/Output/Figure2/ ... Figure10/   PRIMARY
  inferred-zero/Output/Supplementary/           PRIMARY
  inferred-zero/Output/Intermediate/            PRIMARY
  inferred-zero/Output/Logs/                    PRIMARY
  original/Output/Figure2/ ... Figure10/        SENSITIVITY
  original/Output/Supplementary/                SENSITIVITY
  original/Output/Intermediate/                 SENSITIVITY
  original/Output/Logs/                         SENSITIVITY
  Comparison/city/
  Comparison/grid/
~~~

The controller isolates each run in a project workspace. Its Data/Valid
folder is separate from tracked reference data, and the output
workspaces are not overwritten by switching source variants.
The two versions use the original measurement projection, the
same Monte Carlo settings and a **common fixed extra observation
variance** estimated by the Original-input five-block CV
(`Code/cv_observation_variance.R`). The controller fits this value
automatically after the shared LME/prior are ready, or reuses the estimate
if source and workflow fingerprints match. The generated fit is stored at
`Output/Reconstructions/original/Output/Intermediate/ObservationVariance/original/selected_u.csv`.
For deliberate sensitivity checks, use `--u=VALUE` to bypass CV and
supply an independently justified nonnegative constant.
This isolates the effect of the inferred-zero input; the Original
fit is not misrepresented as an independently fitted Augmented u.
The one-time LME/prior source files are fingerprinted to prevent
accidental reuse after the corresponding inputs or code change.

To perform stages separately:

~~~bash
Rscript --vanilla run_all.R --prepare-shared
Rscript --vanilla run_all.R --variant=inferred-zero
Rscript --vanilla run_all.R --variant=original
Rscript --vanilla run_all.R --compare-only
~~~

Set `REPRODUCTION_MODE=full` for manuscript results; smoke mode is only
a pipeline test and uses `u=0` by default to avoid the costly five-block
CV. In full mode, no separate CV command is required.
Use `--rebuild-shared` deliberately if the independent LME/prior inputs
or code change, and then rerun *both* documentary versions. Existing
snapshots are protected against silent replacement.

The coverage/population and interval-compatibility diagnostics still
describe recorded REACHES temperature events; inferred-zero site-years
should not be misinterpreted as directly observed temperature events.
Figure 2 and Supplementary Figure S2 distinguish inferred zero
site-years with outlined markers. The Original descriptive count
plots use recorded documentary entries, while the inferred-zero
version uses aggregated site-years; use the explicit site-year audit
tables for a denominator-matched comparison of reporting density.
The supplementary spatial comparison aligns geographic locations
even when Figure 6 cluster memberships change between versions.

The legacy `run_all.R --only=<step>` entry point still supports
individual workflow steps, which write to root-level `Output/`.
For both complete versions and their separate figure outputs, use
`run_all.R --all`.

## Master Workflow

The master workflow is defined in:

```text
run_all.R
```

It executes each script in a separate R session, stops when a step fails, verifies key outputs, and records the execution status and runtime of each step.

A full two-version manuscript reproduction is run with:

```bash
REPRODUCTION_MODE=full \
Rscript --vanilla run_all.R --all
```

A faster smoke test is available with:

```bash
REPRODUCTION_MODE=smoke \
Rscript --vanilla run_all.R --all
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
Output/Reconstructions/<version>/Output/Logs/run_manifest.csv
```

Individual script logs are written under:

```text
Output/Reconstructions/<version>/Output/Logs/
```

During the one-time shared preparation, `Prior.R` prints its CV
combination progress **live in the terminal** and simultaneously saves it
to `Output/Reconstructions/shared/Output/Logs/*_prior.log`. The other
steps continue writing their detailed output to separate log files.
If a verified shared prior is reused, the Prior step is skipped and no
new CV progress is printed.

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
- converts the prepared temperatures from degrees Celsius to Kelvin to preserve the existing downstream data interface;
- calculates the LME ensemble mean;
- interpolates the LME fields to the nominal GHCN station coordinates for Beijing, Shanghai, and Hong Kong; and
- saves the annual and city-specific LME products.

Analysis scripts convert the stored Kelvin temperatures back to degrees Celsius where required.

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

This script applies WorldClim-based climatological bias corrections to historical CESM-LME temperature fields.

The original WorldClim and modern LME data are processed separately using:

```text
Data/Get_data/GetWorldClim.R
```

This optional upstream preparation generates:

```text
Data/WorldClim/worldclim_offsets_all_269.csv
```

The prepared input contains additive corrections for 266 native LME grid locations and three validation stations.

For each location $s$, the correction is defined by:

$$
d(s) = \overline{W}(s) - \overline{L}_{1970:2000}(s),
$$

where $\overline{W}(s)$ is the WorldClim climatological mean and $\overline{L}_{1970:2000}(s)$ is the corresponding mean of the 13 LME members over 1970–2000.

#### Historical LME correction

The historical LME annual temperatures are first prepared by:

```text
Code/DataPreparation/prepare_lme_annual.R
```

The three-city historical LME series have already undergone thin-plate spline interpolation.

The WorldClim correction is applied afterward as a constant location-specific additive shift:

$$ T_{\mathrm{corrected}}(s,t,j) = T_{\mathrm{original}}(s,t,j) + d(s), $$

where $t$ denotes calendar year and $j$ denotes ensemble member.

The same correction is used for every historical year and ensemble member.

The original uncorrected LME products are preserved.

#### Execution

Run from the repository root:

```bash
Rscript --vanilla Code/DataPreparation/prepare_worldclim_bias_correction.R
```

The script reads the prepared WorldClim offsets and historical LME RDS files. It does not read the original NetCDF or GeoTIFF archives.

#### Outputs

Corrected products are saved under:

```text
Output/Intermediate/LME/WorldClim_Corrected/
```

Important outputs include:

```text
lme_city3_annual_1368_1911.csv
lme_city3_annual_1368_1911.rds
lme_land_grid_corrected_annual_1368_1911.rds
worldclim_grid_correction_status.csv
worldclim_correction_diagnostics.csv
```

The corrected city CSV preserves the `temperature_kelvin` column expected by the downstream analysis scripts.

The corrected city data are used in:

- `Figure3.R`;
- `Figure7-8.R`;
- `Prior.R`; and
- `Figure9d.R`.


#### Spline Consistency and Reconstruction Diagnostics

The historical and modern LME interpolation procedures use identical thin-plate regression-spline settings:

- Spline basis: `bs = "tp"`;
- Basis dimension: `k = 266`;
- Fixed smoothing parameter: `sp = 1e-8`;
- Coordinates: untransformed longitude and latitude in decimal degrees;
- Fitting weights: equal weights for all 266 native LME grid locations.

The spline calculations are implemented directly in the two preparation scripts without requiring an additional shared R file.

All 266 native grid locations are retained for spline fitting. Native-grid reconstruction errors are evaluated only at locations with valid WorldClim land references.

The historical diagnostics include:

```text
Output/Intermediate/LME/WorldClim_Corrected/
    lme_spline_native_grid_errors_historical.csv
    lme_spline_specification_historical.csv
    lme_spline_audit_summary.csv
```

The script also checks whether the common spline settings reproduce the previously generated historical city series and verifies consistency with the modern LME interpolation specification.

These diagnostics do not change the WorldClim bias-correction formula.


#### Spatial correction and missing WorldClim values

WorldClim provides direct reference values only for locations covered by its land-temperature dataset.

The prepared correction table retains all 266 LME grid locations and three validation stations, including locations without valid WorldClim values.

Corrections are applied only where a direct WorldClim value is available. Missing offshore corrections are not automatically replaced by inverse-distance weighting (IDW).

The corrected land-grid RDS is a separate subset product and cannot directly replace the complete 266-location historical LME archive used by `Assimilation_grid.R`.

The treatment of offshore locations in spatial data assimilation remains separate from the three-city validation workflow.

#### Upstream reproducibility

The optional upstream procedure, original data sources, download instructions, and required external files are documented in:

```text
Data/Get_data/GetWorldClim.R
Data/Get_data/README.md
```

The default reproducibility workflow does not require the original large NetCDF or WorldClim GeoTIFF files.



### Supporting documentary-data diagnostics (optional)

The exploratory precipitation-event analysis and the independent
inferred-zero provenance check are kept in
[SupportingAnalyses/README.md](SupportingAnalyses/README.md).
They do not enter the main manuscript figure-generation pipeline.

### Preparing the documentary temperature input

`Code/DataPreparation/prepare_precipitation_augmented_reaches.R` constructs an optional site-year temperature-index dataset from the annual temperature and precipitation worksheets. It preserves the recorded temperature category wherever one exists and assigns a working index of zero only to site-years with an annual precipitation record and no recorded temperature category. This is a documentary-availability assumption, not a direct temperature observation.

The procedure keeps source flags, original text/event identifiers, and precipitation event codes. Annual precipitation records are used as evidence of documentary availability regardless of their rain/snow/hail or dry-weather classification. The original spreadsheets remain unchanged.

Default inputs:

```text
Data/temperature index value.v1.xlsx
Data/precipitation index value.xlsx
```

Output:

```text
Output/Intermediate/REACHES/temperature_site_year_augmented.csv
```

Additional audit and yearly/block summaries are written to the same directory. The source workbooks use longitude and latitude rounded to eight decimal places for site-year matching to avoid spurious floating-point duplicates; observed temperature categories remain unchanged. Recreate the augmented CSV after changing the preparation script.

For the **canonical manuscript reproduction**, run the versioned controller
from the repository root:

```bash
REPRODUCTION_MODE=full Rscript --vanilla run_all.R --all
```

This automatically prepares the precipitation-augmented site-year input,
uses it for the `inferred-zero` primary reconstruction, and uses the
recorded-temperature input for the `original` sensitivity reconstruction.
The resulting figures, tables, and intermediate products are isolated under
`Output/Reconstructions/inferred-zero/` and
`Output/Reconstructions/original/`, respectively; they do **not** overwrite
each other.

For **manual, legacy single-version steps only**, you may instead prepare
the augmented input and explicitly select it for calibration:

```bash
Rscript --vanilla Code/DataPreparation/prepare_precipitation_augmented_reaches.R
REACHES_TEMPERATURE_INPUT=augmented Rscript --vanilla Code/prepare_calibration.R
REACHES_TEMPERATURE_INPUT=augmented Rscript --vanilla Code/Get_tempe_all_data.R
```

These standalone commands write to root-level `Output/` paths, so do not
use them to compare the two final versions. The precipitation workbook path
can be overridden with `REACHES_PRECIPITATION_FILE`. The versioned controller
selects `REACHES_TEMPERATURE_INPUT` for each run automatically.

### Additional measurement variance

`Code/cv_observation_variance.R` is a **required upstream estimation
step** for the default full two-version reconstruction. It estimates
the nonnegative additional measurement variance `u` using five temporal
holdout blocks of the **original recorded-temperature input only**.
The root `run_all.R` automatically runs this estimator after preparing
the shared LME/prior, before launching either documentary-input
reconstruction. It reuses the saved estimate only when input-data,
CV-script and shared-prior fingerprints and CV settings are compatible.

The estimate and its diagnostic score curves are written under:

```text
Output/Reconstructions/original/Output/Intermediate/ObservationVariance/original/
  selected_u.csv
  cv_score_curve.csv
  cv_score_by_fold.csv
  cv_score_by_city.csv
  cv_score_curve.png
```

Both manuscript versions apply the same Original-input estimate to isolate
the effect of inferred-zero documentary observations. The measurement
variance loader and mathematical helpers remain in `Code/Figure9d.R`,
which is used by spatial assimilation. Additional numerical testing,
independent CV validation and precipitation-input diagnostics are **optional**
and documented in [SupportingAnalyses/README.md](SupportingAnalyses/README.md).

For a deliberate fixed-`u` sensitivity run, `run_all.R --all --u=VALUE`
bypasses CV; the supplied value must be independently justified. Smoke
mode defaults to `u=0` and is not suitable for manuscript figures.

### Comparing reconstructions with and without inferred zeros

The root [run_all.R](../run_all.R) now controls both versions:
`inferred-zero` is the primary reconstruction and `original` is
the temperature-record-only sensitivity analysis. It generates and
preserves separate figures, supplementary panels, tables and
intermediate results, while LME/WorldClim preparation and Prior.R
run only once.

```bash
REPRODUCTION_MODE=full Rscript --vanilla run_all.R --all
```

The controller automatically calls the four necessary
`Code/Analysis/` scripts and compares the two versions.
For optional city-only quick checks, see
[SupportingAnalyses/README.md](SupportingAnalyses/README.md).

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

The spatial assimilation uses the fixed penalty values documented in the script. By default, `location_scope = "figure6"` processes the locations selected for Figure 6; setting `location_scope = "all_lme"` processes all native LME grid locations.

#### Measurement-Equation Projection

For each location, quantile mapping converts the kriged REACHES index to the LME temperature scale. A joint Monte Carlo calculation samples from the empirical distribution of annual kriging predictions and their associated prediction uncertainties. For each informative year, the corresponding kriging uncertainty determines a simulated pair of latent and kriged indices. Applying the fitted quantile map to both samples yields the intercept, slope, and projection-residual variance in the working measurement equation.

The measurement Monte Carlo sample size defaults to 100,000 draws with a fixed random seed. The location-specific prior fitting uses a maximum of 150 iterations. The diagnostics file records the Monte Carlo configuration, measurement coefficients, prior convergence information, posterior displacement relative to LME, and numerical checks. Checkpoints are accepted only when their stored settings match the current configuration.

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

#### Measurement-Equation Projection

The measurement-equation coefficients are estimated through a joint Monte Carlo projection using the empirical distribution of the kriged REACHES indices and their prediction standard deviations. Each draw samples a latent index from the empirical reference mixture, then generates a corresponding kriged index using the year-specific information fraction and kriging uncertainty. The fitted quantile map transforms both indices to the temperature scale.

For each informative year, a linear projection of the simulated calibrated proxy on the simulated latent temperature provides the intercept, slope, and residual variance used in the Kalman measurement equation. The default sample size is 100,000 draws with a fixed random seed. The assimilation metadata also records the reference mean, reference sample size, Monte Carlo variance, posterior-minus-LME displacement, and the frequency with which the posterior falls below both contemporaneous inputs.

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

The public entry point is the repository-root `run_all.R --all`.
The controller calls the internal `Code/Workflow/run_single_version.R`
worker to execute each version in its own workspace; it does **not**
estimate the independent LME or penalized prior twice.

```text
run_all.R --all
    |
    +-- Shared preparation (once)
    |     prepare_lme_annual.R
    |       -> prepare_worldclim_bias_correction.R
    |       -> Prior.R
    |
    +-- Code/cv_observation_variance.R
    |     Five-block CV on Original temperature records;
    |     obtain one common observation variance u
    |
    +-- Inferred-zero reconstruction
    |     prepare_precipitation_augmented_reaches.R
    |       -> prepare_calibration.R
    |       -> Get_tempe_all_data.R
    |       -> manuscript figures, assimilation, tables and validation
    |
    +-- Original reconstruction
    |     prepare_calibration.R
    |       -> Get_tempe_all_data.R
    |       -> manuscript figures, assimilation, tables and validation
    |
    +-- Snapshot and compare both completed reconstructions
```

Within each version, `Code/Workflow/run_single_version.R` manages
the detailed step ordering, logs, expected-output checks and
`run_manifest.csv`. The shared prior and LME intermediates are
reused by both versions. The shared Original-input CV estimate is
also applied to both versions.

The supplementary spatial analysis depends on the main spatial
reconstruction:

```text
Figure6.R
    -> Assimilation_grid.R
    -> Supplementary/FigureS1.R
```

The WorldClim-based correction uses the prepared offsets stored in
`Data/WorldClim/worldclim_offsets_all_269.csv`. The original
CESM-LME NetCDF and WorldClim GeoTIFF archives are not required
for this workflow. Their optional upstream preparation is documented
in `Data/Get_data/GetWorldClim.R`.

## Generated Files

The current two-version workflow writes generated files under:

```text
Output/Reconstructions/shared/Output/
Output/Reconstructions/inferred-zero/Output/
Output/Reconstructions/original/Output/
Output/Reconstructions/Comparison/
```

The `Output/` directory is excluded from version control and is
created automatically when the workflow is run. Documentation for
individual scripts uses workspace-relative `Output/...` paths;
legacy individual-stage execution may instead write to the
repository-root `Output/` directories.

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
