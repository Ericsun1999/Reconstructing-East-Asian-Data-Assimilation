# Supporting analyses (optional)

This directory preserves **method checks, validation diagnostics, and
exploratory analyses** that are scientifically useful but are **not
required** to regenerate the main manuscript figures, supplementary
figures, or manuscript tables. The five scripts here are not executed
by the standard `run_all.R --all` workflow.

**Important:** The five-block estimator
[`Code/cv_observation_variance.R`](../cv_observation_variance.R)
is a required *upstream* step for the default full manuscript model.
It has been moved out of this optional directory and is executed
automatically by the versioned `run_all.R` after shared LME/prior
preparation. There is no need to run that CV estimator by hand.

## Independent validation and diagnostics

| File | What it checks | When to use |
|---|---|---|
| `validate_observation_variance_cv.R` | Independently recomputes CV score components, checks fold partitioning and source MD5, and diagnoses low-information city-years | After an automatic full CV fit |
| `test_observation_variance.R` | Fast numerical checks of annual state variance, additive observation variance and both measurement-projection implementations | During methodological development |
| `diagnose_precipitation_augmentation.R` | Checks inferred-zero site-year provenance, duplicate coordinates, unchanged recorded temperature indices and before/after coverage | To audit the documentary-availability assumption |

The mathematical helpers for loading the additional variance,
computing unconditional annual prior variance, and evaluating Gaussian
predictive log score remain in `Code/Figure9d.R`. Optional test
scripts extract these definitions without executing the entire figure
workflow.

The numerical test can run from the repository root:

```bash
Rscript --vanilla Code/SupportingAnalyses/test_observation_variance.R
```

After completing `run_all.R --all`, optional validation can be run in
the separate version workspaces:

```bash
# From the repository root:
cd Output/Reconstructions/original
Rscript --vanilla Code/SupportingAnalyses/validate_observation_variance_cv.R

# Independent diagnostic of the inferred-zero input:
cd ../inferred-zero
Rscript --vanilla Code/SupportingAnalyses/diagnose_precipitation_augmentation.R
```

The CV validator requires the full five-block fit generated in the
`original` workspace; the augmentation diagnosis requires the
prepared `inferred-zero` site-year data. These checks are optional.
Any alternative measurement-projection CV experiments should be
run separately and not substituted into the main reconstruction.

## Exploratory and convenience studies

| File | Purpose |
|---|---|
| `kriging_city_only.R` | Quickly estimates only the Beijing, Shanghai and Hong Kong kriging series for a preliminary check; does not replace full spatial kriging |
| `precip_temperature_joint_distribution.R` | Studies observed temperature-index frequencies alongside preliminary precipitation-event categories; does not impute temperature observations |

Run separately if needed:

```bash
# From the repository root, with standalone analysis inputs available:
Rscript --vanilla Code/SupportingAnalyses/precip_temperature_joint_distribution.R

# From a completed original or inferred-zero project workspace:
cd Output/Reconstructions/inferred-zero
Rscript --vanilla Code/SupportingAnalyses/kriging_city_only.R
```

The city-only kriging utility writes ordinary generated REACHES paths.
Do not use its partial grid products in the complete spatial
reproduction; that requires `Code/Get_tempe_all_data.R`.

## What belongs to the primary reproduction

- `Code/DataPreparation/prepare_precipitation_augmented_reaches.R`
  constructs the inferred-zero working input and retains its
  provenance flags.
- `Code/Analysis/coverage_population_GAM.R` and
  `Code/Analysis/interval_compatibility.R` generate manuscript
  table/validation products.
- `Code/Analysis/snapshot_augmentation_reconstruction.R` and
  `Code/Analysis/compare_augmentation_reconstructions.R` archive and
  compare the two completed versions.
- `Code/cv_observation_variance.R` estimates the required
  Original-input observation variance automatically during full
  manuscript reproduction.
- `Code/Figure9d.R` defines the measurement-variance loader.
- `run_all.R` at the repository root is the public, complete
  dual-version entry point; the internal single-version executor is
  `Code/Workflow/run_single_version.R`.

Inferred-zero site-years are **assumed normal-temperature indices** at
documented site-years, not independently observed temperature records.
Their uncertainty is assessed through the Original-versus-inferred-zero
sensitivity comparison.