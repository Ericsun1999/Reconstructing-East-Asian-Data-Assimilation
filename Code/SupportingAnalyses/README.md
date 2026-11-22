# Supporting analyses (optional)

This directory preserves **method checks, sensitivity diagnostics, and
exploratory analyses** developed alongside the manuscript. They are
scientifically useful and remain documented, but **none of these six
scripts is required to regenerate Figures 2–10, Supplementary Figures
S1–S5, or the principal manuscript tables**.

The complete publication workflow is driven by the repository-root
`run_all.R`. Scripts here are *not* called automatically by that
workflow. Do not mistake their output for the primary manuscript
reconstruction.

## Validation and method diagnostics

| File | What it checks | When to use |
|---|---|---|
| `cv_observation_variance.R` | Five-block temporal CV for the additional Celsius-scale measurement variance, plus fold-specific kriging and predictive-log-score curves | To re-estimate or investigate the supplementary variance choice |
| `validate_observation_variance_cv.R` | Independently recomputes CV score components, checks fold partitioning and source MD5, and diagnoses low-information city-years | After a complete CV fit |
| `test_observation_variance.R` | Fast numerical checks of annual state variance, additive observation variance and both measurement-projection implementations | During methodological development |
| `diagnose_precipitation_augmentation.R` | Checks inferred-zero site-year provenance, duplicate coordinates, unchanged recorded temperature indices and before/after coverage | To audit the documentary-availability assumption |

These checks are not a replacement for the manuscript's primary
reconstruction, and the complete manuscript workflow does not need to
refit the additional variance at every run.

The standalone scripts, when run from the **repository root**, are:

```bash
Rscript --vanilla Code/SupportingAnalyses/test_observation_variance.R
Rscript --vanilla Code/SupportingAnalyses/cv_observation_variance.R
Rscript --vanilla Code/SupportingAnalyses/validate_observation_variance_cv.R
Rscript --vanilla Code/SupportingAnalyses/diagnose_precipitation_augmentation.R
```

The CV scripts need independently generated LME/Prior inputs; the
augmentation diagnosis needs the generated site-year input.
Re-estimating a new variance with another measurement equation must
not be silently mixed with the original-model variance used in the
main two-version comparison. The historical original-model estimate
can be applied with a **common fixed value** to both documentary
versions to isolate the data-input sensitivity.

The mathematical helper functions for loading the additional variance,
computing the unconditional annual prior variance and evaluating
Gaussian predictive log score are defined in `Code/Figure9d.R`.
These diagnostic scripts extract the named definitions without
executing the full Figure 9(d) workflow.

## Exploratory and convenience studies

| File | Purpose |
|---|---|
| `kriging_city_only.R` | Quickly estimates only the Beijing, Shanghai and Hong Kong kriging series for a preliminary check; does not replace full spatial kriging |
| `precip_temperature_joint_distribution.R` | Studies observed temperature-index frequencies alongside preliminary precipitation-event categories; does not impute temperature observations |

Run separately if needed:

```bash
Rscript --vanilla Code/SupportingAnalyses/kriging_city_only.R
Rscript --vanilla Code/SupportingAnalyses/precip_temperature_joint_distribution.R
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
- `Code/Figure9d.R` defines the measurement-variance loader.
- `run_all.R` at the repository root is the public, complete
  dual-version entry point; the internal single-version executor is
  `Code/Workflow/run_single_version.R`.

Inferred-zero site-years are **assumed normal-temperature indices** at
documented site-years, not independently observed temperature records.
Their uncertainty is assessed through the Original-versus-inferred-zero
sensitivity comparison.