# Oligodendroglioma 180-day landmark prediction model: analysis code

This directory contains cleaned, path-portable versions of the final R analysis code supporting the manuscript and deployed calculator.

## Scope

The public code covers the final analysis from the licensed molecular cohort through:

1. treatment-timing / immortal-time sensitivity analysis;
2. predictor and functional-form specification checks;
3. expanded biopsy coding and final 180-day landmark cohort construction;
4. multiple imputation of tumor size;
5. final Cox model fitting and proportional-hazards diagnostics;
6. 20 x 1,000 bootstrap internal validation;
7. uniform shrinkage, baseline-hazard re-estimation, calibration, and deployment.

The frozen model is:

`Surv(post-landmark overall survival, death) ~ age + molecular grade + tumor size + Charlson-Deyo score + index tissue procedure + sex`

with prediction horizons of 12, 24, and 36 months after a 180-day landmark following the index tissue procedure.

## Data availability

The National Cancer Database Participant User File is governed by its data-use agreement and cannot be redistributed here. Investigators with authorized NCDB access must supply their own local 2024 Brain PUF and the molecular cohort described in the manuscript.

Before running the scripts, set:

```r
Sys.setenv(
  OLIGO_PROJECT_DIR = "/path/to/oligodendroglioma_project",
  NCDB_BRAIN_PUF_CSV = "/path/to/NCDBPUF_Brain.0.2024.0.csv"
)
```

The expected molecular-cohort file is:

`<OLIGO_PROJECT_DIR>/data/processed/01_oligo_molecular_all_eligible_years.rds`

## Files

- `00_OLIGO_LANDMARK_config_helpers.R` - configuration, registry harmonization, and performance helpers.
- `04_OLIGO_TIMEVARYING_treatment_sensitivity.R` - naive versus time-varying treatment sensitivity.
- `07_OLIGO_LANDMARK_stage3_final_specification.R` - predictor and functional-form checks.
- `09_OLIGO_LANDMARK_final_cohort_imputation.R` - final biopsy coding, landmark cohort, and multiple imputation.
- `08_OLIGO_LANDMARK_final_model_validation.R` - final model fitting, PH diagnostics, and 20,000-fit bootstrap internal validation.
- `10_OLIGO_LANDMARK_deployment_reference.R` - uniform shrinkage, baseline-hazard re-estimation, and reference prediction implementation.

The exact numerical deployment implementation used by the web calculator is also embedded in `shiny/app.R`.

## Reproducibility note

Because the licensed NCDB PUF itself cannot be redistributed, these scripts are not a one-click public reproduction from raw data. They are provided to make the cohort logic, model specification, validation, and deployment procedures transparent and auditable.

This is a prognostic model. Procedure and treatment coefficients in sensitivity analyses must not be interpreted as causal treatment effects. Independent external validation is required before routine clinical use.
