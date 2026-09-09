############################################################
# 10_OLIGO_LANDMARK_deployment_reference.R
# Cleaned public version of shrinkage / baseline re-estimation / prediction.
############################################################

SCRIPT_DIR <- dirname(normalizePath(sys.frame(1)$ofile, winslash = "/", mustWork = FALSE))
source(file.path(SCRIPT_DIR, "00_OLIGO_LANDMARK_config_helpers.R"))

stage3b_dir <- file.path(OUTPUT_ROOT, "07B_stage3B_expanded_biopsy")
stage4_dir <- file.path(OUTPUT_ROOT, "08_stage4_final_bootstrap")
out_dir <- file.path(OUTPUT_ROOT, "09_stage5_deployment_calibration")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

imp <- readRDS(file.path(stage3b_dir, "06_mice_expanded_biopsy_cohort.rds"))
beta_tbl <- readr::read_csv(file.path(stage4_dir, "02_pooled_final_coefficients_beta.csv"), show_col_types = FALSE)
shrink_tbl <- readr::read_csv(file.path(stage4_dir, "14_provisional_shrinkage_factor.csv"), show_col_types = FALSE)

shrinkage <- as.numeric(shrink_tbl$provisional_uniform_shrinkage_factor[[1]])
beta_raw <- setNames(beta_tbl$estimate, beta_tbl$term)
beta_deploy <- beta_raw * shrinkage

expected_terms <- c(
  "age", "molecular_gradeGrade 3", "tumor_size_mm", "cdcc1", "cdcc2+",
  "procedure_groupSubtotal/partial resection", "procedure_groupGross-total resection", "sexMale"
)
stop_if(!setequal(names(beta_deploy), expected_terms), "Coefficient terms do not match frozen model.")
beta_deploy <- beta_deploy[expected_terms]

completed_sets <- lapply(seq_len(imp$m), function(i) {
  d <- mice::complete(imp, i) %>% droplevels()
  d$procedure_group <- factor(as.character(d$procedure_group),
                              levels = c("Biopsy/local excision", "Subtotal/partial resection", "Gross-total resection"))
  d
})

make_X <- function(d) {
  mm <- model.matrix(~ age + molecular_grade + tumor_size_mm + cdcc + procedure_group + sex, data = d)
  mm <- mm[, colnames(mm) != "(Intercept)", drop = FALSE]
  mm[, expected_terms, drop = FALSE]
}

reference <- data.frame(
  age = 45,
  molecular_grade = factor("Grade 2", levels = levels(completed_sets[[1]]$molecular_grade)),
  tumor_size_mm = 50,
  cdcc = factor("0", levels = levels(completed_sets[[1]]$cdcc)),
  procedure_group = factor("Biopsy/local excision", levels = levels(completed_sets[[1]]$procedure_group)),
  sex = factor("Female", levels = levels(completed_sets[[1]]$sex))
)
Xref <- as.numeric(make_X(reference)[1, ])
names(Xref) <- expected_terms

baseline_by_imp <- lapply(seq_along(completed_sets), function(i) {
  d <- completed_sets[[i]]
  X <- make_X(d)
  lp_centered <- as.numeric((X - matrix(Xref, nrow(X), length(Xref), byrow = TRUE)) %*% beta_deploy)
  fit_offset <- survival::coxph(survival::Surv(landmark_time_months, event) ~ offset(lp_centered),
                                data = d, ties = "efron")
  bh <- survival::basehaz(fit_offset, centered = FALSE)
  data.frame(imputation = i, time_months_after_landmark = bh$time,
             cumulative_baseline_hazard_reference = bh$hazard)
})

grid <- sort(unique(unlist(lapply(baseline_by_imp, `[[`, "time_months_after_landmark"))))
Hmat <- sapply(baseline_by_imp, function(b) {
  approx(c(0, b$time_months_after_landmark),
         c(0, b$cumulative_baseline_hazard_reference),
         xout = grid, method = "constant", f = 0, rule = 2)$y
})
H0 <- rowMeans(Hmat)
baseline_curve <- tibble(
  time_months_after_landmark = grid,
  cumulative_baseline_hazard_reference = H0,
  baseline_survival_reference = exp(-H0)
)
readr::write_csv(baseline_curve, file.path(out_dir, "02_deployment_baseline_survival_full_curve.csv"))

horizons <- c(12, 24, 36)
Hh <- approx(c(0, grid), c(0, H0), xout = horizons, method = "constant", f = 0, rule = 2)$y
baseline_horizons <- tibble(
  horizon_months = horizons,
  cumulative_baseline_hazard_reference = Hh,
  baseline_survival_reference = exp(-Hh)
)
readr::write_csv(baseline_horizons, file.path(out_dir, "03_deployment_baseline_survival_horizons.csv"))

deployment <- tibble(term = expected_terms, development_beta = beta_raw[expected_terms],
                     shrinkage_factor = shrinkage, deployment_beta = beta_deploy)
readr::write_csv(deployment, file.path(out_dir, "04_FINAL_deployment_coefficients.csv"))

predict_oligo_survival <- function(age, grade, tumor_size_mm, cdcc, procedure, sex, horizon_months) {
  lp <- beta_deploy["age"] * (age - 45) +
    beta_deploy["molecular_gradeGrade 3"] * as.numeric(grade == "Grade 3") +
    beta_deploy["tumor_size_mm"] * (tumor_size_mm - 50) +
    beta_deploy["cdcc1"] * as.numeric(cdcc == "1") +
    beta_deploy["cdcc2+"] * as.numeric(cdcc == "2+") +
    beta_deploy["procedure_groupSubtotal/partial resection"] * as.numeric(procedure == "Subtotal/partial resection") +
    beta_deploy["procedure_groupGross-total resection"] * as.numeric(procedure == "Gross-total resection") +
    beta_deploy["sexMale"] * as.numeric(sex == "Male")
  H <- approx(c(0, grid), c(0, H0), xout = horizon_months, method = "constant", f = 0, rule = 2)$y
  exp(-H) ^ exp(lp)
}

stop_if(abs(predict_oligo_survival(45, "Grade 2", 50, "0", "Biopsy/local excision", "Female", 12) - 0.9811932058) > 1e-5,
        "Reference 12-month prediction does not reproduce the frozen deployment model.")

readr::write_csv(
  tibble(status = "COMPLETE", shrinkage_factor = shrinkage,
         reference_survival_12m = predict_oligo_survival(45, "Grade 2", 50, "0", "Biopsy/local excision", "Female", 12),
         reference_survival_24m = predict_oligo_survival(45, "Grade 2", 50, "0", "Biopsy/local excision", "Female", 24),
         reference_survival_36m = predict_oligo_survival(45, "Grade 2", 50, "0", "Biopsy/local excision", "Female", 36)),
  file.path(out_dir, "99_STAGE5_COMPLETE.csv")
)
