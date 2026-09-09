############################################################
# 08_OLIGO_LANDMARK_final_model_validation.R
# Cleaned public version of the final model fitting / validation logic.
############################################################

SCRIPT_DIR <- dirname(normalizePath(sys.frame(1)$ofile, winslash = "/", mustWork = FALSE))
source(file.path(SCRIPT_DIR, "00_OLIGO_LANDMARK_config_helpers.R"))

if (!requireNamespace("mice", quietly = TRUE)) stop("Package 'mice' is required.", call. = FALSE)

stage3b_dir <- file.path(OUTPUT_ROOT, "07B_stage3B_expanded_biopsy")
out_dir <- file.path(OUTPUT_ROOT, "08_stage4_final_bootstrap")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

imp <- readRDS(file.path(stage3b_dir, "06_mice_expanded_biopsy_cohort.rds"))
stop_if(!inherits(imp, "mids"), "Stage 3B object is not a mids object.")
M <- imp$m
completed_sets <- lapply(seq_len(M), function(i) {
  d <- mice::complete(imp, i) %>% droplevels()
  d$procedure_group <- factor(as.character(d$procedure_group),
                              levels = c("Biopsy/local excision", "Subtotal/partial resection", "Gross-total resection"))
  d
})

final_formula <- survival::Surv(landmark_time_months, event) ~
  age + molecular_grade + tumor_size_mm + cdcc + procedure_group + sex

fits <- lapply(completed_sets, function(d) survival::coxph(final_formula, data = d, ties = "efron",
                                                          x = TRUE, y = TRUE, model = TRUE))
mira <- mice::as.mira(fits)
saveRDS(mira, file.path(out_dir, "01_final_model_fits_mira.rds"))

pooled_beta <- summary(mice::pool(mira), conf.int = TRUE, exponentiate = FALSE) %>% as_tibble()
pooled_hr <- summary(mice::pool(mira), conf.int = TRUE, exponentiate = TRUE) %>% as_tibble()
readr::write_csv(pooled_beta, file.path(out_dir, "02_pooled_final_coefficients_beta.csv"))
readr::write_csv(pooled_hr, file.path(out_dir, "03_pooled_final_coefficients_HR.csv"))

readr::write_csv(
  tibble(analysis_n = nrow(completed_sets[[1]]), deaths = sum(completed_sets[[1]]$event == 1L),
         fitted_parameters = length(coef(fits[[1]])),
         events_per_parameter = sum(completed_sets[[1]]$event == 1L) / length(coef(fits[[1]]))),
  file.path(out_dir, "04_events_per_parameter.csv")
)

ph_all <- bind_rows(lapply(seq_len(M), function(i) {
  z <- survival::cox.zph(fits[[i]], terms = TRUE, global = TRUE)
  as.data.frame(z$table) %>% tibble::rownames_to_column("term") %>% as_tibble() %>%
    transmute(imputation = i, term, chisq = chisq, df = df, p = p)
}))
readr::write_csv(ph_all, file.path(out_dir, "05_PH_tests_all_imputations.csv"))
readr::write_csv(
  ph_all %>% group_by(term) %>% summarise(imputations = n(), median_p = median(p, na.rm = TRUE),
                                           min_p = min(p, na.rm = TRUE), max_p = max(p, na.rm = TRUE),
                                           n_p_lt_0_05 = sum(p < 0.05, na.rm = TRUE), .groups = "drop"),
  file.path(out_dir, "06_PH_tests_summary.csv")
)

horizons <- PREDICTION_HORIZONS_MONTHS

apparent_one <- function(d, fit) {
  lp <- as.numeric(predict(fit, newdata = d, type = "lp", reference = "zero"))
  out <- list(c = harrell_c(d, lp), slope = 1)
  for (h in horizons) {
    s <- pred_survival_at(fit, d, h)
    out[[paste0("auc_", h)]] <- ipcw_auc(d, -s, h)
    out[[paste0("brier_", h)]] <- ipcw_brier(d, s, h)
  }
  unlist(out)
}

bootstrap_one <- function(d, B, seed) {
  n <- nrow(d)
  app_fit <- survival::coxph(final_formula, data = d, ties = "efron", x = TRUE, y = TRUE, model = TRUE)
  app <- apparent_one(d, app_fit)

  set.seed(seed)
  boot <- lapply(seq_len(B), function(b) {
    idx <- sample.int(n, n, replace = TRUE)
    db <- d[idx, , drop = FALSE]
    fb <- try(survival::coxph(final_formula, data = db, ties = "efron", x = TRUE, y = TRUE, model = TRUE), silent = TRUE)
    if (inherits(fb, "try-error") || anyNA(coef(fb))) return(rep(NA_real_, 2 + 2 * length(horizons)))

    lp_train <- as.numeric(predict(fb, newdata = db, type = "lp", reference = "zero"))
    lp_test <- as.numeric(predict(fb, newdata = d, type = "lp", reference = "zero"))
    vals <- c(
      c_opt = harrell_c(db, lp_train) - harrell_c(d, lp_test),
      slope_opt = 1 - calibration_slope(d, lp_test)
    )
    for (h in horizons) {
      s_train <- pred_survival_at(fb, db, h)
      s_test <- pred_survival_at(fb, d, h)
      vals[paste0("auc_opt_", h)] <- ipcw_auc(db, -s_train, h) - ipcw_auc(d, -s_test, h)
      vals[paste0("brier_opt_", h)] <- ipcw_brier(db, s_train, h) - ipcw_brier(d, s_test, h)
    }
    vals
  })
  boot <- do.call(rbind, boot)
  list(apparent = app, optimism = colMeans(boot, na.rm = TRUE), boot = boot)
}

validated <- lapply(seq_len(M), function(i) bootstrap_one(completed_sets[[i]], BOOTSTRAP_B, BOOTSTRAP_SEED + i))
app_mat <- do.call(rbind, lapply(validated, `[[`, "apparent"))
opt_mat <- do.call(rbind, lapply(validated, `[[`, "optimism"))

corrected <- colMeans(app_mat, na.rm = TRUE)
corrected["c"] <- mean(app_mat[, "c"], na.rm = TRUE) - mean(opt_mat[, "c_opt"], na.rm = TRUE)
corrected["slope"] <- 1 - mean(opt_mat[, "slope_opt"], na.rm = TRUE)
for (h in horizons) {
  corrected[paste0("auc_", h)] <- mean(app_mat[, paste0("auc_", h)], na.rm = TRUE) -
    mean(opt_mat[, paste0("auc_opt_", h)], na.rm = TRUE)
  corrected[paste0("brier_", h)] <- mean(app_mat[, paste0("brier_", h)], na.rm = TRUE) -
    mean(opt_mat[, paste0("brier_opt_", h)], na.rm = TRUE)
}

readr::write_csv(
  tibble(metric = names(corrected), optimism_corrected_estimate = as.numeric(corrected)),
  file.path(out_dir, "11_FINAL_optimism_corrected_performance.csv")
)
readr::write_csv(
  tibble(provisional_uniform_shrinkage_factor = unname(corrected["slope"])),
  file.path(out_dir, "14_provisional_shrinkage_factor.csv")
)
readr::write_csv(
  tibble(status = "COMPLETE", imputations = M, bootstrap_per_imputation = BOOTSTRAP_B,
         total_requested_fits = M * BOOTSTRAP_B),
  file.path(out_dir, "99_STAGE4_COMPLETE.csv")
)
