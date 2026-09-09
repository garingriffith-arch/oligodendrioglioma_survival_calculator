############################################################
# 07_OLIGO_LANDMARK_stage3_final_specification.R
#
# STAGE 3: FINAL MODEL-SPECIFICATION AUDIT
############################################################

SCRIPT_DIR <- dirname(normalizePath(sys.frame(1)$ofile, winslash = "/", mustWork = FALSE))
source(file.path(SCRIPT_DIR, "00_OLIGO_LANDMARK_config_helpers.R"))

for (pkg in c("mice", "survival", "broom", "splines")) {
  if (!requireNamespace(pkg, quietly = TRUE)) stop("Package '", pkg, "' is required.", call. = FALSE)
}

L <- 180L
in_path <- file.path(PROCESSED_DIR, sprintf("01_landmark_%03dd_all.rds", L))
out_dir <- file.path(OUTPUT_ROOT, "07_stage3_final_specification")
if (dir.exists(out_dir)) unlink(out_dir, recursive = TRUE, force = TRUE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
stopifnot(file.exists(in_path))
dat_all <- readRDS(in_path)

dat_all <- dat_all %>%
  mutate(
    procedure_group = factor(as.character(procedure_group), levels = c("Subtotal/partial resection", "Gross-total resection", "Biopsy/local excision")),
    sex = droplevels(factor(sex)),
    molecular_grade = droplevels(factor(molecular_grade)),
    cdcc = droplevels(factor(cdcc)),
    primary_site_group = droplevels(factor(primary_site_group))
  )

required_core <- c("landmark_time_months", "event", "age", "sex", "molecular_grade", "cdcc", "primary_site_group", "procedure_group")
missing_required <- setdiff(c(required_core, "tumor_size_mm"), names(dat_all))
if (length(missing_required) > 0L) stop("Required Stage 3 variables missing: ", paste(missing_required, collapse = ", "), call. = FALSE)

core_all <- dat_all %>% filter(if_all(all_of(required_core), ~ !is.na(.x))) %>% droplevels()

flow <- tibble(
  step = c("180-day landmark eligible with required non-size core predictors", "Observed tumor size (complete-case reference only)", "Tumor size missing and retained for MI"),
  n = c(nrow(core_all), sum(!is.na(core_all$tumor_size_mm)), sum(is.na(core_all$tumor_size_mm))),
  deaths = c(sum(core_all$event == 1L), sum(core_all$event == 1L & !is.na(core_all$tumor_size_mm)), sum(core_all$event == 1L & is.na(core_all$tumor_size_mm)))
)
readr::write_csv(flow, file.path(out_dir, "00_stage3_analysis_flow.csv"))

parameter_count <- tibble(
  model = c("M1 lean", "M2 + sex", "M3 + primary site"),
  nominal_parameters = c(7L, 8L, 11L),
  deaths = sum(core_all$event == 1L)
) %>% mutate(events_per_parameter = deaths / nominal_parameters)
readr::write_csv(parameter_count, file.path(out_dir, "00A_events_per_parameter.csv"))

add_nelson_aalen <- function(d) {
  fit0 <- survival::coxph(survival::Surv(landmark_time_months, event) ~ 1, data = d, ties = "efron")
  bh <- survival::basehaz(fit0, centered = FALSE)
  idx <- findInterval(d$landmark_time_months, bh$time)
  H <- numeric(nrow(d)); H[idx > 0] <- bh$hazard[idx[idx > 0]]
  d$nelson_aalen <- H
  d
}

core_all <- add_nelson_aalen(core_all)
if (!"dx_year_final" %in% names(core_all)) core_all$dx_year_final <- NA_real_
mi_vars <- c("landmark_time_months", "event", "age", "sex", "molecular_grade", "tumor_size_mm", "cdcc", "primary_site_group", "procedure_group", "dx_year_final", "nelson_aalen")
mi_dat <- core_all %>% select(all_of(mi_vars))

ini <- mice::mice(mi_dat, maxit = 0, printFlag = FALSE)
meth <- ini$method; pred <- ini$predictorMatrix; meth[] <- ""; meth["tumor_size_mm"] <- "pmm"; pred[,] <- 0
size_predictors <- c("landmark_time_months", "event", "age", "sex", "molecular_grade", "cdcc", "primary_site_group", "procedure_group", "dx_year_final", "nelson_aalen")
is_usable_predictor <- function(x) { xx <- x[!is.na(x)]; length(xx) > 1L && length(unique(xx)) > 1L }
present_preds <- size_predictors[size_predictors %in% names(mi_dat)]
usable_size_predictors <- present_preds[vapply(mi_dat[present_preds], is_usable_predictor, logical(1))]
pred["tumor_size_mm", usable_size_predictors] <- 1
pred["tumor_size_mm", "tumor_size_mm"] <- 0

readr::write_csv(tibble(predictor = size_predictors, present = size_predictors %in% names(mi_dat), used_for_tumor_size_imputation = size_predictors %in% usable_size_predictors), file.path(out_dir, "01_MI_predictor_audit.csv"))

set.seed(MI_SEED + 100L)
imp <- mice::mice(mi_dat, m = MI_M, maxit = MI_MAXIT, method = meth, predictorMatrix = pred, seed = MI_SEED + 100L, printFlag = TRUE)
saveRDS(imp, file.path(out_dir, "02_mice_full_core_cohort.rds"))
if (!is.null(imp$loggedEvents) && nrow(imp$loggedEvents) > 0L) readr::write_csv(as_tibble(imp$loggedEvents), file.path(out_dir, "02A_mice_logged_events.csv")) else readr::write_csv(tibble(message = "No logged MICE events"), file.path(out_dir, "02A_mice_logged_events.csv"))

completed_sets <- lapply(seq_len(MI_M), function(i) {
  d <- mice::complete(imp, action = i); d <- droplevels(d)
  d$procedure_group <- factor(as.character(d$procedure_group), levels = c("Subtotal/partial resection", "Gross-total resection", "Biopsy/local excision"))
  d
})

qc <- bind_rows(lapply(seq_along(completed_sets), function(i) {
  d <- completed_sets[[i]]
  tibble(imputation = i, n = nrow(d), deaths = sum(d$event == 1L), tumor_size_missing_n = sum(is.na(d$tumor_size_mm)))
}))
readr::write_csv(qc, file.path(out_dir, "03_completed_dataset_QC.csv"))
if (any(qc$n != nrow(core_all)) || any(qc$tumor_size_missing_n != 0L)) stop("Completed-dataset QC failed.", call. = FALSE)

f_m1 <- survival::Surv(landmark_time_months, event) ~ age + molecular_grade + tumor_size_mm + cdcc + procedure_group
f_m2 <- update(f_m1, . ~ . + sex)
f_m3 <- update(f_m2, . ~ . + primary_site_group)
fit_one <- function(d, formula) survival::coxph(formula, data = d, ties = "efron", x = TRUE, y = TRUE, model = TRUE)

fit_lists <- list(
  M1_lean = lapply(completed_sets, fit_one, formula = f_m1),
  M2_plus_sex = lapply(completed_sets, fit_one, formula = f_m2),
  M3_plus_site = lapply(completed_sets, fit_one, formula = f_m3)
)
miras <- lapply(fit_lists, mice::as.mira)
saveRDS(miras, file.path(out_dir, "04_candidate_model_fits_mira.rds"))

pool_hr <- function(mira_fit, model_name) {
  summary(mice::pool(mira_fit), conf.int = TRUE, exponentiate = TRUE) %>% as_tibble() %>% mutate(model = model_name, p_formatted = p_format(p.value)) %>% relocate(model)
}
pooled <- bind_rows(pool_hr(miras$M1_lean, "M1 lean"), pool_hr(miras$M2_plus_sex, "M2 + sex"), pool_hr(miras$M3_plus_site, "M3 + primary site"))
readr::write_csv(pooled, file.path(out_dir, "05_pooled_candidate_coefficients.csv"))

capture.output(mice::D1(miras$M2_plus_sex, miras$M1_lean), file = file.path(out_dir, "06_D1_increment_sex.txt"))
capture.output(mice::D1(miras$M3_plus_site, miras$M2_plus_sex), file = file.path(out_dir, "07_D1_increment_primary_site.txt"))

c_rows <- bind_rows(lapply(seq_len(MI_M), function(i) {
  d <- completed_sets[[i]]
  bind_rows(lapply(names(fit_lists), function(nm) {
    lp <- stats::predict(fit_lists[[nm]][[i]], newdata = d, type = "lp")
    tibble(imputation = i, model = nm, c_index = harrell_c(d, lp))
  }))
}))
readr::write_csv(c_rows, file.path(out_dir, "08_apparent_C_by_imputation.csv"))

c_wide <- c_rows %>% select(imputation, model, c_index) %>% tidyr::pivot_wider(names_from = model, values_from = c_index) %>% mutate(delta_C_sex = M2_plus_sex - M1_lean, delta_C_site = M3_plus_site - M2_plus_sex)
readr::write_csv(c_wide, file.path(out_dir, "08A_incremental_C_by_imputation.csv"))
readr::write_csv(tibble(comparison = c("M2 + sex vs M1 lean", "M3 + site vs M2 + sex"), mean_delta_C = c(mean(c_wide$delta_C_sex), mean(c_wide$delta_C_site)), sd_delta_C = c(sd(c_wide$delta_C_sex), sd(c_wide$delta_C_site)), min_delta_C = c(min(c_wide$delta_C_sex), min(c_wide$delta_C_site)), max_delta_C = c(max(c_wide$delta_C_sex), max(c_wide$delta_C_site))), file.path(out_dir, "09_incremental_predictor_summary.csv"))

ph_all <- bind_rows(lapply(names(fit_lists), function(nm) bind_rows(lapply(seq_len(MI_M), function(i) {
  z <- survival::cox.zph(fit_lists[[nm]][[i]], terms = TRUE, global = TRUE)
  as.data.frame(z$table) %>% tibble::rownames_to_column("term") %>% as_tibble() %>% transmute(model = nm, imputation = i, term, chisq = chisq, df = df, p = p)
}))))
readr::write_csv(ph_all, file.path(out_dir, "10_PH_tests_all_candidates.csv"))
readr::write_csv(ph_all %>% group_by(model, term) %>% summarise(imputations = n(), median_p = median(p, na.rm = TRUE), min_p = min(p, na.rm = TRUE), max_p = max(p, na.rm = TRUE), proportion_p_lt_0_05 = mean(p < 0.05, na.rm = TRUE), .groups = "drop"), file.path(out_dir, "11_PH_tests_summary.csv"))

f_age_spline <- survival::Surv(landmark_time_months, event) ~ splines::ns(age, df = 3) + molecular_grade + tumor_size_mm + cdcc + procedure_group + sex
f_size_spline <- survival::Surv(landmark_time_months, event) ~ age + molecular_grade + splines::ns(tumor_size_mm, df = 3) + cdcc + procedure_group + sex
f_both_spline <- survival::Surv(landmark_time_months, event) ~ splines::ns(age, df = 3) + molecular_grade + splines::ns(tumor_size_mm, df = 3) + cdcc + procedure_group + sex
spline_lists <- list(M2_linear = fit_lists$M2_plus_sex, age_spline = lapply(completed_sets, fit_one, formula = f_age_spline), size_spline = lapply(completed_sets, fit_one, formula = f_size_spline), age_and_size_spline = lapply(completed_sets, fit_one, formula = f_both_spline))

spline_perf <- bind_rows(lapply(names(spline_lists), function(nm) bind_rows(lapply(seq_len(MI_M), function(i) {
  fit <- spline_lists[[nm]][[i]]; d <- completed_sets[[i]]; lp <- stats::predict(fit, newdata = d, type = "lp")
  tibble(imputation = i, model = nm, c_index = harrell_c(d, lp), AIC = stats::AIC(fit))
}))))
readr::write_csv(spline_perf, file.path(out_dir, "12_continuous_form_performance_by_imputation.csv"))

spline_wide_c <- spline_perf %>% select(imputation, model, c_index) %>% tidyr::pivot_wider(names_from = model, values_from = c_index)
spline_wide_aic <- spline_perf %>% select(imputation, model, AIC) %>% tidyr::pivot_wider(names_from = model, values_from = AIC)
continuous_summary <- tibble(
  model = c("age_spline", "size_spline", "age_and_size_spline"),
  mean_delta_C_vs_linear = c(mean(spline_wide_c$age_spline - spline_wide_c$M2_linear), mean(spline_wide_c$size_spline - spline_wide_c$M2_linear), mean(spline_wide_c$age_and_size_spline - spline_wide_c$M2_linear)),
  mean_delta_AIC_vs_linear = c(mean(spline_wide_aic$age_spline - spline_wide_aic$M2_linear), mean(spline_wide_aic$size_spline - spline_wide_aic$M2_linear), mean(spline_wide_aic$age_and_size_spline - spline_wide_aic$M2_linear)),
  proportion_AIC_lower_than_linear = c(mean(spline_wide_aic$age_spline < spline_wide_aic$M2_linear), mean(spline_wide_aic$size_spline < spline_wide_aic$M2_linear), mean(spline_wide_aic$age_and_size_spline < spline_wide_aic$M2_linear))
)
readr::write_csv(continuous_summary, file.path(out_dir, "13_continuous_form_summary.csv"))

cc <- core_all %>% filter(!is.na(tumor_size_mm)) %>% droplevels()
cc_fit <- survival::coxph(f_m3, data = cc, ties = "efron", x = TRUE, y = TRUE, model = TRUE)
readr::write_csv(broom::tidy(cc_fit, exponentiate = TRUE, conf.int = TRUE) %>% mutate(p_formatted = p_format(p.value)), file.path(out_dir, "14_complete_case_full_model_reference.csv"))

readr::write_csv(tibble(completed = TRUE, timestamp = as.character(Sys.time()), landmark_days = L, m = MI_M, analysis_n = nrow(core_all), deaths = sum(core_all$event == 1L), tumor_size_missing_n = sum(is.na(core_all$tumor_size_mm))), file.path(out_dir, "99_STAGE3_COMPLETE.csv"))
