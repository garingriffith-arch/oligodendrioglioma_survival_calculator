############################################################
# 00_OLIGO_LANDMARK_config_helpers.R
# Oligodendroglioma prognostic-model rebuild
# - postoperative landmark analysis
# - time-varying RT/chemotherapy sensitivity
# - bootstrap internal validation on the full cohort
############################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(survival)
  library(broom)
  library(data.table)
  library(parallel)
})

# ==========================================================
# USER CONFIGURATION
# ==========================================================
# Paths are intentionally supplied at run time because the NCDB PUF cannot be
# redistributed. Example before running:
#   Sys.setenv(
#     OLIGO_PROJECT_DIR = "/path/to/oligodendroglioma_project",
#     NCDB_BRAIN_PUF_CSV = "/path/to/NCDBPUF_Brain.0.2024.0.csv"
#   )
project_dir <- Sys.getenv("OLIGO_PROJECT_DIR", unset = "")
RAW_NCDB_CSV <- Sys.getenv("NCDB_BRAIN_PUF_CSV", unset = "")

if (!nzchar(project_dir)) {
  stop("Set the OLIGO_PROJECT_DIR environment variable before running the analysis.", call. = FALSE)
}
if (!nzchar(RAW_NCDB_CSV) || !file.exists(RAW_NCDB_CSV)) {
  stop(
    paste0(
      "Set NCDB_BRAIN_PUF_CSV to the local 2024 NCDB Brain PUF CSV before running.\n",
      "Current value: ", RAW_NCDB_CSV
    ),
    call. = FALSE
  )
}

# Contemporary molecular cohort created from the licensed NCDB Brain PUF,
# before the prior complete-case restriction.
MOLECULAR_RDS <- file.path(
  project_dir, "data", "processed", "01_oligo_molecular_all_eligible_years.rds"
)

PIPELINE_ROOT <- file.path(project_dir, "oligo_landmark_v2")
PROCESSED_DIR <- file.path(PIPELINE_ROOT, "data", "processed")
OUTPUT_ROOT   <- file.path(PIPELINE_ROOT, "outputs")
dir.create(PROCESSED_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(OUTPUT_ROOT, recursive = TRUE, showWarnings = FALSE)

# Primary landmark is provisional. The sweep script will tell us whether 90 vs 120 days
# captures treatment timing better. Do not change until reviewing 01_landmark_sweep_summary.csv.
PRIMARY_LANDMARK_DAYS <- 180L
LANDMARK_DAYS <- c(60L, 90L, 120L, 180L)

# Internal validation
BOOTSTRAP_B <- 1000L
BOOTSTRAP_SEED <- 20260907L
MAX_WORKERS <- 10L

# Fixed prediction horizons AFTER the landmark
PREDICTION_HORIZONS_MONTHS <- c(12, 24, 36)

DAYS_PER_MONTH <- 365.25 / 12
MIN_POSITIVE_MONTHS <- 0.5 / DAYS_PER_MONTH

# ==========================================================
# STAGE 2: MULTIPLE IMPUTATION / MODEL ADJUDICATION
# ==========================================================
MI_M <- 20L
MI_MAXIT <- 20L
MI_SEED <- 20260908L
STAGE2_BOOTSTRAP_B <- 300L
STAGE2_BOOTSTRAP_SEED <- 20260909L

# ==========================================================
# BASIC HELPERS
# ==========================================================
stop_if <- function(cond, msg) if (isTRUE(cond)) stop(msg, call. = FALSE)

p_format <- function(p) {
  dplyr::case_when(
    is.na(p) ~ NA_character_,
    p < 0.001 ~ "<0.001",
    p > 0.9 ~ ">0.9",
    TRUE ~ sprintf("%.3f", p)
  )
}

safe_num <- function(x) suppressWarnings(as.numeric(as.character(x)))
safe_chr <- function(x) stringr::str_trim(as.character(x))

first_existing <- function(nms, candidates, required = FALSE, label = NULL) {
  hit <- candidates[candidates %in% nms]
  if (length(hit) > 0L) return(hit[[1L]])
  if (required) {
    stop(
      paste0(
        "Required field not found", if (!is.null(label)) paste0(" for ", label) else "",
        ". Tried: ", paste(candidates, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  NA_character_
}

NCDB_ALIASES <- list(
  id = c("PUF_CASE_ID", "id"),
  dx_year = c("YEAR_OF_DIAGNOSIS"),
  dx_surg_started = c("DX_SURG_STARTED_DAYS"),
  dx_defsurg_started = c("DX_DEFSURG_STARTED_DAYS"),
  dx_rad_started = c("DX_RAD_STARTED_DAYS"),
  dx_chemo_started = c("DX_CHEMO_STARTED_DAYS"),
  surgery_code_legacy = c("RX_SUMM_SURG_PRIM_SITE"),
  surgery_code_2023 = c("RX_SUMM_SURG_PRIM_SITE_2023"),
  reason_no_surgery = c("REASON_FOR_NO_SURGERY"),
  reason_no_rad = c("REASON_FOR_NO_RADIATION"),
  chemo_summary = c("RX_SUMM_CHEMO"),
  mort30 = c("PUF_30_DAY_MORT_CD"),
  mort90 = c("PUF_90_DAY_MORT_CD")
)

harmonize_brain_surgery_code <- function(dx_year, legacy_code, code_2023) {
  y <- safe_num(dx_year)
  old <- str_pad(str_replace_all(str_to_upper(safe_chr(legacy_code)), "[^0-9]", ""), 2, pad = "0")
  new <- str_to_upper(safe_chr(code_2023))
  new <- ifelse(new == "", NA_character_, new)
  case_when(
    is.finite(y) & y >= 2023 & !is.na(new) ~ new,
    is.finite(y) & y < 2023 & !is.na(old) ~ old,
    !is.na(new) ~ new,
    !is.na(old) ~ old,
    TRUE ~ NA_character_
  )
}

procedure_group_from_brain_code <- function(code) {
  z <- str_to_upper(safe_chr(code))
  case_when(
    z %in% c("20", "A200") ~ "Biopsy/local excision",
    z %in% c("21", "40", "A210", "A400") ~ "Subtotal/partial resection",
    z %in% c("30", "55", "A300", "A550") ~ "Gross-total resection",
    z %in% c("90", "A900") ~ "Surgery NOS",
    z %in% c("00", "A000") ~ "No primary-site surgery",
    z %in% c("10", "A100") ~ "Tumor destruction",
    z %in% c("22", "A220") ~ "Spinal cord/nerve resection",
    z %in% c("99", "A990") ~ "Unknown surgery",
    TRUE ~ NA_character_
  )
}

derive_radiation_received <- function(rad_start_days, reason_no_rad = NA) {
  rr <- safe_chr(reason_no_rad)
  case_when(
    is.finite(rad_start_days) & rad_start_days >= 0 ~ 1L,
    rr == "0" ~ 1L,
    rr %in% c("1", "2", "5", "6", "7") ~ 0L,
    rr %in% c("8", "9") ~ NA_integer_,
    TRUE ~ NA_integer_
  )
}

derive_chemo_received <- function(chemo_start_days, chemo_summary = NA) {
  cs <- safe_num(chemo_summary)
  case_when(
    is.finite(chemo_start_days) & chemo_start_days >= 0 ~ 1L,
    is.finite(cs) & cs %in% c(1, 2, 3) ~ 1L,
    is.finite(cs) & cs %in% c(0, 82, 85, 86, 87) ~ 0L,
    is.finite(cs) & cs %in% c(88, 99) ~ NA_integer_,
    TRUE ~ NA_integer_
  )
}

status_by_landmark <- function(received, rel_start_days, landmark_days) {
  case_when(
    received == 0L ~ "No",
    is.finite(rel_start_days) & rel_start_days <= landmark_days ~ "Yes",
    received == 1L & is.finite(rel_start_days) & rel_start_days > landmark_days ~ "No",
    TRUE ~ NA_character_
  )
}

make_tx_group <- function(rt, chemo) {
  case_when(
    rt == "No"  & chemo == "No"  ~ "Neither",
    rt == "Yes" & chemo == "No"  ~ "RT only",
    rt == "No"  & chemo == "Yes" ~ "Chemotherapy only",
    rt == "Yes" & chemo == "Yes" ~ "RT + chemotherapy",
    TRUE ~ NA_character_
  )
}

make_tx_sequence <- function(rt_rel, chemo_rel, rt_ever, chemo_ever) {
  case_when(
    rt_ever == 0L & chemo_ever == 0L ~ "Neither",
    rt_ever == 1L & chemo_ever == 0L ~ "RT only",
    rt_ever == 0L & chemo_ever == 1L ~ "Chemotherapy only",
    rt_ever == 1L & chemo_ever == 1L & is.finite(rt_rel) & is.finite(chemo_rel) & rt_rel < chemo_rel ~ "RT before chemotherapy",
    rt_ever == 1L & chemo_ever == 1L & is.finite(rt_rel) & is.finite(chemo_rel) & chemo_rel < rt_rel ~ "Chemotherapy before RT",
    rt_ever == 1L & chemo_ever == 1L & is.finite(rt_rel) & is.finite(chemo_rel) & rt_rel == chemo_rel ~ "Same recorded start day",
    rt_ever == 1L & chemo_ever == 1L ~ "Both; sequence unknown",
    TRUE ~ NA_character_
  )
}

harrell_c <- function(d, lp, time_col = "landmark_time_months", event_col = "event") {
  t <- d[[time_col]]
  e <- d[[event_col]]
  ok <- is.finite(lp) & is.finite(t) & !is.na(e)
  if (sum(ok) < 5L || sum(e[ok] == 1L) < 2L) return(NA_real_)
  cc <- survival::concordance(survival::Surv(t[ok], e[ok]) ~ lp[ok], reverse = TRUE)
  as.numeric(cc$concordance)
}

calibration_slope <- function(d, lp, time_col = "landmark_time_months", event_col = "event") {
  t <- d[[time_col]]
  e <- d[[event_col]]
  ok <- is.finite(lp) & is.finite(t) & !is.na(e)
  if (sum(ok) < 10L || stats::sd(lp[ok]) <= 0) return(NA_real_)
  dd <- data.frame(t = t[ok], e = e[ok], lp = lp[ok])
  z <- try(survival::coxph(survival::Surv(t, e) ~ lp, data = dd, ties = "efron"), silent = TRUE)
  if (inherits(z, "try-error")) return(NA_real_)
  as.numeric(stats::coef(z)[[1L]])
}

km_censor_surv_function <- function(d, time_col = "landmark_time_months", event_col = "event") {
  t <- d[[time_col]]
  e <- d[[event_col]]
  fit <- survival::survfit(survival::Surv(t, 1 - e) ~ 1)
  function(x) {
    if (length(x) == 0L) return(numeric(0))
    out <- vapply(x, function(xx) {
      ss <- summary(fit, times = max(0, xx), extend = TRUE)
      if (length(ss$surv) == 0L) 1 else as.numeric(ss$surv[[1L]])
    }, numeric(1))
    pmax(out, 1e-6)
  }
}

ipcw_auc <- function(d, marker, horizon, time_col = "landmark_time_months", event_col = "event") {
  t <- d[[time_col]]
  e <- d[[event_col]]
  ok <- is.finite(marker) & is.finite(t) & !is.na(e)
  t <- t[ok]; e <- e[ok]; marker <- marker[ok]
  cases <- which(t <= horizon & e == 1L)
  controls <- which(t > horizon)
  if (length(cases) < 2L || length(controls) < 2L) return(NA_real_)
  dd <- data.frame(tt = t, ee = e)
  names(dd) <- c(time_col, event_col)
  G <- km_censor_surv_function(dd, time_col, event_col)
  wc <- 1 / G(pmax(t[cases] - 1e-8, 0))
  wctrl <- rep(1 / G(horizon), length(controls))
  denom <- sum(wc) * sum(wctrl)
  if (!is.finite(denom) || denom <= 0) return(NA_real_)
  numer <- 0
  for (ii in seq_along(cases)) {
    cmp <- (marker[cases[ii]] > marker[controls]) + 0.5 * (marker[cases[ii]] == marker[controls])
    numer <- numer + wc[ii] * sum(wctrl * cmp)
  }
  numer / denom
}

ipcw_brier <- function(d, pred_surv, horizon, time_col = "landmark_time_months", event_col = "event") {
  t <- d[[time_col]]
  e <- d[[event_col]]
  ok <- is.finite(pred_surv) & is.finite(t) & !is.na(e)
  t <- t[ok]; e <- e[ok]; pred_surv <- pred_surv[ok]
  if (length(t) == 0L) return(NA_real_)
  dd <- data.frame(tt = t, ee = e)
  names(dd) <- c(time_col, event_col)
  G <- km_censor_surv_function(dd, time_col, event_col)
  y <- as.numeric(t > horizon)
  w <- rep(0, length(t))
  ev <- t <= horizon & e == 1L
  alive <- t > horizon
  if (any(ev)) w[ev] <- 1 / G(pmax(t[ev] - 1e-8, 0))
  if (any(alive)) w[alive] <- 1 / G(horizon)
  mean(w * (y - pred_surv)^2, na.rm = TRUE)
}

baseline_survival_at <- function(fit, horizon) {
  bh <- survival::basehaz(fit, centered = FALSE)
  z <- bh[bh$time <= horizon, , drop = FALSE]
  if (nrow(z) == 0L) return(1)
  exp(-tail(z$hazard, 1L))
}

pred_survival_at <- function(fit, newdata, horizon) {
  s0 <- baseline_survival_at(fit, horizon)
  lp <- as.numeric(predict(fit, newdata = newdata, type = "lp", reference = "zero"))
  s0 ^ exp(lp)
}

bootstrap_optimism_model <- function(d, formula_obj, B = BOOTSTRAP_B, seed = BOOTSTRAP_SEED,
                                     workers = MAX_WORKERS, time_col = "landmark_time_months") {
  fit0 <- survival::coxph(formula_obj, data = d, ties = "efron", x = TRUE, y = TRUE, model = TRUE)
  if (anyNA(coef(fit0))) stop("Non-estimable coefficient(s) in primary model.", call. = FALSE)
  lp0 <- as.numeric(predict(fit0, newdata = d, type = "lp", reference = "zero"))
  c0 <- harrell_c(d, lp0, time_col = time_col)

  one_boot <- function(b, d, formula_obj, time_col) {
    idx <- sample.int(nrow(d), nrow(d), replace = TRUE)
    db <- d[idx, , drop = FALSE]
    fb <- try(survival::coxph(formula_obj, data = db, ties = "efron", x = TRUE, y = TRUE, model = TRUE), silent = TRUE)
    if (inherits(fb, "try-error") || anyNA(coef(fb))) return(c(copt = NA_real_, sopt = NA_real_))
    lpb <- try(as.numeric(predict(fb, newdata = db, type = "lp", reference = "zero")), silent = TRUE)
    lpt <- try(as.numeric(predict(fb, newdata = d, type = "lp", reference = "zero")), silent = TRUE)
    if (inherits(lpb, "try-error") || inherits(lpt, "try-error")) return(c(copt = NA_real_, sopt = NA_real_))
    copt <- harrell_c(db, lpb, time_col = time_col) - harrell_c(d, lpt, time_col = time_col)
    test_slope <- calibration_slope(d, lpt, time_col = time_col)
    c(copt = copt, sopt = 1 - test_slope)
  }

  detected <- suppressWarnings(parallel::detectCores(logical = TRUE))
  if (!is.finite(detected)) detected <- 2L
  workers <- max(1L, min(as.integer(workers), as.integer(detected) - 1L, as.integer(B)))

  if (workers > 1L) {
    cl <- parallel::makeCluster(workers)
    on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
    parallel::clusterEvalQ(cl, suppressPackageStartupMessages(library(survival)))
    parallel::clusterExport(cl, varlist = c("harrell_c", "calibration_slope"), envir = environment())
    parallel::clusterSetRNGStream(cl, iseed = seed)
    z <- parallel::parLapply(cl, seq_len(B), function(b, d, formula_obj, time_col, one_boot) one_boot(b, d, formula_obj, time_col), d = d, formula_obj = formula_obj, time_col = time_col, one_boot = one_boot)
    parallel::stopCluster(cl)
    on.exit(NULL, add = FALSE)
  } else {
    set.seed(seed)
    z <- lapply(seq_len(B), function(b) one_boot(b, d, formula_obj, time_col))
  }

  zz <- do.call(rbind, z)
  copt <- mean(zz[, "copt"], na.rm = TRUE)
  sopt <- mean(zz[, "sopt"], na.rm = TRUE)

  list(
    fit = fit0,
    summary = tibble(
      n = nrow(d),
      deaths = sum(d$event == 1L),
      fitted_parameters = length(coef(fit0)),
      events_per_parameter = sum(d$event == 1L) / length(coef(fit0)),
      apparent_c_index = c0,
      mean_c_optimism = copt,
      optimism_corrected_c_index = c0 - copt,
      optimism_corrected_calibration_slope = 1 - sopt,
      bootstrap_requested = B,
      bootstrap_successful = sum(is.finite(zz[, "copt"]) & is.finite(zz[, "sopt"]))
    ),
    iterations = as.data.frame(zz)
  )
}

make_tv_intervals_one <- function(id, fu_days, event, rt_time, chemo_time) {
  if (!is.finite(fu_days) || fu_days <= 0) return(NULL)
  rt_eff <- if (is.finite(rt_time)) max(0, rt_time) else Inf
  ch_eff <- if (is.finite(chemo_time)) max(0, chemo_time) else Inf
  cuts <- sort(unique(c(0, rt_eff[rt_eff > 0 & rt_eff < fu_days], ch_eff[ch_eff > 0 & ch_eff < fu_days], fu_days)))
  if (length(cuts) < 2L) return(NULL)
  out <- tibble(
    id = id,
    tstart_days = head(cuts, -1),
    tstop_days = tail(cuts, -1)
  ) %>%
    mutate(
      rt_td = as.integer(tstart_days >= rt_eff),
      chemo_td = as.integer(tstart_days >= ch_eff),
      event_tv = as.integer(event == 1L & abs(tstop_days - fu_days) < 1e-8)
    )
  out
}
