############################################################
# 04_OLIGO_TIMEVARYING_treatment_sensitivity.R
# Demonstrate/repair immortal-time bias for adjuvant RT and chemotherapy.
# Compares:
#   A) naive ever-treated exposure assigned at procedure time
#   B) true time-varying exposure switching on at recorded treatment start
############################################################

SCRIPT_DIR <- dirname(normalizePath(sys.frame(1)$ofile, winslash = "/", mustWork = FALSE))
source(file.path(SCRIPT_DIR, "00_OLIGO_LANDMARK_config_helpers.R"))

in_path <- file.path(PROCESSED_DIR, "01_timing_eligible_postprocedure.rds")
out_dir <- file.path(OUTPUT_ROOT, "04_timevarying")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
stopifnot(file.exists(in_path))

d0 <- readRDS(in_path) %>%
  mutate(
    procedure_group = factor(procedure_group, levels = c("Subtotal/partial resection", "Gross-total resection", "Biopsy/local excision")),
    postop_time_months = pmax(postop_followup_days / DAYS_PER_MONTH, MIN_POSITIVE_MONTHS),
    rt_ever = factor(case_when(rad_received_ever == 0L ~ "No", rad_received_ever == 1L ~ "Yes", TRUE ~ NA_character_), levels = c("No", "Yes")),
    chemo_ever = factor(case_when(chemo_received_ever == 0L ~ "No", chemo_received_ever == 1L ~ "Yes", TRUE ~ NA_character_), levels = c("No", "Yes"))
  )

required_baseline <- c("age", "sex", "molecular_grade", "tumor_size_mm", "cdcc", "primary_site_group", "procedure_group")

tv_eligible <- d0 %>%
  filter(
    if_all(all_of(required_baseline), ~ !is.na(.x)),
    !is.na(rt_ever), !is.na(chemo_ever),
    !(rad_received_ever == 1L & !is.finite(rt_rel_days)),
    !(chemo_received_ever == 1L & !is.finite(chemo_rel_days)),
    is.finite(postop_followup_days), postop_followup_days > 0
  ) %>% droplevels()

qc <- tibble(
  metric = c(
    "Timing-eligible base N",
    "Time-varying common-cohort N",
    "Deaths",
    "Excluded: RT received but start missing",
    "Excluded: chemotherapy received but start missing",
    "Included RT starts before procedure",
    "Included chemotherapy starts before procedure"
  ),
  value = c(
    nrow(d0),
    nrow(tv_eligible),
    sum(tv_eligible$event == 1L),
    sum(d0$rad_received_ever == 1L & !is.finite(d0$rt_rel_days), na.rm = TRUE),
    sum(d0$chemo_received_ever == 1L & !is.finite(d0$chemo_rel_days), na.rm = TRUE),
    sum(is.finite(tv_eligible$rt_rel_days) & tv_eligible$rt_rel_days < 0),
    sum(is.finite(tv_eligible$chemo_rel_days) & tv_eligible$chemo_rel_days < 0)
  )
)
write_csv(qc, file.path(out_dir, "01_timevarying_cohort_QC.csv"))

naive_formula <- survival::Surv(postop_time_months, event) ~
  age + sex + molecular_grade + tumor_size_mm + cdcc + primary_site_group + procedure_group + rt_ever + chemo_ever

fit_naive <- survival::coxph(
  naive_formula, data = tv_eligible, ties = "efron", x = TRUE, y = TRUE, model = TRUE
)
saveRDS(fit_naive, file.path(out_dir, "02_naive_ever_treated_model_BIAS_BENCHMARK.rds"))

intervals <- lapply(seq_len(nrow(tv_eligible)), function(i) {
  r <- tv_eligible[i, ]
  rt_time <- if (r$rad_received_ever == 1L) r$rt_rel_days else Inf
  chemo_time <- if (r$chemo_received_ever == 1L) r$chemo_rel_days else Inf
  make_tv_intervals_one(
    id = r$id,
    fu_days = r$postop_followup_days,
    event = r$event,
    rt_time = rt_time,
    chemo_time = chemo_time
  )
})

tv_long <- bind_rows(intervals) %>%
  left_join(tv_eligible %>% select(id, all_of(required_baseline)), by = "id") %>%
  droplevels()

stop_if(nrow(tv_long) == 0L, "No time-varying intervals were created.")
stop_if(any(tv_long$tstop_days <= tv_long$tstart_days), "Invalid time-varying interval detected.")

write_csv(
  tv_long %>% summarise(
    patients = n_distinct(id),
    intervals = n(),
    deaths = sum(event_tv),
    rt_exposed_intervals = sum(rt_td == 1L),
    chemo_exposed_intervals = sum(chemo_td == 1L)
  ),
  file.path(out_dir, "03_timevarying_interval_summary.csv")
)
saveRDS(tv_long, file.path(PROCESSED_DIR, "04_timevarying_long.rds"))

tv_formula <- survival::Surv(tstart_days, tstop_days, event_tv) ~
  age + sex + molecular_grade + tumor_size_mm + cdcc + primary_site_group + procedure_group + rt_td + chemo_td + cluster(id)

fit_tv <- survival::coxph(
  tv_formula, data = tv_long, ties = "efron", x = TRUE, y = TRUE, model = TRUE
)
saveRDS(fit_tv, file.path(out_dir, "04_timevarying_treatment_model.rds"))

naive_coef <- broom::tidy(fit_naive, conf.int = TRUE) %>%
  filter(term %in% c("rt_everYes", "chemo_everYes")) %>%
  mutate(
    exposure = recode(term, "rt_everYes" = "Radiotherapy", "chemo_everYes" = "Chemotherapy"),
    model = "Naive ever-treated at baseline",
    HR = exp(estimate), HR_low = exp(conf.low), HR_high = exp(conf.high),
    p_formatted = p_format(p.value)
  ) %>% select(exposure, model, term, HR, HR_low, HR_high, p.value, p_formatted)

tv_coef <- broom::tidy(fit_tv, conf.int = TRUE) %>%
  filter(term %in% c("rt_td", "chemo_td")) %>%
  mutate(
    exposure = recode(term, "rt_td" = "Radiotherapy", "chemo_td" = "Chemotherapy"),
    model = "Time-varying exposure",
    HR = exp(estimate), HR_low = exp(conf.low), HR_high = exp(conf.high),
    p_formatted = p_format(p.value)
  ) %>% select(exposure, model, term, HR, HR_low, HR_high, p.value, p_formatted)

comparison <- bind_rows(naive_coef, tv_coef)
write_csv(comparison, file.path(out_dir, "05_naive_vs_timevarying_treatment_HR.csv"))

tv_eligible_postonly <- tv_eligible %>%
  filter((!is.finite(rt_rel_days) | rt_rel_days >= 0), (!is.finite(chemo_rel_days) | chemo_rel_days >= 0))

if (nrow(tv_eligible_postonly) >= 100L && sum(tv_eligible_postonly$event == 1L) >= 20L) {
  intervals2 <- lapply(seq_len(nrow(tv_eligible_postonly)), function(i) {
    r <- tv_eligible_postonly[i, ]
    make_tv_intervals_one(
      id = r$id,
      fu_days = r$postop_followup_days,
      event = r$event,
      rt_time = if (r$rad_received_ever == 1L) r$rt_rel_days else Inf,
      chemo_time = if (r$chemo_received_ever == 1L) r$chemo_rel_days else Inf
    )
  })
  tv_long2 <- bind_rows(intervals2) %>%
    left_join(tv_eligible_postonly %>% select(id, all_of(required_baseline)), by = "id") %>% droplevels()
  fit_tv_postonly <- survival::coxph(tv_formula, data = tv_long2, ties = "efron", x = TRUE, y = TRUE, model = TRUE)
  write_csv(
    broom::tidy(fit_tv_postonly, conf.int = TRUE) %>%
      mutate(HR = exp(estimate), HR_low = exp(conf.low), HR_high = exp(conf.high), p_formatted = p_format(p.value)),
    file.path(out_dir, "06_timevarying_excluding_preprocedure_treatment.csv")
  )
}

cat("\n=== TIME-VARYING TREATMENT SENSITIVITY COMPLETE ===\n")
print(qc)
cat("\nNaive vs time-varying treatment coefficients:\n")
print(comparison)
