############################################################
# 09_OLIGO_LANDMARK_final_cohort_imputation.R
# Cleaned public version of the final Stage 3B cohort / MICE logic.
############################################################

SCRIPT_DIR <- dirname(normalizePath(sys.frame(1)$ofile, winslash = "/", mustWork = FALSE))
source(file.path(SCRIPT_DIR, "00_OLIGO_LANDMARK_config_helpers.R"))

if (!requireNamespace("mice", quietly = TRUE)) stop("Package 'mice' is required.", call. = FALSE)

L <- 180L
out_dir <- file.path(OUTPUT_ROOT, "07B_stage3B_expanded_biopsy")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

molecular <- readRDS(MOLECULAR_RDS)
id_col <- first_existing(names(molecular), c("id", "PUF_CASE_ID"), required = TRUE, label = "patient ID")
molecular <- molecular %>% mutate(id = as.character(.data[[id_col]]))

required_processed <- c("time_months", "event", "age", "sex", "molecular_grade", "tumor_size_mm", "cdcc")
missing_processed <- setdiff(required_processed, names(molecular))
if (length(missing_processed)) stop("Molecular cohort is missing: ", paste(missing_processed, collapse = ", "), call. = FALSE)

raw_names <- names(data.table::fread(RAW_NCDB_CSV, nrows = 0, showProgress = FALSE))
aliases <- list(
  id = c("PUF_CASE_ID", "id"),
  dx_year = c("YEAR_OF_DIAGNOSIS"),
  dx_surg_started = c("DX_SURG_STARTED_DAYS"),
  dx_defsurg_started = c("DX_DEFSURG_STARTED_DAYS"),
  surgery_code_legacy = c("RX_SUMM_SURG_PRIM_SITE"),
  surgery_code_2023 = c("RX_SUMM_SURG_PRIM_SITE_2023"),
  dxstg_proc = c("RX_SUMM_DXSTG_PROC", "RX_SUMM_DX_STG_PROC"),
  dxstg_proc_days = c("DX_STAGING_PROC_DAYS")
)
raw_map <- lapply(names(aliases), function(k) first_existing(raw_names, aliases[[k]], required = k == "id", label = k))
names(raw_map) <- names(aliases)
selected <- unique(na.omit(unlist(raw_map, use.names = FALSE)))
raw <- data.table::fread(RAW_NCDB_CSV, select = selected, colClasses = "character", showProgress = TRUE) %>% as_tibble()

get_raw <- function(k) {
  nm <- raw_map[[k]]
  if (is.na(nm) || !nm %in% names(raw)) rep(NA_character_, nrow(raw)) else raw[[nm]]
}
norm2 <- function(x) {
  z <- safe_chr(x); z[z == ""] <- NA_character_
  n <- suppressWarnings(as.integer(z))
  ifelse(!is.na(n), sprintf("%02d", n), z)
}

raw_clean <- tibble(
  id = as.character(get_raw("id")),
  dx_year_raw = safe_num(get_raw("dx_year")),
  first_surgery_days_raw = safe_num(get_raw("dx_surg_started")),
  definitive_surgery_days_raw = safe_num(get_raw("dx_defsurg_started")),
  surgery_code_legacy_raw = safe_chr(get_raw("surgery_code_legacy")),
  surgery_code_2023_raw = safe_chr(get_raw("surgery_code_2023")),
  dxstg_proc_code = norm2(get_raw("dxstg_proc")),
  dxstg_proc_days_raw = safe_num(get_raw("dxstg_proc_days"))
) %>% distinct(id, .keep_all = TRUE)

existing_dx_year <- if ("dx_year" %in% names(molecular)) safe_num(molecular$dx_year) else rep(NA_real_, nrow(molecular))

dat <- molecular %>%
  left_join(raw_clean, by = "id") %>%
  mutate(
    dx_year_final = coalesce(dx_year_raw, existing_dx_year),
    surgery_code_harmonized = harmonize_brain_surgery_code(dx_year_final, surgery_code_legacy_raw, surgery_code_2023_raw),
    surgery_group_raw = procedure_group_from_brain_code(surgery_code_harmonized),
    surgery_index_days = coalesce(definitive_surgery_days_raw, first_surgery_days_raw),
    diagnostic_primary_site_biopsy_only =
      surgery_group_raw == "No primary-site surgery" &
      dxstg_proc_code == "02" & is.finite(dxstg_proc_days_raw) & dxstg_proc_days_raw >= 0,
    surgery_coded_index_valid =
      surgery_group_raw %in% c("Biopsy/local excision", "Subtotal/partial resection", "Gross-total resection") &
      is.finite(surgery_index_days) & surgery_index_days >= 0,
    procedure_group = case_when(
      surgery_group_raw == "Gross-total resection" & surgery_coded_index_valid ~ "Gross-total resection",
      surgery_group_raw == "Subtotal/partial resection" & surgery_coded_index_valid ~ "Subtotal/partial resection",
      surgery_group_raw == "Biopsy/local excision" & surgery_coded_index_valid ~ "Biopsy/local excision",
      diagnostic_primary_site_biopsy_only ~ "Biopsy/local excision",
      TRUE ~ NA_character_
    ),
    biopsy_mechanism = case_when(
      surgery_group_raw == "Biopsy/local excision" & surgery_coded_index_valid ~ "Surgery-coded biopsy/local excision",
      diagnostic_primary_site_biopsy_only ~ "Diagnostic primary-site biopsy only",
      TRUE ~ NA_character_
    ),
    index_tissue_procedure_days = case_when(
      surgery_coded_index_valid ~ surgery_index_days,
      diagnostic_primary_site_biopsy_only ~ dxstg_proc_days_raw,
      TRUE ~ NA_real_
    ),
    followup_days_from_diagnosis = safe_num(time_months) * DAYS_PER_MONTH,
    postprocedure_followup_days = followup_days_from_diagnosis - index_tissue_procedure_days,
    landmark180_eligible = !is.na(procedure_group) & is.finite(postprocedure_followup_days) & postprocedure_followup_days > L,
    landmark_time_months = (postprocedure_followup_days - L) / DAYS_PER_MONTH
  )

expanded <- dat %>%
  filter(landmark180_eligible) %>%
  mutate(
    procedure_group = factor(procedure_group, levels = c("Biopsy/local excision", "Subtotal/partial resection", "Gross-total resection")),
    sex = factor(sex),
    molecular_grade = factor(molecular_grade),
    cdcc = factor(cdcc)
  ) %>%
  filter(!is.na(age), !is.na(sex), !is.na(molecular_grade), !is.na(cdcc), !is.na(procedure_group), !is.na(event), is.finite(landmark_time_months)) %>%
  droplevels()

stop_if(nrow(expanded) != 5064L, paste0("Frozen-cohort QC failed: expected N=5064, found ", nrow(expanded)))
stop_if(sum(expanded$event == 1L) != 487L, paste0("Frozen-cohort QC failed: expected 487 deaths, found ", sum(expanded$event == 1L)))
stop_if(sum(expanded$procedure_group == "Biopsy/local excision") != 849L, "Biopsy/local excision count mismatch.")
stop_if(sum(expanded$procedure_group == "Subtotal/partial resection") != 1736L, "STR/partial count mismatch.")
stop_if(sum(expanded$procedure_group == "Gross-total resection") != 2479L, "GTR count mismatch.")

fit0 <- survival::coxph(survival::Surv(landmark_time_months, event) ~ 1, data = expanded, ties = "efron")
bh <- survival::basehaz(fit0, centered = FALSE)
idx <- findInterval(expanded$landmark_time_months, bh$time)
expanded$nelson_aalen <- 0
expanded$nelson_aalen[idx > 0] <- bh$hazard[idx[idx > 0]]

mi_dat <- expanded %>%
  select(landmark_time_months, event, age, sex, molecular_grade, tumor_size_mm, cdcc,
         procedure_group, dx_year_final, nelson_aalen, biopsy_mechanism)

ini <- mice::mice(mi_dat, maxit = 0, printFlag = FALSE)
meth <- ini$method
pred <- ini$predictorMatrix
meth[] <- ""
meth["tumor_size_mm"] <- "pmm"
pred[,] <- 0

candidate_predictors <- c("landmark_time_months", "event", "age", "sex", "molecular_grade",
                          "cdcc", "procedure_group", "dx_year_final", "nelson_aalen")
usable <- candidate_predictors[vapply(mi_dat[candidate_predictors], function(x) {
  z <- x[!is.na(x)]
  length(z) > 1L && length(unique(z)) > 1L
}, logical(1))]
pred["tumor_size_mm", usable] <- 1

set.seed(MI_SEED)
imp <- mice::mice(mi_dat, m = 20, maxit = 20, method = meth, predictorMatrix = pred,
                  seed = MI_SEED, printFlag = TRUE)

saveRDS(imp, file.path(out_dir, "06_mice_expanded_biopsy_cohort.rds"))
readr::write_csv(
  tibble(status = "COMPLETE", analysis_n = nrow(expanded), deaths = sum(expanded$event == 1L),
         biopsy_n = sum(expanded$procedure_group == "Biopsy/local excision"),
         STR_n = sum(expanded$procedure_group == "Subtotal/partial resection"),
         GTR_n = sum(expanded$procedure_group == "Gross-total resection")),
  file.path(out_dir, "99_STAGE3B_COMPLETE.csv")
)
