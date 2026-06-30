suppressPackageStartupMessages({
  library(shiny)
  library(survival)
  library(bslib)
  library(ggplot2)
  library(rms)
})

# ============================================================
# Oligodendroglioma Overall Survival Estimator
# Shiny app for NCDB oligodendroglioma calculator
#
# Expected project structure:
#   TBD calc/
#   ├── data/
#   │   └── processed/
#   │       ├── oligodendroglioma_model_objects.rds
#   │       └── oligodendroglioma_pipeline/
#   │           ├── 01_outputs/03_final_model_dataset_split.rds
#   │           └── 02_outputs/fit_final_training_cox.rds
#   └── shiny/
#       ├── app.R
#       └── www/ohsu_logo.png
#
# The preferred model object is:
#   obj$cph_fit
#   obj$df2
#
# If oligodendroglioma_model_objects.rds is not present, this app will
# fall back to the pipeline outputs above.
# ============================================================

# ----------------------------
# 1. Load model object
# ----------------------------
model_object_path <- file.path("..", "data", "processed", "oligodendroglioma_model_objects.rds")
fit_fallback_path <- file.path("..", "data", "processed", "oligodendroglioma_pipeline", "02_outputs", "fit_final_training_cox.rds")
df_fallback_path  <- file.path("..", "data", "processed", "oligodendroglioma_pipeline", "01_outputs", "03_final_model_dataset_split.rds")

if (file.exists(model_object_path)) {
  obj <- readRDS(model_object_path)
  fit <- obj$cph_fit
  df_ref <- as.data.frame(obj$df2)
} else if (file.exists(fit_fallback_path) && file.exists(df_fallback_path)) {
  fit <- readRDS(fit_fallback_path)
  df_ref <- as.data.frame(readRDS(df_fallback_path))
} else {
  stop(
    paste0(
      "Could not find the model object or fallback model files.\n\n",
      "Expected one of these setups:\n",
      "1) ", normalizePath(model_object_path, winslash = "/", mustWork = FALSE), "\n",
      "or\n",
      "2) ", normalizePath(fit_fallback_path, winslash = "/", mustWork = FALSE), "\n",
      "   ", normalizePath(df_fallback_path, winslash = "/", mustWork = FALSE), "\n\n",
      "Make sure the app folder is inside the project folder and that the data/processed folder is one level above shiny/."
    ),
    call. = FALSE
  )
}

# ----------------------------
# 2. Required model variables
# ----------------------------
required_vars <- c(
  "age_years",
  "tumor_size_harmonized_mm",
  "bmm_grade_group",
  "surgery_extent",
  "radiation_status",
  "chemo_status",
  "sex_cat",
  "race_cat",
  "ethnicity_cat",
  "insurance_cat",
  "income_quartile",
  "education_quartile",
  "charlson_deyo_cat"
)

missing_vars <- setdiff(required_vars, names(df_ref))
if (length(missing_vars) > 0) {
  stop(
    paste0(
      "The reference dataset is missing required variables:\n",
      paste(missing_vars, collapse = ", "),
      "\n\nMake sure df_ref/obj$df2 is the final model dataset from 03_final_model_dataset_split.rds."
    ),
    call. = FALSE
  )
}

# ----------------------------
# 3. Restore factor levels exactly as used in the model pipeline
# ----------------------------
df_ref$bmm_grade_group <- factor(
  df_ref$bmm_grade_group,
  levels = c(
    "Oligodendroglioma, IDH-mutant/1p19q-codeleted",
    "Anaplastic/grade 3 oligodendroglioma, IDH-mutant/1p19q-codeleted"
  )
)
df_ref$surgery_extent <- factor(df_ref$surgery_extent, levels = c("No surgery", "STR", "GTR"))
df_ref$radiation_status <- factor(df_ref$radiation_status, levels = c("No radiation", "Radiation given"))
df_ref$chemo_status <- factor(df_ref$chemo_status, levels = c("No chemotherapy", "Chemotherapy given"))
df_ref$sex_cat <- factor(df_ref$sex_cat, levels = c("Female", "Male"))
df_ref$race_cat <- factor(df_ref$race_cat, levels = c("White", "Black", "Other"))
df_ref$ethnicity_cat <- factor(df_ref$ethnicity_cat, levels = c("Non-Hispanic", "Hispanic"))
df_ref$insurance_cat <- factor(df_ref$insurance_cat, levels = c("Private", "Medicare", "Medicaid", "Not insured"))
df_ref$income_quartile <- factor(df_ref$income_quartile, levels = c("Q4 highest income", "Q3", "Q2", "Q1 lowest income"))
df_ref$education_quartile <- factor(df_ref$education_quartile, levels = c("Q4 highest education", "Q3", "Q2", "Q1 lowest education"))
df_ref$charlson_deyo_cat <- factor(df_ref$charlson_deyo_cat, levels = c("0", "1", "2+"))

model_n <- nrow(df_ref)
model_deaths <- if ("event_death" %in% names(df_ref)) sum(df_ref$event_death == 1, na.rm = TRUE) else NA_integer_

# ----------------------------
# 4. Utility functions
# ----------------------------
safe_factor <- function(val, levels) {
  if (is.null(levels) || length(levels) == 0) return(factor(val))
  if (!val %in% levels) val <- levels[1]
  factor(val, levels = levels)
}

clamp_num <- function(x, lower, upper, default) {
  x <- suppressWarnings(as.numeric(x))
  if (!is.finite(x)) x <- default
  min(max(x, lower), upper)
}

mode_value <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA_character_)
  tab <- sort(table(as.character(x)), decreasing = TRUE)
  names(tab)[1]
}

surv_at <- function(sf, t) {
  s <- summary(sf, times = t, extend = TRUE)$surv
  as.numeric(s[1])
}

median_stats <- function(sf) {
  out <- list(med = NA_real_, lower = NA_real_, upper = NA_real_)

  q <- tryCatch(
    quantile(sf, probs = 0.5, conf.int = TRUE),
    error = function(e) NULL
  )

  if (!is.null(q)) {
    out$med <- suppressWarnings(as.numeric(q$quantile[1]))
    if (!is.null(q$lower)) out$lower <- suppressWarnings(as.numeric(q$lower[1]))
    if (!is.null(q$upper)) out$upper <- suppressWarnings(as.numeric(q$upper[1]))
    return(out)
  }

  if (!is.null(sf$time) && !is.null(sf$surv) && length(sf$time) > 0) {
    idx <- which(sf$surv <= 0.5)[1]
    if (!is.na(idx)) out$med <- as.numeric(sf$time[idx])
  }

  out
}

# ----------------------------
# 5. User-facing labels
# ----------------------------
grade_choices <- c(
  "Grade 2 oligodendroglioma, IDH-mutant and 1p/19q-codeleted" =
    "Oligodendroglioma, IDH-mutant/1p19q-codeleted",
  "Grade 3/anaplastic oligodendroglioma, IDH-mutant and 1p/19q-codeleted" =
    "Anaplastic/grade 3 oligodendroglioma, IDH-mutant/1p19q-codeleted"
)

surgery_choices <- c(
  "No surgery" = "No surgery",
  "Subtotal resection (STR)" = "STR",
  "Gross total resection (GTR)" = "GTR"
)

radiation_choices <- c(
  "No radiation" = "No radiation",
  "Radiation given" = "Radiation given"
)

chemo_choices <- c(
  "No chemotherapy" = "No chemotherapy",
  "Chemotherapy given" = "Chemotherapy given"
)

sex_choices <- c(
  "Female" = "Female",
  "Male" = "Male"
)

race_choices <- c(
  "White" = "White",
  "Black" = "Black",
  "Other race" = "Other"
)

ethnicity_choices <- c(
  "Non-Hispanic" = "Non-Hispanic",
  "Hispanic" = "Hispanic"
)

insurance_choices <- c(
  "Private insurance" = "Private",
  "Medicare" = "Medicare",
  "Medicaid" = "Medicaid",
  "Not insured" = "Not insured"
)

income_choices <- c(
  "Q4, highest area-level income" = "Q4 highest income",
  "Q3" = "Q3",
  "Q2" = "Q2",
  "Q1, lowest area-level income" = "Q1 lowest income"
)

education_choices <- c(
  "Q4, highest area-level education" = "Q4 highest education",
  "Q3" = "Q3",
  "Q2" = "Q2",
  "Q1, lowest area-level education" = "Q1 lowest education"
)

charlson_choices <- c(
  "0" = "0",
  "1" = "1",
  "2 or more" = "2+"
)

# Keep only choices that exist in df_ref.
keep_existing_choices <- function(choices, x) {
  lv <- levels(x)
  choices[choices %in% lv]
}

grade_choices <- keep_existing_choices(grade_choices, df_ref$bmm_grade_group)
surgery_choices <- keep_existing_choices(surgery_choices, df_ref$surgery_extent)
radiation_choices <- keep_existing_choices(radiation_choices, df_ref$radiation_status)
chemo_choices <- keep_existing_choices(chemo_choices, df_ref$chemo_status)
sex_choices <- keep_existing_choices(sex_choices, df_ref$sex_cat)
race_choices <- keep_existing_choices(race_choices, df_ref$race_cat)
ethnicity_choices <- keep_existing_choices(ethnicity_choices, df_ref$ethnicity_cat)
insurance_choices <- keep_existing_choices(insurance_choices, df_ref$insurance_cat)
income_choices <- keep_existing_choices(income_choices, df_ref$income_quartile)
education_choices <- keep_existing_choices(education_choices, df_ref$education_quartile)
charlson_choices <- keep_existing_choices(charlson_choices, df_ref$charlson_deyo_cat)

# Defaults.
age_min <- 18
age_max <- 90
age_default <- round(median(df_ref$age_years, na.rm = TRUE))
age_default <- clamp_num(age_default, age_min, age_max, 50)

tumor_min <- 1
tumor_max <- 200
tumor_default <- round(median(df_ref$tumor_size_harmonized_mm, na.rm = TRUE))
tumor_default <- clamp_num(tumor_default, tumor_min, tumor_max, 40)

select_default <- function(x, choices) {
  mv <- mode_value(x)
  if (!is.na(mv) && mv %in% choices) return(mv)
  unname(choices[1])
}

grade_default <- select_default(df_ref$bmm_grade_group, grade_choices)
surgery_default <- select_default(df_ref$surgery_extent, surgery_choices)
radiation_default <- select_default(df_ref$radiation_status, radiation_choices)
chemo_default <- select_default(df_ref$chemo_status, chemo_choices)
sex_default <- select_default(df_ref$sex_cat, sex_choices)
race_default <- select_default(df_ref$race_cat, race_choices)
ethnicity_default <- select_default(df_ref$ethnicity_cat, ethnicity_choices)
insurance_default <- select_default(df_ref$insurance_cat, insurance_choices)
income_default <- select_default(df_ref$income_quartile, income_choices)
education_default <- select_default(df_ref$education_quartile, education_choices)
charlson_default <- select_default(df_ref$charlson_deyo_cat, charlson_choices)

# Logo helper.
logo_ui <- if (file.exists(file.path("www", "ohsu_logo.png"))) {
  img(src = "ohsu_logo.png", class = "ohsu-logo")
} else {
  div("OHSU", class = "ohsu-logo-fallback")
}

# ----------------------------
# 6. UI
# ----------------------------
ui <- page_fluid(
  theme = bs_theme(
    version = 5,
    bootswatch = "flatly",
    base_font = font_google("Inter"),
    heading_font = font_google("Inter"),
    primary = "#1f4e79",
    bg = "#f4f7fb",
    fg = "#243447"
  ),

  tags$head(
    tags$meta(name = "viewport", content = "width=device-width, initial-scale=1"),
    tags$style(HTML("
      :root {
        --page-max: 1320px;
        --card-radius: 24px;
        --shadow-soft: 0 8px 28px rgba(31, 52, 73, 0.07);
        --border-soft: #e7edf5;
        --text-main: #243447;
        --text-muted: #5b6b7f;
        --bg-soft: #f4f7fb;
        --accent: #1f4e79;
      }

      body {
        background: var(--bg-soft);
      }

      .app-container {
        max-width: var(--page-max);
        margin: 0 auto;
        padding: 24px 22px 36px 22px;
      }

      .app-header {
        background: #ffffff;
        border-radius: 28px;
        padding: clamp(18px, 2.2vw, 30px);
        margin-bottom: 24px;
        box-shadow: var(--shadow-soft);
        border: 1px solid var(--border-soft);
      }

      .header-grid {
        display: grid;
        grid-template-columns: minmax(70px, 96px) 1fr;
        gap: 20px;
        align-items: center;
      }

      .logo-wrap {
        display: flex;
        align-items: center;
        justify-content: center;
      }

      .ohsu-logo {
        width: clamp(58px, 6vw, 92px);
        height: auto;
        display: block;
      }

      .ohsu-logo-fallback {
        width: 88px;
        height: 88px;
        border-radius: 22px;
        display: flex;
        align-items: center;
        justify-content: center;
        background: #1f4e79;
        color: white;
        font-weight: 900;
        letter-spacing: 0.06em;
      }

      .header-title {
        margin: 0 0 8px 0;
        font-weight: 800;
        line-height: 1.04;
        font-size: clamp(2rem, 3.7vw, 3.4rem);
        color: var(--text-main);
        max-width: 980px;
      }

      .ohsu-subtitle {
        color: var(--text-muted);
        margin: 0 0 3px 0;
        font-size: 1.05rem;
      }

      .ohsu-dept {
        color: #738396;
        margin: 0;
        font-size: 0.98rem;
      }

      .input-card, .metric-card, .plot-card, .detail-card {
        background: #ffffff;
        border: 1px solid var(--border-soft) !important;
        border-radius: var(--card-radius) !important;
        box-shadow: var(--shadow-soft);
      }

      .metric-card .card-body,
      .plot-card .card-body,
      .detail-card .card-body {
        padding: 22px;
      }

      .input-card .card-body {
        padding: 18px 18px 16px 18px;
      }

      .sticky-panel {
        position: sticky;
        top: 24px;
        max-height: calc(100vh - 48px);
        overflow-y: auto;
        padding-right: 4px;
        scrollbar-width: thin;
      }

      .sticky-panel::-webkit-scrollbar {
        width: 8px;
      }

      .sticky-panel::-webkit-scrollbar-thumb {
        background: #cbd6e2;
        border-radius: 999px;
      }

      .section-title {
        font-weight: 800;
        color: var(--text-main);
        margin-bottom: 14px;
        line-height: 1.06;
        font-size: clamp(1.55rem, 2vw, 2rem);
      }

      .plot-title {
        font-weight: 800;
        color: var(--text-main);
        margin-bottom: 10px;
        font-size: 1.15rem;
      }

      .form-label {
        font-weight: 650;
        color: #2f4257;
        margin-bottom: 5px;
        font-size: 0.97rem;
      }

      .shiny-input-container {
        margin-bottom: 10px;
      }

      .form-control, .form-select {
        border-radius: 14px !important;
        border: 1px solid #d4dde8 !important;
        min-height: 44px;
        box-shadow: none !important;
      }

      .form-control:focus, .form-select:focus {
        border-color: #8db1d5 !important;
        box-shadow: 0 0 0 0.16rem rgba(31, 78, 121, 0.10) !important;
      }

      .btn-primary {
        background-color: #245789 !important;
        border-color: #245789 !important;
        border-radius: 14px !important;
        font-weight: 750;
        min-height: 46px;
        margin-top: 6px;
      }

      .btn-primary:hover {
        background-color: #1d476f !important;
        border-color: #1d476f !important;
      }

      .metric-grid {
        display: grid;
        grid-template-columns: repeat(3, minmax(0, 1fr));
        gap: 16px;
        margin-bottom: 18px;
      }

      .metric-card {
        min-height: 126px;
      }

      .metric-value {
        font-size: clamp(1.8rem, 2.4vw, 2.45rem);
        line-height: 1;
        font-weight: 800;
        color: var(--accent);
        margin-bottom: 10px;
      }

      .metric-label {
        font-size: 0.96rem;
        color: var(--text-muted);
        line-height: 1.35;
      }

      .detail-card h3 {
        font-size: 1.08rem;
        font-weight: 750;
        color: var(--text-main);
        margin-top: 0;
        margin-bottom: 0.8rem;
      }

      .detail-card ul {
        margin-bottom: 0;
        padding-left: 1.15rem;
      }

      .detail-card li {
        color: #425466;
        margin-bottom: 0.48rem;
        line-height: 1.5;
      }

      .block-gap {
        height: 18px;
      }

      .checkbox {
        margin-top: 6px;
        margin-bottom: 0;
      }

      .plot-card .shiny-plot-output {
        margin-top: 2px;
      }

      .detail-grid {
        display: grid;
        grid-template-columns: repeat(2, minmax(0, 1fr));
        gap: 26px 38px;
      }

      .detail-section {
        min-width: 0;
      }

      .input-help {
        color: #6b7b8f;
        font-size: 0.84rem;
        margin-top: -5px;
        margin-bottom: 10px;
        line-height: 1.32;
      }

      @media (max-width: 1199px) {
        .metric-grid {
          grid-template-columns: repeat(2, minmax(0, 1fr));
        }
        .sticky-panel {
          position: static;
          max-height: none;
          overflow-y: visible;
          padding-right: 0;
        }
      }

      @media (max-width: 767px) {
        .app-container {
          padding: 18px 14px 28px 14px;
        }
        .header-grid {
          grid-template-columns: 1fr;
          gap: 14px;
          text-align: center;
        }
        .logo-wrap {
          justify-content: center;
        }
        .metric-grid,
        .detail-grid {
          grid-template-columns: 1fr;
        }
        .input-card .card-body,
        .metric-card .card-body,
        .plot-card .card-body,
        .detail-card .card-body {
          padding: 18px;
        }
        .plot-card .shiny-plot-output {
          height: 420px !important;
        }
      }
    "))
  ),

  div(
    class = "app-container",

    div(
      class = "app-header",
      div(
        class = "header-grid",
        div(class = "logo-wrap", logo_ui),
        div(
          h1("Oligodendroglioma Overall Survival Estimator", class = "header-title"),
          p("Oregon Health & Science University", class = "ohsu-subtitle"),
          p("Department of Neurological Surgery", class = "ohsu-dept")
        )
      )
    ),

    layout_columns(
      col_widths = c(4, 8),

      div(
        class = "sticky-panel",
        card(
          class = "input-card",
          card_body(
            h2("Patient characteristics", class = "section-title"),

            numericInput("age", "Age at diagnosis (years)", value = age_default, min = age_min, max = age_max, step = 1),

            selectInput(
              "sex",
              "Sex",
              choices = sex_choices,
              selected = sex_default,
              selectize = FALSE
            ),

            selectInput(
              "race",
              "Race",
              choices = race_choices,
              selected = race_default,
              selectize = FALSE
            ),

            selectInput(
              "ethnicity",
              "Ethnicity",
              choices = ethnicity_choices,
              selected = ethnicity_default,
              selectize = FALSE
            ),

            selectInput(
              "insurance",
              "Insurance",
              choices = insurance_choices,
              selected = insurance_default,
              selectize = FALSE
            ),

            selectInput(
              "income",
              "Area-level income quartile",
              choices = income_choices,
              selected = income_default,
              selectize = FALSE
            ),

            selectInput(
              "education",
              "Area-level education quartile",
              choices = education_choices,
              selected = education_default,
              selectize = FALSE
            ),

            selectInput(
              "cdcc",
              "Charlson-Deyo comorbidity score",
              choices = charlson_choices,
              selected = charlson_default,
              selectize = FALSE
            ),

            numericInput(
              "tsize_mm",
              "Tumor size (mm)",
              value = tumor_default,
              min = tumor_min,
              max = tumor_max,
              step = 1
            ),

            selectInput(
              "grade",
              "Tumor grade",
              choices = grade_choices,
              selected = grade_default,
              selectize = FALSE
            ),

            selectInput(
              "surgery",
              "Surgery extent",
              choices = surgery_choices,
              selected = surgery_default,
              selectize = FALSE
            ),

            selectInput(
              "radiation",
              "Radiation therapy",
              choices = radiation_choices,
              selected = radiation_default,
              selectize = FALSE
            ),

            selectInput(
              "chemo",
              "Chemotherapy",
              choices = chemo_choices,
              selected = chemo_default,
              selectize = FALSE
            ),

            actionButton("calc", "Estimate survival", class = "btn-primary w-100"),
            checkboxInput("show_ci", "Show confidence bands", value = FALSE)
          )
        )
      ),

      div(
        div(
          class = "metric-grid",

          card(
            class = "metric-card",
            card_body(
              div(textOutput("s3"), class = "metric-value"),
              div("3-year overall survival", class = "metric-label")
            )
          ),

          card(
            class = "metric-card",
            card_body(
              div(textOutput("s5"), class = "metric-value"),
              div("5-year overall survival", class = "metric-label")
            )
          ),

          card(
            class = "metric-card",
            card_body(
              div(textOutput("med"), class = "metric-value"),
              div(
                tags$span(
                  "Median predicted survival",
                  class = "metric-label",
                  title = "‘Not reached’ means the predicted survival curve does not fall below 50% within available follow-up, so the median time cannot be estimated."
                )
              )
            )
          )
        ),

        div(class = "block-gap"),

        card(
          class = "plot-card",
          card_body(
            h2("Estimated overall survival", class = "plot-title"),
            plotOutput("survplot", height = "560px")
          )
        )
      )
    ),

    div(class = "block-gap"),

    card(
      class = "detail-card",
      card_body(
        h2("Model details, analysis summary, and intended use", class = "section-title"),

        div(
          class = "detail-grid",

          div(
            class = "detail-section",
            h3("Model cohort and intended use"),
            tags$ul(
              tags$li(paste0("Model cohort: n = ", format(model_n, big.mark = ","), " adults with oligodendroglioma, IDH-mutant and 1p/19q-codeleted.")),
              tags$li("Intended use: this tool provides diagnosis-time, population-level overall survival estimates derived from the National Cancer Database."),
              tags$li("It is intended to support clinician-patient discussion and risk stratification and does not replace individualized clinical judgment."),
              tags$li("Predictions should be interpreted in the context of imaging, pathology, surgical planning, and multidisciplinary evaluation.")
            )
          ),

          div(
            class = "detail-section",
            h3("Cohort and variables"),
            tags$ul(
              tags$li("Data source: National Cancer Database (NCDB) Brain Participant User File."),
              tags$li("Study population: adults diagnosed from 2018-2023 with primary brain-site oligodendroglioma histology and NCDB brain molecular marker codes consistent with IDH-mutant, 1p/19q-codeleted disease."),
              tags$li("Outcome: overall survival, measured in months from diagnosis."),
              tags$li("Predictors included in the model: age, tumor size, grade group, surgery extent, radiation therapy, chemotherapy, sex, race, Hispanic ethnicity, insurance, area-level income, area-level education, and Charlson-Deyo comorbidity score.")
            )
          ),

          div(
            class = "detail-section",
            h3("Statistical analysis"),
            tags$ul(
              tags$li("Model type: multivariable Cox proportional hazards regression fit in the derivation cohort."),
              tags$li("Displayed outputs include predicted overall survival at 3 and 5 years, the estimated survival curve through 60 months, and median predicted survival when estimable."),
              tags$li("Input ranges are restricted to clinically plausible values and to avoid registry unknown codes being entered as numeric values."),
              tags$li("Internal validation and held-out validation metrics should match the final manuscript analysis.")
            )
          ),

          div(
            class = "detail-section",
            h3("Performance and interpretation"),
            tags$ul(
              tags$li("Model performance: training C-index 0.779, validation C-index 0.753, and overall C-index 0.771. Bootstrap optimism-corrected training C-index was 0.746."),
              tags$li("Time-dependent validation AUC/Brier estimates were 0.883/0.048 at 12 months, 0.726/0.104 at 36 months, and 0.727/0.134 at 60 months."),
              tags$li("Race was modeled as White, Black, and Other race. Other race includes NCDB race categories other than White or Black, including American Indian/Alaska Native, Asian/Pacific Islander, and Other race, NOS."),
              tags$li("MGMT promoter status is intentionally not included in the primary oligodendroglioma calculator."),
              tags$li("The calculator is based on registry data and does not incorporate postoperative residual tumor volume, recurrence, performance status, detailed radiation dose, treatment sequencing, or longitudinal treatment changes."),
              tags$li("External validation is still needed before broad clinical application.")
            )
          )
        )
      )
    )
  )
)

# ----------------------------
# 7. Server
# ----------------------------
server <- function(input, output, session) {

  observe({
    current_age <- suppressWarnings(as.numeric(input$age))
    if (!is.na(current_age) && current_age > age_max) {
      updateNumericInput(session, "age", value = age_max)
    }
    if (!is.na(current_age) && current_age < age_min) {
      updateNumericInput(session, "age", value = age_min)
    }
  })

  observe({
    current_val <- suppressWarnings(as.numeric(input$tsize_mm))
    if (!is.na(current_val) && current_val > tumor_max) {
      updateNumericInput(session, "tsize_mm", value = tumor_max)
    }
    if (!is.na(current_val) && current_val < tumor_min) {
      updateNumericInput(session, "tsize_mm", value = tumor_min)
    }
  })

  newdata <- eventReactive(input$calc, {
    age_val <- clamp_num(input$age, age_min, age_max, age_default)
    tumor_val <- clamp_num(input$tsize_mm, tumor_min, tumor_max, tumor_default)

    data.frame(
      age_years = age_val,
      tumor_size_harmonized_mm = tumor_val,
      bmm_grade_group = safe_factor(input$grade, levels(df_ref$bmm_grade_group)),
      surgery_extent = safe_factor(input$surgery, levels(df_ref$surgery_extent)),
      radiation_status = safe_factor(input$radiation, levels(df_ref$radiation_status)),
      chemo_status = safe_factor(input$chemo, levels(df_ref$chemo_status)),
      sex_cat = safe_factor(input$sex, levels(df_ref$sex_cat)),
      race_cat = safe_factor(input$race, levels(df_ref$race_cat)),
      ethnicity_cat = safe_factor(input$ethnicity, levels(df_ref$ethnicity_cat)),
      insurance_cat = safe_factor(input$insurance, levels(df_ref$insurance_cat)),
      income_quartile = safe_factor(input$income, levels(df_ref$income_quartile)),
      education_quartile = safe_factor(input$education, levels(df_ref$education_quartile)),
      charlson_deyo_cat = safe_factor(input$cdcc, levels(df_ref$charlson_deyo_cat)),
      check.names = FALSE
    )
  }, ignoreNULL = FALSE)

  surv_obj <- eventReactive(input$calc, {
    req(newdata())
    survfit(fit, newdata = newdata())
  }, ignoreNULL = FALSE)

  output$s3 <- renderText({
    req(surv_obj())
    sprintf("%.1f%%", 100 * surv_at(surv_obj(), 36))
  })

  output$s5 <- renderText({
    req(surv_obj())
    sprintf("%.1f%%", 100 * surv_at(surv_obj(), 60))
  })

  output$med <- renderText({
    req(surv_obj())
    ms <- median_stats(surv_obj())

    if (!is.finite(ms$med)) return("Not reached")

    med_round <- as.integer(round(ms$med))

    if (is.finite(ms$lower) && is.finite(ms$upper)) {
      pm <- as.integer(round((ms$upper - ms$lower) / 2))
      if (is.finite(pm) && pm > 0) {
        return(sprintf("%d \u00B1 %d months", med_round, pm))
      }
    }

    sprintf("%d months", med_round)
  })

  output$survplot <- renderPlot({
    req(surv_obj())
    sf <- surv_obj()

    df_plot <- data.frame(
      time = sf$time,
      surv = sf$surv
    )

    if (!is.null(sf$lower) && !is.null(sf$upper)) {
      df_plot$lower <- sf$lower
      df_plot$upper <- sf$upper
    } else {
      df_plot$lower <- NA_real_
      df_plot$upper <- NA_real_
    }

    df_plot <- df_plot[df_plot$time <= 60, , drop = FALSE]

    if (nrow(df_plot) == 0) {
      df_plot <- data.frame(
        time = c(0, 60),
        surv = c(1, 1),
        lower = c(1, 1),
        upper = c(1, 1)
      )
    } else if (min(df_plot$time) > 0) {
      df_plot <- rbind(
        data.frame(time = 0, surv = 1, lower = 1, upper = 1),
        df_plot
      )
    }

    pts <- c(36, 60)
    s_main <- summary(sf, times = pts, extend = TRUE)

    pts_df <- data.frame(
      time = s_main$time,
      surv = s_main$surv
    )

    guide_df <- data.frame(
      time = pts
    )

    ggplot(df_plot, aes(x = time, y = surv)) +
      {
        if (isTRUE(input$show_ci)) {
          geom_ribbon(
            aes(ymin = lower, ymax = upper),
            fill = "#8fb3d9",
            alpha = 0.22
          )
        }
      } +
      geom_vline(
        data = guide_df,
        aes(xintercept = time),
        linetype = "dashed",
        linewidth = 0.55,
        color = "#cbd6e2"
      ) +
      geom_step(
        color = "#1f6feb",
        linewidth = 1.5,
        direction = "hv"
      ) +
      geom_point(
        data = pts_df,
        aes(x = time, y = surv),
        inherit.aes = FALSE,
        color = "#1f6feb",
        size = 3.2
      ) +
      scale_x_continuous(
        limits = c(0, 61),
        breaks = seq(0, 60, by = 12),
        expand = expansion(mult = c(0.01, 0.02))
      ) +
      scale_y_continuous(
        limits = c(0, 1.04),
        breaks = seq(0, 1, by = 0.2),
        labels = function(x) sprintf("%.1f", x),
        expand = expansion(mult = c(0.01, 0.02))
      ) +
      labs(
        x = "Months",
        y = "Overall survival"
      ) +
      theme_minimal(base_size = 14) +
      theme(
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "#e7edf5", linewidth = 0.7),
        axis.title = element_text(color = "#2f4257", face = "bold"),
        axis.text = element_text(color = "#425466"),
        plot.background = element_rect(fill = "#ffffff", color = NA),
        panel.background = element_rect(fill = "#ffffff", color = NA),
        plot.margin = margin(10, 10, 8, 8)
      )
  }, res = 120)
}

shinyApp(ui, server)
