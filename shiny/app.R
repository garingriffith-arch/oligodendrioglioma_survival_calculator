suppressPackageStartupMessages({
  library(shiny)
  library(bslib)
  library(ggplot2)
})

# ============================================================
# Molecularly Defined Oligodendroglioma 180-Day Landmark
# Overall Survival Estimator
#
# Final deployment model from the 2024 NCDB Brain PUF analysis.
# The visual design intentionally mirrors the current OHSU GBM
# calculator; only disease/model-specific content is changed.
#
# Preferred deployment structure:
#   shiny/
#   ├── app.R
#   ├── data/
#   │   └── processed/
#   │       ├── oligodendroglioma_model_objects.rds
#   │       ├── oligo_deployment_coefficients.csv        # fallback
#   │       └── oligo_deployment_baseline_survival_full_curve.csv # fallback
#   └── www/
#       └── ohsu_logo.png
# ============================================================

# ----------------------------
# 1. Load final deployment model
# ----------------------------
MODEL_VERSION <- "oligo_landmark_stage5_2026-09-07"

model_object_candidates <- c(
  file.path("data", "processed", "oligodendroglioma_model_objects.rds"),
  "oligodendroglioma_model_objects.rds",
  file.path("..", "data", "processed", "oligodendroglioma_model_objects.rds")
)

model_object_path <- model_object_candidates[file.exists(model_object_candidates)][1]

obj <- NULL
model_source <- NULL

if (!is.na(model_object_path)) {
  obj_try <- tryCatch(readRDS(model_object_path), error = function(e) NULL)

  if (!is.null(obj_try) &&
      identical(obj_try$model_version, MODEL_VERSION) &&
      !is.null(obj_try$deployment_coefficients) &&
      !is.null(obj_try$baseline_curve)) {
    obj <- obj_try
    model_source <- "rds"
  }
}

# Safe fallback to the exact frozen Stage 5 deployment CSVs.
# This allows deployment before the optional consolidated RDS is built.
if (is.null(obj)) {
  coeff_candidates <- c(
    file.path("data", "processed", "oligo_deployment_coefficients.csv"),
    file.path("model", "oligo_deployment_coefficients.csv"),
    "oligo_deployment_coefficients.csv"
  )
  baseline_candidates <- c(
    file.path("data", "processed", "oligo_deployment_baseline_survival_full_curve.csv"),
    file.path("model", "oligo_deployment_baseline_survival_full_curve.csv"),
    "oligo_deployment_baseline_survival_full_curve.csv"
  )

  coeff_path <- coeff_candidates[file.exists(coeff_candidates)][1]
  baseline_path <- baseline_candidates[file.exists(baseline_candidates)][1]

  if (is.na(coeff_path) || is.na(baseline_path)) {
    stop(
      paste0(
        "Could not find the final oligodendroglioma deployment model.\n\n",
        "Preferred file:\n",
        "  data/processed/oligodendroglioma_model_objects.rds\n\n",
        "Alternatively provide both frozen Stage 5 CSV files:\n",
        "  data/processed/oligo_deployment_coefficients.csv\n",
        "  data/processed/oligo_deployment_baseline_survival_full_curve.csv"
      ),
      call. = FALSE
    )
  }

  obj <- list(
    model_version = MODEL_VERSION,
    deployment_coefficients = read.csv(coeff_path, stringsAsFactors = FALSE, check.names = FALSE),
    baseline_curve = read.csv(baseline_path, stringsAsFactors = FALSE, check.names = FALSE),
    horizons = c(12L, 24L, 36L),
    model_n = 5064L,
    model_deaths = 487L,
    uniform_shrinkage_factor = 0.9760689638084762
  )
  model_source <- "csv"
}

`%||%` <- function(a, b) if (!is.null(a)) a else b

coef_df <- as.data.frame(obj$deployment_coefficients)
baseline_curve <- as.data.frame(obj$baseline_curve)
horizons <- as.integer(obj$horizons %||% c(12L, 24L, 36L))
horizons <- sort(unique(horizons[is.finite(horizons)]))

required_coef_cols <- c("term", "deployment_beta")
if (!all(required_coef_cols %in% names(coef_df))) {
  stop("Deployment coefficient artifact is missing required columns: term and deployment_beta.", call. = FALSE)
}

required_baseline_cols <- c(
  "time_months_after_landmark",
  "cumulative_baseline_hazard_reference"
)
if (!all(required_baseline_cols %in% names(baseline_curve))) {
  stop(
    "Baseline-survival artifact is missing required time/hazard columns.",
    call. = FALSE
  )
}

coef_df$deployment_beta <- suppressWarnings(as.numeric(coef_df$deployment_beta))
coef_lookup <- stats::setNames(coef_df$deployment_beta, coef_df$term)

required_terms <- c(
  "age",
  "molecular_gradeGrade 3",
  "tumor_size_mm",
  "cdcc1",
  "cdcc2+",
  "procedure_groupSubtotal/partial resection",
  "procedure_groupGross-total resection",
  "sexMale"
)

missing_terms <- setdiff(required_terms, names(coef_lookup))
if (length(missing_terms) > 0) {
  stop(
    paste0(
      "Deployment coefficient artifact is missing required terms: ",
      paste(missing_terms, collapse = ", ")
    ),
    call. = FALSE
  )
}

baseline_curve$time_months_after_landmark <- suppressWarnings(
  as.numeric(baseline_curve$time_months_after_landmark)
)
baseline_curve$cumulative_baseline_hazard_reference <- suppressWarnings(
  as.numeric(baseline_curve$cumulative_baseline_hazard_reference)
)

baseline_curve <- baseline_curve[
  is.finite(baseline_curve$time_months_after_landmark) &
    is.finite(baseline_curve$cumulative_baseline_hazard_reference),
  ,
  drop = FALSE
]
baseline_curve <- baseline_curve[
  order(baseline_curve$time_months_after_landmark),
  ,
  drop = FALSE
]

if (nrow(baseline_curve) == 0) {
  stop("Baseline-survival artifact contains no usable rows.", call. = FALSE)
}

# Guarantee an explicit time-zero baseline.
if (baseline_curve$time_months_after_landmark[1] > 0) {
  baseline_curve <- rbind(
    data.frame(
      time_months_after_landmark = 0,
      cumulative_baseline_hazard_reference = 0
    ),
    baseline_curve[, required_baseline_cols, drop = FALSE]
  )
}

if (!identical(horizons, c(12L, 24L, 36L))) {
  stop(
    paste0(
      "The oligodendroglioma deployment model must use 12-, 24-, and 36-month horizons after the landmark. Found: ",
      paste(horizons, collapse = ", "), "."
    ),
    call. = FALSE
  )
}

model_n <- as.integer(obj$model_n %||% 5064L)
model_deaths <- as.integer(obj$model_deaths %||% 487L)
shrinkage_factor <- as.numeric(obj$uniform_shrinkage_factor %||% 0.9760689638084762)

# Frozen internal-validation values from the final analysis.
corrected_c <- 0.7289661234802536
corrected_auc <- c(`12` = 0.7338511997892629,
                   `24` = 0.7455903899504883,
                   `36` = 0.7415886914652117)
corrected_brier <- c(`12` = 0.02508667534348479,
                     `24` = 0.05002006728228309,
                     `36` = 0.07273857614532989)
pre_shrinkage_corrected_slope <- 0.9760689638084762

# ----------------------------
# 2. Prediction utilities
# ----------------------------
clamp_num <- function(x, lower, upper, default) {
  x <- suppressWarnings(as.numeric(x))
  if (!is.finite(x)) x <- default
  min(max(x, lower), upper)
}

fmt_pct <- function(x) {
  ifelse(is.na(x), "—", sprintf("%.1f%%", 100 * x))
}

baseline_hazard_at <- function(times) {
  times <- as.numeric(times)
  idx <- findInterval(times, baseline_curve$time_months_after_landmark)
  out <- numeric(length(times))
  use <- idx > 0
  out[use] <- baseline_curve$cumulative_baseline_hazard_reference[idx[use]]
  out
}

linear_predictor <- function(age, molecular_grade, tumor_size_mm, cdcc, procedure, sex) {
  as.numeric(
    coef_lookup[["age"]] * (age - 45) +
      coef_lookup[["molecular_gradeGrade 3"]] * as.numeric(molecular_grade == "Grade 3") +
      coef_lookup[["tumor_size_mm"]] * (tumor_size_mm - 50) +
      coef_lookup[["cdcc1"]] * as.numeric(cdcc == "1") +
      coef_lookup[["cdcc2+"]] * as.numeric(cdcc == "2+") +
      coef_lookup[["procedure_groupSubtotal/partial resection"]] * as.numeric(procedure == "Subtotal/partial resection") +
      coef_lookup[["procedure_groupGross-total resection"]] * as.numeric(procedure == "Gross-total resection") +
      coef_lookup[["sexMale"]] * as.numeric(sex == "Male")
  )
}

predict_survival_at <- function(age, molecular_grade, tumor_size_mm, cdcc, procedure, sex, times) {
  lp <- linear_predictor(age, molecular_grade, tumor_size_mm, cdcc, procedure, sex)
  h0 <- baseline_hazard_at(times)
  surv <- exp(-h0 * exp(lp))
  pmin(pmax(as.numeric(surv), 0), 1)
}

# ----------------------------
# 3. Defaults and choices
# ----------------------------
age_min <- 18
age_max <- 90
age_default <- 45

tumor_min <- 1
tumor_max <- 200
tumor_default <- 50

sex_choices <- c("Female" = "Female", "Male" = "Male")
grade_choices <- c("Grade 2" = "Grade 2", "Grade 3" = "Grade 3")
charlson_choices <- c("0" = "0", "1" = "1", "2 or more" = "2+")
procedure_choices <- c(
  "Biopsy / local excision" = "Biopsy/local excision",
  "Subtotal / partial resection" = "Subtotal/partial resection",
  "Gross-total resection" = "Gross-total resection"
)

logo_ui <- if (file.exists(file.path("www", "ohsu_logo.png"))) {
  img(src = "ohsu_logo.png", class = "ohsu-logo")
} else {
  div("OHSU", class = "ohsu-logo-fallback")
}

model_status_ui <- if (!identical(obj$model_version, MODEL_VERSION)) {
  div(
    class = "model-warning",
    tags$strong("Model-object update required. "),
    "The loaded oligodendroglioma model artifact does not match the final landmark deployment model."
  )
} else {
  NULL
}

# ----------------------------
# 4. UI
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
      body { background: var(--bg-soft); }
      .app-container { max-width: var(--page-max); margin: 0 auto; padding: 24px 22px 36px 22px; }
      .app-header { background: #ffffff; border-radius: 28px; padding: clamp(18px, 2.2vw, 30px); margin-bottom: 18px; box-shadow: var(--shadow-soft); border: 1px solid var(--border-soft); }
      .header-grid { display: grid; grid-template-columns: minmax(70px, 96px) 1fr; gap: 20px; align-items: center; }
      .logo-wrap { display: flex; align-items: center; justify-content: center; }
      .ohsu-logo { width: clamp(58px, 6vw, 92px); height: auto; display: block; }
      .ohsu-logo-fallback { width: 88px; height: 88px; border-radius: 22px; display: flex; align-items: center; justify-content: center; background: #1f4e79; color: white; font-weight: 900; letter-spacing: 0.06em; }
      .header-title { margin: 0 0 8px 0; font-weight: 800; line-height: 1.04; font-size: clamp(1.9rem, 3.4vw, 3.2rem); color: var(--text-main); max-width: 1000px; }
      .ohsu-subtitle { color: var(--text-muted); margin: 0 0 3px 0; font-size: 1.05rem; }
      .ohsu-dept { color: #738396; margin: 0; font-size: 0.98rem; }
      .model-warning { max-width: var(--page-max); margin: 0 auto 18px auto; padding: 14px 18px; border-radius: 14px; border: 1px solid #e6b800; background: #fff8d8; color: #594600; }
      .input-card, .metric-card, .plot-card, .detail-card { background: #ffffff; border: 1px solid var(--border-soft) !important; border-radius: var(--card-radius) !important; box-shadow: var(--shadow-soft); }
      .metric-card .card-body, .plot-card .card-body, .detail-card .card-body { padding: 22px; }
      .input-card .card-body { padding: 18px 18px 16px 18px; }
      .sticky-panel { position: sticky; top: 24px; max-height: calc(100vh - 48px); overflow-y: auto; padding-right: 4px; scrollbar-width: thin; }
      .section-title { font-weight: 800; color: var(--text-main); margin-bottom: 14px; line-height: 1.06; font-size: clamp(1.55rem, 2vw, 2rem); }
      .plot-title { font-weight: 800; color: var(--text-main); margin-bottom: 10px; font-size: 1.15rem; }
      .form-label { font-weight: 650; color: #2f4257; margin-bottom: 5px; font-size: 0.97rem; }
      .shiny-input-container { margin-bottom: 10px; }
      .form-control, .form-select { border-radius: 14px !important; border: 1px solid #d4dde8 !important; min-height: 44px; box-shadow: none !important; }
      .btn-primary { background-color: #245789 !important; border-color: #245789 !important; border-radius: 14px !important; font-weight: 750; min-height: 46px; margin-top: 6px; }
      .input-note { margin: 12px 2px 2px 2px; color: #65758a; font-size: 0.86rem; line-height: 1.4; }
      .metric-grid { display: grid; grid-template-columns: repeat(3, minmax(0, 1fr)); gap: 14px; margin-bottom: 18px; }
      .metric-card { min-height: 118px; }
      .metric-value { font-size: clamp(1.45rem, 2vw, 2.05rem); line-height: 1; font-weight: 800; color: var(--accent); margin-bottom: 10px; }
      .metric-label { font-size: 0.92rem; color: var(--text-muted); line-height: 1.35; }
      .detail-card h3 { font-size: 1.08rem; font-weight: 750; color: var(--text-main); margin-top: 0; margin-bottom: 0.8rem; }
      .detail-card ul { margin-bottom: 0; padding-left: 1.15rem; }
      .detail-card li { color: #425466; margin-bottom: 0.48rem; line-height: 1.5; }
      .block-gap { height: 18px; }
      .detail-grid { display: grid; grid-template-columns: repeat(2, minmax(0, 1fr)); gap: 26px 38px; }
      @media (max-width: 1199px) {
        .metric-grid { grid-template-columns: repeat(2, minmax(0, 1fr)); }
        .sticky-panel { position: static; max-height: none; overflow-y: visible; padding-right: 0; }
      }
      @media (max-width: 767px) {
        .app-container { padding: 18px 14px 28px 14px; }
        .header-grid, .metric-grid, .detail-grid { grid-template-columns: 1fr; }
        .header-grid { text-align: center; }
        .plot-card .shiny-plot-output { height: 420px !important; }
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

    model_status_ui,

    layout_columns(
      col_widths = c(4, 8),

      div(
        class = "sticky-panel",
        card(
          class = "input-card",
          card_body(
            h2("Patient and tumor characteristics", class = "section-title"),
            numericInput(
              "age",
              "Age at diagnosis (years)",
              value = age_default,
              min = age_min,
              max = age_max,
              step = 1
            ),
            selectInput("sex", "Sex", choices = sex_choices, selected = "Female", selectize = FALSE),
            selectInput("grade", "Molecular grade", choices = grade_choices, selected = "Grade 2", selectize = FALSE),
            selectInput("cdcc", "Charlson-Deyo comorbidity score", choices = charlson_choices, selected = "0", selectize = FALSE),
            numericInput(
              "tsize_mm",
              "Tumor size (mm)",
              value = tumor_default,
              min = tumor_min,
              max = tumor_max,
              step = 1
            ),
            selectInput(
              "procedure",
              "Index tissue procedure",
              choices = procedure_choices,
              selected = "Gross-total resection",
              selectize = FALSE
            ),
            actionButton("calc", "Estimate overall survival", class = "btn-primary w-100"),
            div(
              class = "input-note",
              "Predictions are conditional on being alive 180 days after the index tissue procedure. Procedure categories are prognostic descriptors and must not be interpreted as causal treatment effects."
            )
          )
        )
      ),

      div(
        div(
          class = "metric-grid",
          lapply(horizons, function(h) {
            card(
              class = "metric-card",
              card_body(
                div(textOutput(paste0("surv_", h)), class = "metric-value"),
                div(paste0(h, "-month overall survival"), class = "metric-label")
              )
            )
          })
        ),
        div(class = "block-gap"),
        card(
          class = "plot-card",
          card_body(
            h2("Estimated overall survival after the 180-day landmark", class = "plot-title"),
            plotOutput("survplot", height = "560px")
          )
        )
      )
    ),

    div(class = "block-gap"),

    card(
      class = "detail-card",
      card_body(
        h2("Model details, validation, and intended use", class = "section-title"),
        div(
          class = "detail-grid",
          div(
            h3("Cohort and intended use"),
            tags$ul(
              tags$li(
                paste0(
                  "The final landmark cohort included ",
                  format(model_n, big.mark = ","),
                  " adults with molecularly defined oligodendroglioma who were alive and under follow-up 180 days after their index tissue procedure; ",
                  format(model_deaths, big.mark = ","),
                  " subsequent deaths were observed."
                )
              ),
              tags$li("The calculator estimates all-cause overall survival at 12, 24, and 36 months after the 180-day landmark, not from diagnosis."),
              tags$li("It is intended to support prognostic counseling and risk communication, not to replace multidisciplinary clinical judgment.")
            )
          ),
          div(
            h3("Data and predictors"),
            tags$ul(
              tags$li("Data source: 2024 National Cancer Database Brain Participant User File; diagnoses from 2018 through 2023."),
              tags$li("Study population: adults with molecularly defined IDH-mutant, 1p/19q-codeleted oligodendroglioma, CNS WHO grade 2 or 3."),
              tags$li("Predictors are age, sex, molecular grade, tumor size, Charlson-Deyo comorbidity score, and index tissue procedure."),
              tags$li("Tumor-size missingness was handled using multiple imputation. Race, ethnicity, and area-level socioeconomic measures were not included in the prediction equation.")
            )
          ),
          div(
            h3("Internal validation"),
            tags$ul(
              tags$li("The final Cox model used the full development cohort with 20 multiply imputed datasets and bootstrap internal validation."),
              tags$li(paste0("Optimism-corrected Harrell C was ", sprintf("%.3f", corrected_c), ".")),
              tags$li(
                paste0(
                  "Optimism-corrected time-dependent AUCs at 12, 24, and 36 months were ",
                  paste(sprintf("%.3f", corrected_auc), collapse = ", "),
                  "; corresponding Brier scores were ",
                  paste(sprintf("%.3f", corrected_brier), collapse = ", "), "."
                )
              ),
              tags$li(
                paste0(
                  "The optimism-corrected calibration slope was ",
                  sprintf("%.3f", pre_shrinkage_corrected_slope),
                  ", which was applied as a uniform shrinkage factor before re-estimating the deployment baseline hazard."
                )
              )
            )
          ),
          div(
            h3("Interpretation and limitations"),
            tags$ul(
              tags$li("The procedure coefficient is prognostic and must not be interpreted as the causal survival benefit of subtotal or gross-total resection."),
              tags$li("The registry does not capture several clinically important factors, including performance status, neurologic deficits, postoperative residual tumor volume, recurrence, and longitudinal treatment changes."),
              tags$li("The development cohort predates routine use of vorasidenib; absolute risk may require recalibration as contemporary treatment patterns evolve, particularly for grade 2 disease."),
              tags$li("Independent external validation and, where necessary, recalibration are required before routine clinical implementation.")
            )
          )
        )
      )
    )
  )
)

# ----------------------------
# 5. Server
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
    current_size <- suppressWarnings(as.numeric(input$tsize_mm))
    if (!is.na(current_size) && current_size > tumor_max) {
      updateNumericInput(session, "tsize_mm", value = tumor_max)
    }
    if (!is.na(current_size) && current_size < tumor_min) {
      updateNumericInput(session, "tsize_mm", value = tumor_min)
    }
  })

  patient_values <- eventReactive(input$calc, {
    list(
      age = clamp_num(input$age, age_min, age_max, age_default),
      sex = as.character(input$sex),
      molecular_grade = as.character(input$grade),
      tumor_size_mm = clamp_num(input$tsize_mm, tumor_min, tumor_max, tumor_default),
      cdcc = as.character(input$cdcc),
      procedure = as.character(input$procedure)
    )
  }, ignoreNULL = FALSE)

  horizon_data <- eventReactive(input$calc, {
    p <- patient_values()
    surv <- predict_survival_at(
      age = p$age,
      molecular_grade = p$molecular_grade,
      tumor_size_mm = p$tumor_size_mm,
      cdcc = p$cdcc,
      procedure = p$procedure,
      sex = p$sex,
      times = horizons
    )

    if (any(!is.finite(surv))) {
      showNotification("Unable to calculate survival for the selected values.", type = "error")
    }

    data.frame(
      horizon_months = horizons,
      survival = as.numeric(surv)
    )
  }, ignoreNULL = FALSE)

  curve_data <- eventReactive(input$calc, {
    p <- patient_values()
    time_grid <- seq(0, 36, by = 0.25)
    surv <- predict_survival_at(
      age = p$age,
      molecular_grade = p$molecular_grade,
      tumor_size_mm = p$tumor_size_mm,
      cdcc = p$cdcc,
      procedure = p$procedure,
      sex = p$sex,
      times = time_grid
    )
    data.frame(month = time_grid, survival = surv)
  }, ignoreNULL = FALSE)

  lapply(horizons, function(h) {
    output[[paste0("surv_", h)]] <- renderText({
      d <- horizon_data()
      fmt_pct(d$survival[d$horizon_months == h][1])
    })
  })

  output$survplot <- renderPlot({
    d <- curve_data()
    h <- horizon_data()

    ggplot(d, aes(x = month, y = survival)) +
      geom_line(linewidth = 1.4, color = "#1f6feb") +
      geom_point(
        data = data.frame(month = h$horizon_months, survival = h$survival),
        aes(x = month, y = survival),
        inherit.aes = FALSE,
        size = 3.2,
        color = "#1f6feb"
      ) +
      scale_x_continuous(
        breaks = c(0, horizons),
        limits = c(0, 36)
      ) +
      scale_y_continuous(
        limits = c(0, 1),
        breaks = seq(0, 1, by = 0.2),
        labels = function(x) paste0(round(100 * x), "%")
      ) +
      labs(
        x = "Months after 180-day landmark",
        y = "Predicted overall survival"
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
