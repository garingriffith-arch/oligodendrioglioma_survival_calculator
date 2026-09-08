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
# Deployment note:
# - The final Stage 5 deployment coefficients and reference
#   baseline survivals at 12, 24, and 36 months are embedded
#   directly below.
# - No external .rds or model CSV is required for the app to run.
# - This avoids deployment failures caused by omitted model files.
# ============================================================

# ----------------------------
# 1. Frozen final deployment model
# ----------------------------
MODEL_VERSION <- "oligo_landmark_stage5_2026-09-07"

horizons <- c(12L, 24L, 36L)
model_n <- 5064L
model_deaths <- 487L
shrinkage_factor <- 0.9760689638084762

# Uniformly shrunken Stage 5 deployment coefficients.
BETA <- c(
  age = 0.05095601523983125,
  grade3 = 0.6024565026135468,
  tumor_size_mm = 0.002646459974816687,
  cdcc1 = 0.24724068144979028,
  cdcc2plus = 0.5127431011630753,
  str_partial = -0.2693157192625999,
  gtr = -0.4122570643722357,
  male = 0.1671893429159114
)

# Reference baseline survival after re-estimation of the baseline
# hazard with the shrunken linear predictor. Reference profile:
# age 45, Grade 2, tumor size 50 mm, CDCC 0,
# biopsy/local excision, female.
BASELINE_SURVIVAL <- c(
  `12` = 0.9811932057756714,
  `24` = 0.9589375516251721,
  `36` = 0.9363517030566854
)

# Frozen internal-validation values from the final analysis.
corrected_c <- 0.7289661234802536
corrected_auc <- c(
  `12` = 0.7338511997892629,
  `24` = 0.7455903899504883,
  `36` = 0.7415886914652117
)
corrected_brier <- c(
  `12` = 0.02508667534348479,
  `24` = 0.05002006728228309,
  `36` = 0.07273857614532989
)
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

linear_predictor <- function(age, molecular_grade, tumor_size_mm, cdcc, procedure, sex) {
  as.numeric(
    BETA[["age"]] * (age - 45) +
      BETA[["grade3"]] * as.numeric(molecular_grade == "Grade 3") +
      BETA[["tumor_size_mm"]] * (tumor_size_mm - 50) +
      BETA[["cdcc1"]] * as.numeric(cdcc == "1") +
      BETA[["cdcc2plus"]] * as.numeric(cdcc == "2+") +
      BETA[["str_partial"]] * as.numeric(procedure == "Subtotal/partial resection") +
      BETA[["gtr"]] * as.numeric(procedure == "Gross-total resection") +
      BETA[["male"]] * as.numeric(sex == "Male")
  )
}

predict_survival_at <- function(age, molecular_grade, tumor_size_mm, cdcc, procedure, sex, times) {
  times <- suppressWarnings(as.integer(times))
  if (any(!times %in% horizons)) {
    stop("Predictions are available only at 12, 24, and 36 months after the landmark.", call. = FALSE)
  }

  lp <- linear_predictor(
    age = age,
    molecular_grade = molecular_grade,
    tumor_size_mm = tumor_size_mm,
    cdcc = cdcc,
    procedure = procedure,
    sex = sex
  )

  s0 <- unname(BASELINE_SURVIVAL[as.character(times)])
  surv <- s0 ^ exp(lp)
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

model_status_ui <- NULL


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
            h2("Estimated fixed-horizon overall survival after the 180-day landmark", class = "plot-title"),
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
              tags$li("The final Cox model used the full development cohort with 20 multiply imputed datasets and bootstrap internal validation; the deployed calculator uses the uniformly shrunken final coefficients and re-estimated reference baseline survival."),
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
    h <- horizon_data()
    data.frame(
      month = c(0L, h$horizon_months),
      survival = c(1, h$survival)
    )
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
