# Shiny/app.R
# ------------------------------------------------------------
# Oligodendroglioma Overall Survival Estimator (12-, 24-, and 36-month horizons)
# OHSU | Department of Neurological Surgery
#
# Final deployment model:
# - 180-day landmark after index tissue procedure
# - Uniformly shrunken coefficients from Stage 5
# - Re-estimated pooled baseline hazard from Stage 5
#
# REQUIRED:
# - Place logo at: shiny/www/ohsu_logo.png
# - Place final deployment CSVs in: model/
# ------------------------------------------------------------

suppressPackageStartupMessages({
  library(shiny)
  library(survival)
  library(bslib)
})

# -----------------------------
# Robust project-root discovery (restart-safe)
# -----------------------------
find_project_root <- function() {
  wd <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)

  candidates <- c(
    file.path(wd, "model"),
    file.path(wd, "..", "model")
  )

  if (any(dir.exists(candidates))) {
    if (dir.exists(file.path(wd, "model"))) return(wd)
    if (dir.exists(file.path(wd, "..", "model"))) return(normalizePath(file.path(wd, "..")))
  }

  cur <- wd
  for (i in 1:6) {
    cur <- dirname(cur)
    if (dir.exists(file.path(cur, "model"))) return(cur)
  }

  stop("Could not locate project root containing a model/ directory from working directory: ", wd)
}

project_root <- find_project_root()

# -----------------------------
# Load final Stage 5 deployment artifacts
# -----------------------------
BASELINE_PATHS <- c(
  file.path(project_root, "model", "oligo_deployment_baseline_survival_full_curve.csv"),
  file.path(project_root, "model", "02_deployment_baseline_survival_full_curve.csv")
)

COEF_PATHS <- c(
  file.path(project_root, "model", "oligo_deployment_coefficients.csv"),
  file.path(project_root, "model", "04_FINAL_deployment_coefficients.csv")
)

baseline_path <- BASELINE_PATHS[file.exists(BASELINE_PATHS)][1]
coef_path <- COEF_PATHS[file.exists(COEF_PATHS)][1]

if (length(baseline_path) == 0L || is.na(baseline_path) || !file.exists(baseline_path)) {
  stop("No oligodendroglioma baseline-survival deployment file found. Expected one of:\n",
       paste(BASELINE_PATHS, collapse = "\n"))
}

if (length(coef_path) == 0L || is.na(coef_path) || !file.exists(coef_path)) {
  stop("No oligodendroglioma deployment-coefficient file found. Expected one of:\n",
       paste(COEF_PATHS, collapse = "\n"))
}

baseline_curve <- utils::read.csv(
  baseline_path,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

coef_tbl <- utils::read.csv(
  coef_path,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

required_baseline_cols <- c(
  "time_months_after_landmark",
  "cumulative_baseline_hazard_reference"
)

if (!all(required_baseline_cols %in% names(baseline_curve))) {
  stop(
    "Baseline deployment file is missing required columns: ",
    paste(setdiff(required_baseline_cols, names(baseline_curve)), collapse = ", ")
  )
}

required_coef_cols <- c("term", "deployment_beta")

if (!all(required_coef_cols %in% names(coef_tbl))) {
  stop(
    "Coefficient deployment file is missing required columns: ",
    paste(setdiff(required_coef_cols, names(coef_tbl)), collapse = ", ")
  )
}

baseline_curve$time_months_after_landmark <- suppressWarnings(
  as.numeric(baseline_curve$time_months_after_landmark)
)
baseline_curve$cumulative_baseline_hazard_reference <- suppressWarnings(
  as.numeric(baseline_curve$cumulative_baseline_hazard_reference)
)
coef_tbl$deployment_beta <- suppressWarnings(as.numeric(coef_tbl$deployment_beta))

if (any(!is.finite(baseline_curve$time_months_after_landmark)) ||
    any(!is.finite(baseline_curve$cumulative_baseline_hazard_reference))) {
  stop("Baseline deployment file contains non-finite time or hazard values.")
}

baseline_curve <- baseline_curve[
  order(baseline_curve$time_months_after_landmark),
  ,
  drop = FALSE
]

if (max(baseline_curve$time_months_after_landmark) < 36) {
  stop("Baseline deployment curve does not extend through 36 months after the landmark.")
}

get_deployment_beta <- function(term_name) {
  x <- coef_tbl$deployment_beta[coef_tbl$term == term_name]
  if (length(x) != 1L || !is.finite(x)) {
    stop("Could not uniquely resolve deployment beta for term: ", term_name)
  }
  as.numeric(x)
}

# The Stage 5 deployment model was centered at age 45 years and tumor size 50 mm.
REFERENCE_AGE <- 45
REFERENCE_SIZE_MM <- 50

BETA <- c(
  age = get_deployment_beta("age"),
  grade3 = get_deployment_beta("molecular_gradeGrade 3"),
  tumor_size_mm = get_deployment_beta("tumor_size_mm"),
  cdcc1 = get_deployment_beta("cdcc1"),
  cdcc2plus = get_deployment_beta("cdcc2+"),
  str_partial = get_deployment_beta("procedure_groupSubtotal/partial resection"),
  gtr = get_deployment_beta("procedure_groupGross-total resection"),
  male = get_deployment_beta("sexMale")
)

if (any(!is.finite(BETA))) {
  stop("One or more oligodendroglioma deployment coefficients are non-finite.")
}

# -----------------------------
# Final Stage 5 prediction function
# -----------------------------
predict_oligo_survival <- function(
    age,
    grade,
    tumor_size_mm,
    cdcc,
    procedure,
    sex,
    horizons = c(12, 24, 36)
) {
  age <- suppressWarnings(as.numeric(age))
  tumor_size_mm <- suppressWarnings(as.numeric(tumor_size_mm))
  horizons <- suppressWarnings(as.numeric(horizons))

  if (length(age) != 1L || !is.finite(age)) {
    stop("Age must be a finite numeric value.")
  }
  if (age < 18 || age > 100) {
    stop("Age must be between 18 and 100 years.")
  }
  if (length(tumor_size_mm) != 1L || !is.finite(tumor_size_mm)) {
    stop("Tumor size must be a finite numeric value.")
  }
  if (tumor_size_mm <= 0 || tumor_size_mm > 300) {
    stop("Tumor size must be greater than 0 and no more than 300 mm.")
  }
  if (any(!is.finite(horizons)) || any(horizons < 0)) {
    stop("Invalid prediction horizon.")
  }

  grade <- match.arg(as.character(grade), c("Grade 2", "Grade 3"))
  cdcc <- match.arg(as.character(cdcc), c("0", "1", "2+"))
  procedure <- match.arg(
    as.character(procedure),
    c("Biopsy/local excision", "Subtotal/partial resection", "Gross-total resection")
  )
  sex <- match.arg(as.character(sex), c("Female", "Male"))

  lp <- BETA[["age"]] * (age - REFERENCE_AGE) +
    BETA[["grade3"]] * as.numeric(grade == "Grade 3") +
    BETA[["tumor_size_mm"]] * (tumor_size_mm - REFERENCE_SIZE_MM) +
    BETA[["cdcc1"]] * as.numeric(cdcc == "1") +
    BETA[["cdcc2plus"]] * as.numeric(cdcc == "2+") +
    BETA[["str_partial"]] * as.numeric(procedure == "Subtotal/partial resection") +
    BETA[["gtr"]] * as.numeric(procedure == "Gross-total resection") +
    BETA[["male"]] * as.numeric(sex == "Male")

  H0 <- vapply(horizons, function(t) {
    idx <- findInterval(t, baseline_curve$time_months_after_landmark)
    if (idx == 0L) {
      0
    } else {
      baseline_curve$cumulative_baseline_hazard_reference[[idx]]
    }
  }, numeric(1))

  surv <- exp(-H0 * exp(lp))

  if (any(!is.finite(surv)) || any(surv < 0) || any(surv > 1)) {
    stop("Prediction produced an invalid survival probability.")
  }

  data.frame(
    horizon_months = horizons,
    survival = surv
  )
}

# -----------------------------
# Build newdata row matching the final model
# -----------------------------
make_newdata <- function(input) {
  list(
    age = as.numeric(input$age),
    grade = input$grade,
    tumor_size_mm = as.numeric(input$size),
    cdcc = input$cdcc,
    procedure = input$procedure,
    sex = input$sex
  )
}

# -----------------------------
# Survival curve (0–36 months after landmark)
# -----------------------------
predict_curve_0_36 <- function(newdata, t_max = 36, step = 0.25) {
  t_grid <- seq(0, t_max, by = step)

  p <- predict_oligo_survival(
    age = newdata$age,
    grade = newdata$grade,
    tumor_size_mm = newdata$tumor_size_mm,
    cdcc = newdata$cdcc,
    procedure = newdata$procedure,
    sex = newdata$sex,
    horizons = t_grid
  )

  data.frame(
    Month = p$horizon_months,
    Survival = p$survival
  )
}

# -----------------------------
# Choices from final model
# -----------------------------
sex_levels <- c("Female", "Male")
grade_levels <- c("Grade 2", "Grade 3")
cdcc_levels <- c("0", "1", "2+")
procedure_levels <- c(
  "Biopsy/local excision",
  "Subtotal/partial resection",
  "Gross-total resection"
)

# -----------------------------
# UI
# -----------------------------
ui <- fluidPage(
  theme = bs_theme(
    version = 5,
    bootswatch = "flatly",
    base_font = font_google("Roboto"),
    heading_font = font_google("Roboto")
  ),

  # Header: logo + left-aligned text block right next to it
  fluidRow(
    style = "
      background-color: #f8f9fa;
      padding: 18px 15px;
      border-bottom: 1px solid #ddd;
      display: flex;
      align-items: flex-start;
      gap: 10px;
    ",
    column(
      width = 2,
      tags$img(
        src = "ohsu_logo.png",
        style = "height: 95px; width: auto; display: block;"
      )
    ),
    column(
      width = 10,
      tags$div(
        style = "display: flex; flex-direction: column; justify-content: flex-start; align-items: flex-start; text-align: left;",
        tags$h3("Oligodendroglioma Overall Survival Estimator (12-, 24-, and 36-month)", style = "margin: 0 0 6px 0;"),
        tags$h5("Oregon Health & Science University", style = "color: #555; margin: 0;"),
        tags$h6("Department of Neurological Surgery", style = "color: #666; margin: 2px 0 0 0;")
      )
    )
  ),

  tags$div(style = "height: 10px;"),

  sidebarLayout(
    sidebarPanel(
      numericInput("age", "Age (years)", value = 45, min = 18, max = 100),

      selectInput("sex", "Sex", choices = sex_levels, selected = sex_levels[1]),
      selectInput("grade", "Molecular grade", choices = grade_levels, selected = grade_levels[1]),
      selectInput("cdcc", "Charlson–Deyo comorbidity", choices = cdcc_levels, selected = cdcc_levels[1]),
      numericInput("size", "Tumor size (mm)", value = 50, min = 1, max = 300, step = 1),

      tags$hr(),

      selectInput(
        "procedure",
        "Index tissue procedure",
        choices = procedure_levels,
        selected = procedure_levels[1]
      ),

      tags$hr(),

      actionButton("calc", "Calculate", class = "btn-primary"),
      br(), br(),
      downloadButton("download_pred", "Download 12/24/36 predictions (CSV)")
    ),

    mainPanel(
      tags$h4("Predicted overall survival probabilities"),
      tableOutput("pred_table"),
      br(),

      tags$h4("Predicted survival curve (0–36 months after landmark)"),
      plotOutput("surv_plot", height = "320px"),
      br(),

      tags$h4("Disclaimer"),
      tags$p(
        paste0(
          "This tool estimates conditional overall survival for adults with molecularly defined oligodendroglioma (IDH-mutant and 1p/19q-codeleted) who remain alive 180 days after their index tissue procedure. ",
          "The 12-, 24-, and 36-month predictions are measured from that 180-day landmark, not from diagnosis or surgery. ",
          "Predictions are based on National Cancer Database cases diagnosed from 2018 through 2023 and use age, sex, molecular grade, tumor size, Charlson–Deyo comorbidity, and index tissue procedure. ",
          "The model is prognostic and should not be interpreted as estimating the causal benefit of biopsy, resection, radiation, chemotherapy, or other therapy. ",
          "It has undergone internal validation but has not yet undergone independent external validation, and the development cohort predates routine use of vorasidenib. ",
          "This tool is intended to support research and clinical discussion and is not a substitute for individualized clinical judgment."
        )
      )
    )
  )
)

# -----------------------------
# Server
# -----------------------------
server <- function(input, output, session) {

  pred_store <- reactiveVal(NULL)
  curve_store <- reactiveVal(NULL)

  observeEvent(input$calc, {
    tryCatch({
      nd <- make_newdata(input)

      pred_raw <- predict_oligo_survival(
        age = nd$age,
        grade = nd$grade,
        tumor_size_mm = nd$tumor_size_mm,
        cdcc = nd$cdcc,
        procedure = nd$procedure,
        sex = nd$sex,
        horizons = c(12, 24, 36)
      )

      pred_tbl <- data.frame(
        Horizon_months_after_landmark = pred_raw$horizon_months,
        Survival_probability = round(pred_raw$survival, 3)
      )

      curve_df <- predict_curve_0_36(nd, t_max = 36, step = 0.25)

      pred_store(pred_tbl)
      curve_store(curve_df)

      output$pred_table <- renderTable(
        pred_tbl,
        striped = TRUE,
        spacing = "s",
        digits = 3
      )

      output$surv_plot <- renderPlot({
        plot(
          curve_df$Month,
          curve_df$Survival,
          type = "l",
          lwd = 2,
          col = "#2C7FB8",
          xlab = "Months after 180-day landmark",
          ylab = "Survival probability",
          ylim = c(0, 1),
          xlim = c(0, 36)
        )
        grid()
        abline(v = c(12, 24, 36), lty = 3, col = "gray50")
        points(
          pred_tbl$Horizon_months_after_landmark,
          pred_tbl$Survival_probability,
          pch = 16,
          col = "#D95F02"
        )
        legend(
          "topright",
          legend = c("Predicted survival", "12/24/36 mo points"),
          lty = c(1, NA),
          pch = c(NA, 16),
          col = c("#2C7FB8", "#D95F02"),
          bty = "n"
        )
      })

    }, error = function(e) {
      pred_store(NULL)
      curve_store(NULL)
      showNotification(
        paste0("Prediction error: ", conditionMessage(e)),
        type = "error",
        duration = NULL
      )
    })
  })

  output$download_pred <- downloadHandler(
    filename = function() {
      paste0("oligodendroglioma_predictions_12_24_36_", Sys.Date(), ".csv")
    },
    content = function(file) {
      pred <- pred_store()
      if (is.null(pred)) {
        pred <- data.frame(
          Horizon_months_after_landmark = c(12, 24, 36),
          Survival_probability = NA_real_
        )
      }
      write.csv(pred, file, row.names = FALSE)
    }
  )
}

shinyApp(ui, server)
