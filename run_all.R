# ============================================================
# Master reproducibility workflow
#
# Usage from the repository root:
#
#   Rscript --vanilla run_all.R
#   Rscript --vanilla run_all.R --clean
#   Rscript --vanilla run_all.R --only=figure2
#   Rscript --vanilla run_all.R --only=figure4
#   Rscript --vanilla run_all.R --only=kriging
#
# Reproducible input boundary:
#   - The public workflow begins with Data/LME data/a1.csv.gz,
#     ..., a13.csv.gz.
#   - Data/Get_data/GetLME.R documents the optional conversion
#     from the original NetCDF files, but is not run here because
#     those upstream NetCDF files are not redistributed.
#
# Design:
#   - Every analysis script is launched in a separate clean R
#     session using Rscript --vanilla.
#   - The workflow stops immediately when a script fails.
#   - Standard output and errors are written to Output/Logs/.
#   - A machine-readable run manifest is written after every step.
#
# IMPORTANT:
#   Review the `workflow` object below and make its paths match
#   the final repository. Add any computational figure/table
#   scripts that are not yet listed.
# ============================================================

# ------------------------------------------------------------
# 1. Locate repository root without depending on the here package
# ------------------------------------------------------------

command_line <- commandArgs(
  trailingOnly = FALSE
)

file_argument <- grep(
  "^--file=",
  command_line,
  value = TRUE
)

if (length(file_argument) != 1L) {
  stop(
    "run_all.R must be executed with Rscript."
  )
}

run_all_file <- sub(
  "^--file=",
  "",
  file_argument
)

repository_root <- dirname(
  normalizePath(
    run_all_file,
    winslash = "/",
    mustWork = TRUE
  )
)

setwd(
  repository_root
)

options(
  stringsAsFactors = FALSE
)

Sys.setenv(
  TZ = "UTC"
)

# ------------------------------------------------------------
# 2. Parse command-line options
# ------------------------------------------------------------

arguments <- commandArgs(
  trailingOnly = TRUE
)

clean_run <- "--clean" %in% arguments

only_argument <- grep(
  "^--only=",
  arguments,
  value = TRUE
)

only_step <- if (length(only_argument) == 1L) {
  sub(
    "^--only=",
    "",
    only_argument
  )
} else {
  ""
}

unknown_arguments <- arguments[
  !arguments %in% "--clean" &
    !grepl(
      "^--only=",
      arguments
    )
]

if (length(unknown_arguments) > 0L) {
  stop(
    "Unknown argument(s): ",
    paste(
      unknown_arguments,
      collapse = ", "
    )
  )
}

# ------------------------------------------------------------
# 3. Optional cleanup
# ------------------------------------------------------------

remove_generated_validation_files <- function() {

  validation_directory <- file.path(
    repository_root,
    "Data",
    "Valid"
  )

  if (!dir.exists(
    validation_directory
  )) {
    return(
      invisible(
        NULL
      )
    )
  }

  generated_patterns <- c(
    "^assimilation_(Beijing|Shanghai|HongKong)\\.csv$",
    "^assimilated_posterior_lme_grid_(mean|variance)\\.csv$"
  )

  validation_files <- list.files(
    validation_directory,
    full.names = TRUE
  )

  generated_files <- validation_files[
    vapply(
      basename(
        validation_files
      ),
      function(filename) {
        any(
          grepl(
            generated_patterns,
            filename
          )
        )
      },
      logical(1)
    )
  ]

  if (length(generated_files) > 0L) {
    unlink(
      generated_files,
      force = TRUE
    )
  }

  invisible(
    NULL
  )
}

if (clean_run) {

  output_directory <- file.path(
    repository_root,
    "Output"
  )

  if (dir.exists(
    output_directory
  )) {
    message(
      "Removing generated Output directory: ",
      output_directory
    )

    unlink(
      output_directory,
      recursive = TRUE,
      force = TRUE
    )
  }

  remove_generated_validation_files()
}

log_directory <- file.path(
  repository_root,
  "Output",
  "Logs"
)

dir.create(
  log_directory,
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# 4. Workflow definition
#
# Each step contains:
#   id       : compact command-line identifier
#   script   : path relative to repository root
#   expected : key outputs that must exist after successful run
#
# Keep this list in dependency order.
# ------------------------------------------------------------

workflow <- list(
  list(
    id = "figure2",
    script = file.path(
      "Code",
      "Figure2.R"
    ),
    expected = c(
      file.path(
        "Output",
        "Figure2",
        "Figure2(a).png"
      ),
      file.path(
        "Output",
        "Figure2",
        "Figure2(b).png"
      ),
      file.path(
        "Output",
        "Figure2",
        "Figure2(c).png"
      )
    )
  ),
  list(
    id = "prepare_lme",
    script = file.path(
      "Code",
      "DataPreparation",
      "prepare_lme_annual.R"
    ),
    expected = c(
      file.path(
        "Output",
        "Intermediate",
        "LME",
        "lme_annual_1368_1911.rds"
      ),
      file.path(
        "Output",
        "Intermediate",
        "LME",
        "lme_city3_annual_1368_1911.csv"
      )
    )
  ),
  list(
    id = "calibration",
    script = file.path(
      "Code",
      "prepare_calibration.R"
    ),
    expected = c(
      file.path(
        "Output",
        "Intermediate",
        "calibration_parameters.rds"
      )
    )
  ),
  list(
    id = "figure4",
    script = file.path(
      "Code",
      "Figure4.R"
    ),
    expected = c(
      file.path(
        "Output",
        "Figure4",
        "Figure4(a).png"
      ),
      file.path(
        "Output",
        "Figure4",
        "Figure4(b).png"
      ),
      file.path(
        "Output",
        "Figure4",
        "Figure4(c).png"
      )
    )
  ),
  list(
    id = "kriging",
    script = file.path(
      "Code",
      "Get_tempe_all_data.R"
    ),
    expected = c(
      file.path(
        "Output",
        "Intermediate",
        "REACHES",
        "reaches_kriging_lme_grid_mean.csv"
      ),
      file.path(
        "Output",
        "Intermediate",
        "REACHES",
        "reaches_kriging_city3_mean.csv"
      ),
      file.path(
        "Output",
        "Intermediate",
        "REACHES",
        "reaches_kriging_city3_sd.csv"
      )
    )
  ),
  list(
    id = "figure3",
    script = file.path(
      "Code",
      "Figure3.R"
    ),
    expected = character(
      0
    )
  ),
  list(
    id = "figure5",
    script = file.path(
      "Code",
      "Figure5.R"
    ),
    expected = character(
      0
    )
  ),
  list(
    id = "figure6",
    script = file.path(
      "Code",
      "Figure6.R"
    ),
    expected = c(
      file.path(
        "Output",
        "Figure6",
        "Figure6_cluster_assignments.csv"
      )
    )
  ),
  list(
    id = "figure7e",
    script = file.path(
      "Code",
      "Figure7e.R"
    ),
    expected = character(
      0
    )
  ),
  list(
    id = "figure7_8",
    script = file.path(
      "Code",
      "Figure7-8.R"
    ),
    expected = character(
      0
    )
  ),
  list(
    id = "prior",
    script = file.path(
      "Code",
      "Prior.R"
    ),
    expected = c(
      file.path(
        "Output",
        "Intermediate",
        "Prior",
        "mtB.csv"
      ),
      file.path(
        "Output",
        "Intermediate",
        "Prior",
        "muB.csv"
      ),
      file.path(
        "Output",
        "Intermediate",
        "Prior",
        "rtB.csv"
      ),
      file.path(
        "Output",
        "Intermediate",
        "Prior",
        "mtS.csv"
      ),
      file.path(
        "Output",
        "Intermediate",
        "Prior",
        "mtH.csv"
      )
    )
  ),
  list(
    id = "figure9abc",
    script = file.path(
      "Code",
      "Figure9abc.R"
    ),
    expected = character(
      0
    )
  ),
  list(
    id = "figure9d",
    script = file.path(
      "Code",
      "Figure9d.R"
    ),
    expected = c(
      file.path(
        "Data",
        "Valid",
        "assimilation_Beijing.csv"
      ),
      file.path(
        "Data",
        "Valid",
        "assimilation_Shanghai.csv"
      ),
      file.path(
        "Data",
        "Valid",
        "assimilation_HongKong.csv"
      ),
      file.path(
        "Output",
        "Intermediate",
        "Assimilation",
        "assimilation_metrics.csv"
      )
    )
  ),
  list(
    id = "figure10",
    script = file.path(
      "Code",
      "Figure10.R"
    ),
    expected = c(
      file.path(
        "Output",
        "Figure10",
        "Figure10_metrics.csv"
      )
    )
  ),
  list(
    id = "figureS1",
    script = file.path(
      "Code",
      "Supplementary",
      "FigureS1.R"
    ),
    expected = character(
      0
    )
  ),
  list(
    id = "figureS2",
    script = file.path(
      "Code",
      "Supplementary",
      "FigureS2.R"
    ),
    expected = c(
      file.path(
        "Output",
        "Supplementary",
        "FigureS2_plot_data.csv"
      )
    )
  ),
  list(
    id = "grid_assimilation",
    script = file.path(
      "Code",
      "Assimilation_grid.R"
    ),
    expected = c(
      file.path(
        "Output",
        "Intermediate",
        "Assimilation",
        "assimilated_posterior_lme_grid_mean.csv"
      )
    )
  ),
  list(
    id = "figureS5",
    script = file.path(
      "Code",
      "Supplementary",
      "FigureS5.R"
    ),
    expected = c(
      file.path(
        "Output",
        "Intermediate",
        "FigureS5",
        "FigureS5_panel_files.csv"
      )
    )
  ),
  list(
    id = "coverage_population",
    script = file.path(
      "Code",
      "Analysis",
      "coverage_population_GAM.R"
    ),
    expected = c(
      file.path(
        "Output",
        "Tables",
        "Table2_coverage_GAM.csv"
      ),
      file.path(
        "Output",
        "Tables",
        "TableS2_M4_parametric_coefficients.csv"
      )
    )
  ),
  list(
    id = "internal_validation",
    script = file.path(
      "Code",
      "Analysis",
      "interval_compatibility.R"
    ),
    expected = c(
      file.path(
        "Output",
        "Validation",
        "interval_compatibility_summary.csv"
      )
    )
  )
)

step_ids <- vapply(
  workflow,
  function(step) {
    step$id
  },
  character(1)
)

if (
  nzchar(
    only_step
  ) &&
    !only_step %in%
      step_ids
) {
  stop(
    "Unknown workflow step '",
    only_step,
    "'. Available steps:\n  ",
    paste(
      step_ids,
      collapse = "\n  "
    )
  )
}

if (nzchar(
  only_step
)) {
  workflow <- workflow[
    step_ids ==
      only_step
  ]
}

# ------------------------------------------------------------
# 5. Preflight checks
# ------------------------------------------------------------

# The reproducible LME workflow starts from the 13 prepared,
# gzip-compressed member files. The original NetCDF files are
# outside the distributed workflow.
required_lme_member_files <- file.path(
  repository_root,
  "Data",
  "LME data",
  paste0(
    "a",
    1:13,
    ".csv.gz"
  )
)

missing_lme_member_files <- required_lme_member_files[
  !file.exists(
    required_lme_member_files
  )
]

if (length(
  missing_lme_member_files
) > 0L) {
  stop(
    "The reproducible workflow begins with the prepared LME ",
    "member files a1.csv.gz through a13.csv.gz. The following ",
    "files are missing:\n  ",
    paste(
      missing_lme_member_files,
      collapse = "\n  "
    ),
    "\nSee Data/Get_data/README.md for provenance and the ",
    "optional NetCDF-to-CSV preparation step."
  )
}

required_core_input_files <- c(
  file.path(
    repository_root,
    "Data",
    "temperature index value.v1.xlsx"
  ),
  file.path(
    repository_root,
    "Data",
    "GHCNv4.xlsx"
  )
)

missing_core_input_files <- required_core_input_files[
  !file.exists(
    required_core_input_files
  )
]

if (length(
  missing_core_input_files
) > 0L) {
  stop(
    "Required core input file(s) are missing:\n  ",
    paste(
      missing_core_input_files,
      collapse = "\n  "
    )
  )
}

missing_scripts <- vapply(
  workflow,
  function(step) {
    !file.exists(
      file.path(
        repository_root,
        step$script
      )
    )
  },
  logical(1)
)

if (any(
  missing_scripts
)) {
  stop(
    "The following workflow scripts do not exist:\n  ",
    paste(
      vapply(
        workflow[
          missing_scripts
        ],
        function(step) {
          step$script
        },
        character(1)
      ),
      collapse = "\n  "
    ),
    "\nUpdate the workflow paths in run_all.R."
  )
}

rscript_executable <- Sys.which(
  "Rscript"
)

if (!nzchar(
  rscript_executable
)) {
  stop(
    "Rscript was not found on PATH."
  )
}

# Report suspicious local paths without preventing the run.
code_files <- list.files(
  file.path(
    repository_root,
    "Code"
  ),
  pattern = "\\.[Rr]$",
  recursive = TRUE,
  full.names = TRUE
)

active_code_files <- code_files[
  !grepl(
    "/(Archive|archive|Old|old)/",
    gsub(
      "\\\\",
      "/",
      code_files
    )
  )
]

suspicious_patterns <- c(
  "~/Downloads",
  "setwd\\(",
  "/Users/",
  "[A-Za-z]:\\\\"
)

suspicious_hits <- list()

for (code_file in active_code_files) {

  code_lines <- readLines(
    code_file,
    warn = FALSE
  )

  matched_lines <- which(
    vapply(
      code_lines,
      function(code_line) {
        any(
          vapply(
            suspicious_patterns,
            function(pattern) {
              grepl(
                pattern,
                code_line
              )
            },
            logical(1)
          )
        )
      },
      logical(1)
    )
  )

  if (length(
    matched_lines
  ) > 0L) {
    suspicious_hits[[code_file]] <- matched_lines
  }
}

if (length(
  suspicious_hits
) > 0L) {
  warning(
    "Potential hard-coded local paths remain in active Code/ files:\n",
    paste(
      vapply(
        names(
          suspicious_hits
        ),
        function(code_file) {
          paste0(
            normalizePath(
              code_file,
              winslash = "/",
              mustWork = FALSE
            ),
            ":",
            paste(
              suspicious_hits[[code_file]],
              collapse = ","
            )
          )
        },
        character(1)
      ),
      collapse = "\n"
    )
  )
}

# ------------------------------------------------------------
# 6. Execute each step in a separate clean R process
# ------------------------------------------------------------

run_manifest_file <- file.path(
  log_directory,
  "run_manifest.csv"
)

run_records <- vector(
  "list",
  length(
    workflow
  )
)

write_run_manifest <- function(
    records) {

  completed_records <- records[
    !vapply(
      records,
      is.null,
      logical(1)
    )
  ]

  if (length(
    completed_records
  ) == 0L) {
    return(
      invisible(
        NULL
      )
    )
  }

  manifest <- do.call(
    rbind,
    completed_records
  )

  utils::write.csv(
    manifest,
    run_manifest_file,
    row.names = FALSE,
    na = ""
  )

  invisible(
    manifest
  )
}

for (step_index in seq_along(
  workflow
)) {

  step <- workflow[[step_index]]

  script_file <- normalizePath(
    file.path(
      repository_root,
      step$script
    ),
    winslash = "/",
    mustWork = TRUE
  )

  log_file <- file.path(
    log_directory,
    paste0(
      sprintf(
        "%02d",
        step_index
      ),
      "_",
      step$id,
      ".log"
    )
  )

  message(
    "\n[",
    step_index,
    "/",
    length(
      workflow
    ),
    "] Running ",
    step$id,
    "\n    Script: ",
    step$script,
    "\n    Log:    ",
    log_file
  )

  start_time <- Sys.time()

  status <- system2(
    command = rscript_executable,
    args = c(
      "--vanilla",
      shQuote(
        script_file
      )
    ),
    stdout = log_file,
    stderr = log_file,
    wait = TRUE
  )

  end_time <- Sys.time()

  elapsed_seconds <- as.numeric(
    difftime(
      end_time,
      start_time,
      units = "secs"
    )
  )

  expected_files <- file.path(
    repository_root,
    step$expected
  )

  missing_outputs <- expected_files[
    !file.exists(
      expected_files
    )
  ]

  step_success <- identical(
    status,
    0L
  ) &&
    length(
      missing_outputs
    ) ==
      0L

  run_records[[step_index]] <- data.frame(
    step_number = step_index,
    step_id = step$id,
    script = step$script,
    start_time = format(
      start_time,
      "%Y-%m-%d %H:%M:%S %Z"
    ),
    end_time = format(
      end_time,
      "%Y-%m-%d %H:%M:%S %Z"
    ),
    elapsed_seconds =
      elapsed_seconds,
    exit_status = status,
    success = step_success,
    missing_outputs = paste(
      sub(
        paste0(
          "^",
          gsub(
            "([][{}()+*^$|\\\\?.])",
            "\\\\\\1",
            repository_root
          ),
          "/?"
        ),
        "",
        missing_outputs
      ),
      collapse = "; "
    ),
    log_file = file.path(
      "Output",
      "Logs",
      basename(
        log_file
      )
    ),
    stringsAsFactors = FALSE
  )

  write_run_manifest(
    run_records
  )

  if (!step_success) {

    log_lines <- if (file.exists(
      log_file
    )) {
      readLines(
        log_file,
        warn = FALSE
      )
    } else {
      character(
        0
      )
    }

    log_tail <- tail(
      log_lines,
      50L
    )

    message(
      "\nLast lines of the failed log:\n",
      paste(
        log_tail,
        collapse = "\n"
      )
    )

    if (status != 0L) {
      stop(
        "Workflow step '",
        step$id,
        "' failed with exit status ",
        status,
        ". See ",
        log_file,
        "."
      )
    }

    stop(
      "Workflow step '",
      step$id,
      "' completed but did not create:\n  ",
      paste(
        missing_outputs,
        collapse = "\n  "
      )
    )
  }

  message(
    "Completed ",
    step$id,
    " in ",
    round(
      elapsed_seconds /
        60,
      2
    ),
    " minutes."
  )
}

final_manifest <- write_run_manifest(
  run_records
)

message(
  "\nReproducibility workflow completed successfully."
)

message(
  "Run manifest: ",
  run_manifest_file
)

message(
  "Logs: ",
  log_directory
)
