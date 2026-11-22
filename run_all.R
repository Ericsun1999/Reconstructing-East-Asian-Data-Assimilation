here::i_am("run_all.R")

# Legacy single-step compatibility. The default entry point now runs
# both documentary input versions, while --only and --clean remain
# available to reproduce a specific historical workflow stage.
legacy_arguments <- commandArgs(trailingOnly = TRUE)
if (any(grepl("^--only=", legacy_arguments)) ||
    "--clean" %in% legacy_arguments ||
    "--skip-shared" %in% legacy_arguments) {
  source(here::here("Code", "Workflow", "run_single_version.R"),
         local = FALSE)
  quit(save = "no", status = 0L)
}

# Independent, fully versioned reproductions:
#   inferred-zero = primary reconstruction (Original + documented zeros)
#   original      = sensitivity reconstruction (temperature records only)
#
# All figure, supplementary, calibration, kriging, posterior, validation
# and log outputs remain in separate workspaces under:
#   Output/Reconstructions/inferred-zero/Output/
#   Output/Reconstructions/original/Output/
#
# The independent LME products and penalized prior are prepared ONCE:
#   Output/Reconstructions/shared/Output/Intermediate/{LME,Prior}/
# Both versioned workspaces link those directories read-only in intent.
#
# Usage (from the real repository root):
#   Rscript --vanilla run_all.R --all
#   Rscript --vanilla run_all.R --prepare-shared
#   Rscript --vanilla run_all.R --variant=inferred-zero
#   Rscript --vanilla run_all.R --variant=original
#   Rscript --vanilla run_all.R --compare-only
#   Rscript --vanilla run_all.R --all --u=0
#
# By default u comes from the original blocked-CV selected_u.csv.
# This controller does not change the statistical methods.

arguments <- commandArgs(trailingOnly = TRUE)
choices <- c("--all", "--prepare-shared", "--compare-only",
             "--rebuild-shared")
variant_arg <- grep("^--variant=", arguments, value = TRUE)
u_arg <- grep("^--u=", arguments, value = TRUE)
unknown <- arguments[
  !arguments %in% choices &
  !grepl("^--variant=", arguments) &
  !grepl("^--u=", arguments)
]
if (length(unknown) ||
    length(variant_arg) > 1L || length(u_arg) > 1L) {
  stop("Unknown or repeated arguments: ",
       paste(c(unknown, variant_arg[-1L], u_arg[-1L]),
             collapse = ", "))
}
actions <- c(
  all = "--all" %in% arguments,
  prepare_shared = "--prepare-shared" %in% arguments,
  compare_only = "--compare-only" %in% arguments,
  variant = length(variant_arg) == 1L
)
if (sum(actions) > 1L) {
  stop("Choose one action: --all, --prepare-shared, --variant=..., ",
       "or --compare-only.")
}
if (!any(actions)) actions["all"] <- TRUE
requested_variant <- if (actions["variant"]) {
  sub("^--variant=", "", variant_arg)
} else {
  ""
}
if (actions["variant"] &&
    !requested_variant %in% c("inferred-zero", "original")) {
  stop("--variant must be inferred-zero or original.")
}
if ("--rebuild-shared" %in% arguments && actions["compare_only"]) {
  stop("--rebuild-shared is not meaningful with --compare-only.")
}

repository_root <- here::here()
reconstruction_root <- file.path(
  repository_root, "Output", "Reconstructions"
)
rscript <- Sys.which("Rscript")
if (!nzchar(rscript)) stop("Rscript was not found on PATH.")
reproduction_mode <- Sys.getenv("REPRODUCTION_MODE", unset = "full")
if (!reproduction_mode %in% c("full", "smoke")) {
  stop("REPRODUCTION_MODE must be full or smoke.")
}
if (nzchar(Sys.getenv("WORLDCLIM_OFFSETS_FILE", unset = ""))) {
  stop("The versioned common run uses the repository's prepared ",
       "WorldClim offsets. Unset WORLDCLIM_OFFSETS_FILE so the ",
       "shared-input fingerprint matches the actual data.")
}
# Calendar holdout settings are for cross-validation, not a full
# reconstruction. Never allow inherited terminal exports to silently
# restrict the historical source years or reroute calibration output.
Sys.unsetenv(c(
  "REACHES_CV_HOLDOUT_START", "REACHES_CV_HOLDOUT_END",
  "REACHES_CV_CALIBRATION_OUTPUT_DIR", "REACHES_CV_ONLY_FOLD",
  "REACHES_CV_REUSE_CITY_KRIGING"
))
if (actions["all"] || actions["variant"] ||
    actions["prepare_shared"]) {
  dir.create(reconstruction_root, recursive = TRUE, showWarnings = FALSE)
}

# Resolve the selected Original-input variance estimate once. Both documentary inputs
# deliberately use the SAME observation error setting.
fixed_u <- NA_real_
variance_source <- ""
if (actions["all"] || actions["variant"]) {
  u_option <- if (length(u_arg)) {
    sub("^--u=", "", u_arg[1L])
  } else {
    "from-cv"
  }
  if (identical(u_option, "from-cv")) {
    cv_file <- file.path(
      repository_root, "Output", "Intermediate",
      "ObservationVariance", "original", "selected_u.csv"
    )
    if (!file.exists(cv_file)) {
      stop("Missing original blocked-CV variance: ", cv_file,
           "\nUse --u=<nonnegative-number> only if you have ",
           "independently verified the intended original-model u.")
    }
    chosen <- read.csv(cv_file, check.names = FALSE)
    if (nrow(chosen) != 1L || !"selected_u" %in% names(chosen)) {
      stop("The original CV estimate is incomplete.")
    }
    if ("temperature_input_mode" %in% names(chosen) &&
        !identical(as.character(chosen$temperature_input_mode[1L]),
                   "original")) {
      stop("The selected variance is not from Original observations.")
    }
    if ("measurement_projection_method" %in% names(chosen) &&
        !is.na(chosen$measurement_projection_method[1L]) &&
        !identical(
          as.character(chosen$measurement_projection_method[1L]),
          "empirical_reference_joint_mc_v1"
        )) {
      stop("The saved u was fitted under another measurement ",
           "projection; do not use it with the original model.")
    }
    fixed_u <- as.numeric(chosen$selected_u[1L])
    variance_source <- normalizePath(cv_file, winslash = "/")
  } else {
    fixed_u <- suppressWarnings(as.numeric(u_option))
    variance_source <- paste0("explicit:", u_option)
  }
  if (length(fixed_u) != 1L ||
      !is.finite(fixed_u) || fixed_u < 0) {
    stop("The common additional observation variance must be ",
         "a finite nonnegative number.")
  }
}

is_link <- function(path) {
  link <- Sys.readlink(path)
  length(link) == 1L && !is.na(link) && nzchar(link)
}
make_link <- function(target, path) {
  if (!file.exists(target) && !dir.exists(target)) {
    stop("Link source does not exist: ", target)
  }
  if (file.exists(path) || dir.exists(path) || is_link(path)) {
    if (is_link(path) &&
        identical(normalizePath(path, winslash = "/"),
                  normalizePath(target, winslash = "/"))) {
      return(invisible(NULL))
    }
    stop("Refusing to overwrite nonmatching workspace entry: ", path)
  }
  if (!file.symlink(target, path)) {
    stop("Could not link workspace source ", target, " to ", path)
  }
  invisible(NULL)
}

# A local Data/Valid folder prevents any run from mutating the
# Git-tracked validation files in the actual repository. All other
# Data entries are input-only links to the exact same source files.
prepare_workspace <- function(name) {
  workspace <- file.path(reconstruction_root, name)
  dir.create(workspace, recursive = TRUE, showWarnings = FALSE)
  file.create(file.path(workspace, ".here"))
  make_link(file.path(repository_root, "Code"),
            file.path(workspace, "Code"))
  make_link(file.path(repository_root, "renv"),
            file.path(workspace, "renv"))
  make_link(file.path(repository_root, "renv.lock"),
            file.path(workspace, "renv.lock"))
  profile <- file.path(repository_root, ".Rprofile")
  if (file.exists(profile)) {
    # Copy rather than link a writable startup file into the source tree.
    profile_copy <- file.path(workspace, ".Rprofile")
    if (!file.copy(profile, profile_copy, overwrite = TRUE)) {
      stop("Unable to copy renv startup profile to ", workspace)
    }
  }
  run_from <- file.path(repository_root, "Code", "Workflow",
                        "run_single_version.R")
  run_to <- file.path(workspace, "run_all.R")
  if (!file.copy(run_from, run_to, overwrite = TRUE)) {
    stop("Could not copy run_all.R into ", workspace)
  }

  data_dir <- file.path(workspace, "Data")
  dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
  original_data <- file.path(repository_root, "Data")
  entries <- list.files(original_data, all.files = TRUE,
                        no.. = TRUE, full.names = TRUE)
  for (entry in entries) {
    if (basename(entry) != "Valid") {
      make_link(entry, file.path(data_dir, basename(entry)))
    }
  }
  dir.create(file.path(data_dir, "Valid"), showWarnings = FALSE)

  out_dir <- file.path(workspace, "Output")
  dir.create(file.path(out_dir, "Intermediate"),
             recursive = TRUE, showWarnings = FALSE)
  if (name != "shared") {
    shared_intermediate <- file.path(
      reconstruction_root, "shared", "Output", "Intermediate"
    )
    for (item in c("LME", "Prior")) {
      make_link(
        file.path(shared_intermediate, item),
        file.path(out_dir, "Intermediate", item)
      )
    }
  }
  workspace
}

common_outputs <- c(
  file.path("LME", "lme_annual_1368_1911.rds"),
  file.path("LME", "lme_city3_annual_1368_1911.csv"),
  file.path("LME", "lme_ensemble_mean_1368_1911.csv"),
  file.path("LME", "WorldClim_Corrected",
            "lme_city3_annual_1368_1911.csv"),
  file.path("LME", "WorldClim_Corrected",
            "lme_land_grid_corrected_annual_1368_1911.rds"),
  file.path("Prior", c(
    "mtB.csv", "muB.csv", "rtB.csv",
    "mtS.csv", "muS.csv", "rtS.csv",
    "mtH.csv", "muH.csv", "rtH.csv"
  ))
)
common_root <- file.path(
  reconstruction_root, "shared", "Output", "Intermediate"
)
common_manifest <- file.path(
  reconstruction_root, "shared", "shared_inputs_manifest.rds"
)

# The inputs and source code that determine the independent LME/prior
# outputs are fingerprinted. Changing them requires an explicit
# --rebuild-shared, preventing silent reuse of stale prior parameters.
dependency_files <- c(
  file.path(repository_root, "Code", "Prior.R"),
  file.path(repository_root, "Code", "DataPreparation",
            "prepare_lme_annual.R"),
  file.path(repository_root, "Code", "DataPreparation",
            "prepare_worldclim_bias_correction.R"),
  file.path(repository_root, "Data", "WorldClim",
            "worldclim_offsets_all_269.csv"),
  file.path(repository_root, "Data", "LME data",
            paste0("a", 1:13, ".csv.gz")),
  file.path(repository_root, "renv.lock")
)
missing_dependencies <- dependency_files[!file.exists(dependency_files)]
if ((actions["all"] || actions["variant"] ||
     actions["prepare_shared"]) && length(missing_dependencies)) {
  stop("Missing input(s) for shared LME/prior generation:\n",
       paste(missing_dependencies, collapse = "\n"))
}
common_fingerprint <- function() {
  list(
    reproduction_mode = reproduction_mode,
    inputs = unname(tools::md5sum(dependency_files)),
    paths = dependency_files
  )
}

run_in_workspace <- function(workspace, args) {
  previous <- getwd()
  on.exit(setwd(previous), add = TRUE)
  setwd(workspace)
  status <- system2(
    rscript,
    args = c("--vanilla", shQuote(file.path(workspace, "run_all.R")),
             args),
    stdout = "", stderr = "", wait = TRUE
  )
  if (!identical(as.integer(status), 0L)) {
    stop("Workflow failed in ", workspace, " (status ", status, "). ",
         "Read the step logs under Output/Logs in that workspace.")
  }
  invisible(NULL)
}
run_script_in_workspace <- function(workspace, script, args) {
  previous <- getwd()
  on.exit(setwd(previous), add = TRUE)
  setwd(workspace)
  status <- system2(
    rscript,
    args = c("--vanilla", shQuote(file.path(workspace, script)),
             args),
    stdout = "", stderr = "", wait = TRUE
  )
  if (!identical(as.integer(status), 0L)) {
    stop("Script failed: ", script, " in ", workspace)
  }
  invisible(NULL)
}

prepare_shared <- function(rebuild = FALSE) {
  present <- all(file.exists(file.path(common_root, common_outputs)))
  previous <- if (file.exists(common_manifest)) {
    readRDS(common_manifest)
  } else {
    NULL
  }
  current <- common_fingerprint()
  matching <- present &&
    is.list(previous) &&
    identical(previous$fingerprint, current) &&
    identical(
      previous$output_md5,
      unname(tools::md5sum(
        file.path(common_root, common_outputs)
      ))
    )
  if (matching && !rebuild) {
    message("Reusing the verified independent LME and prior (once).")
    return(invisible(NULL))
  }
  if (present && !rebuild) {
    stop("Existing shared LME/prior outputs do not match the ",
         "current input and code fingerprint. Use --rebuild-shared ",
         "deliberately, then rerun BOTH documentary-input versions.")
  }

  workspace <- prepare_workspace("shared")
  if (rebuild) {
    for (part in c("LME", "Prior")) {
      old <- file.path(common_root, part)
      if (is_link(old)) {
        stop("Shared output unexpectedly contains a symlink: ", old)
      }
      if (dir.exists(old)) unlink(old, recursive = TRUE)
    }
    if (file.exists(common_manifest)) file.remove(common_manifest)
  }
  Sys.setenv(
    REPRODUCTION_MODE = reproduction_mode,
    REACHES_VERSIONED_RUN = "shared",
    REACHES_TEMPERATURE_INPUT = "original",
    REACHES_EXTRA_VARIANCE_FIXED = "",
    REACHES_EXTRA_VARIANCE_FILE = ""
  )
  for (step in c("prepare_lme", "correct_lme_city", "prior")) {
    message("\nShared preparation: ", step)
    run_in_workspace(workspace, paste0("--only=", step))
  }
  missing <- file.path(common_root, common_outputs)[
    !file.exists(file.path(common_root, common_outputs))
  ]
  if (length(missing)) {
    stop("Shared LME/prior preparation was incomplete:\n",
         paste(missing, collapse = "\n"))
  }
  saveRDS(
    list(fingerprint = current,
         created = Sys.time(),
         output_md5 = unname(tools::md5sum(
           file.path(common_root, common_outputs)
         ))),
    common_manifest
  )
  message("Shared prior has been generated ONCE.")
  invisible(NULL)
}

check_shared <- function() {
  if (!file.exists(common_manifest) ||
      !all(file.exists(file.path(common_root, common_outputs)))) {
    stop("Shared LME and prior missing; use --prepare-shared or --all.")
  }
  previous <- readRDS(common_manifest)
  if (!identical(previous$fingerprint, common_fingerprint()) ||
      !identical(previous$output_md5,
                 unname(tools::md5sum(
                   file.path(common_root, common_outputs)
                 )))) {
    stop("Shared prior/LME inputs changed. Use --rebuild-shared ",
         "and rerun both versions.")
  }
}

run_variant <- function(name) {
  check_shared()
  temperature_mode <- if (name == "inferred-zero") "augmented" else "original"
  workspace <- prepare_workspace(name)
  for (stage in c("city", "grid")) {
    existing <- file.path(
      workspace, "Output", "Intermediate",
      "AugmentationSensitivity", temperature_mode, stage,
      "snapshot_manifest.csv"
    )
    if (file.exists(existing)) {
      stop("A previous completed snapshot exists at: ", existing,
           ". Archive or explicitly remove that previous version ",
           "before starting another full run.")
    }
  }
  Sys.setenv(
    REPRODUCTION_MODE = reproduction_mode,
    REACHES_VERSIONED_RUN = name,
    REACHES_TEMPERATURE_INPUT = temperature_mode,
    REACHES_WRITE_VALIDATION_COPY = "0",
    REACHES_MEASUREMENT_PROJECTION_METHOD =
      "empirical_reference_joint_mc_v1",
    REACHES_EXTRA_VARIANCE_FIXED = sprintf("%.17g", fixed_u),
    REACHES_EXTRA_VARIANCE_FILE = ""
  )
  message("\nStarting ", name, " full reconstruction ",
          "(shared LME/prior reused; original measurement equation).")
  run_in_workspace(workspace, "--skip-shared")

  manifest_file <- file.path(workspace, "version_manifest.csv")
  result <- data.frame(
    variant = name,
    documentary_input = temperature_mode,
    role = if (name == "inferred-zero") {
      "primary_inferred_zero"
    } else {
      "original_record_sensitivity"
    },
    shared_prior_dir = file.path(common_root, "Prior"),
    shared_lme_dir = file.path(common_root, "LME"),
    shared_manifest_md5 = unname(tools::md5sum(common_manifest)),
    fixed_u = fixed_u,
    u_source = variance_source,
    projection_method = "empirical_reference_joint_mc_v1",
    reproduction_mode = reproduction_mode,
    completed_at = as.character(Sys.time())
  )
  write.csv(result, manifest_file, row.names = FALSE)
  message("Reconstruction complete: ", workspace, "/Output")

  # Full-grid snapshot also contains all three city reconstructions.
  for (stage in c("city", "grid")) {
    saved <- file.path(
      workspace, "Output", "Intermediate",
      "AugmentationSensitivity", temperature_mode, stage,
      "snapshot_manifest.csv"
    )
    if (file.exists(saved)) {
      # Re-running should never silently overwrite an older snapshot.
      stop("Existing ", name, " snapshot blocks replacement: ",
           saved, ". Archive the old snapshot first.")
    }
    run_script_in_workspace(
      workspace,
      "Code/Analysis/snapshot_augmentation_reconstruction.R",
      c(temperature_mode, stage)
    )
  }
  invisible(NULL)
}

compare_versions <- function() {
  manifests <- lapply(c("inferred-zero", "original"), function(name) {
    p <- file.path(reconstruction_root, name, "version_manifest.csv")
    if (!file.exists(p)) stop("Missing completed version manifest: ", p)
    m <- read.csv(p, stringsAsFactors = FALSE)
    if (nrow(m) != 1L ||
        !identical(as.character(m$variant[1L]), name)) {
      stop("Invalid version metadata for ", name)
    }
    m
  })
  a <- manifests[[1L]]
  b <- manifests[[2L]]
  same <- c("shared_manifest_md5", "fixed_u", "projection_method",
            "reproduction_mode")
  if (!all(same %in% names(a)) || !all(same %in% names(b)) ||
      !identical(as.list(a[1L, same]), as.list(b[1L, same])) ||
      !identical(as.character(a$role[1L]), "primary_inferred_zero") ||
      !identical(as.character(b$role[1L]),
                 "original_record_sensitivity")) {
    stop("Version manifests are not compatible; the two runs must ",
         "use the same LME/prior, fixed u, projection and mode.")
  }
  Sys.setenv(REACHES_COMPARISON_VERSION_ROOT = reconstruction_root)
  for (stage in c("city", "grid")) {
    message("\nComparing ", stage, " reconstructions ...")
    previous <- getwd()
    status <- system2(
      rscript,
      args = c(
        "--vanilla",
        shQuote(file.path(
          repository_root,
          "Code", "Analysis", "compare_augmentation_reconstructions.R"
        )),
        stage
      ),
      stdout = "", stderr = "", wait = TRUE
    )
    if (!identical(as.integer(status), 0L)) {
      stop("Comparison failed for ", stage,
           "; inspect the versioned snapshot manifests.")
    }
  }
  message("Main figures: ", file.path(
    reconstruction_root, "inferred-zero", "Output"
  ))
  message("Original sensitivity figures: ", file.path(
    reconstruction_root, "original", "Output"
  ))
  message("Comparisons: ", file.path(
    reconstruction_root, "Comparison"
  ))
}

if (actions["compare_only"]) {
  compare_versions()
} else {
  prepare_shared(rebuild = "--rebuild-shared" %in% arguments)
  if (actions["all"]) {
    run_variant("inferred-zero")
    run_variant("original")
    compare_versions()
  } else if (actions["variant"]) {
    run_variant(requested_variant)
  }
}
