# ==============================================================
# project_paths.R -- robust project-root and path handling
# ============================================================== 

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L) return(y)
  if (length(x) == 1L && is.na(x)) return(y)
  x
}

normalize_path2 <- function(path, mustWork = FALSE) {
  normalizePath(path, winslash = "/", mustWork = mustWork)
}

parse_named_args <- function(args = commandArgs(trailingOnly = TRUE)) {
  out <- list()
  if (length(args) == 0L) return(out)
  named <- grep("^--[A-Za-z0-9_.-]+=", args, value = TRUE)
  if (length(named) == 0L) return(out)
  for (x in named) {
    key <- sub("^--([^=]+)=.*$", "\\1", x)
    val <- sub("^--[^=]+=", "", x)
    out[[key]] <- val
  }
  out
}

get_script_path <- function() {
  args_full <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args_full, value = TRUE)
  if (length(file_arg) >= 1L) {
    return(normalize_path2(sub("^--file=", "", file_arg[1L]), mustWork = FALSE))
  }

  ofile <- tryCatch(sys.frames()[[1]]$ofile, error = function(e) NULL)
  if (!is.null(ofile) && nzchar(ofile)) {
    return(normalize_path2(ofile, mustWork = FALSE))
  }

  if (requireNamespace("rstudioapi", quietly = TRUE)) {
    p1 <- tryCatch(rstudioapi::getSourceEditorContext()$path, error = function(e) "")
    if (nzchar(p1)) return(normalize_path2(p1, mustWork = FALSE))
    p2 <- tryCatch(rstudioapi::getActiveDocumentContext()$path, error = function(e) "")
    if (nzchar(p2)) return(normalize_path2(p2, mustWork = FALSE))
  }

  NULL
}

find_project_root <- function(start = NULL,
                              markers = c("run_null_streamingcurve_simulation.R", file.path("R", "project_paths.R"))) {
  here <- start %||% get_script_path() %||% getwd()
  here <- normalize_path2(here, mustWork = FALSE)
  if (!dir.exists(here)) here <- dirname(here)

  cur <- here
  repeat {
    ok <- all(file.exists(file.path(cur, markers)))
    if (ok) return(cur)
    parent <- dirname(cur)
    if (identical(parent, cur)) break
    cur <- parent
  }

  stop(
    "Could not locate the simulation project root. Run the script from inside the project folder, ",
    "or pass --root=/path/to/functional-representation-monitoring/simulations/contaminated_training."
  )
}

resolve_project_root <- function(default_start = NULL) {
  args <- parse_named_args()
  env_root <- Sys.getenv("STREAMING_CURVE_SIM_ROOT", unset = "")

  if (!is.null(args$root) && nzchar(args$root)) {
    return(normalize_path2(args$root, mustWork = FALSE))
  }
  if (nzchar(env_root)) {
    return(normalize_path2(env_root, mustWork = FALSE))
  }

  find_project_root(start = default_start)
}

project_path <- function(root, ...) file.path(root, ...)

ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

backup_file_if_exists <- function(path) {
  if (!file.exists(path)) return(invisible(NULL))
  stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  bak <- paste0(path, ".bak_", stamp)
  file.copy(path, bak, overwrite = TRUE)
  invisible(bak)
}
