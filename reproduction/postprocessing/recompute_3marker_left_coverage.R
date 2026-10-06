# Recompute the 3-marker one-sided upper-endpoint coverage from saved interval objects.
#
# The historical April 2026 summary helper counted an infinite upper endpoint as
# noncoverage. For an interval [landmark, Inf], every finite event time is covered.
# This script corrects only that summary metric; no PRC model is refitted.
#
# Example:
# Rscript reproduction/postprocessing/recompute_3marker_left_coverage.R \
#   --root /path/to/landmark_longitudinal_two_scenarios_N300_1000_3000 \
#   --output corrected_3marker_left_coverage.csv

get_arg <- function(flag, default = NULL) {
  args <- commandArgs(trailingOnly = TRUE)
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) return(default)
  args[i + 1L]
}

root <- get_arg("--root", Sys.getenv("SIM3_ROOT", ""))
out_file <- get_arg("--output", "corrected_3marker_left_coverage.csv")

if (!nzchar(root) || !dir.exists(root)) {
  stop("Provide the original 3-marker results root with --root or SIM3_ROOT.")
}

collect_files <- function(path) {
  if (!dir.exists(path)) return(character(0))
  list.files(path, pattern = "\\.RData$", recursive = TRUE, full.names = TRUE)
}

files <- c(
  collect_files(file.path(root, "Naive", "Results")),
  collect_files(file.path(root, "Qin", "Results"))
)

if (!length(files)) stop("No saved result .RData files found.")

rows <- vector("list", length(files))
k <- 1L

for (fp in files) {
  e <- new.env(parent = emptyenv())
  ok <- tryCatch({ load(fp, envir = e); TRUE }, error = function(err) FALSE)
  if (!ok || !exists("result_row", envir = e, inherits = FALSE) ||
      !exists("interval_df", envir = e, inherits = FALSE)) next

  rr <- as.data.frame(e$result_row)
  int <- as.data.frame(e$interval_df)
  needed_rr <- c("method", "scenario", "N", "cens", "lmk", "rep")
  if (!all(needed_rr %in% names(rr))) next
  if (!all(c("true_time", "upper_left") %in% names(int))) next

  # Keep the final design cells only.
  if (!(rr$scenario[1] %in% c("weakNPH", "strongNPH"))) next
  if (!(as.integer(rr$N[1]) %in% c(300L, 1000L, 1500L))) next
  if (!(as.integer(rr$cens[1]) %in% c(20L, 35L, 50L))) next
  if (!(as.numeric(rr$lmk[1]) %in% c(2, 4))) next

  method <- tolower(as.character(rr$method[1]))
  if (method %in% c("qin", "conformal") && "B" %in% names(rr)) {
    if (as.integer(rr$B[1]) != 500L) next
  }

  tt <- as.numeric(int$true_time)
  up <- as.numeric(int$upper_left)
  corrected <- mean(tt <= up, na.rm = TRUE)
  historical <- if ("cov_left" %in% names(rr)) as.numeric(rr$cov_left[1]) else NA_real_

  rows[[k]] <- data.frame(
    file = fp,
    method = method,
    scenario = as.character(rr$scenario[1]),
    N = as.integer(rr$N[1]),
    cens = as.integer(rr$cens[1]),
    landmark = as.numeric(rr$lmk[1]),
    rep = as.integer(rr$rep[1]),
    B = if ("B" %in% names(rr)) as.integer(rr$B[1]) else NA_integer_,
    historical_left_cov = historical,
    corrected_left_cov = corrected,
    stringsAsFactors = FALSE
  )
  k <- k + 1L
}

rows <- rows[seq_len(k - 1L)]
if (!length(rows)) stop("No compatible final result objects found.")
d <- do.call(rbind, rows)

# Resolve duplicate files for the same design replicate by keeping the newest file.
d$mtime <- file.info(d$file)$mtime
key <- interaction(d$method, d$scenario, d$N, d$cens, d$landmark, d$rep, drop = TRUE)
keep <- unlist(lapply(split(seq_len(nrow(d)), key), function(ii) ii[which.max(d$mtime[ii])]), use.names = FALSE)
d <- d[sort(keep), , drop = FALSE]

agg <- aggregate(
  cbind(historical_left_cov, corrected_left_cov) ~ method + scenario + N + cens + landmark,
  data = d,
  FUN = mean,
  na.rm = TRUE
)
counts <- aggregate(rep ~ method + scenario + N + cens + landmark, data = d, FUN = length)
names(counts)[names(counts) == "rep"] <- "n_reps"
agg <- merge(agg, counts, by = c("method", "scenario", "N", "cens", "landmark"), all.x = TRUE)
agg$change <- agg$corrected_left_cov - agg$historical_left_cov
agg <- agg[order(agg$scenario, agg$method, agg$cens, agg$landmark, agg$N), ]

write.csv(agg, out_file, row.names = FALSE)
print(agg)
cat("\nSaved corrected grouped coverage to:", normalizePath(out_file, mustWork = FALSE), "\n")
cat("No model fitting was repeated. Only the stored interval endpoints were re-summarised.\n")
