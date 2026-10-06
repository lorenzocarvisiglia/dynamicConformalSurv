# Re-aggregate the existing ADNI interval CSV files without refitting models.
#
# Corrections relative to an old development aggregator:
# 1) baseline-Cox files are restricted to method == "conformal" instead of
#    pooling naive and conformal baseline rows under one label;
# 2) both finite-only and eta-truncated upper lengths are reported explicitly.
#
# Example:
# Rscript reproduction/postprocessing/reaggregate_adni_existing.R \
#   --root /path/to/adni_qin_cv --output adni_cv_table_corrected.csv

if (!requireNamespace("data.table", quietly = TRUE)) {
  stop("Package 'data.table' is required.")
}

get_arg <- function(flag, default = NULL) {
  args <- commandArgs(trailingOnly = TRUE)
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) return(default)
  args[i + 1L]
}

root <- get_arg("--root", Sys.getenv("ADNIROOT", ""))
out_file <- get_arg("--output", "adni_cv_table_corrected.csv")
if (!nzchar(root) || !dir.exists(root)) stop("Provide the ADNI analysis root with --root or ADNIROOT.")

files <- list.files(
  file.path(root, "output"),
  pattern = "intervals.*fold.*lmk.*B0500\\.csv$",
  recursive = TRUE,
  full.names = TRUE
)
if (!length(files)) stop("No B=500 ADNI interval CSV files found.")

file_tables <- setNames(lapply(files, data.table::fread), files)
x <- data.table::rbindlist(file_tables, fill = TRUE, idcol = "source_file")
if ("alpha" %in% names(x)) x <- x[abs(alpha - 0.10) < 1e-12]

needed <- c("fold", "method", "landmark", "obs_time", "event", "lower_right", "upper_left", "eta")
missing <- setdiff(needed, names(x))
if (length(missing)) stop("Missing columns: ", paste(missing, collapse = ", "))

x[, source_is_baseline := grepl("baseline", source_file, ignore.case = TRUE)]

# The baseline worker saved both naive and conformal rows in each baseline file.
# The thesis comparison requires conformal baseline Cox only.
x <- x[!source_is_baseline | tolower(method) == "conformal"]

x[, method_label := NA_character_]
x[source_is_baseline, method_label := "Conformal baseline Cox"]
x[!source_is_baseline & tolower(method) == "conformal", method_label := "Conformal PRC"]
x[!source_is_baseline & tolower(method) == "naive", method_label := "Naive PRC"]
x <- x[method_label %in% c("Naive PRC", "Conformal PRC", "Conformal baseline Cox")]

x[, lmk := if ("lmk_value" %in% names(x)) lmk_value else landmark]
x[, upper_inf := is.na(upper_left) | is.infinite(upper_left)]
x[, upper_length_finite := upper_left - lmk]
x[!is.finite(upper_length_finite) | upper_length_finite < 0, upper_length_finite := NA_real_]
x[, upper_left_trunc := data.table::fifelse(upper_inf, eta, pmin(upper_left, eta))]
x[, upper_length_truncated := upper_left_trunc - lmk]
x[!is.finite(upper_length_truncated) | upper_length_truncated < 0, upper_length_truncated := NA_real_]

ev <- x[event == 1]
ev[, upper_cov_ok := data.table::fifelse(upper_inf, TRUE, obs_time <= upper_left)]
ev[, lower_trunc_cov_ok := data.table::fifelse(
  is.na(lower_right) | is.infinite(lower_right) | !is.finite(eta),
  FALSE,
  lower_right <= obs_time & obs_time <= eta
)]

main <- x[, list(
  n_test = .N,
  upper_length_finite = mean(upper_length_finite, na.rm = TRUE),
  upper_length_truncated = mean(upper_length_truncated, na.rm = TRUE),
  inf_upper = mean(upper_inf)
), by = list(landmark, method_label)]

events <- ev[, list(
  n_events = .N,
  upper_cov = mean(upper_cov_ok, na.rm = TRUE),
  lower_trunc_cov = mean(lower_trunc_cov_ok, na.rm = TRUE)
), by = list(landmark, method_label)]

ans <- merge(main, events, by = c("landmark", "method_label"), all.x = TRUE)
ans[, order_method := match(method_label, c("Naive PRC", "Conformal PRC", "Conformal baseline Cox"))]
data.table::setorder(ans, landmark, order_method)
ans[, order_method := NULL]

data.table::fwrite(ans, out_file)
print(ans)
cat("\nSaved corrected ADNI aggregation to:", normalizePath(out_file, mustWork = FALSE), "\n")
cat("Use upper_length_finite to reproduce the historical length convention, or upper_length_truncated if the table note explicitly defines truncation at eta.\n")
