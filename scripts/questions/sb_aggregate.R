# =====================================================================================================================
# STEP 1 USE OF THIS SCRIPT: pool all step-1 data sets into the tables of README.md
#   Rscript scripts/questions/sb_aggregate.R --question SB1_one_record   (INLA's default prior, 85 data sets)
#   Rscript scripts/questions/sb_aggregate.R --question SB1_pc           (the same data refitted with PC priors, 80)
# Its output is results/step1_default_prior/ and results/step1_pc_prior/ of this repository.
# A fit counts as "collapsed" when sigma2_g < 0.001: it then sits at the bottom of INLA's default prior (about 1e-4)
# and says nothing about the genetic variance. In step 1 no pe model is fitted, so the WAIC reference is GHS (the true
# model); the report's generic header line "truth = GPHS" applies to step 3, where a pe exists.
# =====================================================================================================================

# Stage SB-AGGREGATE — genetic parameters of study SB (one breed; docs/questions/SB_single_breed.md), pooled.
#   Rscript scripts/questions/sb_aggregate.R --question SB_single_breed
# Reads runs_q/<question>/<cell>_rNN/{output/q_vars.csv, output/q_metrics.csv, truth/sim_settings.csv,
# truth/true_hyperparameters.csv}. A replicate counts only when q_extract.R has written its files (= all models fitted).
# Writes output/questions/<question>/{sb_by_rep.csv, sb_summary.csv, report.md}.
#   h2 = s2g / S and repeatability = (s2g + s2pe) / S, with S = the sum of the model's variance components (field: SPDE
#   variance). Truth: the same formulas with the true values.

args <- commandArgs(trailingOnly = TRUE)
Q <- args[match("--question", args) + 1]
if (is.na(Q)) stop("--question <name>")
root <- file.path("runs_q", Q)
out <- file.path("output/questions", Q)
dir.create(out, recursive = TRUE, showWarnings = FALSE)

runs <- list.dirs(root, recursive = FALSE)
runs <- runs[file.exists(file.path(runs, "output/q_vars.csv")) & file.exists(file.path(runs, "output/q_metrics.csv"))]
n_unfinished <- length(list.dirs(root, recursive = FALSE)) - length(runs)
if (!length(runs)) stop("No finished replicate under ", root)

term_key <- c("animal" = "s2g", "perm" = "s2pe", "herd" = "s2h", "field" = "s2s", "Gaussian observations" = "s2e")
rows <- lapply(runs, function(r) {
  v <- read.csv(file.path(r, "output/q_vars.csv"), stringsAsFactors = FALSE)
  m <- read.csv(file.path(r, "output/q_metrics.csv"), stringsAsFactors = FALSE)
  st <- read.csv(file.path(r, "truth/sim_settings.csv"), stringsAsFactors = FALSE)
  hp <- read.csv(file.path(r, "truth/true_hyperparameters.csv"))
  st <- setNames(st$value, st$setting); tv <- setNames(hp$value, hp$parameter)
  # Reference model for ΔWAIC / Δmlik: GPHS (the truth model) when it was fitted, else GHS (ladder steps 1 and 2).
  ref <- if ("GPHS" %in% v$model) "GPHS" else "GHS"
  ic <- function(mo, k) { x <- m$value[m$model == mo & m$metric == k & m$group == ""]; if (length(x)) x[1] else NA }
  do.call(rbind, lapply(unique(v$model), function(mo) {
    vm <- v[v$model == mo, ]
    get <- function(k, col = "mean") { i <- which(vm$what == "variance" & unname(term_key[vm$term]) %in% k); if (length(i)) vm[[col]][i[1]] else NA }
    rng <- vm$mean[grepl("^Range for field", vm$term)]
    met <- function(k, g = "") { x <- m$value[m$model == mo & m$metric == k & m$group == g]; if (length(x)) x[1] else NA }
    vals <- sapply(c("s2g", "s2pe", "s2h", "s2s", "s2e"), get)
    S <- sum(vals, na.rm = TRUE)
    data.frame(cell = sub("_r[0-9]+$", "", basename(r)), rep = as.integer(sub(".*_r", "", basename(r))), model = mo,
               ncow = as.numeric(st[["ncow"]]), p1 = as.numeric(st[["p1"]]), p1_realised = as.numeric(st[["p1_realised"]]),
               cows_single = as.numeric(st[["cows_single_share"]]), n_herd = as.numeric(st[["n_herd"]]),
               n_rec = as.numeric(st[["n_rec"]]), var_s_herds = tv[["realised_var_s_herds"]],
               s2g = vals[["s2g"]], s2pe = vals[["s2pe"]], s2h = vals[["s2h"]], s2s = vals[["s2s"]], s2e = vals[["s2e"]],
               range = if (length(rng)) rng[1] else NA,
               s2g_lo = get("s2g", "q025"), s2g_hi = get("s2g", "q975"),
               s2pe_lo = get("s2pe", "q025"), s2pe_hi = get("s2pe", "q975"),
               s2h_lo = get("s2h", "q025"), s2h_hi = get("s2h", "q975"),
               h2 = vals[["s2g"]] / S, rep_cow = sum(vals[c("s2g", "s2pe")], na.rm = TRUE) / S,
               t_s2g = tv[["sigma2_g"]], t_s2pe = tv[["sigma2_pe"]], t_s2h = tv[["sigma2_h"]], t_s2s = tv[["sigma2_s"]],
               t_s2e = tv[["sigma2_e"]], t_range = tv[["range_km"]],
               r_ebv = met("r_ebv_a0", "cows"), slope = met("slope_a0_on_ebv", "cows"),
               ref = ref, dwaic = ic(mo, "waic") - ic(ref, "waic"), dmlik = ic(mo, "mlik") - ic(ref, "mlik"),
               minutes = met("minutes"),
               # σ̂²g < 0.001: the fit sits at the prior's lower end (σ²g "collapsed"); its EBVs are ≈ 0
               collapsed = !is.na(vals[["s2g"]]) && vals[["s2g"]] < 1e-3)
  }))
})
d <- do.call(rbind, rows)
St <- d$t_s2g + d$t_s2pe + d$t_s2h + d$t_s2s + d$t_s2e
d$t_h2 <- d$t_s2g / St
d$t_rep <- (d$t_s2g + d$t_s2pe) / St
write.csv(d, file.path(out, "sb_by_rep.csv"), row.names = FALSE)

model_order <- c("G", "GH", "GS", "GHS", "GP", "GPH", "GPS", "GPHS")
models <- intersect(model_order, unique(d$model))
cells <- unique(d$cell[order(d$ncow, d$p1)])
ms <- function(x, dig = 3) { x <- x[!is.na(x)]; if (!length(x)) "–" else if (length(x) < 2) formatC(mean(x), format = "f", digits = dig) else
  sprintf(paste0("%.", dig, "f (%.", dig, "f)"), mean(x), sd(x)) }
nrep_of <- function(cl) length(unique(d$rep[d$cell == cl]))
table_of <- function(col, dig = 3, f = identity) {
  hdr <- c(paste0("| cell | reps | ", paste(models, collapse = " | "), " |"), paste0("|", strrep("---|", length(models) + 2)))
  body <- vapply(cells, function(cl) paste0("| ", cl, " | ", nrep_of(cl), " | ",
    paste(vapply(models, function(mo) ms(f(d[d$cell == cl & d$model == mo, ])[[col]], dig), ""), collapse = " | "), " |"), "")
  c(hdr, body)
}
cover <- function(x, lo, hi, t) if (all(is.na(x))) NA else mean(t >= lo & t <= hi, na.rm = TRUE)

summ <- do.call(rbind, lapply(cells, function(cl) do.call(rbind, lapply(models, function(mo) {
  x <- d[d$cell == cl & d$model == mo, ]
  if (!nrow(x)) return(NULL)
  data.frame(cell = cl, model = mo, n_rep = nrow(x), ncow = x$ncow[1], p1 = x$p1[1],
             s2g = mean(x$s2g), s2g_sd = sd(x$s2g), s2g_bias_pct = 100 * (mean(x$s2g) / x$t_s2g[1] - 1),
             s2g_cov = cover(x$s2g, x$s2g_lo, x$s2g_hi, x$t_s2g),
             s2pe = mean(x$s2pe), s2pe_sd = sd(x$s2pe), s2pe_cov = cover(x$s2pe, x$s2pe_lo, x$s2pe_hi, x$t_s2pe),
             s2h = mean(x$s2h), s2h_sd = sd(x$s2h), s2h_cov = cover(x$s2h, x$s2h_lo, x$s2h_hi, x$t_s2h),
             s2s = mean(x$s2s), range = mean(x$range), s2e = mean(x$s2e),
             h2 = mean(x$h2), h2_sd = sd(x$h2), rep_cow = mean(x$rep_cow),
             r_ebv = mean(x$r_ebv), n_collapsed = sum(x$collapsed),
             slope = if (all(x$collapsed)) NA else mean(x$slope[!x$collapsed]), dwaic = mean(x$dwaic), dmlik = mean(x$dmlik), minutes = mean(x$minutes),
             cor_pe_h = if (nrow(x) > 2 && all(!is.na(x$s2pe)) && all(!is.na(x$s2h))) cor(x$s2pe, x$s2h) else NA)
}))))
write.csv(summ, file.path(out, "sb_summary.csv"), row.names = FALSE)

tr <- d[1, ]
ref_lab <- paste(unique(d$ref), collapse = " / ")
design <- aggregate(cbind(n_herd, p1_realised, cows_single, n_rec, var_s_herds) ~ cell, data = d[d$model == models[1], ], FUN = mean)
design <- design[match(cells, design$cell), ]
lines <- c(paste0("# ", Q, " — genetic parameters (pooled)"), "",
  sprintf("Generated by `scripts/questions/sb_aggregate.R` on %s from %d finished replicates%s. SIMULATED data.",
          format(Sys.time(), "%Y-%m-%d %H:%M"), length(runs),
          if (n_unfinished > 0) sprintf(" (%d replicate folders not finished yet: INTERIM)", n_unfinished) else ""),
  "Entries: mean (sd over replicates). Models: Isidore's G, GH, GS, GHS and his pe series GP, GPH, GPS, GPHS (truth = GPHS).", "",
  sprintf("**Truth:** σ²g %.3f, σ²pe %.3f, σ²h %.3f, σ²s %.3f (SPDE), range %g km, σ²e %.3f → h² %.3f, repeatability %.3f.",
          tr$t_s2g, tr$t_s2pe, tr$t_s2h, tr$t_s2s, tr$t_range, tr$t_s2e, tr$t_h2, tr$t_rep), "",
  "## Design as realised", "", "| cell | herds | single-cow herds | cows in single-cow herds | records | field variance across herds |",
  "|---|---|---|---|---|---|",
  sprintf("| %s | %.0f | %.1f %% | %.1f %% | %.0f | %.3f |", design$cell, design$n_herd, 100 * design$p1_realised,
          100 * design$cows_single, design$n_rec, design$var_s_herds), "",
  sprintf("## Genetic variance σ²g (truth %.3f)", tr$t_s2g), "", table_of("s2g"), "",
  "Relative bias of σ²g (%), and share of 95 % intervals covering the truth:", "",
  c(paste0("| cell | ", paste(models, collapse = " | "), " |"), paste0("|", strrep("---|", length(models) + 1)),
    vapply(cells, function(cl) paste0("| ", cl, " | ", paste(vapply(models, function(mo) {
      s <- summ[summ$cell == cl & summ$model == mo, ]; if (!nrow(s)) "–" else sprintf("%+.0f %% (cov %.0f %%)", s$s2g_bias_pct, 100 * s$s2g_cov) }, ""),
      collapse = " | "), " |"), "")), "",
  sprintf("## Heritability h² (truth %.3f)", tr$t_h2), "", table_of("h2"), "",
  sprintf("## Permanent environment σ²pe (truth %.3f)", tr$t_s2pe), "", table_of("s2pe"), "",
  sprintf("## Herd σ²h (truth %.3f)", tr$t_s2h), "", table_of("s2h"), "",
  "GPHS: correlation of σ̂²pe and σ̂²h over replicates (E3: negative = the two trade off):", "",
  c("| cell | cor(σ̂²pe, σ̂²h) | cov σ²pe | cov σ²h |", "|---|---|---|---|",
    vapply(cells, function(cl) { s <- summ[summ$cell == cl & summ$model == "GPHS", ]
      if (!nrow(s)) sprintf("| %s | – | – | – |", cl) else sprintf("| %s | %s | %.0f %% | %.0f %% |", cl,
        if (is.na(s$cor_pe_h)) "–" else sprintf("%.2f", s$cor_pe_h), 100 * s$s2pe_cov, 100 * s$s2h_cov) }, "")), "",
  sprintf("## Repeatability (σ²g + σ²pe) / Σ (truth %.3f)", tr$t_rep), "", table_of("rep_cow"), "",
  sprintf("## Spatial variance σ²s, SPDE parameter (truth %.3f)", tr$t_s2s), "", table_of("s2s"), "",
  sprintf("## Range, km (truth %g)", tr$t_range), "", table_of("range", 1), "",
  sprintf("## Residual σ²e (truth %.3f)", tr$t_s2e), "", table_of("s2e"), "",
  "## EBV accuracy of the cows, r(EBV, a)", "", table_of("r_ebv"), "",
  "## Fits where σ̂²g collapsed to ≈ 0 (σ̂²g < 0.001: the fit sits at the prior's lower end; its EBVs are ≈ 0)", "",
  c(paste0("| cell | ", paste(models, collapse = " | "), " |"), paste0("|", strrep("---|", length(models) + 1)),
    vapply(cells, function(cl) paste0("| ", cl, " | ", paste(vapply(models, function(mo) {
      x <- d[d$cell == cl & d$model == mo, ]; if (!nrow(x)) "–" else sprintf("%d of %d", sum(x$collapsed), nrow(x)) }, ""),
      collapse = " | "), " |"), "")), "",
  "## EBV dispersion: slope of a on EBV (1 = right scale; < 1 = over-dispersed; > 1 = over-shrunk). Collapsed fits left out", "",
  table_of("slope", f = function(x) x[!x$collapsed, ]), "",
  sprintf("## ΔWAIC against %s (positive = worse than %s)", ref_lab, ref_lab), "", table_of("dwaic", 1), "",
  sprintf("## Δ log marginal likelihood against %s (positive = better than %s)", ref_lab, ref_lab), "", table_of("dmlik", 1), "",
  "## Fit time, minutes", "", table_of("minutes", 0), "")

# ---- cross-validation (q_validate.R --cv kfold): pooled rows, every cow predicted once ------------------
cv_done <- vapply(runs, function(r) {
  lg <- file.path(r, "output/q_validate_kfold_log.txt")
  file.exists(lg) && any(grepl("Q-VALIDATE finished", readLines(lg, warn = FALSE)))
}, logical(1))
if (any(cv_done)) {
  cvd <- do.call(rbind, lapply(runs[cv_done], function(r) {
    x <- read.csv(file.path(r, "output/q_validate_kfold.csv"), stringsAsFactors = FALSE)
    x <- x[x$mask == "pooled", ]
    x$cell <- sub("_r[0-9]+$", "", basename(r)); x$rep <- as.integer(sub(".*_r", "", basename(r))); x
  }))
  write.csv(cvd, file.path(out, "sb_cv_by_rep.csv"), row.names = FALSE)
  cv_models <- intersect(model_order, unique(cvd$model))
  cv_cells <- cells[cells %in% cvd$cell]
  tab_cv <- function(metric, grp = "all") {
    hdr <- c(paste0("| cell | reps | ", paste(cv_models, collapse = " | "), " |"), paste0("|", strrep("---|", length(cv_models) + 2)))
    body <- vapply(cv_cells, function(cl) paste0("| ", cl, " | ", length(unique(cvd$rep[cvd$cell == cl])), " | ",
      paste(vapply(cv_models, function(mo) ms(cvd[[metric]][cvd$cell == cl & cvd$model == mo & cvd$group == grp]), ""), collapse = " | "),
      " |"), "")
    c(hdr, body)
  }
  lines <- c(lines, "# Cross-validation (5-fold)", "",
    sprintf("%d replicates validated. The cows are split at random into 5 folds; each fold (all records of its cows) is masked once and predicted from a fit to the other 4, so every cow is predicted once. Pooled over all cows.", sum(cv_done)), "",
    "## True accuracy of the PREDICTED EBVs of the masked cows (no own records), r(EBV, a)", "", tab_cv("r_ebv_a"), "",
    "## For contrast: r(EBV, a) of the training cows (with own records) in the same fits, mean over folds", "", tab_cv("r_ebv_a_train"), "",
    "## Houaga's accuracy cor(y, ŷ) of the masked records", "", tab_cv("cor_y_yhat"), "",
    "## Its genetic part cor(y, fixed + EBV)", "", tab_cv("cor_y_gen"), "",
    "## Its environmental part cor(y, herd + field)", "", tab_cv("cor_y_env"), "",
    "## cor(y, ŷ) for masked cows whose herd still has training records", "", tab_cv("cor_y_yhat", "herd_seen"), "",
    "## cor(y, ŷ) for masked cows whose herd has NO training records (e.g. single-cow herds)", "", tab_cv("cor_y_yhat", "herd_unseen"), "",
    "## Best possible genetic prediction cor(y, Xb + a) of the masked records", "", tab_cv("ceil_gen"), "")
}
writeLines(lines, file.path(out, "report.md"), useBytes = TRUE)
cat("wrote", file.path(out, "report.md"), "from", length(runs), "replicates\n")
