# =====================================================================================================================
# STEP 1 USE OF THIS SCRIPT: score the fitted models against the simulated truth for one data set
#   Rscript scripts/questions/q_extract.R --run runs_q/SB1_one_record/n1000_p70_r01
#
# WHY: only a simulation knows the true breeding values and variances, so only here can we measure
#   - the bias of each variance component (estimate / truth - 1) and whether its 95 % interval covers the truth;
#   - accuracy r(EBV, a): correlation of the estimated with the TRUE breeding values of the cows (the ranking);
#   - the slope of a on EBV: the scale of the EBVs (1 = right; above 1 = EBVs shrunk too much, sigma2_g too low;
#     below 1 = EBVs too spread, sigma2_g too high);
#   - WAIC and the marginal likelihood of each model (the true model, GHS, should be best).
# Step 1 uses the cow-level metrics and the variance table. The bull, G x E and breed metrics below belong to other
# studies of the project and are empty or unused in step 1 (in study SB the bulls are not genotyped).
# Writes output/q_vars.csv and output/q_metrics.csv, which sb_aggregate.R pools.
# =====================================================================================================================

# Stage Q-EXTRACT — decision-oriented metrics for one question-study replicate.
#   Rscript scripts/questions/q_extract.R --run runs_q/Q1Q2/L1s_D4_r01
# Reads truth/ and output/q_fit_summaries.rds; writes output/q_metrics.csv (long format: model, metric, value)
# and output/q_vars.csv (variance components vs truth).
#
# Metrics (see docs/expansion_plan_realistic_questions.md, Section 1):
#   cows      r(EBV, a0), slope a0 on EBV, bias trend over farm-environment quintiles
#   proven AI bulls / test bulls   Spearman(EBV, a0), top-k overlap, selection loss,
#             "environment leakage" = slope of (EBV - a0) on the mean environment of the bull's daughters
#   G x E     accuracy of a(E = -1) and a(E = +1) for bulls and cows, re-ranking, rg / slope variance
#   breed x environment   cor(estimated breed slope at herds, true breed slope)
#   space     r(s_hat, s) within regions, spread of regional offsets
#   fit       WAIC, minutes

source("scripts/sim/_paths.R")
args <- commandArgs(trailingOnly = TRUE)
S <- readRDS(file.path(out_dir, "q_fit_summaries.rds"))
ta <- read.csv(file.path(truth_dir, "true_animals.csv"), stringsAsFactors = FALSE)
th <- read.csv(file.path(truth_dir, "true_herds.csv"), stringsAsFactors = FALSE)
tb <- read.csv(file.path(truth_dir, "true_bulls.csv"), stringsAsFactors = FALSE)
hp <- read.csv(file.path(truth_dir, "true_hyperparameters.csv"))
tv <- setNames(hp$value, hp$parameter)
settings <- read.csv(file.path(truth_dir, "sim_settings.csv"), stringsAsFactors = FALSE)
an <- read.csv(file.path(data_dir, "animals.csv"), stringsAsFactors = FALSE)

rows <- list()
add <- function(model, metric, value, group = "") rows[[length(rows) + 1]] <<- data.frame(model = model, group = group, metric = metric, value = as.numeric(value))
spear <- function(x, y) suppressWarnings(cor(x, y, method = "spearman"))
topk <- function(est, tru, k) length(intersect(order(-est)[seq_len(k)], order(-tru)[seq_len(k)])) / k
sel_loss <- function(est, tru, k) mean(sort(tru, decreasing = TRUE)[seq_len(k)]) - mean(tru[order(-est)[seq_len(k)]])
leak <- function(err, env) if (sum(is.finite(env)) > 3) unname(coef(lm(err ~ env))[2]) else NA

vars_out <- list()
for (m in names(S)) {
  s <- S[[m]]
  ids <- s$animal_ids
  # animal effect: ungrouped (n) or grouped by environment class (3n, class-major)
  n_a <- length(ids)
  grouped <- nrow(s$animal) == 3 * n_a
  ebv <- if (grouped) matrix(s$animal$mean, ncol = 3) else matrix(s$animal$mean, ncol = 1)
  rownames(ebv) <- ids
  ebv_avg <- rowMeans(ebv)                          # average over classes = EBV for the average environment
  a0 <- ta$a0_true[match(ids, ta$animal_id)]
  role <- ta$role[match(ids, ta$animal_id)]
  cow <- role == "cow"

  add(m, "minutes", s$seconds / 60)
  add(m, "waic", s$waic)
  add(m, "mlik", s$mlik)
  add(m, "r_ebv_a0", cor(ebv_avg[cow], a0[cow]), "cows")
  add(m, "slope_a0_on_ebv", coef(lm(a0[cow] ~ ebv_avg[cow]))[2], "cows")
  # bias trend over farm-environment quintiles (cows)
  fe <- ta$farm_E[match(ids[cow], ta$animal_id)]
  q5 <- cut(fe, quantile(fe, 0:5 / 5), include.lowest = TRUE, labels = FALSE)
  bias_q <- tapply(ebv_avg[cow] - a0[cow], q5, mean)
  add(m, "bias_trend_per_quintile", coef(lm(bias_q ~ seq_along(bias_q)))[2], "cows")
  add(m, "bias_top_minus_bottom_quintile", bias_q[5] - bias_q[1], "cows")

  for (grp in c("AI_proven", "AI_test")) {
    b <- tb[tb$type == grp & tb$n_daughters > 0, ]
    if (nrow(b) < 5) next
    e <- ebv_avg[match(b$id, ids)]
    k <- if (grp == "AI_proven") 10 else 5
    add(m, "spearman", spear(e, b$a0_true), grp)
    add(m, "pearson", cor(e, b$a0_true), grp)
    add(m, paste0("top", k, "_overlap"), topk(e, b$a0_true, k), grp)
    add(m, paste0("selection_loss_top", k), sel_loss(e, b$a0_true, k), grp)
    err <- e - b$a0_true
    add(m, "env_leakage_slope_E", leak(err - mean(err), b$mean_E_daughters), grp)
    add(m, "env_leakage_slope_driver", leak(err - mean(err), b$mean_zdriver_daughters), grp)
    # Partial leakage: slope of the EBV error on the daughters' environment HOLDING TRUE MERIT FIXED. When good bulls
    # go to good farms (lambda > 0), environment and merit are correlated, and plain shrinkage of EBVs (error ~ -merit)
    # would show up as a negative slope on environment even without any leakage.
    pl <- function(env) if (sd(env, na.rm = TRUE) > 0.05) unname(coef(lm(err ~ env + b$a0_true))[2]) else NA
    add(m, "env_leakage_partial_E", pl(b$mean_E_daughters), grp)
    add(m, "env_leakage_partial_driver", pl(b$mean_zdriver_daughters), grp)
    add(m, "sd_mean_E_daughters", sd(b$mean_E_daughters), grp)
    add(m, "cor_merit_env_daughters", cor(b$a0_true, b$mean_zdriver_daughters), grp)
    add(m, "slope_a0_on_ebv", coef(lm(b$a0_true ~ e))[2], grp)
  }

  # G x E: environment-specific breeding values a(E) = a0 + a1 E at E = -1 and +1 (bulls), and at the cow's own farm
  if (any(ta$a1_true != 0)) {
    a1 <- ta$a1_true[match(ids, ta$animal_id)]
    a_lo <- a0 - a1; a_hi <- a0 + a1
    e_lo <- e_hi <- ebv_avg
    if (!is.null(s$slope)) {                       # sire random regression: bull slope = half the bull's own slope
      sl <- setNames(s$slope$mean, s$slope_ids)
      ib <- match(s$slope_ids, ids)
      e_lo[ib] <- ebv_avg[ib] - 2 * sl; e_hi[ib] <- ebv_avg[ib] + 2 * sl
      add(m, "r_slope_bulls", cor(sl, a1[ib]), "AI_bulls")
    } else if (grouped) {
      e_lo <- ebv[, 1]; e_hi <- ebv[, 3]
    } else if (!is.null(s$gk)) {                   # sire x location kernel: mean cell effect in low / high E_hat tertile
      gc <- s$gk_cells; gc$eff <- s$gk
      cut3 <- quantile(s$Ehat_herd, c(1/3, 2/3))
      for (b in unique(gc$sire)) {
        g <- gc[gc$sire == b, ]
        ii <- match(b, ids)
        lo <- g$eff[g$Ehat <= cut3[1]]; hi <- g$eff[g$Ehat >= cut3[2]]
        if (length(lo)) e_lo[ii] <- ebv_avg[ii] + 2 * mean(lo)
        if (length(hi)) e_hi[ii] <- ebv_avg[ii] + 2 * mean(hi)
      }
    }
    for (grp in c("AI_proven", "AI_test")) {
      b <- tb[tb$type == grp & tb$n_daughters > 0, ]
      if (nrow(b) < 5) next
      ii <- match(b$id, ids)
      add(m, "r_low_env", cor(e_lo[ii], a_lo[ii]), grp)
      add(m, "r_high_env", cor(e_hi[ii], a_hi[ii]), grp)
      add(m, "true_rerank_spearman_low_high", spear(a_lo[ii], a_hi[ii]), grp)
      add(m, "est_rerank_spearman_low_high", spear(e_lo[ii], e_hi[ii]), grp)
      k <- min(10, floor(length(ii) / 3))
      add(m, "topk_overlap_low_env", topk(e_lo[ii], a_lo[ii], k), grp)
      add(m, "topk_overlap_high_env", topk(e_hi[ii], a_hi[ii], k), grp)
    }
    # cows: breeding value at their own farm environment
    cid <- ids[cow]
    a_own <- ta$a_own_env[match(cid, ta$animal_id)]
    cherd <- ta$herd[match(cid, ta$animal_id)]
    own <- ebv_avg[cow]
    if (!is.null(s$slope)) {
      csire <- an$sire[match(cid, an$animal_id)]
      has <- !is.na(csire) & csire %in% s$slope_ids
      own[has] <- own[has] + sl[csire[has]] * s$Ehat_herd[match(cherd[has], s$herd_ids)]
    } else if (grouped) {
      own <- ebv[cbind(which(cow), s$eclass_herd[match(cherd, s$herd_ids)])]
    } else if (!is.null(s$gk)) {
      csire <- an$sire[match(cid, an$animal_id)]
      key <- paste(csire, cherd); ck <- paste(s$gk_cells$sire, s$gk_cells$herd)
      has <- key %in% ck
      own[has] <- own[has] + s$gk[match(key[has], ck)]
    }
    add(m, "r_own_env_bv", cor(own, a_own), "cows")
  }

  # spatial effects at herds
  if (!is.null(s$field_herds)) {
    f <- if (is.matrix(s$field_herds)) rowMeans(s$field_herds) else s$field_herds
    st <- th$s_true[match(s$herd_ids, th$herd)]
    rg <- th$region[match(s$herd_ids, th$herd)]
    add(m, "r_shat_s", cor(f, st), "herds")
    add(m, "r_shat_s_within_region", cor(f - ave(f, rg), st - ave(st, rg)), "herds")
    add(m, "var_shat_over_var_s", var(f) / var(st), "herds")          # S3: how much of the spatial variance is captured
    off <- tapply(f - st, rg, mean)
    add(m, "region_offset_spread", diff(range(off)), "herds")
    if (is.matrix(s$field_herds) && all(paste0("Y", 1:4) %in% names(th))) {
      sty <- as.matrix(th[match(s$herd_ids, th$herd), paste0("Y", 1:4)])
      add(m, "r_shat_s_by_year", cor(as.numeric(s$field_herds), as.numeric(sty)), "herds")
    }
  }
  # S3: the TOTAL farm environment (herd + ward + spatial) against its estimate (herd + field), whatever its shape
  # [[ ]] not $: s$herd would partially match s$herd_ids in models without a herd effect (smoke test 4 Nov 2026)
  if (!is.null(s[["herd"]]) || !is.null(s[["field_herds"]])) {
    env_hat <- rep(0, length(s$herd_ids))
    if (!is.null(s[["herd"]])) env_hat <- env_hat + s[["herd"]]
    if (!is.null(s[["field_herds"]])) env_hat <- env_hat + (if (is.matrix(s[["field_herds"]])) rowMeans(s[["field_herds"]]) else s[["field_herds"]])
    thm <- th[match(s$herd_ids, th$herd), ]
    add(m, "r_env_hat_env", cor(env_hat, thm$h_true + thm$ward_true + thm$s_true), "herds")
  }
  if (!is.null(s$bfield_herds)) {
    bt <- th$breed_slope_true[match(s$herd_ids, th$herd)]
    if (sd(bt) > 0) add(m, "r_breed_slope", cor(s$bfield_herds, bt), "herds")
    add(m, "sd_breed_slope_hat", sd(s$bfield_herds), "herds")
  }

  # variance components vs truth
  v <- s$vars
  truth_of <- function(term) switch(term,
    "animal" = tv[["sigma2_g"]], "perm" = tv[["sigma2_pe"]], "herd" = tv[["sigma2_h"]],
    "Gaussian observations" = tv[["sigma2_e"]], "field" = tv[["sigma2_s"]], "slope" = tv[["sigma2_slope"]],
    NA_real_)
  v$truth <- vapply(seq_len(nrow(v)), function(i) if (v$what[i] == "variance") truth_of(v$term[i]) else
    if (grepl("^Range for field", v$term[i])) tv[["range_km"]] else
      if (grepl("^GroupRho for animal", v$term[i])) NA_real_ else NA_real_, numeric(1))
  v$model <- m
  vars_out[[m]] <- v
}
# S1 / S4: paired comparisons with the free-range GPHS of the same replicate (WAIC and marginal-likelihood
# differences; agreement of EBVs). Positive dwaic = worse than GPHS; positive dmlik = better than GPHS.
if ("GPHS" %in% names(S)) {
  ref <- S[["GPHS"]]
  ref_ebv <- setNames(ref$animal$mean, ref$animal_ids)
  cow_ids <- intersect(ta$animal_id[ta$role == "cow"], ref$animal_ids)          # (test-region fits hold a subset)
  prov_ids <- intersect(tb$id[tb$type == "AI_proven" & tb$n_daughters > 0], ref$animal_ids)
  for (m in setdiff(names(S), "GPHS")) {
    s <- S[[m]]
    add(m, "dwaic_vs_GPHS", s$waic - ref$waic)
    add(m, "dmlik_vs_GPHS", s$mlik - ref$mlik)
    if (nrow(s$animal) == length(s$animal_ids)) {
      e_m <- setNames(s$animal$mean, s$animal_ids)
      add(m, "r_ebv_vs_GPHS", cor(e_m[cow_ids], ref_ebv[cow_ids]), "cows")
      if (length(prov_ids) >= 5) add(m, "r_ebv_vs_GPHS", cor(e_m[prov_ids], ref_ebv[prov_ids]), "AI_proven")
    }
  }
}
met <- do.call(rbind, rows)
sw <- as.data.frame(as.list(setNames(settings$value, settings$setting)), stringsAsFactors = FALSE)
met <- cbind(sw[rep(1, nrow(met)), , drop = FALSE], met, row.names = NULL)
write.csv(met, file.path(out_dir, "q_metrics.csv"), row.names = FALSE)
vv <- do.call(rbind, vars_out)
write.csv(vv, file.path(out_dir, "q_vars.csv"), row.names = FALSE)
cat("Q-EXTRACT finished:", nrow(met), "metrics,", nrow(vv), "variance rows ->", out_dir, "\n")
