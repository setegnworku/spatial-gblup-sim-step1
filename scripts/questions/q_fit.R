# =====================================================================================================================
# STEP 1 USE OF THIS SCRIPT: fit G, GH, GS and GHS to one simulated data set
#   Rscript scripts/questions/q_fit.R --run runs_q/SB1_one_record/n1000_p70_r01 --qmodels G,GH,GS,GHS --priors isidore
#
# WHY THESE FOUR MODELS: they are the models of Houaga et al. (2026). In step 1 the true model is GHS (one record per
#   cow, no pe). G, GH and GS each leave out a true term on purpose, to show where that term's variance goes.
#     G    animal                     GH   animal + herd
#     GS   animal + field             GHS  animal + herd + field   (= truth in step 1)
#   "animal" = genomic breeding value with covariance sigma2_g * G (VanRaden G of the cows, G* = 0.99 G + 0.01 I).
#
# PRIORS (step 1 was fitted twice, on the same data sets):
#   --priors isidore  INLA's DEFAULT prior, loggamma(1, 5e-5) on each precision (= inverse-gamma(1, 5e-5) on the
#                     variance). It is what Houaga's script used because it sets no prior; the flag is named after that
#                     script, it is not a prior anyone chose. On the variance scale it puts 61 % of its probability
#                     below 0.0001 and 99.9 % below 0.06.
#   --priors pc       penalised-complexity (PC) priors (Simpson, Rue, Riebler, Martins and Sorbye 2017, Statistical
#                     Science 32: 1-28): an exponential prior on the standard deviation that shrinks towards "no effect"
#                     at a constant rate, set by P(sd > U) = alpha. Here P(sd > 0.5) = 0.5 for animal and herd and
#                     P(sd > 1) = 0.5 for the residual. It puts 1.4 % of its probability below a variance of 0.0001.
#   In both, the field keeps Houaga's PC-Matern prior: P(sd > 0.5) = 0.5, P(range < 50 km) = 0.8.
#   Integration: empirical Bayes (int.strategy = "eb"). On one data set it agreed with INLA's default ("auto") to three
#   decimals, so EB does not cause the collapses reported in README.md.
#
# EXPECTED (written before the results)
#   7 Nov 2026, default prior: GHS unbiased for sigma2_g; G and GS put the environment into sigma2_g; GH puts the
#                field into the herd; WAIC prefers GHS; GHS accuracy about 0.35-0.45.
#   9 Nov 2026, PC refit: PC1 the collapses vanish; PC2 GHS at 1,000 cows rises to about 0.05-0.06; PC3 GHS at 2,000
#                cows lands within 1 SE of 0.06 if the prior caused its low bias; PC4 WAIC better, accuracy unchanged.
#
# The model library below also serves later steps and other studies of the project; step 1 uses only G, GH, GS and GHS.
# =====================================================================================================================

# Stage Q-FIT — fit a chosen set of models to one question-study replicate and save compact summaries.
#
#   Rscript scripts/questions/q_fit.R --run runs_q/Q1Q2/L1s_D4_r01 --qmodels GH,GHS,GPHS,GPHS_sg
#   options: --test-region R4  (fit one region only, for quick local tests)   --keep-fits
#
# All models use the SAME priors (PC priors: P(sd > 0.5) = 0.5 for animal, pe, herd and slope effects;
# P(sd > 1) = 0.5 residual; Houaga's PC-Matern for fields) and empirical-Bayes integration, so that the
# differences between models come from their structure. Fixed part = Houaga's full fixed effects.
#
# Model library (terms besides the fixed part):
#   G         animal                                (Houaga's base model; S2)
#   GS        animal + field                        (Houaga; S2)
#   GH        animal + herd                         (no spatial term)
#   GHS       animal + herd + field                 (Houaga's full model, no pe)
#   GP, GPH, GPS   animal + pe [+ herd | + field]   (Isidore's pe series; study SB)
#   GPHS      animal + pe + herd + field            (reference model of the question studies)
#   GPHS_sg   GPHS + sire-group fixed effect        (genetic groups for bull origin; Q1)
#   GPHS_SVC  GPHS + exotic-proportion x second field  (spatially varying breed effect; Q3a)
#   GPHS_RN   GPHS + AI-sire slope x E_hat: sire random regression on the GPHS spatial environment E_hat = z(s_hat),
#             slopes correlated through G among AI bulls (two-step reaction norm; Q3b)
#   GPHS_ME   GPHS with the animal effect split into 3 environment classes (tertiles of E_hat), G x exchangeable (rg; Q3b)
#   GPHS_GK   GPHS + sire x herd effect with covariance G_bulls o M, M = Matern correlation of the herds at the GPHS
#             range: a "sire x location" kernel at fixed range (Q3b)
# Why the G x E terms sit at the SIRE level: a slope or G o M effect per COW is confounded with her herd and pe
# effects (72 % single-cow herds; most cows nearly unrelated). Found in the smoke test on 3 Oct 2026: a cow-level
# slope took 0.41 of variance and pushed herd and pe to ~0. AI bulls have many daughters across environments.
#   GPHS_ST   animal + pe + herd + field grouped by year, AR(1) across years    (space-time; Q5)
#   GPHS_R<km>  GPHS with the field's range FIXED at <km> (e.g. GPHS_R16.5): pcmatern prior.range = c(km, NA);
#             the spatial variance is still estimated. The "frozen kappa" of the BLUPF90/Mix99 route (S4)
# RN, ME and GK need GPHS first (E_hat, classes, range): it is fitted first automatically.
# --priors isidore (study SB): INLA's default loggamma(1, 5e-5) on every precision, as in Isidore's script; the field
# keeps his PC-Matern. Default --priors pc (all earlier studies).

source("scripts/sim/_fit_setup.R")    # data, Ginv, mesh, matern, fixed_formula (reads only data/)
suppressPackageStartupMessages(library(Matrix))
args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(flag, default = NULL) if (flag %in% args) args[match(flag, args) + 1] else default
qmodels <- strsplit(arg_value("--qmodels", "GH,GHS,GPHS"), ",")[[1]]
test_region <- arg_value("--test-region")
keep_fits <- "--keep-fits" %in% args
priors <- arg_value("--priors", "pc")
if (!priors %in% c("pc", "isidore")) stop("--priors must be pc or isidore")
inla.setOption(num.threads = Sys.getenv("INLA_THREADS", "4:1"))
fit_dir <- file.path(out_dir, "q_fits")
dir.create(fit_dir, recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(out_dir, "q_fit_log.txt")
cat("", file = log_file, append = FALSE)
say <- function(...) { txt <- paste0(...); cat(txt, "\n"); cat(txt, "\n", file = log_file, append = TRUE) }

G_full <- readRDS(file.path(data_dir, "G.rds"))
animals$sire_group[is.na(animals$sire_group)] <- "none"

# ---- optional: one region only (local smoke tests) ----------------------------------------------
if (!is.null(test_region)) {
  herds <- herds[herds$region == test_region, ]
  dat <- dat[dat$herd %in% herds$herd, ]
  keep_ids <- c(unique(dat$animal_id), animals$animal_id[animals$role != "cow"])
  G_full <- G_full[keep_ids, keep_ids]
  Ginv <- solve(G_full); Ginv <- (Ginv + t(Ginv)) / 2; dimnames(Ginv) <- dimnames(G_full)
  n_animal <- nrow(Ginv); n_herd <- nrow(herds)
  dat$cowI <- match(dat$animal_id, rownames(Ginv)); dat$herdI <- match(dat$herd, herds$herd)
  loc_herd <- as.matrix(herds[, c("x_km", "y_km")])
  mesh <- build_mesh(loc_herd, mesh_type_of(data_dir))
  matern <- inla.spde2.pcmatern(mesh, prior.sigma = c(sqrt(0.25), 0.5), prior.range = c(50, 0.8))
  for (v in c("tyrmn", "cyrsn", "lacgr")) dat[[v]] <- droplevels(dat[[v]])
  say("TEST MODE: region ", test_region, " only (", nrow(herds), " herds, ", nrow(dat), " records)")
}

# ---- extra covariates for the new models --------------------------------------------------------
dat$sgrp <- factor(animals$sire_group[match(dat$animal_id, animals$animal_id)])
p_cls <- animals$exotic_class_prop[match(dat$animal_id, animals$animal_id)]
dat$p_c <- p_cls - mean(p_cls)
dat$ygrp <- if ("ygrp" %in% names(dat)) as.integer(dat$ygrp) else 1L
# AI sires (genotyped, in G): index among the bulls; records of village-sired cows get weight 0 in sire terms
bull_ids <- animals$animal_id[animals$role != "cow" & animals$animal_id %in% rownames(Ginv)]
n_bull <- length(bull_ids)
Gb_inv <- NULL                                    # no genotyped bulls (study SB: G holds the cows only, as Isidore's)
if (n_bull > 0) { Gb_inv <- solve(G_full[bull_ids, bull_ids]); Gb_inv <- (Gb_inv + t(Gb_inv)) / 2 }
cow_sire <- animals$sire[match(dat$animal_id, animals$animal_id)]
dat$w_sire <- as.numeric(!is.na(cow_sire) & cow_sire %in% bull_ids)
dat$sireJ <- ifelse(dat$w_sire > 0, match(cow_sire, bull_ids), 1L)
B_herd_fit <- fm_basis(mesh, as.matrix(herds[, c("x_km", "y_km")]))
if (is.list(B_herd_fit)) B_herd_fit <- B_herd_fit$A

pc <- list(prec = list(prior = "pc.prec", param = c(0.5, 0.5)))
pc_slope <- list(prec = list(prior = "pc.prec", param = c(0.25, 0.5)))
family_control <- list(hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.5))))
if (priors == "isidore") { pc <- NULL; family_control <- list() }   # hyper = NULL -> INLA default loggamma(1, 5e-5)
matern2 <- inla.spde2.pcmatern(mesh, prior.sigma = c(sqrt(0.25), 0.5), prior.range = c(50, 0.8))
fx_txt <- function(sg = FALSE) {
  f <- paste(deparse(fixed_formula), collapse = "")
  if (sg) f <- paste(f, "+ sgrp")
  paste0("fixed_effects(main = ", f, ", model = \"fixed\")")
}
T_animal <- "animal(main = cowI, model = \"generic\", Cmatrix = Ginv, mapper = bru_mapper_index(n_animal), hyper = pc)"
T_perm <- "perm(main = cowI, model = \"iid\", mapper = bru_mapper_index(n_animal), hyper = pc)"
T_herd <- "herd(main = herdI, model = \"iid\", mapper = bru_mapper_index(n_herd), hyper = pc)"
T_field <- "field(main = cbind(x_km, y_km), model = matern)"
cmp <- function(terms) as.formula(paste("~ -1 +", paste(terms, collapse = " + ")), env = globalenv())

fit_one <- function(components, data = dat) {
  bru(components,
      bru_obs(formula = y ~ ., family = "gaussian", data = data, control.family = family_control),
      options = list(control.inla = list(int.strategy = "eb"),
                     control.compute = list(waic = TRUE, dic = TRUE), verbose = FALSE))
}

# Posterior summaries of variances from the hyperparameter marginals
var_summary <- function(fit) {
  mh <- fit$marginals.hyperpar
  do.call(rbind, lapply(names(mh), function(nm) tryCatch({
    m <- mh[[nm]]
    if (grepl("^Precision", nm)) {
      v <- inla.tmarginal(function(x) 1 / x, m)
      data.frame(term = sub("^Precision for (the )?", "", nm), what = "variance", mean = inla.emarginal(identity, v),
                 q025 = inla.qmarginal(0.025, v), q975 = inla.qmarginal(0.975, v))
    } else if (grepl("^Stdev", nm)) {
      v <- inla.tmarginal(function(x) x^2, m)
      data.frame(term = sub("^Stdev for ", "", nm), what = "variance", mean = inla.emarginal(identity, v),
                 q025 = inla.qmarginal(0.025, v), q975 = inla.qmarginal(0.975, v))
    } else {
      data.frame(term = nm, what = "other", mean = inla.emarginal(identity, m),
                 q025 = inla.qmarginal(0.025, m), q975 = inla.qmarginal(0.975, m))
    }
  }, error = function(e) {
    # degenerate marginal (e.g. the breed-field variance when its true variance is 0, Q3a B0_D3_r03 failed here
    # on 5 Oct 2026): keep the row with NA instead of losing the whole run
    say("   var_summary: marginal '", nm, "' not summarisable (", conditionMessage(e), ")")
    data.frame(term = sub("^(Precision|Stdev) for (the )?", "", nm), what = "unsummarised", mean = NA_real_, q025 = NA_real_, q975 = NA_real_)
  })))
}
herd_field <- function(fit, name = "field") {
  rf <- fit$summary.random[[name]]
  if (is.null(rf)) return(NULL)
  if (nrow(rf) == mesh$n) return(as.numeric(B_herd_fit %*% rf$mean))
  # grouped field (space-time): mesh$n x ngroup, column-major by group
  matrix(as.numeric(B_herd_fit %*% matrix(rf$mean, nrow = mesh$n)), nrow = nrow(herds))
}
summarise_fit <- function(fit, name, secs) {
  ra <- fit$summary.random$animal
  out <- list(model = name, seconds = secs, waic = fit$waic$waic, dic = fit$dic$dic, mlik = fit$mlik[1, 1],
              vars = var_summary(fit), animal_ids = rownames(Ginv), herd_ids = herds$herd)
  out$animal <- data.frame(ID = ra$ID, mean = ra$mean, sd = ra$sd)
  if (!is.null(fit$summary.random$perm)) out$perm <- fit$summary.random$perm$mean
  if (!is.null(fit$summary.random$herd)) out$herd <- fit$summary.random$herd$mean
  if (!is.null(fit$summary.random$slope)) out$slope <- data.frame(mean = fit$summary.random$slope$mean, sd = fit$summary.random$slope$sd)
  if (!is.null(fit$summary.random$gk)) out$gk <- fit$summary.random$gk$mean
  out$field_herds <- herd_field(fit, "field")
  out$bfield_herds <- herd_field(fit, "bfield")
  fe <- fit$summary.random$fixed_effects
  out$fixed <- data.frame(ID = as.character(fe$ID), mean = fe$mean, sd = fe$sd)
  out
}

# ---- which models, in which order ---------------------------------------------------------------
needs_gphs <- any(c("GPHS_RN", "GPHS_ME", "GPHS_GK") %in% qmodels)
order_models <- unique(c(if (needs_gphs) "GPHS", qmodels))
say("=== Q-FIT ", run_dir, " | models: ", paste(order_models, collapse = ", "), " | records ", nrow(dat),
    " | animals in G ", n_animal, " | herds ", n_herd, " | mesh ", mesh$n, " | priors ", priors,
    " | fixed ", paste(deparse(fixed_formula), collapse = ""))
results <- list()
for (m in order_models) {
  t0 <- Sys.time()
  if (m == "G") components <- cmp(c(fx_txt(), T_animal))
  if (m == "GS") components <- cmp(c(fx_txt(), T_animal, T_field))
  if (grepl("^GPHS_R[0-9.]+$", m)) {
    r_fix <- as.numeric(sub("^GPHS_R", "", m))
    matern_fix <- inla.spde2.pcmatern(mesh, prior.sigma = c(sqrt(0.25), 0.5), prior.range = c(r_fix, NA))   # NA = range fixed
    components <- cmp(c(fx_txt(), T_animal, T_perm, T_herd, "field(main = cbind(x_km, y_km), model = matern_fix)"))
  }
  if (m == "GP") components <- cmp(c(fx_txt(), T_animal, T_perm))
  if (m == "GPH") components <- cmp(c(fx_txt(), T_animal, T_perm, T_herd))
  if (m == "GPS") components <- cmp(c(fx_txt(), T_animal, T_perm, T_field))
  if (m == "GH") components <- cmp(c(fx_txt(), T_animal, T_herd))
  if (m == "GHS") components <- cmp(c(fx_txt(), T_animal, T_herd, T_field))
  if (m == "GPHS") components <- cmp(c(fx_txt(), T_animal, T_perm, T_herd, T_field))
  if (m == "GPHS_sg") components <- cmp(c(fx_txt(sg = TRUE), T_animal, T_perm, T_herd, T_field))
  if (m == "GPHS_SVC") components <- cmp(c(fx_txt(), T_animal, T_perm, T_herd, T_field,
                                           "bfield(main = cbind(x_km, y_km), weights = p_c, model = matern2)"))
  if (m == "GPHS_ST") components <- cmp(c(fx_txt(), T_animal, T_perm, T_herd,
                                          "field(main = cbind(x_km, y_km), model = matern, group = ygrp, ngroup = 4, control.group = list(model = \"ar1\"))"))
  if (m %in% c("GPHS_RN", "GPHS_ME", "GPHS_GK")) {
    ref <- results[["GPHS"]]
    # Environment descriptor = the estimated SPATIAL effect only. Using h_hat + s_hat was circular: with 72 %
    # single-cow herds h_hat is essentially the cow's own residual, and the slope then absorbed the herd
    # variance (smoke test 3 Oct 2026: herd 0.20 -> 0.01). s_hat is estimated from many neighbouring herds.
    E_raw <- ref$field_herds
    E_hat <- (E_raw - mean(E_raw)) / sd(E_raw)
    dat$Ehat <- E_hat[dat$herdI]
    dat$eclass <- as.integer(cut(dat$Ehat, quantile(E_hat, c(0, 1/3, 2/3, 1)), include.lowest = TRUE))
  }
  if (m == "GPHS_RN") components <- cmp(c(fx_txt(), T_animal, T_perm, T_herd, T_field,
    "slope(main = sireJ, weights = w_sire * Ehat, model = \"generic\", Cmatrix = Gb_inv, mapper = bru_mapper_index(n_bull), hyper = pc_slope)"))
  if (m == "GPHS_ME") components <- cmp(c(fx_txt(),
    "animal(main = cowI, model = \"generic\", Cmatrix = Ginv, mapper = bru_mapper_index(n_animal), hyper = pc, group = eclass, ngroup = 3, control.group = list(model = \"exchangeable\"))",
    T_perm, T_herd, T_field))
  if (m == "GPHS_GK") {
    rng <- ref$vars$mean[grepl("^Range", ref$vars$term)][1]          # GPHS posterior mean range (km)
    ai_rec <- dat$w_sire > 0
    cell_key <- paste(cow_sire[ai_rec], dat$herd[ai_rec])
    cells <- unique(cell_key)
    cell_sire <- sub(" .*", "", cells); cell_herd <- sub("^[^ ]+ ", "", cells)
    xy <- as.matrix(herds[match(cell_herd, herds$herd), c("x_km", "y_km")])
    D <- as.matrix(dist(xy)); kap <- sqrt(8) / rng
    M <- ifelse(D == 0, 1, (kap * D) * besselK(pmax(kap * D, 1e-12), 1))           # Matern nu = 1 correlation
    GK <- G_full[cell_sire, cell_sire] * M + diag(1e-4, length(cells))
    GKinv <- solve(GK); GKinv <- (GKinv + t(GKinv)) / 2
    dat$cellK <- 1L; dat$cellK[ai_rec] <- match(cell_key, cells)
    n_gk <- length(cells)
    components <- cmp(c(fx_txt(), T_animal, T_perm, T_herd, T_field,
      "gk(main = cellK, weights = w_sire, model = \"generic\", Cmatrix = GKinv, mapper = bru_mapper_index(n_gk), hyper = pc)"))
    say(sprintf("   sire x location kernel: %d sire-herd cells (%d AI bulls), range fixed at %.1f km (GPHS posterior mean)",
                n_gk, length(unique(cell_sire)), rng))
  }
  say("--- fitting ", m)
  fit <- fit_one(components, dat)
  secs <- as.numeric(Sys.time() - t0, units = "secs")
  s <- summarise_fit(fit, m, secs)
  if (m == "GPHS_ME") s$eclass_herd <- tapply(dat$eclass, dat$herdI, function(x) x[1])
  if (m %in% c("GPHS_RN", "GPHS_ME", "GPHS_GK")) s$Ehat_herd <- E_hat
  if (m == "GPHS_RN") s$slope_ids <- bull_ids
  if (m == "GPHS_GK") s$gk_cells <- data.frame(sire = cell_sire, herd = cell_herd, Ehat = E_hat[match(cell_herd, herds$herd)])
  results[[m]] <- s
  if (keep_fits) saveRDS(fit, file.path(fit_dir, paste0("fit_", m, ".rds")))
  say(sprintf("   %s: %.1f min | WAIC %.1f | %s", m, secs / 60, s$waic,
              paste(sprintf("%s %.3f", s$vars$term, s$vars$mean), collapse = "; ")))
  saveRDS(results, file.path(out_dir, "q_fit_summaries.rds"))      # saved after every model
  rm(fit); gc(verbose = FALSE)
}
say("Q-FIT finished: ", file.path(out_dir, "q_fit_summaries.rds"))
