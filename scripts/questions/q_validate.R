# =====================================================================================================================
# STEP 1 USE OF THIS SCRIPT: 5-fold cross-validation (replicates 1-5 of the 1,000- and 2,000-cow cells; 40 data sets)
#   Rscript scripts/questions/q_validate.R --run runs_q/SB1_one_record/n1000_p70_r01 --cv kfold --qmodels G,GH,GS,GHS --priors isidore
#
# WHY: a breeding programme needs EBVs of animals WITHOUT their own records (young heifers, cows in new herds). The
#   cows are split at random into 5 groups; each group's records are hidden in turn and predicted from a fit to the
#   other four, so every cow is predicted once, from her relatives only. We also compute Houaga's accuracy,
#   cor(y, yhat), and split it into its genetic and environmental parts, which only a simulation can do.
#   Test cows are also split by their herd: "seen" (other cows of her herd are in training) and "unseen" (none are,
#   e.g. every single-cow herd).
#
# EXPECTED (written 7 November 2026, before any result)
#   L7  The masked cows' EBVs are far less accurate than the training cows' (guess 0.15-0.30 against 0.35-0.6);
#       accuracy rises with the number of cows.
#   L8  cor(y, yhat) rises by +0.2 to +0.4 from G to GHS, almost all of it environment. For unseen herds GH is no
#       better than G, and only the field helps (GS about equal to GHS).
#
# The folds depend only on the replicate (seed 20261101 + replicate), so steps 1-3 of the ladder hide the same cows.
# The class, forward and region options below belong to other studies of the project.
# =====================================================================================================================

# Stage Q-VALIDATE — Houaga's validation (mask -> refit -> predict) on one simulated replicate, with every
# prediction split into its parts. Only a simulation can do the split, because it needs the true parts. Study S2.
#
#   Rscript scripts/questions/q_validate.R --run runs_q/S2_validation/g07_pe05_class_r01 --cv class --qmodels G,GH,GS,GHS
#   --cv class    mask each exotic-genome class (dgrp 1-4) in turn: Houaga Table 5, "by genome class"
#   --cv forward  mask the young cows: first lactation, first record in calving-year index >= 5 (about 110 cows).
#                 The backbone has no birth years; this stands in for Houaga's cows born 2016-17 (146 cows).
#                 Reported overall, by class and by region: Houaga Table 6
#   --cv region   mask each region in turn
#   --cv kfold    k-fold cross-validation (--folds 5): the cows are split at random into k folds; each fold (all records
#                 of its cows) is masked once and predicted from the other k - 1, so EVERY cow is predicted once.
#                 The split depends only on the replicate number, so the steps of study SB (which share their cows)
#                 use the same folds. Rows per fold, plus POOLED rows (mask "pooled") over all cows: groups
#                 "all", "herd_seen" (the masked cow's herd still has training records) and "herd_unseen" (it has
#                 none, e.g. a single-cow herd)
#   options: --test-region R4 (one region only, for smoke tests), --samples 200 (posterior samples per prediction)
#            --priors pc|isidore (as q_fit.R; default pc)
#
# For every model x masked set x group it writes output/q_validate_<cv>.csv with
#   cor_y_yhat     Houaga's accuracy: cor(y, yhat), yhat = fixed + animal (+ perm) (+ herd) (+ field)
#   cor_y_gen      cor(y, fixed + animal): the genetic part of the prediction (plus fixed effects)
#   cor_y_env      cor(y, herd + field): the environmental part
#   cor_y_fixed    cor(y, fixed): fixed effects alone (reference)
#   cor_y_animal   cor(y, animal): EBVs alone
#   r_ebv_a        cor(EBV, true a0) over the masked cows (cow level)
#   r_shat_s       cor(s_hat, true s) over the masked cows' herds (herd level)
#   ceil_gen       cor(y, true Xb + true a0): best possible genetic + fixed prediction
#   ceil_all       cor(y, y - e): best possible prediction of any kind
#   r_ebv_a_train  cor(EBV, true a0) over the TRAINING cows of the same fit (cows with own records), for contrast
# Same priors and integration as q_fit.R (PC priors, Houaga's PC-Matern, empirical Bayes); config = TRUE for predict().

source("scripts/sim/_fit_setup.R")    # dat, herds, animals, Ginv, mesh, matern, fixed_formula (reads only data/)
suppressPackageStartupMessages(library(Matrix))
args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(flag, default = NULL) if (flag %in% args) args[match(flag, args) + 1] else default
cv <- arg_value("--cv", "class")
if (!cv %in% c("class", "forward", "region", "kfold")) stop("--cv class|forward|region|kfold")
n_folds <- as.integer(arg_value("--folds", "5"))
priors <- arg_value("--priors", "pc")
if (!priors %in% c("pc", "isidore")) stop("--priors must be pc or isidore")
qmodels <- strsplit(arg_value("--qmodels", "G,GH,GS,GHS"), ",")[[1]]
n_samp <- as.integer(arg_value("--samples", "200"))
test_region <- arg_value("--test-region")
inla.setOption(num.threads = Sys.getenv("INLA_THREADS", "4:1"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(out_dir, paste0("q_validate_", cv, "_log.txt"))
cat("", file = log_file)
say <- function(...) { txt <- paste0(...); cat(txt, "\n"); cat(txt, "\n", file = log_file, append = TRUE) }

# ---- truth (read only to score the predictions, after fitting) ---------------------------------
f_tr <- file.path(truth_dir, "true_records.csv")
if (!file.exists(f_tr)) stop("Missing ", f_tr, ": simulate with q_simulate.R of 4 Nov 2026 or later")
tr <- read.csv(f_tr, stringsAsFactors = FALSE)
ta <- read.csv(file.path(truth_dir, "true_animals.csv"), stringsAsFactors = FALSE)
th <- read.csv(file.path(truth_dir, "true_herds.csv"), stringsAsFactors = FALSE)

# ---- optional: one region only (smoke tests), as in q_fit.R -------------------------------------
if (!is.null(test_region)) {
  G_full <- readRDS(file.path(data_dir, "G.rds"))
  herds <- herds[herds$region == test_region, ]
  dat <- dat[dat$herd %in% herds$herd, ]
  keep_ids <- c(unique(dat$animal_id), animals$animal_id[animals$role != "cow"])
  G_full <- G_full[keep_ids, keep_ids]
  Ginv <- solve(G_full); Ginv <- (Ginv + t(Ginv)) / 2; dimnames(Ginv) <- dimnames(G_full)
  n_animal <- nrow(Ginv); n_herd <- nrow(herds)
  dat$cowI <- match(dat$animal_id, rownames(Ginv)); dat$herdI <- match(dat$herd, herds$herd)
  mesh <- build_mesh(as.matrix(herds[, c("x_km", "y_km")]), mesh_type_of(data_dir))
  matern <- inla.spde2.pcmatern(mesh, prior.sigma = c(sqrt(0.25), 0.5), prior.range = c(50, 0.8))
  for (v in c("tyrmn", "cyrsn", "lacgr")) dat[[v]] <- droplevels(dat[[v]])
  dat$region <- droplevels(dat$region)
  say("TEST MODE: region ", test_region, " only (", nrow(herds), " herds, ", nrow(dat), " records)")
}

# ---- models: same terms and priors as q_fit.R ---------------------------------------------------
pc <- list(prec = list(prior = "pc.prec", param = c(0.5, 0.5)))
family_control <- list(hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.5))))
if (priors == "isidore") { pc <- NULL; family_control <- list() }   # INLA default loggamma(1, 5e-5), as Isidore
fx <- paste0("fixed_effects(main = ", paste(deparse(fixed_formula), collapse = ""), ", model = \"fixed\")")
T_animal <- "animal(main = cowI, model = \"generic\", Cmatrix = Ginv, mapper = bru_mapper_index(n_animal), hyper = pc)"
T_perm <- "perm(main = cowI, model = \"iid\", mapper = bru_mapper_index(n_animal), hyper = pc)"
T_herd <- "herd(main = herdI, model = \"iid\", mapper = bru_mapper_index(n_herd), hyper = pc)"
T_field <- "field(main = cbind(x_km, y_km), model = matern)"
lib <- list(G = c(fx, T_animal), GH = c(fx, T_animal, T_herd), GS = c(fx, T_animal, T_field),
            GHS = c(fx, T_animal, T_herd, T_field), GPHS = c(fx, T_animal, T_perm, T_herd, T_field),
            GP = c(fx, T_animal, T_perm), GPH = c(fx, T_animal, T_perm, T_herd), GPS = c(fx, T_animal, T_perm, T_field))
parts_of <- list(G = c("fixed_effects", "animal"), GH = c("fixed_effects", "animal", "herd"),
                 GS = c("fixed_effects", "animal", "field"), GHS = c("fixed_effects", "animal", "herd", "field"),
                 GPHS = c("fixed_effects", "animal", "perm", "herd", "field"), GP = c("fixed_effects", "animal", "perm"),
                 GPH = c("fixed_effects", "animal", "perm", "herd"), GPS = c("fixed_effects", "animal", "perm", "field"))
bad <- setdiff(qmodels, names(lib)); if (length(bad)) stop("q_validate.R has no model ", paste(bad, collapse = ", "))
fit_masked <- function(m, data) {
  cmpf <- as.formula(paste("~ -1 +", paste(lib[[m]], collapse = " + ")), env = globalenv())
  bru(cmpf,
      bru_obs(formula = y ~ ., family = "gaussian", data = data, control.family = family_control),
      options = list(control.inla = list(int.strategy = "eb"),
                     control.compute = list(config = TRUE), verbose = FALSE))
}

# ---- what is masked -------------------------------------------------------------------------------
dat$year_n <- as.integer(dat$year)
dat$lac_n <- as.integer(as.character(dat$lacgr))
if (cv == "class") {
  masks <- lapply(setNames(levels(droplevels(dat$dgrp)), paste0("class", levels(droplevels(dat$dgrp)))), function(k) dat$dgrp == k)
} else if (cv == "region") {
  masks <- lapply(setNames(levels(dat$region), levels(dat$region)), function(k) dat$region == k)
} else if (cv == "kfold") {
  rep_no <- suppressWarnings(as.integer(regmatches(run_dir, regexpr("[0-9]+$", run_dir))))
  if (!length(rep_no) || is.na(rep_no)) rep_no <- 0L
  cows_all <- sort(unique(dat$animal_id))
  set.seed(20261101L + rep_no)
  fold_of <- setNames(sample(rep_len(seq_len(n_folds), length(cows_all))), cows_all)
  masks <- lapply(setNames(seq_len(n_folds), paste0("fold", seq_len(n_folds))), function(k) dat$animal_id %in% cows_all[fold_of == k])
  say(sprintf("%d-fold cross-validation: %d cows, %s cows per fold; every cow is predicted once", n_folds, length(cows_all),
              paste(as.numeric(table(fold_of)), collapse = "/")))
} else {
  first_year <- tapply(dat$year_n, dat$animal_id, min)
  first_lac <- tapply(dat$lac_n, dat$animal_id, min)
  young <- names(first_year)[first_year >= 5 & first_lac == 1]
  masks <- list(young = dat$animal_id %in% young)
  say(sprintf("forward validation: %d young cows (first lactation, first record in year index >= 5), %d records",
              length(young), sum(masks$young)))
}

# ---- fit, predict, split ----------------------------------------------------------------------------
dat$rec_id <- as.character(dat$rec_id)
tr <- tr[match(dat$rec_id, tr$rec_id), ]
stopifnot(!anyNA(tr$rec_id))
B_all <- fm_basis(mesh, as.matrix(herds[, c("x_km", "y_km")])); if (is.list(B_all)) B_all <- B_all$A
out_csv <- file.path(out_dir, paste0("q_validate_", cv, ".csv"))
res <- list()
pool_rec <- list(); pool_cow <- list(); pool_fold <- list()          # k-fold: predictions of every fold
cr <- function(a, b) if (length(b) > 2 && sd(b) > 0) cor(a, b) else NA_real_
say("=== Q-VALIDATE ", run_dir, " | cv ", cv, " | models ", paste(qmodels, collapse = ","), " | masks ", paste(names(masks), collapse = ","),
    " | priors ", priors)
for (mk in names(masks)) {
  idx_test <- which(masks[[mk]])
  test <- dat[idx_test, ]
  tt <- tr[idx_test, ]
  dat_na <- dat; dat_na$y[idx_test] <- NA
  for (m in qmodels) {
    t0 <- Sys.time()
    fit <- tryCatch(fit_masked(m, dat_na), error = function(e) { say("   FIT FAILED ", m, " / ", mk, ": ", conditionMessage(e)); NULL })
    if (is.null(fit) || is.null(fit$summary.random$animal)) next
    pr <- sapply(parts_of[[m]], function(p)
      predict(fit, newdata = test, formula = as.formula(paste("~", p)), n.samples = n_samp, seed = 1L)$mean)
    pr <- matrix(pr, nrow = nrow(test), dimnames = list(NULL, parts_of[[m]]))
    part <- function(p) if (p %in% colnames(pr)) pr[, p] else rep(0, nrow(test))
    yhat <- rowSums(pr)
    gen <- part("fixed_effects") + part("animal")
    env <- part("herd") + part("field")
    has_env <- any(c("herd", "field") %in% colnames(pr))
    ebv_all <- setNames(fit$summary.random$animal$mean, rownames(Ginv))
    s_hat_herd <- if (!is.null(fit$summary.random$field)) setNames(as.numeric(B_all %*% fit$summary.random$field$mean), herds$herd) else NULL
    groups <- list(all = seq_len(nrow(test)))
    if (cv == "forward") {
      for (k in sort(unique(as.character(test$dgrp)))) groups[[paste0("class", k)]] <- which(as.character(test$dgrp) == k)
      for (k in sort(unique(as.character(test$region)))) groups[[k]] <- which(as.character(test$region) == k)
    }
    seen <- test$herd %in% dat_na$herd[!is.na(dat_na$y)]   # does the masked cow's herd still have training records?
    if (cv == "kfold") { groups$herd_seen <- which(seen); groups$herd_unseen <- which(!seen) }
    train_cows <- setdiff(unique(dat_na$animal_id[!is.na(dat_na$y)]), unique(test$animal_id))
    r_train <- cor(ebv_all[train_cows], ta$a0_true[match(train_cows, ta$animal_id)])
    mins <- as.numeric(Sys.time() - t0, units = "mins")
    if (cv == "kfold") {                                      # keep every prediction for the pooled rows
      mcows <- unique(test$animal_id)
      pool_rec[[length(pool_rec) + 1]] <- data.frame(model = m, idx = idx_test, yhat = yhat, gen = gen, env = env,
                                                     fixed = part("fixed_effects"), animal = part("animal"), seen = seen)
      pool_cow[[length(pool_cow) + 1]] <- data.frame(model = m, cow = mcows, ebv = unname(ebv_all[mcows]))
      pool_fold[[length(pool_fold) + 1]] <- data.frame(model = m, r_train = r_train, minutes = mins, has_env = has_env)
    }
    first_row <- length(res) + 1                                    # the "all" group, always written
    for (g in names(groups)) {
      ii <- groups[[g]]
      if (length(ii) < 10) next
      cows <- unique(test$animal_id[ii]); hh <- unique(test$herd[ii])
      res[[length(res) + 1]] <- data.frame(
        cv = cv, mask = mk, group = g, model = m, n_records = length(ii), n_cows = length(cows), n_herds = length(hh),
        cor_y_yhat = cr(test$y[ii], yhat[ii]), cor_y_gen = cr(test$y[ii], gen[ii]),
        cor_y_env = if (has_env) cr(test$y[ii], env[ii]) else NA_real_,
        cor_y_fixed = cr(test$y[ii], part("fixed_effects")[ii]), cor_y_animal = cr(test$y[ii], part("animal")[ii]),
        r_ebv_a = if (length(cows) > 5) cor(ebv_all[cows], ta$a0_true[match(cows, ta$animal_id)]) else NA_real_,
        r_shat_s = if (!is.null(s_hat_herd) && length(hh) > 5) cor(s_hat_herd[hh], th$s_true[match(hh, th$herd)]) else NA_real_,
        ceil_gen = cr(test$y[ii], tt$Xb_true[ii] + tt$a0[ii]),
        ceil_all = cr(test$y[ii], test$y[ii] - tt$e[ii]),
        r_ebv_a_train = r_train,
        minutes = as.numeric(Sys.time() - t0, units = "mins"))
    }
    r_all <- res[[first_row]]
    say(sprintf("%-4s mask %-7s %5.1f min | cor(y,yhat) %.3f | gen %.3f | env %s | r(EBV,a) %.3f", m, mk,
                as.numeric(Sys.time() - t0, units = "mins"), r_all$cor_y_yhat, r_all$cor_y_gen,
                ifelse(is.na(r_all$cor_y_env), "  -  ", sprintf("%.3f", r_all$cor_y_env)), r_all$r_ebv_a))
    write.csv(do.call(rbind, res), out_csv, row.names = FALSE)      # saved after every fit
    rm(fit); gc(verbose = FALSE)
  }
}
# ---- k-fold: pooled accuracy over all cows (each predicted once, from the fit without her fold) ------------
if (cv == "kfold" && length(pool_rec)) {
  PR <- do.call(rbind, pool_rec); PC <- do.call(rbind, pool_cow); PF <- do.call(rbind, pool_fold)
  for (m in unique(PR$model)) {
    nf <- sum(PF$model == m)
    if (nf < n_folds) { say("   pooled: ", m, " has ", nf, " of ", n_folds, " folds; no pooled rows"); next }
    p <- PR[PR$model == m, ]
    y <- dat$y[p$idx]; tp <- tr[p$idx, ]
    ebv_m <- setNames(PC$ebv[PC$model == m], PC$cow[PC$model == m])
    env_on <- any(PF$has_env[PF$model == m])
    for (g in c("all", "herd_seen", "herd_unseen")) {
      ii <- switch(g, all = seq_len(nrow(p)), herd_seen = which(p$seen), herd_unseen = which(!p$seen))
      if (length(ii) < 10) next
      cows <- unique(dat$animal_id[p$idx[ii]])
      res[[length(res) + 1]] <- data.frame(
        cv = cv, mask = "pooled", group = g, model = m, n_records = length(ii), n_cows = length(cows),
        n_herds = length(unique(dat$herd[p$idx[ii]])),
        cor_y_yhat = cr(y[ii], p$yhat[ii]), cor_y_gen = cr(y[ii], p$gen[ii]),
        cor_y_env = if (env_on) cr(y[ii], p$env[ii]) else NA_real_,
        cor_y_fixed = cr(y[ii], p$fixed[ii]), cor_y_animal = cr(y[ii], p$animal[ii]),
        r_ebv_a = cor(ebv_m[cows], ta$a0_true[match(cows, ta$animal_id)]),
        r_shat_s = NA_real_,
        ceil_gen = cr(y[ii], tp$Xb_true[ii] + tp$a0[ii]), ceil_all = cr(y[ii], y[ii] - tp$e[ii]),
        r_ebv_a_train = mean(PF$r_train[PF$model == m]), minutes = sum(PF$minutes[PF$model == m]))
    }
    r_all <- res[[which(vapply(res, function(r) r$mask == "pooled" && r$group == "all" && r$model == m, logical(1)))]]
    say(sprintf("%-4s POOLED over %d folds | cor(y,yhat) %.3f | gen %.3f | r(EBV,a) masked %.3f vs training %.3f", m, nf,
                r_all$cor_y_yhat, r_all$cor_y_gen, r_all$r_ebv_a, r_all$r_ebv_a_train))
  }
  write.csv(do.call(rbind, res), out_csv, row.names = FALSE)
}
say("Q-VALIDATE finished: ", out_csv)
