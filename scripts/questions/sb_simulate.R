# =====================================================================================================================
# STUDY SB, STEP 1: ONE RECORD PER COW  (this repository holds step 1 only; README.md has the full story and results)
# ---------------------------------------------------------------------------------------------------------------------
# WHY WE SIMULATE
#   Houaga et al. (2026, Genetics Selection Evolution) fitted genomic models with a herd effect and an SPDE spatial
#   field (G, GH, GS, GHS) to test-day milk records of smallholder dairy cows in Tanzania. In real data the true
#   genetic variance and the true breeding values are unknown, so nobody can check whether the estimates are right.
#   A simulation puts KNOWN values into data built like the real data (the same kind of herd locations, herd sizes
#   and recording pattern), fits the same models, and compares every estimate with the truth.
#   Questions for step 1:
#     1. Can the models tell genetics, herd and location apart when most herds have ONE cow?
#     2. How many cows are needed (1,000, 2,000 or 4,000)?
#     3. Does the share of single-cow herds (50, 60, 70, 80 %) matter?
#     4. Where does a term go when a model leaves it out?
#     5. Does the spatial field improve the breeding values, also of cows with no records?
#     6. Does the prior on the variances matter (INLA's default prior against a PC prior)?
#
# WHY A LADDER, AND WHY STEP 1 FIRST
#   The real problem mixes several things at once: repeated records, a cow permanent environment (pe), crossbreds.
#   The ladder adds ONE thing per step on the SAME simulated cows (same seeds, same founders, same field):
#     step 1  one record per cow, no pe            <- this repository
#     step 2  about ten records per cow, no pe
#     step 3  about ten records per cow, true pe
#   Step 1 is the baseline. With one record per cow nothing repeats, so only four kinds of sharing are in play:
#   genetics (shared by relatives), herd (shared by herd-mates), field (shared by neighbours), residual (shared by nobody).
#
# STEP 1 SETTINGS (passed by cluster/tasks_SB1*.txt):  --records one --spe 0
#   --records one  each cow keeps ONE random test-day of the record pattern she was given (same cows, herds and field
#                  as the later steps);
#   --spe 0        no permanent environment (pe_true is all 0).
#   Truth: sigma2_g 0.06, sigma2_h 0.15, sigma2_s 0.41 (SPDE parameter; about 0.35 across the herds), range 33 km,
#          sigma2_e 0.29  ->  h2 = 0.06 / 0.91 = 0.066.
#   Design: 1,000 and 2,000 cows x 50 / 60 / 70 / 80 % single-cow herds, 10 replicates each (80 data sets), plus
#           4,000 cows at 70 % (5 replicates): 85 data sets. Replicate r uses the same founders and field in every cell.
#
# EXPECTED RESULTS (written 7 November 2026, before any step-1 result)
#   L1  GHS, the true model, estimates sigma2_g without bias (mean within +-15 % of 0.06), but imprecisely.
#   L2  Herd and residual trade off, because only multi-cow herds tell them apart; more so at 80 % single-cow herds.
#   L3  A model without a term misplaces it: G and GS put herd and field into sigma2_g; GH puts the field into the herd.
#   L4  Accuracy of the recorded cows' EBVs, r(EBV, a), of GHS is about 0.35-0.45.
#   (L7 and L8, for cross-validation, are in q_validate.R.)
#
# INPUT NOT IN THIS REPOSITORY
#   data/backbone/backbone_tz.rds is built from the real field data (herd coordinates, record patterns, fixed-effect
#   sizes) and is not shared. README.md ("Data not included") describes its structure. Comments below that cite
#   docs/... refer to project notes that are not in this repository; README.md holds their step-1 content.
# =====================================================================================================================

# Stage SB-SIM — ONE breed, a chosen number of cows and share of single-cow herds, repeated records with a
# permanent environment. Truth = GPHS:  y = Xb + a + pe + herd + field + e.
# Pre-registered note: docs/questions/SB_single_breed.md.
#
#   Rscript scripts/questions/sb_simulate.R --run runs_q/SB_single_breed/n2000_p70_r01 --rep 1 --ncow 2000 --p1 0.7
#   options: --sg 0.06 --spe 0.04 --sh 0.15 --ss 0.41 --range 33 --se 0.29   (true variances, range in km)
#            --records all|one   all = the template cow's full record pattern (default); one = ONE random test-day
#                                of it (step 1 of the ladder; the same cows, herds and field as with all records)
#            --dps 10  daughters per sire on average: n_sire = ncow / dps (default 10 = every earlier SB run, unchanged;
#                      P9 relatives study, 13 Nov 2026: 5, 2, 1)
#            --quick   (small genome, for local tests)
#
# What comes from the real AADGG backbone (data/backbone/backbone_tz.rds, anonymised; no milk record is used):
#   * herd LOCATIONS: resampled real herd locations (+ 1 km jitter); herds beyond the 1,386 real ones are extra
#     farms in the same villages;
#   * sizes of the multi-cow herds: drawn from the real multi-cow herds (2: 299, 3: 71, 4: 13, 5: 5, 9: 1);
#   * RECORD PATTERN of each cow: copied from a random real cow (number of test-days, calendar, lactation group,
#     age, Legendre covariates), and the real fixed-effect sizes WITHOUT the breed-class part.
# The share of single-cow herds (--p1) and the number of cows (--ncow) are the design factors.
# Writes the same data/ and truth/ files as q_simulate.R, so q_fit.R and q_extract.R run unchanged.
# SIMULATED PHENOTYPES — not AADGG milk records.

source("scripts/sim/_paths.R")      # data_dir, truth_dir, out_dir from --run
source("scripts/sim/_mesh.R")
if (is.null(run_dir)) stop("Give --run runs_q/SB_single_breed/<cell>_rNN")
suppressPackageStartupMessages({ library(AlphaSimR); library(Matrix); library(fmesher); library(INLA) })
inla.setOption(num.threads = "1:1")

arg_value <- function(flag, default = NULL) if (flag %in% cli_args) cli_args[match(flag, cli_args) + 1] else default
rep_id <- as.integer(arg_value("--rep", "1"))
n_cow <- as.integer(arg_value("--ncow", "2000"))
p1 <- as.numeric(arg_value("--p1", "0.7"))
quick <- "--quick" %in% cli_args
records_mode <- arg_value("--records", "all")
if (!records_mode %in% c("all", "one")) stop("--records must be all or one")
sigma2_g <- as.numeric(arg_value("--sg", "0.06"))
sigma2_pe <- as.numeric(arg_value("--spe", "0.04"))
sigma2_h <- as.numeric(arg_value("--sh", "0.15"))
sigma2_s <- as.numeric(arg_value("--ss", "0.41"))
range_km <- as.numeric(arg_value("--range", "33"))
sigma2_e <- as.numeric(arg_value("--se", "0.29"))
dps <- as.numeric(arg_value("--dps", "10"))
if (!(p1 > 0 && p1 < 1)) stop("--p1 must be between 0 and 1")
if (!(dps >= 1)) stop("--dps must be at least 1")
seed <- 20261001L + 1000L * rep_id            # as every other study: replicate r = the same field geography
nthr <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))

for (d in c(data_dir, truth_dir, out_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(out_dir, "00_sim_summary.txt")
cat("", file = log_file)
say <- function(...) { txt <- paste0(...); cat(txt, "\n"); cat(txt, "\n", file = log_file, append = TRUE) }
zscore <- function(x) (x - mean(x)) / sd(x)
say(sprintf("=== SB-SIM | rep %d | seed %d | cows %d | single-cow herds %.0f %% | records %s | sg %g spe %g sh %g ss %g range %g se %g%s ===",
            rep_id, seed, n_cow, 100 * p1, records_mode, sigma2_g, sigma2_pe, sigma2_h, sigma2_s, range_km, sigma2_e,
            if (quick) " | QUICK" else ""))

bb <- readRDS("data/backbone/backbone_tz.rds")
real_h <- bb$herds

# ---- 1. herd structure: exactly n_cow cows, share p1 of the herds with one cow ---------------------
tail_sizes <- as.integer(table(bb$cows$herd))
tail_sizes <- tail_sizes[tail_sizes >= 2]                     # real multi-cow herds (mean 2.31)
set.seed(seed + 202L)
n_multi <- round(n_cow * (1 - p1) / (p1 + (1 - p1) * mean(tail_sizes)))
size_multi <- sample(tail_sizes, n_multi, replace = TRUE)
n_single <- n_cow - sum(size_multi)
if (n_single < 0) stop("too few cows for this herd structure")
sizes <- sample(c(rep(1L, n_single), size_multi))
n_herd <- length(sizes)

# ---- 2. herd locations: resampled real ones -------------------------------------------------------
set.seed(seed + 203L)
nr <- nrow(real_h)
idx <- if (n_herd <= nr) sample.int(nr, n_herd) else c(sample.int(nr), sample.int(nr, n_herd - nr, replace = TRUE))
herds <- data.frame(herd = sprintf("H%05d", seq_len(n_herd)), region = real_h$region[idx], ward = real_h$ward[idx],
                    x_km = round(real_h$x_km[idx] + rnorm(n_herd, 0, 1), 4),
                    y_km = round(real_h$y_km[idx] + rnorm(n_herd, 0, 1), 4), n_cows = sizes)
say(sprintf("-- herds %d | single-cow herds %d (%.1f %% of herds, %.1f %% of cows) | herd sizes 2/3/4/5+: %s | real locations reused: %s",
            n_herd, sum(sizes == 1), 100 * mean(sizes == 1), 100 * sum(sizes == 1) / n_cow,
            paste(c(sum(sizes == 2), sum(sizes == 3), sum(sizes == 4), sum(sizes >= 5)), collapse = "/"),
            if (n_herd <= nr) "no repeats" else sprintf("all %d + %d repeats", nr, n_herd - nr)))

# ---- 3. spatial field: on the mesh of the REAL herd locations, the same in every cell of a replicate --------
# (same mesh and seed as all other question studies; the analyst's mesh is rebuilt from the simulated herds by
# _fit_setup.R, so it is close to, but not exactly, the mesh the truth was drawn on)
nu <- 1; kappa <- sqrt(8 * nu) / range_km; tau <- 1 / (sqrt(4 * pi) * kappa * sqrt(sigma2_s))
mesh_truth <- build_mesh(as.matrix(real_h[, c("x_km", "y_km")]), "backbone")
fem <- fm_fem(mesh_truth, order = 2)
Q_w <- as(forceSymmetric(tau^2 * (kappa^4 * fem$c0 + 2 * kappa^2 * fem$g1 + fem$g2)), "CsparseMatrix")
w_true <- as.numeric(inla.qsample(n = 1, Q = Q_w, seed = seed + 4L, num.threads = "1:1"))
B <- fm_basis(mesh_truth, as.matrix(herds[, c("x_km", "y_km")])); if (is.list(B)) B <- B$A
s_herd <- setNames(as.numeric(B %*% w_true), herds$herd)
mesh_fit <- build_mesh(as.matrix(herds[, c("x_km", "y_km")]), "backbone")
set.seed(seed + 210L)
h_true <- setNames(rnorm(n_herd, 0, sqrt(sigma2_h)), herds$herd)
say(sprintf("-- field: truth mesh %d nodes (real herds), analyst mesh %d nodes | realised var(s) across herds %.3f | var(h) %.3f",
            mesh_truth$n, mesh_fit$n, var(s_herd), var(h_true)))

# ---- 4. one breed: founders, population, sires, dams, cows -----------------------------------------
gen <- if (quick) list(nChr = 5, segSites = 700, nQtl = 100, nSnp = 500) else list(nChr = 29, segSites = 2000, nQtl = 200, nSnp = 1500)
founder_cache <- sprintf("runs/founders/founderSB%s_rep%02d.rds", if (quick) "quick" else "", rep_id)
if (file.exists(founder_cache)) {
  founder <- readRDS(founder_cache); founder_src <- founder_cache
} else {
  t0 <- Sys.time()
  set.seed(seed + 206L)
  founder <- runMacs(nInd = 1000, nChr = gen$nChr, segSites = gen$segSites, species = "CATTLE", nThreads = nthr)
  dir.create(dirname(founder_cache), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(founder_cache, ".tmp", Sys.getpid())          # atomic write: tasks of other cells may read it
  saveRDS(founder, tmp); file.rename(tmp, founder_cache)
  founder_src <- sprintf("%s (new, %.1f min)", founder_cache, as.numeric(Sys.time() - t0, units = "mins"))
}
say("-- founders: ", founder_src)
set.seed(seed + 261L)
SP <- SimParam$new(founder)
SP$nThreads <- nthr
SP$setSexes("yes_sys")
SP$restrSegSites(minQtlPerChr = gen$nQtl, minSnpPerChr = gen$nSnp, overlap = FALSE)
SP$addTraitA(nQtlPerChr = gen$nQtl, mean = 0, var = 1)
SP$addSnpChip(nSnpPerChr = gen$nSnp)
base <- newPop(founder, simParam = SP)
P1 <- randCross(base, nCrosses = 3000, simParam = SP)
P2 <- randCross(P1, nCrosses = 10000, simParam = SP)         # the breed today (same in every cell of a replicate)

set.seed(seed + 262L)
n_sire <- max(20L, as.integer(round(n_cow / dps)))       # dps = 10 (default) reproduces every earlier SB run
males <- P2[P2@sex == "M"]; fems <- P2[P2@sex == "F"]
if (fems@nInd < n_cow) stop("not enough females in the population for ", n_cow, " cows")
sires <- males[sample.int(males@nInd, n_sire)]
dams <- fems[sample.int(fems@nInd, n_cow)]
sire_of <- sample.int(n_sire, n_cow, replace = TRUE, prob = rlnorm(n_sire, 0, 1))   # uneven AI-bull use
calves <- makeCross2(females = dams, males = sires, crossPlan = cbind(seq_len(n_cow), sire_of), simParam = SP)
calves@sex[] <- "F"
cow_ids <- sprintf("C%05d", seq_len(n_cow))
dps <- tabulate(sire_of, n_sire)
say(sprintf("-- one breed: %d sires (not genotyped) | daughters per sire: median %g, mean %.1f, max %d | %d sires unused",
            n_sire, median(dps[dps > 0]), mean(dps[dps > 0]), max(dps), sum(dps == 0)))

# true breeding values from the QTL, scaled to sigma2_g (gv: fixed QTL scale, so sires are on the same scale)
g_cow <- gv(calves)[, 1]
k_g <- sqrt(sigma2_g / var(g_cow))
a_true <- setNames(k_g * (g_cow - mean(g_cow)), cow_ids)
a_sire <- k_g * (gv(sires)[, 1] - mean(g_cow))
set.seed(seed + 281L)
pe_true <- setNames(rnorm(n_cow) * sqrt(sigma2_pe), cow_ids)

# G of the cows (VanRaden), G* = 0.99 G + 0.01 I, as Isidore
M_snp <- pullSnpGeno(calves, simParam = SP)
p_snp <- colMeans(M_snp) / 2
keep <- pmin(p_snp, 1 - p_snp) >= 0.01
Zc <- sweep(M_snp[, keep], 2, 2 * p_snp[keep]); rm(M_snp)
G_raw <- tcrossprod(Zc) / (2 * sum(p_snp[keep] * (1 - p_snp[keep]))); rm(Zc)
dimnames(G_raw) <- list(cow_ids, cow_ids)
G <- 0.99 * G_raw + 0.01 * diag(n_cow); G <- (G + t(G)) / 2
Ginv <- solve(G); Ginv <- (Ginv + t(Ginv)) / 2; dimnames(Ginv) <- dimnames(G)
off <- G_raw[upper.tri(G_raw)]
say(sprintf("-- G: %d cows | SNPs %d (MAF >= 0.01) | QTL %d | diag mean %.3f | off-diagonal sd %.4f | max|G Ginv - I| %.1e",
            n_cow, sum(keep), gen$nChr * gen$nQtl, mean(diag(G_raw)), sd(off), max(abs(G %*% Ginv - diag(n_cow)))))

# ---- 5. cows to herds (at random: no link between genetics and environment) -------------------------
set.seed(seed + 263L)
cow_herd <- setNames(sample(rep(herds$herd, herds$n_cows)), cow_ids)

# ---- 6. records: the record pattern of a random real cow ---------------------------------------------
set.seed(seed + 204L)
rec_by_cow <- split(seq_len(nrow(bb$records)), bb$records$animal_id)
tmpl <- sample(names(rec_by_cow), n_cow, replace = TRUE)
rows <- rec_by_cow[tmpl]
if (records_mode == "one") {                                   # one random test-day of the template cow
  set.seed(seed + 205L)
  rows <- lapply(rows, function(r) r[sample.int(length(r), 1)])
}
nrec <- lengths(rows)
records <- bb$records[unlist(rows), c("tyrmn", "cyrsn", "lacgr", "ageZ", "leg1", "leg2", "year", "dgrp", "Xb_true")]
b_dg <- c(0, bb$b[c("dgrp2", "dgrp3", "dgrp4")]); b_dg[is.na(b_dg)] <- 0
records$Xb_true <- records$Xb_true - b_dg[records$dgrp]       # one breed: remove the breed-class part
records$dgrp <- NULL
records$animal_id <- rep(cow_ids, nrec)
records$herd <- unname(cow_herd[records$animal_id])
records$rec_id <- sprintf("R%06d", seq_len(nrow(records)))
records$ygrp <- pmin(pmax(records$year, 3L), 6L) - 2L
set.seed(seed + 211L)
e_true <- rnorm(nrow(records), 0, sqrt(sigma2_e))
recA <- records$animal_id; recH <- records$herd
records$y <- as.numeric(records$Xb_true + a_true[recA] + pe_true[recA] + h_true[recH] + s_herd[recH] + e_true)
say(sprintf("-- records %d | per cow: mean %.1f, median %g, range %d-%d | var: Xb %.3f, a %.3f, pe %.3f, y %.3f",
            nrow(records), mean(nrec), median(nrec), min(nrec), max(nrec), var(records$Xb_true),
            var(a_true), var(pe_true), var(records$y)))

# ---- 7. write data (analyst) and truth --------------------------------------------------------------
write.csv(herds[, c("herd", "region", "ward", "x_km", "y_km")], file.path(data_dir, "herds.csv"), row.names = FALSE)
write.csv(data.frame(animal_id = cow_ids, herd = unname(cow_herd), dgrp = 1L, sex = "F", role = "cow", phenotyped = TRUE,
                     sire_group = "none", exotic_class_prop = 1, sire = NA),
          file.path(data_dir, "animals.csv"), row.names = FALSE)
write.csv(records[, c("rec_id", "animal_id", "herd", "tyrmn", "cyrsn", "lacgr", "ageZ", "leg1", "leg2", "year", "ygrp", "y")],
          file.path(data_dir, "records.csv"), row.names = FALSE)
saveRDS(G, file.path(data_dir, "G.rds"))
saveRDS(Ginv, file.path(data_dir, "Ginv.rds"))
writeLines("backbone", file.path(data_dir, "mesh_type.txt"))

farm_E <- setNames(zscore(h_true + s_herd), herds$herd)
settings <- data.frame(setting = c("rep", "seed", "ncow", "p1", "p1_realised", "cows_single_share", "n_herd", "n_sire", "n_rec",
                                   "records", "quick", "lambda", "sigma2_g", "sigma2_pe", "env_shape"),
                       value = c(rep_id, seed, n_cow, p1, round(mean(sizes == 1), 4), round(sum(sizes == 1) / n_cow, 4), n_herd,
                                 n_sire, nrow(records), records_mode, quick, 0, sigma2_g, sigma2_pe, "spde"))
write.csv(settings, file.path(truth_dir, "sim_settings.csv"), row.names = FALSE)
write.csv(data.frame(parameter = c("sigma2_g", "sigma2_pe", "sigma2_h", "sigma2_ward", "sigma2_s", "range_km", "sigma2_e", "sigma2_slope",
                                   "gxe_breed", "realised_var_s_herds", "realised_var_s_within_region", "realised_var_a0_cows",
                                   "realised_var_pe_cows", "realised_var_h", "realised_var_y"),
                     value = c(sigma2_g, sigma2_pe, sigma2_h, 0, sigma2_s, range_km, sigma2_e, 0, 0, var(s_herd),
                               sum((s_herd - ave(s_herd, herds$region))^2) / (n_herd - length(unique(herds$region))),
                               var(a_true), var(pe_true), var(h_true), var(records$y))),
          file.path(truth_dir, "true_hyperparameters.csv"), row.names = FALSE)
write.csv(data.frame(animal_id = cow_ids, a0_true = unname(a_true), a1_true = 0, role = "cow", herd = unname(cow_herd),
                     pe_true = unname(pe_true), farm_E = unname(farm_E[cow_herd]), a_own_env = unname(a_true), exotic_prop = 1,
                     sire = sprintf("S%04d", sire_of), template_cow = tmpl, n_records = nrec),
          file.path(truth_dir, "true_animals.csv"), row.names = FALSE)
write.csv(data.frame(herd = herds$herd, region = herds$region, ward = herds$ward, n_cows = herds$n_cows, h_true = unname(h_true),
                     ward_true = 0, s_true = unname(s_herd), farm_E = unname(farm_E), breed_slope_true = 0),
          file.path(truth_dir, "true_herds.csv"), row.names = FALSE)
# sires are not genotyped and not in G (as in Isidore's data); type "sire" keeps them out of the AI-bull metrics
write.csv(data.frame(id = sprintf("S%04d", seq_len(n_sire)), type = "sire", a0_true = a_sire, a1_true = 0, n_daughters = dps),
          file.path(truth_dir, "true_bulls.csv"), row.names = FALSE)
write.csv(data.frame(rec_id = records$rec_id, Xb_true = records$Xb_true, a0 = unname(a_true[recA]), pe = unname(pe_true[recA]),
                     h = unname(h_true[recH]), ward = 0, s = unname(s_herd[recH]), gxe = 0, e = e_true),
          file.path(truth_dir, "true_records.csv"), row.names = FALSE)
saveRDS(list(mesh_truth = mesh_truth, w_true = w_true, env_shape = "spde", s_herd = s_herd), file.path(truth_dir, "field_truth.rds"))
asr_sha <- tryCatch(packageDescription("AlphaSimR")$RemoteSha, error = function(e) NULL)
writeLines(c(paste("seed", seed), paste("founders:", founder_src),
             paste("AlphaSimR", packageVersion("AlphaSimR"), if (!is.null(asr_sha)) substr(asr_sha, 1, 10) else "",
                   "| INLA", packageVersion("INLA"), "| fmesher", packageVersion("fmesher"))),
           file.path(truth_dir, "versions.txt"))
say("SB-SIM finished: ", run_dir)
