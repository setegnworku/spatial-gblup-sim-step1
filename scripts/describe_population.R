# =====================================================================================================================
# DESCRIBE THE SIMULATED POPULATION OF STEP 1 (docs/POPULATION.md). Run from the repository root after all step-1 data
# sets exist in runs_q/SB1_one_record/ (with founder caches in runs/founders/):  Rscript scripts/describe_population.R
# Part A reads every data set: herds, bulls and daughters, half-sisters and where they live, distances, G, true field.
# Part B rebuilds the pedigree (founders -> generation 1 -> generation 2 -> cows) of replicates 1-3 with the seeds of
# sb_simulate.R, checks it against the saved truth (same bull for every cow, same true breeding values), and computes
# pedigree relationships. Output: summary tables only (no coordinates) in tmp/popdesc/; copies are in results/population/.
# Note: in A_sire_classes.csv a class with no bull has cows = NA; read it as 0 (docs/POPULATION.md does).
# =====================================================================================================================
# Describe the simulated population of study SB, step 1 (SIMULATED): families, genome, geography, placement in herds
# and in space. Run from ~/Spacial_Cl on ilri-hpc. Writes only summary tables (no coordinates) to tmp/popdesc/.
#   Part A: every saved step-1 data set (runs_q/SB1_one_record/*): herds, bulls, half-sisters, distances, G.
#   Part B: rebuild the pedigree (founders -> generation 1 -> generation 2 -> cows) of replicates 1-3 from the cached
#           founders with the seeds of sb_simulate.R, check it against the saved truth, then pedigree relationships.
suppressPackageStartupMessages(library(Matrix))
root <- "runs_q/SB1_one_record"
out <- "tmp/popdesc"; dir.create(out, recursive = TRUE, showWarnings = FALSE)
runs <- sort(list.dirs(root, recursive = FALSE))
cat("runs:", length(runs), "\n")
cbin <- function(d, same) {   # cow-pair distance classes
  b <- cut(d, c(-Inf, 10, 33, 100, Inf), labels = c("0-10 km", "10-33 km", "33-100 km", "100+ km"))
  b <- as.character(b); b[same] <- "same herd"; factor(b, levels = c("same herd", "0-10 km", "10-33 km", "33-100 km", "100+ km"))
}
hbreaks <- c(0, 5, 10, 20, 33, 66, 100, 200, Inf)
A_rows <- list(); A_sire <- list(); A_field <- list(); A_gdist <- list(); A_hsdist <- list()

for (r in runs) {
  cell <- basename(r)
  hd <- read.csv(file.path(r, "data/herds.csv"), stringsAsFactors = FALSE)
  th <- read.csv(file.path(r, "truth/true_herds.csv"), stringsAsFactors = FALSE)
  ta <- read.csv(file.path(r, "truth/true_animals.csv"), stringsAsFactors = FALSE)
  tb <- read.csv(file.path(r, "truth/true_bulls.csv"), stringsAsFactors = FALSE)
  ncow <- nrow(ta); nh <- nrow(hd)
  th <- th[match(hd$herd, th$herd), ]
  sz <- th$n_cows
  xy <- as.matrix(hd[, c("x_km", "y_km")])
  Dh <- as.matrix(dist(xy)); diag(Dh) <- Inf
  nn <- apply(Dh, 1, min)
  w10 <- rowSums(Dh <= 10); w33 <- rowSums(Dh <= 33)
  cows10 <- as.numeric((Dh <= 10) %*% sz); cows33 <- as.numeric((Dh <= 33) %*% sz)   # cows in OTHER herds nearby
  hi <- match(ta$herd, hd$herd)
  # field correlogram between herds (true field values)
  s <- th$s_true; s <- (s - mean(s)) / sd(s)
  ut <- which(upper.tri(Dh)); ij <- arrayInd(ut, dim(Dh))
  hb <- cut(Dh[ut], hbreaks, right = FALSE)
  pr <- s[ij[, 1]] * s[ij[, 2]]
  A_field[[cell]] <- data.frame(cell = cell, bin = levels(hb), pairs = as.integer(table(hb)),
                                corr = as.numeric(tapply(pr, hb, mean)))
  rm(pr, hb)
  # sires and daughters
  d <- tb$n_daughters; stopifnot(sum(d) == ncow)
  sire_i <- match(ta$sire, tb$id)
  hs_per_cow <- d[sire_i] - 1
  cls <- cut(d, c(-1, 0, 1, 4, 9, 19, 49, Inf), labels = c("0", "1", "2-4", "5-9", "10-19", "20-49", "50+"))
  A_sire[[cell]] <- data.frame(cell = cell, class = levels(cls), sires = as.integer(table(cls)),
                               cows = as.integer(tapply(d, cls, sum)))
  # half-sister pairs: same herd / distance class
  hsd <- list()
  for (k in which(d >= 2)) {
    cw <- which(sire_i == k); h <- hi[cw]
    pp <- combn(length(cw), 2)
    same <- h[pp[1, ]] == h[pp[2, ]]
    dd <- sqrt(rowSums((xy[h[pp[1, ]], , drop = FALSE] - xy[h[pp[2, ]], , drop = FALSE])^2))
    hsd[[length(hsd) + 1]] <- data.frame(a = cw[pp[1, ]], b = cw[pp[2, ]], same = same, d = dd)
  }
  hsd <- do.call(rbind, hsd)
  hsd$bin <- cbin(hsd$d, hsd$same)
  per_cow <- function(sel) tabulate(c(hsd$a[sel], hsd$b[sel]), ncow)
  hs_same <- per_cow(hsd$same); hs_10 <- per_cow(hsd$d <= 10 | hsd$same); hs_33 <- per_cow(hsd$d <= 33 | hsd$same)
  # all cow pairs by distance class, from herd-level counts
  wpair <- sz[ij[, 1]] * sz[ij[, 2]]
  allb <- cbin(Dh[ut], rep(FALSE, length(ut)))
  all_by <- tapply(wpair, allb, sum); all_by[is.na(all_by)] <- 0
  all_by["same herd"] <- sum(choose(sz, 2))
  hs_by <- table(hsd$bin)
  A_hsdist[[cell]] <- data.frame(cell = cell, bin = names(hs_by), hs_pairs = as.integer(hs_by),
                                 all_pairs = as.numeric(all_by[names(hs_by)]))
  # herd-mates that are half-sisters
  hm_pairs <- sum(choose(sz, 2))
  # G (raw VanRaden; G.rds = 0.99 G + 0.01 I), where it was kept
  gstat <- list(g_diag = NA, g_off_mean = NA, g_off_sd = NA, g_gt01 = NA, g_gt02 = NA, g_hs_mean = NA, g_hs_sd = NA,
                g_hs_min = NA, g_hs_max = NA, g_nonhs_mean = NA, g_nonhs_sd = NA)
  fG <- file.path(r, "data/G.rds")
  if (file.exists(fG)) {
    G <- readRDS(fG); G <- G[ta$animal_id, ta$animal_id]
    Graw_d <- (diag(G) - 0.01) / 0.99
    utc <- which(upper.tri(G)); go <- G[utc] / 0.99
    ijc <- arrayInd(utc, dim(G))
    ishs <- sire_i[ijc[, 1]] == sire_i[ijc[, 2]]
    dc <- sqrt(rowSums((xy[hi[ijc[, 1]], ] - xy[hi[ijc[, 2]], ])^2))
    gb <- cbin(dc, hi[ijc[, 1]] == hi[ijc[, 2]])
    A_gdist[[cell]] <- data.frame(cell = cell, bin = levels(gb), pairs = as.integer(table(gb)),
                                  g_mean = as.numeric(tapply(go, gb, mean)),
                                  hs_share = as.numeric(tapply(ishs, gb, mean)))
    gstat <- list(g_diag = mean(Graw_d), g_off_mean = mean(go), g_off_sd = sd(go), g_gt01 = mean(go > 0.1),
                  g_gt02 = mean(go > 0.2), g_hs_mean = mean(go[ishs]), g_hs_sd = sd(go[ishs]), g_hs_min = min(go[ishs]),
                  g_hs_max = max(go[ishs]), g_nonhs_mean = mean(go[!ishs]), g_nonhs_sd = sd(go[!ishs]))
    rm(G, go, ijc, ishs, dc, gb, utc); gc(verbose = FALSE)
  }
  log <- readLines(file.path(r, "output/00_sim_summary.txt"))
  gl <- grep("^-- G:", log, value = TRUE)
  nsnp <- as.integer(sub(".*SNPs ([0-9]+).*", "\\1", gl))
  A_rows[[cell]] <- data.frame(
    cell = cell, ncow = ncow, rep = as.integer(sub(".*_r", "", cell)), p1 = as.numeric(sub(".*_p([0-9]+)_.*", "\\1", cell)) / 100,
    n_herd = nh, herd_single = mean(sz == 1), cows_single = sum(sz[sz == 1]) / ncow,
    herd_2 = sum(sz == 2), herd_3 = sum(sz == 3), herd_4 = sum(sz == 4), herd_5p = sum(sz >= 5), herd_max = max(sz),
    n_region = length(unique(hd$region)), n_ward = length(unique(hd$ward)),
    ext_x = diff(range(xy[, 1])), ext_y = diff(range(xy[, 2])),
    nn_med = median(nn), nn_q10 = unname(quantile(nn, 0.1)), nn_q90 = unname(quantile(nn, 0.9)),
    herds10_med = median(w10), herds33_med = median(w33),
    cows10_mean = mean(cows10[hi] + sz[hi] - 1), cows33_mean = mean(cows33[hi] + sz[hi] - 1),
    n_sire = length(d), sires_used = sum(d > 0), dps_mean = mean(d[d > 0]), dps_med = median(d[d > 0]),
    dps_q90 = unname(quantile(d[d > 0], 0.9)), dps_max = max(d), top10_share = sum(sort(d, decreasing = TRUE)[1:10]) / ncow,
    hs_pairs = nrow(hsd), hs_per_cow_mean = mean(hs_per_cow), hs_per_cow_med = median(hs_per_cow), cows_no_hs = mean(hs_per_cow == 0),
    hs_same_herd_cows = mean(hs_same > 0), hs_10_cows = mean(hs_10 > 0), hs_33_cows = mean(hs_33 > 0),
    hs_10_mean = mean(hs_10), hs_33_mean = mean(hs_33),
    hm_pairs = hm_pairs, hm_hs_pairs = sum(hsd$same), hs_d_med = median(hsd$d[!hsd$same]),
    nsnp = nsnp, gstat)
  cat(cell, "done\n")
}
wcsv <- function(x, f) write.csv(do.call(rbind, x), file.path(out, f), row.names = FALSE)
wcsv(A_rows, "A_by_dataset.csv"); wcsv(A_sire, "A_sire_classes.csv"); wcsv(A_field, "A_field_corr.csv")
wcsv(A_gdist, "A_G_by_distance.csv"); wcsv(A_hsdist, "A_halfsis_by_distance.csv")
cat("PART A finished\n")

# ---- Part B: pedigree of replicates 1-3 ------------------------------------------------------------------------------
suppressPackageStartupMessages(library(AlphaSimR))
B_rows <- list(); B_cls <- list(); B_off <- list()
tally <- function(x, who, what) { t <- as.integer(table(x)); data.frame(who = who, what = what, mean = mean(t), sd = sd(t),
                                                                         min = min(t), max = max(t), n_parents = length(t)) }
for (rep_id in 1:3) {
  seed <- 20261001L + 1000L * rep_id
  founder <- readRDS(sprintf("runs/founders/founderSB_rep%02d.rds", rep_id))
  set.seed(seed + 261L)
  SP <- SimParam$new(founder)
  SP$nThreads <- 1L
  SP$setSexes("yes_sys")
  SP$restrSegSites(minQtlPerChr = 200, minSnpPerChr = 1500, overlap = FALSE)
  SP$addTraitA(nQtlPerChr = 200, mean = 0, var = 1)
  SP$addSnpChip(nSnpPerChr = 1500)
  base <- newPop(founder, simParam = SP)
  P1 <- randCross(base, nCrosses = 3000, simParam = SP)
  P2 <- randCross(P1, nCrosses = 10000, simParam = SP)
  B_off[[length(B_off) + 1]] <- rbind(
    tally(P1@father, "founder sires", sprintf("rep%d: offspring in generation 1", rep_id)),
    tally(P1@mother, "founder dams", sprintf("rep%d: offspring in generation 1", rep_id)),
    tally(P2@father, "generation-1 sires", sprintf("rep%d: offspring in generation 2", rep_id)),
    tally(P2@mother, "generation-1 dams", sprintf("rep%d: offspring in generation 2", rep_id)))
  B_off[[length(B_off)]]$founders_unused <- c(sum(!base@id[base@sex == "M"] %in% P1@father), sum(!base@id[base@sex == "F"] %in% P1@mother),
                                              sum(!P1@id[P1@sex == "M"] %in% P2@father), sum(!P1@id[P1@sex == "F"] %in% P2@mother))
  # generation-1 relationships (founders unrelated in the pedigree)
  for (n_cow in c(1000, 2000, 4000)) {
    set.seed(seed + 262L)
    n_sire <- max(20L, as.integer(round(n_cow / 10)))
    males <- P2[P2@sex == "M"]; fems <- P2[P2@sex == "F"]
    sires <- males[sample.int(males@nInd, n_sire)]
    dams <- fems[sample.int(fems@nInd, n_cow)]
    sire_of <- sample.int(n_sire, n_cow, replace = TRUE, prob = rlnorm(n_sire, 0, 1))
    calves <- makeCross2(females = dams, males = sires, crossPlan = cbind(seq_len(n_cow), sire_of), simParam = SP)
    g_cow <- gv(calves)[, 1]; a_true <- sqrt(0.06 / var(g_cow)) * (g_cow - mean(g_cow))
    chk_dir <- sprintf("%s/n%d_p70_r%02d", root, n_cow, rep_id)
    ta <- read.csv(file.path(chk_dir, "truth/true_animals.csv"), stringsAsFactors = FALSE)
    same_sire <- all(ta$sire == sprintf("S%04d", sire_of)); max_da <- max(abs(ta$a0_true - a_true))
    # pedigree A: generation 1 (needed animals) -> generation 2 (sires and dams) -> cows
    S2_f <- c(sires@father, dams@father); S2_m <- c(sires@mother, dams@mother)
    U1 <- unique(c(S2_f, S2_m))
    p1i <- match(U1, P1@id); f1 <- P1@father[p1i]; m1 <- P1@mother[p1i]
    fo <- unique(c(f1, m1))
    K1 <- sparseMatrix(i = c(seq_along(U1), seq_along(U1)), j = match(c(f1, m1), fo), x = 1, dims = c(length(U1), length(fo)))
    A1 <- as.matrix(0.25 * tcrossprod(K1)); diag(A1) <- 1
    n2 <- n_sire + n_cow
    K2 <- sparseMatrix(i = c(1:n2, 1:n2), j = match(c(S2_f, S2_m), U1), x = 1, dims = c(n2, length(U1)))
    A2 <- as.matrix(0.25 * K2 %*% A1 %*% t(K2))
    f2 <- match(S2_f, U1); m2 <- match(S2_m, U1)
    diag(A2) <- 1 + 0.5 * A1[cbind(f2, m2)]
    K3 <- sparseMatrix(i = c(1:n_cow, 1:n_cow), j = c(sire_of, n_sire + seq_len(n_cow)), x = 1, dims = c(n_cow, n2))
    A3 <- as.matrix(0.25 * K3 %*% A2 %*% t(K3))
    Fcow <- 0.5 * A2[cbind(sire_of, n_sire + seq_len(n_cow))]
    As <- A2[1:n_sire, 1:n_sire]; Ad <- A2[n_sire + 1:n_cow, n_sire + 1:n_cow]; Asd <- A2[1:n_sire, n_sire + 1:n_cow]
    ut_s <- As[upper.tri(As)]; ut_d <- Ad[upper.tri(Ad)]
    sire_pairs <- c(full = sum(ut_s >= 0.5 - 1e-9), half = sum(ut_s >= 0.25 - 1e-9 & ut_s < 0.5 - 1e-9), cousin = sum(ut_s > 0 & ut_s < 0.25 - 1e-9))
    dam_pairs <- c(full = sum(ut_d >= 0.5 - 1e-9), half = sum(ut_d >= 0.25 - 1e-9 & ut_d < 0.5 - 1e-9), cousin = sum(ut_d > 0 & ut_d < 0.25 - 1e-9))
    utc <- which(upper.tri(A3)); a3 <- A3[utc]
    cl <- cut(a3, c(-1, 1e-12, 0.015625 - 1e-9, 0.0625 - 1e-9, 0.25 - 1e-9, Inf),
              labels = c("unrelated (0)", "distant (< 1/64)", "second-cousin level (1/64 to < 1/16)",
                         "cousin level (1/16 to < 1/4)", "half-sisters or closer (>= 1/4)"))
    fG <- file.path(chk_dir, "data/G.rds"); gm <- rep(NA, nlevels(cl)); gs <- gm
    if (file.exists(fG)) {
      G <- readRDS(fG)[ta$animal_id, ta$animal_id]; go <- G[utc] / 0.99
      gm <- as.numeric(tapply(go, cl, mean)); gs <- as.numeric(tapply(go, cl, sd))
      cor_ga <- cor(go, a3); rm(G, go)
    } else cor_ga <- NA
    B_cls[[length(B_cls) + 1]] <- data.frame(rep = rep_id, ncow = n_cow, class = levels(cl), pairs = as.numeric(table(cl)),
                                             per_cow = 2 * as.numeric(table(cl)) / n_cow, g_mean = gm, g_sd = gs)
    sire_par <- unique(sires@father); dam_par <- unique(dams@father)
    B_rows[[length(B_rows) + 1]] <- data.frame(
      rep = rep_id, ncow = n_cow, check_same_sires = same_sire, check_max_diff_a = max_da,
      founders = base@nInd, founders_m = sum(base@sex == "M"), gen1 = P1@nInd, gen1_m = sum(P1@sex == "M"),
      gen2 = P2@nInd, gen2_m = sum(P2@sex == "M"), n_sire = n_sire, n_dam = n_cow,
      sires_distinct_fathers = length(unique(sires@father)), sires_distinct_mothers = length(unique(sires@mother)),
      dams_distinct_fathers = length(unique(dams@father)), dams_distinct_mothers = length(unique(dams@mother)),
      sire_pairs_full = sire_pairs[["full"]], sire_pairs_half = sire_pairs[["half"]], sire_pairs_cousin = sire_pairs[["cousin"]],
      sire_pairs_total = length(ut_s),
      dam_pairs_full = dam_pairs[["full"]], dam_pairs_half = dam_pairs[["half"]], dam_pairs_cousin = dam_pairs[["cousin"]],
      dam_pairs_total = length(ut_d),
      sire_dam_related_share = mean(Asd > 0), cow_F_mean = mean(Fcow), cow_F_share_pos = mean(Fcow > 0), cow_F_max = max(Fcow),
      dams_repeated = sum(duplicated(dams@id)),
      cor_G_A = cor_ga)
    cat(sprintf("rep %d N %d: sires reproduced %s, max |a diff| %.2e\n", rep_id, n_cow, same_sire, max_da))
    rm(A1, A2, A3, a3, cl, K1, K2, K3); gc(verbose = FALSE)
  }
  rm(base, P1, P2, SP, founder); gc(verbose = FALSE)
}
wcsv(B_rows, "B_pedigree.csv"); wcsv(B_cls, "B_relationship_classes.csv"); wcsv(B_off, "B_offspring_per_parent.csv")
cat("PART B finished\n")
