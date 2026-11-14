# =====================================================================================================================
# SYNTHETIC STAND-IN FOR THE BACKBONE (data/backbone/backbone_tz.rds)
#   Rscript scripts/make_synthetic_backbone.R                     -> synthetic/backbone_synthetic.rds
#   Rscript scripts/make_synthetic_backbone.R <out.rds>
#
# WHY: the simulator (scripts/questions/sb_simulate.R) borrows the STRUCTURE of the real field data: herd locations,
# sizes of the multi-cow herds, and each cow's recording pattern (calendar, lactation, age, lactation curve) with the
# true sizes of the fixed effects. That file is built from the real data and is not shared. This script makes a file
# with the SAME FORMAT from invented values, so that every step of the pipeline can be run and checked by anyone.
#
# WHAT IS INVENTED (nothing below comes from the real data, except the herd-size counts of the real data set of
# Houaga et al. 2026, which the paper itself reports only as 1,386 herds, more than 71 % of them with one cow):
#   herds    1,386 herds in 156 wards of 4 regions on a km grid of about 700 x 700 km; wards scattered around their
#            region centre (sd 55 km), herds clustered around their ward centre (sd 3 km), so most herds have a
#            neighbour within about a kilometre;
#   cows     1,894 cows; herd sizes as in the real data set (997 herds with 1 cow; 299 with 2, 71 with 3, 13 with
#            4, 5 with 5, 1 with 9), assigned to herds at random; one breed class;
#   records  3-37 monthly test-days per cow (mean about 10) over 61 calendar months: first lactation group 1/2/3 with
#            shares 0.54/0.34/0.12, days in milk 5-305, a 2-month dry period between lactations, age at test,
#            calving year-season; Legendre covariates leg1, leg2 of days in milk;
#   Xb_true  fixed-effect part of each record from invented effects: test month (sd 0.11), calving season (sd 0.30),
#            age and lactation-curve slopes within lactation group; centred, variance about 0.15.
# The numbers of the README come from the real backbone; with this stand-in the pipeline runs the same way, but the
# herd geography and record patterns differ, so the results will be similar in kind, not identical.
# =====================================================================================================================
args <- commandArgs(trailingOnly = TRUE)
out <- if (length(args)) args[1] else "synthetic/backbone_synthetic.rds"
set.seed(20261114L)

# ---- herds: 4 regions, 156 wards, 1,386 herds --------------------------------------------------------------------
n_herd <- 1386L; n_ward <- 156L
reg_c <- rbind(c(150, 560), c(520, 600), c(230, 180), c(560, 230))         # invented region centres (km)
ward_reg <- rep(1:4, times = c(48, 42, 36, 30))
ward_c <- reg_c[ward_reg, ] + matrix(rnorm(2 * n_ward, 0, 55), ncol = 2)
w_size <- as.integer(rmultinom(1, n_herd - n_ward, prob = rgamma(n_ward, 0.8))) + 1L   # every ward has a herd
herd_ward <- rep(seq_len(n_ward), w_size)
xy <- ward_c[herd_ward, ] + matrix(rnorm(2 * n_herd, 0, 3), ncol = 2)
herds <- data.frame(herd = sprintf("H%04d", seq_len(n_herd)), region = paste0("R", ward_reg[herd_ward]),
                    x_km = round(xy[, 1] - min(xy[, 1]) + 5, 3), y_km = round(xy[, 2] - min(xy[, 2]) + 5, 3),
                    ward = sprintf("W%03d", herd_ward), stringsAsFactors = FALSE)

# ---- cows: herd-size counts of the real data set, assigned to herds at random -------------------------------------
sizes <- sample(c(rep(1L, 997), rep(2L, 299), rep(3L, 71), rep(4L, 13), rep(5L, 5), 9L))
stopifnot(length(sizes) == n_herd, sum(sizes) == 1894L)
cows <- data.frame(animal_id = sprintf("C%04d", seq_len(sum(sizes))), herd = rep(herds$herd, sizes), dgrp = 1L,
                   stringsAsFactors = FALSE)
n_cow <- nrow(cows)

# ---- records: monthly test-days with lactations ------------------------------------------------------------------
n_month <- 61L
n_rec <- pmin(3L + rnbinom(n_cow, size = 2.2, mu = 7.2), 37L)
L0 <- sample(1:3, n_cow, replace = TRUE, prob = c(0.54, 0.34, 0.12))
age0 <- c(30, 45, 62)[L0] + rnorm(n_cow, 0, c(4, 5, 8)[L0])                 # age (months) at the first test
rows <- vector("list", n_cow)
for (i in seq_len(n_cow)) {
  dim <- sample(5:60, 1); lac <- L0[i]; rel_m <- 0L; recs <- NULL
  for (k in seq_len(n_rec[i])) {
    if (k > 1) {
      dim <- dim + 30.4; rel_m <- rel_m + 1L
      if (dim > 305) { dim <- sample(5:35, 1); lac <- min(lac + 1L, 3L); rel_m <- rel_m + 2L }   # dry period, calving
    }
    recs <- rbind(recs, c(rel_m, round(dim), lac))
  }
  span <- recs[nrow(recs), 1]
  m0 <- sample.int(max(1L, n_month - span), 1)
  m <- m0 + recs[, 1]
  keep <- m <= n_month
  rows[[i]] <- data.frame(animal_id = cows$animal_id[i], m = m[keep], dim = recs[keep, 2], lacgr = as.integer(recs[keep, 3]),
                          age = age0[i] + recs[keep, 1])
}
r <- do.call(rbind, rows)
r$herd <- cows$herd[match(r$animal_id, cows$animal_id)]
r$tyrmn <- sprintf("T%02d", r$m)
calv <- r$m - round(r$dim / 30.4)                                            # calendar month of calving
r$cyrsn <- sprintf("S%02d", (calv + 2L) %/% 3L + 1L)
r$year <- pmin(7L, 2L + as.integer(ceiling(r$m / 12)))
x <- pmin(pmax(2 * (r$dim - 5) / 300 - 1, -1), 1)
r$leg1 <- sqrt(3 / 2) * x
r$leg2 <- sqrt(5 / 2) * (3 * x^2 - 1) / 2
r$ageZ <- (r$age - mean(r$age)) / sd(r$age)
r$dgrp <- 1L

# ---- fixed effects (invented) ------------------------------------------------------------------------------------
lv_t <- sort(unique(r$tyrmn)); lv_c <- sort(unique(r$cyrsn))
b_t <- setNames(c(0, rnorm(length(lv_t) - 1, 0, 0.11)), lv_t)
b_c <- setNames(c(0, rnorm(length(lv_c) - 1, 0, 0.30)), lv_c)
b_age <- c(0.10, 0.05, 0.02); b_l1 <- c(-0.25, -0.30, -0.30); b_l2 <- c(-0.10, -0.10, -0.10)
xb <- b_t[r$tyrmn] + b_c[r$cyrsn] + b_age[r$lacgr] * r$ageZ + b_l1[r$lacgr] * r$leg1 + b_l2[r$lacgr] * r$leg2
r$Xb_true <- as.numeric(xb - mean(xb))
r$rec_id <- sprintf("R%05d", seq_len(nrow(r)))
records <- r[, c("animal_id", "herd", "tyrmn", "cyrsn", "dgrp", "lacgr", "ageZ", "leg1", "leg2", "year", "Xb_true", "rec_id")]
b <- c("(Intercept)" = 0, setNames(b_t[-1], paste0("tyrmn", lv_t[-1])), setNames(b_c[-1], paste0("cyrsn", lv_c[-1])),
       dgrp2 = 0, dgrp3 = 0, dgrp4 = 0, setNames(b_age, paste0("ageZ:lacgr", 1:3)),
       setNames(b_l1, paste0("lacgr", 1:3, ":leg1")), setNames(b_l2, paste0("lacgr", 1:3, ":leg2")))

bb <- list(herds = herds, cows = cows, records = records, b = b,
           fixed_formula = ~tyrmn + cyrsn + dgrp + ageZ:lacgr + leg1:lacgr + leg2:lacgr,
           meta = list(source = "SYNTHETIC stand-in: invented herd locations and records (no real data)",
                       crs = "none (km grid)", jitter_sd_km = 0, built = format(Sys.Date()),
                       n_herd = nrow(herds), n_cow = n_cow, n_rec = nrow(records)))
dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
saveRDS(bb, out)
nn <- apply(as.matrix(dist(herds[, c("x_km", "y_km")])) + diag(Inf, n_herd), 1, min)
nrc <- table(records$animal_id)
cat(sprintf("synthetic backbone -> %s\n  herds %d (single-cow %d) in %d wards, %d regions; extent %.0f x %.0f km; nearest herd median %.1f km\n  cows %d; records %d (per cow mean %.1f, median %g, range %d-%d); tyrmn %d levels, cyrsn %d levels; var(Xb) %.3f\n",
            out, n_herd, sum(sizes == 1), length(unique(herds$ward)), 4, diff(range(herds$x_km)), diff(range(herds$y_km)), median(nn),
            n_cow, nrow(records), mean(nrc), median(nrc), min(nrc), max(nrc), length(lv_t), length(lv_c), var(records$Xb_true)))
