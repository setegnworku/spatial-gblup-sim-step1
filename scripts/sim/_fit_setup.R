# =====================================================================================================================
# STEP 1 NOTE: shared set-up for q_fit.R and q_validate.R. It reads ONLY the analyst's data/ folder (never truth/), so
# the models see what a real analyst would see. In study SB (one breed) the fixed part is Houaga's without the breed
# class: tyrmn + cyrsn + ageZ:lacgr + leg1:lacgr + leg2:lacgr (test year-month, calving year-season, age and the
# lactation curve within lactation group). The mesh is rebuilt from the simulated herds with Houaga's recipe
# (_mesh.R), and the field gets his PC-Matern prior. The model list, fit_bru() and the --pc-priors switch further down
# belong to the project's older pipeline: q_fit.R and q_validate.R define the step-1 models and priors themselves.
# =====================================================================================================================

# Shared setup for the fitting scripts (01_fit, 03_validate): data, mesh, Matern, components.
# Reads ONLY the data folder (data/sim*/ or runs/<tag>/data/) — never truth/. Spec: houaga_gblup_spatial_sim_notes.md, Parts 4 and 6.
#
# Two inlabru >= 2.13 pitfalls handled here (note Part 6):
#   1. A factor `cow` passed to model = "generic" is silently disconnected from the
#      data (all EBVs = 0, only a "'-' not meaningful for factors" warning). We pass
#      an integer index into Ginv instead.
#   2. inlabru adds an automatic Intercept; fixed_effects(main = ~ dgrp + tyrmn)
#      already has one. The component formula starts with -1 to drop the extra one.

source("scripts/sim/_paths.R")   # data_dir, truth_dir, out_dir, model_set, pc_priors
need <- c(file.path(data_dir, "herds.csv"), file.path(data_dir, "animals.csv"), file.path(data_dir, "records.csv"), file.path(data_dir, "Ginv.rds"))
if (!all(file.exists(need))) stop("Missing ", paste(need[!file.exists(need)], collapse = ", "),
                                  ". Run scripts/sim/00_simulate_houaga_ghs.R (or scripts/real/10_simulate_backbone.R) first.")

suppressPackageStartupMessages({
  library(INLA)
  library(inlabru)
  library(fmesher)
  library(Matrix)
})

herds <- read.csv(file.path(data_dir, "herds.csv"), colClasses = c(herd = "character", region = "character"))
animals <- read.csv(file.path(data_dir, "animals.csv"), colClasses = c(animal_id = "character", herd = "character"))
records <- read.csv(file.path(data_dir, "records.csv"), colClasses = c(rec_id = "character", animal_id = "character",
                                                            herd = "character"))
Ginv <- readRDS(file.path(data_dir, "Ginv.rds"))
n_animal <- nrow(Ginv)
n_herd <- nrow(herds)

dat <- records
dat$region <- factor(herds$region[match(dat$herd, herds$herd)])
dat$x_km <- herds$x_km[match(dat$herd, herds$herd)]
dat$y_km <- herds$y_km[match(dat$herd, herds$herd)]
dat$dgrp <- factor(animals$dgrp[match(dat$animal_id, animals$animal_id)], levels = 1:4)
dat$tyrmn <- factor(dat$tyrmn, levels = sort(unique(dat$tyrmn)))
# Real-data backbone runs carry Houaga's full set of covariates; use his full fixed part then.
full_fixed <- all(c("cyrsn", "lacgr", "ageZ", "leg1", "leg2") %in% names(dat))
if (full_fixed) {
  dat$cyrsn <- factor(dat$cyrsn, levels = sort(unique(dat$cyrsn)))
  dat$lacgr <- factor(dat$lacgr, levels = sort(unique(dat$lacgr)))
}
one_breed <- length(unique(dat$dgrp)) == 1      # study SB (one breed): no breed-class effect
fixed_formula <- if (full_fixed) {
  if (one_breed) ~ tyrmn + cyrsn + ageZ:lacgr + leg1:lacgr + leg2:lacgr else
    ~ tyrmn + cyrsn + dgrp + ageZ:lacgr + leg1:lacgr + leg2:lacgr        # Houaga 1.Scripts_Spatial_Modelling.R
} else {
  if (one_breed) ~ tyrmn else ~ dgrp + tyrmn
}
# Integer indices: position of the cow in Ginv, position of the herd in herds.csv.
# Genotyped animals without records (the 20 sires) remain as indices with no data,
# exactly like Houaga's unused GRMInv factor levels — they still get EBVs through G.
dat$cowI <- match(dat$animal_id, rownames(Ginv))
dat$herdI <- match(dat$herd, herds$herd)
stopifnot(!anyNA(dat$cowI), !anyNA(dat$herdI), !anyNA(dat$dgrp), !anyNA(dat$x_km),
          isSymmetric(unname(Ginv)))

# Mesh: built by the same function as the simulation (note Part 4.10). PC priors exactly as Houaga.
source("scripts/sim/_mesh.R")
loc_herd <- as.matrix(herds[, c("x_km", "y_km")])
mesh <- build_mesh(loc_herd, mesh_type_of(data_dir))
matern <- inla.spde2.pcmatern(mesh,
                              prior.sigma = c(sqrt(0.25), 0.5),   # P(sigma > 0.5) = 0.5
                              prior.range = c(50, 0.8))           # P(range < 50 km) = 0.8

# Variance priors. Default (as Houaga): INLA's loggamma(1, 5e-5) on every precision.
# --pc-priors (Selle et al. 2020 style): PC priors, P(sd > 0.5) = 0.5 for animal, permanent
# environment and herd, P(sd > 1) = 0.5 for the residual — weakly informative, and unlike the
# default they do not pull a weakly identified variance towards zero.
hyper_re <- list(prec = list(prior = "pc.prec", param = c(0.5, 0.5)))
hyper_resid <- list(prec = list(prior = "pc.prec", param = c(1, 0.5)))
pc_txt <- if (pc_priors) ", hyper = hyper_re" else ""

term_txt <- c(
  fixed = paste0("fixed_effects(main = ", paste(deparse(fixed_formula), collapse = ""), ", model = \"fixed\")"),
  animal = paste0("animal(main = cowI, model = \"generic\", Cmatrix = Ginv, mapper = bru_mapper_index(n_animal)", pc_txt, ")"),
  perm = paste0("perm(main = cowI, model = \"iid\", mapper = bru_mapper_index(n_animal)", pc_txt, ")"),
  herd = paste0("herd(main = herdI, model = \"iid\", mapper = bru_mapper_index(n_herd)", pc_txt, ")"),
  field = "field(main = cbind(x_km, y_km), model = matern)"
)
model_terms <- if (model_set == "gp") {
  list(GP = c("fixed", "animal", "perm"), GPH = c("fixed", "animal", "perm", "herd"),
       GPS = c("fixed", "animal", "perm", "field"), GPHS = c("fixed", "animal", "perm", "herd", "field"))
} else {
  list(G = c("fixed", "animal"), GH = c("fixed", "animal", "herd"),
       GS = c("fixed", "animal", "field"), GHS = c("fixed", "animal", "herd", "field"))
}
components <- lapply(model_terms, function(tt)
  as.formula(paste("~ -1 +", paste(term_txt[tt], collapse = " + ")), env = globalenv()))

fit_bru <- function(cmp, data, int_strategy) {
  obs <- if (pc_priors) {
    bru_obs(formula = y ~ ., family = "gaussian", data = data,
            control.family = list(hyper = hyper_resid))
  } else {
    bru_obs(formula = y ~ ., family = "gaussian", data = data)   # Houaga: like(family = "Gaussian", ...)
  }
  bru(cmp, obs,
      options = list(
        control.inla = list(int.strategy = int_strategy),
        control.compute = list(waic = TRUE, dic = TRUE, config = TRUE),
        verbose = FALSE
      ))
}
