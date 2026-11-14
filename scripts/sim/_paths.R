# =====================================================================================================================
# STEP 1 NOTE: every step-1 script is called with --run runs_q/<study>/<cell>_rNN. That folder then holds one data set:
#   data/    what the analyst sees (herds, animals, records, G, G inverse);
#   truth/   the simulated true values (never read by the fitting code);
#   output/  fits, metrics and logs.
# Run every script from the repository root. The --scale, --models and --pc-priors options below are not used in step 1.
# =====================================================================================================================

# Scenario paths, shared by every stage.
#
# Synthetic map (scripts/sim/00): `--scale k` multiplies every count (herds, cows, sires,
# dams, breeding populations, records) by k on the same 80 x 80 km map, mesh and founder
# haplotypes. k = 1 is the specified size and uses the base folders:
#   k = 1  -> data/sim,      truth,      output/sim
#   k = 10 -> data/sim_x10,  truth_x10,  output/sim_x10   (paper-size: 1,200 herds, 1,800 cows)
#
# Real-data backbone runs (scripts/real/10): `--run runs/<tag>` puts everything of one
# replicate under that folder:  runs/<tag>/data,  runs/<tag>/truth,  runs/<tag>/output
#
# Fit-stage variants write next to the default output so they sit side by side:
#   `--models gp`   Houaga's GP series (adds a permanent-environment effect)  -> <output>_gp
#   `--pc-priors`   PC priors on every precision                               -> <output>_pc
#   both                                                                       -> <output>_gp_pc

if (!dir.exists("scripts/sim")) stop("Run from the project root (folder that contains scripts/sim/).")

cli_args <- commandArgs(trailingOnly = TRUE)
scale_k <- 1L
if ("--scale" %in% cli_args) scale_k <- as.integer(cli_args[match("--scale", cli_args) + 1])
if (is.na(scale_k) || scale_k < 1) stop("--scale must be followed by a whole number >= 1")

scenario_suffix <- if (scale_k == 1) "" else paste0("_x", scale_k)
data_dir <- paste0("data/sim", scenario_suffix)
truth_dir <- paste0("truth", scenario_suffix)
out_dir <- paste0("output/sim", scenario_suffix)
scale_arg <- if (scale_k == 1) "" else paste(" --scale", scale_k)

run_dir <- NULL
if ("--run" %in% cli_args) {
  run_dir <- sub("/+$", "", cli_args[match("--run", cli_args) + 1])
  if (is.na(run_dir) || startsWith(run_dir, "--")) stop("--run must be followed by a folder, e.g. runs/bb_r01")
  data_dir <- file.path(run_dir, "data")
  truth_dir <- file.path(run_dir, "truth")
  out_dir <- file.path(run_dir, "output")
  scale_arg <- paste(" --run", run_dir)
}

# `--models gp` fits Houaga's GP series (G + permanent environment) instead of G/GH/GS/GHS
# and writes to <output>_gp.
model_set <- if ("--models" %in% cli_args) cli_args[match("--models", cli_args) + 1] else "houaga"
if (!model_set %in% c("houaga", "gp")) stop("--models must be houaga or gp")
if (model_set == "gp") {
  out_dir <- paste0(out_dir, "_gp")
  scale_arg <- paste0(scale_arg, " --models gp")
}

pc_priors <- "--pc-priors" %in% cli_args
if (pc_priors) {
  out_dir <- paste0(out_dir, "_pc")
  scale_arg <- paste0(scale_arg, " --pc-priors")
}

# Founder haplotypes are shared by all scales (runMacs is not reproducible; note Part 2 Step 6)
founder_file <- "truth/founder_haplotypes_macs.rds"
