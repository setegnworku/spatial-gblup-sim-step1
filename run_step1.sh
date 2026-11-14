#!/usr/bin/env bash
# =====================================================================================================================
# Run step 1 on ONE machine (no Slurm), from the repository root. Each data set goes through every step:
#   simulate (genome -> population -> G -> herds -> field -> one record per cow)  ->  fit G, GH, GS, GHS  ->
#   score against the truth  ->  5-fold cross-validation (replicates <= CV_MAXREP)  ->  pool all finished data sets.
# Steps already done are skipped, so the script can be re-run after an interruption.
#
#   bash run_step1.sh                                   # one data set: 1,000 cows, 70 % single-cow herds, replicate 1
#   CELLS="n1000_p50 n1000_p70" REPS="1 2" bash run_step1.sh
#   CELLS="n1000_p50 n1000_p60 n1000_p70 n1000_p80 n2000_p50 n2000_p60 n2000_p70 n2000_p80" REPS="$(seq 1 10)" \
#     bash run_step1.sh                                 # the 80 data sets of the paper (days on one machine)
#   CELLS="n4000_p70" REPS="$(seq 1 5)" CV=no bash run_step1.sh          # the 4,000-cow data sets (no CV)
#   STUDY=SB1_pc PRIORS=pc CV=no CELLS="..." REPS="..." bash run_step1.sh   # PC-prior refit of the SAME data sets
#
# Needs: founders in runs/founders/ (cp founders/*.rds runs/founders/) and a backbone in data/backbone/backbone_tz.rds
# (the real one, or the synthetic stand-in: cp synthetic/backbone_synthetic.rds data/backbone/backbone_tz.rds).
# Environment: RSCRIPT (default Rscript), INLA_THREADS (default 4:1).
# =====================================================================================================================
set -euo pipefail
RSCRIPT=${RSCRIPT:-Rscript}
STUDY=${STUDY:-SB1_one_record}
PRIORS=${PRIORS:-isidore}            # isidore = INLA's default prior (as Houaga et al.); pc = PC priors
CV=${CV:-yes}; CV_MAXREP=${CV_MAXREP:-5}
CELLS=${CELLS:-n1000_p70}; REPS=${REPS:-1}
export INLA_THREADS=${INLA_THREADS:-4:1}
[ -d scripts/questions ] || { echo "run from the repository root"; exit 1; }
[ -s data/backbone/backbone_tz.rds ] || { echo "missing data/backbone/backbone_tz.rds (see REPRODUCE.md, step 2)"; exit 1; }
mkdir -p runs/founders
step() { local t0=$SECONDS; echo "[$(date '+%F %T')] $*"; "$@"; echo "    done in $((SECONDS - t0)) s"; }

for cell in $CELLS; do
  ncow=${cell#n}; ncow=${ncow%%_*}; p=${cell##*_p}
  for rep in $REPS; do
    RUN=runs_q/$STUDY/${cell}_r$(printf %02d "$rep")
    SRC=runs_q/SB1_one_record/${cell}_r$(printf %02d "$rep")
    if [ ! -s "$RUN/truth/true_animals.csv" ]; then
      if [ "$STUDY" != "SB1_one_record" ] && [ -s "$SRC/truth/true_animals.csv" ]; then
        mkdir -p "$RUN" && cp -r "$SRC/data" "$SRC/truth" "$RUN/"          # refit the same data set
        echo "[$(date '+%F %T')] copied data/ and truth/ from $SRC"
      else
        step "$RSCRIPT" scripts/questions/sb_simulate.R --run "$RUN" --rep "$rep" --ncow "$ncow" --p1 "0.$p" --records one --spe 0
      fi
    fi
    if [ ! -s "$RUN/output/q_vars.csv" ]; then
      step "$RSCRIPT" scripts/questions/q_fit.R --run "$RUN" --qmodels G,GH,GS,GHS --priors "$PRIORS"
      step "$RSCRIPT" scripts/questions/q_extract.R --run "$RUN"
    fi
    if [ "$CV" = yes ] && [ "$rep" -le "$CV_MAXREP" ] && ! grep -qs "Q-VALIDATE finished" "$RUN/output/q_validate_kfold_log.txt"; then
      step "$RSCRIPT" scripts/questions/q_validate.R --run "$RUN" --cv kfold --qmodels G,GH,GS,GHS --priors "$PRIORS"
    fi
  done
done
step "$RSCRIPT" scripts/questions/sb_aggregate.R --question "$STUDY"
echo "pooled results: output/questions/$STUDY/  (compare with results/ of this repository)"
