#!/bin/bash -l
# ---------------------------------------------------------------------------------------------------------------------
# STEP 1 (this repository). One array task = one data set: simulate -> fit G, GH, GS, GHS -> extract -> (replicates
# 1-5 only) 5-fold cross-validation. Edit the "cd" line below to the folder where this repository is cloned, and
# create logs/ first. Submitted from the repository root:
#   INLA's default prior, 80 data sets at 1,000 and 2,000 cows (CV in replicates 1-5):
#     SIMSCRIPT=scripts/questions/sb_simulate.R FITFLAGS="--priors isidore" VCV=kfold VMODELS=G,GH,GS,GHS VREPS=5 \
#       sbatch --export=ALL,TASKS=cluster/tasks_SB1.txt --array=1-80%4 cluster/q_task.sbatch
#   4,000 cows (5 data sets, no CV; each needs many hours and more memory; run on a second node with the same versions):
#     SIMSCRIPT=scripts/questions/sb_simulate.R FITFLAGS="--priors isidore" \
#       sbatch --export=ALL,TASKS=cluster/tasks_SB1big.txt --array=1-5%2 --mem=64G --time=7-00:00:00 cluster/q_task.sbatch
#   PC-prior refit of the same 80 data sets (copy data/ and truth/ first; see README.md, "How to run"):
#     SIMSCRIPT=scripts/questions/sb_simulate.R FITFLAGS="--priors pc" \
#       sbatch --export=ALL,TASKS=cluster/tasks_SB1pc.txt --array=1-80%4 cluster/q_task.sbatch
# ---------------------------------------------------------------------------------------------------------------------
# One question-study task = simulate -> fit -> extract for one cell x replicate.
# The task table (TASKS=cluster/questions/tasks_<Q>.txt) has lines:  question | cell | rep | sim flags | models
# Everything lands in ~/Spacial_Cl/runs_q/<question>/<cell>_rNN/ ; each step is appended to STEPS.md there.
#   sbatch --export=ALL,TASKS=cluster/questions/tasks_Q1Q2.txt --array=1-24%4 cluster/questions/q_task.sbatch
# Optional: SIMSCRIPT (default scripts/questions/q_simulate.R; study SB: scripts/questions/sb_simulate.R) and
# FITFLAGS (extra q_fit.R / q_validate.R flags; study SB: "--priors isidore").
# Optional cross-validation after the fit: VCV (q_validate.R --cv, e.g. kfold) with VMODELS (default: MODELS) and
# VREPS (only replicates <= VREPS are validated; default all). A run whose fit is already extracted (output/q_vars.csv)
# is not refitted, so a table can add validation to finished runs.
#SBATCH --job-name=q
#SBATCH --partition=batch
#SBATCH --cpus-per-task=5
#SBATCH --mem=24G
#SBATCH --time=2-00:00:00
#SBATCH --output=logs/q_%x_%A_%a.out
set -uo pipefail
module load R/4.5
export R_LIBS_USER=$HOME/R/x86_64-pc-linux-gnu-library/4.5
export INLA_THREADS="${SLURM_CPUS_PER_TASK}:1"
cd "$HOME/Spacial_Cl"
LINE=$(grep -v '^#' "$TASKS" | grep -v '^\s*$' | sed -n "${SLURM_ARRAY_TASK_ID}p")
IFS='|' read -r Q CELL REP SIMFLAGS MODELS <<< "$LINE"
Q=$(echo $Q); CELL=$(echo $CELL); REP=$(echo $REP); SIMFLAGS=$(echo $SIMFLAGS); MODELS=$(echo $MODELS)
RUN="runs_q/${Q}/${CELL}_r$(printf %02d "$REP")"
mkdir -p "$RUN"
STEPS="$RUN/STEPS.md"
step() {   # step <name> <command...> : run, time and record in STEPS.md
  local name=$1; shift
  local t0=$(date +%s)
  echo "- $(date '+%F %T') **${name}** started on $(hostname) (job ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}): \`$*\`" >> "$STEPS"
  "$@"; local rc=$?
  echo "  - finished rc=${rc} after $(( $(date +%s) - t0 )) s" >> "$STEPS"
  return $rc
}
[ -s "$STEPS" ] || printf "# %s / %s / replicate %s\n\nSimulation flags: \`%s\`  \nModels: \`%s\`\n\n## Steps\n\n" "$Q" "$CELL" "$REP" "$SIMFLAGS" "$MODELS" > "$STEPS"
echo "task ${SLURM_ARRAY_TASK_ID}: $RUN | sim: $SIMFLAGS | models: $MODELS | $(date)"
if [ ! -s "$RUN/truth/true_animals.csv" ]; then
  step simulate Rscript "${SIMSCRIPT:-scripts/questions/q_simulate.R}" --run "$RUN" --rep "$REP" $SIMFLAGS || exit 1
fi
if [ ! -s "$RUN/output/q_vars.csv" ]; then
  step fit Rscript scripts/questions/q_fit.R --run "$RUN" --qmodels "$MODELS" ${FITFLAGS:-} || exit 2
  step extract Rscript scripts/questions/q_extract.R --run "$RUN" || exit 3
fi
if [ -n "${VCV:-}" ] && [ "$REP" -le "${VREPS:-999}" ] && ! grep -qs "Q-VALIDATE finished" "$RUN/output/q_validate_${VCV}_log.txt"; then
  step validate Rscript scripts/questions/q_validate.R --run "$RUN" --cv "$VCV" --qmodels "${VMODELS:-$MODELS}" ${FITFLAGS:-} || exit 4
fi
echo "done $(date)"
