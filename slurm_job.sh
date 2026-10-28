#!/bin/bash
# Driver job: runs snakemake with the SLURM profile, which submits one SLURM
# job per rule instance. Submit from the repo root with: sbatch slurm_job.sh

#SBATCH --job-name=wf-eIF-deltaTE
#SBATCH --output=logs/slurm_%j.out
#SBATCH --error=logs/slurm_%j.err
#SBATCH --time=2-00:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate snake

# Untracked machine-specific paths, layered over config/config.yaml (README.md);
# then $RUN_CONFIG if set (another dataset's config, e.g.
# RUN_CONFIG=config/local/R002/config.yaml), in the same --configfile: a second
# --configfile replaces the first
CONFIGS=()
[[ -f config/local/config.yaml ]] && CONFIGS+=(config/local/config.yaml)
[[ -n "$RUN_CONFIG" ]] && CONFIGS+=("$RUN_CONFIG")
LOCAL_CONFIG=()
(( ${#CONFIGS[@]} )) && LOCAL_CONFIG=(--configfile "${CONFIGS[@]}")

snakemake --workflow-profile workflow/profiles/slurm -s workflow/Snakefile "$@" "${LOCAL_CONFIG[@]}"
