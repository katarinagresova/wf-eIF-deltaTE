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

# Untracked machine-specific paths, layered over config/config.yaml (README.md)
LOCAL_CONFIG=()
[[ -f config/local/config.yaml ]] && LOCAL_CONFIG=(--configfile config/local/config.yaml)

snakemake --workflow-profile workflow/profiles/slurm -s workflow/Snakefile "$@" "${LOCAL_CONFIG[@]}"
