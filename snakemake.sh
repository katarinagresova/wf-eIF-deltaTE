#!/bin/bash
# Local run: ./snakemake.sh [snakemake args], e.g. ./snakemake.sh -n for a dry run.
set -eo pipefail

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate snake

# Untracked machine-specific paths, layered over config/config.yaml (README.md)
LOCAL_CONFIG=()
[[ -f config/local/config.yaml ]] && LOCAL_CONFIG=(--configfile config/local/config.yaml)

snakemake --workflow-profile workflow/profiles/default -s workflow/Snakefile "$@" "${LOCAL_CONFIG[@]}"
