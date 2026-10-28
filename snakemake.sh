#!/bin/bash
# Local run: ./snakemake.sh [snakemake args], e.g. ./snakemake.sh -n for a dry run.
set -eo pipefail

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

snakemake --workflow-profile workflow/profiles/default -s workflow/Snakefile "$@" "${LOCAL_CONFIG[@]}"
