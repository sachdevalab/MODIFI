#!/bin/bash
# Set up annotation tools for the Klebsiella-ECE deep dive.
# Creates a conda env with mob_suite + AMRFinderPlus + bakta, and downloads their DBs to borg
# (NEVER /home, which is a small SSD). Run in background and watch the log.
set -euo pipefail

export HF_HUB_ENABLE_XET=0
CONDA=/shared/software/miniconda3/latest
DBROOT=/home/shuaiw/borg/revision/kp_eces/dbs
# The shared pkgs cache is read-only for us; use a writable cache under borg.
export CONDA_PKGS_DIRS=/home/shuaiw/borg/revision/kp_eces/conda_pkgs
mkdir -p "$DBROOT" "$CONDA_PKGS_DIRS"

echo "[setup] $(date) creating env kp_eces (mob_suite + amrfinderplus)"
# NOTE: do NOT list seqkit explicitly - it is pulled transitively and an explicit + transitive
# request triggers a conda 'record already exists in the prefix' rollback. seqkit lives in env tldr.
# --override-channels avoids the implicit defaults/repo.anaconda.com channel (DNS-flaky, unneeded).
"$CONDA/bin/mamba" create -y -n kp_eces --override-channels -c conda-forge -c bioconda \
    mob_suite ncbi-amrfinderplus blast

ENV=/home/shuaiw/miniconda3/envs/kp_eces

echo "[setup] $(date) mob_suite DB -> $DBROOT/mob_suite"
mkdir -p "$DBROOT/mob_suite"
"$ENV/bin/mob_init" --database_directory "$DBROOT/mob_suite" || \
    echo "[setup] mob_init returned nonzero (may already be initialised) - continuing"

echo "[setup] $(date) AMRFinderPlus DB -> $DBROOT/amrfinder"
mkdir -p "$DBROOT/amrfinder"
"$ENV/bin/amrfinder" -u -d "$DBROOT/amrfinder" || \
    "$ENV/bin/amrfinder_update" -d "$DBROOT/amrfinder" || \
    echo "[setup] amrfinder update returned nonzero - check log"

echo "[setup] $(date) creating env bakta_env + light DB"
"$CONDA/bin/mamba" create -y -n bakta_env --override-channels -c conda-forge -c bioconda bakta
BENV=/home/shuaiw/miniconda3/envs/bakta_env
mkdir -p "$DBROOT/bakta"
"$BENV/bin/bakta_db" download --output "$DBROOT/bakta" --type light || \
    echo "[setup] bakta_db download returned nonzero - check log"

echo "[setup] $(date) DONE. Tool locations:"
ls -la "$ENV/bin/mob_typer" "$ENV/bin/amrfinder" "$BENV/bin/bakta" 2>&1 || true
echo "[setup] DB dirs:"; du -sh "$DBROOT"/* 2>/dev/null || true
