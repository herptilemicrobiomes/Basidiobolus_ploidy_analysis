#!/usr/bin/env bash
# run_nquire.sh — sbatch launcher for the nQuire Nextflow head process.
#
# Submit:  sbatch run_nquire.sh
# Test:    nextflow run pipeline/ploidy/nquire.nf \
#              -c nextflow.config -profile test -stub-run --n_test 2

#SBATCH -p epyc
#SBATCH -N 1
#SBATCH -n 2
#SBATCH --mem 8G
#SBATCH -t 2-00:00:00
#SBATCH --job-name nf-nquire
#SBATCH -o logs/slurm/nf_nquire_%j.out
#SBATCH -e logs/slurm/nf_nquire_%j.err

set -euo pipefail

mkdir -p logs/slurm logs/nextflow

source /etc/profile.d/modules.sh 2>/dev/null || true
module load nextflow

nextflow run pipeline/ploidy/nquire.nf \
    -c nextflow.config \
    -profile nquire \
    -resume \
    "$@"
