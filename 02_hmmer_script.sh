#!/bin/bash
# Run as a standalone script:  bash 02_hmmer_script.sh            (loops over all proteomes)
# Run as a SLURM array job:    add your own #SBATCH header incl. --array=0-<N-1> (N = number of proteomes)
#                              and submit with: sbatch 02_hmmer_script.sh

# module load hmmer   # uncomment on a cluster that uses modules

THREADS=${SLURM_CPUS_PER_TASK:-4}

HMM=PF00319.hmm   # MADS-box Pfam HMM; repeat with the K-box Pfam HMM

proteomes=(
species1_proteome.fasta
species2_proteome.fasta
species3_proteome.fasta
)

run_hmmer() {
    proteome=$1
    prefix=$(basename "$proteome" .fasta)
    hmmsearch \
        --cpu "$THREADS" \
        --tblout "${prefix}_hmmer.txt" \
        -E 0.1 \
        "$HMM" \
        "$proteome"
}

if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    run_hmmer "${proteomes[$SLURM_ARRAY_TASK_ID]}"   # SLURM: one proteome per array task
else
    for p in "${proteomes[@]}"; do run_hmmer "$p"; done   # standalone: all proteomes
fi
