#!/bin/bash
# Run as a standalone script:  bash 03_Interproscan.sh            (loops over all FASTA files)
# Run as a SLURM array job:    add your own #SBATCH header incl. --array=0-<N-1> (N = number of FASTA files)
#                              and submit with: sbatch 03_Interproscan.sh

# module load interproscan   # uncomment on a cluster that uses modules

THREADS=${SLURM_CPUS_PER_TASK:-8}

# directory containing the protein FASTA files
cd /path/to/interpro_dir || exit 1

# Collect FASTA files
fastas=(*.fasta)

run_interpro() {
    fasta=$1
    base=$(basename "$fasta" .fasta)
    interproscan.sh \
        --cpu "$THREADS" \
        -t p \
        -i "$fasta" \
        -o "${base}_interpro.tsv" \
        -f TSV \
        -appl Pfam,SMART,CDD
}

if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    run_interpro "${fastas[$SLURM_ARRAY_TASK_ID]}"   # SLURM: one FASTA per array task
else
    for f in "${fastas[@]}"; do run_interpro "$f"; done   # standalone: all FASTA files
fi
