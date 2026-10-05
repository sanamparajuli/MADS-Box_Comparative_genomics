#!/bin/bash
# Run as a standalone script:  bash 04_mafft_alignment.sh
# Run as a SLURM script:       add your own #SBATCH header (cpus, mem, time, output/error)
#                              and submit with: sbatch 04_mafft_alignment.sh

set -euo pipefail

# module load mafft   # uncomment on a cluster that uses modules

THREADS=${SLURM_CPUS_PER_TASK:-8}

INDIR=/path/to/fasta_dir      # directory with the *.fasta files to align
OUTDIR=/path/to/alignments    # output directory
mkdir -p "${OUTDIR}"

for fa in "${INDIR}"/*.fasta; do
    # Skip unwanted file
    [[ "$(basename "${fa}")" == "unwanted.fasta" ]] && continue

    base=$(basename "${fa}" .fasta)
    out="${OUTDIR}/${base}.aln.fasta"

    echo "[$(date)] Aligning ${fa}"
    mafft --auto --thread "${THREADS}" --reorder "${fa}" > "${out}"
done
