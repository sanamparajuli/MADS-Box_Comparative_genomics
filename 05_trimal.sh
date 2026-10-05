#!/bin/bash
# Run as a standalone script:  bash 05_trimal.sh   (after: conda activate trimal)
# Run as a SLURM script:       add your own #SBATCH header and 'module load trimal' (or conda activate trimal),
#                              then submit with: sbatch 05_trimal.sh
# Use an appropriate outgroup as needed.

INDIR=/path/to/alignments       # directory with *.aln.fasta files
OUTDIR=$INDIR/gt09
mkdir -p "$OUTDIR"

for aln in "$INDIR"/*.aln.fasta; do
    base=$(basename "$aln" .aln.fasta)
    trimal -in "$aln" -out "$OUTDIR/${base}.aln.trimal.gt09.fasta" -gt 0.9
done
