#!/bin/bash
# Run as a standalone script:  bash 06_iqtree_phylogeny.sh        (loops over all alignments)
# Run as a SLURM array job:    add your own #SBATCH header incl. --array=0-<N-1> (N = number of alignments)
#                              and submit with: sbatch 06_iqtree_phylogeny.sh

# module load iqtree   # uncomment on a cluster that uses modules

cd /path/to/trimmed_alignments || exit 1

alignments=( *.gt09.fasta )

run_iqtree() {
    aln_file=$1
    prefix=$(basename "$aln_file" aln.gt09.fasta)
    iqtree3 -s "$aln_file" \
            -m TEST \
            -bb 1000 \
            -nt AUTO \
            -pre "${prefix}_iqtree"
}

if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    run_iqtree "${alignments[$SLURM_ARRAY_TASK_ID]}"   # SLURM: one alignment per array task
else
    for a in "${alignments[@]}"; do run_iqtree "$a"; done   # standalone: all alignments
fi
