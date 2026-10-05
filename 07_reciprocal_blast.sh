#!/bin/bash
# Run as a standalone script:  bash 07_reciprocal_blast.sh        (loops over all targets)
# Run as a SLURM array job:    add your own #SBATCH header incl. --array=0-<N-1> (N = number of targets)
#                              and submit with: sbatch 07_reciprocal_blast.sh
# Requires BLAST databases (makeblastdb -dbtype prot) for the reference and every target, in the working dir.

# module load ncbi-blast   # uncomment on a cluster that uses modules

cd /path/to/working_dir || exit 1

THREADS=${SLURM_CPUS_PER_TASK:-8}

REF=reference_typeII.fasta
REF_DB=reference_typeII

targets=(
species1_typeII.fasta
species2_typeII.fasta
species3_typeII.fasta
)

run_rbh() {
    TARGET=$1
    PREFIX=$(basename "$TARGET" .fasta)
    TARGET_DB=$PREFIX

    mkdir -p rbh
    pushd rbh > /dev/null || exit 1

    # Reference -> target
    blastp \
      -query ../$REF \
      -db ../$TARGET_DB \
      -out reference_vs_${PREFIX}.blastp.tsv \
      -evalue 1e-5 \
      -max_target_seqs 5 \
      -num_threads "$THREADS" \
      -outfmt 6

    # Target -> reference
    blastp \
      -query ../$TARGET \
      -db ../$REF_DB \
      -out ${PREFIX}_vs_reference.blastp.tsv \
      -evalue 1e-5 \
      -max_target_seqs 5 \
      -num_threads "$THREADS" \
      -outfmt 6

    popd > /dev/null
}

if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    run_rbh "${targets[$SLURM_ARRAY_TASK_ID]}"   # SLURM: one target per array task
else
    for t in "${targets[@]}"; do run_rbh "$t"; done   # standalone: all targets
fi
