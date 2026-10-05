#!/bin/bash
# Run as a standalone script:  bash 01_Braker_genome_annotation.sh
# Run as a SLURM script:       add your own #SBATCH header (partition, cpus, mem, time, output/error)
#                              and submit with: sbatch 01_Braker_genome_annotation.sh
# Note: BRAKER is memory-intensive and slow; allow plenty of RAM and time.

# module load singularity   # uncomment on a cluster that uses modules

THREADS=${SLURM_CPUS_PER_TASK:-8}   # uses SLURM allocation if present, otherwise 8

BASE=/path/to/genome_dir            # directory containing your genome and inputs

# input paths
GENOME=$BASE/your_genome.fa
PROT=$BASE/protein_evidence.fasta   # e.g. OrthoDB Viridiplantae proteins
AUGCFG=$BASE/augustus_config
SIF=$BASE/braker3.sif               # download this sif image first
WORKDIR=$BASE/braker

mkdir -p "$WORKDIR"

singularity exec \
  --bind "$BASE":/mnt \
  --bind "$AUGCFG":/opt/Augustus/config \
  "$SIF" \
  braker.pl \
    --genome=/mnt/your_genome.fa \
    --prot_seq=/mnt/protein_evidence.fasta \
    --softmasking \
    --species=<your_species> --useexisting \
    --gff3 \
    --threads="$THREADS" \
    --workingdir=/mnt/braker \
    --AUGUSTUS_CONFIG_PATH=/opt/Augustus/config
# --prot_seq     : protein evidence
# --useexisting  : use existing Augustus species config path
