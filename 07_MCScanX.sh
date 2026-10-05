#steps for MCScanX
#1 first rename the proteome sequences to add <species_name> in the front of each protein id
# then parse the gff files for MCScanX (see MCScanX github for guidelines)
# cat all renamed.fasta files as all.fasta
# do all Vs all blast --> all.blast (with diamond or ncbi blast)

diamond makedb --in all.fasta -d all
diamond blastp -d all -q all.fasta -o all.blast --outfmt 6 --max-target-seqs 5 --evalue 1e-10

#use the all.* files for MCScanX
