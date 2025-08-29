#!/bin/bash

HUMAN_FASTA="Homo_sapiens.GRCh38.dna.primary_assembly.fa"

cd /data/cephfs-1/work/groups/buchauer/h5n1risk/genomes/individualGenomes

for VIRUS_FILE in *.fasta; do
    # Define output file name
    OUTPUT_FILE="../combinedGenomes/${VIRUS_FILE%.fasta}_GRCh38.fa"

    echo "Combining $VIRUS_FILE with $HUMAN_FASTA -> $OUTPUT_FILE"

    # Concatenate the virus FASTA and human FASTA, then compress
    cat "$VIRUS_FILE" "$HUMAN_FASTA" > "$OUTPUT_FILE"
done

