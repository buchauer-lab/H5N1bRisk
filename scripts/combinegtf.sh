#!/bin/bash

HUMAN_GTF="Homo_sapiens.GRCh38.113.gtf"

cd /data/cephfs-1/work/groups/buchauer/h5n1risk/genomes/individualGenomes

# Loop through all GTF files
for VIRUS_FILE in *vir.gtf; do
    # Skip the human GTF file
    [[ "$VIRUS_FILE" == "$HUMAN_GTF" ]] && continue

    OUTPUT_FILE="../combinedGenomes/${VIRUS_FILE%.gtf}_GRCh38.gtf.gz"

    echo "Combining: $VIRUS_FILE + $HUMAN_GTF -> $OUTPUT_FILE"

    # Combine the human GTF with the virus GTF (excluding virus GTF headers), then compress
    cat "$HUMAN_GTF" <(grep -v '^#' "$VIRUS_FILE") > "$OUTPUT_FILE"
done
