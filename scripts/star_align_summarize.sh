#!/bin/bash

ROOTPATH="/data/cephfs-1/work/groups/buchauer/users/${USER}/h5n1risk/"
ALIGNMENTPATH="${ROOTPATH}/alignment/"
OUTFILE="$ALIGNMENTPATH/star_summary.csv"

# Output CSV header
echo -e "Sample,InputReads,AvgReadLen,UniquelyMappedPct,MultiMappedPct,UnmappedTooShortPct,MismatchRatePct" > $OUTFILE

# Loop through all subdirectories in ALIGNMENTPATH
for sample_dir in "${ALIGNMENTPATH}"*/; do
    sample=$(basename "$sample_dir")
    log_file="${sample_dir}/${sample}_Log.final.out"
    
    # Extract metrics from the log file
    if [[ -f "$log_file" ]]; then
        input_reads=$(grep "Number of input reads" "$log_file" | awk '{print $NF}')
        avg_read_len=$(grep "Average input read length" "$log_file" | awk '{print $NF}')
        uniq_pct=$(grep "Uniquely mapped reads %" "$log_file" | awk '{print $NF}' | sed 's/%//')
        multi_pct=$(grep "% of reads mapped to multiple loci" "$log_file" | awk '{print $NF}' | sed 's/%//')
        short_pct=$(grep "% of reads unmapped: too short" "$log_file" | awk '{print $NF}' | sed 's/%//')
        mismatch_pct=$(grep "Mismatch rate per base" "$log_file" | awk '{print $NF}' | sed 's/%//')

        # Append to CSV
        echo -e "${sample},${input_reads},${avg_read_len},${uniq_pct},${multi_pct},${short_pct},${mismatch_pct}" >> $OUTFILE
    else
        echo "Warning: Log file not found for sample $sample" >&2
    fi
done
