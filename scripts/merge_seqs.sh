#!/bin/bash

RAW_DIR="/data/cephfs-1/work/groups/buchauer/users/"${USER}"/h5n1risk/raw_run1/"
RE_RAW_DIR="/data/cephfs-1/work/groups/buchauer/users/"${USER}"/h5n1risk/raw_run2/"
MERGED_DIR="/data/cephfs-1/work/groups/buchauer/users/"${USER}"/h5n1risk/raw/"

mkdir -p "$MERGED_DIR"

for file in "$RAW_DIR"/*.fastq.gz; do

    filename=$(basename "$file")
    refile="$RE_RAW_DIR/$filename"
    mergedfile="$MERGED_DIR/$filename"

    # Check if the matching file exists in re_raw
    if [[ -f "$refile" ]]; then
        # Merge both files and write to merged directory
        cat "$file" "$refile" > "$mergedfile"
    else
        echo "Warning: Matching file not found for $filename in $RE_RAW_DIR"
        cp "$file" "$mergedfile"
    fi
done