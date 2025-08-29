#!/bin/bash

FASTQPATH="/data/cephfs-1/work/groups/buchauer/users/${USER}/h5n1risk/raw/"
QCPATH="${FASTQPATH}QC"
mkdir -p ${QCPATH}

fastqc -o $QCPATH -t 12 ${FASTQPATH}*.fastq.gz

multiqc -o $QCPATH $QCPATH