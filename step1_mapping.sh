````
bwa index ref.fa
samtools faidx ref.fa

#!/bin/bash
set -euo pipefail

REF="ref.fa"
RAW_DIR="raw_reads"
OUT_DIR="bam"
THREADS=8

mkdir -p $OUT_DIR

for R1 in ${RAW_DIR}/*_R1.fastq.gz
do
    SAMPLE=$(basename $R1 _R1.fastq.gz)
    R2=${RAW_DIR}/${SAMPLE}_R2.fastq.gz

    echo "Mapping $SAMPLE"

    bwa mem -t $THREADS \
        -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tLB:lib1" \
        $REF $R1 $R2 | \
    samblaster --excludeDups --addMateTags --maxSplitCount 2 --minNonOverlap 20 | \
    samtools view -bS - | \
    samtools sort -@ $THREADS -o ${OUT_DIR}/${SAMPLE}.bam

    samtools index ${OUT_DIR}/${SAMPLE}.bam
done
````
