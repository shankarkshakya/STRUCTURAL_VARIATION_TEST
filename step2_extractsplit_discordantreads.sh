````
#!/bin/bash
set -euo pipefail

BAM_DIR="bam"
OUT_DIR="lumpy"

mkdir -p $OUT_DIR/tmp

for BAM in $BAM_DIR/*.bam
do
    SAMPLE=$(basename $BAM .bam)

    echo "Processing $SAMPLE"

    # Discordant reads
    samtools view -b -F 1294 $BAM > $OUT_DIR/tmp/${SAMPLE}.discordants.bam

    # Split reads
    samtools view -h $BAM | \
        extractSplitReads_BwaMem -i stdin | \
        samtools view -Sb - > $OUT_DIR/tmp/${SAMPLE}.splitters.bam

    # Sort + index
    samtools sort -o $OUT_DIR/tmp/${SAMPLE}.discordants.sorted.bam $OUT_DIR/tmp/${SAMPLE}.discordants.bam
    samtools sort -o $OUT_DIR/tmp/${SAMPLE}.splitters.sorted.bam $OUT_DIR/tmp/${SAMPLE}.splitters.bam

    samtools index $OUT_DIR/tmp/${SAMPLE}.discordants.sorted.bam
    samtools index $OUT_DIR/tmp/${SAMPLE}.splitters.sorted.bam
done
````
