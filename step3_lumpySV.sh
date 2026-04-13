````
#!/bin/bash
set -euo pipefail

BAM_DIR="bam"
OUT_DIR="lumpy"

mkdir -p $OUT_DIR/results

for BAM in $BAM_DIR/*.bam
do
    SAMPLE=$(basename $BAM .bam)

    lumpyexpress \
        -B $BAM \
        -S $OUT_DIR/tmp/${SAMPLE}.splitters.sorted.bam \
        -D $OUT_DIR/tmp/${SAMPLE}.discordants.sorted.bam \
        -o $OUT_DIR/results/${SAMPLE}.vcf
done
````
