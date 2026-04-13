#!/bin/bash
set -euo pipefail

VCF_DIR="lumpy/results"

for VCF in $VCF_DIR/*.vcf
do
    SAMPLE=$(basename $VCF .vcf)

    bcftools view -i 'SVTYPE!="BND" && INFO/SU>=5' $VCF \
    > $VCF_DIR/${SAMPLE}.filtered.vcf
done
