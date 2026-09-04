## run boc_loop.sh to get BOC values, make list of breadth4x.tsv file and run filter-merge-boc.sh script, this will  filter BOC values and return merged_tsv file. run unique.sh on the merged tsv This file now can be processed in R.

##  compute BOC for window to calculate SV
````
#!/bin/bash

# ============================================================
# Breadth-of-coverage SV screening
#
# Calculates, for every 10-kb window:
#   1. Number of bases with >=4x coverage
#   2. Breadth of coverage = bases >=4x / window size
#
# Usage:
#   ./breadth_coverage.sh sample.bam reference.fa output_dir
# ============================================================

set -euo pipefail

BAM=$1
REF=$2
OUT=$3

SAMPLE=$(basename "$BAM" .bam)

mkdir -p "$OUT"

# ------------------------------------------------------------
# 1. Reference index
# ------------------------------------------------------------

if [ ! -f "${REF}.fai" ]; then
    samtools faidx "$REF"
fi

# ------------------------------------------------------------
# 2. Genome sizes
# ------------------------------------------------------------

cut -f1,2 "${REF}.fai" > "$OUT/genome.txt"

# ------------------------------------------------------------
# 3. Make 10-kb windows
# ------------------------------------------------------------

bedtools makewindows \
    -g "$OUT/genome.txt" \
    -w 10000 \
    > "$OUT/windows.bed"
    # ------------------------------------------------------------
# 4. Calculate per-base coverage
# ------------------------------------------------------------

mosdepth \
    --threads 16 \
    --by "$OUT/windows.bed" \
    --thresholds 4 \
    "$OUT/$SAMPLE" \
    "$BAM"
# ------------------------------------------------------------
# 5. Calculate breadth >=4x for each 10-kb window
# ------------------------------------------------------------

samtools depth \
    -aa \
    -b "$OUT/windows.bed" \
    "$BAM" \
    > "$OUT/${SAMPLE}.depth.txt"

# ------------------------------------------------------------
# 6. Aggregate per-base depth into windows
# ------------------------------------------------------------

awk '
BEGIN {
    OFS="\t"
    print "Sample","Contig","Start","End","Window_bp","Bases_ge_4x","Breadth_4x"
}

{
    contig=$1
    pos=$2
    depth=$3

    # Identify 10-kb window
    start=int((pos-1)/10000)*10000
    end=start+10000

    key=contig ":" start ":" end

    total[key]++

    if (depth >= 4)
        covered[key]++
}

END {
    for (key in total) {

        split(key,a,":")

        contig=a[1]
        start=a[2]
        end=a[3]

        n=total[key]
        c=covered[key]+0

        breadth=c/n

        print "'"$SAMPLE"'",contig,start,end,n,c,breadth
    }
}
' "$OUT/${SAMPLE}.depth.txt" \
| sort -k2,2 -k3,3n \
> "$OUT/${SAMPLE}.breadth4x.tsv"
echo ""
echo "Output:"
echo "$OUT/${SAMPLE}.breadth4x.tsv"

echo ""
echo "First 10 windows:"
head "$OUT/${SAMPLE}.breadth4x.tsv"
````

## loop throguh boc.sh scrupt
````
for bam in /mnt/etoposide-marigold/SORTED_BAMS/*.sorted.bam; do
    sample=$(basename "$bam" .sorted.bam)
    echo "Starting $sample"

    ./boc.sh \
        "$bam" \
        /mnt/etoposide-marigold/REF/OG99_ragtag_min10K.fasta \
        "$sample" \
        > "${sample}.log" 2>&1 &
done

wait
echo "All coverage jobs finished."
````



## filter merge BOC output
./script.sh listoftsv merged-output.tsv
````
#!/bin/bash

LIST="$1"
OUTPUT="$2"

echo -e "Sample\tContig\tStart\tEnd\tSize_bp\tWindows_merged" > "$OUTPUT"

while read -r TSV; do

    awk -F'\t' '
    BEGIN { OFS="\t" }

    NR==1 { next }

    # Breadth_4x < 0.05
    $7 >= 0.05 { next }

    {
        if (!active) {
            sample=$1
            contig=$2
            region_start=$3
            region_end=$4
            n=1
            active=1
            next
        }

        # Adjacent qualifying window
        if ($1 == sample &&
            $2 == contig &&
            $3 == region_end) {

            region_end=$4
            n++

        } else {

            print sample,contig,
                  region_start,region_end,
                  region_end-region_start,n

            sample=$1
            contig=$2
            region_start=$3
            region_end=$4
            n=1
        }
    }

    END {
        if (active) {
            print sample,contig,
                  region_start,region_end,
                  region_end-region_start,n
        }
    }
    ' "$TSV" >> "$OUTPUT"

done < "$LIST"

````

## unique.sh
````
#!/bin/bash

INPUT="$1"
OUTPUT="$2"

# Create proper BED file:
# Contig  Start  End  Sample
awk -F'\t' '
NR > 1 {
    print $2, $3, $4, $1
}' OFS="\t" "$INPUT" > all_sv.bed


# Find SVs overlapping an SV from another sample
bedtools intersect \
    -a all_sv.bed \
    -b all_sv.bed \
    -wa -wb |
awk -F'\t' '
BEGIN { OFS="\t" }

{
    contig1=$1
    start1=$2
    end1=$3
    sample1=$4

    contig2=$5
    start2=$6
    end2=$7
    sample2=$8

    # Same sample = not evidence of sharing
    if (sample1 == sample2)
        next

    # Different contigs cannot overlap
    if (contig1 != contig2)
        next

    overlap_start = (start1 > start2 ? start1 : start2)
    overlap_end   = (end1 < end2 ? end1 : end2)

    if (overlap_end > overlap_start) {

        overlap = overlap_end - overlap_start
        len1 = end1 - start1
        len2 = end2 - start2

        # 50% reciprocal overlap
        if (overlap / len1 >= 0.50 &&
            overlap / len2 >= 0.50) {

            print sample1, contig1, start1, end1
        }
    }
}
' > shared_sv.tmp


# Identify SVs that are NOT shared with another sample
awk -F'\t' '
BEGIN { OFS="\t" }

FNR==NR {
    shared[$1 FS $2 FS $3 FS $4] = 1
    next
}

NR==1 {
    print $0, "Unique"
    next
}

{
    key=$1 FS $2 FS $3 FS $4

    if (key in shared)
        print $0, "FALSE"
    else
        print $0, "TRUE"
}
' shared_sv.tmp "$INPUT" > "$OUTPUT"


rm -f all_sv.bed shared_sv.tmp


````
