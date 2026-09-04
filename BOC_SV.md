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

## Now filter the BOC table based on BOC cutoff

````
awk -F'\t' 'BEGIN{OFS="\t"}
NR==1 || $7 < 0.05
' MAR_Et3_01_S84.sorted.breadth4x.tsv | less -S
````

## filter BOC output
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
