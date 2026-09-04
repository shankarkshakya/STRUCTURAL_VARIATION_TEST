## first calucalte normalize cov for window

````
#!/bin/bash

# ============================================================
# Coverage-based SV detection
# 10-kb windows + normalized coverage
#
# Usage:
#   ./coverage_10kb.sh sample.bam reference.fa output_dir
#
# Output:
#   sample.regions.bed.gz
#   sample.coverage.tsv
#   windows.bed
#   genome.txt
#
# Columns:
#   Sample
#   Contig
#   Start
#   End
#   Mean_Coverage
#   Normalized_Coverage
# ============================================================

set -euo pipefail

BAM=$1
REF=$2
OUT=$3

SAMPLE=$(basename "$BAM" .bam)

mkdir -p "$OUT"

echo "=========================================="
echo "Sample:    $SAMPLE"
echo "BAM:       $BAM"
echo "Reference: $REF"
echo "Output:    $OUT"
echo "=========================================="

# ------------------------------------------------------------
# 1. Reference index
# ------------------------------------------------------------

if [ ! -f "${REF}.fai" ]; then
    echo "Creating reference index..."
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

echo "Created $(wc -l < "$OUT/windows.bed") windows."

# ------------------------------------------------------------
# 4. Calculate coverage
# ------------------------------------------------------------

mosdepth \
    --threads 16 \
    --by "$OUT/windows.bed" \
    "$OUT/$SAMPLE" \
    "$BAM"

# ------------------------------------------------------------
# 5. Calculate genome-wide median coverage
# ------------------------------------------------------------
#
# mosdepth regions file:
#   column 1 = contig
#   column 2 = start
#   column 3 = end
#   column 4 = mean coverage
#
# We calculate the median of all window mean coverages.
# ------------------------------------------------------------
MEDIAN=$(zcat "$OUT/${SAMPLE}.regions.bed.gz" \
    | awk '{print $4}' \
    | sort -n \
    | awk '
        {
            a[NR]=$1
        }
        END {
            if (NR == 0) {
                exit 1
            }
            if (NR % 2 == 1) {
                print a[(NR+1)/2]
            } else {
                print (a[NR/2] + a[NR/2+1]) / 2
            }
        }
    ')

echo ""
echo "Genome-wide median coverage: $MEDIAN"

# ------------------------------------------------------------
# 6. Create coverage + normalized coverage table
# ------------------------------------------------------------

zcat "$OUT/${SAMPLE}.regions.bed.gz" \
    | awk -v sample="$SAMPLE" -v median="$MEDIAN" '
        BEGIN {
            OFS="\t"
            print "Sample","Contig","Start","End","Mean_Coverage","Normalized_Coverage"
        }
        {
            normalized = $4 / median
            print sample,$1,$2,$3,$4,normalized
        }
    ' \
    > "$OUT/${SAMPLE}.coverage.tsv"

# ------------------------------------------------------------
# 7. Display results
# ------------------------------------------------------------

echo ""
echo "=========================================="
echo "Coverage results"
echo "=========================================="

head -n 11 "$OUT/${SAMPLE}.coverage.tsv"

echo ""
echo "Output:"
echo "$OUT/${SAMPLE}.coverage.tsv"

echo ""
echo "DONE"

````
## next step is to take tsv file from above and identify SV

````
#!/bin/bash

# ============================================================
# Identify SVs from normalized coverage
#
# Input TSV columns:
# Sample  Contig  Start  End  Mean_Coverage  Normalized_Coverage
#
# Usage:
#   ./identify_coverage_SV.sh sample.coverage.tsv
#
# Default:
#   DEL: normalized coverage < 0.65
#   DUP: normalized coverage > 1.35
#   Minimum SV size: 20 kb
#
# Output:
#   sample.DEL.SV.tsv
#   sample.DUP.SV.tsv
#   sample.coverage.SV.tsv
# ============================================================

set -euo pipefail

INPUT=$1

# ------------------------------------------------------------
# Parameters
# ------------------------------------------------------------

DEL_THRESHOLD=0.5
DUP_THRESHOLD=1.35
MIN_SV_SIZE=50000

PREFIX="${INPUT%.tsv}"

# ------------------------------------------------------------
# Check input
# ------------------------------------------------------------

if [ ! -f "$INPUT" ]; then
    echo "ERROR: Input file not found: $INPUT"
    exit 1
fi

# ------------------------------------------------------------
# Identify abnormal windows
# ------------------------------------------------------------


awk -v del="$DEL_THRESHOLD" \
    -v dup="$DUP_THRESHOLD" '
BEGIN {
    FS=OFS="\t"
}

NR == 1 {
    print > "'"$PREFIX"'.abnormal.windows.tsv"
    next
}

{
    norm=$6

    if (norm < del) {
        print $0, "DEL" >> "'"$PREFIX"'.abnormal.windows.tsv"
    }
    else if (norm > dup) {
        print $0, "DUP" >> "'"$PREFIX"'.abnormal.windows.tsv"
    }
}
' "$INPUT"

# ------------------------------------------------------------
# Sort abnormal windows
# ------------------------------------------------------------

{
    head -n 1 "$PREFIX.abnormal.windows.tsv"

    tail -n +2 "$PREFIX.abnormal.windows.tsv" \
        | sort -k2,2 -k3,3n

} > "$PREFIX.abnormal.sorted.tsv"

# ------------------------------------------------------------
# Merge consecutive/adjacent abnormal windows
# ------------------------------------------------------------

awk '
BEGIN {
    FS=OFS="\t"
}

NR == 1 {
    print "Sample","SV_Type","Contig","Start","End",
          "Size_bp","Windows","Mean_Normalized_Coverage"
    next
}

{
    sample=$1
    contig=$2
    start=$3
    end=$4
    norm=$6
    type=$7

    # Start first SV
    if (!active) {
        prev_contig=contig
        prev_type=type
        sv_start=start
        sv_end=end
        sum_norm=norm
        n=1
        active=1
        next
    }

    # Same chromosome + same SV type + adjacent windows
    if (contig == prev_contig &&
        type == prev_type &&
        start <= sv_end) {

        sv_end=end
        sum_norm += norm
        n++

    } else {

        size=sv_end-sv_start
        mean_norm=sum_norm/n

        if (size >= 20000) {
            print sample,prev_type,prev_contig,
                  sv_start,sv_end,size,n,mean_norm
        }

        prev_contig=contig
        prev_type=type
        sv_start=start
        sv_end=end
        sum_norm=norm
        n=1
    }
}

END {
    if (active) {

        size=sv_end-sv_start
        mean_norm=sum_norm/n

        if (size >= 20000) {
            print sample,prev_type,prev_contig,
                  sv_start,sv_end,size,n,mean_norm
        }
    }
}
' "$PREFIX.abnormal.sorted.tsv" \
> "$PREFIX.coverage.SV.tsv"

# ------------------------------------------------------------
# Separate DEL and DUP
# ------------------------------------------------------------

head -n 1 "$PREFIX.coverage.SV.tsv" \
> "$PREFIX.DEL.SV.tsv"

awk -F'\t' '$2=="DEL"' "$PREFIX.coverage.SV.tsv" \
>> "$PREFIX.DEL.SV.tsv"


head -n 1 "$PREFIX.coverage.SV.tsv" \
> "$PREFIX.DUP.SV.tsv"

awk -F'\t' '$2=="DUP"' "$PREFIX.coverage.SV.tsv" \
>> "$PREFIX.DUP.SV.tsv"

# ------------------------------------------------------------
# Summary
# ------------------------------------------------------------


echo ""
echo "=========================================="
echo "Coverage-based SV detection"
echo "=========================================="

echo "Input:           $INPUT"
echo "DEL threshold:   < $DEL_THRESHOLD"
echo "DUP threshold:   > $DUP_THRESHOLD"
echo "Minimum SV size: $MIN_SV_SIZE bp"

echo ""
echo "Candidate deletions:"
awk 'NR > 1 {print}' "$PREFIX.DEL.SV.tsv" | wc -l

echo "Candidate duplications:"
awk 'NR > 1 {print}' "$PREFIX.DUP.SV.tsv" | wc -l

echo ""
echo "Output:"
echo "  $PREFIX.DEL.SV.tsv"
echo "  $PREFIX.DUP.SV.tsv"
echo "  $PREFIX.coverage.SV.tsv"

echo ""
echo "DONE"







````
