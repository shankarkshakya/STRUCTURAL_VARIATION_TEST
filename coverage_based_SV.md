````
once the bam is produced from the step1 mapping. this can be used further with bedtools coverage to identify SV based on coverage.
The following script will need a 100jb_windows.bed file. This can be generated from ref.fa file.

samtools faidx Col-CEN_v1.2.fasta
 cut -f1,2 Col-CEN_v1.2.fasta.fai > genome.txt
bedtools makewindows -g genome.txt -w 100000 > 100kb_windows.bed
./bedtools_cov.sh

````

````
#!/bin/bash
set -euo pipefail

# =========================
# CONFIG
# =========================
BAM_DIR="bam"
WINDOWS="100kb_windows.bed"

COV_DIR="coverage"
NORM_DIR="normalized"
RESULTS_DIR="results"

mkdir -p $COV_DIR $NORM_DIR $RESULTS_DIR

# =========================
# STEP 1: COVERAGE (auto-detect BAMs)
# =========================
echo "Calculating coverage..."

for BAM in $BAM_DIR/*.bam
do
    SAMPLE=$(basename "$BAM" .bam)
    echo "Processing $SAMPLE"

    # Coverage using properly paired reads (streaming, no temp BAM)
    samtools view -b -f 2 -F 4 "$BAM" | \
    bedtools coverage \
        -a "$WINDOWS" \
        -b stdin \
        -counts \
    > $COV_DIR/${SAMPLE}.cov

done

# =========================
# STEP 2: NORMALIZE (RPM)
# =========================
echo "Normalizing coverage..."

for BAM in $BAM_DIR/*.bam
do
    SAMPLE=$(basename "$BAM" .bam)

    TOTAL_READS=$(samtools view -c -f 2 -F 4 "$BAM")

    awk -v total=$TOTAL_READS '{
        rpm = ($4 / total) * 1000000;
        print $1"\t"$2"\t"$3"\t"rpm
    }' $COV_DIR/${SAMPLE}.cov > $NORM_DIR/${SAMPLE}.rpm

done


# =========================
# STEP 3: MEDIAN PER WINDOW
# =========================
echo "Computing median across samples..."

paste $NORM_DIR/*.rpm | \
awk '{
    chr=$1; start=$2; end=$3;

    n=0;
    for(i=4; i<=NF; i+=4){
        vals[n++]=$i;
    }

    asort(vals);

    if(n%2==1){
        median=vals[int(n/2)];
    } else {
        median=(vals[n/2-1]+vals[n/2])/2;
    }

    print chr"\t"start"\t"end"\t"median;
}' > $RESULTS_DIR/median_windows.txt


# =========================
# STEP 4: DOSAGE + SV CALLING
# =========================
echo "Calling structural variants..."

for BAM in $BAM_DIR/*.bam
do
    SAMPLE=$(basename "$BAM" .bam)

    paste $NORM_DIR/${SAMPLE}.rpm $RESULTS_DIR/median_windows.txt | \
    awk '{
        chr=$1; start=$2; end=$3;
        rpm=$4;
        median=$8;

        if(median==0){
            dosage="NA";
            status="NA";
        } else {
            dosage=(rpm/median)*2;

            if(dosage >= 4){
                status="HOMO_DUP";
            } else if(dosage >= 3){
                status="HET_DUP";
            } else if(dosage <= 0){
                status="HOMO_DEL";
            } else if(dosage <= 1){
                status="HET_DEL";
            } else {
                status="NORMAL";
            }
        }

        print chr"\t"start"\t"end"\t"rpm"\t"median"\t"dosage"\t"status;
    }' > $RESULTS_DIR/${SAMPLE}.sv_calls.txt

done

echo "Pipeline complete!"
````
