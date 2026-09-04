#!/bin/bash
REF=genome.fa
BAM=dedup_bam
OUT=sv_results
mkdir -p $OUT/{manta,delly,smoove,merge}
# Manta
for b in $BAM/*.dedup.bam; do s=$(basename $b .dedup.bam); configManta.py --bam $b --referenceFasta $REF --runDir $OUT/manta/$s && python $OUT/manta/$s/runWorkflow.py && cp $OUT/manta/$s/results/variants/diploidSV.vcf.gz $OUT/manta/${s}.vcf.gz; done
# Delly
for b in $BAM/*.dedup.bam; do s=$(basename $b .dedup.bam); delly call -g $REF -o $OUT/delly/${s}.bcf $b && bcftools view $OUT/delly/${s}.bcf > $OUT/delly/${s}.vcf; done
# Smoove
for b in $BAM/*.dedup.bam; do s=$(basename $b .dedup.bam); smoove call --outdir $OUT/smoove --name $s --fasta $REF --genotype $b; done
# Per-sample merge
for s in $(ls $BAM/*.dedup.bam | sed 's/.*\///;s/\.dedup\.bam//'); do
    echo -e "$OUT/manta/${s}.vcf.gz\n$OUT/smoove/${s}.smoove.genotyped.vcf\n$OUT/delly/${s}.vcf" > $OUT/merge/list_$s.txt
    SURVIVOR merge $OUT/merge/list_$s.txt 1000 2 1 0 0 50 $OUT/merge/${s}_merged.vcf
done
# Population-level merge
ls $OUT/merge/*_merged.vcf > $OUT/merge/all.txt
SURVIVOR merge $OUT/merge/all.txt 1000 1 1 0 0 50 $OUT/merge/population.vcf
