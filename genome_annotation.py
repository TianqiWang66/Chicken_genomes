#!/bin/bash
GENOME=genome.fasta                        
PREFIX=mydb                 
trf $GENOME 2 7 7 80 10 50 2000 -d -h
misa.pl $GENOME
BuildDatabase -name $PREFIX $GENOME
RepeatModeler -database $PREFIX -threads $THREADS -LTRStruct
LTR_FINDER -threads $THREADS -harvest_out -size 1000000 -time 300 $GENOME

gt suffixerator -db $GENOME -indexname genome_idx -tis -suf -lcp -des -ssp -sds -dna
gt ltrharvest -index genome_idx \
    -minlenltr 100 -maxlenltr 7000 -mintsd 4 -maxtsd 6 \
    -motif TGCA -motifmis 1 -similar 85 -vic 10 -seed 20 \
    -seqids yes -out ltrharvest.out

LTR_retriever -genome $GENOME \
    -inharvest ltrharvest.out \
    -infinder ${GENOME}.finder.combine \   
    -threads $THREADS -noanno

cat ${PREFIX}-families.fa LTR_retriever_result.fa > denovo_lib_raw.fa
teclass -i denovo_lib_raw.fa -o denovo_lib_classified.fa
RepeatMasker -lib denovo_lib_classified.fa \
    -noLowSimple -pvalue 0.0001 -threads $THREADS $GENOME
