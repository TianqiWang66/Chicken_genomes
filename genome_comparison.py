#!/bin/bash
REF="hapBfinal.fa"     
QRY="hapAfinal.fa"      
PREFIX="out"                  
NUCMER_OPTS="--maxmatch -l 100 -c 500"
FILTER_OPTS="-1 -l 4000"             
COORDS_OPTS="-THrd"               

nucmer $NUCMER_OPTS $REF $QRY -prefix $PREFIX
delta-filter $FILTER_OPTS ${PREFIX}.delta > ${PREFIX}.filtered.delta
show-coords $COORDS_OPTS ${PREFIX}.filtered.delta > ${PREFIX}.coords
syri -c ${PREFIX}.coords -d ${PREFIX}.filtered.delta -r $REF -q $QRY --prefix $PREFIX
plotsr \
    --sr ${PREFIX}syri.out \
    --genomes <(echo -e "${REF}\n${QRY}") \   # 临时构建 genomes.txt
    -f 20 -W 10 -H 35 \
    -o ${PREFIX}_plot.pdf
