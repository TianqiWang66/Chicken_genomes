#!/bin/bash
GENOME="/path/to/genome.fasta"            
HAPMER1="/path/to/genomeZ1.hapmer.meryl"      
HAPMER2="/path/to/genomeZ2.hapmer.meryl"       
OUT="out"                                    
BUSCO_LINEAGE="busco_downloads/lineages/aves_odb10"
THREADS=10

source ~/miniconda3/etc/profile.d/conda.sh     #

conda activate merqury
$CONDA_PREFIX/share/merqury/trio/phase_block.sh \
    "$GENOME" "$HAPMER1" "$HAPMER2" "$OUT"
$CONDA_PREFIX/share/merqury/trio/switch_error.sh \
    "${OUT}.sort.bed" "$OUT" 1 100 F
conda deactivate

conda activate busco
busco -i "$GENOME" -c $THREADS -o final63 -m geno -l "$BUSCO_LINEAGE"
quast "$GENOME" -o quast_output
conda deactivate
