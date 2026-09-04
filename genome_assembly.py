###Genome assembly
#!/bin/bash
yak count -k31 -b37 -t16 -o pat.yak paternal.fq.gz
yak count -k31 -b37 -t16 -o mat.yak maternal.fq.gz
hifiasm -o genome.asm -t32 --ul ul.fq.gz -1 pat.yak -2 mat.yak HiFi-reads.fq.gz


THREADS=32
ENZYME=GATC           
CHROM_NUM=39         
bwa index contig.fasta
bwa mem -5SP -t $THREADS contig.fasta R1.fq.gz R2.fq.gz | samtools view -bS - > hic.bam
samtools sort -n -@ $THREADS hic.bam -o hic.sortn.bam
allhic extract hic.sortn.bam contig.fasta --RE $ENZYME
allhic optimize contig.clm $CHROM_NUM
allhic build contig.fasta contig.agp
echo "ALLHiC done. Check contig.chromosome.fasta and contig.agp"
