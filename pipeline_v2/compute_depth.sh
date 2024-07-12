#!/bin/bash
faidx $1 -i chromsizes >  chrom.sizes
awk '/^chr[0-9XY]*\t/ {printf("%s\t0\t%s\n",$1,$2);}' $2 > bed_split.bed
samtools view -L bed_split.bed -o out.bam $3
cat chrom.sizes|sort -V > sizes.genome.sort
bedtools coverage -a ~/datos_exomas/coverage/hg38/xgen-exome.bed  -b out.bam -g sizes.genome.sort -sorted -hist 
