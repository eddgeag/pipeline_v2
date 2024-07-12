#!/bin/bash
echo "compute faidx" 
faidx ~/datos_exomas/datos_gatk/referencia/ensembl/referencia_filtered.fasta -i chromsizes >  chrom.sizes
echo "encontrando tamanos de cromosomas canonicos" 
awk '/^chr[0-9XY]*\t/ {printf("%s\t0\t%s\n",$1,$2);}' ~/datos_exomas/datos_gatk/referencia/ensembl/referencia_filtered.fasta.fai > bed_split.bed
echo "cortar el bam" 
samtools view -L bed_split.bed -o out.bam $1
echo "sortear el bam a 16 gb" 
samtools sort out.bam -@ 32 -o myfile_sorted.bam
echo "sortear los chrom.sizes" 
cat chrom.sizes|sort -V > sizes.genome.sort
echo "compuar cobertura" 
bedtools coverage -a ~/datos_exomas/coverage/hg38/xgen-exome.bed  -b myfile_sorted.bam -g sizes.genome.sort -sorted -hist 
