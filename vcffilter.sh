#!/bin/bash

# will make very simple filter of SNPs only
# outputs vcf

SAMPLE=$1
VCF=/home/ivieira/MappingBias/vcfs/$SAMPLE.unflt.bcf
BAM=/home/ivieira/MappingBias/bam/$SAMPLE.realigned.bam
OUTFOLDER=/home/ivieira/MappingBias/vcfs

module load samtools/1.17
module load bcftools/1.17

MAXDEP=`samtools mpileup -q 20 -Q 20  $BAM  | cut -f 4 \
    | perl -ne ' $tota+=$_;$n++;print 2.5*$tota/$n,"\n";' | tail -n 1 | cut -d'.' -f1`

echo $MAXDEP

bcftools filter -O v  -ssimpleFlt1 -e 'QUAL <= 10 || INFO/DP > '$MAXDEP' || FMT/DP<8 || ( GT="0/1" && DP4[0]+DP4[1] < 2  ) ' $VCF  | grep -v simpleFlt1 | grep -v INDEL > $OUTFOLDER/$SAMPLE.flt1.snps.vcf
