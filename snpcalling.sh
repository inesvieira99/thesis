#!/bin/bash


module load bcftools/1.17
module load samtools/1.17
 
SAMPLE=$1
BAM=$SAMPLE.realigned.bam

REFERENCE=/home/ivieira/MappingBias/genome/Lmutabilis.mollA12.v1.0.asm.p.fa
INPUTFOLDER=/home/ivieira/MappingBias/bam
OUTPUTFOLDER=/home/ivieira/MappingBias/vcfs


MINDEP=8

#======================================================================================

MAXDEP=`samtools mpileup -q 20 -Q 20 $INPUTFOLDER/$BAM | cut -f 4 \
| perl -ne ' $tota+=$_;$n++;print 2.5*$tota/$n,"\n";' | tail -n 1 | cut
-d'.' -f1`

bcftools mpileup -Ov -q 20 -Q 20 -a FORMAT/DP -f $REFERENCE $INPUTFOLDER/
$BAM | \
bcftools call -mO b -g $MINDEP,$MAXDEP | \
bcftools filter -O v -g3 -G10 -s filter311023 -e "QUAL < 20 || FMT/DP >
$MAXDEP || FMT/DP < 8 || GT='0/1'" | \
grep -v filter311023 > $OUTPUTFOLDER/$SAMPLE.flt.vcf


#bcftools mpileup -Ou -q 20 -Q 20  -a FORMAT/DP -f $REFERENCE $INPUTFOLDER/$BAM | \
  #bcftools call -mO b -g $MINDEP > $OUTPUTFOLDER/$SAMPLE.unflt.bcf

#bcftools index $OUTPUTFOLDER/$SAMPLE.unflt.bcf
