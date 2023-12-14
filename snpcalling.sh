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

bcftools mpileup -Ou -q 20 -Q 20  -a FORMAT/DP -f $REFERENCE $INPUTFOLDER/$BAM | \
  bcftools call -mO b -g $MINDEP > $OUTPUTFOLDER/$SAMPLE.unflt.bcf

bcftools index $OUTPUTFOLDER/$SAMPLE.unflt.bcf
