#!/bin/bash
SAMPLE=$1
NT=$2

[ $# != 2 ] && echo "Need 2 args: sampleName NT " && exit


# Will map trimmed data to reference. Includes gatk3 indel realignment step 
# Expects input and output folders as below

INPUTFOLDER=/home/ivieira/MappingBias/trimmed
OUTPUTFOLDER=/home/ivieira/MappingBias/mapk15
REFERENCE=/home/ivieira/MappingBias/genome/Lmutabilis.mollA12.v1.0.asm.p.fa
FREADS=$INPUTFOLDER/$SAMPLE.1_val_1.fq.gz
RREADS=$INPUTFOLDER/$SAMPLE.2_val_2.fq.gz

echo "Starting mapping for $SAMPLE at" `date`


[ ! -e $FREADS ] && echo "ERROR: cant find F trimmed file: $FREADS"  && exit
[ ! -e $RREADS ] && echo "ERROR: cant find R trimmed file: $RREADS"  && exit
[ ! -e $REFERENCE ] && echo "ERROR: cant find reference file: $REFERENCE"  && return 1
[ -e $OUTPUTFOLDER/$SAMPLE.realigned.bam ] && echo "ERROR: output file $SAMPLE.realigned.bam exists (check $OUTPUTFOLDER)"  && exit


RFOLDER=`mktemp -d tempXXXXXXX`
cd $RFOLDER

# =====================================================================================
module load bwa/0.7.17
module load gatk/4.4.0.0
module load samtools/1.17
conda activate gatk

#======================================================================================
# MAPPING

bwa mem -M -t $NT $REFERENCE $FREADS $RREADS | \
samtools fixmate -u -m - - | \
samtools sort -u -@ $NT - | \
samtools markdup -@ $NT --reference $REFERENCE - $SAMPLE.bam
samtools index $SAMPLE.bam


# add read groups
gatk AddOrReplaceReadGroups \
    -I $SAMPLE.bam \
    -O $SAMPLE.RGOK.bam \
    -ID 1 -LB 1 -PL illumina -PU 1 -SM "$SAMPLE"
samtools index $SAMPLE.RGOK.bam
samtools flagstat $SAMPLE.RGOK.bam > $OUTPUTFOLDER/$SAMPLE.stats  



#======================================================================================
## GATK INDEL REALIGNMENT - requires gatk 3

conda deactivate
module purge
module load gatk/3.8.1
module load samtools

gatk3 -T RealignerTargetCreator \
    -R $REFERENCE -I $SAMPLE.RGOK.bam \
    -o $SAMPLE.realigner.intervals

gatk3 -T IndelRealigner -maxReads 100000 \
    -R $REFERENCE -I $SAMPLE.RGOK.bam \
    -targetIntervals $SAMPLE.realigner.intervals \
    -o $SAMPLE.realigned.bam
samtools index $SAMPLE.realigned.bam

mv $SAMPLE.realigned.bam $OUTPUTFOLDER
mv $SAMPLE.realigned.bam.bai $OUTPUTFOLDER
#;cd ..
#rm -r $RFOLDER
