#!/bin/bash

# Will run trim_galore and fastqc on both raw and trimmed files
# Expects raw files are on rawdata folder, and only 2 raw data files exist for each sample (F & R)


SAMPLE=$1
RAWFOLDER=/home/rawdata/Lupinus/Resequencing
OUTPUTFOLDER=/home/ivieira/MappingBias/trimmed


FREADS=$RAWFOLDER/$1/$1.1.fq.gz
RREADS=$RAWFOLDER/$1/$1.2.fq.gz

[ ! -e $FREADS ] && echo "ERROR: cant find F raw file: $FREADS"  && exit
[ ! -e $RREADS ] && echo "ERROR: cant find R raw file: $RREADS"  && exit
[ -e $OUTPUTFOLDER/$SAMPLE.1_val_1.fq.gz ] && echo "ERROR: output file F exists (check $OUTPUTFOLDER)"  && exit
[ -e $OUTPUTFOLDER/$SAMPLE.2_val_2.fq.gz ] && echo "ERROR: output file F exists (check $OUTPUTFOLDER)"  && exit


module load trimgalore/0.6.10
module load fastqc/0.12.1


fastqc $FREADS -o /home/ivieira/MappingBias/rawdata
fastqc $RREADS -o /home/ivieira/MappingBias/rawdata

trim_galore -q 20 --length 50 --fastqc --output_dir $OUTPUTFOLDER \
            --paired $FREADS $RREADS