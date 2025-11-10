#obtain reads with a quality score above 20

count=$(samtools view -F 4 "$input_bam" | awk '$5 >= 20' | wc -l)


#total number of reads that passed quality filters and successfully aligned to the reference genome and the total number of aligned bases with a quality score of at least Q20

gatk CollectAlignmentSummaryMetrics \ R=/home/ivieira/MappingBias/genome/Lmutabilis.mollA12.v1.0.asm.p.fa \
I=/home/ivieira/MappingBias/bam/$SAMPLE.realigned.bam \ O=/home/ivieira/MappingBias/bam/$SAMPLE.txt


#depth,  number of sites that meet the depth criteria, and the average coverage at those sites.

samtools depth -q 20 -Q 20 $SAMPLE.realigned.bam > $SAMPLE.depth.txt
awk '$3 >= 8 && $3 <= *max depth* { sum += $3; count++ } END { if (count) print "Count:", count, "Mean:", sum / count }' Lup2108.depth.txt

#heterozygity and missingness

vcf2fas -reference /home/ivieira/MappingBias/genome/Lmutabilis.mollA12.v1
.0.asm.p.fa -vcfs /home/ivieira/MappingBias/bam/samples.txt -gf GT - strictIUPAC 1
for i in {1..24} do
cleanfasta2 -infile /home/ivieira/MappingBias/bam/scaffold_$i.fas -
maxMissing 100 -verbose 1 -format stats -gff /home/ivieira/MappingBias/ genome/Lmutabilis.mollA12.v1.0.asm.p.EDTA.TEanno.gff3 -scaffold scaffold_$i -include 0 -feature all -outfile clean_scaffold$i.fas
done

