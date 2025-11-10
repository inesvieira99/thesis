
#separating the mapped and unmapped reads

samtools view -f 0x4 -F 0x800 Lup2164.realigned.bam > unmapped2164.txt 
samtools view -F 0x4 -F 0x800 Lup2164.realigned.bam > mapped2164.txt

SAMPLE=$1
# Shuffle the lines and keep the first 10,000 lines
shuf mapped$SAMPLE.txt | head -n 10000 > shuffled_and_sampled$SAMPLE.txt
# Extract the 10th column and create a new file
awk '{print $10}' shuffled_and_sampled$SAMPLE.txt > $SAMPLE.
extracted_column.txt
# Combine line numbers and contents from the 10th column
awk '{print ">line" NR "\n" $0}' $SAMPLE.extracted_column.txt > $SAMPLE. fasta.fa

#Blast
for SAMPLE in "$@"; do
echo "Running blast for $SAMPLE"
blastn -query "/home/ivieira/MappingBias/bam/${SAMPLE}.fasta.fa" -db "/ home/bnevado/software/blast/db/nt" -max_target_seqs 5 -outfmt "6 ' qseqid qaccver saccver pident length staxid ssciname sskingdom" > "
Lup${SAMPLE}.blast.txt"
echo "Blast for $SAMPLE is done"
done