
#adds polimorfism and divergence
python3 poli.py --input scaffold_1.1.fa --divergence $div --polimorfism $poli --outdiv ${div}_div.fa --outpoli ${div}_${poli}_div.fa

#combines the diploid sequences generated in the previous step
python3 iupac2.py --input1 ${div}_div.fa --input2 ${div}_${poli}_div.fa -- output ${div}_genome.fa


#add indels
/home/bnevado/temp/ines/add_indels/v2/add_indels2 -infile ${div}_genome.fa -outfile out${div}_${indel}_sim.fa -pIndel $indel -normMean 0 -normVar
7.7

#converts IUPAC codes into diploid FASTA sequences
diploid2randomphase -in ~/MappingBias/sim_indels/171024/out${div}_${indel} _sim.fa -out ${div}_${indel}_sim.fa
 
  
#art simulation
art_illumina -p -sam -i ${div}_${indel}_sim.fa -l 150 -f 8 -m 400 -s 10 -o ${div}_${indel}indel -ss HS25
 