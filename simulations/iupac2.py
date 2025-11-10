mport sys
import getopt
from Bio import SeqIO, Seq

# IUPAC ambiguity code for consensus sequence generation
iupac_code = {
    ('A', 'A'): 'A', ('C', 'C'): 'C', ('G', 'G'): 'G', ('T', 'T'): 'T',
    ('A', 'G'): 'R', ('G', 'A'): 'R', ('C', 'T'): 'Y', ('T', 'C'): 'Y',
    ('G', 'C'): 'S', ('C', 'G'): 'S', ('A', 'T'): 'W', ('T', 'A'): 'W',
    ('G', 'T'): 'K', ('T', 'G'): 'K', ('A', 'C'): 'M', ('C', 'A'): 'M',
    ('C', 'G', 'T'): 'B', ('A', 'G', 'T'): 'D', ('A', 'C', 'T'): 'H', ('A', 'C', 'G'): 'V',
    ('A', 'C', 'G', 'T'): 'N'
}

# Function to compare two sequences and create consensus
def compare_sequences(seq1, seq2):
    # Check if the sequences have the same length
    if len(seq1) != len(seq2):
        raise ValueError(f"Sequences are not the same length: {len(seq1)} vs {len(seq2)}")

    consensus_seq = []
    for base1, base2 in zip(seq1, seq2):
        if base1 == base2:
            consensus_seq.append(base1)
        else:
            # Handle cases with two different bases using IUPAC codes
            consensus_seq.append(iupac_code.get(tuple(sorted(set([base1, base2]))), 'N'))

    return ''.join(consensus_seq)

# Function to compare FASTA files and generate consensus
def compare_fasta_files(file1, file2, output_file):
    try:
        records1 = SeqIO.to_dict(SeqIO.parse(file1, "fasta"))
        records2 = SeqIO.to_dict(SeqIO.parse(file2, "fasta"))
    except Exception as e:
        print(f"Error reading files: {e}")
        sys.exit(1)

    consensus_records = []

    for record_id in records1:
        if record_id in records2:
            seq1 = str(records1[record_id].seq)
            seq2 = str(records2[record_id].seq)
            
            # Compare sequences and ensure they are the same length
            try:
                consensus_sequence = compare_sequences(seq1, seq2)
            except ValueError as ve:
                print(f"Error comparing sequences for {record_id}: {ve}")
                sys.exit(1)

            consensus_record = records1[record_id]
            # Set the consensus sequence without any newlines
            consensus_record.seq = Seq.Seq(consensus_sequence)
            consensus_records.append(consensus_record)
        else:
            print(f"Record {record_id} is only found in the first file and will be skipped.")

    for record_id in records2:
        if record_id not in records1:
            print(f"Record {record_id} is only found in the second file and will be skipped.")

    # Write the consensus sequences to the output file as a single continuous line
    with open(output_file, "w") as output_handle:
        for record in consensus_records:
            # Write the header and the sequence in one line
            output_handle.write(f">{record.id}\n{str(record.seq)}\n")

    print(f"Consensus sequence has been saved to {output_file}")

def main(argv):
    file1 = ''
    file2 = ''
    output_file = ''

    try:
        # Define expected options with short flags
        opts, args = getopt.getopt(argv, "", ["input1=", "input2=", "output="])
        
        # Process options
        for opt, arg in opts:
            if opt == "--input1":
                file1 = arg
            elif opt == "--input2":
                file2 = arg
            elif opt == "--output":
                output_file = arg

        # Check that all required arguments are provided
        if not file1 or not file2 or not output_file:
            raise getopt.GetoptError("All three arguments (-input1, -input2, -output) are required")

    except getopt.GetoptError as err:
        print(f"Error: {err}")
        print('Usage: python3 iupac.py -input1 file1.fa -input2 file2.fa -output output.fa')
        sys.exit(2)

    # Call the function to compare files and generate consensus
    compare_fasta_files(file1, file2, output_file)
    
if __name__ == "__main__":
    # Call the main function with the provided arguments (excluding the script name)
    main(sys.argv[1:])