import vcf
from Bio import SeqIO


def read_vcf(vcf_file):
    records = []
    vcf_reader = vcf.Reader(open(vcf_file, 'r'))

    for record in vcf_reader:
        if record.samples[0]['GT'] == '1/1':
            records.append({
                'CHROM': record.CHROM,
                'POS': record.POS,
                'REF': record.REF,
                'ALT': str(record.ALT[0])
            })
    return records


def update_fasta(fasta_file, vcf_records, output_file):
    updated_sequences = []

    # Load the original FASTA file
    sequences = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))

    # Update the sequences based on VCF records
    for vcf_record in vcf_records:
        chrom = vcf_record['CHROM']
        pos = vcf_record['POS']
        ref = vcf_record['REF']
        alt = str(vcf_record['ALT'])

        if chrom in sequences:
            original_sequence = sequences[chrom].seq

            # Check if the position in the sequence matches the VCF record
            if len(original_sequence) >= pos and original_sequence[pos - 1] == ref:
                # Update the sequence
                updated_sequence = original_sequence[:pos -
                                                     1] + alt + original_sequence[pos:]
                sequences[chrom].seq = updated_sequence

    # Save the updated sequences to a new FASTA file
    with open(output_file, 'w') as output_handle:
        SeqIO.write(sequences.values(), output_handle, 'fasta')


sample_file = input("Sample: ")
vcf_file = sample_file + ".flt.vcf"
fasta_file = input("Fasta: ")
output_file = sample_file + "_11.fa"
vcf_records = read_vcf(vcf_file)

update_fasta(fasta_file, vcf_records, output_file)
