import vcf
import xlsxwriter

filename = input("Sample: ")
vcf_file = vcf.Reader(open(filename + '.flt2.vcf', 'r'))
samples_info = []

workbook = xlsxwriter.Workbook(filename + '.xlsx')
worksheet = workbook.add_worksheet()

header_format = workbook.add_format(
    {'bold': True, 'align': 'center', 'valign': 'vcenter'})

worksheet.write('A1', 'Scaffold', header_format)
worksheet.write('B1', 'POS', header_format)
worksheet.write('C1', 'GT', header_format)
worksheet.write('D1', 'DP', header_format)
worksheet.write('E1', 'DP4', header_format)
worksheet.write('F1', 'DP4 Total', header_format)
worksheet.write('G1', 'Final Sum', header_format)

row = 1

for record in vcf_file:
    for sample in record.samples:
        dp4 = record.INFO.get('DP4')
        if sample['GT'] == '0/1':
            samples_info.append({
                'Scaffold': record.CHROM,
                'POS': record.POS,
                'GT': sample["GT"],
                'DP': sample["DP"],
                'DP4': dp4,
                'DP4 Total': sum(dp4),
                'Final Sum': (dp4[0] + dp4[1]) / sum(dp4)
            })

for sample_info in samples_info:
    worksheet.write(row, 0, sample_info['Scaffold'])
    worksheet.write(row, 1, sample_info['POS'])
    worksheet.write(row, 2, sample_info['GT'])
    worksheet.write(row, 3, sample_info['DP'])
    worksheet.write(row, 4, str(sample_info['DP4']))
    worksheet.write(row, 5, sample_info['DP4 Total'])
    worksheet.write(row, 6, sample_info['Final Sum'])
    row += 1

workbook.close()
