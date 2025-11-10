install.packages("readxl")
install.packages("vioplot")

library(readxl)
library(vioplot)

create_violin_plot <- function(file_path, label) {
  data <- read_excel(file_path)
  return(data$`Final Sum`)
}

#file_paths <- c( "Lup2163.xlsx", "Lup2189.xlsx", "Lup2221.xlsx", "Lup2241.xlsx", "Lup2223.xlsx")
#plot_labels <- c( "Lup2163", "Lup2189", "Lup2221", "Lup2241", "Lup2223") 
#plot_colors <- c("lightblue")

file_paths <- c( "Lup2108.xlsx", "Lup2163.xlsx", "Lup2164.xlsx", "Lup2189.xlsx", "Lup2191.xlsx", "Lup2202.xlsx","Lup2221.xlsx", "Lup2223.xlsx", "Lup2241.xlsx", "Lup2312.xlsx", "Lup2313.xlsx", "Lup2323.xlsx", "Lup2370.xlsx", "Lup2378.xlsx")
plot_labels <- c( "Lup2108", "Lup2163", "Lup2164", "Lup2189", "Lup2191", "Lup2202","Lup2221","Lup2223", "Lup2241", "Lup2312", "Lup2313", "Lup2323", "Lup2370", "Lup2378")
plot_colors <- c("blue")

data_list <- lapply(file_paths, create_violin_plot)

pdf("Violinplott.pdf", width = 1 * length(data_list))

vioplot(data_list, names = plot_labels, col = plot_colors, horizontal = FALSE)

abline(h = 0.5, col = "black", lty = 2)

dev.off()


#png("vioplot.png") 
#vioplot(data_list, names = plot_labels, col = plot_colors, border = plot_borders, horizontal = FALSE)

#abline(h = 0.5, col = "black", lty = 2)  # You can customize the color and line type

#dev.off()

-------- vioplot
import vcf

filename = input("Sample: ")
vcf_file = vcf.Reader(open(filename+'.flt1.snps.vcf', 'r'))
samples_info = []

with open(filename+'.txt', 'w') as output_file:
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
        output_file.write(
            f"Scaffold: {sample_info['Scaffold']}, POS: {sample_info['POS']}, GT: {sample_info['GT']}, DP4: {sample_info['DP4']}, DP4 Total: {sample_info['DP4 Total']}, Valor: {sample_info['Final Sum']}\n")

output_file.close()

------

import vcf
import xlsxwriter

filename = input("Sample: ")
vcf_file = vcf.Reader(open(filename + '.flt1.snps.vcf', 'r'))
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
                'Final Sum': dp4[0] + dp4[1] / sum(dp4)
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



