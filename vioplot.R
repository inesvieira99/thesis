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
