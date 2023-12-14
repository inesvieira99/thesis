library(ggplot2)
library(readxl)
library(dplyr)

# List of file names
file_names <- c("Lup2108.xlsx", "Lup2163.xlsx", "Lup2164.xlsx", "Lup2189.xlsx", "Lup2191.xlsx" , "Lup2202.xlsx" , "Lup2221.xlsx" , "Lup2223.xlsx", "Lup2241.xlsx", "Lup2312.xlsx","Lup2313.xlsx" , "Lup2323.xlsx" , "Lup2370.xlsx" , "Lup2370.xlsx")

# Initialize an empty list to store data frames
data_list <- list()

# Read and combine data from the XLSX files
for (file in file_names) {
  data <- read_xlsx(file)
  data_list <- append(data_list, list(data))
}

# Combine all data frames into a single data frame
combined_data <- bind_rows(data_list, .id = "Sample")

# Create a vector of custom legend labels in the desired order
custom_legend_labels <- c("Lup2108", "Lup2163", "Lup2164", "Lup2189", "Lup2191" , "Lup2202" , "Lup2221", "Lup2223", "Lup2241", "Lup2312","Lup2313" , "Lup2323" , "Lup2370" , "Lup2378")

# Create a color palette corresponding to your samples
custom_colors <- c("red", "blue", "green", "purple", "orange", "pink", "violet" , "yellow" , "black", "brown", "turquoise","sienna", "gray", "gold")

# Update the ggplot code with scale_color_manual to customize the legend
ggplot(combined_data, aes(x = GCContent, y = Count, color = Sample, group = Sample)) +
  geom_line() +
  labs(title = "GC Content vs. Count", x = "GC Content", y = "Count") +
  scale_color_manual(values = custom_colors, labels = custom_legend_labels)

png(plot.png)
