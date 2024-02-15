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


----------
#Script alterado para correção de cores 
library(ggplot2)
library(readxl)
library(dplyr)

# List of file names
file_names <- c("Lup2108.xlsx", "Lup2163.xlsx", "Lup2164.xlsx", "Lup2189.xlsx","Lup2191.xlsx","Lup2202.xlsx","Lup2221.xlsx","Lup2223.xlsx","Lup2241.xlsx","Lup2312.xlsx","Lup2313.xlsx","Lup2323.xlsx","Lup2370.xlsx","Lup2378.xlsx","Lup2383.xlsx","Lup2395.xlsx","Lup2453.xlsx","Lup2456.xlsx")

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
custom_labels <- c("Lup2108","Lup2163","Lup2164","Lup2189","Lup2191","Lup2202","Lup2221","Lup2223","Lup2241","Lup2312","Lup2313","Lup2323","Lup2370","Lup2378","Lup2383","Lup2395", "Lup2453","Lup2456")

# Create a color palette corresponding to your samples
custom_colors <- rainbow(length(custom_labels))

# Update the ggplot code with scale_color_manual to customize the legend
ggplot(combined_data,
       aes(
         x = GCContent,
         y = Count,
         color = Sample,
       )) +
  geom_line() +
  labs(title = "GC Content vs. Count", x = "GC Content", y = "Count") +
  scale_color_manual(values = custom_colors, labels = custom_labels)
  
  ggsave("greyplot.pdf")

-------
#Script alterado para ter linhas a cinza e apenas higlighted as samples que quero 

library(ggplot2)
library(readxl)
library(dplyr)

# List of file names
file_names <- c("Lup2108.xlsx", "Lup2163.xlsx", "Lup2164.xlsx", "Lup2189.xlsx","Lup2191.xlsx","Lup2202.xlsx","Lup2221.xlsx","Lup2223.xlsx","Lup2241.xlsx","Lup2312.xlsx","Lup2313.xlsx","Lup2323.xlsx","Lup2370.xlsx","Lup2378.xlsx","Lup2383.xlsx","Lup2395.xlsx","Lup2453.xlsx","Lup2456.xlsx")

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
custom_labels <- c("Lup2108","Lup2163","Lup2164","Lup2189","Lup2191","Lup2202","Lup2221","Lup2223","Lup2241","Lup2312","Lup2313","Lup2323","Lup2370","Lup2378","Lup2383","Lup2395", "Lup2453","Lup2456")

custom_colors <- c('grey',"blue", "red", "green",'grey','grey','grey','grey','grey','grey','grey','grey','grey','grey','grey','grey','grey','grey')

custom_size <- c(0.5, 2, 2, 2, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5)

# Update the ggplot code with scale_color_manual to customize the legend
ggplot(combined_data, aes(x = GCContent, y = Count, color = Sample, size = Sample)) +
  geom_line() +
  labs(title = "GC Content vs. Count", x = "GC Content", y = "Count") +
  scale_color_manual(values = custom_colors, labels = custom_labels) +
  scale_size_manual(values = custom_size, labels = custom_labels)

ggsave("greyplot.pdf")




