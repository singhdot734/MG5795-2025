# R - class2
```
setwd("/fs/ess/PAS3124/MOLGEN_5795_OSU/materials")

# Read the bed file with gene IDs, gene lengths and strand
data1 <- read.table("pc_genes_len.bed", header = FALSE, sep = "\t")

# Understanding the data frame
# Dimension of the data frame
dim(data1)

# Summary of the data frame
summary(data1)

# Structure of the data frame
str(data1)

# Column 2 is numeric; basic arithmetic functions on $V2
mean(data1$V2)
min(data1$V2)
max(data1$V2)

# Look at head and tail of the data frame
head(data1)
tail(data1)

# Get specific values from the data frame
# Remember, data frame is two-dimensional (rows and columns).
# Therefore, to select a specific value we can use [ ] (brackets).

# Value in row 1, column 1 (row is on the left, column on right)
data1[1,1]

data1[5,2]

# All values in row 1
data1[1,]

# Column 1 values
data1[,1]
data1$V1

# All column 1 and 2 values
data1[1:2]

# Plotting the data in the data frame
# Basic histogram
hist(data1$V2)

# histogram of log10 transformed gene lengths
hist(log10(data1$V2))

# more customized histogram of log10 transformed gene lengths
hist(log10(data1$V2),
     main = "hg38 Gene Lengths",
     xlab = "Gene Length (log10(bp))",
     ylab = "Gene Count",
     col = "lightblue",
     border = "black")

# even more customized histogram of log10 transformed gene lengths
hist(log10(data1$V2),
     main = "hg38 Gene Lengths",
     xlab = "Gene Length (log10(bp))",
     ylab = "Gene Count",
     breaks = 20,
     xlim = c(1,7),
     col = "blue",
     border = "black")

# Plot smooth density curve
plot(density(data1$V2), 
     main = "hg38 Gene Lengths", 
     xlab = "Gene Length (log10(bp))",
     col = "red", 
     lwd = 2)

# Plot smooth density curve
plot(density(log10(data1$V2)), 
     main = "hg38 Gene Lengths", 
     xlab = "Gene Length (log10(bp))",
     col = "red", 
     lwd = 2)

# Histogram version 2 - split genes by strand to plot gene lengths from two strands
# First, lets find all rows with + strand - outputs a logical vector
data1$V3 == "+"

# filtering for plus strand rows - pass the logical vector to return only TRUE rows 
data1[data1$V3 == "+",]

# Read the data and give column names
data2 <- read.table("pc_genes_len.bed", header = FALSE, sep = "\t",
                   col.names = c("ENSG_ID", "Gene_Length", "Strand"))

#Add a new column Log10_Length that has log10-transformed gene length
data2$Log10_Length <- log10(data2$Gene_Length)

# Split data by strand
plus_strand <- subset(data2, Strand == "+")
minus_strand <- subset(data2, Strand == "-")

# Histogram for + strand
hist(plus_strand$Log10_Length,
     main = "Gene Lengths on + Strand",
     xlab = "log 10(Gene Length, bp)",
     ylab = "Gene counts",
     xlim = c(1,8),
     col = "gray",
     border = "black",
     breaks = 20)

# Histogram for - strand
hist(minus_strand$Log10_Length,
     main = "Gene Lengths on - Strand",
     xlab = "log10(Gene Length, bp)",
     ylab = "Gene counts",
     xlim = c(1,8),
     col = "black",
     border = "black",
     breaks = 20)

# Plot density for + strand
plot(density(plus_strand$Log10_Length), 
     col = "blue", 
     lwd = 2,
     main = "hg38 Gene Length",
     xlab = "log10(Gene Length)")

# Add density for - strand
lines(density(minus_strand$Log10_Length), col = "red", lwd = 2)

# Add legend
legend("topright", 
       legend = c("+ strand", "- strand"),
       col = c("blue", "red"), 
       lwd = 4)

library(ggplot2)

# Create histogram
ggplot(data2, aes(x = log10(Gene_Length), fill = Strand)) +
  geom_histogram(position = "identity", alpha = 0.6, bins = 30) +
  labs(title = "Histogram of Gene Lengths by Strand",
       x = "Gene Length",
       y = "Count") +
  theme_minimal() +
  scale_fill_manual(values = c("blue", "red"))

```
# R - class 3

```
# Class 3
setwd("~/materials")
library(ggplot2) # libraries (this one and any other) have to be loaded in every session

# Read the data and give column names
data2 <- read.table("pc_genes_len.bed", header = FALSE, sep = "\t",
                    col.names = c("ENSG_ID", "Gene_Length", "Strand"))

# ggplot usage
# ggplot(data = <DATA>, mapping = aes(<MAPPINGS>)) +  <GEOM_FUNCTION>()

# Create histogram
ggplot(data2, aes(x = log10(Gene_Length), fill = Strand)) + # map x variable to the aesthetics function of ggplot
  geom_histogram()

# Modify histogram to not stack the bars
ggplot(data2, aes(x = log10(Gene_Length), fill = Strand)) + # assign x variable to the aesthetics function of ggplot
  geom_histogram(position = "identity")

# Make histograms more transparent
ggplot(data2, aes(x = log10(Gene_Length), fill = Strand)) + # assign x variable to the aesthetics function of ggplot
  geom_histogram(position = "identity", alpha = 0.6)

# Add labels to histogram and custom color
ggplot(data2, aes(x = log10(Gene_Length), fill = Strand)) + #assign x variable to the aesthetics function of ggplot
  geom_histogram(position = "identity", alpha = 0.6) +
  labs(title = "Histogram of Gene Lengths by Strand",
       x = "Gene Length",
       y = "Count") +
  scale_fill_manual(values = c("purple", "darkgreen")) # add color to two strand vectors

# Create density plot
ggplot(data2, aes(x = log10(Gene_Length), fill = Strand)) +
  geom_density(alpha = 0.1) +
  labs(title = "Histogram of Gene Lengths by Strand",
       x = "Gene Length",
       y = "Count") +
  scale_fill_manual(values = c("blue", "red"))

# Load dplyr and readr packages to use tidy data formats

library(dplyr)
library(readr)
library(tidyr)

# Read data using read_csv()
data3 <- read_csv("Horste_S1_mod.csv")

# glimpse() to see the structure of data
glimpse(data3)

# mode = dbl signifies "double-precision" numeric vector

# Introducing Pipe - %>% - and count()
data3 %>% count(Loc) # this works like "cut -f n | sort | uniq -c" on command line

# Changing mode of a column using mutate()
data3_clean <- data3 %>%
  mutate(
    mRNA_length_mane = as.numeric(mRNA_length_mane) # change mode
  )

data3_clean %>% count(mRNA_length_mane)

# Changing mode of additional numeric columns
data3_clean <- data3 %>%
  mutate(
    mRNA_length_mane = as.numeric(mRNA_length_mane), # change mode
    Protein_length_mane = as.numeric(Protein_length_mane), # change mode
    Exon_number_mane = as.numeric(Exon_number_mane), # change mode
    av_CDS_exon_length_mane = as.numeric(av_CDS_exon_length_mane) # change mode
  )

glimpse(data3_clean)

# How many variables are for mRNA localization score (Loc)?
data3_clean %>% count(Loc)

# Let's say we want to test if there is a relationship between mRNA length and protein length.

# Create scatter plot
ggplot(data3_clean, aes(x = mRNA_length_mane, y = Protein_length_mane)) + # map the x and y variables to ggplot aesthetics
  geom_point() # makes a basic scatter plot
    
# Customize the scatter plot
ggplot(data3_clean, aes(x = mRNA_length_mane, y = Protein_length_mane)) +
  geom_point(color = "blue", alpha = 0.6) +
  labs(
    title = "Scatter Plot of mRNA Length vs Protein length",
    x = "mRNA length",
    y = "Protein length"
  )

# Create scatter plot with log10
ggplot(data3_clean, aes(x = log10(mRNA_length_mane), y = log10(Protein_length_mane))) +
  geom_point(color = "lightblue", alpha = 0.6) +
  labs(
    title = "Scatter Plot of mRNA Length vs Protein length",
    x = "mRNA length",
    y = "Protein length"
  )

# Create scatter plot with a linear regression line
ggplot(data3_clean, aes(x = mRNA_length_mane, y = Protein_length_mane)) +
  geom_point(color = "red", alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE) + # Add linear regression line (lm) with confidence interval (se = TRUE)
  labs(
    title = "Scatter Plot of mRNA Length vs Protein length",
    x = "mRNA Length",
    y = "Protein length"
  )

# Let's see if there is any relationship between exon number and mRNA localization
ggplot(data3_clean, aes(x = Loc, y = Exon_number_mane)) +
  geom_boxplot()

# Customize the box plot
ggplot(data3_clean, aes(x = Loc, y = Exon_number_mane)) +
  geom_boxplot(fill = "skyblue", color = "darkblue") +
  labs(title = "Exon number vs mRNA localization",
       x = "Location",
       y = "Exon Number") +
  theme_minimal() # this will make a boxplot (or any other geom object without background)


# Select certain columns
select(data3_clean, mRNA_length_mane, Protein_length_mane) #select(which_df, col_1, col_2)

# filter rows
filter(data3_clean, Loc == "ER") #filter(which_df, col == row_value)

```
