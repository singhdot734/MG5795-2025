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
