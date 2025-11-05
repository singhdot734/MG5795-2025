# RNA-Seq analysis
As all tools needed from OSC for this analysis are available only in the Ascend cluster, all command line tasks should be done on this cluster (Trimmomatic is available only on this cluster). To connect to this cluster from the OnDemand dashboard, click the arrow to the right of "Open in terminal" button and choose "Ascend" from the list of three clusters.
 
## Step I: Downloading data
1. Load sratoolkit `module load sratoolkit/3.0.2`
2. The following command will download data from a SRA run ID into the same directory where the command will be issued from: `fasterq-dump --split-files <SRR#>`
3. Replace `<SRR#>` with the SRA Run ID to be downloaded.
4. The --split-files option ensures the R1 and R2 files are separate for paired end data.

## Step II: Trimming adapters
1. Adapaters to be trimmed will need to be provided in fasta format.
2. Create this file in the OnDemand page in the classroom account as follows:
3. +New file > In the pop-up window, give it a name "adapters.fa".
4. To the right of the new file, press the drop-down button to edit it.
5. In the new tab that opens, paste the following:
```
#adapters.fa
>IlluminaUniversal
AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
```
6. If there are any other adapters used in samples, add their name and sequence on new lines in fasta format.
7. Save the file.
8. Load trimmomatic in your environment: `module load trimmomatic/0.38`
9. The following trimmomatic script can be used to trim adapters these adapters from paired-end samples.

```
  trimmomatic PE \
  SRR#_R1.fastq SRR#_R2.fastq \ 
  SRR#_paired_F.fastq SRR#_unpaired_F.fastq SRR#_paired_R.fastq SRR#_unpaired_R.fastq \
  ILLUMINACLIP:adapters.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:18 
```
10. Replace `SRR#` with file name of the input and output fasta file names.
11. As no path is given for each input and output file, the above script assumes input files (line 2 above) are in the same directory where the above command is issued. Output fastq files are on line 3.
12. The `adapters.fa` file also has to be in this directory.
13. Output files will also be in this directory.
14. See trimmomatic manual for more details for options, as needed: https://github.com/usadellab/Trimmomatic

## Step III: Alignment
1. Following items will be needed to do the alignment:
   a. Genome index (STAR index), which is in this directory: /fs/ess/PAS3124/data/STAR_index/GRCh38_gencode48_index_STAR_v2_7_11b
   b. Trimmed fastq files
   c. A shell bash script to submit/queue the alignment job OSC compute nodes. This script is below.
2. To create a SBATCH script in OnDemand, do as follows:
3. +New file > In the pop-up window, give it a name "<my_star_alignment_script.sh>".
4. To the right of the new file, press the drop-down button to edit it.
5. In the new tab that opens, paste the following:
```
#!/bin/bash
#SBATCH --account=PAS3124
#SBATCH --time=02:00:00
#SBATCH --mem=35G
#SBATCH --cpus-per-task=8
#SBATCH -o run_star_<sample_id>.log

module load star/2.7.10b

genome_index=/fs/ess/PAS3124/data/STAR_index/GRCh38_gencode48_index_STAR_v2_7_11b
R1_fq=/fs/ess/PAS3124/MOLGEN_5795_OSU/path_to/<sample_id>_paired_F.fastq
R2_fq=/fs/ess/PAS3124/MOLGEN_5795_OSU/path_to/<sample_id>_paired_R.fastq
out_prefix=<sample_id>_

STAR \
    --readFilesIn $R1_fq $R2_fq \
    --genomeDir $genome_index \
    --outFileNamePrefix $out_prefix \
    --runThreadN 8 \
    --quantMode GeneCounts
```
6. In the above file, absolute paths are provided for the variables `genome_index`, `R1_fq` and `R2_fq`. Replace these paths as needed where your files are stored. Also make sure to update <sample_id> where ever needed in the script.
7. Run the job submission script `sbatch <path_to/my_star_alignment_script.sh>`
8. To check the progress of your job: `squeue -u <your username>` or `squeue --job <JOBID>`
9. Output files will come into the directory where the sbatch command is run from.
10. For more information on STAR alignment, consult the manual: https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf

## Step IV: Creating bam files and visualization
1. Use samtools to convert sam file into sorted bam file.
```
module load samtools/1.21
samtools view -b -o <sample_id>_unsorted.bam <sample_id>_Aligned.out.sam
samtools sort -o <sample_id>_sorted.bam <sample_id>_unsorted.bam
```
2. Subset the sorted bam file, in this case to capture only chr20 alignments, and save as sam file.
3. Convert back to bam, sort and index. 
``` 
samtools view -h <sample_id>_sorted.bam chr20 > <sample_id>_sorted.chr20.sam
samtools view -b <sample_id>_sorted.chr20.sam > <sample_id>_sorted.chr20.bam
samtools sort -o <sample_id>_sorted.chr20.bam <sample_id>_sorted.chr20.bam # this step may not be needed.
samtools index <sample_id>_sorted.chr20.bam
```
4. The bam and corresponding bam.bai can be uploaded to IGV.

## Step V: Creating gene counts table
1. To perform differential expression analysis, first step is to create a gene counts table. As the dataset we are working with is likely unstranded, we will use the unstranded counts in the `_ReadsPerGene.out.tab` files from the `--quantMode GeneCounts` output from STAR alignment.
2. Start an R Studio app in the OSC classroom.
3. The following script will make a counts table in a format that can be used for differential expression analysis using DESeq2.
   
```
setwd("~/materials")

# Make a counts table

#read in the counts for control siRNA 1
pc3_siC_1_count_data = read.table("larp4B_star_counts/pc3_siC_1_align_ReadsPerGene.out.tab",
                                  sep = "\t",
                                  skip = 4,
                                  col.names = c("gene_id", "unstranded", "forward_stranded", "reverse_stranded"))
#read in the counts for control siRNA 2
pc3_siC_2_count_data = read.table("larp4B_star_counts/pc3_siC_2_align_ReadsPerGene.out.tab",
                                  sep = "\t",
                                  skip = 4,
                                  col.names = c("gene_id", "unstranded", "forward_stranded", "reverse_stranded"))
#read in the counts for LARP4B siRNA 1
pc3_si4B_1_count_data = read.table("larp4B_star_counts/pc3_si4B_1_align_ReadsPerGene.out.tab",
                                   sep = "\t",
                                   skip = 4,
                                   col.names = c("gene_id", "unstranded", "forward_stranded", "reverse_stranded"))
#read in the counts for LARP4B siRNA 2
pc3_si4B_2_count_data = read.table("larp4B_star_counts/pc3_si4B_2_align_ReadsPerGene.out.tab",
                                   sep = "\t",
                                   skip = 4,
                                   col.names = c("gene_id", "unstranded", "forward_stranded", "reverse_stranded"))

#check the dimensions of the data frames we have just created
dim(pc3_siC_1_count_data)
dim(pc3_siC_2_count_data)
dim(pc3_si4B_1_count_data)
dim(pc3_si4B_2_count_data)

#check that the gene ids appear in the same order in each file, all or any options 
stopifnot(all(pc3_siC_1_count_data$gene_id == pc3_siC_2_count_data$gene_id))
stopifnot(all(pc3_siC_2_count_data$gene_id == pc3_si4B_1_count_data$gene_id))
stopifnot(all(pc3_si4B_1_count_data$gene_id == pc3_si4B_2_count_data$gene_id))
stopifnot(all(pc3_si4B_2_count_data$gene_id == pc3_siC_1_count_data$gene_id))

#create count data frame for the control samples
count_table = data.frame("pc3_siC_1" = pc3_siC_1_count_data$unstranded,
                         "pc3_siC_2" = pc3_siC_2_count_data$unstranded,
                         "pc3_si4B_1" = pc3_si4B_1_count_data$unstranded,
                         "pc3_si4B_2" = pc3_si4B_2_count_data$unstranded)

#get a vector of gene ids to add as row names
gene_ids = pc3_siC_1_count_data$gene_id 
rownames(count_table) <- gene_ids

# write count table as csv file
write.csv(count_table,
          "pc3_siC_si4B_count_table.csv")
```

4. Next, we will compare counts in the two replicates for each knockdown. We expect the replicate samples to have similar counts.
5. These steps will use the `count_table` made in the previous steps. So, these should follow in the same environment as steps 1-3.
   
```
library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)

#identify genes with sufficient counts to analyse
count_at_least_10 = count_table >= 10
num_samples_with_enough_counts_each_gene = rowSums(count_at_least_10)
gene_passes_filter = num_samples_with_enough_counts_each_gene >= 3
count_table_filtered = count_table[gene_passes_filter, ]

# compare control replicates - scatter plot 
ggplot(count_table_filtered, aes(x = log10(pc3_siC_1), y = log10(pc3_siC_2))) + # map the x and y variables to ggplot aesthetics
  geom_point(color = "blue", alpha = 0.6) +
  labs(
    title = "PC3 siNC1 vs siNC2",
    x = "Control siRNA 1",
    y = "Control siRNA 2"
  )

# compare LARP4B kd replicates - scatter plot 
ggplot(count_table_filtered, aes(x = log10(pc3_si4B_1), y = log10(pc3_si4B_2))) + # map the x and y variables to ggplot aesthetics
  geom_point(color = "red", alpha = 0.6) +
  labs(
    title = "PC3 siLARP4B_1 vs siLARP4B_2",
    x = "LARP4B kd 1",
    y = "LARP4B kd 2"
  )
```

## Step VI: DESeq2
1. The steps below will run DESeq2 algorithm for differential expression analysis between the control and LARP4B knockdown samples.
2. In a new R script file, paste the following code: 

```
# DESeq2 analysis

setwd("~/materials")

library("DESeq2")
library("biomaRt")

#read in count table
count_table = read.csv("pc3_siC_si4B_count_table.csv",
                       row.names = 1)

#identify genes with sufficient counts to analyse
count_at_least_10 = count_table >= 10
num_samples_with_enough_counts_each_gene = rowSums(count_at_least_10)
gene_passes_filter = num_samples_with_enough_counts_each_gene >= 3
count_table_filtered = count_table[gene_passes_filter, ]

#define metadata
metadata = data.frame("condition" = c("ctrl", "ctrl", "kd", "kd"))
rownames(metadata) = colnames(count_table_filtered)

# turn condition column of metadata into a factor
metadata$condition = factor(metadata$condition)

#create DESeqDataSet Object
DE_dataset = DESeqDataSetFromMatrix(countData = count_table_filtered,
                                    colData = metadata,
                                    design = ~ condition)

# Make a PCA plot
# first do variance stabilizing transformation
vst = varianceStabilizingTransformation(DE_dataset)

# Now plot PCA
plotPCA(vst)

# Set control knockdown as reference level
DE_dataset$condition <- relevel(DE_dataset$condition, "ctrl")

#run the differential expression analysis
DE_dataset = DESeq(DE_dataset)

# plot MA plot
plotMA(DE_dataset,
       alpha = 0.1,
       main = "",
       xlab = "mean of normalized counts",
       colNonSig = "gray60",
       colSig = "blue",
       colLine = "grey40",
       returnData = FALSE,
       MLE = FALSE
       )

#get the results for knockdown vs control comparison
condition_results = results(DE_dataset, 
                            contrast = c("condition","kd", "ctrl"),
                            independentFiltering = FALSE)

#connect to ensembl specifying v114 (gencode v48) annotation
ensembl_114 = useEnsembl(biomart = "genes", 
                         dataset = "hsapiens_gene_ensembl",
                         version = 114)

#get a data frame mapping ensembl gene ids to HGNC symbols (gene names)
gene_id_table = getBM(mart = ensembl_114,
                      attributes = c("ensembl_gene_id", "hgnc_symbol"))

#remove duplicates for very rare cases an ensembl id maps to multiple symbols
gene_id_table = gene_id_table[!duplicated(gene_id_table$ensembl_gene_id), ]

#remove gene version from ensembl id in results table
rownames(condition_results) = gsub("\\.[0-9]{1,3}", "", rownames(condition_results))

#add HGNC symbol to our results table
condition_results_with_names = merge(gene_id_table, as.data.frame(condition_results), all.y = TRUE, by.x = "ensembl_gene_id", by.y = "row.names")

#save the results to a file
write.csv(condition_results_with_names, "siC_vs_siLARP4B_results.csv", row.names = FALSE)
```
3. Next, the DESeq2 output results csv will be used to understand the type of effect LARP4B has on gene expression.
4. The block of code below will generate a volcano plot of differentially expressed genes and will output lists of such genes for GO term analysis.
```
# Analysis of DESeq2 results

setwd("~/materials")

library(ggplot2)
library(dplyr)
library(readr)
library(tidyr)

LARP4B_DEseq <- read_csv("siC_vs_siLARP4B_results.csv")

glimpse(LARP4B_DEseq)

# Select certain columns and get summary
log2fc <- select(LARP4B_DEseq, log2FoldChange) 
summary(log2fc)

# Create a new column for 2-fold significance
LARP4B_DEseq$sig_2f <- with(LARP4B_DEseq, ifelse(log2FoldChange > 1 & padj < 0.05, "Up",
                                                       ifelse(log2FoldChange < -1 & padj < 0.05, "Down", "Not significant")))

# Volcano Plot 1 (2-fold change)
ggplot(LARP4B_DEseq, aes(x = log2FoldChange, y = -log10(padj), color = sig_2f)) +
  geom_point(alpha = 0.6) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey") +
  scale_color_manual(values = c("Up" = "red", "Down" = "blue", "Not significant" = "black")) +
  labs(title = "LARP4B vs Control KD (2-fold change)", x = "Log2 Fold Change", y = "-log10(Adjusted P-value)") +
  theme_minimal()

# Next we will make a list of genes that go significantly up and down. But first remove rows where hgnc_symbol is NA
LARP4B_DEseq <- subset(LARP4B_DEseq, !is.na(hgnc_symbol))

# Subset for Up-regulated genes
up_genes <- subset(LARP4B_DEseq, significance == "Up")

# Subset for Down-regulated genes
down_genes <- subset(LARP4B_DEseq, significance == "Down")

# Write to CSV (only hgnc_symbol column)
write.csv(up_genes$hgnc_symbol, "up_genes.csv", row.names = FALSE)
write.csv(down_genes$hgnc_symbol, "down_genes.csv", row.names = FALSE)
write.csv(LARP4B_DEseq$hgnc_symbol, "all_genes.csv", row.names = FALSE)
```
