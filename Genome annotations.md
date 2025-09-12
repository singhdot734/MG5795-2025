# Class 5
## Visualizing genome annotations in NCBI RefSeq, Ensembl and UCSC
1. Check annotations for the gene TBX1 in the following databases:
  NCBI RefSeq: `https://www.ncbi.nlm.nih.gov/gene/6899`
  Ensembl: `https://www.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=ENSG00000184058;r=22:19756703-19783593`
  UCSC: `https://genome.ucsc.edu/` (search for TBX1, in Gene and Gene Predictions group, turn on Gencode V48, RefSeq and MANE tracks)

## Download GTF files from RefSeq and Ensembl
1. In your personal directory, make a new directory called `hg38_gtfs`. Change directory into this new directory.
2. Use `curl` to download GTF files from the following links. Direct output to a file using `>`
3. RefSeq GTF: `https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.ncbiRefSeq.gtf.gz`
4. Ensembl GTF: `https://ftp.ensembl.org/pub/release-115/gtf/homo_sapiens/Homo_sapiens.GRCh38.115.gtf.gz`
5. This is how the command will look for RefSeq `curl -s https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.ncbiRefSeq.gtf.gz > hg38.ncbiRefSeq.gtf.gz`
6. Both are gzipped files. Use `gunzip` to unzip the files. `-c` option will be needed if you want to save the original file. Don't forget to direct output to a file using `>`. You can use the same name as the downloaded file except for .gz.
7. This is how the command will look for the RefSeq GTF: `gunzip -c hg38.ncbiRefSeq.gtf.gz > hg38.ncbiRefSeq.gtf`

### Exercise 1: Once both files are downloaded and unzipped, do `pwd` and `ls` and take a screenshot to submit your progress.

# Class 6
## Searching for patterns in GTF files
1. How many lines are in each GTF file?
2. Take a look at each GTF file using `head`.
3. Too many lines can be hard to read. Try to print only one line on the screen using `head`, `tail` etc.
4. Use `grep` to look for TBX1 gene in each GTF file. Option `-w` may be needed to avoid matches such as TBX11 in the output. Instead of option `-w`, word boundry can be defined using `\b`, like `grep "TBX1\b"`
5. Pipe output into `head -2` to limit to the first two lines. Otherwise there are several lines containing the entry "TBX1". We will talk why.
6. Look at the individual tab separated fields of one line and cross check with GTF/GFF file format description in lecture slides to understand what is listed in each field of the GTF file.
7. How many lines have the entry "TBX1" in each GTF file?

## How many chromosomes are represented in the GTF files?
1. Before we go further, lets see what chromosomes are listed in GTF files.
2. We will introduce a new command to select certain columns. Meet `cut`, which can select a particular column. For example, `cut -f 2` will select column (field) 2.
3. We will first use `cat` to open data stream from a gtf file and then pipe it into `cut` to select column 1, and then sort and collapse column 1 contents. 
4. This is how it can be done `cat hg38.ncbiRefSeq.gtf | cut -f 1 | sort | uniq -c`.
5. What did you get?
6. Hmmm...there are lot more entries than just chr 1-22 plus chr X and Y. There are alternate assembly patches for each chromosome.
7. So, if we want to count chromosomes, we will have to refine our search patterns. For example, if we want to use a combination of `grep` for chr22 and `wc -l` to see how many features are on chr22, we will have to restrict counting to pattern "chr22", either using `-w` option or define boundaries on either side "\bchr22\b". This will force exact match to "chr22".
8. So, go ahead and try to count all the rows/lines that contain chr22. This will do it: `grep -w 'chr22' hg38.ncbiRefSeq.gtf | wc -l`
9. Using the above steps, similarly investigate all chromosome entries in the Ensembl GTF file.

## Counting number of genes and other features in GTF files
1. Column 3 of GTF file lists what gene feature is represented in a row. If we count how many rows have a particular feature, we can get an idea about their counts in a genome.
2. `grep`, `sort`, `uniq` can be used for collapsing and counting feautres in column 3 but pattern matching is harder in GTF file as other columns (e.g., column 9) can have same patterns. The solution is to first select certain columns and then do pattern matching and counting.
3. We will use `cat` to open data stream from a file. Then we will select particular column using `cut`.
4. To summarize all gene features present in RefSeq GTF file, we can do this: `cat hg38.ncbiRefSeq.gtf | cut -f 3 | sort | uniq`
5. These features can be counted by adding flag `-c` to `uniq`: `cat hg38.ncbiRefSeq.gtf | cut -f 3 | sort | uniq -c`
6. To limit our summary or counting to a particular chr, say chr22, we can also add column 1 and then `grep` for chr22. Modify the command in step 5: `cat hg38.ncbiRefSeq.gtf | cut -f 1,3 | grep -w 'chr22' | sort | uniq -c`
7. So, how many transcripts on chr22?

### Exercise 1: Ensembl GTF file column 3 also has entries for "gene" (missing in RefSeq GTF). Write a command to summarize and count all features including genes in the entire Ensembl GTF file and on chr22. Note the difference in how chromosome numbers are listed in column 1 of the Ensembl GTF file. So, how many genes are in the Ensembl GTF? And how many on chr22?

## Meet AWK to extract specific information from large text (GTF) files
1. Very often in genomic data analysis, you will need to extract rows or columns from a text file like a GTF file to make a new text file. One of the most common text file needed during analysis is a BED file. We will learn some tricks of a simple yet powerful new programming language called `awk`.
2. To make a custom BED file of all genes in Ensembl GTF file,  this awk one-liner will do the trick: `awk -F'\t' '{if ($3=="gene") {print $1,$4-1,$5,$3,$6,$7}}' Homo_sapiens.GRCh38.115.gtf | head`
3. First use `head` to view the top 10 lines of the newly generated BED file.
4. If the command works, send output to a file called `hg38_genes.bed`.
5. Column 4 has the same value in every row and is not very useful. It will be better to have gene_id in this column but that requires learning some more `awk` tricks.

## So what are all the types of genes annotated in the human genome
1. Ensembl GTF file column 9 has a tag called "gene_biotype", which specifies various types of gene classes, e.g., protein coding, miRNA, etc.
2. If we can select all the rows that have column 3 = gene and then extract gene_biotype value into a new column.
3. To do this, we want to learn how to extract specific information from rows and columns and save it in a new column.
4. Here is an `awk` command that will do what we said in step 2. Note that output is piped into `head` so that we can make sure the command works and produces an output we intend to: `awk -F'\t' '$3 == "gene" {match($9, /gene_biotype "([^"]+)";/); gene_biotype = substr($9, RSTART + 14, RLENGTH - 15); print $0 "\t" gene_biotype}' Homo_sapiens.GRCh38.115.gtf | head`
5. We will break this up into pieces and examine what it does and why.
6. Note that the gene_biotype has been added into a new column #10.

### Exercise 2: Pipe the output of command in step 4 into `cut` to isolate column #10, and then use `sort` and `uniq` to collapse and count the number of gene_biotype features. Paste the screenshot and output of the final command and submit.



