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

## Searching for patterns in GTF files
1. How many lines are in each GTF file?
2. Take a look at each GTF file using `head`.
3. Too many lines can be hard to read. Try to print only one line on the screen using `head`, `tail` etc.
4. Use `grep` to look for TBX1 gene in each GTF file. Option `-w` may be needed to avoid matches such as TBX11 in the output. Instead of option `-w`, word boundry can be defined using `\b`, like `grep "TBX1\b"`
5. Pipe output into `head -2` to limit to the first two lines. Otherwise there are several lines containing the entry "TBX1". We will talk why.
6. Look at the individual tab separated fields of one line and cross check with GTF/GFF file format description in lecture slides to understand what is listed in each field of the GTF file.
7. How many lines have the entry "TBX1" in each GTF file?
8. Use a combination of `grep` for chr22 and `wc -l` to see how many features are on chr22?
9. There are alternate assembly patches for each chromosome, which will also be counted. To restrict counting to pattern "chr22", either use `-w` option or define boundaries on either side "\bchr22\b".

## Counting number of genes and other features in GTF files
1. Column 3 of GTF file lists what gene feature is represented in a row. This can be used to count (or select) lines with particular features.
2. `grep`, `sort`, `uniq` can be used for collapsing and counting feautres in column 3 but pattern matching is harder in GTF file as other columns (e.g., column 9) can have same patterns. The solution is to first select certain columns and then do pattern matching and counting.
3. We will use `cat` to open data stream from a file. Then we will select particular column using `cut`. For example, `cut -f 2` will select column (field) 2.
4. To summarize all gene features present in RefSeq GTF file, we can do this: `cat hg38.ncbiRefSeq.gtf | cut -f 3 | sort | uniq`
5. These features can be counted by adding flag `-c` to `uniq`: `cat hg38.ncbiRefSeq.gtf | cut -f 3 | sort | uniq -c`
6. To limit our summary or counting to a particular chr, say chr22, we can also add column 1 and then `grep` for chr22. Modify the command in step 5: `cat hg38.ncbiRefSeq.gtf | cut -f 1,3 | grep -w `chr22` | sort | uniq -c`
7. So, how many transcripts on chr22?

### Exercise 2: Emsembl GTF file column 3 also has entries for "gene" (missing in RefSeq GTF). Write a command to summarize and count all features including genes in the entire Ensembl GTF file and on chr22. Note the difference in how chromosome numbers are listed in column 1 of the Ensembl GTF file. So, how many genes in the Ensembl GTF? And how many on chr22?
