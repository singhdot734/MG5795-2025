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
8. 
9. Can you use a combination of `grep` for chr22 and `wc -l` to see how many features are on chr22?
10. Is there a quick way to get an estimate of how many protein coding transcripts are on chr22? Think about what is a feature unique to only protein coding genes/transcripts that is listed in the third field of the GTF file. Then you can combine multiple `grep` commands to look for particular strings and pipe it into line count.
11. So how many (estimated) protein coding transcripts are on chr22? Take a screenshot of your code and answer from line 9 to submit one of the two activity 3 answers.
