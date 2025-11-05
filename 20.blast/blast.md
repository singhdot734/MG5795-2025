## Set up your environment ##
### create and actvate a conda environment called "mapping"
### you can name it what ever you like, just something simple and easy to remember
```
module load miniconda3/24.1.2-py310  
conda create -n mapping -y  
conda activate mapping  
```

### install blast to the mapping environment
```
conda install -c bioconda blast  
blastn -h
blastn -help
```

### create a folder called "blast" in your directory
### enter your course folder; name.# is your OSU ID
```
cd ~/MG5795/class/name.#/MG5795-2025/20.blast
```

### create a softlink of the genome
```
ln -s ~/MG5795/class/data/Rice_AGIS1.fasta  
less Rice_AGIS1.fasta  
```

### Explore blast output formats
```
less blast_test1.fa  
blastn -query blast_test1.fa -subject Rice_AGIS1.fasta -outfmt 0  
blastn -query blast_test1.fa -subject Rice_AGIS1.fasta -outfmt 6  
blastn -query blast_test1.fa -subject Rice_AGIS1.fasta -outfmt 7  
blastn -query blast_test1.fa -subject Rice_AGIS1.fasta -outfmt '6 qseqid sseqid pident length qlen qcovs evalue'
```

### run blast with many sequences
```
less B73_chr2_5M.Zm00001eb.1.genes.CDS.fa  
blastn -query B73_chr2_5M.Zm00001eb.1.genes.CDS.fa -subject Rice_MSU7.fasta -outfmt '6 qseqid sseqid pident length qlen qcovs evalue' > B73_chr2_5M.Zm00001eb.1.genes.CDS.fa.out  
```

### filter blast results
```
awk '{if ($6 > 80 && $3 > 70) print $0}' B73_chr2_5M.Zm00001eb.1.genes.CDS.fa.out | less  
```


