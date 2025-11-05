### Explore blast output formats
```
blastn -query blast_test1.fa -subject Rice_MSU7.fasta -outfmt 0  
blastn -query blast_test1.fa -subject Rice_MSU7.fasta -outfmt 6  
blastn -query blast_test1.fa -subject Rice_MSU7.fasta -outfmt 7  
blastn -query blast_test1.fa -subject Rice_MSU7.fasta -outfmt '6 qseqid sseqid pident length qlen qcovs evalue'
```

### run blast with many sequences
```
blastn -query B73_chr2_5M.Zm00001eb.1.genes.CDS.fa -subject Rice_MSU7.fasta -outfmt '6 qseqid sseqid pident length qlen qcovs evalue' > B73_chr2_5M.Zm00001eb.1.genes.CDS.fa.out  
```

### filter blast results
```
awk '{if ($6 > 80 && $3 > 70) print $0}' B73_chr2_5M.Zm00001eb.1.genes.CDS.fa.out | less  
```

