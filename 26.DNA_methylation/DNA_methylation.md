# alignment and call CpG
```
pbmm2 index B73_chr2_5M.fa B73_chr2_5M.mmi  
pbmm2 align B73_chr2_5M.fa B73_chr2_5M.hifi_reads.bam B73_chr2_5M.hifi_reads.aln.bam --sort -j 4 -J 2 &  
aligned_bam_to_cpg_scores --bam B73_chr2_5M.hifi_reads.aln.bam --output-prefix B73_chr2_5M.hifi_reads.CpG --threads 4  
```

