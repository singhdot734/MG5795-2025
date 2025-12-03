# PacBio CpG calling tutorial
This exercise extracts CpG methylation probability from PacBio HiFi data and performs analyses
Sequences used in this tutorial were derived from:   
1. SRR15447420: [Chen et al. (2023)](https://www.nature.com/articles/s41588-023-01419-6)   
2. SRR11606869: [Hon et al. (2020)](https://www.nature.com/articles/s41597-020-00743-4)   
3. Mo17 assembly: [Chen et al. (2023)](https://www.nature.com/articles/s41588-023-01419-6)   
4. B73 assembly: [Hufford et al. (2021)](https://www.science.org/doi/10.1126/science.abg5289)   

Tools used in this tutorial:   
1. For alignment: [pbmm2](https://github.com/PacificBiosciences/pbmm2)
2. For methylation calling: [pb-CpG-tools](https://github.com/PacificBiosciences/pb-CpG-tools)


## Part 0: Setup

### Request compute resources.
```
sinteractive -A PAS3124 -c 48 -t 1:00:00
```

### Update your course folder
```
cd ~/MG5645/name.#/MG5795-2025
git pull
```

### Activate the `mapping` environment and install tools
```
cd 26.DNA_methylation
module load miniconda3/24.1.2-py310
conda activate mapping
conda install pbmm2 pb-cpg-tools -c bioconda -y
```

## Part 1: Start our analysis 
First, let check out the PacBio BAM file that contains 5mC information
```bash
samtools view B73_chr2_5M.hifi_reads.bam | less
```
Do you know what does each field means? Try to copy the information of a single read and ask AI.


### Align PacBio long reads to the reference
`pbmm2` is a PacBio wrapped version of `minimap2`. In `pbmm2 align`, we use `--preset CCS` to use alignment settings optimized for PacBio HiFi reads.
```bash
module load miniconda3/24.1.2-py310
conda activate mapping
pbmm2 index B73_chr2_5M.fa B73_chr2_5M.mmi  
pbmm2 align B73_chr2_5M.fa B73_chr2_5M.hifi_reads.bam B73_chr2_5M.hifi_reads.aln.bam \
	--sort --num-threads 30 -J 2 --preset CCS

# check out the aligned BAM file.
samtools view B73_chr2_5M.hifi_reads.aln.bam | less
```

### Call CpG methylation probabilities
Signals from a single read will be noisy. We can increase the methylation calling confidence by the accumulation of signals from multiple reads.
```bash
# check out available parameters
aligned_bam_to_cpg_scores --help

# call CpG with selected parameters
aligned_bam_to_cpg_scores --bam B73_chr2_5M.hifi_reads.aln.bam \
	--modsites-mode reference --ref B73_chr2_5M.fa \
	--pileup-mode model --min-mapq 60 --min-coverage 4 \
	--output-prefix B73_chr2_5M.hifi_reads.CpG --threads 30

```
Can you explain each of these parameters?

## Part 2: Check out the CpG calling results
### Check out the methylation calling
Can you explain each field in your data?
```bash
zcat B73_chr2_5M.hifi_reads.CpG.combined.bed.gz | less
```

### Download these files and load onto your IGV
>B73_chr2_5M.fa                      
B73_chr2_5M.EDTA.TEanno.gff3   
B73_chr2_5M.Zm00001eb.1.genes.gff3   
B73_chr2_5M.hifi_reads.aln.bam    
B73_chr2_5M.hifi_reads.aln.bam.bai   
B73_chr2_5M.hifi_reads.CpG.combined.bw          


