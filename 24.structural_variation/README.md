# Structural variants calling

This tutorial steps through some basic tasks in alignment and structural variant calling using some maize data sets.
This tutorial is modified from the [Krumlov_teaching_SV](https://github.com/fritzsedlazeck/teaching_material/blob/main/Krumlov/Krumlov_teaching_SV.md).

## Goals of this module
The goal of this module is to get you familiarized with structural variant identification across assembly, short read based mapping and long read based mapping.          
For Structural Variation (SV) detection we will be using multiple methods (Assemblytics, minimap2, Sniffles) and then later compare them across to obtain more insights about advantages and disadvantages across the different the differnet approaches.

### IMPORTANT NOTES
All the instructions below requires you to think along. You are responsible for a clear data structure on your workspace and to find the individual files that you need. I am happy to assist but please take a moment to look around (`ll`) or think where a certain file could be.

### The main steps in this Module are:
1. Assembly based SV detection (using [Assemblytics](http://assemblytics.com/))
2. Long read based mapping based SV detection (using [Sniffles](https://github.com/fritzsedlazeck/Sniffles))
3. SV comparison (using [SURVIVOR](https://github.com/fritzsedlazeck/SURVIVOR))

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

### Install packages

We're going to use a bunch of fun tools for working with genomic data. You can install these tools via conda. Let's create a new environment called "sv" because the packages used in today's analysis are in conflict with the "mapping" environment we previously created.

```bash
module load miniconda3/24.1.2-py310
conda create -n sv -c bioconda -c conda-forge survivor "sniffles>2" assemblytics minimap2 bcftools samtools
conda activate sv
```


## Part 1: Illumina based small InDels detection 
Let's first look at the small Insertions and Deletions (InDels) detected through aligning illumina short reads. These Indels are usually shorter than 50 bp due to read length limitations and alignment uncertainties.

Identify short Indels from the short-read-based VCF that was created from the SNP lecture.
```bash
bcftools view -v indels Zea.chr2_5M.filtered.annotated.vcf > B73_chr2_5M.Illumina.InDel.vcf

# check it out
less -S B73_chr2_5M.Illumina.InDel.vcf

# how many SV we could identify
grep -vc '#' B73_chr2_5M.Illumina.InDel.vcf
```


## Part 2: Assembly based SV detection 
As discussed in the lecture, assembly based SV detection is quite comprehensive. Nevertheless, before we start, we need to align a new assembly to an existing e.g. reference genome to identify structural differences. For this we will be using MUMmer (`nucmer`). Alternatively, you could also utilize `Dipcall` as a program especially if you have a phased assembly at hand. MUMmer is a commonly used package that allows you to rapidly compare two sequences together and includes multiple packages for summary reports over the aligned sequences. This includes variant reporting, coordinate reporting or even the generation of dot plots. As shown in the lecture dot plot is very helpful and intuitive method to compare sequences for us.

To begin to call variants we need to first compare the reference to our assembly. Now initiate the alignment. It will create an alignment file named `out.delta` for our next step.
```
nucmer -maxmatch -l 100 -c 500 B73_chr2_5M.fa Mo17_chr2_5M.fa
```

Use Assemblytics to identify SV candidates with minimum variant size of 50 bp and maximum variant size of 50 kb.
```bash
Assemblytics out.delta B73-Mo17_chr2_5M 100 50 50000
```

Let's check out the summary of SVs detected by Assemblytics:
```bash
less B73-Mo17_chr2_5M.Assemblytics_structural_variants.bed
less B73-Mo17_chr2_5M.Assemblytics_structural_variants.summary
```

To convert the Assemblytics SV into VCF format we will need SURVIVOR. Minimum variant size is set at 10 bp.
```
SURVIVOR convertAssemblytics B73-Mo17_chr2_5M.Assemblytics_structural_variants.bed 10 B73_chr2_5M.WGA.SV.vcf
less -S B73_chr2_5M.WGA.SV.vcf
```

Calculate the average length of SVs
```bash
# SV length without filtering
grep -v '#' B73-Mo17_chr2_5M.Assemblytics_structural_variants.bed | awk '{LEN += $5} END {print LEN/NR}'

# SV length with filtering
grep -v '#' B73_chr2_5M.WGA.SV.vcf | sed 's/.*SVLEN=//; s/;PE.*//' | awk '{LEN+=$0} END {print LEN/NR}'
```


## Part 3: Long read based SV detection 
Finally we are ready for detecting SVs using PacBio long reads. First, let's align the long reads to the reference.

```bash
# Map Mo17 long reads to the B73 reference sequence
minimap2 -ax map-hifi -t 48 -R '@RG\tID:Mo17\tSM:Mo17' B73_chr2_5M.fa Mo17_chr2_5M.fastq.gz | samtools sort - > B73_chr2_5M.Mo17.bam

# Map B73 long reads to the B73 reference sequence
minimap2 -ax map-hifi -t 48 -R '@RG\tID:B73\tSM:B73' B73_chr2_5M.fa B73_chr2_5M.fastq.gz | samtools sort - > B73_chr2_5M.B73.bam
```

Using Sniffles v2 to identify SVs that are >= 10 bp with at least 3 reads supporting. Filtering out unreliable read alignments with mapping quality 40 help to control the calling confidence.
Here we perform multi-sample SV calling
```bash
# index both bam files with one command
samtools index -M B73_chr2_5M.Mo17.bam B73_chr2_5M.B73.bam

# Create .snf for each sample
sniffles --input B73_chr2_5M.B73.bam --snf B73_chr2_5M.B73.snf
sniffles --input B73_chr2_5M.Mo17.bam --snf B73_chr2_5M.Mo17.snf

# Joint SV calling
sniffles --input B73_chr2_5M.Mo17.snf B73_chr2_5M.B73.snf \
	--vcf B73_chr2_5M.PacBio.SV.vcf \
	--reference B73_chr2_5M.fa \
	--minsvlen 10 --minsupport 3 --mapq 40
```


Next we can inspect the file with e.g.:
```
less -S B73_chr2_5M.PacBio.SV.vcf
grep -vc '#' B73_chr2_5M.PacBio.SV.vcf
```

How many SV did you detect? 


## Part 4: Structural Variant comparison
Now that we generated Assembly, Illumina, and PacBio based SV calls, it is time to compare them. One tool that you can use for this very easily is SURVIVOR. SURVIVOR is a very simple method to compare SV but also includes multiple other methods that might be useful.

For SURVIVOR, we want to use the merge option. Before doing this, the merge option requires a file including all paths and VCF files that you want to compare. Thus, we generate the file like this:
```
ls B73_chr2_5M.Illumina.InDel.vcf > vcf_files
ls B73_chr2_5M.WGA.SV.vcf >> vcf_files
ls B73_chr2_5M.PacBio.SV.vcf >> vcf_files
```

Next we can initiate the compare with 100bp wobble and requiring that we are only merging with SV type agreement. Furthermore, we will only take variants into account with 1bp+. 
```bash
SURVIVOR merge vcf_files 100 1 1 0 0 0 B73_chr2_5M.merged.SV.vcf
```
Lets check how good/bad the overlap is:
```bash
perl -ne 'print "$1\n" if /SUPP_VEC=([^,;]+)/' B73_chr2_5M.merged.SV.vcf | sort | uniq -c 
```

As you can see you will get the pattern and the number of times the pattern occurs. The first number is the number of times it can be observed in the VCF file. The 2nd number is the pattern (0 or 1 depending if it was observed or not) 

