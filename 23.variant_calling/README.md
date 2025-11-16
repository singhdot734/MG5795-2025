# NGS alignment and variant calling

This tutorial steps through some basic tasks in alignment and variant calling using a handful of Illumina sequencing data sets. 
This tutorial is modified from the [alignment-and-variant-calling-tutorial](https://github.com/ekg/alignment-and-variant-calling-tutorial/tree/master).

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

We're going to use a bunch of fun tools for working with genomic data:

1. [minimap2](https://github.com/lh3/minimap2)
2. [samtools](https://github.com/samtools/samtools)
3. [freebayes](https://github.com/ekg/freebayes)
4. [vcflib](https://github.com/ekg/vcflib/)
5. [sambamba](https://github.com/lomereiter/sambamba)
6. [vt](https://github.com/atks/vt)
7. [vcftools](https://vcftools.github.io/index.html)
8. [bcftools](https://github.com/samtools/bcftools)
9. [snpEff](https://pcingola.github.io/SnpEff/)

You can install these tools via conda:

```bash
module load miniconda3/24.1.2-py310
conda activate mapping
conda install -c bioconda -c conda-forge minimap2 samtools freebayes vcflib sambamba vcftools tectonic pdflatex bcftools vt matplotlib snpeff
```

## Part 1: Aligning maize data with `minimap2`

### Setting up our reference indexes

#### FASTA file index

First, we'll want to allow tools (such as our variant caller) to quickly access certain regions in the reference. This is done using the samtools `.fai` FASTA index format, which records the lengths of the various sequences in the reference and their offsets from the beginning of the file.

```bash
cd 23.variant_calling
samtools faidx B73_chr2_5M.fa
```

Now it's possible to quickly obtain any part of the B73_chr2_5M.fa reference sequence. For instance, we can get the 200bp from position 1000000 to 1000200. We'll use a special format to describe the target region `[chr]:[start]-[end]`.

```bash
samtools faidx B73_chr2_5M.fa B73_chr2:1000000-1000200
```

We get back a small FASTA-format file describing the region:

```text
>B73_chr2:1000000-1000200
ACTGGGGCAAGGGGTGAAAGCGTTGAACTCTCGAACAACCAGAACGCTGACATTTTTGGCCTTTTAAAAGAATTCCACTATCATATCCAGAAATATACCAATGGGAGCTAACGCATAAGAAACGAAGAATTAAGACAGCGTAAAAACAACCGCTTGCTGTACAGGTCACCTGCCACCAGCAAAGACAGACATGTACATCAG
```


### Aligning our data against the B73_chr2_5M.fa reference

Here's an outline of the steps we'll follow to align our B73 reads against the B73 reference:

1. use minimap2 to generate SAM records for each read
2. convert the output to BAM
3. sort the output
4. mark PCR duplicates that result from exact duplication of a template during amplification

We could do the steps one-by-one, generating an intermediate file for each step.
However, this isn't really necessary unless we want to debug the process, and it will make a lot of excess files which will do nothing but confuse us when we come to work with the data later.
Thankfully, it's easy to use [unix pipes](https://en.wikiepdia.org/wiki/Pipeline_%28Unix%29) to stream most of these tools together (see this [nice thread about piping bwa and samtools together on biostar](https://www.biostars.org/p/43677/) for a discussion of the benefits and possible drawbacks of this).

You can now run the alignment using a piped approach. _Replace `$threads` with the number of CPUs you would like to use for alignment._ 
We will use `minimap2` to align the reads to the reference genome. In the case of these short reads, we'd use it as follows.

```bash
minimap2 -ax sr -t 10 -R '@RG\tID:SRR30894834\tSM:SRR30894834' \
    B73_chr2_5M.fa SRR30894834.chr2_1M.fastq.gz \
    | samtools view -b - > SRR30894834.raw.bam
sambamba sort SRR30894834.raw.bam
sambamba markdup SRR30894834.raw.sorted.bam SRR30894834.raw.sorted.dedup.bam
```

Breaking it down by line:

- *alignment with minimap2*: `minimap2 -t $threads -R '@RG\tID:SRR30894834\tSM:SRR30894834'` --- this says "align using so many threads" and also "give the reads the read group SRR30894834 (the NCBI data ID for B73) and the sample name SRR30894834". We also tell it we're working with short read data using `-ax sr`.
- *reference and FASTQs* `B73_chr2_5M.fa SRR30894834.chr2_1M.fastq` --- this just specifies the base reference file name and the input alignment files.
- *conversion to BAM*: `samtools view -b -` --- this reads SAM from stdin (the `-` specifier in place of the file name indicates this) and converts to BAM.
- *sorting the BAM file*: `sambamba sort SRR30894834.bam` --- sort the BAM file, writing it to `.sorted.bam`.
- *marking PCR duplicates*: `sambamba markdup SRR30894834.sorted.bam SRR30894834.sorted.dedup.bam` --- this marks reads which appear to be redundant PCR duplicates based on their read mapping position. It [uses the same criteria for marking duplicates as picard](http://lomereiter.github.io/sambamba/docs/sambamba-markdup.html).


Now, run the same alignment process for the other three samples. Make sure to specify a different sample name via the `-R '@RG...` flag incantation to specify the identity of the data in the BAM file header and in the alignment records themselves. Run the following commands line by line:

```bash
for i in SRR27984607 SRR30894789 SRR30894891; do 
	minimap2 -ax sr -t 10 -R "@RG\tID:$i\tSM:$i" B73_chr2_5M.fa $i.chr2_1M.fastq.gz \
		| samtools view -b - > $i.raw.bam & \
done

for i in SRR27984607 SRR30894789 SRR30894891; do sambamba sort $i.raw.bam & done

for i in SRR27984607 SRR30894789 SRR30894891; do sambamba markdup $i.raw.sorted.bam $i.raw.sorted.dedup.bam & done

```

## Part 2: Calling variants

Now that we have our alignments sorted, we can quickly determine variation against the reference by scanning through them using a variant caller.
There are many options, including [samtools mpileup](http://samtools.sourceforge.net/samtools.shtml), [platypus](http://www.well.ox.ac.uk/platypus), and the [GATK](https://www.broadinstitute.org/gatk/).

For this tutorial, we'll keep things simple and use [freebayes](https://github.com/ekg/freebayes). It has a number of advantages in this context (bacterial genomes), such as long-term support for haploid (and polyploid) genomes. However, the best reason to use it is that it's very easy to set up and run, and it produces a very well-annotated VCF output that is suitable for immediate downstream filtering.

### Joint calling with `freebayes`

It's quite easy to use `freebayes` provided you have your BAM file completed. We use `--ploidy 2` to indicate that the sample should be genotyped as diploid.

We can put the samples together if we want to find differences between them. Calling them jointly can help if we have a population of samples to use to help remove calls from paralogous regions. The Bayesian model in freebayes combines the data likelihoods from sequencing data with an estimate of the probability of observing a given set of genotypes under assumptions of neutral evolution and a [panmictic](https://en.wikipedia.org/wiki/Panmixia) population. For instance, [it would be very unusual to find a locus at which all the samples are heterozygous](https://en.wikipedia.org/wiki/Hardy%E2%80%93Weinberg_principle). It also helps improve statistics about observational biases (like strand bias, read placement bias, and allele balance in heterozygotes) by bringing more data into the algorithm.

Another reason to call them jointly is to make sure we have a genotype for each one at every locus where a non-reference allele passes the caller's thresholds in each sample. We would run a joint call by dropping in all BAMs on the command line to freebayes:

```bash
freebayes -f B73_chr2_5M.fa --ploidy 2 *sorted.dedup.bam > Zea.chr2_5M.vcf
```

### Hard filtering strategies

The initial variant calling is usually noisy. We can apply some quick filters to increase the calling validity. 

```bash
# a basic filter to remove low-quality sites
vcffilter -f 'QUAL > 10' Zea.chr2_5M.vcf | vt peek -

# scaling quality by depth is like requiring that the additional log-unit contribution
# of each read is at least N
vcffilter -f 'QUAL / AO > 10' Zea.chr2_5M.vcf | vt peek -

# Apply stricter filters:
# 	Keep only biallelic sites;
#	Require at least 2 alt reads (AO >= 2)
# 	Require quality-per-alt-read (QUAL/AO > 10)
# 	Require no missing genotypes (all samples called)
NS=4 # four samples
bcftools view -m2 -M2 -i "INFO/AO>=2 && QUAL/INFO/AO>10 && COUNT(GT!='mis')==$NS" -o Zea.chr2_5M.filtered.vcf Zea.chr2_5M.vcf
```

### Annotate SNP effects
```bash
# Create a genome folder in snpEff/data
mkdir -p snpEff/data/MyGenome
cp B73_chr2_5M.fa snpEff/data/MyGenome/sequences.fa
cp ../21.gene_annotation/B73_chr2_5M.Zm00001eb.1.genes.gff3 snpEff/data/MyGenome/genes.gff

# Create and add entry to snpEff.config
echo "MyGenome.genome : MyGenome" > snpEff/snpEff.config

# Build database for this genome
snpEff build -gff3 -noCheckCds -noCheckProtein -v MyGenome -c snpEff/snpEff.config

# Annotate your SNP VCF
snpEff ann MyGenome Zea.chr2_5M.filtered.vcf -c snpEff/snpEff.config > Zea.chr2_5M.filtered.annotated.vcf

# Check out the results
less snpEff_genes.txt

# You need to download snpEff_summary.html and open with your browser for the annotation summary
```

### Output variant table
```bash
# Compress the VCF file and index it using bcftools
bgzip Zea.chr2_5M.filtered.vcf
bcftools index Zea.chr2_5M.filtered.vcf.gz

# Convert the VCF file to a tabulated variant table
bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/DP\t[%SAMPLE\t%GT\t%DP\t]\n' Zea.chr2_5M.filtered.vcf.gz -o Zea.chr2_5M.filtered.txt
```

### Visualize results
```bash
# Summarize the filter VCF file. Outputs are stored in the vcfstats_plots folder.
bcftools stats -s - Zea.chr2_5M.filtered.vcf.gz > stats.txt
plot-vcfstats -p vcfstats_plots stats.txt

```

