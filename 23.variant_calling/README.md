# NGS alignment and variant calling

This tutorial steps through some basic tasks in alignment and variant calling using a handful of Illumina sequencing data sets. 
This tutorial is modified from the [alignment-and-variant-calling-tutorial](https://github.com/ekg/alignment-and-variant-calling-tutorial/tree/master).

## Part 0: Setup

We're going to use a bunch of fun tools for working with genomic data:

1. [minimap2](https://github.com/lh3/minimap2)
2. [samtools](https://github.com/samtools/samtools)
3. [freebayes](https://github.com/ekg/freebayes)
4. [vcflib](https://github.com/ekg/vcflib/)
5. [sambamba](https://github.com/lomereiter/sambamba)
6. [vg](https://github.com/vgteam/vg)
7. [vcftools](https://vcftools.github.io/index.html)

In most cases, you can install these tools via conda:

```bash
conda install -c bioconda -c conda-forge bwa 
```

Otherwise, let's assume you're in an environment where you've already got them available.

## Part 1: Aligning maize data with `minimap2`

### Setting up our reference indexes

#### FASTA file index

First, we'll want to allow tools (such as our variant caller) to quickly access certain regions in the reference. This is done using the samtools `.fai` FASTA index format, which records the lengths of the various sequences in the reference and their offsets from the beginning of the file.

```bash
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
    B73_chr2_5M.fa SRR30894834.chr2_1M.fastq \
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
	minimap2 -ax sr -t 10 -R "@RG\tID:$i\tSM:$i" B73_chr2_5M.fa $i.chr2_1M.fastq \
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
freebayes -f B73_chr2_5M.fa --ploidy 2 SRR30894789.sorted.dedup.bam SRR30894834.sorted.dedup.bam SRR30894891.sorted.dedup.bam SRR27984607.sorted.dedup.bam > Zea.chr2_5M.vcf
```

### Hard filtering strategies

The initial variant calling is usually noisy. We can apply some quick filters to increase the calling validity. 

```bash
# a basic filter to remove low-quality sites
vcffilter -f 'QUAL > 10' Zea.chr2_5M.vcf | vt peek -

# scaling quality by depth is like requiring that the additional log-unit contribution
# of each read is at least N
vcffilter -f 'QUAL / AO > 10' Zea.chr2_5M.vcf | vt peek -
```

Another rule of thumb is restricting variant coverage lower than 3 times the median coverage to control for false alignment.

```bash
VCF=Zea.chr2_5M.vcf

# median of site-level DP (use INFO/DP) to define the cap
MED=$(bcftools query -f '%INFO/DP\n' "$VCF" | awk 'NF' | sort -n | \
      awk '{a[NR]=$1} END{ if(NR%2){print a[(NR+1)/2]} else {print (a[NR/2]+a[NR/2+1])/2} }')
MAX=$(awk -v m="$MED" 'BEGIN{printf "%.0f", 3*m}')

# set FORMAT/DP outside [, MAX] to missing; keep others
bcftools filter -e "FORMAT/DP>${MAX}" -S . -Ou "$VCF" \
  | bcftools view -g ^miss -o Zea.chr2_5M.dpflt.vcf
```

### Visualize results
```bash
bcftools stats -s 

Zea.chr2_5M.dpflt.vcfvcftools --vcf Zea.chr2_5M.dpflt.vcf --max-missing 0.75 --recode --recode-INFO-all \
  --out Zea.chr2_5M.miss25 > stats.txt
plot-vcfstats -p vcfstats_plots stats.txt
```
