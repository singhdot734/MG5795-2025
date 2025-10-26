# RNA-Seq analysis
 
## Step 1: Downloading data
1. Load sratoolkit
2. `module load sratoolkit/3.0.2`
3. The following command will download data from a SRA run ID into the same directory where the command will be issued from:
4. `fasterq-dump --split-files <SRR#>`
5. Replace `<SRR#>` with the SRA Run ID to be downloaded.
6. The --split-files option ensures the R1 and R2 files are separate for paired end data.

## Step 2: Trimming adapters
Adapaters to be trimmed will need to be provided in fasta format. 
Create this file in the OnDemand page in the classroom account as follows: 
+New file > In the pop-up window, give it a name "adapters.fa".
To the right of the new file, press the drop-down button to edit it.
In the new tab that opens, paste the following:
```
#adapters.fa
>IlluminaUniversal
AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
```
If there are any other adapters used in samples, add their name and sequence on new lines in fasta format.
Save the file.

The following trimmomatic script can be used to trim adapters these adapters from paired-end samples.

```
  trimmomatic PE \
  SRR#_R1.fastq SRR#_R2.fastq \ #these are the two paired-end input fastq files 
  SRR#_paired_F.fastq SRR#_unpaired_F.fastq SRR#_paired_R.fastq SRR#_unpaired_R.fastq \ #these are four ouput fastq files
  ILLUMINACLIP:adapters.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:18 
```
Replace `SRR#` with file name of the input and output fasta file names. 
As no path is given for each input and output file, the above script assumes input files are in the same directory where the above command is issued. 
The `adapters.fa` file also has to be in this directory.
Output files will also be in this directory.

## Step 3: Alignment
Following items will be needed to do the alignment:
a. Genome index (STAR index), which is in this directory: /fs/ess/PAS3124/MOLGEN_5795_OSU/materials/STAR_index/GRCh38_gencode48_index_STAR_v2_7_11b
b. Trimmed fastq files
c. A shell bash script to submit/queue the alignment job OSC compute nodes. This script is below.
To create a SBATCH script in OnDemand, do as follows: 
+New file > In the pop-up window, give it a name "<my_star_alignment_script.sh>".
To the right of the new file, press the drop-down button to edit it.
In the new tab that opens, paste the following:
```
#!/bin/bash
#SBATCH --account=PAS3124
#SBATCH --time=02:00:00
#SBATCH --mem=35G
#SBATCH --cpus-per-task=8
#SBATCH -o run_star_<sample_id>.log

module load star/2.7.10b

genome_index=/fs/ess/PAS3124/MOLGEN_5795_OSU/materials/STAR_index/GRCh38_gencode48_index_STAR_v2_7_11b
R1_fq=/fs/ess/PAS3124/MOLGEN_5795_OSU/materials/GSE222467/trimmed_data/pc3_siC_1_100k_pf.fastq
R2_fq=/fs/ess/PAS3124/MOLGEN_5795_OSU/materials/GSE222467/trimmed_data/pc3_siC_1_100k_pr.fastq
out_prefix=pc3_siC_1_100k_align_

STAR \
    --readFilesIn $R1_fq $R2_fq \
    --genomeDir $genome_index \
    --outFileNamePrefix $out_prefix \
    --runThreadN 8 \
    --quantMode GeneCounts
```
In the above file, absolute paths are provided for the variables `genome_index`, `R1_fq` and `R2_fq`. Replace these paths as needed where your files are stored.
Run the job submission script `sbatch <my_star_alignment_script.sh>`
To check the progress of your job: `squeue -u <your username>` or `squeue --job <JOBID>`



