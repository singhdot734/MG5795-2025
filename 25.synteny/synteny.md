## This pipeline uses a reference proteome to liftover gene annotations and perform syntenic analysis using McScan.
Assemblies used in this tutorial were derived from:  
1. Arabidopsis thaliana TAIR10: [GCF_000001735.4](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001735.4/)  
2. Arabidopsis halleri: [GCA_964271285.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_964271285.1/)  
3. Capsella bursa-pastoris: [GCA_036452645.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_036452645.1/)  

Tools used in this tutorial:
1. For liftover: [miniprot](https://github.com/lh3/miniprot)  
2. For synteny: [mcscan](https://github.com/tanghaibao/jcvi/wiki/Mcscan-(python-version))


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

### Create and activate the environment.
```
cd 25.synteny
module load miniconda3/24.1.2-py310
conda env create -f synteny.yml 
conda activate synteny
```

### Clone the repo used for syntenic analysis. 
```
git clone https://github.com/cwb14/synLTR.git
```

## Part 1: Synteny analysis

### Get three Arabidopsis genomes. 
```
# Arabidopsis thaliana.
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/735/GCF_000001735.4_TAIR10.1/GCF_000001735.4_TAIR10.1_genomic.fna.gz && gunzip GCF_000001735.4_TAIR10.1_genomic.fna.gz && mv GCF_000001735.4_TAIR10.1_genomic.fna Athal.fa

# Arabidopsis halleri.
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/964/271/285/GCA_964271285.1_ddAraHall_1.3/GCA_964271285.1_ddAraHall_1.3_genomic.fna.gz && gunzip GCA_964271285.1_ddAraHall_1.3_genomic.fna.gz && mv GCA_964271285.1_ddAraHall_1.3_genomic.fna Ahall.fa

# Capsella bursa-pastoris.
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/036/452/645/GCA_036452645.1_ASM3645264v1/GCA_036452645.1_ASM3645264v1_genomic.fna.gz && gunzip GCA_036452645.1_ASM3645264v1_genomic.fna.gz && mv GCA_036452645.1_ASM3645264v1_genomic.fna Cburs.fa
```

### Get the refrence protein file. 
```
# Arabidopsis thaliana.
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/735/GCF_000001735.4_TAIR10.1/GCF_000001735.4_TAIR10.1_protein.faa.gz && gunzip GCF_000001735.4_TAIR10.1_protein.faa.gz && mv GCF_000001735.4_TAIR10.1_protein.faa Athal.pep
```

### Run the syntenic analysis.
```
# activate your conda env if you have not done so
module load miniconda3/24.1.2-py310
conda activate synteny

# Get colinear 1:1 orthologous blocks. 
python synLTR/module1.py --genomes *fa --threads 48 --dir_name results --protein_fa Athal.pep --miniprot_outn 5 --script_dir synLTR/module1/
```

### View results.
```
# Line 1.
cat results/Ahall.Athal.anchors.coords.polished2.consolidated | head -1
Ahall_chr1:2844..49097	Athal_chr1:387478..422154	-
```
*Arabidopsis thaliana* `Athal_chr1:387478-422154` is syntenic to *Arabidopsis halleri* at `Ahall_chr1:2844-49097` on the reverse strand.  

### Download syntenic dotplots from server to personal computer. 
```
ll ./results/dotplots/.pdf   
```

### NCBI houses pre-computed syntenic blocks!
https://www.ncbi.nlm.nih.gov/cgv/browse/GCA_905216605.1/GCF_000001735.4/48945/3702
https://www.ncbi.nlm.nih.gov/cgv/plot/GCA_905216605.1/GCF_000001735.4/48945/3702  
The idiogram and the dotplot are two ways to view synteny. 


## Part 2: coordinate operations

### Install BEDtools
```bash
conda activate synteny
conda install -c bioconda bedtools
```

### Find overlaps between genes and SVs
```bash
bedtools intersect -a <(grep CDS B73_chr2_5M.Zm00001eb.1.genes.gff3) -b B73_chr2_5M.PacBio.SV.vcf -wo|less
```

