## This pipeline uses a reference proteome to liftover gene annotations and perform syntenic analysis using McScan.
For liftover: https://github.com/lh3/miniprot  
For synteny: https://github.com/tanghaibao/jcvi/wiki/Mcscan-(python-version)


### Initiate interactive job.
```
sinteractive -A <PAS2387> -c 48 -t 0:30:00
```

### Create and activate the environment.
```
mamba env create -f Genomic-Data-Analysis/synteny/synteny.yml 
conda activate synteny
```
### Clone the repo used for syntenic analysis. 
```
git clone https://github.com/cwb14/synLTR.git
```
### Get the genomes. 
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
# Get colinear 1:1 orthologous blocks. 
python synLTR/module1.py --genomes *fa --threads 48 --dir_name results --protein_fa Athal.pep --miniprot_outn 5 --script_dir synLTR/module1/
```

# View results.
```
# Line 1.
cat results/Ahall.Athal.anchors.coords.polished2.consolidated | head -1
Ahall_chr1:2844..49097	Athal_chr1:387478..422154	-
```
*Arabidopsis thaliana* `Athal_chr1:387478-422154` is syntenic to *Arabidopsis halleri* at `Ahall_chr1:2844-49097` on the reverse strand.  


### Download syntenic dotplots from server to personal computer. 
```
scp -r username@remote.server.address:/path/to/remote/file/results/dotplots/\*pdf /path/to/local/destination
```

### NCBI houses pre-computed syntenic blocks!
https://www.ncbi.nlm.nih.gov/cgv/browse/GCA_905216605.1/GCF_000001735.4/48945/3702
https://www.ncbi.nlm.nih.gov/cgv/plot/GCA_905216605.1/GCF_000001735.4/48945/3702  
The idiogram and the dotplot are two ways to view synteny. 
