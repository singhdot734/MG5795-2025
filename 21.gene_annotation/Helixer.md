## Using Deep Learning for ab initio gene predictions

#### Helixer GitHub ###
This is the original Helixer that includes model training and prediction  
https://github.com/usadellab/Helixer  
https://github.com/gglyptodon/helixer-docker  

#### HelixerLite GitHub ###
This is a lightweight and easy-to-install Helixer that only includes model prediction  
https://github.com/nextgenusfs/helixerlite  

#### Install HelixerLite ###
```
# set up the base environment
module load miniconda3/24.1.2-py310 
conda env create -f helixer.yml
conda activate helixer

# install with either way
python -m pip install helixerlite
pip install -r helixer_pip.yml

# test if you have helixerlite installed
helixerlite -h
```

### Annotate Maize sequences with different models ###
#### 1. request compute resources ####
```
sinteractive -A PAS3124 -n 48 -t 1:00:00
```

#### 2. update your course folder
```
cd ~/MG5645/name.#/MG5795-2025  
git pull  
```

#### 3. use HelixerLite to predict genes using varying models ####
```
conda activate helixer
cd 21.gene_annotation
nohup helixerlite --fasta B73_chr2_5M.fa --lineage land_plant --out B73_chr2_5M.land_plant.output.gff3 -c 30 &  
nohup helixerlite --fasta B73_chr2_5M.fa --lineage fungi --out B73_chr2_5M.fungi.output.gff3 -c 30 &  
nohup helixerlite --fasta B73_chr2_5M.fa --lineage vertebrate --out B73_chr2_5M.vertebrate.output.gff3 -c 30 &  
nohup helixerlite --fasta B73_chr2_5M.fa --lineage invertebrate --out B73_chr2_5M.invertebrate.output.gff3 -c 30 &
```

### Download and Install IGV on your own computer ###
https://igv.org/download/html/oldtempfixForDownload.html
