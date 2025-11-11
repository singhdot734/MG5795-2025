### Request compute resources.
```
sinteractive -A PAS3124 -c 48 -t 1:00:00
```

### Update your course folder
```
cd ~/MG5645/name.#/MG5795-2025  
git pull  
```

### Install EDTA.
Run them line by line
```
cd 22.TE_annotation  
git clone https://github.com/oushujun/EDTA.git  
module load miniconda3/24.1.2-py310

# Install EDTA using either approach
# Approach 1
conda env create -f ./EDTA/EDTA_2.2.x.yml 

# Approach 2
conda create -n EDTA -y  
conda activate EDTA  
conda install -c conda-forge -c bioconda edta -y

```

### Run EDTA.
```
perl ./EDTA/EDTA.pl -h  
cd ./EDTA/test  
perl ../EDTA.pl --genome ./genome.fa --cds ./genome.cds.fa --curatedlib ../database/rice7.0.0.liban --exclude genome.exclude.bed --overwrite 1 --sensitive 1 --anno 1 --threads 48
```
