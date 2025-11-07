### Request compute resources.
```
sinteractive -A PAS3124 -c 48 -t 2:00:00
```

### Install EDTA.
```
module load miniconda3/24.1.2-py310
conda create -n EDTA -y
conda activate EDTA
conda install -c conda-forge -c bioconda edta
```

### Update your course folder
```
cd ~/MG5645/name.#/MG5795-2025  
git pull  
```

### Run EDTA.
```
cd 22.TE_annotation  
git clone https://github.com/oushujun/EDTA.git  
./EDTA/EDTA.pl -h  
cd ./EDTA/test  
perl ../EDTA.pl --genome genome.fa --cds genome.cds.fa --curatedlib ../database/rice7.0.0.liban --exclude genome.exclude.bed --overwrite 1 --sensitive 1 --anno 1 --threads 48  
```
