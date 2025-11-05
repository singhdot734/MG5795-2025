### Request compute resources.
```
sinteractive -A <PAS2387> -c 48 -t 1:00:00
```

### Install EDTA.
```
conda create -n EDTA
conda activate EDTA
mamba install -c conda-forge -c bioconda edta
```

### Run EDTA.
```
git clone https://github.com/oushujun/EDTA.git
cd ./EDTA/test
perl ../EDTA.pl --genome genome.fa --cds genome.cds.fa --curatedlib ../database/rice7.0.0.liban --exclude genome.exclude.bed --overwrite 1 --sensitive 1 --anno 1 --threads 48
```
