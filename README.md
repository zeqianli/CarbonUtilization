# Predict microbial metabolic traits from genomes

Zeqian Li, Ahmed Selim, Seppe Kuehn

Last update: Mar 24th, 2024

## Files

- `snakemake`: snakemake pipeline for bioinformatic pipeline and fba simulation. See `snakemake/README`.
- `organize_data.ipynb`: organize the three datasets
- `figures.ipynb`: Jupyter notebook to generate figures
    - `figures_*.ipynb:` Jupyter notebook to generate specific panels
- `figures_SI.ipynb`: code to generate SI figures
- `workflows/gapfill`: snakemake pipelines to study the effects of gap-filled carbons. 

## Data

Data files are deposited to the [Open Science Framework (OSF)](https://osf.io/jwkr7/)

## Setup

1. If you use conda environments or pyenv, create and switch to a new environment. 
2. Install packages
```
pip install -r requirements.txt
```
3. Download data
```
mkdir data && cd data
osf -p jwkr7 clone
mv -v data/jwkr7/osfstorage/* .
rm -rf jwkr7
```

`figure_data.zip` contains pre-trained data saved in pickle. 

## Citation

[Li, Zeqian, Ahmed Selim, and Seppe Kuehn. "Statistical prediction of microbial metabolic traits from genomes." PLOS Computational Biology 19, no. 12 (2023): e1011705.](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1011705)

## Contact 

Zeqian Li (zeqianli.chicago@gmail.com)

