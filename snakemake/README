# FBA + annotaiton snakemake pipeline. 

## Installation

1. Installing [snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) and two packages: `click`, `biopython`. Note that installing using mamba is strongly recommended. 

code:
```
mamba create -c conda-forge -c bioconda -n snakemake snakemake click biopython
```

2. Create these conda environments and install the packages:
-  "barrnap": [barrnap](https://github.com/tseemann/barrnap) 
-  "prodigal": [prodigal](https://github.com/hyattpd/Prodigal)
-  "kofamscan": [kofamscan](https://github.com/takaram/kofam_scan)
    - Also download the KEGG database. 
- "fba": carveme, cplex 
    - [CPLEX academic license](https://www.ibm.com/academic/home) is required. Installing via conda does not work. 
        - Note that cplex installs to a root directory by default. I strongly discourage that. Also when setting up cplex's Python API, the default method (running the `setup.py` script) destroyed all my conda environments. I recomment simply changing the `PYTHONPATH` environment variable. 
        - [TODO] After installation, change the code in `rules/build_models.smk` to your own `PYTHONPATH`.

(Ideally, snakemake can install conda environments from .yml files. But I don't want to do that.)

## Run the pipeline

1. Change configurations in `config/config.yml`. See commments in the file. 

2. Dry run for testing: 

```
snakemake --cores all -n
```

3. Start the actual run:

```
snakemake --cores all --conda-frontend conda --use-conda
```
