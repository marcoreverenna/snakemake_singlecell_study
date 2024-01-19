## Imputation pipeline on single-cell SNP array data

<p align="left"> 
    <br> Snakemake workflow for integrating imputation algorithm and processing PGT embryos (SNP array data)
</p>

---

## Table of contents
- [About](#about)
- [Getting Started](#getting_started)
- [Prerequisites and installing](#prerequisites_and_installing)
- [Repository structure](#repository_structure)
- [References](#references)
- [Authors](#authors)
- [Acknowledgments](#acknowledgement)

## About <a name = "about"></a>
The objective of this project is to restore unsatisfactory genotypes using an imputation algorithm. Analyzing single-cell SNP array data poses challenges due to the amplification of DNA, which is necessary because of the limited genetic material available. After quality control, unsuitable genotypes are detected and eliminated. Subsequently, they are imputed again to enhance the accuracy of downstream analysis.

## Getting Started <a name = "getting_started"></a>
This repository run a workflow using [Genomedk](https://genome.au.dk) HPC on Leuven PGDXXX family, integrating an imputation algorithm to enrich the genomic information on our SNP array single-cell data.

## Prerequisites and installing <a name = "prerequisites_and_installing"></a>
This workflow is set up to be executed on genomedk cluster. Therefore, the only prerequisite is to have a genomedk login and to be a member of the meiomap project group on genomedk. Genomedk uses package manager called [Conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html "Conda") and [Git](https://github.com/git-guides/install-git "Git") to download the entire project. Date are confidential; reference panel and map are in the folder on the server but downloadable from [Beagle5.4](https://faculty.washington.edu/browning/beagle/beagle.html).

1. Install dependencies using conda and activate the environment
```
conda env create --name imp_proj --file imp_env.yaml

source activate imp_proj
```
2. Download data for imputation
```
sh get_data.sh
```
3. Process parquet file
```
python processing_parquet.py
```
4. Execute workflow using snakemake command tool
```
snakemake --snakefile Snakefile --cores 10
```

## Repository structure <a name = "repository_structure"></a>
The table below provides an overview of the key files and directories in this repository, along with a brief description of each.
|File  |Description            |
|:----:|-----------------------|
|[Snakefile](Snakefile)|Snakemake main file to run the workflow.|
|[beagle.22Jul22.46e.jar](beagle.22Jul22.46e.jar)|Imputation program developded by Brian Browning.|
|[imp_env.yaml](imp_env.yaml)|Conda requirements to create the environment.|

## References <a name = "references"></a>
- [GenomeDK](https://genome.au.dk/) - Server Environment
- [Beagle5.4](https://faculty.washington.edu/browning/beagle/beagle.html) - Imputation algorithm
- [Snakemake](https://snakemake.readthedocs.io/en/stable/) - Snakemake workflow
  
## Authors <a name = "authors"></a>
- [@marcoreverenna](https://github.com/marcoreverenna) -
- [@ivanvogel](https://github.com/puko818) -
  
## Acknowledgements <a name = "acknowledgement"></a>
I would like to extend my heartfelt gratitude to [KU](https://www.ku.dk/english/) and [CCS](https://ccs.ku.dk)(Center for Chromosome Stability) for providing the essential resources and support that have been fundamental in the development and success of [Eva Hoffmann group](https://icmm.ku.dk/english/research-groups/hoffmann-group/) projects.
