## Imputation pipeline on single-cell SNP array data

<div align="center">

[![Status](https://img.shields.io/badge/status-active-success.svg)]()
[![License](https://img.shields.io/badge/license-MIT-blue.svg)](/LICENSE)

</div>

---

<p align="center"> 
    <br> Snakemake workflow for integrating imputation algorithm and processing PGT embryos (SNP array)
</p>

## 📝 Table of Contents

- [About](#about)
- [Getting Started](#getting_started)
- [Prerequisites and installing](#prerequisites_and_installing)
- [Deployment](#deployment)
- [Usage](#usage)
- [Built Using](#built_using)
- [Authors](#authors)
- [Acknowledgments](#acknowledgement)

## 🧐 About <a name = "about"></a>
The objective of this project is to restore unsatisfactory genotypes using an imputation algorithm. Analyzing single-cell SNP array data poses challenges due to the amplification of DNA, which is necessary because of the limited genetic material available. After quality control, unsuitable genotypes are detected and eliminated. Subsequently, they are imputed again to enhance the accuracy of downstream analysis.

## 🏁 Getting Started <a name = "getting_started"></a>
These instructions will enable you to have a copy of the project up and running on your genomedk profile for development, testing and playback of the workflow for all Leuven PGD036 family.

### 🔧 Prerequisites and installing <a name = "prerequisites_and_installing"></a>
This workflow is set up to be executed on genomedk cluster. Therefore, the only prerequisite is to have a genomedk login and to be a member of the meiomap project group on genomedk. Genomedk uses package manager called [Conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html "Conda") and [Git](https://github.com/git-guides/install-git "Git") to download the entire project.

1. Clone workflow into working directory
```
git clone https://github.com/marcoreverenna/snakemake_tutorial path/to/workdir
cd path/to/workdir
```
2. install dependencies into isolated environment
```
conda env create -n myworkflow --file environment.yaml
```
4. activate environment
```
source activate myworkflow
```
5. execute workflow, i.e.
```
snakemake --snakefile snakemake_one_sample.smk --cores 10
```

## 🚀 Deployment <a name = "deployment"></a>
...............................................

## 🎈 Usage <a name="usage"></a>
...............................................
Add notes about how to use the system.

Add additional notes about how to deploy this on a live system.
## ⛏️ Built Using <a name = "built_using"></a>
- [GenomeDK](https://genome.au.dk/) - Server Environment
- [Beagle5.4](https://faculty.washington.edu/browning/beagle/beagle.html) - Imputation algorithm
## ✍️ Authors <a name = "authors"></a>
- [@marcoreverenna](https://github.com/marcoreverenna) -
- [@ivanvogel](https://github.com/puko818) -
## 🎉 Acknowledgements <a name = "acknowledgement"></a>
...............................................
