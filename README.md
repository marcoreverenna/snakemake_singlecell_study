<p align="center">
  <a href="" rel="noopener">
 <img width=200px height=200px src="https://i.imgur.com/6wj0hh6.jpg" alt="Project logo"></a>
</p>

<h3 align="center">EmbryoHMM</h3>

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
- [Deployment](#deployment)
- [Usage](#usage)
- [Built Using](#built_using)
- [TODO](../TODO.md)
- [Contributing](../CONTRIBUTING.md)
- [Authors](#authors)
- [Acknowledgments](#acknowledgement)

## 🧐 About <a name = "about"></a>



## 🏁 Getting Started <a name = "getting_started"></a>

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See [deployment](#deployment) for notes on how to deploy the project on a live system.

### Prerequisites


[Git](https://github.com/git-guides/install-git "Git") and [Conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html "Conda") should be installed prior to installation of this workflow. 


### Installing

- Clone workflow into working directory

```
git clone https://github.com/marcoreverenna/snakemake_tutorial path/to/workdir
cd path/to/workdir
```


- install dependencies into isolated environment
conda env create -n myworkflow --file environment.yaml

- activate environment
source activate myworkflow

- execute workflow, i.e.
```
snakemake -c1 "data/crossovers/PGD043_chrom:2_mat_crossovers.csv"
```

## 🔧 Running the tests <a name = "tests"></a>

Explain how to run the automated tests for this system.

### Break down into end to end tests

Explain what these tests test and why

```
Give an example
```

### And coding style tests

Explain what these tests test and why

```
Give an example
```

## 🎈 Usage <a name="usage"></a>

Add notes about how to use the system.

## 🚀 Deployment <a name = "deployment"></a>

Add additional notes about how to deploy this on a live system.

## ⛏️ Built Using <a name = "built_using"></a>
- [GenomeDK](https://genome.au.dk/) - Server Environment
- [Beagle5.4](https://faculty.washington.edu/browning/beagle/beagle.html) - Imputation algorithm
## ✍️ Authors <a name = "authors"></a>
- [@marcoreverenna](https://github.com/marcoreverenna) -
- [@ivanvogel](https://github.com/puko818) -
## 🎉 Acknowledgements <a name = "acknowledgement"></a>

- A
- B
- C
