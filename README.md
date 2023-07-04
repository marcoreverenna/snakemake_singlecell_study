<p align="center">
  <a href="" rel="noopener">
 <img width=200px height=200px src="[https://i.imgur.com/6wj0hh6.jpg](https://www.bing.com/images/search?view=detailV2&ccid=LYo6pFwa&id=E704250C7552C84EA19A624E0980256E2022A4F7&thid=OIP.LYo6pFwaa8LnUqM9LrwisQAAAA&mediaurl=https%3a%2f%2fd33wubrfki0l68.cloudfront.net%2f8d66ec8e891f44d59bc223211a41625c11dd1ddb%2f77472%2fimages%2fsnakemake-logo.svg&cdnurl=https%3a%2f%2fth.bing.com%2fth%2fid%2fR.2d8a3aa45c1a6bc2e752a33d2ebc22b1%3frik%3d96QiIG4lgAlOYg%26pid%3dImgRaw%26r%3d0&exph=200&expw=200&q=pipeline+snakemake&simid=607998813628562687&FORM=IRPRST&ck=38BB189A4F6FD8587841CECF200DD61B&selectedIndex=32)" alt="Project logo"></a>
</p>

<h3 align="center">EmbryoHMM</h3>

<div align="center">

[![Status](https://img.shields.io/badge/status-active-success.svg)]()
[![License](https://img.shields.io/badge/license-MIT-blue.svg)](/LICENSE)

</div>

---

<p align="center"> 
    <br> Snakemake workflow for processing PGT embryos (SNP array)
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
git clone https://github.com/Meiomap/ParallelPhasing path/to/workdir
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

- [MongoDB](https://www.mongodb.com/) - Database
- [Express](https://expressjs.com/) - Server Framework
- [VueJs](https://vuejs.org/) - Web Framework
- [NodeJs](https://nodejs.org/en/) - Server Environment

## ✍️ Authors <a name = "authors"></a>

- [@kylelobo](https://github.com/kylelobo) - Idea & Initial work

See also the list of [contributors](https://github.com/kylelobo/The-Documentation-Compendium/contributors) who participated in this project.

## 🎉 Acknowledgements <a name = "acknowledgement"></a>

- Hat tip to anyone whose code was used
- Inspiration
- References
