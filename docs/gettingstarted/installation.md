# Installation
Conda is the simplest way to install Panaroo and all its dependencies.

### Conda
Although Panaroo supports the latest version of python, currently Bioconda does not support versions >= 3.10 for all packages. An older version of python can be installed by creating an environment as

```
conda create -n panaroo python=3.9
conda activate panaroo
```

Panaroo can then be installed by running the following within a conda environment

```
conda install -c conda-forge -c bioconda -c defaults 'panaroo>=1.3'
```

alternatively it is often faster to use the [mamba](https://github.com/mamba-org/mamba) solver. This can be installed by running

```
conda install mamba -c conda-forge
```

Panaroo can then be installed using

```
mamba install -c conda-forge -c bioconda -c defaults 'panaroo>=1.3'
```

### Manual

Install cd-hit using either conda or by following these [instructions](https://github.com/weizhongli/cdhit).

If you would like to build multiple sequence alignments install the corresponding aligner (MAFFT by default). MAFFT can be installed using conda or by following these [instructions](https://mafft.cbrc.jp/alignment/software/source.html)

Alternatively cd-hit and MAFFT can be installed using apt-get on some linux systems.

```
sudo apt-get install cd-hit

sudo apt-get install mafft
```

Download or clone the repository and install by running

```
pip3 install git+https://github.com/gtonkinhill/panaroo
```
        
or alternatively

```
git clone https://github.com/gtonkinhill/panaroo
cd panaroo
python3 setup.py install
```

To install locally, instead run

```
python3 setup.py install --user
```

Panaroo relies on a number of dependencies that must be installed separately.

#### Dependencies
Required:
* biopython
* numpy
* networkx
* gffutils
* edlib
* joblib
* tdqm
* cd-hit

Optional:
* prokka
* prank
* mafft
* clustal
* mash
