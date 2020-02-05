# Installation

Conda is the simplest way to install Panaroo and all its dependencies.

### Conda

```
conda install -c bioconda panaroo
```

For some systems this can install a single threaded version of the cd-hit dependency. In that case the following version can sometimes resolve the issue.

```
conda install -c bioconda/label/cf201901 panaroo
```

### Manual

Install cd-hit using either conda or by following the instructions at https://github.com/weizhongli/cdhit

If you would like to build multiple sequence alignments install the corresponding aligner (MAFFT by default). MAFFT can be installed using conda or by following the instructions at https://mafft.cbrc.jp/alignment/software/source.html

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

Panaroo relies on a number of dependencies that must be installed seperately.

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
