# Installation

Conda is the simplest way to install Panaroo and all its dependencies.

### Conda

```
conda install -c bioconda panaroo
```

### Manual

Install cd-hit using either conda or by following the instructions at https://github.com/weizhongli/cdhit

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
