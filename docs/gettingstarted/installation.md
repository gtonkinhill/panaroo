# Installation

Conda is the simplest way to install Panaurus and all its dependencies.

### Conda

### PyPi

### Manual

Download or clone the repository and install by running
```
git clone https://github.com/gtonkinhill/panaroo
cd panaroo
python3 setup.py install
```

To install locally, instead run

```
python3 setup.py install --user
```

Panaurus relies on a number of dependencies that must be installed seperately.

#### Dependencies
Required:
* biopython
* numpy
* networkx
* gffutils
* joblib
* tdqm
* cd-hit

Optional:
* prokka
* prank
* mafft
* clustal
