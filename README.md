# panaroo

[![Travis-CI Build Status](https://travis-ci.com/gtonkinhill/panaroo.svg?token=sFfw198f325BoK14Aor2&branch=master)](https://travis-ci.com/gtonkinhill/panaroo)

An updated pipeline for pan-genome investigation

<p align="center">
<img src="https://github.com/gtonkinhill/panaroo/blob/master/panaroo.png" alt="alt text" width="500">
</p>

## Installation

Prior to python installation [cd-hit](http://weizhongli-lab.org/cd-hit/) is the only external dependency that is required.

```
python setup.py install
```

Then run `panaroo` or `run_prokka`.

If cloning the repository, instead use `python panaroo-runner.py` or `python prokka-runner.py`.

### Dependencies
Required:
* biopython
* numpy
* networkx
* gffutils
* joblib
* edlib
* tdqm
* cd-hit

Optional:
* prodigal
* prokka
* prank
* mafft
* clustal

## Basic usage

To regenerate gene annotations from sequence assemblies using a consistent training model:

```
run_prokka -i *.gff -o reannotated
```

Using these GFFs, or alternatively those from Prokka:

```
panaroo --verbose -i reannotated/*.gff -o results
```
