# Project template
Creator, creation date

## Description

## Project directory description
A project directory structure based on https://drivendata.github.io/cookiecutter-data-science/

## install bioconda
### http://ddocent.com//bioconda/
conda config --add channels r
conda config --add channels defaults
conda config --add channels conda-forge
conda config --add channels bioconda

```
├── LICENSE
├── README.md          <- The top-level README for developers using this project.
├── data
│   ├── derived        <- The final, canonical data sets for modeling.
│   └── raw-data       <- The original, immutable data dump.
│
├── docs               <- A default Sphinx project; see sphinx-doc.org for details
│
├── notebooks          <- Jupyter notebooks. Naming convention is a number (for ordering),
│                         the creator's initials, and a short `-` delimited description, e.g.
│                         `1.0-jqp-initial-data-exploration`.
│
├── references         <- Data dictionaries, manuals, and all other explanatory materials.
│
├── reports            <- Generated analysis as HTML, PDF, LaTeX, etc.
│   └── figures        <- Generated graphics and figures to be used in reporting
│
├── requirements.txt   <- The requirements file for reproducing the analysis environment, e.g.
│                         generated with `pip freeze > requirements.txt`
│
├── requirements_r.txt <- The requirements file for installing non-python packages.
│                         Can be installed with `conda install --file requirements_r.txt`
│
└── src                <- Source code for use in this project.
```
