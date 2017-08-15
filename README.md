# pheno-calour
[Calour](https://github.com/amnona/Calour) interface to [Phenotype database](https://doi.org/10.6084/m9.figshare.4272392)

Based on:

[Hiding in Plain Sight: Mining Bacterial Species Records for Phenotypic Trait Information](http://msphere.asm.org/content/2/4/e00237-17)

Albert Barberán, Hildamarie Caceres Velazquez, Stuart Jones, Noah Fierer

DOI: 10.1128/mSphere.00237-17

## Installation:
1. Install Calour:

Follow instructions [Here](https://github.com/amnona/Calour)

2. Install pheno-calour:

```
pip install git+git://github.com/amnona/pheno-calour
```

## Usage:
Easiest way to use calour is via the [EZCalour](https://github.com/amnona/EZCalour) GUI.

Alternatively, Calour can be used via Jupyter notebook/ipython.

In order to obtain phenotype annotations for V4 region (forward primer 515F) [deblurred](https://github.com/biocore/deblur) sequences, add:

databases = ['phenotype', ...] to the plot / plot_sort functions


If phenotype data is available for the selected sequence, it will be shown in the details box.

double clicking on the annotation will open the paper from where the phenotypes were obtained.
