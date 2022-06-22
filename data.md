---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13.2
    jupytext_version: 1.13.8
kernelspec:
  display_name: Python 3
  language: python
  name: pgip
---


(sec_data)=

# Data

```{admonition} FIXME
:class: warning
- Describe data sets used in course
```

```{note}
Datasets are presumably of three types:

1. simulated data
2. experimental data
3. precomputed results

Datasets 1 and 2 will be used for analyses. Datasets of type 3 are
results that take a long time to compute (e.g. on a HPC) and will be
needed to enable compilation of jupyter book examples

Only small datasets are committed to this repo. In the case we need
larger datasets, setup an independent data repo.
```

Data sets reside in the `data` subdirectory. Simulated data sets have
been generated with the script `data/make_test_data.py`.

(sec_data_simulated)=

## Simulated data

The script `data/make_test_data.py` has subcommands to simulate
genealogies, mutations, and read data from a demographic model.
Demographic models are defined using the

(sec_data_simulated_ooa)=

### Out of Africa (ooa)

The out of Africa dataset is based on the
{cite}`gutenkunst_InferringJointDemographic_2009` demographic model
(see also [msprime out of Africa
example](https://tskit.dev/msprime/docs/stable/demography.html#population-tree)
and [the stdpopsim out of Africa model
specification](https://popsim-consortium.github.io/stdpopsim-docs/stable/catalog.html#sec_catalog_homsap_models_outofafrica_3g09))


```{code-cell} pgip
:"tags": ["remove-output", "hide-input"]
from myst_nb import glue
import demes
import demesdraw
ooa = demes.load("data/ooa/ooa.demes.yaml")
ax = demesdraw.tubes(ooa, num_lines_per_migration=6, log_time=True, seed=32)
glue("fig-ooa-demesdraw", ax.figure, display=False)
```
```{glue:figure} fig-ooa-demesdraw
:name: fig-ooa-demesdraw

[demesdraw](https://grahamgower.github.io/demesdraw/latest/quickstart.html) illustration of out of Africa model.
```


### Out of Africa with outgroups (ooa-outgroups)

The out of Africa with outgroups is an extension of the previous data
set. In addition to simulating human data, the model has been expanded
to include three outgroups chimpanzee, gorilla and orangutan (Figure
{numref}`fig-ooa-with-outgroups-demesdraw`).


```{code-cell} pgip
:"tags": ["remove-output", "hide-input"]
from myst_nb import glue
import demes
import demesdraw
ooa = demes.load("data/ooa-outgroups/ooa_with_outgroups.demes.yaml")
ax = demesdraw.tubes(ooa, num_lines_per_migration=6, log_time=True, seed=32)
glue("fig-ooa-with-outgroups-demesdraw", ax.figure, display=False)
```
```{glue:figure} fig-ooa-with-outgroups-demesdraw
:name: fig-ooa-with-outgroups-demesdraw

[demesdraw](https://grahamgower.github.io/demesdraw/latest/quickstart.html) illustration of out of Africa model, including outgroups chimpanzee, gorilla, and orangutan.
```
