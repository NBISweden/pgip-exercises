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


```{code-cell} pgip
:"tags": ["hide-output"]
from bokeh.plotting import figure, show, output_notebook
from pgip import utils
output_notebook()
```


```{note}
- We assume that the input data has already been filtered and pre-processed to obtain a high-quality variant call set. See [](sec_variantfiltering).
- sgkit pca method does not support missing data - add example?
```

(sec_ld_pruning_pca)=

# Principal component analysis

Variant data is high-dimensional and some sort of dimensionality
reduction is often employed to provide a condensed overview of the
data and as a means to identify potential population structure.
Principal component analysis (PCA) is one common technique usually
applied to genetic variation data (e.g.
{cite}`patterson_PopulationStructureEigenanalysis_2006`).

For large data sets, it is common practice to first prune the data to
remove correlation between data points. In this exercise, you will
prune the data using either [plink](sec_pca_with_plink) or
[sgkit](sec_pca_with_sgkit) code snippets.

## Data preparation

```{admonition} FIXME
:class: warning
sgkit is likely too cumbersome? Keep as reference for now, use EIGENSOFT or plink for pruning and pca.

```

```{code-cell} pgip
import os
import sgkit as sg
from sgkit.io import vcf as sgvcf
```
Then make output directories and convert vcf to the zarr format.

```{code-cell} pgip
outdir = "exercises/ld_pruning"
os.makedirs(outdir, exist_ok=True)
vcf = "data/ooa/ooa.gatk.vcf.gz"
zarr = os.path.join(outdir, "ooa.gatk.zarr")
sgvcf.vcf_to_zarr(vcf, zarr)
```


We load the data and remember the original chunk size (needed
internally for sgkit):

```{code-cell} pgip
from myst_nb import glue
ds = sg.load_dataset(zarr)
original_chunk_size = ds.chunks["variants"][0]
glue("nvar", len(ds.variants))
```


## Linkage disequilibrium structure

Before pruning the data, we make a plot of the LD structure. We will
use `python` for plotting and first load the necessary libraries:

## Linkage disequilibrium pruning

Depending on the number of variants, it may be necessary to work on a
subset of variants. Here with only {glue:}`nvar` variants we skip this
step.

For pruning, we first need to calculate the dosage:
```{code-cell} pgip
ds['dosage'] = ds['call_genotype'].sum(dim='ploidy')
```
and then divide into windows by variants:
```{code-cell} pgip
# ds = sg.window_by_variant(ds, size=200)
# print(ds)
```

## PCA

Before proceeding with the pca, we need to compute quality control
variant statistics from genotype calls followed by calculation of and
filtering on allele frequencies.

<!-- See https://github.com/pystatgen/sgkit/issues/752 -->
<!-- All the pca code is based on that example -->

```{code-cell} pgip
ds = sg.variant_stats(ds)

ds = ds.assign(**ds[["variant_allele_frequency"]].compute()).pipe(
lambda ds: ds.sel(variants=(ds.variant_allele_frequency[:, 1] > 0.01))
)
ds = ds.chunk(original_chunk_size)
```

Now we can perform pca.

```{code-cell} pgip
# Transform genotype call into number of non-reference alleles
ds_pca = sg.stats.pca.count_call_alternate_alleles(ds)
# Ensure data has positive variance across samples
variant_mask = ((ds_pca.call_alternate_allele_count < 0).any(dim="samples")) | (
ds_pca.call_alternate_allele_count.std(dim="samples") <= 0.0
)

ds_pca = ds_pca.sel(variants=~variant_mask)
ds_pca["call_alternate_allele_count"] = ds_pca.call_alternate_allele_count.chunk(
(None, -1)
)
ds_pca = sg.pca(ds_pca, n_components=10)
```


We convert the data structure to a pandas dataframe and add various
metadata, such as sample and population names.

```{code-cell} pgip
explained = ds_pca.sample_pca_explained_variance_ratio.values * 100
df = (
ds_pca.sample_pca_projection.to_dataframe()
.reset_index()
.pivot(index="samples", columns="components")
)
df.index.names = ["sample"]
df.columns = [f"PC{i+1}" for _, i in df.columns]
df.reset_index(inplace=True)
df["sample"] = ds_pca.sample_id[df["sample"].values].values
df["population"] = df["sample"].str.slice(0, 3)
colors = dict(zip(set(df["population"]), utils._get_palette(n=3)))
df["color"] = [colors[x] for x in df["population"]]
```

For the purpose of this exercise, we want to also investigate the
effect of sequencing coverage.

```{code-cell} pgip
import io
import pandas as pd
import subprocess
output = subprocess.check_output("samtools depth -a data/ooa/*.bam", shell=True)
coverage = pd.read_table(io.StringIO(output.decode()), header=None).mean(axis=0)
df["coverage"] = coverage[1:]
```

Finally, we can plot the pca with bokeh.
```{code-cell} pgip
from pgip.plotting import bokeh_plot_pca
p = bokeh_plot_pca(df, explained, ncomp=3, ncols=2, width=400, height=400)
show(p)
```
