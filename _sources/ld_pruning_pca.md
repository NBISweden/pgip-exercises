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
from bokeh.plotting import figure, show, output_notebook
output_notebook()
```


```{note}
- We assume that the input data has already been filtered and pre-processed to obtain a high-quality variant call set. See [](sec_variantfiltering).
- sgkit pca method does not support missing data - add example?
```

(sec_ld_pruning_pca)=

# Ld pruning and principal component analysis

Variant data is high-dimensional and some sort of dimensionality
reduction is often employed to provide a condensed overview of the
data and as a means to identify potential population structure.
Principal component analysis (PCA) is one common technique usually
applied to genetic variation data (e.g.
{cite}`patterson_PopulationStructureEigenanalysis_2006`).

For large data sets, it is common practice to first prune the data to
remove correlation between data points. In this exercise, you will
prune the data using either `plink` or `sgkit` code snippets.

## Data preparation

```{admonition} FIXME
:class: warning
sgkit is likely too cumbersome to explain to include. Keep as reference, use EIGENSOFT or plink for pruning and pca.

```

Before pruning the data, we make a plot of the LD structure. We will
use `python` for plotting and first load the necessary libraries:


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
from bokeh.palettes import Set3
from bokeh.plotting import figure


def _get_palette(cmap=Set3[12], n=12, start=0, end=1):
    import matplotlib
    import numpy as np

    linspace = np.linspace(start, end, n)
    cmap = matplotlib.colors.LinearSegmentedColormap.from_list("customcmap", cmap)
    palette = cmap(linspace)
    hex_palette = [matplotlib.colors.rgb2hex(c) for c in palette]
    return hex_palette

```



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

and plot the first two components

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
colors = dict(zip(set(df["population"]), _get_palette(n=3)))
df["color"] = [colors[x] for x in df["population"]]
```


using specialized functions
```{code-cell} pgip
# NB: https://speciationgenomics.github.io/pca/ calculates explained
# variance from sums of plink eigenvals, although not all components
# have been calculated
import re
from bokeh.io import output_file
from bokeh.layouts import gridplot

from bokeh.models import ColumnDataSource
import itertools
import math

def bokeh_plot_pca_coords(df, explained, *, pc1=1, pc2=2, **kw):
    source = ColumnDataSource(df)
    x = explained[pc1 - 1]
    y = explained[pc2 - 1]
    xlab = f"PC{pc1} ({x:.1f}%)"
    ylab = f"PC{pc2} ({y:.1f}%)"

    metadata_columns = [x for x in df.columns if not re.match("^(PC[0-9]+|color)$", x)]
    tooltips = [(x, f"@{x}") for x in metadata_columns]
    p = figure(
        x_axis_label=xlab,
        y_axis_label=ylab,
        tooltips=tooltips,
        title=f"PC{pc1} vs PC{pc2}",
        **kw,
    )
    p.circle(
        x=f"PC{pc1}",
        y=f"PC{pc2}",
        source=source,
        color="color",
        size=15,
        alpha=0.8,
        line_color="black",
        legend_group="population",
    )
    p.add_layout(p.legend[0], "right")
    return p


def bokeh_plot_pca(df, eigenvals, ncomp=6, filename=None, **kw):
    pairs = list(itertools.combinations(range(ncomp), 2))
    n = len(pairs)
    ncols = math.floor(math.sqrt(n))
    plots = []
    for (i, j) in pairs:
        p = bokeh_plot_pca_coords(df, eigenvals, pc1=i + 1, pc2=j + 1, **kw)
        plots.append(p)
    gp = gridplot(plots, ncols=ncols)
    if filename is not None:
        output_file(filename)
        show(gp)
    else:
        return gp

```


```{code-cell} pgip
p = bokeh_plot_pca(df, explained)
show(p)
```
