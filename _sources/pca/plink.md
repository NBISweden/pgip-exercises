---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13.2
    jupytext_version: 1.13.8
kernelspec:
  display_name: Python 3
  language: ipython3
  name: pgip
---


```{code-cell} pgip
:"tags": ["hide-output", "hide-input"]
from bokeh.plotting import figure, show, output_notebook
from pgip import utils
import subprocess
output_notebook()
```


(sec_pca_with_plink)=

# PCA with plink

```{code-cell} pgip
:"tags": ["hide-input"]
subprocess.run(["plink"])
```


