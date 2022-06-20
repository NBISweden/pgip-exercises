---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13.2
    jupytext_version: 1.13.8
kernelspec:
  display_name: Bash
  language: bash
  name: bash
---


```{code-cell} bash
:"tags": ["remove-input", "remove-output"]
bind 'set enable-bracketed-paste off'
```


(sec_fstatistics)=

# f-statistics #

# Exercises #

````{tab-set-code}


```{code-block} shell
:name: Dsuite
pwd
```

```{code-block} ipython3
:name: sgkit
import sgkit
```

```{code} ipython3
:name: msprime
import msprime
```

````

```{code-cell} bash
pwd
```

NOTE: Dsuite uses allele frequency estimates. Start example by
calculating stuff by hand.


Want also to relate the allele frequencies to branch lengths; cf
Peters
