---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.12
    jupytext_version: 1.9.1
kernelspec:
  display_name: Bash
  language: bash
  name: bash
---


```{code-cell} bash
:"tags": ["remove-output"]
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

````

```{code-cell} bash
pwd
```


NOTE: Dsuite uses allele frequency estimates. Start example by
calculating stuff by hand.


Want also to relate the allele frequencies to branch lengths; cf
Peters
