---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.12
    jupytext_version: 1.9.1
kernelspec:
  display_name: Python 3
  language: python
  name: pgip
---

(sec_intro)=

# Population Genomics in Practice exercises

Welcome to the Population Genomics in Practice exercises homepage!

(sec_exercises_about)=

## About the exercises

As the focus of the course is on hands-on work, the topics have been
designed to cover the fundamental analyses that are common in many
population genomics studies. 

(sec_exercises_about_manuscript_route)=

### The manuscript route

In principle, you could imagine the course structure to follow that of
a manuscript (e.g {cite}`fuller_PopulationGeneticsCoral_2020`).

High-throughput DNA sequencing has now made it possible to generate
whole-genome resequencing data for multiple individuals and
populations, and a first step is to map sequence data to a reference,
perform variant calling and variant filtering.

Once a high-quality variant set has been obtained, a common task is to
describe variation, either in terms of summary statistics such as
nucleotide diversity ($\pi$) or site-frequency spectra (sfs), or as
descriptions of population structure in terms of admixture or pca
plots.

Genetic diversity is also affected by population history and
demographic processes such as population expansion, bottlenecks,
migration events and hybridizations.

Finally, it is often of interest to identify adaptive traits, to which
end selection tests and scans can be performed. The tests are designed
to detect signals of selection, either via direct selection on loci,
or by looking at haplotype structures to detect linked selection.

(sec_exercises_about_baseline)=

### The baseline model

Much of what has been described in [The manuscript
route](sec_exercises_about_manuscript_route) has recently been treated
in an article on statistical inference in population genomics
{cite}`johri_StatisticalInferencePopulation_2021`. In it, the authors
point out that whereas historically theoretical advances outpaced data
production, that is no longer true due to the advent of
next-generation sequencing. In particular, they caution researchers to
attach too much faith to a test that explains the data well, as there
are many alternative hypotheses with equal explanatory power, but with
drastically different conclusions. At the very least, a population
genomics study should aim at first generating a baseline model
consisting of all or several of the following components:

1. mutation
2. recombination
3. reassortment
4. gene conversion
5. purifying selection acting on functional regions and its effects on
   linked variants (background selection)
6. genetic drift with demographic history and geographic structure

The exercises are designed to address many of the points above, and to
highlight cases where competing hypotheses may actually explain data
to equal degrees.

(sec_exercises_about_additional_material)=

### Additional material

```{admonition} FIXME
:class: warning
Rewrite/restructure and add more additional topics
```


```{note}
NB: with time we could use repo to collect all kinds of exercises
(even those that don't fit current curriculum), providing the
possibility to create optional course plans
``` 

A five day course cannot but provide an overview of possible analyses
and topics. Among the things we won't have time to cover are

1. machine learning and AI for population genomics
2. recombination landscape estimation
3. gene conversion
4. experimental design
5. coalescent and Wright-Fisher simulations
6. ABC methods
7. ancestral recombination graphs and tree sequence inference
8. ancient DNA
9. spatial genetics
10. ...


