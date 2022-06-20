---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13.2
    jupytext_version: 1.13.8
kernelspec:
  display_name: bash
  language: bash
  name: bash
---

(sec_nucleotide_diversity)=

# Nucleotide diversity
The nucleotide diversity is defined as the average number of pairwise
differences per site for a number of sequences
{cite}`nei_MolecularEvolutionPhylogenetics_2000`, p. 251:
[@chan_GenomeWideFineScaleRecombination_2012]
{cite}`chan_GenomeWideFineScaleRecombination_2012`
$$
 \pi = \sum_{ij} x_i x_j \pi_{ij}
$$

Here $x_i$ is the population frequency of sequence $i$. It is unclear
whether it means the number of occurrences or the fraction of sequence
$i$ in a population of $n$ samples.

Since $\pi_{ii} = 0$ and $\pi_{ij} = \pi_{ji}$, the equation can be
rewritten as

$$
 \pi = \sum_{ij} x_i x_j \pi_{ij} = \sum_{i<j} x_i x_j \pi_{ij} + \sum_{j<i} x_j x_i \pi_{ji} + \sum_{i} x_i x_i \pi_{ii} = 2\sum_{ij} x_i x_j \pi_{ij}
$$

The (second) summation has $n(n-1)/2$ terms.

The nucleotide diversity can also be expressed as a function of the
allele frequency spectrum (afs) $p = (p_1,p_2,...,p_q)$, where
$p_i$ is the frequency of sites with $i$ minor alleles, and $q$ is the
number of possible non-monomorphic allele configurations. For diploid
organisms, $q = 2k - 1$. $p$ is the *unfolded* afs. There is a
*folded* version $p^* = (p^*_1, p^*_2, p^*_{\lceil(q+1)/2\rceil})$,
where $p^*_i = p_i + p_{q-i}$ for $i \neq \lceil(q+1)/2\rceil$ and
$p^*_i = p_i$ otherwise.

Now, for a given site $s$ with $i$ minor alleles, there are $i(q-i)$
different combinations where $\pi_{ij} \neq 0$. Moreover, $x_i = x_j =
 1$, such that

$$
 \pi_s = \sum_{i<j} x_i x_j \pi_{ij} = i(q - i)
$$

Generalizing to $S = \sum_{i=1}^{q} p_i$ sites, there are $p_i$ sites
with allele configuration $i, (q-i)$, such that

$$
 \pi = \sum_{i=1}^{q} p_i i (q-i)
$$

Normalizing by the number of terms for each site, $n(n-1)/2$ gives the
average pairwise nucleotide diversity per site.

Letting $k = \lceil(q+1)/2\rceil$, we can further rewrite $\pi$ as

$$
 \pi = \sum_{i=1}^{q} p_i i (q-i) = p_k k (q-k) +  \sum_{i=1}^{k - 1} p_i i (q-i) + p_{q-i} (q-i) i = \sum_{i=1}^{k} p^*_i i (k - i)
$$
 
