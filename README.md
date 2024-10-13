# Sorted Sandpiles and Configurations

## Sorted Sandpile

The class `SortedSandpile` implements a modified version of the Sandpile model.

Given a sandpile $S$ on a graph $G = (V,E)$, we fix a permutation subgroup $\Gamma < \operatorname{Aut}(G)$. In this implementation we assume that $\Gamma$ consists of all possible permutations between subsets of vertices, i.e. $\Gamma \simeq S_{A_1} \times \dots \times S_{A_r}$ where $V = A_1 \sqcup \dots \sqcup A_r$.

The tuple $(S,\Gamma)$ forms a so-called sorted sandpile.


## Sorted Configurations

The class `SandpileSortConfig` is a configuration for a sorted sandpile. It is the analogous of `SandpileConfig` for the class `Sandpile`, with more efficient methods for computing "sorted" recurrent configurations.


# Methods

## Methods for Sorted Sandpile

TBD

## Methods for Sorted Configurations

TBD