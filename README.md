# Sorted Sandpiles and Configurations

## Sorted Sandpile

The class `SortedSandpile` implements a modified version of the Sandpile model.

Given a sandpile $S$ on a graph $G = (V,E)$, we fix a permutation subgroup $\Gamma < \text{Aut}(G)$. In this implementation we assume that $\Gamma$ consists of all possible permutations between subsets of vertices, i.e. $\Gamma \simeq S_{A_1} \times \dots \times S_{A_r}$ where $V = A_1 \sqcup \dots \sqcup A_r$.

The tuple $(S,\Gamma)$ forms a so-called sorted sandpile.


## Sorted Configurations

The class `SandpileSortConfig` is a configuration for a sorted sandpile. It is the analogous of `SandpileConfig` for the class `Sandpile`, with more efficient methods for computing "sorted" recurrent configurations.


# Methods

## Methods for Sorted Sandpile

The following list contains all methods implemented for `SortedSandpile` class.
-   `__init__(graph, sink, permut, opt=[], sort_rec = [])`
    This function initialize the class.
<br>
-   `nonsink_vertices()`
    The function returns a list of all non-sink vertices of the sorted sandpile.
<br>
-   `sink()`
    The function returns the sink vertex of the sorted sandpile.
<br>
-   `sandpile()`
    The function returns a `Sandpile` class with the non-sorted structure.
<br>
-   `reduced_laplacian()`
    The function returns the reduced laplacian matrix of the underlying sandpile structure.
<br>
-   `simple_recurrents(verbose = True)`
    The function computes a list with all simple recurrents for the sandpile structure. It calls the function `Sandpile.recurrents()`.
    Details:
    - INPUT : the arguments are
        - `verbose = True`  : if it is `True` (default) then it returns a list of dictionaries, if `False` a list of lists.
    - OUTPUT: depending on `verbose`, it returns a list of dictionaries (`vertex:value`) or of lists.
<br>
-   `sorted_recurrents(option = 0)`
    The function computes a list of all sorted recurrents. Depending on the option, a different algorithm is implemented:
    - `option = 0`  : computes all simple recurrents and then deletes equivalent configurations.
    - `option = 1`  : modified version of `Sandpile.recurrents()` which doesn't append configurations equivalent to those already in the lists.
    - `option = 2`  : modified version of the previous one. The lists are divided up by degree so linear search is more efficient.

    Details:
    - INPUT : the arguments are
        - `option = 0`  : the argument (default is `0`) specify which algorithm to implement for the computation (usually higher value is more efficient with large orbits).
    - OUTPUT: a list of dictionaries with all sorted configurations.
<br>
-   `q_Polynomial(opt = 2, override = False)`
    The function computes the q-polynomial for the sorted sandpile. It's the sum over all sorted configurations $c \in \text{SortRec}$ of $q^level(c)$.
    
    Details:
    - INPUT : the arguments are
        - `opt = 2`         : if specified (default is `2`) it chooses the algorithm to compute sorted recurrents.
        - 'override = False': if specified (default is `False`) and is `True` the sorted recurrents are re-computed even if the list is already non-empty.
    - OUTPUT: the q-polynomial, contained in `FractionField(QQ['q'])`.
<br>
-   `qt_Polynomial(ordered = [], opt = 2, override = False)`
    The function computes the (q,t)-polynomial for the sorted sandpile. It's the sum over all sorted configurations $c \in \text{SortRec}$ of $q^{\text{level}(c)}t^{\text{delay}(c)}$.
    
    Details:
    - INPUT : the arguments are
        - `ordered = []`    : if specified (default is `[]`) it defines the reading order for computing delay.
        - `opt = 2`         : if specified (default is `2`) it chooses the algorithm to compute sorted recurrents.
        - 'override = False': if specified (default is `False`) and is `True` the sorted recurrents are re-computed even if the list is already non-empty.
    - OUTPUT: the q-polynomial, contained in `FractionField(QQ['q'])`.
<br>
-   `associated_ring`
<br>
-   `symm_poly`
<br>
-   `sortrec_ideal`
<br>
-   `associated_homog_ring`
<br>
-   `sorterec_homog_ideal`
<br>
-   `show`
<br>
-   `export`
<br>
-   `save`
<br>
-   `load`
<br>

## Methods for Sorted Configurations

TBD