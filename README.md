# Sorted Sandpiles and Configurations

In this repository you can find the file `sorted_sandpile.py` that contains the implementation of the new sorted sandpile model.

In the following lines, we redact a more detailed description of the classes and methods for sorted sandpiles and their configurations.

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

-   `nonsink_vertices()`
    The function returns a list of all non-sink vertices of the sorted sandpile.

-   `sink()`
    The function returns the sink vertex of the sorted sandpile.

-   `sandpile()`
    The function returns a `Sandpile` class with the non-sorted structure.

-   `reduced_laplacian()`
    The function returns the reduced laplacian matrix of the underlying sandpile structure.

-   `simple_recurrents(verbose = True)`
    The function computes a list with all simple recurrents for the sandpile structure. It calls the function `Sandpile.recurrents()`.
    Details:
    - INPUT : the arguments are
        - `verbose = True`  : if it is `True` (default) then it returns a list of dictionaries, if `False` a list of lists.
    - OUTPUT: depending on `verbose`, it returns a list of dictionaries (`vertex:value`) or of lists.

-   `sorted_recurrents(option = 0)`
    The function computes a list of all sorted recurrents. Depending on the option, a different algorithm is implemented:
    - `option = 0`  : computes all simple recurrents and then deletes equivalent configurations.
    - `option = 1`  : modified version of `Sandpile.recurrents()` which doesn't append configurations equivalent to those already in the lists.
    - `option = 2`  : modified version of the previous one. The lists are divided up by degree so linear search is more efficient.

    Details:
    - INPUT : the arguments are
        - `option = 0`  : the argument (default is `0`) specify which algorithm to implement for the computation (usually higher value is more efficient with large orbits).
    - OUTPUT: a list of dictionaries with all sorted configurations.

-   `q_Polynomial(opt = 2, override = False)`
    The function computes the q-polynomial for the sorted sandpile. It's the sum over all sorted configurations $c \in \text{SortRec}$ of $q^level(c)$.
    
    Details:
    - INPUT : the arguments are
        - `opt = 2`         : if specified (default is `2`) it chooses the algorithm to compute sorted recurrents.
        - `override = False`: if specified (default is `False`) and is `True` the sorted recurrents are re-computed even if the list is already non-empty.
    - OUTPUT: the q-polynomial, contained in `FractionField(QQ['q'])`.

-   `qt_Polynomial(ordered = [], opt = 2, override = False)`
    The function computes the (q,t)-polynomial for the sorted sandpile. It's the sum over all sorted configurations $c \in \text{SortRec}$ of $q^{\text{level}(c)}t^{\text{delay}(c)}$.
    
    Details:
    - INPUT : the arguments are
        - `ordered = []`    : if specified it defines the reading order for computing delay.
        - `opt = 2`         : if specified (default is `2`) it chooses the algorithm to compute sorted recurrents.
        - `override = False`: if specified (default is `False`) and is `True` the sorted recurrents are re-computed even if the list is already non-empty.
    - OUTPUT: the q-polynomial, contained in `FractionField(QQ['q'])`.

-   `associated_ring(coeff_ring, order = [], homog = False)`
    The function returns a tuple with a polynomial ring and it's toppling ideal.

    Details:
    - INPUT : the arguments are
        - `coeff_ring`      : it is the coefficient ring for the polynomials.
        - `order = []`      : if specified it is the order in which variables $x_i$ are associated to vertices.
        - `homog = False`   : if specified (default is `False`) and `True`, we associate a variable to the sink.
    - OUTPUT: a 2-tuple consisting of the polynomial ring and the ideal associated to the topplings in the sandpile.

-   `symm_poly(coeff_ring, config, order = [])`
    The function returns the simmetric polynomial where each monomial is associated to one of the equivalent configurations under the action that "sorts" the sandpile.

    Details:
    - INPUT : the arguments are
        - `coeff_ring`  : it is the coefficient ring for the polynomials.
        - `config`      : the configuration to which the polynomial is associated.
        - `order = []`  : if specified it is the order in which variables $x_i$ are associated to vertices.
    - OUTPUT: the polynomial symmetric over the vertices associated to the "sorting" of the sandpile.

-   `sortrec_ideal(coeff_ring, sorted = True, opt = 2, override = False, order = [], homog = False)`
    The function computes the ideal in the associated ring containing the polynomials linked to (sorted) configurations.
    _This function is still under construction, since there is no meaningful way of defining the association between configurations and polynomials._

    Details:
    - INPUT : the arguments are
        - `coeff_ring`      : it is the coefficient ring for the polynomials.
        - `sorted = True`   : if specified (default is `True`) and `True` the ring and the ideal are constructed considering the sorted sandpile and not just the sandpile.
        - `opt = 2`         : if specified (default is `2`) it chooses which algorithm to use when calling `SortedSandpile.sorted_recurrents()`.
        - `override = False`: if specified (default is `False`) and `True` the function computes the (sorted recurrents) even if the list `SortedSandpile.sorted_rec` is non-empty.
        - `order = []`      : if specified it is the order in which variables $x_i$ are associated to vertices.
        - `homog = False`   : if specified (default is `False`) and `True` we associate a variable to the sink.
    - OUTPUT: a 3-tuple consisting of the polynomial ring, the toppling ideal and the ideal associated to (sorted) recurrents.

-   `show()`
    This function displays a representation of the sorted sandpile. Depending on the specific options of the class object, the drawing may differ:
    - `default`             : it calls the function `Sandpile.show()` ignoring the "sorted" structure.
    - `clique-indep`        : if the sorted sandpile was created calling `CliqueIndependent_SortedSandpile()`, the representation displays vertices on a wheel (the sink is at the center) and the vertices are placed as defined in the article [[1](https://arxiv.org/abs/2401.06488)]. Colors are used to distinguish different orbits under the action on the sandpile.
    - `gen-clique-indep`    : it displays the cliques and independents as nodes on a graph where edges indicate the relations between components. The color coding is used to differ cliques (in blue) and independents (in red).
    - `mul-clique-indep`    : it is similar to `default` but edges are labeled with their multiplicities.
    - `mulgen-clique-indep` : it is similar to `gen-clique-indep` but edges are labeled with their multiplicities.

    Details:
    - INPUT : the arguments are
        - `default = False` : if specified (default is `False`) and `True`, the function follows the `default` setting even if specific options for the sorted sandpile are present.
    - OUTPUT: none.

-   `export(saveopt = 0, opt = 2)`
    This function returns the critical information of the sorted sandpile in a format that can be saved using pickle.
    If `saveopt = 0` (which is the only implemented option), the function returns a list with:
        -   saveopt
        -   a dictionary of the underlying graph
        -   the sink of the sandpile
        -   the permutation group
        -   the specific options for the sandpile
        -   the order of vertices for...
        -   ...the list of sorted recurrents
        -   the qt-polynomial associated to the sorted sandpile.
    
    Details:
    - INPUT : the arguments are
        - `saveopt = 0` : if specified (default is `0`) it changes the exporting format, according to the description above.
        - `opt = 2`     : if specified (default is `2`) it indicates with which algorithm to compute the sorted recurrents.
    - OUTPUT: a list with all information explained above.

-   `save(namefile, saveopt = 0, opt = 2)`
    This function saves the critical information on the sandpile obtained via `SortedSandpile.extract()` and saves it to the specified location.

    Details:
    - INPUT : the arguments are 
        - `namefile`    : it indicates where to save the file. If `namefile = "PATH/NAME.EXT"` then the informations are saved in `current_dir/PATH/NAME.EXT`.
        - `saveopt = 0` : it is the option `saveopt` passed to `SortedSandpile.export()` function to decide the format for critical information.
        - `opt = 2`     : it is the option `opt` passed to `SortedSandpile.export()` function to compute the list of sorted recurrents.
    - OUTPUT: none.

-   `load(namefile)`
    This function returns the sorted sandpile saved in a given file location.

    Details:
    - INPUT : the arguments are
        - `namefile`    : it indicates which file to read. If `namefile = "PATH/NAME.EXT"` then the function search a file located at `current_dir/PATH/NAME.EXT`.
    - OUTPUT: a class `SortedSandpile` obtained from the given file.


## Methods for Sorted Configurations

The following list contains all methods for class `SandpileSortConfig`.
-   `__init__(sandp, conf, permut, sort = True, verts = [])`
    Initialization of the class SandpileSortConfig.
    
    Details:
    -   INPUT: the possible arguments are:
        -   `sandp`         : an object from `SortedSandpile` class.
        -   `conf`          : it can be a list, dictionary or a `SandpileConfig` object.
        -   `permut`        : this is a list of sub-list. Each sub-list correspond to a subset of all vertices equivalent by sorting.
        -   `sort = True`   : if True, the configuration is sorted using the convention of being increasing on permut sublists.
        -   `verts = []`    : if non-empty this list modify the default naming of vertices.
    -   OUTPUT: a `SandpileSortConfig` implementation with the given input information.

-   `sandpile_struct()`
    Returns the underlying sandpile structure.

-   `sort(change = True)`
    

## Specific Families of Sandpiles

-   `CliqueIndependent_SortedSandpile(mu, nu)`
    This function constructs the original sorted sandpile $G(\mu;\nu)$.
    Details:
    -   INPUT: the arguments are:
        -   `mu`    : a list representing a composition, associated to clique sets.
        -   `nu`    : a list representing a composition, associated to independent sets.
    -   OUTPUT: a `SortedSandpile` obtained using $G(\mu;\nu)$.

-   `General_CliqueIndependent_SortedSandpile(cells_graph, card_cell, order_cells = [])`
    This function constructs a sorted sandpile given a "clique-independent" structure.
    For every vertex `v` of the graph `cells_graph`, the construction considers `card_cell[v]` and
    -   if `card_cell[v] > 0` it creates a clique set of `card_cell[v]` vertices.
    -   if `card_cell[v] < 0` it creates an independent set of `abs(card_cell[v])` vertices.

    For every edge in `cells_graph`, an edge is drawn between all vertices on the two clique/independent sets linked by such edge.
    Finally, add a sink connected to every vertex.
    Details:
    -   INPUT:
        -   `cells_graph`       : a graph structure.
        -   `card_cell`         : a dictionary associating vertices of `cells_graph` to integer values.
        -   `order_cells = []`  : a list containing the reading order of cells in the graph.
    -   OUTPUT: a `SortedSandpile` described by the input structure.

-   `Multi_CliqueIndependent_SortedSandpile(mu, nu, kmul, hmul = -1, sinkmul = 1)`
    This function constructs the sorted sandpile $G(\mu;\nu)$ where edges inside cliques have multeplicity `kmul` and between them multeplicity `hmul`.
    Details:
    -   INPUT: the arguments are:
        -   `mu`        : a list representing a composition, associated to clique sets.
        -   `nu`        : a list representing a composition, associated to independent sets.
        -   `kmul`      : multeplicity of edges inside cliques.
        -   `hmul`      : multeplicity of edges between cliques.
        -   `sinkmul`   : multeplicity of edges between sink and all other vertices.
    -   OUTPUT: a `SortedSandpile` obtained using $G(\mu;\nu)$.

-   `MultiGeneral_CliqueIndependent_SortedSandpile(cells_graph, card_cell, multi_sink = 1, multiedge_cell = {}, order_cells = [])`
    This function constructs a sorted sandpile given a "clique-independent" structure.
    For every vertex `v` of the graph `cells_graph`, the construction considers `card_cell[v]` and
    -   if `card_cell[v] > 0` it creates a clique set of `card_cell[v]` vertices and edges of multeplicity `multiedge_cell[v]`.
    -   if `card_cell[v] < 0` it creates an independent set of `abs(card_cell[v])` vertices.

    For every edge in `cells_graph`, an edge is drawn between all vertices on the two clique/independent sets linked by such edge, with its corresponding multeplicity.
    Finally, add a sink connected to every vertex with multeplicity `multi_sink`.
    Details:
    -   INPUT:
        -   `cells_graph`           : a graph structure.
        -   `card_cell`             : a dictionary associating vertices of `cells_graph` to integer values.
        -   `multi_sink = 1`        : an integer, the multeplicity of edges from sink.
        -   `multiedge_cell = {}`   : list of multeplicities in each clique set.
        -   `order_cells = []`      : a list containing the reading order of cells in the graph.
    -   OUTPUT: a `SortedSandpile` described by the input structure.