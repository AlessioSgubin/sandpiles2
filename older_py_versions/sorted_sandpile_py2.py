# This version is equivalent to the previous one, but written for Python versions 3.4 or lower

import itertools
import copy
import numpy as np

# Functions not supported in Python 3.4 or less...

def merge_two_dicts(x, y):
    z = x.copy()   # start with keys and values of x
    z.update(y)    # modifies z with keys and values of y
    return z



                    #####################################################
                    #           Sorted Sandpile Configuration           #
                    #####################################################

class SandpileSortConfig():

    def __init__(self, sandp, conf, permut, sort = True, verts = []):      ## Defines the class SandpileSortConfig
        r"""
            Definition of the class SandpileSortConfig.
            
            The arguments given are a sandpile and a partition of the vertex set (without sink).
        """

        ### TODO: check if the arguments given correspond to a good SandpileSortConfig

        ### Handle the type of configuration given
        if isinstance(conf, SandpileConfig):  #type: ignore
            config = deepcopy({v:conf[v] for v in sandp.nonsink_vertices()})  #type: ignore
        elif isinstance(conf, list):
            vert = sandp.nonsink_vertices()
            config = deepcopy({vert[i]:conf[i] for i in range(len(vert))})    #type: ignore
        elif isinstance(conf, dict):
            config = deepcopy(conf)                                           #type: ignore
        else:
            print(type(conf))
            raise TypeError

        self.sandpile_struct = sandp                                          # Define the sandpile_struct as sandp
        self.perm_group = permut                                              # Define the permutation group on the sandpile
        if verts == []:                                                       # Define the nonsink vertex set
            self.vertices = sandp.nonsink_vertices()
        else:
            self.vertices = verts
        self.sink = self.sandpile_struct.sink()                               # Define the sink vertex (efficiency!)

        if sort:
            # Reorder the configuration based on the permutation group
            result = []
            pos = {}
            index = 0
            for part in self.perm_group:
                for i in part:
                    pos = merge_two_dicts(pos, {i:index})
                    index += 1
                part_conf = [config[i] for i in part]     # Get the partial configuration associated with an orbit
                part_conf.sort()                          # Sort it
                result = result + part_conf               # Modify the configuration
            sort_conf = {i:result[pos[i]] for i in self.vertices}
        else:
            sort_conf = config

        self.sandpile_config = SandpileConfig(sandp, sort_conf) # type: ignore  # Define the configuration for the sorted sandpile
        

    def __repr__(self):                             ## Returns a description of the class Sandpile2
        return "A sorted configuration for a sandpile with vertices {} and sink {}".format(self.sandpile_struct.vertices(), self.sink)
    

    def sandpile_struct(self):
        r"""
            Returns the underlying sandpile structure.
        """
        return self.sandpile_config.sandpile()

    
    def sort(self, change = True):
        r"""
            Rearrange the configuration values by increasing order on each orbit.

            If change is True, it modifies the configuration, if it is False, it just returns the sorted one.
        """
        result = []
        pos = {}
        index = 0

        for part in self.perm_group:
            for i in part:
                pos = merge_two_dicts(pos, {i:index})
                index += 1
            part_conf = [self.sandpile_config[i] for i in part]     # Get the partial configuration associated with an orbit
            part_conf.sort()                          # Sort it
            result = result + part_conf               # Modify the configuration
        sort_conf = {i:result[pos[i]] for i in self.vertices}
        if change:
            self.sandpile_config = SandpileConfig(self.sandpile_struct, sort_conf) # type: ignore
        return sort_conf
    

    def sort2(self, change = True):
        r"""
            Rearrange the configuration values by increasing order on each orbit.

            If change is True, it modifies the configuration, if it is False, it just returns the sorted one.
        """
        temp_conf = deepcopy(self.sandpile_config)  #type: ignore
        for part in self.perm_group:
            if len(part) > 1:
                temp = copy.copy([temp_conf[i] for i in part])
                temp.sort()
                for i in range(len(part)):
                    temp_conf[part[i]] = temp[i]

        if change:
            self.sandpile_config = temp_conf # type: ignore
        return temp_conf
    

    def __deepcopy__(self, memo):
        r"""
            Overrides the deepcopy method for dict.
        """
        c = SandpileSortConfig(self.sandpile_struct, dict(self.sandpile_config), self.perm_group)
        c.sandpile_config.__dict__.update(self.sandpile_config.__dict__)
        return c


    def __eq__(self, other):                        ## Compares two sorted configurations to see if they are equivalent under the permutation group
        r"""
            Checks if the two sorded configurations are the same or not.
        """
        return self.sandpile_config == other.sandpile_config

    
    def is_recurrent(self):                         ## Check if the configuration is recurrent
        r"""
            Check if the configuration is recurrent.
            
            This function just calls the function is_recurrent() from the class SandpileConfig. 
        """
        return self.sandpile_config.is_recurrent()
    

    def is_stable(self):                            ## Check if the configuration is stable
        r"""
            Check if the configuration is stable.

            This function just calls the function is_stable() from the class SandpileConfig.
        """
        return self.sandpile_config.is_stable()


    def topple_sink(self, sorting = True):                          ## Topples the sink vertex
        r"""
            Topple the sink and change the configuration stored.
        """
        #for vert in self.sandpile_struct.neighbors(self.sandpile_struct.sink()):
        #    self.sandpile_config[vert] += 1
        c = dict(self.sandpile_config)
        for e in self.sandpile_struct.outgoing_edge_iterator(self.sink):
            c[e[1]] += e[2]
        self.sandpile_config = SandpileConfig(self.sandpile_struct, c)    #type: ignore
        if sorting:
            self.sort()
        return self.sandpile_config
    

    def topple_vertex(self, vert, sorting = True):                  ## Topples a vertex of the configuration
        r"""
            Topple a given vertex of the sandpile. This changes the configuration stored.
            
            This function just calls the function fire_vertex() from the class SandpileConfig and applies it to self.sandpile_config. 
        """
        self.sandpile_config = self.sandpile_config.fire_vertex(vert)
        if sorting:
            self.sort()
        return self.sandpile_config
    

    def topple_sequence(self, list, sorting = True):                ## Topples a sequence of vertices
        r"""
            Topple a list of vertices of the sandpile. This changes the configuration stored.
            
            This function calls iteratively the function fire_vertex() from the class SandpileConfig and applies it to self.sandpile_config. 
        """
        for vert in list:
            self.sandpile_config = self.sandpile_config.fire_vertex(vert)
        if sorting:
            self.sort()
        return self.sandpile_config
    

    def __invert__(self, sorting = True):                           ## Topplings until stable
        r"""
            Modifies the Sorted Configuration to get the stable configuration.

            The function calls the SandpileConfig.__invert__ function, then applies the sorting.
        """
        self.sandpile_config = ~self.sandpile_config
        if sorting:
            self.sort()
        return self.sandpile_config


    def level(self):                                ## Returns the level statistic of the configuration
        r"""
            Returns the level of the configuration.
        """
        not_inc_sink = len(self.sandpile_struct.to_undirected().edges()) - self.sandpile_struct.out_degree(self.sink)
        return (- not_inc_sink + self.sandpile_config.deg())
    

    def delay(self, order = [], check_rec = True):  ## Returns the delay statistic of the configuration
        r"""
            Given a reading order of nonsink vertices, this function computes the configuration's delay.

            - order         : the order for reading vertices. If undefined, the decreasing order on vertices is assumed.
            - check_rec     : this option can be used to override the is_recurrent() call.
        """
        if (not self.is_recurrent()) and check_rec:           # Check if the configuration is recurrent
            print(self.sandpile_config)
            print(self.sandpile_config.is_recurrent())
            self.sandpile_struct.show()
            raise Exception("The sorted configuration is not recurrent, hence delay is not defined.")

        nonsink_vert = self.vertices

        if order == []:                                     # No order has been assigned: take decreasing order
            order = copy.copy(nonsink_vert)
            order.sort()
            order.reverse()
        
        delay = 0                                                   # The delay value
        loop_count = 0                                              # The loop count
        not_toppled = copy.copy(nonsink_vert)                                  # A list with the vertices still to be toppled
        self.topple_sink(sorting = False)                                  # Start by toppling the sink
        while len(not_toppled) > 0:
            for ind in order:                               # Search for non-toppled vertices in the right order
                value = self.sandpile_config[ind]
                if (self.sandpile_struct.out_degree(ind) <= value) and (ind in not_toppled):
                                                            # Can be toppled!
                    delay += loop_count                             # Raise the delay by suitable amount
                    self.topple_vertex(ind, sorting = False)        # Topple the vertex without sorting it back!
                    not_toppled.remove(ind)                         # Remove the vertex from not-toppled
            loop_count += 1
        self.sort()
        return delay
    



                    #####################################################
                    #                   Sorted Sandpile                 #
                    #####################################################

class SortedSandpile():

    def __init__(self, graph, sink, permut, opt=[]):    ## Initialize the class SortedSandpile
        r"""
            Definition of the class SortedSandpile.
            
            The arguments given are a sandpile and a partition of the vertex set (without sink).
        """

        ### TODO: check if the arguments given correspond to a good SandpileSortConfig
        
        self.sandpile_struct = Sandpile(graph, sink) # type: ignore                         # Define the Sandpile structure
        self.vertices = self.sandpile_struct.nonsink_vertices()                             # Define the nonsink vertices
        self.perm_group = permut                                                            # Define the permutation group on the graph vertices
        self.specific_opt = opt                                                             # Define the specific options, to optimize in specific cases
        self.sorted_rec = []                                                                # Eventually store the sorted recurrent configurations if computed


    def __repr__(self):                                 ## Description of SortedSandpile class
        return "A sorted sandpile on vertices {} and sink {}.".format(self.sandpile_struct.vertices(), self.sink)
    

    def _max_stable(self):
        r"""
            Returns the maximal stable sorted configuration of the sandpile
        """
        return SandpileSortConfig(self.sandpile_struct, {v:self.sandpile_struct.out_degree(v)-1 for v in self.vertices}, self.perm_group)


    def simple_recurrents(self):                        ## Computes the simple recurrents
        r"""
            Returns the list of recurrent configurations ignoring the action of perm_group.

            It's a call of the recurrents() function from Sandpile class.
        """
        return self.sandpile_struct.recurrents()
    

    def sorted_recurrents(self, option = 0):            ## Computes a list of all sorted recurrent configurations
        r"""
            Computes the sorted recurrent configurations.

            - option        : based on the value the search algorithm changes
                = 0     : calls the recurrent() function from Sandpile class [default]
                = 1     : is a modified version of recurrent() function that ignores the same orbits all together
        """
        if option == 0:                                                         ## OPTION 0: computes all recurrents and then quotient by perm_group action
            simpl_rec = self.simple_recurrents()        # Compute all the recurrent configurations
            for i in range(len(simpl_rec)):             # Reorder elements in each configuration
                temp = SandpileSortConfig(self.sandpile_struct, simpl_rec[i], self.perm_group)
                X = temp.sort()
                simpl_rec[i] = [X[j] for j in list(X)]  # Change dictionary in list
            simpl_rec.sort()                            # Remove duplicates
            simpl_rec = list(simpl_rec for simpl_rec,_ in itertools.groupby(simpl_rec))
            # Now back to a dictionary...
            self.sorted_rec = []
            for i in range(len(simpl_rec)):
                sort_dict = {self.vertices[j]:simpl_rec[i][j] for j in range(len(simpl_rec[i]))}
                self.sorted_rec = self.sorted_rec + [sort_dict]
            return self.sorted_rec

        elif option == 1:                                                         ## OPTION 1: modified version of 0, doesn't append the equivalent configurations
            sorted_temp = []                    # Empty sorted_rec list         
            active = [self._max_stable()]       # Start just with the maximal stable
            while active:                       # While the active list is non-empty...
                c = active.pop()
                sorted_temp.append(c)
                for v in self.vertices:
                    cnext = deepcopy(c) # type: ignore                   # Deepcopy the configuration
                    cnext.sandpile_config[v] += 1                        # Add 1 to a vertex in the configuration
                    cnext.sandpile_config = ~cnext.sandpile_config       # Stabilize the new configuration
                    cnext.sort()
                    if (cnext not in active) and (cnext not in sorted_temp) and (cnext != c):            # If it is still to be discovered and not repeating...
                        active.insert(0,cnext)                           # Instead of appending, put at the start of list!
                # Now convert all SandpileSortConfig to dictionaries...
                self.sorted_rec = [x.sandpile_config for x in sorted_temp]
            return self.sorted_rec

        else:
            raise Exception("The given option is not valid.")
    

    def q_Polynomial(self, ordered = [], opt = 1):      ## Computes the q,t polynomial on (level, delay)
        r"""
            Returns the q - polynomial corresponding to the sorted sandpile's recurrent configurations.

            - order     : if specified, it fixes the reading order for the delay statistic.
            - opt       : if specified, it uses a different computing algorithm
        """
        R = FractionField(QQ['q'])      # type: ignore
        q = R.gen()
        poly = 0*q                        # Define the polynomial as 0

        if self.sorted_rec == []:           # If sorted recurrents still to be computed...
            self.sorted_recurrents(option=opt)        # ...compute them!
        
        # TODO: be sure that the delay doesn't depend on the order in the same orbit...

        if ordered == []:                   # If there is no explicit order check for a specific one
            if self.specific_opt[0] == "clique-indep":      # The reading order that defines delay...
                ordered = self.specific_opt[2]
                
        for config in self.sorted_rec:      # Compute the polynomial
            sortedconfig = SandpileSortConfig(self.sandpile_struct, config, self.perm_group, sort = False, verts = self.vertices)
            q_exp = sortedconfig.level()
            #print("Configurazione {} con livello {} e delay {}".format(sortedconfig.sandpile_config, q_exp, t_exp))
            poly = poly + (q**q_exp)
        return poly


    def qt_Polynomial(self, ordered = [], opt = 1):     ## Computes the q,t polynomial on (level, delay)
        r"""
            Returns the q,t - polynomial corresponding to the sorted sandpile's recurrent configurations.

            - order     : if specified, it fixes the reading order for the delay statistic.
            - opt       : if specified, it uses a different computing algorithm
        """
        R = FractionField(QQ['q, t'])      # type: ignore
        q,t = R.gens()
        poly = 0*q*t                        # Define the polynomial as 0

        if self.sorted_rec == []:           # If sorted recurrents still to be computed...
            self.sorted_recurrents(option=opt)        # ...compute them!
        
        # TODO: be sure that the delay doesn't depend on the order in the same orbit...

        if ordered == []:                   # If there is no explicit order check for a specific one
            if self.specific_opt[0] == "clique-indep":      # The reading order that defines delay...
                ordered = self.specific_opt[2]
                
        for config in self.sorted_rec:      # Compute the polynomial
            sortedconfig = SandpileSortConfig(self.sandpile_struct, config, self.perm_group, sort = False, verts = self.vertices)
            q_exp = sortedconfig.level()
            t_exp = sortedconfig.delay(order = ordered, check_rec=False)
            #print("Configurazione {} con livello {} e delay {}".format(sortedconfig.sandpile_config, q_exp, t_exp))
            poly = poly + (q**q_exp) * (t**t_exp)
        return poly
    
    
    def show(self):                                     ## Function that displays the Sorted Sandpile
        r"""
            This function plots the sorted sandpile.

            If specific_opt is non-empty, the show function uses specific parameters:
            - default                   : call to the Sandpile.show() function
            - clique-independent graph  : call to the Sandpile.show() function fixing the position of the vertices,
                                          the sink in the center and the nonsink in a circle.
        """
        if self.specific_opt == []:                         # Default
            self.sandpile_struct.show()
        elif self.specific_opt[0] == "clique-indep":        # Clique-independent graph
            [mu, nu] = self.specific_opt[1]
            mu_num = sum(mu)
            nu_num = sum(nu)
            # Define the position of each vertex
            positions = merge_two_dicts({0:(0,0)} , {i+1:(-np.sin(2*np.pi*i/(mu_num + nu_num)), np.cos(2*np.pi*i/(mu_num + nu_num)))  for i in range(mu_num + nu_num)} )
            # Define the color of each part
            palette = rainbow(len(mu) + len(nu) + 1) # type: ignore
            col = merge_two_dicts({palette[0]:[0]} , {palette[i+1]:self.perm_group[i] for i in range(len(nu) + len(mu))})

            G = Graph(self.sandpile_struct.dict()).to_undirected() # type: ignore
            
            G.show(pos = positions, vertex_colors = col)
        elif self.specific_opt[0] == "gen-clique-indep":    # Generalized clique-independent graph
            vertex_set = self.specific_opt[1].vertices()
            cliq = []
            indep = []
            for vert in vertex_set:
                if self.specific_opt[2][vert] > 0:
                    cliq.append(vert)
                else:
                    indep.append(vert)
            col = {'blue':cliq, 'red':indep}
            self.specific_opt[1].show(vertex_colors = col)
            



                    ########################################################
                    #                   Specific Sandpiles                 #
                    ########################################################


def CliqueIndependent_SortedSandpile(mu, nu):           ## Specific type of Sandpile
    r"""
        Construction of a Sorted Sandpile on the clique-independent graph given by parameters mu and nu.
    """
    mu_num = sum(mu)                                    # Number of independent vertices
    nu_num = sum(nu)                                    # Number of vertices in cliques
    d = {0 : [i+1 for i in range(mu_num + nu_num)]}     # Initialize the dictionary that will define the graph, link the sink to each vertex
    perm_group = []                                     # Initialize the permutation group acting on the graph

    part_first = 1                                      # Keeps track of first vertex of current part
    for part_nu in nu:                                                  # Add edges for independent vertices. For each part...
        for i in range(part_nu):                                        # ...for each vertex in the part...
            d[part_first + i] = [vert for vert in range(mu_num + nu_num + 1) if (vert < part_first) or (vert >= part_first + part_nu)]
                                                                        # ...add all edges except for other vertices in part_nu    
        perm_group.append([part_first+j for j in range(part_nu)])       # Add the permutation orbit for the nu_part
        part_first += part_nu
    mu_rev = copy.copy(mu)              # We need the reversed partition mu...
    mu_rev.reverse()
    for part_mu in mu_rev:
        for i in range(part_mu):
            d[part_first + i] = [vert for vert in range(mu_num + nu_num + 1) if vert != part_first + i]
        perm_group.append([part_first+j for j in range(part_mu)])       # Add the permutation orbit for the mu_part
        part_first += part_mu

    ordered = []
    for i in range(len(mu)):
        temp = copy.copy(perm_group[len(perm_group)-i-1])
        temp.sort()
        ordered = ordered + temp
    for j in range(len(nu)):
        temp = copy.copy(perm_group[len(nu)-j-1])
        temp.sort()
        temp.reverse()
        ordered = ordered + temp
    
    G = Graph(d)            # type: ignore      # Define the clique-independent graph
    specif_opt = ["clique-indep", [mu, nu], ordered]
    '''
        Specific Options:
        1-  list of the defining partitions
        2-  reading order
    '''
    S = SortedSandpile(G, 0, perm_group, specif_opt)
    return S

def General_CliqueIndependent_SortedSandpile(cells_graph, card_cell, order_cells = []):
    r"""
        Construction of a Sorted Sandpile on the generalized clique-independent graph. Arguments are a graph and a dictionary with values for vertices, such that:
            - nodes:    represents a clique or an independent component with |card_cell[node]| vertices.
                        If the sign of card_cell[node] is:
                            - positive:     we have a clique.
                            - negative:     we have an independent.
            - edges:    each edge correspond to all possible edges between the two cells connected.
    """
    cell_list = cells_graph.vertices()
    
    if order_cells == []:
        order_cells = copy.copy(cell_list)
    
    num_vert = sum([abs(card_cell[i]) for i in cell_list])      # Set the total number of vertices for sandpile

    dict_cell = {}                                              # Dictionary "cell -> list vertices"
    index_next_cell = 1
    for cell in cell_list:
        dict_cell = merge_two_dicts(dict_cell, {cell:[index_next_cell + j for j in range(abs(card_cell[cell]))]})
        index_next_cell += abs(card_cell[cell])
    
    edges_set = {0:[i+1 for i in range(num_vert)]}              # Set of edges, sink to everyone
    perm_group = []                                             # Set the permutation group

    now = 1
    index_next_cell = 1
    for cell in cell_list:                              # Add the cell
        perm_group.append(dict_cell[cell])

        for vert in range(abs(card_cell[cell])):
            if card_cell[cell] > 0:                         # Add for clique
                next_neigh = [0] + [index_next_cell + j for j in range(abs(card_cell[cell])) if index_next_cell + j != now]
            else:                                           # Add nothing for independent
                next_neigh = [0]
            # Add complete edges with other cells...
            for neigh_cell in cells_graph.neighbors(cell):
                next_neigh = next_neigh + dict_cell[neigh_cell]

            # Add to the dictionary
            edges_set = merge_two_dicts(edges_set, {now:next_neigh})
            now += 1

        index_next_cell += abs(card_cell[cell])
    
    G = Graph(edges_set)    #type: ignore
    spec_opt = ["gen-clique-indep", cells_graph, card_cell, []]
    '''
        Specific Options:
        1-  Cell graph
        2-  Cell cardinality dictionary
        3-  Reading order (TODO)
    '''
    S = SortedSandpile(G, 0, perm_group, spec_opt)
    return S
