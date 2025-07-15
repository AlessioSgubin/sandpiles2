import itertools
import copy
import pickle
import numpy as np


                    #####################################################
                    #                Auxiliary Functions                #
                    #####################################################

def is_symmetric(poly):         # Check if polynomial is symmetric in q and t
    R = FractionField(QQ['q','t','x'])  # type: ignore
    q,t,x = R.gens()
    poly1 = R(poly)     # Casting in new ring
    poly2 = poly1
    poly2(q=x)
    poly2(t=q)
    poly2(x=t)
    return (poly == poly2)

def is_increasing_list(lst):
    r"""
        Checks if a list is weakly increasing.
    """
    i = 0
    check = True
    while i < len(lst)-1 and check:
        check = (lst[i]<=lst[i+1])
        i += 1
    return check


                    #####################################################
                    #           Sorted Sandpile Configuration           #
                    #####################################################


class SandpileSortConfig():

    def __init__(self, sandp, conf, permut, sort = True, verts = []):      ## Defines the class SandpileSortConfig
        r"""
            Initialization of the class SandpileSortConfig.
            
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
                    pos = pos | {i:index}
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
                pos = pos | {i:index}
                index += 1
            part_conf = [self.sandpile_config[i] for i in part]     # Get the partial configuration associated with an orbit
            part_conf.sort()                          # Sort it
            result = result + part_conf               # Modify the configuration
        sort_conf = {i:result[pos[i]] for i in self.vertices}
        if change:
            self.sandpile_config = SandpileConfig(self.sandpile_struct, sort_conf) # type: ignore
        return sort_conf
    

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


    def is_vertex_stable(self, v):                  ## Check if a vertex of the configuration is stable
        r"""
            Check if the vertex of the configuration is stable.
        """
        neigh = self.sandpile_struct.edges(vertices=v)
        outdeg_v = sum(x[2] for x in neigh)
        return self.sandpile_config[v] < outdeg_v


    def is_sorted(self):                            ## Check if the configuration has already been sorted
        r"""
            Check if the configuration has already been sorted.
        """
        for perm in self.perm_group:
            conf = [self.sandpile_config[i] for i in perm]
        return is_increasing_list(conf)


    def deg(self):                                  ## Returns the degree of the configuration
        r"""
            Returns the degree of the configuration.
        """
        return self.sandpile_config.deg()


    def topple_sink(self, sorting = True):          ## Topples the sink vertex
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
    

    def single_topple(self, vert, threshold = 0, sorting = True):   ## Toppling sending one grain for each arc
        r"""
            Topple a vertex by sending a grain to each edge.
        """
        neigh = self.sandpile_struct.edges(vert)
        for v in neigh:
            if (v[1] != self.sink) and (v[2] > threshold):
                self.sandpile_config[v[1]] += 1
                self.sandpile_config[vert] -= 1
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


    def level(self):                                                ## Returns the level statistic of the configuration
        r"""
            Returns the level of the configuration.
        """
        not_inc_sink = sum([v[2] for v in self.sandpile_struct.to_undirected().edges()]) - self.sandpile_struct.out_degree(self.sink)
        return (- not_inc_sink + self.sandpile_config.deg())
    

    def delay(self, order = [], check_rec = True):                  ## Returns the delay statistic of the configuration
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
    

    def k_delay(self, k = -1, order = [], check_rec = True):                ## Returns the NEW delay statistic of the configuration
        r"""
            Given a reading order of nonsink vertices, this function computes the configuration's new conjectured delay.
            - k             : the multeplicity of all edges (except incident to the sink).
            - order         : the order for reading vertices. If undefined, the decreasing order on vertices is assumed.
            - check_rec     : this option can be used to override the is_recurrent() call.
        """
        if (not self.is_recurrent()) and check_rec:           # Check if the configuration is recurrent
            print(self.sandpile_config)
            print(self.sandpile_config.is_recurrent())
            self.sandpile_struct.show()
            raise Exception("The sorted configuration is not recurrent, hence delay is not defined.")

        if k == -1:
            k = max([edge[2] for edge in self.sandpile_struct.edges()])

        nonsink_vert = self.vertices
        n = len(nonsink_vert)

        if order == []:                                     # No order has been assigned: take decreasing order
            order = copy.copy(nonsink_vert)
            order.sort()
            order.reverse()

        toppl = [0]*n                                           # Stores information on how many partial topplings are still needed
        finalv = [k]*n                                          # We exit the loop when toppl == finalv
        delay = 0
        plus = 0
        self.topple_sink(sorting = False)                               # Start by toppling the sink
        while toppl != finalv:                               # Until everything has been toppled k times...
            for i in range(len(order)):
                if (self.sandpile_struct.out_degree(order[i]) <= self.sandpile_config[order[i]]) and (toppl[i] == 0):
                                                                        # Can be toppled for the first time!
                    self.single_topple(order[i], threshold = k-1-toppl[i], sorting = False)
                    toppl[i] += 1
                    delay += plus
                else:
                    if toppl[i] < finalv[i] and toppl[i] > 0:                   # Topple later time...
                        self.single_topple(order[i], threshold = k-1-toppl[i], sorting = False)
                        toppl[i] += 1
            plus += 1
        #print("The configuration {} has k-delay {}".format(self.sandpile_config, delay))
        self.sort()
        return delay
    

    def k_delay_test(self, k = -1, order = [], info = False, check_rec = True):                ## Returns the NEW delay statistic of the configuration
        r"""
            Given a reading order of nonsink vertices, this function computes the configuration's new conjectured delay.
            - k             : the multeplicity of all edges (except incident to the sink).
            - order         : the order for reading vertices. If undefined, the decreasing order on vertices is assumed.
            - check_rec     : this option can be used to override the is_recurrent() call.
        """
        if (not self.is_recurrent()) and check_rec:           # Check if the configuration is recurrent
            print(self.sandpile_config)
            print(self.sandpile_config.is_recurrent())
            self.sandpile_struct.show()
            raise Exception("The sorted configuration is not recurrent, hence delay is not defined.")

        if k == -1:
            k = max([edge[2] for edge in self.sandpile_struct.edges()])

        nonsink_vert = self.vertices
        n = len(nonsink_vert)

        if order == []:                                     # No order has been assigned: take decreasing order
            order = copy.copy(nonsink_vert)
            order.sort()
            order.reverse()

        toppl = [0]*n                                           # Stores information on how many partial topplings are still needed
        finalv = [k]*n                                          # We exit the loop when toppl == finalv
        delay = 0
        plus = 0
        wtopp = []      # Toppling word
        record = []     # Record of iterations
        self.topple_sink(sorting = False)                               # Start by toppling the sink
        latex_table = "0" + "".join([" & {}".format(self.sandpile_config[n - v]) for v in range(n)] + ["\\\\ \n"])
        while toppl != finalv:                               # Until everything has been toppled k times...
            for i in range(len(order)):
                if (self.sandpile_struct.out_degree(order[i]) <= self.sandpile_config[order[i]]) and (toppl[i] == 0):       # LATEX CODE
                                                                        # Can be toppled for the first time!
                    self.single_topple(order[i], threshold = k-1-toppl[i], sorting = False)
                    toppl[i] += 1
                    delay += plus
                    wtopp = wtopp + [order[i]]
                    record = record + [order[i]]
                    latex_table += "{}".format(order[i])
                    for ind in range(n):                                                                                    # LATEX CODE
                        if self.is_vertex_stable(n - ind):
                            latex_table += " & {}".format(self.sandpile_config[n - ind])
                        else:
                            latex_table += " & \\textbf{{ {} }}".format(self.sandpile_config[n - ind])
                    latex_table += "\\\\ \n"
                else:
                    if toppl[i] < finalv[i] and toppl[i] > 0:                   # Topple later time...
                        self.single_topple(order[i], threshold = k-1-toppl[i], sorting = False)
                        toppl[i] += 1
                        wtopp = wtopp + [order[i]]
                        record = record + [order[i]]
                        latex_table += "{}".format(order[i])
                        for ind in range(n):                                                                                    # LATEX CODE
                            if self.is_vertex_stable(n - ind):
                                latex_table += " & {}".format(self.sandpile_config[n - ind])
                            else:
                                latex_table += " & \\textbf{{ {} }}".format(self.sandpile_config[n - ind])
                        latex_table += "\\\\ \n"
                    else:
                        record = record + [-1]
                        #latex_table += "\\textbullet"
                        #for ind in range(n):                                                                                    # LATEX CODE
                        #    if self.is_vertex_stable(n - ind):
                        #        latex_table += " & {}".format(self.sandpile_config[n - ind])
                        #    else:
                        #        latex_table += " & \\textbf{{ {} }}".format(self.sandpile_config[n - ind])
                        #latex_table += "\\\\ \n"
            plus += 1
        self.sort()
        if info:
            return [delay, wtopp, record, latex_table]
        else:
            return delay
  

    def show(self, sink = True, colors = False, heights = False, directed = False):     ## Returns a drawing of the configuration
        r"""
            Returns a drawing of the sorted sandpile configuration.
        """
        # Define the color of each part
        palette = rainbow(len(self.perm_group)+1) # type: ignore
        if sink:
            col = {palette[0]:[self.sink]} | {palette[i+1]:self.perm_group[i] for i in range(len(self.perm_group))}
        else:
            col = {palette[i+1]:self.perm_group[i] for i in range(len(self.perm_group))}
        
        # Call SandpileConfig.show()
        self.sandpile_config.show(sink = sink, colors = colors, heights = heights, directed = directed, vertex_colors = col)
    



                    #####################################################
                    #                   Sorted Sandpile                 #
                    #####################################################

class SortedSandpile():

    def __init__(self, graph, sink, permut, opt=[], sort_rec = []):    ## Initialize the class SortedSandpile
        r"""
            Definition of the class SortedSandpile.
            
            The arguments given are a sandpile and a partition of the vertex set (without sink).
        """

        ### TODO: check if the arguments given correspond to a good SandpileSortConfig
        
        self.sandpile_struct = Sandpile(graph, sink) # type: ignore                         # Define the Sandpile structure
        self.vertices = self.sandpile_struct.nonsink_vertices()                             # Define the nonsink vertices
        self.perm_group = permut                                                            # Define the permutation group on the graph vertices
        self.specific_opt = opt                                                             # Define the specific options, to optimize in specific cases
        self.sorted_rec = sort_rec                                                          # Eventually store the sorted recurrent configurations if computed


    def __repr__(self):                                 ## Description of SortedSandpile class
        return "A sorted sandpile on vertices {} and sink {}.".format(self.sandpile_struct.vertices(), self.sandpile_struct.sink())
    

    def _max_stable(self):                              ## Returns maximal stable configuration
        r"""
            Returns the maximal stable sorted configuration of the sandpile.
        """
        return SandpileSortConfig(self.sandpile_struct, {v:self.sandpile_struct.out_degree(v)-1 for v in self.vertices}, self.perm_group)
    
    
    def nonsink_vertices(self):                         ## Returns non-sink vertices
        r"""
            Returns the non-sink vertices of the sorted sandpile.
        """
        return self.sandpile_struct.nonsink_vertices()
    

    def sink(self):                                     ## Returns sink
        r"""
            Returns the sink of the sorted sandpile.
        """
        return self.sandpile_struct.sink()
    

    def sandpile(self):                                 ## Returns the unsorted sandpile
        r"""
            Returns the sandpile without sorting.
        """
        return self.sandpile_struct

    
    def reduced_laplacian(self):                        ## Reduced Laplacian
        r"""
            Returns the reduced Laplacian matrix of the sandpile.
        """
        return self.sandpile_struct.reduced_laplacian()


    def simple_recurrents(self, verbose = True):                            ## Computes the simple recurrents
        r"""
            Returns the list of recurrent configurations ignoring the action of perm_group.

            It's a call of the recurrents() function from Sandpile class.
        """
        return self.sandpile_struct.recurrents(verbose = verbose)
    

    def sorted_recurrents(self, option = 2):                                ## Computes a list of all sorted recurrent configurations
        r"""
            Computes the sorted recurrent configurations.

            - option        : based on the value the search algorithm changes
                = 0     : calls the recurrent() function from Sandpile class [default].
                = 1     : is a modified version of recurrent() function that ignores the same orbits all together.
                = 2     : modified version of the previous, where the duplicate search is more efficient since we sort by configuration degree.
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
        
        elif option == 2:                                                         ## OPTION 2: modified version of 1 that sorts the recurrent configurations and active by degree (invariant for sorting).
                                                                                  #            The idea is to make the search in active and sorted_temp more efficient.
            max_deg = self._max_stable().deg()
            sorted_temp = []                                                                    # Empty sorted_rec list     
            sorted_dict = {i:[] for i in range(max_deg+1)}                                      # Empty sorted_rec dictionary, keyed by degree   
            active = [self._max_stable()]                                                       # Empty active list
            active_dict = {max_deg:[self._max_stable()]} | {i:[] for i in range(max_deg)}       # Empty active dictionary, keyed by degree
            while active:                       # While the active list is non-empty...
                c = active.pop()
                deg_c = c.deg()
                active_dict[deg_c].pop()
                sorted_temp.append(c)               # Append to sorted list
                sorted_dict[deg_c].insert(0,c)      # Append to sorted list dictionary
                for v in self.vertices:
                    cnext = deepcopy(c) # type: ignore                   # Deepcopy the configuration
                    cnext.sandpile_config[v] += 1                        # Add 1 to a vertex in the configuration
                    cnext.sandpile_config = ~cnext.sandpile_config       # Stabilize the new configuration
                    cnext.sort()
                    deg_cnext = cnext.deg()
                    if (cnext not in active_dict[deg_cnext]) and (cnext not in sorted_dict[deg_cnext]) and (cnext != c):            # If it is still to be discovered and not repeating...
                        active.insert(0,cnext)
                        active_dict[deg_cnext].insert(0,cnext)
            # Now convert all SandpileSortConfig to dictionaries...
            self.sorted_rec = [x.sandpile_config for x in sorted_temp]
            return self.sorted_rec

        else:
            raise Exception("The given option is not valid.")
    

    def q_Polynomial(self, opt = 2, override = False):                      ## Computes the q,t polynomial on (level, delay)
        r"""
            Returns the q - polynomial corresponding to the sorted sandpile's recurrent configurations.

            - order     : if specified, it fixes the reading order for the delay statistic.
            - opt       : if specified, it uses a different computing algorithm.
            - override  _ if specified, the sorted recurrents are computed even if sorted_rec is non-empty.
        """
        R = FractionField(QQ['q'])      # type: ignore
        q = R.gen()
        poly = 0*q                                  # Define the polynomial as 0

        if self.sorted_rec == [] or override:       # If sorted recurrents still to be computed...
            self.sorted_recurrents(option=opt)          # ...compute them!
        
        for config in self.sorted_rec:      # Compute the polynomial
            sortedconfig = SandpileSortConfig(self.sandpile_struct, config, self.perm_group, sort = False, verts = self.vertices)
            q_exp = sortedconfig.level()
            #print("Configurazione {} con livello {} e delay {}".format(sortedconfig.sandpile_config, q_exp, t_exp))
            poly = poly + (q**q_exp)
        return poly


    def qt_Polynomial(self, ordered = [], opt = 2, override = False):       ## Computes the q,t polynomial on (level, delay)
        r"""
            Returns the q,t - polynomial corresponding to the sorted sandpile's recurrent configurations.

            - order     : if specified, it fixes the reading order for the delay statistic.
            - opt       : if specified, it uses a different computing algorithm.
            - override  : if specified, the sorted recurrents are computed even if sorted_rec is non-empty.
        """
        R = FractionField(QQ['q, t'])      # type: ignore
        q,t = R.gens()
        poly = 0*q*t                            # Define the polynomial as 0

        if self.sorted_rec == [] or override:   # If sorted recurrents still to be computed...
            self.sorted_recurrents(option=opt)        # ...compute them!
        
        # TODO: be sure that the delay doesn't depend on the order in the same orbit...

        if ordered == []:                   # If there is no explicit order check for a specific one
            if self.specific_opt[0] == "clique-indep" or self.specific_opt[0] == "mul-clique-indep":          # The reading order that defines delay...
                ordered = self.specific_opt[2]
            elif self.specific_opt[0] == "gen-clique-indep" or self.specific_opt[0] == "mulgen-clique-indep":
                ordered = self.specific_opt[3]
            
        if self.specific_opt[0] == "mul-clique-indep":              # Compute the polynomial with k_delay
            for config in self.sorted_rec:
                sortedconfig = SandpileSortConfig(self.sandpile_struct, config, self.perm_group, sort = False, verts = self.vertices)
                q_exp = sortedconfig.level()
                t_exp = sortedconfig.k_delay(order = ordered, check_rec=False)
                #print("Configurazione {} con livello {} e delay {}".format(sortedconfig.sandpile_config, q_exp, t_exp))
                poly = poly + (q**q_exp) * (t**t_exp)
        else:                                                       # Compute the polynomial with regular delay
            for config in self.sorted_rec:
                sortedconfig = SandpileSortConfig(self.sandpile_struct, config, self.perm_group, sort = False, verts = self.vertices)
                q_exp = sortedconfig.level()
                t_exp = sortedconfig.k_delay(order = ordered, check_rec=False)
                #print("Configurazione {} con livello {} e delay {}".format(sortedconfig.sandpile_config, q_exp, t_exp))
                poly = poly + (q**q_exp) * (t**t_exp)
        
        return poly


    def associated_ring(self, coeff_ring, order = [], homog = False):       ## Compute the associated ring
        r"""
            Compute the ring associated to the sorted sandpile. The possible arguments are:
                - coeff_ring    : specify the coefficient ring for the polynomial ring.
                - order         : the vertex order associated to variables x0, ..., xn (optional).
                - homog         : asks whether to construct the homogeneous ideal or not (optional, default is False).

            The function returns a list with the quotient ring and a tuple with the polynomial ring and the ideal.
        """
        if not homog:                                           ## Not homogeneous
            # Define the multivariate polynomial ring
            R = PolynomialRing(coeff_ring, ['x%s'%v for v in range(len(self.vertices))])    #type: ignore
            t = R.gens()
            if order == []:
                x = {self.vertices[i]:t[i] for i in range(len(self.vertices))}      # Dictionary with variables
            else:
                x = {order[i]:t[i] for i in range(len(self.vertices))}

            # Create a list with generators for the ideal
            ideal_gens = []
            for v in self.vertices:         # Toppling polynomials
                poly = - x[v]**(self.sandpile_struct.out_degree(v)) + R.prod([x[e[1]] for e in self.sandpile_struct.edges(v, labels = False) if e[1] in self.vertices])
                ideal_gens.append(poly)
            I = R.ideal(ideal_gens)         # Define the ideal

            return (R, I)
        else:                                                   ## Homogeneous
            # Define the multivariate polynomial ring
            R = PolynomialRing(coeff_ring, ['x%s'%v for v in range(len(self.vertices)+1)])    #type: ignore
            t = R.gens()
            if order == []:
                order = self.vertices + [self.sink()]
            x = {order[i]:t[i] for i in range(len(self.vertices)+1)}

            # Create a list with generators for the ideal
            ideal_gens = [x[self.sink()]-1]
            for v in self.vertices:         # Toppling polynomials
                poly = - x[v]**(self.sandpile_struct.out_degree(v)) + R.prod([x[e[1]] for e in self.sandpile_struct.edges(v, labels = False)])
                ideal_gens.append(poly)
            I = R.ideal(ideal_gens)         # Define the ideal

            return (R, I)


    def symm_poly(self, coeff_ring, config, order = []):                    ## Given a configuration, this function returns the symmetrized
        r"""
            This function computes the polynomial symmetric with respect to perm_group associated to the configuration config.
        """
        (R,I) = self.associated_ring(coeff_ring, order = order)    # Define the ring
        t = R.gens()
        if order == []:
            order = self.vertices
        x = {order[i]:t[i] for i in range(len(self.vertices))}

        poly = R.prod([x[v]**config[v] for v in self.vertices])
        Sn = SymmetricGroup(range(len(order)))      #type:ignore
        for subgr in self.perm_group:
            new_poly = 0
            indexes = [i for i in range(len(order)) if order[i] in subgr]
            for sigma in Permutations(indexes):     #type: ignore
                ext_dict = {v:v for v in order if v not in subgr} | {subgr[i]:order[sigma[i]] for i in range(len(indexes))}
                ext_perm = [ext_dict[v] for v in order]
                ext_sigma = Sn([ext_perm.index(v) for v in order])   #type: ignore
                new_poly = new_poly + poly(*ext_sigma(t))
            poly = new_poly
        return poly
    

    def sortrec_ideal(self, coeff_ring, sorted = True, opt = 2, override = False, order = [], homog = False):     ## Compute the sorted recurrent ideal
        r"""
            Compute the ideal to sorted recurrents. The possible arguments are:
                - coeff_ring    : specify the coefficient ring for the polynomial ring.
                - sorted        : boolean, it indicates if the considered sandpile is sorted (optional, default is True).
                - opt           : option for computing sorted recurrents (optional).
                - override      : boolean, if true the sorted recurrents are computed even if sorted_rec is non-empty (optional, default is False)
                - order         : the vertex order associated to variables x0, ..., xn (optional).
                - homog         : asks whether to construct the homogeneous ideal or not (optional, default is False).

            The function returns a list with the quotient ring, the polynomial ring and the ideal.
        """
        if not homog:
            if self.sorted_rec == [] or override:                       # Computes the sorted recurrents if necessary
                self.sorted_recurrents(option = opt)

            (R,I) = self.associated_ring(coeff_ring, sorted=sorted, order=order)     # Construct the associated ring
            t = R.gens()
            if order == []:
                order = self.vertices
            x = {order[i]:t[i] for i in range(len(self.vertices))}
            
            ideal_gens = []                                             # Compute ideal generators for R
            for conf in self.sorted_rec:
                if sorted:
                    poly = self.symm_poly(coeff_ring, conf, order=order)
                else:
                    poly = R.prod([x[v]**conf[v] for v in self.vertices])
                ideal_gens.append(poly)
            J = R.ideal(ideal_gens)

            return (R,I,J)
        else:
            if self.sorted_rec == [] or override:                       # Computes the sorted recurrents if necessary
                self.sorted_recurrents(option = opt)

            (R,I) = self.associated_ring(coeff_ring, sorted=sorted, order=order, homog=True)     # Construct the associated ring
            t = R.gens()
            if order == []:
                order = self.vertices + [self.sink()]
            x = {order[i]:t[i] for i in range(len(self.vertices)+1)}

            ideal_gens = []                                             # Compute ideal generators for R
            for conf in self.sorted_rec:
                if sorted:
                    poly = self.symm_poly(coeff_ring, conf, order=[v for v in order if v != self.sink()])
                else:
                    poly = R.prod([x[v]**conf[v] for v in self.vertices])
                ideal_gens.append(poly)
            J = R.ideal(ideal_gens)

            return (R,I,J)
    
    
    def show(self, default = False):                                        ## Function that displays the Sorted Sandpile
        r"""
            This function plots the sorted sandpile. The possible arguments are:
            - basic = False : if True, the representation is the same as Sandpile.show().

            If specific_opt is non-empty, the show function uses specific parameters:
            - default                   : call to the Sandpile.show() function
            - clique-independent graph  : call to the Sandpile.show() function fixing the position of the vertices,
                                          the sink in the center and the nonsink in a circle.
        """
        if self.specific_opt == [] or default:                         # Default
            if self.sandpile_struct.has_multiple_edges():
                self.sandpile_struct.show(edge_labels=True)
            else:
                self.sandpile_struct.show()
        elif self.specific_opt[0] == "clique-indep":        # Clique-independent graph
            [mu, nu] = self.specific_opt[1]
            mu_num = sum(mu)
            nu_num = sum(nu)
            # Define the position of each vertex
            positions = {0:(0,0)} | {i+1:(-np.sin(2*np.pi*i/(mu_num + nu_num)), np.cos(2*np.pi*i/(mu_num + nu_num)))  for i in range(mu_num + nu_num)}
            # Define the color of each part
            palette = rainbow(len(mu) + len(nu) + 1) # type: ignore
            col = {palette[0]:[0]} | {palette[i+1]:self.perm_group[i] for i in range(len(nu) + len(mu))}

            G = Graph(self.sandpile_struct.dict()).to_undirected() # type: ignore
            
            G.show(pos = positions, vertex_colors = col)
        elif self.specific_opt[0] == "gen-clique-indep":        # Generalized clique-independent graph
            vertex_set = self.specific_opt[1].vertices()
            cliq = []
            indep = []
            for vert in vertex_set:
                if self.specific_opt[2][vert] > 0:
                    cliq.append(vert)
                else:
                    indep.append(vert)
            palette = rainbow(2)    #type: ignore
            col = {palette[0]:cliq, palette[1]:indep}    #type: ignore
            self.specific_opt[1].show(vertex_colors = col, vertex_labels=self.specific_opt[2])
        elif self.specific_opt[0] == "mul-clique-indep":            # Multi clique-independent graph
            [mu, nu] = self.specific_opt[1]
            mu_num = sum(mu)
            nu_num = sum(nu)
            # Define the position of each vertex
            positions = {0:(0,0)} | {i+1:(-np.sin(2*np.pi*i/(mu_num + nu_num)), np.cos(2*np.pi*i/(mu_num + nu_num)))  for i in range(mu_num + nu_num)}
            # Define the color of each part
            palette = rainbow(len(mu) + len(nu) + 1) # type: ignore
            col = {palette[0]:[0]} | {palette[i+1]:self.perm_group[i] for i in range(len(nu) + len(mu))}

            G = Graph(self.sandpile_struct.dict()).to_undirected() # type: ignore
            
            G.show(pos = positions, vertex_colors = col, edge_labels = True)
        elif self.specific_opt[0] == "mulgen-clique-indep":    # Multi generalized clique-independent graph
            vertex_set = self.specific_opt[1].vertices()
            cliq = []
            indep = []
            for vert in vertex_set:
                if self.specific_opt[2][vert] > 0:
                    cliq.append(vert)
                else:
                    indep.append(vert)
            palette = rainbow(2)    #type: ignore
            col = {palette[0]:cliq, palette[1]:indep}    #type: ignore
            self.specific_opt[1].show(vertex_colors = col, vertex_labels=self.specific_opt[2], edge_labels=False)
    

    def export(self, saveopt = 0, opt = 2):                 ## Export informations for the Sorted Sandpile
        r"""
            This function returns the critical information of the sorted sandpile in a format that can be saved using pickle.
            If saveopt = 0, the function returns a list with:
                -   saveopt
                -   a dictionary of the underlying graph
                -   the sink of the sandpile
                -   the permutation group
                -   the specific options for the sandpile
                -   the order of vertices for...
                -   ...the list of sorted recurrents
                -   the qt-polynomial associated to the sorted sandpile.
        """
        qt_poly = self.qt_Polynomial(opt = opt)
        key_list = list(self.sorted_rec[0].keys())
        conflist = [[conf[i] for i in key_list] for conf in self.sorted_rec]
        return [saveopt, self.sandpile_struct.dict(), self.sandpile_struct.sink(), self.perm_group, self.specific_opt, key_list, conflist, qt_poly]


    def save(self, namefile, saveopt = 0, opt = 2):         ## Save the information on a file
        r"""
            This function saves the critical information of the sorted sandpile in a namefile
        """
        info = self.export(saveopt = saveopt, opt = opt)
        
        with open(namefile, 'wb') as handle:
            pickle.dump(info, handle)
    

    def load(namefile):                                     ## Load the information for a Sandpile
        r"""
            This function reads a file and returns a sorted sandpile with the information.
        """
        with open(namefile, 'rb') as handle:
            info = pickle.load(handle)

        if info[0] == 0:
            sorted_rec = []
            keys = info[5]
            for conf in info[6]:
                sort_conf = {keys[i]:conf[i] for i in range(len(keys))}
                sorted_rec.append(sort_conf)
            S = SortedSandpile(info[1], info[2], info[3], opt=info[4], sort_rec = sorted_rec)
            return S
        else:
            raise ImportError("The file cannot be read.")


                    ########################################################
                    #                   Specific Sandpiles                 #
                    ########################################################


def CliqueIndependent_SortedSandpile(mu, nu):                                                                                           ## Specific type of Sandpile
    r"""
        Construction of a Sorted Sandpile on the clique-independent graph given by parameters:
            - mu    : partition associated to the number and size of cliques in the graph.
            - nu    : partition associated to the number and size of independents in the graph.
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


def General_CliqueIndependent_SortedSandpile(cells_graph, card_cell, order_cells = []):                                                 ## Specific type of Sandpile
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
        dict_cell = dict_cell | {cell:[index_next_cell + j for j in range(abs(card_cell[cell]))]}
        index_next_cell += abs(card_cell[cell])

    order = []                                                  # Reading order for qt-Polynomials
    for cell in order_cells:                                    
        order = order + dict_cell[cell]
    
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
            edges_set = edges_set | {now:next_neigh}
            now += 1

        index_next_cell += abs(card_cell[cell])

    G = Graph(edges_set)    #type: ignore
    spec_opt = ["gen-clique-indep", cells_graph, card_cell, order]
    '''
        Specific Options:
        1-  Cell graph
        2-  Cell cardinality dictionary
        3-  Reading order (TODO)
    '''
    S = SortedSandpile(G, 0, perm_group, spec_opt)
    return S


def Multi_CliqueIndependent_SortedSandpile(mu, nu, kmul, hmul = -1, sinkmul = 1):                                                                    ## Specific type of Sandpile
    r"""
        Construction of the Sorted Sandpile, given two partitions mu, nu, where edges have multeplicity k in each clique and multeplicity h between components.
    """
    if hmul == -1:     # If third argument is not given, assume it is equal to k
        hmul = kmul

    mu_num = sum(mu)                                    # Number of independent vertices
    nu_num = sum(nu)                                    # Number of vertices in cliques
    d = {0 : [i+1 for i in range(mu_num + nu_num)]*sinkmul}     # Initialize the dictionary that will define the graph, link the sink to each vertex
    perm_group = []                                     # Initialize the permutation group acting on the graph

    part_first = 1                                      # Keeps track of first vertex of current part
    for part_nu in nu:                                          # Add edges for independent vertices. For each part...
        for i in range(part_nu):
            indep = [vert for vert in range(mu_num + nu_num + 1) for mult in range(kmul-1) if (part_first <= vert) and (vert < part_first + part_nu) and (vert != part_first + i)]
            others = [vert for vert in range(mu_num + nu_num + 1) for mult in range(hmul) if ((0 < vert) and (vert < part_first)) or (vert >= part_first + part_nu)]
            d[part_first + i] = [0]*sinkmul + indep + others
                                                                        # ...add all edges except for other vertices in part_nu    
        perm_group.append([part_first+j for j in range(part_nu)])       # Add the permutation orbit for the nu_part
        part_first += part_nu
    mu_rev = copy.copy(mu)                                              # We need the reversed partition mu...
    mu_rev.reverse()
    for part_mu in mu_rev:                                      # Add edges for clique sets. For each part...
        for i in range(part_mu):
            clique = [vert for vert in range(mu_num + nu_num + 1) for mult in range(kmul) if (part_first <= vert) and (vert < part_first + part_mu) and (vert != part_first + i)]
            others = [vert for vert in range(mu_num + nu_num + 1) for mult in range(hmul) if ((0 < vert) and (vert < part_first)) or (vert >= part_first + part_mu)]
            d[part_first + i] = [0]*sinkmul + clique + others
        perm_group.append([part_first+j for j in range(part_mu)])       # Add the permutation orbit for the mu_part
        part_first += part_mu

    ordered = []
    for i in range(len(mu)):
        temp = copy.copy(perm_group[len(perm_group)-i-1])
        temp.sort()
        ordered = ordered + temp
    for j in range(len(nu)):
        temp = copy.copy(perm_group[len(nu)-j-1])
        temp.sort(reverse = True)
        ordered = ordered + temp
    
    G = Graph(d)            # type: ignore      # Define the clique-independent graph
    specif_opt = ["mul-clique-indep", [mu, nu],  ordered, kmul, hmul]
    '''
        Specific Options:
        1-  list of the defining partitions
        2-  reading order
    '''
    S = SortedSandpile(G, 0, perm_group, specif_opt)
    return S


def MultiGeneral_CliqueIndependent_SortedSandpile(cells_graph, card_cell, multi_sink = 1, multiedge_cell = {}, order_cells = []):       ## Specific type of Sandpile
    r"""
        Construction of a Sorted Sandpile on the generalized clique-independent graph. Arguments are a graph and a dictionary with values for vertices, such that:
            - nodes:    represents a clique or an independent component with |card_cell[node]| vertices.
                        If the sign of card_cell[node] is:
                            - positive:     we have a clique.
                            - negative:     we have an independent.
            - edges:    each edge correspond to all possible edges between the two cells connected.
        If multiedge_cell is non-empty, each clique will have edges of multeplicity multiedge_cell[node].
        If cells_graph has multiple edges, the edges between components will be multiple.
    """
    cell_list = cells_graph.vertices()
    
    if order_cells == []:
        order_cells = copy.copy(cell_list)
    
    num_vert = sum([abs(card_cell[i]) for i in cell_list])      # Set the total number of vertices for sandpile

    dict_cell = {}                                              # Dictionary "cell -> list vertices"
    index_next_cell = 1
    for cell in cell_list:
        dict_cell = dict_cell | {cell:[index_next_cell + j for j in range(abs(card_cell[cell]))]}
        index_next_cell += abs(card_cell[cell])
    
    order = []                                                  # Reading order for qt-Polynomials
    for cell in order_cells:                                    
        order = order + dict_cell[cell]

    edges_set = {0:[i+1 for i in range(num_vert) for k in range(multi_sink)]}   # Set of edges, sink to everyone with multeplicity "multi_sink"
    perm_group = []                                                             # Set the permutation group

    now = 1
    index_next_cell = 1
    for cell in cell_list:                              # Add the cell
        perm_group.append(dict_cell[cell])

        for vert in range(abs(card_cell[cell])):
            # Add cell with its own vertices
            if card_cell[cell] > 0:                         # Add for clique
                if cell in multiedge_cell.keys():               # If we have multiple edges...
                    next_neigh = [0 for k in range(multi_sink)] + [index_next_cell + j for j in range(abs(card_cell[cell])) for k in range(multiedge_cell[cell]) if index_next_cell + j != now]
                else:                                           # If not...
                    next_neigh = [0 for k in range(multi_sink)] + [index_next_cell + j for j in range(abs(card_cell[cell])) if index_next_cell + j != now]
            else:                                           # Add nothing for independent
                next_neigh = [0 for k in range(multi_sink)]
            # Add complete with outer edges
            for neigh_cell in cells_graph.neighbors(cell):
                multi = list(cells_graph.edges(labels=False)).count((cell, neigh_cell))
                next_neigh = next_neigh + [x for x in dict_cell[neigh_cell] for k in range(multi)]

            # Add to the dictionary
            edges_set = edges_set | {now:next_neigh}
            now += 1

        index_next_cell += abs(card_cell[cell])
    
    # Conversion to sage-math's format
    multi_edges_set = {}
    for v in range(num_vert + 1):
        vertex_dict = {}
        for w in set(edges_set[v]):
            vertex_dict = vertex_dict | {w:list(edges_set[v]).count(w)}
        multi_edges_set = multi_edges_set | {v:vertex_dict}
    # Construct the graph
    G = Graph(edges_set)    #type: ignore
    spec_opt = ["mulgen-clique-indep", cells_graph, card_cell, order]
    '''
        Specific Options:
        1-  Cell graph
        2-  Cell cardinality dictionary
        3-  Reading order
    '''
    S = SortedSandpile(G, 0, perm_group, spec_opt)
    return S