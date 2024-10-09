import itertools
import copy
import numpy as np

                    #####################################################
                    #           Sorted Sandpile Configuration           #
                    #####################################################

class SandpileSortConfig():

    def __init__(self, sandp, config, permut):      ## Defines the class SandpileSortConfig
        r"""
            Definition of the class SandpileSortConfig.
            
            The arguments given are a sandpile and a partition of the vertex set (without sink).
        """

        ### TODO: check if the arguments given correspond to a good SandpileSortConfig

        self.sandpile_struct = sandp                                                        # Define the sandpile_struct as sandp
        self.sandpile_config = SandpileConfig(self.sandpile_struct, config) # type: ignore  # Define the configuration for 
        self.perm_group = permut                                                            # Define the permutation group on the sandpile


    def __repr__(self):                             ## Returns a description of the class Sandpile2
        return "A sorted configuration for a sandpile with vertices {} and sink {}".format(self.sandpile_struct.vertices(), self.sandpile_struct.sink())
    

    def sort(self, change = True):
        r"""
            Rearrange the configuration values by increasing order on each orbit.

            If change is True, it modifies the configuration, if it is False, it just returns the sorted one.
        """
        result = self.sandpile_config
        for part in self.perm_group:
            part_conf = [result[i] for i in part]     # Get the partial configuration associated with an orbit
            part_conf.sort()                          # Sort it
            for j in range(len(part)):                # Modify the configuration
                result[part[j]] = part_conf[j]
        if change:                                    # Changing the actual configuration
            self.sandpile_config = result
        return result
    

    def __eq__(self, other):                        ## Compares two sorted configurations to see if they are equivalent under the permutation group
        r"""
            Checks if the two sorded configurations are the same or not.
        """
        if (self.sandpile_struct != other.sandpile_struct) or (self.perm_group != other.perm_group):    # Not comparable
            raise Exception("The two sorted configurations are not comparable.")
        elif self.sort(change = False) == other.sort(change = False):                                   # Compare sorted configurations
            return True
        else:
            return False
            
    
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

    def topple_sink(self):                          ## Topples the sink vertex
        r"""
            Topple the sink and change the configuration stored.
        """
        for vert in self.sandpile_struct.neighbors(self.sandpile_struct.sink()):
            self.sandpile_config[vert] += 1
        return self.sandpile_config
    

    def topple_vertex(self, vert):                  ## Topples a vertex of the configuration
        r"""
            Topple a given vertex of the sandpile. This changes the configuration stored.
            
            This function just calls the function fire_vertex() from the class SandpileConfig and applies it to self.sandpile_config. 
        """
        self.sandpile_config = self.sandpile_config.fire_vertex(vert)
        return self.sandpile_config
    

    def topple_sequence(self, list):                ## Topples a sequence of vertices
        r"""
            Topple a list of vertices of the sandpile. This changes the configuration stored.
            
            This function calls iteratively the function fire_vertex() from the class SandpileConfig and applies it to self.sandpile_config. 
        """
        for vert in list:
            self.sandpile_config = self.sandpile_config.fire_vertex(vert)
        return self.sandpile_config
    

    def level(self):                                ## Returns the level statistic of the configuration
        r"""
            Returns the level of the configuration.
        """
        not_inc_sink = len(self.sandpile_struct.to_undirected().edges()) - self.sandpile_struct.out_degree(self.sandpile_struct.sink())
        return (- not_inc_sink + self.sandpile_config.deg())
    

    def delay(self, order = [], check_rec = True):  ## Returns the delay statistic of the configuration
        r"""
            Given a reading order of nonsink vertices, this function computes the configuration's delay.

            - order         : the order for reading vertices. If undefined, the decreasing order on vertices is assumed.
            - check_rec     : this option can be used to override the is_recurrent() call.
        """
        if not self.is_recurrent() and check_rec:           # Check if the configuration is recurrent
            raise Exception("The sorted configuration is not recurrent, hence delay is not defined.")

        nonsink_vert = list(self.sandpile_struct.nonsink_vertices())

        if order == []:                                     # No order has been assigned: take decreasing order
            order = copy.copy(nonsink_vert)
            order.sort()
            order.reverse()
        
        delay = 0                                                   # The delay value
        loop_count = 0                                              # The loop count
        not_toppled = copy.copy(nonsink_vert)                                  # A list with the vertices still to be toppled
        self.topple_sink()                                  # Start by toppling the sink
        while len(not_toppled) > 0:
            for ind in order:                               # Search for non-toppled vertices in the right order
                value = self.sandpile_config[ind]
                if (self.sandpile_struct.out_degree(ind) <= value) and (ind in not_toppled):
                                                            # Can be toppled!
                    delay += loop_count                             # Raise the delay by suitable amount
                    self.topple_vertex(ind)                         # Topple the vertex
                    not_toppled.remove(ind)                         # Remove the vertex from not-toppled
            loop_count += 1
        return delay
    



                    #####################################################
                    #                   Sorted Sandpile                 #
                    #####################################################

class SortedSandpile():

    def __init__(self, graph, sink, permut, opt=[]):            ## Initialize the class SortedSandpile
        r"""
            Definition of the class SortedSandpile.
            
            The arguments given are a sandpile and a partition of the vertex set (without sink).
        """

        ### TODO: check if the arguments given correspond to a good SandpileSortConfig
        
        self.sandpile_struct = Sandpile(graph, sink) # type: ignore                         # Define the Sandpile structure
        self.perm_group = permut                                                            # Define the permutation group on the graph vertices
        self.specific_opt = opt                                                             # Define the specific options, to optimize in specific cases
        self.sorted_rec = []                                                                # Eventually store the sorted recurrent configurations if computed


    def __repr__(self):                             ## Description of SortedSandpile class
        return "A sorted sandpile on vertices {} and sink {}.".format(self.sandpile_struct.vertices(), self.sandpile_struct.sink())
    

    def simple_recurrents(self):                    ## Computes the simple recurrents
        r"""
            Returns the list of recurrent configurations ignoring the action of perm_group.

            It's a call of the recurrents() function from Sandpile class.
        """
        return self.sandpile_struct.recurrents()
    

    def sorted_recurrents(self, option = 0):        ## Computes a list of all sorted recurrent configurations
        r"""
            Computes the sorted recurrent configurations.

            - option        : based on the value the search algorithm changes
                = 0     : calls the recurrent() function from Sandpile class [default]
                = 1     : HOPEFULLY ONE DAY...
        """
        match option:
            
            case 0:                                                         ## OPTION 0
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
                    self.sorted_rec = self.sorted_rec + [{self.sandpile_struct.nonsink_vertices()[j]:simpl_rec[i][j] for j in range(len(simpl_rec[i]))}]
                return self.sorted_rec
            
            case _:
                raise Exception("The given option is not valid.")
            

    def specific_sort_recurrents(self):
        if self.sorted_rec != []:
            if self.specific_opt != []:         # Order the configurations in a specific order
                if self.specific_opt[0] == "clique-indep":          # Clique-independent graphs
                    for conf in range(len(self.sorted_rec)):
                        new_config = {}
                        [mu, nu] = self.specific_opt[1]
                        for j in range(len(nu)):
                            temp = copy.copy([self.sorted_rec[conf][i] for i in self.perm_group[j]])
                            temp.sort()
                            new_config = new_config | {self.perm_group[j][b]:temp[b] for b in range(len(self.perm_group[j]))}
                        for j in range(len(mu)):
                            temp = copy.copy([self.sorted_rec[conf][i] for i in self.perm_group[len(nu)+j]])
                            temp.sort()
                            temp.reverse()
                            new_config = new_config | {self.perm_group[len(nu) + j][b]:temp[b] for b in range(len(self.perm_group[len(nu) + j]))}
                        self.sorted_rec[conf] = new_config
            else:
                raise Exception("Sorted Sandpile has no specific option!")
        else:
            raise Exception("Sorted recurrent configurations not yet computed!")

    
    def qt_Polynomial(self, ordered = []):            ## Computes the q,t polynomial on (level, delay)
        r"""
            Returns the q,t - polynomial corresponding to the sorted sandpile's recurrent configurations.

            - order     : if specified, it fixes the reading order for the delay statistic.
        """
        R = FractionField(QQ['q, t'])      # type: ignore
        q,t = R.gens()
        poly = 0*q*t                        # Define the polynomial as 0

        if self.sorted_rec == []:           # If sorted recurrents still to be computed...
            self.sorted_recurrents()        # ...compute them!

        # TODO: be sure that the delay doesn't depend on the order in the same orbit...
        if self.specific_opt != []:         # If there is a specific option...
            self.specific_sort_recurrents() # ...sort in a particular way each configuration

        for config in self.sorted_rec:      # Compute the polynomial
            sortedconfig = SandpileSortConfig(self.sandpile_struct, config, self.perm_group)
            q_exp = sortedconfig.level()
            t_exp = sortedconfig.delay(order = ordered, check_rec=False)
            poly = poly + (q**q_exp) * (t**t_exp)
        return poly
    
    
    def show(self):                                 ## Function that displays the Sorted Sandpile
        r"""
            This function plots the sorted sandpile.

            If specific_opt is non-empty, the show function uses specific parameters:
            - default                   : call to the Sandpile.show() function
            - clique-independent graph  : call to the Sandpile.show() function fixing the position of the vertices,
                                          the sink in the center and the nonsink in a circle.
        """
        if self.specific_opt == []:                     # Default
            self.sandpile_struct.show()
        elif self.specific_opt[0] == "clique-indep":    # Clique-independent graph
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
            
    

def CliqueIndependent_SortedSandpile(mu, nu):       ## Specific type of Sandpile
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
    
    G = Graph(d)            # type: ignore      # Define the clique-independent graph
    S = SortedSandpile(G, 0, perm_group, ["clique-indep", [mu, nu]])
    return S
