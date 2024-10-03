import itertools

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

        self.sandpile_struct = sandp                                            # Define the sandpile_struct as sandp
        self.sandpile_config = SandpileConfig(self.sandpile_struct, config)     # type: ignore # Define the configuration for 
        self.perm_group = permut                                                # Define the permutation group on the sandpile


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
        not_inc_sink = len(self.sandpile_struct.vertices()) - 1 - self.sandpile_struct.out_degree(self.sandpile_struct.sink())
        return - not_inc_sink + self.sandpile_config.deg()
    

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
            order = nonsink_vert
            order.sort()
            order.reverse()
            print(order)
        
        delay = 0                                                   # The delay value
        loop_count = 0                                              # The loop count
        not_toppled = nonsink_vert                                  # A list with the vertices still to be toppled
        self.topple_sink()                                  # Start by toppling the sink
        while len(not_toppled) > 0:
            for ind in order:                               # Search for non-toppled vertices in the right order
                if (self.sandpile_struct.out_degree(ind) <= self.sandpile_config[ind]) and (ind in not_toppled):
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

    def __init__(self, graph, sink, permut):              ## Initialize the class SortedSandpile
        r"""
            Definition of the class SortedSandpile.
            
            The arguments given are a sandpile and a partition of the vertex set (without sink).
        """

        ### TODO: check if the arguments given correspond to a good SandpileSortConfig
        
        self.sandpile_struct = Sandpile(graph, sink) # type: ignore
        self.perm_group = permut

    def __repr__(self):                             ## Description of SortedSandpile class
        return "A sorted sandpile on vertices {} and sink {}.".format(self.sandpile_struct.vertices(), self.sandpile_struct.sink())
    

    def simple_recurrents(self):                    ## Computes the simple recurrents
        r"""
            Returns the list of recurrent configurations ignoring the action of perm_group.

            It's a call of the recurrents() function from Sandpile class.
        """
        return self.sandpile_struct.recurrents()
    

    def sorted_recurrents(self, option = 0):
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
                    simpl_rec[i] = temp.sort()
                simpl_rec.sort()                            # Remove duplicates
                return list(simpl_rec for simpl_rec,_ in itertools.groupby(simpl_rec))
            
            case _:
                raise Exception("The given option is not valid.")