import copy                 # To deep-copy iterated lists
import itertools            # To remove unhashable duplicates in lists
import numpy as np          # To use arrays

class Sandpile2():          ### Define the class Sandpile2

    def __init__(self, G, sink):            ## Implicit definition of Sandpile2
        r"""
            Creates an object SandPile2, the defining arguments are the graph (directed or undirected) and a sink.
        """
    ### Recording initial data and defining values
        self.graph = G                              # Graph structure
        self.sink = sink                            # Sink
        self.vertices = G.vertices(sort = True)     # Vertex set (sorted, so we have an order on the edges)
        
        ### Checking whether the sink is a vertex, if not declare 0 to be the new sink (adding it if necessary)
        if not (self.sink in self.vertices):
            self.sink = 0
            self.graph.add_vertex(0)
            self.vertices = self.graph.vertices(sort = True)
        
    ### Initializing other data
        self.v_len = len(self.vertices)                                                 # Vertex number
        
        # Since the vertices can be given any name, we want to standardize their numbering for ease of computation
        self.v_vertices = [i for i in range(self.v_len)]                                # Virtual vertices
        self.v_name = {i:self.vertices[i] for i in self.v_vertices}                     # Dictionary that returns the actual vertex of a given index (0 -> sink)
        self.v_index = {self.vertices[i]:i for i in self.v_vertices}                    # Dictionary that returns the index of an actual vertex (sink -> 0)
        # Swap vertex associated to 0 and the sink, so the sink is associated to 0...
        temp_vertex = self.v_index[self.sink]       # ... for v_index
        temp_key = (list(self.v_index.keys())[list(self.v_index.values()).index(0)])
        self.v_index[temp_key] = temp_vertex
        self.v_index[self.sink] = 0
        temp_vertex = self.v_name[0]               # ... for v_name
        temp_key = (list(self.v_name.keys())[list(self.v_name.values()).index(self.sink)])
        self.v_name[temp_key] = temp_vertex
        self.v_name[0] = self.sink
        del temp_vertex, temp_key
        
        # Recovers informations on graph, translating to new numbering scheme
        self.current_conf = [0 for i in range(self.v_len)]                                                               # Sets the current configuration
        self.directed_neigh = [[self.v_index[j] for j in self.graph.neighbors(self.v_name[i])] for i in self.v_vertices]  # Records the outward edges of each vertex, using the new naming scheme
        self.out_degree = [len(self.directed_neigh[i]) for i in self.v_vertices]                                              # Computes the degree of each vertex
        
    
    def __repr__(self):                     ## Returns a description of the class Sandpile2
        return "A sandpile structure. The graph has vertex set {} and sink {}".format(self.vertices, self.sink) 
        

    def show(self, v_indexing=False):                         ## Returns an image of the graph with labeling given by the configuration
        r"""
            Draws the usual representation of the graph. The sink is colored in blue and the others in red.
            The possible options are:
                v_indexing          If set "True" the vertex labels are the virtual ones, where "0" is the sink
        """
        if v_indexing == False:
            self.graph.show(vertex_color = 'lightcoral', vertex_colors = {'lightskyblue':[self.sink]})

        else:
            self.graph.show(vertex_labels = self.v_index, vertex_color = 'lightcoral', vertex_colors = {'lightskyblue':[self.sink]})
        
    
    def topple_vertex(self, vert, virtual = True):
        r"""
            The function does a toppling on a given vertex. Returns the new current configuration.
        """
        if not virtual:         # If the vertex is not virtual, convert it
            vert = self.v_index[vert]
        self.current_conf[vert] -= self.out_degree[vert]
        for i in self.directed_neigh[vert]:
            self.current_conf[i] += 1
        return self.current_conf
    

    def topple_sequence(self, list_vert, virtual = True):
        r"""
            The function does a sequence of topplings on a list of vertices. Returns the new current configuration.
        """
        if not virtual:         # If the vertices are not virtual, convert them
            v_list = [self.v_index[i] for i in list_vert]
            list_vert = v_list
        for vert in list_vert:
            for i in self.directed_neigh[self.v_index[vert]]:
                self.current_conf[i] += 1
        return self.current_conf
    

    def is_stable(self):
        r"""
            Checks if the current configuration is stable.
        """
        for i in self.v_vertices:
            if (i != 0) and (self.current_conf[i] >= self.out_degree[i] or self.current_conf[i]< 0):
                return False
        return True


    def is_recurrent(self, virtual=True):
        r"""
            Checks if the current configuration is recurrent.
            The function returns the boolean value and a good toppling order. If there is no good toppling order, return list "[sink]".
        """
        start_conf = self.current_conf                                                      # The configuration before computations
        
        self.topple_vertex(self.sink)                                                       # Topple the sink first
        branchings = [ [self.current_conf, [0], [i+1 for i in range(self.v_len - 1)], [x for x in self.v_vertices if (self.current_conf[x] >= self.out_degree[x] and x != 0) ] ] ]   
                                                                                            # List of the branching points for the depth-first search of recurrent orders.
                                                                                            # It stores lists containing "[configuration, toppling order, remaining to be toppled, good candidates]"
        depth = 0                                                                           # Stores the depth of the current search
        
        while depth >= 0:                           # Continue the search depth-first until we checked everything
            if depth == self.v_len - 1:                     # All vertices have been toppled
                self.current_conf = start_conf
                if not virtual:                             # Return True and the toppling order
                    return [True, [self.v_name[x] for x in branchings[depth][1]]]
                else:
                    return [True, branchings[depth][1]]               
            elif len(branchings[depth][3]) == 0:            # If there is no more candidate, go back
                if depth == 0:                                      # Check if we have finished candidates
                    self.current_conf = start_conf
                    if not virtual:                         # Return False and the sink
                        return [False, [self.sink]]
                    else:
                        return [False, [0]]
                else:
                    depth -= 1                                      # Decrease the depth
                    branchings[depth][3].pop(0)                     # Remove the node from good candidates (it was the first one of the list)
                    branchings.pop()                                # Remove last element from branchings
            else:                                           # There is a new candidate, grow the branching tree
                new_cand = branchings[depth][3][0]
                self.topple_vertex(new_cand)                            # Topple the first candidate of the list (the good candidates are unstable vertices, so no need to check non-negativity)
                new_toppl_ord = branchings[depth][1]
                new_toppl_ord.append(new_cand)                          # Add the new candidate to the possible toppling order
                new_remaining = branchings[depth][2]
                new_remaining.remove(new_cand)                          # Removing new candidate from remaining vertices
                new_good_cand = [x for x in new_remaining if self.current_conf[x] >= self.out_degree[x]]
                                                                        # List the good next candidates, unstable vertices jet to be toppled
                if new_remaining == new_good_cand:                      # If all remaining are unstable, we are finished
                    self.current_conf = start_conf
                    new_toppl_ord.extend(new_remaining)
                    if not virtual:                                             # Return True and the toppling order
                        return [True, [self.v_name[x] for x in new_toppl_ord]]
                    else:
                        return [True, new_toppl_ord] 
                new_branch = [self.current_conf, new_toppl_ord, new_remaining, new_good_cand]  
                branchings.append(new_branch)                           # Append the new branching point
                depth += 1
        self.current_conf = start_conf
        if not virtual:                                         # Return False and the sink
            return [False, [self.sink]]
        else:
            return [False, [0]]


    def minimal_conf(self, configs):
        r"""
            This function takes a list of configurations and returns the subset of minimal ones with respect to the componentwise integer order.
        """
        ### Delete bad configs
        print(len(configs))
        count = 0
        for i in configs:
            print(count)
            j = 0
            check = False
            while j <= self.v_len and check == False:
                if i[j] >= self.out_degree[j]:
                    configs.pop(count)
                    check = True
                    print("YEP.")
                    count -= 1
                j += 1
            count += 1


        ### Delete duplicates
        print(len(configs))
        configs.sort()
        list(configs for configs,_ in itertools.groupby(configs))
        print(len(configs))

        ### Delete non-minimal
        count = 0
        iterat = 0
        index = 0
        for x in configs:
            iterat += 1
            print("Iterazione {} e rimossi {}".format(iterat, count))
            for y in configs[index+1:]:
                removable = True
                i = 0
                while (removable == True) and (i < len(x)):
                    if x[i] > y[i]:
                        removable = False
                    else:
                        i += 1
                if removable:
                    configs.remove(y)
                    count += 1
            index +=1
        count                               # Display the redundants
        return configs


    def recurrent_conf(self, minimal=False):
        r"""
            Function that computes all recurrent configurations. The possible options are:
                - minimal       : if True it runs a O(n!) algorithm to find just the minimal recurrent configurations
        """
        print("Starting the search for recurrent configurations...")
        start_conf = self.current_conf
        self.current_conf = [0 for i in self.vertices]                                      # The configuration before computations
        recurrent_list = []                                                                 # Initialize list for recurrents

        self.topple_vertex(self.sink)                                                       # Topple the sink first
        branchings = [ [self.current_conf, [0], [i+1 for i in range(self.v_len - 1)], [x for x in self.v_vertices if self.current_conf[x]>0], [0 for i in self.vertices]] ]
                                                                                            # List of the branching points for the depth-first search of recurrent orders.
                                                                                            # It stores lists containing "[how many grains added, toppling order, remaining to be toppled, good candidates, configuration being built]"
        depth = 0                                                                           # Stores the depth of the current search
        try:
            while depth >= 0:                           # Continue the search depth-first until we checked everything
                print("Depth {}".format(depth)) #, branching {} con dati {}".format(depth, len(branchings), branchings[depth]))
                if len(branchings[depth][3]) == 0:          # There is no more good candidate...
                    if depth == self.v_len - 1:                     # ...because we have a good configuration! (USELESS CASE, I THINK)
                        ###print("Good configuration")
                        recurrent_list.append(branchings[depth][4])             # Append the new configuration found
                        depth -= 1                                              # Decrease the depth
                        branchings[depth][3].pop(0)                             # Remove the node from good candidates (it was the first one of the list)
                        branchings[depth][4]
                        branchings.pop()                                        # Remove last element from branchings
                        self.current_conf = copy.copy(branchings[depth][0])     # Restore the configuration of grains added
                    elif depth == 0:                                # ...because we have finished the entire search!
                        ###print("Finished all!")
                        self.current_conf = start_conf                          # Reset the configuration
                        recurrent_list = self.minimal_conf(recurrent_list)      # Remove non-minimal recurrent configurations found
                        return [len(recurrent_list), recurrent_list]            # Return number of configurations and list
                    else:                                           # ...because we are on a dead branch!
                        ###print("Dead branch")
                        depth -= 1                                              # Decrease the depth
                        branchings[depth][3].pop(0)                             # Remove the node from good candidates (it was the first one of the list)
                        branchings.pop()                                        # Remove last element from branchings
                        self.current_conf = copy.copy(branchings[depth][0])     # Restore the configuration of grains added
                else:                                       # There is a new candidate, grow the branching tree
                    new_cand = branchings[depth][3][0]                          # Topple the first candidate of the list
                    ###print("New candidate {}".format(new_cand))
                    self.topple_vertex(new_cand)          
                    new_grainsadd = copy.copy(self.current_conf)
                    new_toppl_ord = copy.copy(branchings[depth][1])             # Add the new candidate to the possible toppling order
                    new_toppl_ord.append(new_cand)                          
                    new_remaining = copy.copy(branchings[depth][2])             # Removing new candidate from remaining vertices
                    new_remaining.remove(new_cand)                          
                    new_builtconf = copy.copy(branchings[depth][4])             # Adding to the configuration being built the minimal value for the new candidate to be toppled
                    new_builtconf[new_cand] = self.out_degree[new_cand] - self.current_conf[new_cand]
                    new_good_cand = [x for x in new_remaining if self.current_conf[x] > 0]
                                                                                # List the good next candidates, unstable vertices jet to be toppled
                    
                    new_branch = [new_grainsadd, new_toppl_ord, new_remaining, new_good_cand, new_builtconf]  
                    branchings.append(new_branch)                               # Append the new branching point
                    depth += 1
        except:
            for i in branchings:
                print(i)
        
        self.current_conf = start_conf
        print("...search has concluded! Now reduce the configurations to minimal ones...")
        recurrent_list = self.minimal_conf(recurrent_list)  # Remove non-minimal recurrent configurations found
        print("...return the results, we found {} minimal configurations".format(len(recurrent_list)))
        return [ len(recurrent_list), recurrent_list]       # Return the number and list of recurrent configurations
        
    
    def over_configs(self, min_conf):
        r"""
            This function returns a list of configurations starting from minimal set.
        """
        all_conf = []
        for v in min_conf:
            if v not in all_conf:
                branchings = [ [v, 0] ]          # Start a recurrent addition of configurations
                depth = 0
                while depth >= 0:
                    return False
                    ########## TO BE FINISHED!