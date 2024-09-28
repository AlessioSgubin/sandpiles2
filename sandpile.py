class Sandpile2():          ### Define the class Sandpile2

    def __init__(self, G, sink, sparse=False):            ## Implicit definition of Sandpile2
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
        self.current_conf = [0 for i in range(self.v_len)]                                                                # Sets the current configuration
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


    def is_recurrent(self):
        r"""
            Checks if the current configuration is recurrent.
            The function returns the boolean value and a good toppling order.
        """
        start_conf = self.current_conf                                                      # The configuration before computations
        
        self.topple_vertex(self.sink)                                                       # Topple the sink first
        branchings = [ [self.current_conf, [0], [i+1 for i in range(self.v_len - 1)], [x for x in self.v_vertices if self.current_conf[x] >= self.out_degree[x] ]] ]   
                                                                                            # List of the branching points for the depth-first search of recurrent orders.
                                                                                            # It stores lists containing "[configuration, toppling order, remaining to be toppled, good candidates]"
        depth = 0                                                                           # Stores the depth of the current search
        
        while depth >= 0:                           # Continue the search depth-first until we checked everything
            if depth == self.v_len - 1:                     # All vertices have been toppled
                self.current_conf = start_conf
                return [True, (branchings[depth])[1]]               # Return True and the toppling order
            elif len(branchings[depth][3]) == 0:            # If there is no more candidate, go back
                if depth == 0:                                      # Check if we have finished candidates
                    self.current_conf = start_conf
                    return [False, 0]
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
                    return [True, (branchings[depth])[1]]               # Return True and the toppling order
                new_branch = [self.current_conf, new_toppl_ord, new_remaining, new_good_cand]  
                branchings.append(new_branch)                           # Append the new branching point
                depth += 1
        self.current_conf = start_conf
        return [False, [0]]

    