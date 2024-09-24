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
        self.v_name[self.sink] = 0
        del temp_vertex, temp_key
        
        self.current_conf = [0 for i in range(self.v_len)]                                      # Sets the current configuration
        self.directed_neigh = [self.v_index[j] for i in range(self.v_vertices) for j in range(self.graph.neighbors(self.v_name[j]))]             # Records the outward edges of each vertex, using the new naming scheme
        self.degree = [len((self.graph).neighbors(self.v_names[i])) for i in self.v_vertices]   # Computes the degree of each vertex
        
         
    
    def __repr__(self):                     ## Returns a description of the class Sandpile2
        return "A sandpile structure. The graph has vertex set {} and sink {}".format(self.vertices, self.sink) 
        
    
    def show(self):                         ## Returns an image of the graph with labeling given by the configuration
        r"""
            Draws the usual representation of the graph.
        """
        self.graph.show()
        
    
