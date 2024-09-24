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
        self.v_names = {i:self.vertices[i] for i in self.v_vertices}                    # Names for virtual vertices
        # Swap vertex associated to 0 and the sink, so the sink is associated to 0
        temp_vertex = self.v_names[0]
        temp_key = (list(self.v_names.keys())[list(self.v_names.values()).index(self.sink)])
        self.v_names[temp_key] = temp_vertex
        self.v_names[0] = self.sink
        del temp_vertex, temp_key
        
        self.current_conf = [0 for i in range(self.v_len)]                                      # Sets the current configuration
        self.degree = [len((self.graph).neighbors(self.v_names[i])) for i in self.v_vertices]   # Computes the degree of each vertex
        
         
    
    def __repr__(self):                     ## Returns a description of the class Sandpile2
        return "A sandpile structure. The graph has vertex set {} and sink {}".format(self.vertices, self.sink) 
        
    
    def show(self):                         ## Returns an image of the graph with labeling given by the configuration
        r"""
            Draws the usual representation of the graph.
        """
        self.graph.show()
        
    
