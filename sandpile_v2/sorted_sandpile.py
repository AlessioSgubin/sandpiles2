import copy

class SandpileSortConfig():

    def __init__(self, sandp, permut):
        r"""
            Definition of the class SandpileSortConfig. The arguments given are a sandpile and a partition of the vertex set.
        """

        

        self.sandpile_struct = sandp        # Define the sandpile_struct as sandp
        self.perm_group = permut            # Define the permutation group on the sandpile

