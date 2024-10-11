#                           GENERALIZED CLIQUE-INDEPENDENT SORTED SANDPILES
#
#   These are some tests for particular structures of generalized clique-independent sorted sandpiles

from sage.combinat.q_analogues import qt_catalan_number
import time

exp_num = 0

if exp_num == 0:        #   Try a chain of independent sets: two parameters (a,n):
                        #       a   : vertices in one independent
                        #       n   : length of a chain
    a = 2
    poly_list = []
    for a in range(4):
        a = -a-1
        for n in range(1):
            n = 2
            print("Computing case a={} and n={}".format(a,n))
            t = time.time()
            ind_dict = {i:[i+1] for i in range(n)}  # Structure of gen independent graph
            conf = {i:a for i in range(n+1)}        # All independents have "a" vertices
            G = Graph(ind_dict)
            S = General_CliqueIndependent_SortedSandpile(G, conf)
            #S.sandpile_struct.show()
            q_poly = S.q_Polynomial()
            poly_list.append(q_poly)
            print("Finished in {} seconds".format(time.time() - t))
    eval_list = [p(1) for p in poly_list]
    print(poly_list)
    print(eval_list)
    print(oeis(eval_list))
