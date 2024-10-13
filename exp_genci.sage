#                           GENERALIZED CLIQUE-INDEPENDENT SORTED SANDPILES
#
#   These are some tests for particular structures of generalized clique-independent sorted sandpiles

from sage.combinat.q_analogues import qt_catalan_number
import time

exp_num = 2

if exp_num == 0:        #   Try a chain of independent/clique sets: two parameters (a,n):
                        #       a   : vertices in one independent/clique
                        #       n   : length of a chain
    a = 2
    poly_list = []
    for a in range(1):
        #a = -a-1
        a = -4
        for n in range(4):
            #n = 2
            print("Computing case a={} and n={}".format(a,n))
            t = time.time()
            ind_dict = {i:[i+1] for i in range(n)}  # Structure of gen independent/clique graph
            conf = {i:a for i in range(n+1)}        # All independents/cliques have "a" vertices
            G = Graph(ind_dict)
            S = General_CliqueIndependent_SortedSandpile(G, conf)
            #S.sandpile_struct.show()
            q_poly = S.q_Polynomial()
            poly_list.append(q_poly)
            print("Time elapsed: {} seconds".format(time.time() - t))
            print("q-polynomial: {}".format(q_poly))
            print("Number sorted recurrents: {}\n".format(q_poly(1)))
    eval_list = [p(1) for p in poly_list]
    print(poly_list)
    print(eval_list)
    print(oeis(eval_list))

if exp_num == 1:        #   Try a loop of independent/clique sets: two parameters (a,n):
                        #       a   : vertices in one independent/clique
                        #       n   : length of the loop
    a = 2
    poly_list = []
    for a in range(1):
        #a = -a-1
        a = 4
        for n in range(3):
            #n = 0
            print("Computing case a={} and n={}".format(a,n))
            t = time.time()
            if n == 0:
                ind_dict = {i:[i+1] for i in range(n)}
            else:
                ind_dict = {n:[0]} | {i:[i+1] for i in range(n)}  # Structure of gen independent/clique graph
            conf = {i:a for i in range(n+1)}                  # All independents/cliques have "a" vertices
            G = Graph(ind_dict)
            S = General_CliqueIndependent_SortedSandpile(G, conf)
            #S.sandpile_struct.show()
            q_poly = S.q_Polynomial()
            poly_list.append(q_poly)
            print("Time elapsed: {} seconds".format(time.time() - t))
            print("q-polynomial: {}".format(q_poly))
            print("Number sorted recurrents: {}\n".format(q_poly(1)))
    eval_list = [p(1) for p in poly_list]
    print(poly_list)
    print(eval_list)
    print(oeis(eval_list))

if exp_num == 2:        #   Try a chain of independent/clique sets: two parameters (a,n):
                        #       a   : vertices in one independent/clique
                        #       n   : length of a chain
    a = 2
    poly_list = []
    for (a,n) in [(-3,2),(-4,1),(-2,4),(3,2),(4,1),(2,4)]:
        print("Computing case a={} and n={}".format(a,n))
        t = time.time()
        ind_dict = {i:[i+1] for i in range(n)}  # Structure of gen independent/clique graph
        conf = {i:a for i in range(n+1)}        # All independents/cliques have "a" vertices
        G = Graph(ind_dict)
        S = General_CliqueIndependent_SortedSandpile(G, conf)
        #S.sandpile_struct.show()
        q_poly = S.q_Polynomial()
        poly_list.append(q_poly)
        print("Time elapsed: {} seconds".format(time.time() - t))
        print("q-polynomial: {}".format(q_poly))
        print("Number sorted recurrents: {}\n".format(q_poly(1)))
    eval_list = [p(1) for p in poly_list]
    print(poly_list)
    print(eval_list)
    print(oeis(eval_list))