#                           GENERALIZED CLIQUE-INDEPENDENT SORTED SANDPILES
#
#   These are some tests for particular structures of generalized clique-independent sorted sandpiles

from sage.combinat.q_analogues import qt_catalan_number
import time

exp_num = 5         #####   EXPERIMENT NUMBER

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

if exp_num == 3:        #   Try the loop of independents (a,n) = (-2,3) which has a 3 min computation time with option 1
    poly_list = []
    for (a,n,option) in [(-2,3,1),(-2,3,2)]:
        print("Computing case a={} and n={} with option={}.".format(a,n,option))
        t = time.time()
        ind_dict = {n:[0]}|{i:[i+1] for i in range(n)}  # Structure of gen independent/clique graph
        conf = {i:a for i in range(n+1)}                # All independents/cliques have "a" vertices
        G = Graph(ind_dict)
        S = General_CliqueIndependent_SortedSandpile(G, conf)
        q_poly = S.q_Polynomial(opt=option)
        poly_list.append(q_poly)
        print("Time elapsed: {} seconds".format(time.time() - t))
        print("(a,n) = ({},{}) has q-polynomial: {}".format(a,n,q_poly))
        print("Number sorted recurrents: {}\n".format(q_poly(1)))
    eval_list = [p(1) for p in poly_list]


if exp_num == 4:        #   Try the loop of independents (a,n) = (-2,3) which has a 3 min computation time with option 1
    poly_list = []
    for (a,n,option) in [(1,1,2),(1,2,2),(1,3,2),(1,4,2),(1,5,2)]:
        print("Computing case a={} and n={} with option={}.".format(a,n,option))
        t = time.time()
        ind_dict = {i:[i+1] for i in range(n)}          # FOR CHAIN
        #ind_dict = {n:[0]}|{i:[i+1] for i in range(n)}  # FOR LOOP
        conf = {i:a for i in range(n+1)}                # All independents/cliques have "a" vertices
        G = Graph(ind_dict)
        S = General_CliqueIndependent_SortedSandpile(G, conf)
        q_poly = S.qt_Polynomial(opt=option)
        poly_list.append(q_poly)
        print("Time elapsed: {} seconds".format(time.time() - t))
        print("(a,n) = ({},{}) has qt-polynomial: {}".format(a,n,q_poly))
        print("(a,n) = ({},{}) has q-polynomial: {}".format(a,n,q_poly(t=1)))
        print("(a,n) = ({},{}) has t-polynomial: {}".format(a,n,q_poly(q=1)))
        print("Sorted recurrents: {}\n".format(q_poly(q=1,t=1)))

if exp_num == 5:
    ind_dict = {0:[1], 1:[2,2]}#, 2:[3,3,3]}, 3:[4,4,4,4], 4:[5,5,5,5,5]}
    conf = {0:-2, 1:-2, 2:-2}# 3:-2}#, 4:1, 5:1}
    G = Graph(ind_dict)
    S = MultiGeneral_CliqueIndependent_SortedSandpile(G,conf)