#   This code computes the sorted recurrent configurations and their q,t polynomials and stores them in a database
#   
#   We are following this naming scheme:
#       -   CHC_XX_YY    : a chain of cliques, sorted sandpile with XX = a and YY = n
#       -   CHI_XX_YY    : a chain of independents, sorted sandpile with XX = a and YY = n
#       -   LPC_XX_YY    : a loop of cliques, sorted sandpile with XX = a and YY = n
#       -   LPI_XX_YY    : a loop of independents, sorted sandpile with XX = a and YY = n

import pickle
import time
import os.path

type_test = "CHC"   # Type of sandpile
verbose = True     # Print stats about the sandpile
override = False    # If file already present, re-compute or not?

if type_test == "CHC":
    for test in [(1,2,2),(2,2,2),(1,3,2),(1,4,2),(1,5,2),(1,6,2),(2,3,2),(3,2,2)]:
        a = test[0]
        n = test[1]
        option = test[2]
        namefile = str("CHC_{0:02d}".format(a) + "_{0:02d}".format(n))
        fpath = str("database/" + namefile + ".ssdp")
        if override or (not os.path.isfile(str(fpath))):
            print(str("Computing the sorted sandpile coded " + namefile))
            t = time.time()
            ind_dict = {i:[i+1] for i in range(n)}          # FOR CHAIN
            conf = {i:a for i in range(n+1)}                # All independents/cliques have "a" vertices
            G = Graph(ind_dict)
            S = General_CliqueIndependent_SortedSandpile(G, conf)

            q_poly = S.qt_Polynomial(opt=option)

            print("Time elapsed: {} seconds".format(time.time() - t))
            if verbose:
                print("(a,n) = ({},{}) has qt-polynomial: {}".format(a,n,q_poly))
                print("(a,n) = ({},{}) has q-polynomial: {}".format(a,n,q_poly(t=1)))
                print("(a,n) = ({},{}) has t-polynomial: {}".format(a,n,q_poly(q=1)))
            print("Sorted recurrents: {}".format(q_poly(q=1,t=1)))

            S.save(fpath)
            print("Saved file {}\n".format(fpath))
        else:
            print("The sorted sandpile {} is already in the database.\n".format(namefile))

#type_test = "CHI"

if type_test == "CHI":
    for test in [(2,1,2),(2,2,2),(2,3,2),(3,1,2),(3,2,2),(4,1,2)]:
        a = test[0]
        n = test[1]
        option = test[2]
        namefile = str("CHI_{0:02d}".format(a) + "_{0:02d}".format(n))
        fpath = str("database/" + namefile + ".ssdp")
        if override or (not os.path.isfile(str(fpath))):
            print(str("Computing the sorted sandpile coded " + namefile))
            t = time.time()
            ind_dict = {i:[i+1] for i in range(n)}          # FOR CHAIN
            conf = {i:-a for i in range(n+1)}                # All independents/cliques have "a" vertices
            G = Graph(ind_dict)
            S = General_CliqueIndependent_SortedSandpile(G, conf)

            q_poly = S.qt_Polynomial(opt=option)

            print("Time elapsed: {} seconds".format(time.time() - t))
            if verbose:
                print("(a,n) = ({},{}) has qt-polynomial: {}".format(a,n,q_poly))
                print("(a,n) = ({},{}) has q-polynomial: {}".format(a,n,q_poly(t=1)))
                print("(a,n) = ({},{}) has t-polynomial: {}".format(a,n,q_poly(q=1)))
            print("Sorted recurrents: {}".format(q_poly(q=1,t=1)))

            S.save(fpath)
            print("Saved file {}\n".format(fpath))
        else:
            print("The sorted sandpile {} is already in the database.\n".format(namefile))

#type_test = "LPC"

if type_test == "LPC":
    for test in [(2,1,2),(2,2,2),(2,3,2),(3,1,2),(4,1,2)]:
        a = test[0]
        n = test[1]
        option = test[2]
        namefile = str("LPC_{0:02d}".format(a) + "_{0:02d}".format(n))
        fpath = str("database/" + namefile + ".ssdp")
        if override or (not os.path.isfile(str(fpath))):
            print(str("Computing the sorted sandpile coded " + namefile))
            t = time.time()
            ind_dict = {n:[0]} | {i:[i+1] for i in range(n)}          # FOR LOOP
            conf = {i:-a for i in range(n+1)}                # All independents/cliques have "a" vertices
            G = Graph(ind_dict)
            S = General_CliqueIndependent_SortedSandpile(G, conf)

            q_poly = S.qt_Polynomial(opt=option)

            print("Time elapsed: {} seconds".format(time.time() - t))
            if verbose:
                print("(a,n) = ({},{}) has qt-polynomial: {}".format(a,n,q_poly))
                print("(a,n) = ({},{}) has q-polynomial: {}".format(a,n,q_poly(t=1)))
                print("(a,n) = ({},{}) has t-polynomial: {}".format(a,n,q_poly(q=1)))
            print("Sorted recurrents: {}".format(q_poly(q=1,t=1)))

            S.save(fpath)
            print("Saved file {}\n".format(fpath))
        else:
            print("The sorted sandpile {} is already in the database.\n".format(namefile))

#type_test = "LPI"

if type_test == "LPI":
    for test in [(2,1,2),(2,2,2),(2,3,2),(3,1,2),(4,1,2),(5,1,2),(3,2,2)]:
        a = test[0]
        n = test[1]
        option = test[2]
        namefile = str("LPI_{0:02d}".format(a) + "_{0:02d}".format(n))
        fpath = str("database/" + namefile + ".ssdp")
        if override or (not os.path.isfile(str(fpath))):
            print(str("Computing the sorted sandpile coded " + namefile))
            t = time.time()
            ind_dict = {n:[0]} | {i:[i+1] for i in range(n)}          # FOR LOOP
            conf = {i:-a for i in range(n+1)}                # All independents/cliques have "a" vertices
            G = Graph(ind_dict)
            S = General_CliqueIndependent_SortedSandpile(G, conf)

            q_poly = S.qt_Polynomial(opt=option)

            print("Time elapsed: {} seconds".format(time.time() - t))
            if verbose:
                print("(a,n) = ({},{}) has qt-polynomial: {}".format(a,n,q_poly))
                print("(a,n) = ({},{}) has q-polynomial: {}".format(a,n,q_poly(t=1)))
                print("(a,n) = ({},{}) has t-polynomial: {}".format(a,n,q_poly(q=1)))
            print("Sorted recurrents: {}".format(q_poly(q=1,t=1)))

            S.save(fpath)
            print("Saved file {}\n".format(fpath))
        else:
            print("The sorted sandpile {} is already in the database.\n".format(namefile))
