### Testing the new class SandpileSortCofig

test_number = 3

if test_number == 0:

    G = graphs.CompleteGraph(5)

    S = Sandpile(G, 0)

    print(S.sink())
    C = SandpileConfig(S, [0,1,2,3])

    print(C.is_recurrent())

elif test_number == 1:

    G = graphs.CompleteGraph(5)
    S = Sandpile(G, 0)

    C1 = SandpileSortConfig(S, [0,1,2,3], [[1],[2],[3,4]])
    C2 = SandpileSortConfig(S, [0,1,3,2], [[1],[2],[3,4]])
    C3 = SandpileSortConfig(S, [0,3,1,2], [[1],[2],[3,4]])
    
    print(C1 == C2)
    print(C1 == C3)

    print(C1.level())
    print(C1.is_recurrent())
    print(C1.delay())

    print(C2.sort())

elif test_number == 2:

    G = graphs.CompleteGraph(6)
    S = SortedSandpile(G, 0, [[1,2,3,4]])

    lis = S.sorted_recurrents()
    #print("In tutto sono {} elementi".format(len(lis)))

    poly = S.qt_Polynomial()
    print("Il polinomio associato al Sorted Sandpile è {}".format(poly))

elif test_number == 3:

    S = CliqueIndependent_SortedSandpile([4,3],[3,2])
    poly = S.qt_Polynomial()
    print("Il polinomio associato al Sorted Sandpile è {}".format(poly))