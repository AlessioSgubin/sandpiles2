### Testing the new class SandpileSortCofig
from sage.combinat.q_analogues import qt_catalan_number


#test_number = 2                     # Calcolo configurazioni Sandpile

if test_number == 0:

    G = graphs.CompleteGraph(5)

    S = Sandpile(G, 0)

    print(S.sink())
    C = SandpileConfig(S, [0,1,2,3])

    print(C.is_recurrent())

elif test_number == 1:              # Controllo della funzionalità di Sorted Sandpile Configurations

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

elif test_number == 2:              # Computare il qt-polynomial di grafi completi...

    G = graphs.CompleteGraph(6)
    S = SortedSandpile(G, 0, [[1,2,3,4,5]])

    lis = S.sorted_recurrents()
    print("In tutto sono {} elementi".format(len(lis)))

    poly = S.qt_Polynomial()
    print("Il polinomio associato al Sorted Sandpile è {}".format(poly))

elif test_number == 3:              # Computare il qt-polynomial di clique-independent
    n = 5
    mu = [n]
    nu = []
    Sym = SymmetricFunctions(FractionField(QQ['q','t']))
    e = Sym.elementary()
    h = Sym.homogeneous()
    left = e([n]).nabla()
    right = e(mu)*h(nu)
    poly2 = left.scalar_qt(right)
    #print("Fattore a sinistra è {}".format(left))
    #print("Fattore a destra è {}".format(right))
    print("Il polinomio ottenuto col prodotto scalare è {}".format(poly2))

    S = CliqueIndependent_SortedSandpile(mu, nu)
    #S.show()
    print("Calcoliamo le configurazioni...")
    S.sorted_recurrents()
    print(S.sorted_rec)
    print("Calcoliamo il polinomio...")
    poly = S.qt_Polynomial()
    print("Il polinomio associato al Sorted Sandpile è {}".format(poly))

    check = (poly == poly2)
    print("I due polinomi sono uguali? {}".format(check))
    check = (poly == qt_catalan_number(n))
    print("I due polinomi sono uguali? {}".format(check))

elif test_number == 4:
    S = CliqueIndependent_SortedSandpile([6],[])
    lis = S.sorted_recurrents()
    print(len(lis))

    G = graphs.CompleteGraph(6)
    T = SortedSandpile(G, 0, [[1,2,3,4,5]])
    lis2 = T.sorted_recurrents()
    print(len(lis2))

elif test_number == 5:
    seq = []
    for i in range(5):
        S = CliqueIndependent_SortedSandpile([i+1],[])
        lis = S.sorted_recurrents()
        seq = seq + [len(lis)]
        print(seq)
    oeis(seq)