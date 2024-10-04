### Testing the new class SandpileSortCofig

test_number = 3                     # Calcolo configurazioni Sandpile

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
    S = SortedSandpile(G, 0, [[1,2,3,4]])

    lis = S.sorted_recurrents()
    #print("In tutto sono {} elementi".format(len(lis)))

    poly = S.qt_Polynomial()
    print("Il polinomio associato al Sorted Sandpile è {}".format(poly))

elif test_number == 3:              # Computare il qt-polynomial di clique-independent
    n = 2 + 1 + 2 + 1
    mu = [2,1]
    nu = [2,1]
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
    S.show()
    print("Calcoliamo le configurazioni...")
    S.sorted_recurrents()
    print("Calcoliamo il polinomio...")
    poly = S.qt_Polynomial()
    print("Il polinomio associato al Sorted Sandpile è {}".format(poly))

    check = (poly == poly2)
    print("I due polinomi sono uguali? {}".format(check))