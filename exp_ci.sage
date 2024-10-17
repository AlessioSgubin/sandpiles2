### Testing the new class SandpileSortCofig
from sage.combinat.q_analogues import qt_catalan_number
import time

# Choose which test to run...
#test_number = 8

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
    mu = [2,1]
    nu = [2]
    Sym = SymmetricFunctions(FractionField(QQ['q','t']))
    #e = SymmetricFunctions(QQ['q','t']).e()
    #h = SymmetricFunctions(QQ['q','t']).h()
    e = Sym.e()
    h = Sym.h()
    left = e([n])
    left = left.nabla()
    right = e(mu)*h(nu)
    poly2 = left.scalar(right)
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

elif test_number == 4:
    S = CliqueIndependent_SortedSandpile([6],[])
    lis = S.sorted_recurrents()
    print(len(lis))

    G = graphs.CompleteGraph(6)
    T = SortedSandpile(G, 0, [[1,2,3,4,5]])
    lis2 = T.sorted_recurrents()
    print(len(lis2))

elif test_number == 5:              # Test per alcune partizioni
    n = 5
    for i in range(n-1):
        # Compute the qt-Poly
        S = CliqueIndependent_SortedSandpile([n-i-1,1],[i])
        poly_sand = S.qt_Polynomial()
        # Compute the scalar product
        Sym = SymmetricFunctions(FractionField(QQ['q','t']))
        e = Sym.e()
        h = Sym.h()
        left = e([n])
        left = left.nabla()
        right = e([n-i-1,1])*h([i])
        poly_Hall = left.scalar(right)

        check = (poly_sand == poly_Hall)
        print("Combinazione {},[{}] funziona? {}".format([n-i-1,1],i,check))

elif test_number == 6:              # Test su tutte le combinazioni per n fissato
    n = 6
    Sym = SymmetricFunctions(FractionField(QQ['q','t']))
    e = Sym.e()
    h = Sym.h()
    timing = []
    for m in range(n+1):
        if m == 0:
            Part1 = [[]]
            Part2 = Partitions(n-m).list()
        elif m == n:
            Part1 = Partitions(m).list()
            Part2 = [[]]
        else:
            Part1 = Partitions(m).list()
            Part2 = Partitions(n-m).list()
        for mu in Part1:
            for nu in Part2:
                if len(mu) + len(nu) <= 3:
                    t = time.time()         # Set timer
                    mu_l = list(mu)
                    nu_l = list(nu)
                    # Compute the qt-Poly
                    S = CliqueIndependent_SortedSandpile(mu_l,nu_l)
                    poly_sand = S.qt_Polynomial(opt = 2)
                    # Compute the scalar product
                    left = e([n])
                    left = left.nabla()
                    right = e(mu_l)*h(nu_l)
                    poly_Hall = left.scalar(right)
                    # Check
                    check = (poly_sand == poly_Hall)
                    print("The q,t-polynomial for mu={} and nu={} works? {}".format(mu,nu,check))
                    timing.append(time.time() - t)
    print("The time elapsed is {}".format(timing))

elif test_number == 7:              # Test per creazione del clique-indep generalizzato
    G = Graph({ 0:[1],  1:[2]})
    conf = {0: 2, 1:2, 3:2}
    S = General_CliqueIndependent_SortedSandpile(G, conf)

elif test_number == 8:              # Test su tutte le combinazioni per n fissato, per edge multipli
    n = 5
    kmul = 2
    hmul = 2
    Sym = SymmetricFunctions(FractionField(QQ['q','t']))
    e = Sym.e()
    h = Sym.h()
    timing = []
    for m in range(n+1):
        if m == 0:
            Part1 = [[]]
            Part2 = Partitions(n-m).list()
        elif m == n:
            Part1 = Partitions(m).list()
            Part2 = [[]]
        else:
            Part1 = Partitions(m).list()
            Part2 = Partitions(n-m).list()
        for mu in Part1:
            for nu in Part2:
                if len(mu) + len(nu) <= 3:
                    t = time.time()         # Set timer
                    mu_l = list(mu)
                    nu_l = list(nu)
                    # Compute the qt-Poly
                    S = Multi_CliqueIndependent_SortedSandpile(mu_l,nu_l,kmul,hmul)
                    poly_sand = S.qt_Polynomial(opt = 2)
                    '''
                    # Compute the scalar product
                    left = e([n])
                    left = left.nabla()
                    right = e(mu_l)*h(nu_l)
                    poly_Hall = left.scalar(right)
                    # Check
                    check = (poly_sand == poly_Hall)
                    print("The q,t-polynomial for mu={} and nu={} works? {}".format(mu,nu,check))
                    '''
                    poly1 = poly_sand(q = 1)
                    poly2 = poly_sand(t = 1)
                    check = (poly1 == poly2)
                    print("Partitions mu = {} and nu = {} is: {}".format(mu, nu, check))
                    print("Polynomials: {}\t\t\t{}\n".format(poly1,poly2))
                    #print("Using partitions mu = {} and nu = {} we get:\n\t{}".format(mu,nu,poly_sand))
                    timing.append(time.time() - t)
    print("The time elapsed is {}".format(timing))

elif test_number == 9:              # Test su tutte le combinazioni per n fissato, per edge multipli
    n = 5
    kmul = 3
    hmul = 3
    R.<q> = PolynomialRing(QQbar)
    T.<t> = PolynomialRing(QQbar)
    Sym = SymmetricFunctions(FractionField(QQ['q','t']))
    e = Sym.e()
    h = Sym.h()
    timing = []
    for m in range(n+1):
        if m == 0:
            Part1 = [[]]
            Part2 = Partitions(n-m).list()
        elif m == n:
            Part1 = Partitions(m).list()
            Part2 = [[]]
        else:
            Part1 = Partitions(m).list()
            Part2 = Partitions(n-m).list()
        for mu in Part1:
            for nu in Part2:
                if len(mu) + len(nu) <= 3:
                    t = time.time()         # Set timer
                    mu_l = list(mu)
                    nu_l = list(nu)
                    # Compute the qt-Poly
                    S = Multi_CliqueIndependent_SortedSandpile(mu_l,nu_l,kmul,hmul)
                    poly_sand = S.qt_Polynomial(opt = 2)
                    '''
                    # Compute the scalar product
                    left = e([n])
                    left = left.nabla()
                    right = e(mu_l)*h(nu_l)
                    poly_Hall = left.scalar(right)
                    # Check
                    check = (poly_sand == poly_Hall)
                    print("The q,t-polynomial for mu={} and nu={} works? {}".format(mu,nu,check))
                    '''
                    poly1 = poly_sand(q = 1)
                    poly2 = poly_sand(t = 1)
                    check = (poly1 == poly2)
                    print("Partitions mu = {} and nu = {} is: {}".format(mu, nu, check))

                    # Compute the degrees
                    deg1 = poly1.degree()
                    deg2 = poly2.degree()
                    int_poly1 = T(poly1)
                    int_poly2 = R(poly2)
                    coeffs1 = int_poly1.coefficients(sparse = False)
                    coeffs2 = int_poly2.coefficients(sparse = False)
                    mindeg1 = min([i for i in range(len(coeffs1)) if coeffs1[i] != 0])
                    mindeg2 = min([i for i in range(len(coeffs2)) if coeffs2[i] != 0])

                    print("Max/min q-degree: {}\t\t Max/min t-degree: {}\n".format((deg2, mindeg2), (deg1, mindeg1)))
                    timing.append(time.time() - t)
    print("The time elapsed is {}".format(timing))