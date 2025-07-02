### This file will test all methods and properties of the sorted sandpile.

### LIST OF TESTS
#   1   : verify the sorted sandpile correspondence for classical case
#   2   : verify the sorted sandpile correspondence for the new case k > 1

load("sorted_sandpile.py")
# Put a list of tests to run. To run all of them, add 0 to the list.
TEST_RUN = [1]

### SYMMETRIC FUNCTIONS SIDE

RR = FractionField(QQ['q','t'])
q,t = RR.gens()
Sym = SymmetricFunctions(RR)
e = Sym.e()
h = Sym.h()

def is_symmetric(poly):         # Check if polynomial is symmetric in q and t
    RR = FractionField(QQ['q','t'])  # type: ignore
    q,t,x = R.gens()
    poly1 = R(poly)     # Casting in new ring
    poly2 = poly1
    poly2(q=x)
    poly2(t=q)
    poly2(x=t)
    return (poly == poly2)


###     TEST    1   ####################################################

if 1 in TEST_RUN or 0 in TEST_RUN:
    print("\nStarting:\t TEST 1 - Clique-Independent sets and Shuffle Theorem")

    n = 5
    kmul = 1
    timing = []
    for mm in range(n+1):       # Run over all partitions that sum to n
        if mm == 0:
            Part1 = [[]]
            Part2 = Compositions(n-mm).list()
        elif mm == n:
            Part1 = Compositions(mm).list()
            Part2 = [[]]
        else:
            Part1 = Compositions(mm).list()
            Part2 = Compositions(n-mm).list()
        for mu in Part1:
            for nu in Part2:
                mu_l = list(mu)
                nu_l = list(nu)
                mu_p = list(mu)
                mu_p.sort(reverse=True)
                nu_p = list(nu)
                nu_p.sort(reverse=True)
                print("Testing compositions mu={}\t and nu={}".format(mu_l,nu_l))

                # Compute the qt-Poly
                S = Multi_CliqueIndependent_SortedSandpile(mu_l,nu_l,kmul)
                poly_sand = S.qt_Polynomial(opt = 2)
                # Compute the scalar product
                left = e([n])
                left = left.nabla(power=kmul)
                right = e(mu_p)*h(nu_p)
                poly_Hall = left.scalar(right)
                # Check
                check = (poly_sand == poly_Hall)
                if not check:
                    print("\tThe q,t-polynomial for mu={} and nu={} does NOT work!".format(mu,nu,check))
                    print("\tPolynomial symmetric functions:\t{}".format(poly_Hall))
                    print("\tPolynomial sandpiles:\t\t{}\n".format(poly_sand))

    print("Completed:\t TEST 1\n")

###     TEST    2   ####################################################

if 2 in TEST_RUN or 0 in TEST_RUN:
    print("\nStarting:\t TEST 2 - Multiplicity 2 Clique-Independent sets and a Rational Shuffle Theorem")

    n = 5
    kmul = 2
    timing = []
    for mm in range(n+1):       # Run over all partitions that sum to n
        if mm == 0:
            Part1 = [[]]
            Part2 = Compositions(n-mm).list()
        elif mm == n:
            Part1 = Compositions(mm).list()
            Part2 = [[]]
        else:
            Part1 = Compositions(mm).list()
            Part2 = Compositions(n-mm).list()
        for mu in Part1:
            for nu in Part2:
                print("Testing compositions mu={}\t and nu={}".format(mu,nu))
                mu_l = list(mu)
                nu_l = list(nu)
                mu_p = list(mu)
                mu_p.sort(reverse=True)
                nu_p = list(nu)
                nu_p.sort(reverse=True)

                # Compute the qt-Poly
                S = Multi_CliqueIndependent_SortedSandpile(mu_l,nu_l,kmul)
                poly_sand = S.qt_Polynomial(opt = 2)
                
                # Compute the scalar product
                left = e([n])
                left = left.nabla(power=kmul)
                right = e(mu_p)*h(nu_p)
                poly_Hall = left.scalar(right)

                # Check
                check = (poly_sand.numerator() == poly_Hall.numerator())
                if not check:
                    print("\tThe q,t-polynomial for mu={} and nu={} does NOT work!".format(mu,nu,check))
                    print("\tPolynomial symmetric functions:\t{}".format(poly_Hall.numerator()))
                    print("\tPolynomial sandpiles:\t\t{}\n".format(poly_sand))

    print("Completed:\t TEST 2\n")

###     TEST    3   ####################################################

print("\nStarting:\t TEST 3")



print("Completed:\t TEST 3\n")