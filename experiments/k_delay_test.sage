#
#       SOME TESTS ON THE NEW DELAY STATISTIC
#
attach("./sorted_sandpile.py")

def is_symmetric(poly):         # Check if polynomial is symmetric in q and t
    R = FractionField(QQ['q','t','x'])  # type: ignore
    q,t,x = R.gens()
    poly1 = R(poly)     # Casting in new ring
    poly2 = poly1
    poly2(q=x)
    poly2(t=q)
    poly2(x=t)
    return (poly == poly2)

mu = [2]
nu = [2,2]
n = sum(mu) + sum(nu)
k = 2

#
#   Compute E_n,k                         ####################################
#

AA=QQ['q', 't', 'u'].fraction_field()
AA.inject_variables(verbose=False)
q = AA.gens()[0]
t = AA.gens()[1]
u = AA.gens()[2]

SymAA = SymmetricFunctions(AA)
SymAA.inject_shorthands(verbose=False)
H = SymAA.macdonald().Ht()
e = SymAA.e()

def C_alpha(alpha, f=SymAA.schur()[0]):         # Zabrocki's operator C_\alpha.
    if len(alpha) == 0:
        return f
    else:
        sf = sum((-q)**(1-alpha[-1]) * q**(-r) * SymAA.homogeneous()[alpha[-1]+r] * f.skew_by(SymAA.homogeneous()[r](SymAA.schur()[1]*(1-q))) for r in range(f.degree()+1))
        return C_alpha(alpha[:-1], sf)

def E(n, k):                                    # The E_{n,k} symmetric functions.
    return sum([C_alpha(alpha, SymAA.schur()[0]) for alpha in Compositions(n) if len(alpha) == k])

#
#   Compute the q,t Hall scalar product   ###################################
#

left = e([n])
left = left.nabla(power=k)
right = e(mu)*h(nu)
poly_symmetric = left.scalar(right)

#
#   Compute delay of sandpile             ###################################
#

S = Multi_CliqueIndependent_SortedSandpile(mu,nu,k)

poly_sandpile = S.qt_Polynomial(opt=2)

print("Are the polynomials equal? {}".format(poly_symmetric == poly_sandpile))