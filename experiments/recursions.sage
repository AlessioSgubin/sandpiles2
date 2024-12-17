#
# Auxiliary functions                       ####################################
#
from sage.combinat.q_analogues import q_multinomial

def kDyckPaths(n,k):                    # kDYCKPATHS: computes the set of all (nk,n)-Dyck paths
    possible_paths = [[1]]
    for i in range(n*(k+1)-1):
        new_possible = []
        for part_path in possible_paths:
            height = sum(part_path)
            distance = sum([1-part_path[j] for j in range(i+1)])
            if height == n:                 # Reached max height, just 0's
                temp = deepcopy(part_path)
                temp.append(0)
                new_possible.append(temp)
            elif k*height < distance + 1:   # Reached the "diagonal", just a 1
                temp = deepcopy(part_path)
                temp.append(1)
                new_possible.append(temp)
            else:
                temp = deepcopy(part_path)
                temp.append(0)
                new_possible.append(temp)
                temp = deepcopy(part_path)
                temp.append(1)
                new_possible.append(temp)
        possible_paths = new_possible
    dyck_k_paths = [[sum([1-w[i+j] for j in range(len(w)-i)]) for i in range(len(w)) if w[i]==1] for w in possible_paths]
    return dyck_k_paths

def area(n,k,w):                        # AREA: computes the area of a given path
    return sum([w[i]-n*k+i*k for i in range(n)])

def pdinv(n,k,w):                       # PDINV: computes the path dinv of a given path
    dinv = 0
    # Read all cells above path
    for i in range(n):                  # Loop on rows
        for j in range(n*k-w[i]):       # Loop on columns
            arm = (n*k - w[i]) - (j+1)
            leg = i - max([h for h in range(n) if n*k-w[h]-1 < j]) - 1
            if (arm <= k*(leg+1)) and (k*leg < arm + 1):
                dinv += 1
    return dinv

def diag_touch(n,k,path,m=0):               # DIAG_TOUCH: computes the times the path touches the diagonal (= modulo m)
    touch = 0
    for i in range(n):  # in each row check area...
        if path[i] == (n-i)*k+m:
            touch += 1
    return touch

def gen_poly(n,k,set_paths,ring):
    q = ring.gens()[0]
    t = ring.gens()[1]
    poly = 0*q*t
    for path in set_paths:
        poly = poly + (q**pdinv(n,k,path))*(t**area(n,k,path))
    return poly


#
# Compute E_n,k                         ####################################
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

def C_alpha(alpha, f=SymAA.schur()[0]):
    # Zabrocki's operator C_\alpha.
    if len(alpha) == 0:
        return f
    else:
        sf = sum((-q)**(1-alpha[-1]) * q**(-r) * SymAA.homogeneous()[alpha[-1]+r] *
                f.skew_by(SymAA.homogeneous()[r](SymAA.schur()[1]*(1-q))) for r in range(f.degree()+1))
        return C_alpha(alpha[:-1], sf)

def E(n, k):
    # The E_{n,k} symmetric functions.
    return sum([C_alpha(alpha, SymAA.schur()[0]) for alpha in Compositions(n) if len(alpha) == k])

#
# Compute all D(n,k)                    ####################################
#
n = 5
k = 2
dyck_k_paths = kDyckPaths(n,k)

#
# Try to see if we get any recurrence...####################################
#
for h in range(n):              # See if we get a relation between En,k and subsets of kDyckPaths
    h += 1
    # Compute poly from dyck paths
    nkh_dyck_poly = 0*q*t*u
    for path in dyck_k_paths:
        if diag_touch(n,k,path) == h:
            nkh_dyck_poly = nkh_dyck_poly + (t**(area(n,k,path))*(q**(pdinv(n,k,path))))
    # Compute poly from symmetric functions
    left = E(n,h)
    left = left.nabla(power = k)
    nkh_Hall_poly = left.scalar(e([n]))

    # Print the results
    #print("For n = {}, k = {} and h = {}: equality holds? {}".format(n,k,h, nkh_dyck_poly==nkh_Hall_poly))
    #print("\t Hall: {}".format(nkh_Hall_poly))
    #print("\t Dyck: {}\n".format(nkh_dyck_poly))

#
# Try to see if we get a hint on the actual recursion?
#
n = 6
k = 2
h0 = 2
h1 = 1
# All (nk,n)-Dyck paths
dyck1 = kDyckPaths(n,k)
dyck1_touch = [path for path in dyck1 if diag_touch(n,k,path) == h0 and diag_touch(n,k,path,1) == h1]
# Associated poly
poly = gen_poly(n,k,dyck1_touch,AA)

# All ((n-h0)k, n-h0)
dyck2 = kDyckPaths(n,k)
poly_rec = 0*q*t
for h2 in range(n-h0-h1+1):
    dyck2_touch = [path for path in dyck2 if diag_touch(n,k,path) == h0+h1 and diag_touch(n,k,path,1) == h2]
    poly_rec = poly_rec + q**(n-h0) * t**(binomial(h0,2)) * q_multinomial([h0-1,h1,h2], q=t) * gen_poly(n,k,dyck2_touch,AA)

print(poly(t=1))
print(poly(q=1))
print(poly_rec(t=1))
print(poly_rec(q=1))