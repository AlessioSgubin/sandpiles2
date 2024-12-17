#
#               SOME OTHER EXPERIMENTS
#
test = 5

###     FUNCTIONS
def k_ladder(n,k):                      # This function computes and returns the k-ladder
    return [k*(n-i) for i in range(n)]

def argmin(iterable):
    return min(enumerate(iterable), key=lambda x: x[1])[0]

#
### STATISTICS FOR PATH
def decomp_path(n,k,w,dyck):
    if k == 1:
        if max(w) <= n:
            return [[w]]
        else:
         return []
    else:
        possible = []
        k_ladd = k_ladder(n,k-1)
        for d in dyck:
            check = [w[i]-d[i]-k_ladd[i] for i in range(n)]
            if min(check) >= 0:
                w_less = [w[i]-d[i] for i in range(n)]
                recurs = decomp_path(n,k-1,w_less,dyck)
                new_poss = [[d] + t for t in recurs]
                possible = possible + new_poss
        return possible

def area(n,k,w):
    return sum([w[i]-n*k+i*k for i in range(n)])

def dinv(w):
    n = len(w)
    pairs = []
    u = [w[i] - n + i for i in range(n)]
    for i in range(n-1):
        for j in range(n-i-1):
            if (u[i] == u[i+j+1] or u[i] == u[i+j+1]+1):
                pairs.append((i,j))
    return len(pairs)

def bounce(w):
    n = len(w)
    b_path = 0
    bounce = 0
    add_b = 0
    u = [w[i] - n + i for i in range(n)]
    for i in range(n):
        if u[i] >= b_path:    # bounce_path continues up
            b_path += 1
        else:                 # bounce_path hits path
            b_path = 1
            add_b += 1
        bounce += add_b
    return bounce

def min_bounce(n,k,w,dyck):
    bounce_list = []
    for decomp in decomp_path(n,k,w,dyck):
        bounce_dec = [bounce(p) for p in decomp]
        bounce_list.append(sum(bounce_dec))
    return min(bounce_list)

def min_bounce_info(n,k,w,dyck):
    bounce_list = []
    decompos = decomp_path(n,k,w,dyck)
    for decomp in decompos:
        bounce_dec = [bounce(p) for p in decomp]
        bounce_list.append(sum(bounce_dec))
    return (min(bounce_list), decompos[argmin(bounce_list)], len(decompos), argmin(bounce_list))

def pdinv(n,k,w):
    dinv = 0
    # Read all cells above path
    for i in range(n):                  # Loop on rows
        for j in range(n*k-w[i]):       # Loop on columns
            arm = (n*k - w[i]) - (j+1)
            leg = i - max([h for h in range(n) if n*k-w[h]-1 < j]) - 1
            if (arm <= k*(leg+1)) and (k*leg < arm + 1):
                dinv += 1
    return dinv

def pdinv_noedge(n,k,w):
    dinv = 0
    # Read all cells above path
    for i in range(n):                  # Loop on rows
        for j in range(n*k-w[i]):       # Loop on columns
            arm = (n*k - w[i]) - (j+1)
            leg = i - max([h for h in range(n) if n*k-w[h]-1 < j]) - 1
            if (arm <= k*(leg+1)-1) and (k*leg < arm + 1):
                dinv += 1
    return dinv

def diag_touch(n,k,path):
    touch = 0
    for i in range(n):  # in each row check area...
        if path[i] == (n-i)*k:
            touch += 1
    return touch

#
### PATH FUNCTIONS
def rows_to_cols(n,k,w):        # Path written as rows to cols
    return [max([i for i in range(n) if n*k <= w[i]+j])+1  for j in range(n*k)]

def cols_to_rows(n,k,w):        # Path written as cols to rows
    return [max([j for j in range(n*k) if w[n*k-j-1]>i])+1 for i in range(n)]

#
### k-DYCK PATHS
def kDyckPaths(n,k):
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


#
###     ORDERING OF TABLEAUX FOR INCREASING DELAY
if test == 0:
    import functools

    def lex_gt(a,b):        # Lexicographic order of lists
        ind = 0
        cont = True
        check = 0
        while ind < len(a) and cont:
            if a[ind] < b[ind]:
                check = 1
                cont = False
            elif a[ind] > b[ind]:
                check = -1
                cont = False
            else:
                ind += 1
        return check

    def delay(w):           # Computes delay on a word
        res = 0
        loops = 0
        prev = 0
        ran = [i for i in range(len(w))]
        ran.reverse()
        for i in ran:
            nex = w.index(i)
            if prev <= nex:     # Still same loop
                res += loops
            else:               # New loop
                loops += 1
                res += loops
            prev = nex
        return res

    mu = [4,2]
    n = sum(mu)
    list_inj = []

    for sigma in Permutations(range(n)):
        # Sort all orbits
        sort_sigma = []
        index = 0
        for part in mu:
            temp = [sigma[index + i] for i in range(part)]
            temp.sort()
            temp.reverse()
            sort_sigma = sort_sigma + temp
            index += part
        # Append if necessary
        if sort_sigma not in list_inj:
            list_inj.append(sort_sigma)

    list_inj.sort(key=functools.cmp_to_key(lex_gt))

    list_delays = [delay(w) for w in list_inj]
    max_del = max(list_delays)
    list_order = []

    for i in range(max_del+1):
        list_order = list_order + [w for w in list_inj if delay(w) == i]

    print(list_order)
    print([delay(w) for w in list_order])

    # Convert all to tableaux
    list_tabl = []
    for w in list_order:
        tabl = []
        ind = 0
        for part in mu:
            tabl.append([w[ind + i] for i in range(part)])
            ind += part
        list_tabl.append(Tableau(tabl))

    for i in range(len(list_tabl)):
        list_tabl[i].pp()
        print("Delay is {}\n".format(delay(list_order[i])))


#
###     MIN DELAY
if test == 1:

    n = 4
    k = 2

    # Compute all D(n,k)
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

    # Compute all D(n,1)
    dyck_paths = [[n]]
    for i in range(n-1):
        new_paths = []
        for w in dyck_paths:
            for j in range(w[i]-n+i+2):
                temp = w + [w[i]-j]
                new_paths.append(temp)
        dyck_paths = new_paths

    R = FractionField(QQ['q, t, x'])
    q,t,x = R.gens()
    p = 0*q*t
    for w in dyck_k_paths:
        term = q**(area(n,k,w)) * t**(min_bounce(n,k,w,dyck_paths))
        p = p + term
    
    p_symm = p(q=x)
    p_symm = p_symm(t=q)
    p_symm = p_symm(x=t)


    # Symmetric functions
    Sym = SymmetricFunctions(FractionField(QQ['q','t']))
    e = Sym.e()
    left = e([n])
    left = left.nabla(power = k)
    right = e([n])
    poly = left.scalar(right)

    print(p_symm(q=1,t=1))


if test == 2:
    n = 4
    k = 2

    # Compute all D(n,k)
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

    '''
    #   ASSOCIATED POLY
    #   Compute the polynomial associated to the Dyck paths
    R = FractionField(QQ['q,t,x'])
    q,t,x = R.gens()

    poly = 0*q*t*x
    for p in dyck_k_paths:
        poly += q**(pdinv(n,k,p)) * t**(area(n,k,p))
    poly2 = poly
    poly2 = poly2(q=x)
    poly2 = poly2(t=q)
    poly2 = poly2(x=t)
    #print(poly)
    #print(poly2)
    #print(poly == poly2)

    #   SYMMETRIC FUNCTIONS
    #   Check that the formula agrees with the symmetric function
    Sym = SymmetricFunctions(FractionField(QQ['q','t']))
    e = Sym.e()
    left = e([n])
    left = left.nabla(power = k)
    right = e([n])
    qpoly = left.scalar(right)
    qpoly2 = qpoly
    qpoly2 = qpoly2(q=x)
    qpoly2 = qpoly2(t=q)
    qpoly2 = qpoly2(x=t)
    #print(qpoly == qpoly2)
    '''
    #   TRY DECOMPOSITION
    #   Check if a formula works
    for w in dyck_k_paths:
        # Decompose into right number of parts
        w_col = rows_to_cols(n,k,w)
        decomp_col = [[w_col[k*i+h] for i in range(n)] for h in range(k)]
        decomp_row = [cols_to_rows(n,1,w) for w in decomp_col]
        # Compute
        val_pdinv = pdinv(n,k,w)
        decomp_dinv = pdinv(n,1,decomp_row[0])
        for i in range(k-1):
            decomp_dinv = decomp_dinv + pdinv_noedge(n,1,decomp_row[i+1])
        if val_pdinv != decomp_dinv:
            print("Path {} has pdinv {} and decomp_dinv {}".format(w, val_pdinv, decomp_dinv))

if test == 3:
    n = 4
    k = 3

    # Compute all D(n,k)
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

    # Compute all D(n,1)
    dyck_paths = [[n]]
    for i in range(n-1):
        new_paths = []
        for w in dyck_paths:
            for j in range(w[i]-n+i+2):
                temp = w + [w[i]-j]
                new_paths.append(temp)
        dyck_paths = new_paths

    for w in dyck_k_paths:
        # Compute infos
        val,dec,a,b = min_bounce_info(n,k,w,dyck_paths)
        # Print infos
        print("Path {} decomposed as {} has minbounce {}. Between {} chose {}".format(w,dec,val, a, b))

if test == 4:
    #
    #       Enumeration of P_w,u
    #
    w = [0,0,1,0,1]
    u = [0,1,1]

    k = 2
    m  = len(w) + len(u)
    good_list = [[w[0]]]

    for lop in range(m-1):
        new_list = []
        for part in good_list:
            pos_w = len([i for i in part if i < k])
            pos_u = len([i for i in part if i >= k])
            if pos_w < len(w):                                                                      # Can still add a letter from w
                new_list.append(part + [w[pos_w]])
            if (pos_u < len(u)) and (part[len(part)-1] >= k):                                       # Last letter in part was from u
                new_list.append(part + [k + u[pos_u]])
            if (pos_u < len(u)) and (part[len(part)-1] < k) and (part[len(part)-1] >= u[pos_u]):    # From w and right descent      
                new_list.append(part + [k + u[pos_u]])
        good_list = new_list

    print(good_list)
    print(len(good_list))


if test == 5:       # Try the decompositions with En,k

    # Compute E_n,k                         ####################################
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

    # Compute all D(n,k)                    ####################################
    n = 5
    k = 2
    dyck_k_paths = kDyckPaths(n,k)

    # Try to see if we get any recurrence   ####################################
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

        print("For n = {}, k = {} and h = {}:".format(n,k,h))
        print("\t Hall: {}".format(nkh_Hall_poly))
        print("\t Dyck: {}\n".format(nkh_dyck_poly))

