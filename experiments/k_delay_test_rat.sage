#
#       SOME TESTS ON THE NEW DELAY STATISTIC
#
attach("./sorted_sandpile.py")


def level_rat(sandp, restr = 1):                                                ## Returns the level statistic of the configuration
    r"""
        Returns the level of the configuration.
    """
    not_inc_sink = sum([v[2] for v in sandp.sandpile_struct.to_undirected().edges()]) - sandp.sandpile_struct.out_degree(sandp.sink)
    # RESTRICTION OF LEVEL...
    return (- not_inc_sink + sandp.sandpile_config.deg())

def k_delay_rat(sandp, order = [], check_rec = True, restr = 1):             ## Returns the NEW delay statistic of the configuration
    r"""
        Given a reading order of nonsink vertices, this function computes the configuration's new conjectured delay.
        - k             : the multeplicity of all edges (except incident to the sink).
        - order         : the order for reading vertices. If undefined, the decreasing order on vertices is assumed.
        - check_rec     : this option can be used to override the is_recurrent() call.
    """
    if (not sandp.is_recurrent()) and check_rec:           # Check if the configuration is recurrent
        print(sandp.sandpile_config)
        print(sandp.sandpile_config.is_recurrent())
        sandp.sandpile_struct.show()
        raise Exception("The sorted configuration is not recurrent, hence delay is not defined.")
    k = max([edge[2] for edge in sandp.sandpile_struct.edges()])
    nonsink_vert = sandp.vertices
    n = len(nonsink_vert)

    if order == []:                                     # No order has been assigned: take decreasing order
        order = copy.copy(nonsink_vert)
        order.sort()
        order.reverse()

    toppl = [0]*n                                           # Stores information on how many partial topplings are still needed
    finalv = [k]*n                                          # We exit the loop when toppl == finalv
    delay = 0
    visit_vertex = 0
    visit_toppl = -1
    plus = 0
    delay_word = []
    sandp.topple_sink(sorting = False)                               # Start by toppling the sink
    while toppl != finalv:                               # Until everything has been toppled k times...
        for i in range(len(order)):
            #i = len(order) - 1 - i
            if (sandp.sandpile_struct.out_degree(order[i]) <= sandp.sandpile_config[order[i]]) and (toppl[i] == 0):
                                                                    # Can be toppled for the first time!
                #print("Toppl word {} and finalv {} for config {}".format(toppl,finalv,sandp.sandpile_config))
                sandp.single_topple(order[i], threshold = toppl[i], sorting = False)
                toppl[i] += 1
                visit_toppl += 1
                delay_word = delay_word + [100*order[i]]

                #incr = ceil(plus/restr)
                #incr = ceil((floor(visit_vertex/restr))/(n))
                #incr = floor(ceil(visit_vertex/restr)/n)
                #incr = ceil(visit_vertex/restr/n)
                #incr = ceil((floor(visit_vertex/n))/restr)
                #incr = floor(visit_toppl/restr)
                #incr = floor(visit_vertex/8)

                print("Toppl word {} and finalv {} for config {} and visit {} with restr {} and n {} and incr {}".format(toppl,finalv,sandp.sandpile_config,visit_vertex, restr, n, incr))
                delay += incr
            else:
                if toppl[i] < finalv[i] and toppl[i] > 0:                   # Topple later time...
                    #print("Toppl word {} and finalv {} for config {}".format(toppl,finalv,sandp.sandpile_config))
                    sandp.single_topple(order[i], threshold = toppl[i], sorting = False)
                    toppl[i] += 1
                    visit_toppl += 1
                    delay_word = delay_word + [order[i]]
            visit_vertex += 1
        plus += 1
    #print("The configuration {} has k-delay {}".format(sandp.sandpile_config, delay))
    print("Delay word: {}".format(delay_word))
    sandp.sort()
    return delay

def qt_Polynomial_rat(sandp, ordered = [], opt = 2, override = False, restr = 1):       ## Computes the q,t polynomial on (level, delay)
    r"""
        Returns the q,t - polynomial corresponding to the sorted sandpile's recurrent configurations.

        - order     : if specified, it fixes the reading order for the delay statistic.
        - opt       : if specified, it uses a different computing algorithm.
        - override  : if specified, the sorted recurrents are computed even if sorted_rec is non-empty.
    """
    R = FractionField(QQ['q, t'])      # type: ignore
    q,t = R.gens()
    poly = 0*q*t                            # Define the polynomial as 0

    if sandp.sorted_rec == [] or override:   # If sorted recurrents still to be computed...
        sandp.sorted_recurrents(option=opt)        # ...compute them!
    
    # TODO: be sure that the delay doesn't depend on the order in the same orbit...

    if ordered == []:                   # If there is no explicit order check for a specific one
        if sandp.specific_opt[0] == "clique-indep" or sandp.specific_opt[0] == "mul-clique-indep":          # The reading order that defines delay...
            ordered = sandp.specific_opt[2]
        elif sandp.specific_opt[0] == "gen-clique-indep" or sandp.specific_opt[0] == "mulgen-clique-indep":
            ordered = sandp.specific_opt[3]
        
    if sandp.specific_opt[0] == "mul-clique-indep":              # Compute the polynomial with k_delay
        # Restrict the sorted recurrents
        k_mul = sandp.specific_opt[3]
        topp_sink = {edge[1]:edge[2] for edge in sandp.sandpile_struct.edges(sandp.sink())}      # Compute addition of the sink grains
        restr_sort_rec = [conf for conf in sandp.sorted_rec if [v for v in sandp.vertices if (conf[v]+k_mul)%restr != 0] == []]
        q_shift = sum([k_mul*i%restr for i in range(len(sandp.vertices))])
        # Loop through all configurations
        for config in restr_sort_rec:
            sortedconfig = SandpileSortConfig(sandp.sandpile_struct, config, sandp.perm_group, sort = False, verts = sandp.vertices)
            q_exp = int((level_rat(sortedconfig, restr = restr) - q_shift)/ restr)
            #print("DELAY")
            #print("Config {}".format(sortedconfig.sandpile_config))
            t_exp = k_delay_rat(sortedconfig, order = ordered, check_rec=True, restr = restr)
            print("Configurazione {} con livello {} e delay {}\n".format(sortedconfig.sandpile_config, q_exp, t_exp))
            poly = poly + (q**q_exp) * (t**t_exp)
    else:                                                       # Compute the polynomial with regular delay
        # DO SAME THING
        for config in sandp.sorted_rec:
            sortedconfig = SandpileSortConfig(sandp.sandpile_struct, config, sandp.perm_group, sort = False, verts = sandp.vertices)
            #print("HERE!")
            q_exp = sortedconfig.level()
            t_exp = sortedconfig.k_delay(order = ordered, check_rec=False)
            #print("Configurazione {} con livello {} e delay {}".format(sortedconfig.sandpile_config, q_exp, t_exp))
            poly = poly + (q**q_exp) * (t**t_exp)
    
    return poly


def visual_degree(p):
    q_max = p.degrees()[0]
    t_max = p.degrees()[1]
    max_deg = max([q_max,t_max])
    # Compute the matrix
    coeff_matrix = [[p.coefficient([q_deg, t_deg]) for t_deg in range(max_deg + 1)] for q_deg in range(max_deg + 1)]
    # Plot the matrix
    return matrix_plot(matrix(coeff_matrix))


def is_symmetric(poly):         # Check if polynomial is symmetric in q and t
    R = FractionField(QQ['q','t','x'])  # type: ignore
    q,t,x = R.gens()
    poly1 = R(poly)     # Casting in new ring
    poly2 = poly1
    poly2 = poly2(q=x)
    poly2 = poly2(t=q)
    poly2 = poly2(x=t)
    return poly1 == poly2

n = 5
m = 7
restr = int(n/gcd(n,m))
k = int(m/gcd(n,m))

mu = [n]
nu = []

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
Q = SymAA.macdonald().Q()
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

S = Multi_CliqueIndependent_SortedSandpile(mu,nu,k,sinkmul=restr)

poly_sandpile = qt_Polynomial_rat(S, opt=2, restr = restr)

print("\nPolynomial")
print(poly_sandpile)
print("\nPolynomial's specializations")
print(poly_sandpile(t=1))
print(poly_sandpile(q=1))
print("\nPolynomial's coefficients")
print(poly_sandpile(t=1).coefficients())
print(poly_sandpile(q=1).coefficients())
print("\nIs it symmetric?")
print(is_symmetric(poly_sandpile))

p_sand = poly_sandpile.numerator()
plotted = visual_degree(p_sand)


#
#   Consider a different sandpile
#

#T = Multi_CliqueIndependent_SortedSandpile(mu,nu,k,hmul)

#poly_sandpile2 = T.qt_Polynomial(opt=2)

#print("Is the polynomial symmetric? {}".format(is_symmetric(poly_sandpile2)))