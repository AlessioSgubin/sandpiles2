### Testing new delay
from sage.combinat.q_analogues import qt_catalan_number

def new_delay(sandconf, n, k):                      # Given a configuration, it computes the delay
    #print("Starting Config is {}".format(sandconf.sandpile_config))
    toppl = [k for i in range(n)]                   # Stores information on how many topplings
    nullv = [0 for i in range(n)]
    order = copy.copy(sandconf.vertices)
    order.sort()
    delay = 0
    plus = 0
    sandconf.topple_sink(sorting = False)
    while toppl != nullv:                           # Until everything has been toppled k times...
        #print(sandconf.sandpile_config)
        for i in range(len(order)):
            if (sandconf.sandpile_struct.out_degree(order[i]) <= sandconf.sandpile_config[order[i]]) and (toppl[i] == k):
                                                            # Can be toppled for the first time!
                sandconf.single_topple(order[i], sorting = False)
                toppl[i] -= 1
                delay += plus
                #print("Topple! Now delay = {} and plus = {}".format(delay, plus))
            else:
                if toppl[i] < k and toppl[i] > 0:
                                                            # Topple later time...
                    sandconf.single_topple(order[i], sorting = False)
                    toppl[i] -= 1
                    #print("Delayed toppling")
        plus += 1
    return delay

def qt_poly(n,k,sandp):
    R = FractionField(QQ['q, t'])      # type: ignore
    q,t = R.gens()
    poly = 0*q*t
    print("Sorted")
    sandp.sorted_recurrents(option = 2)
    print("Computing")
    for conf in sandp.sorted_rec:
        sortedconfig = SandpileSortConfig(sandp.sandpile_struct, conf, sandp.perm_group, sort = False, verts = sandp.vertices)
        q_exp = sortedconfig.level()
        t_exp = new_delay(sortedconfig, n, k)
        poly = poly + (q**q_exp) * (t**t_exp)
    return poly

def is_symmetric(poly):         # Check if polynomial is symmetric in q and t
    R = FractionField(QQ['q','t','x'])  # type: ignore
    q,t,x = R.gens()
    poly1 = R(poly)     # Casting in new ring
    poly2 = poly1
    poly2(q=x)
    poly2(t=q)
    poly2(x=t)
    return (poly == poly2)

n = 6
k = 3

S = Multi_CliqueIndependent_SortedSandpile([n],[],k)        # Define clique with n vertices

poly = qt_poly(n,k,S)
print(poly)
print("The polynomial is symmetric? {}".format(is_symmetric(poly)))