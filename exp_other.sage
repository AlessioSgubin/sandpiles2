#
#               SOME OTHER EXPERIMENTS
#
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
