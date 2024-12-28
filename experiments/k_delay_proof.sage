#
#       A program to test Parking Functions and the bijection with k-Delay
#
TRY = True

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

def kParkingFunctions(n,k):             # kPARKINGFUNCTIONS: computes the set of all (n,k)-parking functions
    parking_functions = []
    ### Loop on k-Dyck paths
    for path in kDyckPaths(n,k):
        print(path)
        possible_labels = [[[],set([i+1 for i in range(n)])]]
        for row in range(n*k):
            #print(possible_labels)
            new_possible_labels = []
            for [label,remaining] in possible_labels:
                s = len([j for j in path if j == n*k - row])
                for sub in Subsets(remaining,s):
                    lab = sorted(sub)
                    new_possible_labels = new_possible_labels + [[(label + lab),remaining-set(sub)]]
            possible_labels = new_possible_labels
        parking_functions = parking_functions + [[path, label] for [label,rem] in possible_labels]
    return parking_functions


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

def tdinv(n,k,w_path,w_label):
    tdinv = 0
    row_cresc = [w_label.index(h+1) for h in range(n)]
    for i in range(n):
        rank_c = n*w_path[row_cresc[i]] + (n*k+1)*row_cresc[i]
        for j in range(n-i-1):
            rank_d = n*w_path[row_cresc[i+j+1]] + (n*k+1)*row_cresc[i+j+1]
            #print("Cell {} has rank {} and cell {} has rank {}".format(i+1,rank_c,i+j+2,rank_d))
            if rank_c < rank_d and rank_d < rank_c + (n*k):
                tdinv += 1
    return tdinv

def maxtdinv(n,k,w_path):
    rank_list = [(n*k+1)*i + n*w_path[i] for i in range(n)]
    order = sorted(range(len(rank_list)), key=lambda k: rank_list[k])
    park = sorted(range(len(order)), key=lambda k: order[k])
    park = [i+1 for i in park]
    return tdinv(n,k,w_path,park)

def dinv(n,k,w_path,w_label):           # DINV: computes the dinv of a parking function
    #print("La statistica pdinv è {}".format(pdinv(n,k,w_path)))
    #print("La statistica maxtdinv è {}".format(maxtdinv(n,k,w_path)))
    #print("La statistica tdinv è {}".format(tdinv(n,k,w_path,w_label)))
    return pdinv(n,k,w_path) + tdinv(n,k,w_path,w_label) - maxtdinv(n,k,w_path)

def delay(n,k,w_path,w_label):          # DELAY: computes the delay of a parking function
    delay = 0
    plus = 0
    current = n+1
    reading_word = []
    buffer = []
    for i in range(n*k):
        # Adding to buffer new labels with k multeplicity
        new = [w_label[j] for j in range(n) if w_path[j] == n*k-i]
        buffer = buffer + [new[floor(j/k)] for j in range(len(new)*k)]
        # Get next step
        if current <= min(buffer):          # Starting new ascending chain
            plus += 1
            current = max(buffer)
        else:                               # Continuing current ascending chain
            current = max([w for w in buffer if w < current])
        
        if current not in reading_word:     # Adding delay if needed...
            delay += plus
        buffer.remove(current)                 # Removing one copy of current from buffer

        # writing down the new letter in reading word
        reading_word = reading_word + [current]
    return delay

def gen_poly(n,k,set_paths,ring):
    q = ring.gens()[0]
    t = ring.gens()[1]
    poly = 0*q*t
    for path in set_paths:
        poly = poly + (q**pdinv(n,k,path))*(t**area(n,k,path))
    return poly
'''
def Phi(n,k,w_path,w_label):
    # This function tests the bijection between (dinv,area) -> (area,delay)

    # Create the path word, ordered by rank
    rank_list = [(n*k+1)*i + n*w_path[i] for i in range(n)]
    path_order = sorted(range(len(rank_list)), key=lambda k: rank_list[k])      # Reading order of the columns
    path_word = [w_label[path_order[i]] for i in range(n)]                      # Labels in order
    #print("The path_word: {}".format(path_word))

    ### Compute the not_tdinv
    # Compute all pairs of cars making tdinv
    tdinv_pairs = []            
    row_cresc = [w_label.index(h+1) for h in range(n)]
    for i in range(n):
        rank_c = n*w_path[row_cresc[i]] + (n*k+1)*row_cresc[i]
        for j in range(n-i-1):
            rank_d = n*w_path[row_cresc[i+j+1]] + (n*k+1)*row_cresc[i+j+1]
            if rank_c < rank_d and rank_d < rank_c + (n*k+1):
                tdinv_pairs = tdinv_pairs + [(i+1,i+j+2)]
    # Compute not_tdinv with previous cars
    not_tdinv = []
    for i in range(n):
        temp = 0
        for j in range(i):
            if ((path_word[i],path_word[j]) not in tdinv_pairs) and ((path_word[j],path_word[i]) not in tdinv_pairs):
                temp += 1
        not_tdinv = not_tdinv + [temp]
    #print("The not_tdinv: {}".format(not_tdinv))

    ### Compute the not_dinvcorr
    # Compute all pairs of cars making dinvcorr
    dinvcorr_pairs = []
    for i in range(n):
        rank_c = n*w_path[i] + (n*k+1)*i
        for j in range(n-i-1):
            rank_d = n*w_path[i+j+1] + (n*k+1)*(i+j+1)
            if rank_d < rank_c + n and rank_c - n*(k-1) < rank_d:
                dinvcorr_pairs = dinvcorr_pairs + [(w_label[i],w_label[i+j+1])]
    # Compute not_dinvcorr with previous cars
    not_dinvcorr = []
    for i in range(n):
        temp = 0
        for j in range(i):
            if ((path_word[i],path_word[j]) not in dinvcorr_pairs) and ((path_word[j],path_word[i]) not in dinvcorr_pairs):
                temp += 1
        not_dinvcorr = not_dinvcorr + [temp]
    #print("The not_dinvcorr: {}".format(not_dinvcorr))

    ### Compute the not_dinv
    not_dinv = [not_tdinv[i] + not_dinvcorr[i] for i in range(n)]
    #print("The not_dinv: {}".format(not_dinv))

    ### Compute the w_path of the image
    not_dinv_sort = sorted(not_dinv, key=lambda k: k)
    area_word = [n*k-not_dinv_sort[i] for i in range(n)]
    
    ### Compute the w_label of the image
    label_word = [path_word[i] for i in sorted(range(n), key=lambda k: not_dinv[k]*n + path_word[k])]
    return [area_word, label_word]
'''
def Psi(n,k,w_path,w_label):
    # This function tests the bijection between (dinv,area) -> (area,delay)

    ### Create the path word, ordered by rank
    rank_list = [(n*k+1)*i + n*w_path[i] for i in range(n)]
    path_order = sorted(range(len(rank_list)), key=lambda k: rank_list[k])      # Reading order of the columns
    path_word = [w_label[path_order[i]] for i in range(n)]                      # Labels in order
    #print("The path_word: {}".format(path_word))

    ### Compute the not_tdinv
    # Compute all pairs of cars making tdinv
    tdinv_pairs = []            
    row_cresc = [w_label.index(h+1) for h in range(n)]
    for i in range(n):
        rank_c = n*w_path[row_cresc[i]] + (n*k+1)*row_cresc[i]
        for j in range(n-i-1):
            rank_d = n*w_path[row_cresc[i+j+1]] + (n*k+1)*row_cresc[i+j+1]
            if rank_c < rank_d and rank_d < rank_c + (n*k+1):
                tdinv_pairs = tdinv_pairs + [(i+1,i+j+2)]

    # Compute not_tdinv with previous cars
    not_tdinv = []
    for i in range(n):
        temp = 0
        for j in range(i):
            if ((path_word[i],path_word[j]) not in tdinv_pairs) and ((path_word[j],path_word[i]) not in tdinv_pairs):
                temp += 1
        not_tdinv = not_tdinv + [temp]
    #print("The not_tdinv: {}".format(not_tdinv))

    ### Compute the not_dinvcorr
    # Compute all pairs of cars making dinvcorr
    dinvcorr_pairs = []
    for i in range(n):
        rank_c = n*w_path[i] + (n*k+1)*i
        for j in range(n-i-1):
            rank_d = n*w_path[i+j+1] + (n*k+1)*(i+j+1)
            for shift in range(k-1):
                if rank_d < rank_c + n*(k-1) - n*shift and rank_c - n - n*shift < rank_d:
                    dinvcorr_pairs = dinvcorr_pairs + [(w_label[i],w_label[i+j+1])]
    # Compute not_dinvcorr with previous cars
    not_dinvcorr = []
    for i in range(n):
        temp = 0
        for j in range(i):
            temp += k - 1 - len([1 for (a,b) in dinvcorr_pairs if (a,b)==(path_word[i],path_word[j]) or (a,b)==(path_word[j],path_word[i])])
        not_dinvcorr = not_dinvcorr + [temp]
    #print("The not_dinvcorr: {}".format(not_dinvcorr))

    ### Compute the not_dinv
    not_dinv = [not_tdinv[i] + not_dinvcorr[i] for i in range(n)]
    #print("The not_dinv: {}".format(not_dinv))

    ### Compute the w_path of the image
    not_dinv_sort = sorted(not_dinv, key=lambda k: k)
    area_word = [n*k-not_dinv_sort[i] for i in range(n)]
    
    ### Compute the w_label of the image
    label_word = [path_word[i] for i in sorted(range(n), key=lambda k: not_dinv[k]*n + path_word[k])]
    return [area_word, label_word]


if TRY:

    ### Now let's test the bijection...
    n = 5
    k = 2
    label1 = [i+1 for i in range(n)]
    for [path1,label1] in kParkingFunctions(n,k):
        area1 = area(n,k,path1)
        dinv1 = dinv(n,k,path1,label1)

        path2,label2 = Psi(n,k,path1,label1)
        area2 = area(n,k,path2)
        delay2 = delay(n,k,path2,label2)

        if (area1 == delay2) and (dinv1 == area2):
            print("Funziona: [{},{}] ---> [{},{}]".format(path1,label1,path2,label2))
        else:
            print("\nNO")
            print("Path1 {} con Label1 {} ha dinv {} e area {}".format(path1,label1,dinv1,area1))
            print("Path2 {} con Label2 {} ha area {} e delay {}\n".format(path2,label2,area2,delay2))