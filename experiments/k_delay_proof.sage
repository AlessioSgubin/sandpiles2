#
#       A program to test Parking Functions and the bijection with k-Delay
#
import copy
import sys

TRY1 = False
TRY2 = False
TRY3 = False
#TRY1 = True
#TRY2 = True
TRY3 = True

def kDyckPaths(n,k):                        # kDYCKPATHS: computes the set of all (nk,n)-Dyck paths
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

def kParkingFunctions(n,k):                 # kPARKINGFUNCTIONS: computes the set of all (n,k)-parking functions
    parking_functions = []
    counter = 0
    # Compute total number of Dyck paths
    totalnum = binomial(n*(k+1)+1, n)/(n*(k+1)+1)
    print("Starting the algorithm to compute all parking functions.")
    ### Loop on k-Dyck paths
    for path in kDyckPaths(n,k):
        counter += 1
        perc = floor(100*counter/totalnum)
        sys.stdout.write('\r')                          # Reset to start of line
        sys.stdout.write("Percentage %3d %%, computing parking functions for path %6d of %6d" % (perc, counter, totalnum))
        sys.stdout.flush()
        #print("Computing {}".format(path))
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
    print("\nCompleted computation!\n")
    return parking_functions

def nkParkingFunctions(n,k):
    parking_functions = []
    counter = 0
    # Compute total number of Dyck paths
    totalnum = binomial(n*(k+1)+1, n)/(n*(k+1)+1)
    print("Starting the algorithm to compute all parking functions.")
    # Compute all possible set partitions
    set_partitions = {}
    for lamb in Partitions(n):
        temp = []
        for possible in OrderedSetPartitions(n,lamb):
            poss = [list(part) for part in possible]
            temp.append(poss)
        lambt = tuple(lamb)
        set_partitions.update({lambt:temp})

    ### Loop on k-Dyck paths
    for path in kDyckPaths(n,k):
        counter += 1
        perc = floor(100*counter/totalnum)
        sys.stdout.write('\r')                          # Reset to start of line
        sys.stdout.write("Percentage %3d %%, computing parking functions for path %6d of %6d" % (perc, counter, totalnum))
        sys.stdout.flush()
        #print("Computing {}".format(path))

        # Compute the partitions we want
        rec_partition = [len([1 for w in path if w==n*k-i]) for i in range(n*k) if (n*k-i) in path]
        #print(rec_partition)
        ord_partition = copy.copy(rec_partition)
        ord_partition.sort()
        ord_partition.reverse()
        #print(ord_partition)
        ord_partitiont = tuple(ord_partition)
        #print(ord_partitiont)
        part_decr_ord = sorted(range(len(rec_partition)), key=lambda i: -n*rec_partition[i])
        read_order = sorted(range(len(part_decr_ord)), key=lambda i: part_decr_ord[i])
        for set_part in set_partitions[ord_partitiont]:
            label = []
            for j in read_order:
                temp2 = list(list(set_part)[j])
                temp2.sort()
                label = label + temp2
            parking_functions.append(copy.copy([path,label]))

    print("\nCompleted computation!\n")
    return parking_functions


def area(n,k,w):                            # AREA: computes the area of a given path
    return sum([w[i]-n*k+i*k for i in range(n)])

def pdinv(n,k,w):                           # PDINV: computes the path dinv of a given path
    dinv = 0
    # Read all cells above path
    for i in range(n):                  # Loop on rows
        for j in range(n*k-w[i]):       # Loop on columns
            arm = (n*k - w[i]) - (j+1)
            leg = i - max([h for h in range(n) if n*k-w[h]-1 < j]) - 1
            if (arm <= k*(leg+1)) and (k*leg < arm + 1):
                dinv += 1
    return dinv

def tdinv(n,k,w_path,w_label):              # TDINV
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

def maxtdinv(n,k,w_path):                   # MAXTDINV
    rank_list = [(n*k+1)*i + n*w_path[i] for i in range(n)]
    order = sorted(range(len(rank_list)), key=lambda k: rank_list[k])
    park = sorted(range(len(order)), key=lambda k: order[k])
    park = [i+1 for i in park]
    return tdinv(n,k,w_path,park)

def dinv(n,k,w_path,w_label):               # DINV: computes the dinv of a parking function
    #print("La statistica pdinv è {}".format(pdinv(n,k,w_path)))
    #print("La statistica maxtdinv è {}".format(maxtdinv(n,k,w_path)))
    #print("La statistica tdinv è {}".format(tdinv(n,k,w_path,w_label)))
    return pdinv(n,k,w_path) + tdinv(n,k,w_path,w_label) - maxtdinv(n,k,w_path)

def delay(n,k,w_path,w_label,contr=False):              # DELAY: computes the delay of a parking function
    delay = 0
    plus = 0
    current = n+1
    reading_word = []
    contributes = []
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
            contributes = contributes + [[current,plus]]
        buffer.remove(current)                 # Removing one copy of current from buffer

        # writing down the new letter in reading word
        reading_word = reading_word + [current]
    if contr:
        return [delay,contributes]
    else:
        return delay

def gen_poly(n,k,set_paths,ring):           # GEN_POLY: polynomial generated by paths
    q = ring.gens()[0]
    t = ring.gens()[1]
    poly = 0*q*t
    for path in set_paths:
        poly = poly + (q**pdinv(n,k,path))*(t**area(n,k,path))
    return poly

def path_word(n,k,w_path,w_label):
    rank_list = [(n*k+1)*i + n*w_path[i] for i in range(n)]
    path_order = sorted(range(len(rank_list)), key=lambda k: rank_list[k])      # Reading order of the columns
    path_word = [w_label[path_order[i]] for i in range(n)]
    return path_word

def Psi(n,k,w_path,w_label, infos = False):                # PSI: this is the bijection
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
    '''
    if not all(a <= b for a, b in zip(not_dinvcorr, not_dinvcorr[1:])):
        print("NO\nNO\nNO")
    '''

    ### Compute the not_dinv
    not_dinv = [not_tdinv[i] + not_dinvcorr[i] for i in range(n)]
    #print("The not_dinv: {}".format(not_dinv))

    ### Compute the w_path of the image
    not_dinv_sort = sorted(not_dinv, key=lambda k: k)
    area_word = [n*k-not_dinv_sort[i] for i in range(n)]
    
    ### Compute the w_label of the image
    label_word = [path_word[i] for i in sorted(range(n), key=lambda k: not_dinv[k]*n + path_word[k])]
    if not infos:
        return [area_word, label_word]
    else:
        return [area_word, label_word, not_tdinv, not_dinvcorr]

def sweep_map(n,k,w_path):
    # Write down the sequence of north (N)/east (E) steps
    steps = ['N']
    for i in range(n-1):
        steps = steps + ['E' for j in range(w_path[i]-w_path[i+1])] + ['N']
    steps = steps + ['E' for j in range(w_path[n-1]+1)]
    #print(steps)
    # Write the rank of each "wand"
    rank_steps = [0]
    current = 0
    for i in range(n*(k+1)):
        if steps[i] == 'N':
            current += n*k+1
        else:
            current -= n
        rank_steps = rank_steps + [current]
    # Write down sequence of steps reordered by rank
    new_order = sorted(range(len(rank_steps)), key=lambda i: rank_steps[i])
    #print(new_order)
    sweep_steps = [steps[i] for i in new_order]
    #print(rank_steps)
    #print(sweep_steps)
    # Compute new w_path
    w_sweep_path = []
    current = n*k
    for w in sweep_steps:
        if w == 'N':
            w_sweep_path = w_sweep_path + [current]
        else:
            current -= 1
    return w_sweep_path

def area_word(n,k,w_path):
    return [w_path[i] - (n-i)*k for i in range(n)]

def tdinv_word(n,k,w_path,w_label,path_word = []):
    ### Create the path word, ordered by rank
    rank_list = [(n*k+1)*i + n*w_path[i] for i in range(n)]
    path_order = sorted(range(len(rank_list)), key=lambda k: rank_list[k])      # Reading order of the columns
    if path_word == []:
        path_word = [w_label[path_order[i]] for i in range(n)]                      # Labels in reading order
    #print("The path_word: {}".format(path_word))

    ### Compute the not_tdinv
    # Compute all pairs of cars making tdinv
    tdinv_pairs = []            
    row_cresc = sorted(range(n), key=lambda i: w_label[i])
    for i in range(n):
        rank_c = n*w_path[row_cresc[i]] + (n*k+1)*row_cresc[i]
        for j in range(n-i-1):
            rank_d = n*w_path[row_cresc[i+j+1]] + (n*k+1)*row_cresc[i+j+1]
            if rank_c < rank_d and rank_d < rank_c + (n*k+1):
                tdinv_pairs = tdinv_pairs + [(i+1,i+j+2)]
    # Compute not_tdinv with previous cars
    tdinv = []
    labels = path_word
    for i in range(n):
        temp = 0
        for j in range(i):
            if ((labels[i],labels[j]) in tdinv_pairs) or ((labels[j],labels[i]) in tdinv_pairs):
                temp += 1
        tdinv = tdinv + [temp]
    return tdinv

def max_tdinv_word(n,k,w_path):
    # Compute w_label for maximum tdinv
    rank_list = [(n*k+1)*i + n*w_path[i] for i in range(n)]
    order = sorted(range(len(rank_list)), key=lambda k: rank_list[k])
    park = sorted(range(len(order)), key=lambda k: order[k])
    w_label_max = [i+1 for i in park]
    tdinv = tdinv_word(n,k,w_path,w_label_max)
    return tdinv

def dinvcorr_word(n,k,w_path,w_label,path_word = []):
    ### Create the path word, ordered by rank
    if path_word == []:
        rank_list = [(n*k+1)*i + n*w_path[i] for i in range(n)]
        path_order = sorted(range(len(rank_list)), key=lambda k: rank_list[k])      # Reading order of the columns
        path_word = [w_label[path_order[i]] for i in range(n)]                      # Labels in reading order
    #print("The path_word: {}".format(path_word))
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
    # Compute dinvcorr with previous cars
    dinvcorr = []
    for i in range(n):
        temp = 0
        for j in range(i):
            temp += len([1 for (a,b) in dinvcorr_pairs if (a,b)==(path_word[i],path_word[j]) or (a,b)==(path_word[j],path_word[i])])
        dinvcorr = dinvcorr + [temp]
    return dinvcorr

def not_dinvcorr_word(n,k,w_path,w_label,path_word = []):
    ### Create the path word, ordered by rank
    if path_word == []:
        rank_list = [(n*k+1)*i + n*w_path[i] for i in range(n)]
        path_order = sorted(range(len(rank_list)), key=lambda k: rank_list[k])      # Reading order of the columns
        path_word = [w_label[path_order[i]] for i in range(n)]                      # Labels in reading order
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
    return not_dinvcorr

def tries_w_reading(n,k,w_path,w_label):
    ### For now just try to get the reading word...
    [d,contrib] = delay(n,k,w_path,w_label,contr=True)
    rank_list = [(n*k+1)*i + n*w_path[i] for i in range(n)]
    rank_order = [w_label[i] for i in sorted(range(n), key=lambda j: rank_list[j])]
    #order = sorted(range(n), key=lambda i: n*contrib[i][1] + rank_order.index(contrib[i][0]))
    order = sorted(range(n), key=lambda i: n*k*contrib[i][1] + w_path[w_label.index(contrib[i][0])])
    w_reading = [contrib[i][0] for i in order]
    return w_reading

def w_reading(n,k,w_path,w_label):
    rank_list = [(n*k+1)*i + n*w_path[i] for i in range(n)]
    path_order = sorted(range(len(rank_list)), key=lambda k: rank_list[k])
    path_word = [w_label[path_order[i]] for i in range(n)]
    return path_word

def Phi_bruteforce(n,k,w_path,w_label):
    # Compute the not_dinv associated to labels
    order = sorted(range(n), key=lambda i: w_label[i])
    not_dinv_labels = [n*k - w_path[i] for i in order]

    # Get the delay word, that will correspond to the area
    val_delay,w_area = delay(n,k,w_path,w_label,contr=True)

    # Gather labels by area...
    labels_by_area = []
    for i in range(n*(k-1)+1):
        labels_by_area = labels_by_area + [ [w_area[j][0] for j in range(n) if w_area[j][1] == i] ]
    
    # List areas by label
    order = sorted(range(n),key=lambda i: w_area[i][0])
    area_by_label = [w_area[i][1] for i in order]
    
    # Try possibilities...
    partial_readings = [[]]
    for ar in range(n*(k-1)+1):
        adding = labels_by_area[ar]
        if adding != []:            # Try possible configurations
            new_possibilities = []
            for pos in Subsets(range(len(partial_readings[0])+len(adding)),len(adding)):
                positions = pos.list()
                for adding_order in Permutations(adding).list():
                    try_readings = copy.copy(partial_readings)
                    for temp in try_readings:
                        try_read = copy.copy(temp)
                        for j in range(len(adding)):            # Shuffle the labels
                            try_read.insert(positions[j],adding_order[j])
                        # Test if it is right
                        n2 = len(try_read)
                        path2 = [n2*(k-i) + area_by_label[try_read[i]-1] for i in range(n2)]
                        not_tdinv2 = tdinv_word(n2,k,path2,try_read)
                        not_dinvcorr2 = dinvcorr_word(n2,k,path2)
                        not_dinv2 = [not_tdinv2[i]+not_dinvcorr2[i] for i in range(n2)]
                        not_dinv = [not_dinv_labels[i-1] for i in try_read]
                        if not_dinv == not_dinv2:
                            new_possibilities = new_possibilities + [try_read]
            if len(new_possibilities) >= 1:
                partial_readings = new_possibilities
            else:
                print("There are {} possibilities:".format(len(new_possibilities)))
                print("Partial reading: {}".format(partial_readings))
                print("Possibilities {}".format(new_possibilities))
                raise Exception("SBAGLIATO!")
    return partial_readings


### Now let's test the bijection...
n = 5
k = 3


if TRY1:
    for [path1,label1] in nkParkingFunctions(n,k):
        area1 = area(n,k,path1)
        dinv1 = dinv(n,k,path1,label1)

        path2,label2,not_tdinv1,not_dinvcorr1 = Psi(n,k,path1,label1, infos=True)
        area2 = area(n,k,path2)

        area1_word = [- k*(n-i) + path1[i] for i in range(n)]
        [delay2,contrib] = delay(n,k,path2,label2,contr=True)
        delay_ord = sorted(range(n), key=lambda l: label1.index(contrib[l][0]))
        delay_word = [contrib[i][1] for i in delay_ord]
        check = (area1_word == delay_word)

        w_read1 = w_reading(n,k,path1,label1)
        w_read2 = tries_w_reading(n,k,path2,label2)
        check2 = (w_read1 == w_read2)

        path3 = sweep_map(n,k,path1)
        difference = [path3[i]-path2[i] for i in range(n)]
        check3 = (sum(difference) == maxtdinv(n,k,path1) - tdinv(n,k,path1,label1))

        difference2 = [not_dinvcorr1[i]-not_tdinv1[i] for i in range(n)]
        check4 = (difference == difference2)

        maxtdinv1 = max_tdinv_word(n,k,path1)
        tdinv1 = tdinv_word(n,k,path1,label1)
        difference3 = [- tdinv1[i] + maxtdinv1[i] for i in range(n)]
        check5 = (difference == difference3)

        '''
        if (area1 == delay2) and (dinv1 == area2):
            print("Funziona: \t[{},{}] \t---> \t[{},{}]".format(path1,label1,path2,label2))
            print("|\tdinv1 = area2\t{}".format(dinv1))
            print("|\tarea1 = delay2\t{}".format(area1))
            print("|\tw_area1 =?= w_delay2: {}\t {}".format(check,area1_word))
            print("|\tRight w_reading: {}\t {}\t{}".format(check2,w_read1,w_read2))
            print("")
        else:
            print("\nNO")
            print("Path1 {} con Label1 {} ha dinv {} e area {}".format(path1,label1,dinv1,area1))
            print("Path2 {} con Label2 {} ha area {} e delay {}\n".format(path2,label2,area2,delay2))
        '''
        '''
        if not check5:
            print("Funziona: \t[{},{}] \t---> \t[{},{}]".format(path1,label1,path2,label2))
            #print("|\tdinv1 = area2\t{}".format(dinv1))
            #print("|\tarea1 = delay2\t{}".format(area1))
            #print("|\tw_area1 =?= w_delay2: {}\t {}".format(check,area1_word))
            #print("|\tRight w_reading: {}\t {}\t{}".format(check2,w_read1,w_read2))
            #print("|\tPaths {}\t {}\t Difference {}".format(path2,path3,[path3[i]-path2[i] for i in range(n)]))
            print("|\tDifferences\t{}\t {}".format(difference, difference3))
            print("|\t\t\t{}\t {} with {}\t{}".format(tdinv1,maxtdinv1,tdinv(n,k,path1,label1),maxtdinv(n,k,path1)))
            print("")
        '''
        print("Funziona: \t[{},{}] \t---> \t[{},{}]".format(path1,label1,path2,label2))
        
        reading_new = Phi_bruteforce(n,k,path2,label2)
        reading_p1 = w_reading(n,k,path1,label1)
        
        if reading_new != reading_p1:
            print("NO! {}\t{}".format(reading_p1,reading_new))
            print(reading_new)
        else:
            print("YES!")

if TRY2:
    nkPF = nkParkingFunctions(n,k)
    nkPF2 = copy.copy(nkPF)
    print("There are {} parking functions".format(len(nkPF)))
    for [path1,label1] in nkPF:
        print("Doing {},{}".format(path1,label1))
        path2,label2 = Psi(n,k,path1,label1)

        if [path2,label2] in nkPF2:
            nkPF2.remove([path2,label2])
        else:
            print("Not bijective...")
    if nkPF2 == []:
        print("The map is bijective!")
    else:
        print("The map is NOT bijective!")

if TRY3:
    nkPF = nkParkingFunctions(n,k)
    for [path1,label1] in nkPF:
        # Compute area word
        area1 = area_word(n,k,path1)
        # Compute path word
        path_word1 = path_word(n,k,path1,label1)
        # Compute dinvcorr word (order: path word)
        not_dinvcorr_temp = not_dinvcorr_word(n,k,path1,label1)
        ordering = [path_word1.index(label1[i]) for i in range(n)]
        not_dinvcorr1 = [not_dinvcorr_temp[i] for i in ordering]

        diff = [i + not_dinvcorr1[i] - area1[i] for i in range(n)]
        print("Difference word of parking function {} \tis {}".format([path1,label1],diff))
