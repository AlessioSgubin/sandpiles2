# This tests the bijection $\Sort_k(\mu,\nu)$ with $\Sort_k(\mu \sqcup \nu, \varnothing)$
load("./sorted_sandpile.py")

# Define the tilde map
def tilde_map(mu,nu,config):
    new_conf = {}       # Define new dictionary
    n = sum(mu) + sum(nu)
    i = 0
    for part_nu in nu:
        for j in range(part_nu):
            i += 1
            new_conf = new_conf | {i:config[i]+j}
    for j in range(n-i):
        new_conf = new_conf | {i+j+1:config[i+j+1]}
    return new_conf         # !!!!!! RETURNS A DICTIONARY !!!!!!

# Construct the two sandpiles
mu = [2]
nu = [4]
nu_rev = [nu[len(nu) - i - 1] for i in range(len(nu))]
k = 2
n = sum(mu) + sum(nu)

S1 = Multi_CliqueIndependent_SortedSandpile(mu, nu, k)
S2 = Multi_CliqueIndependent_SortedSandpile(mu + nu_rev, [], k)

# Compute recurrents of S1 and apply tilde map to configurations
print("Computing recurrents...")
dictrec1 = S1.sorted_recurrents()
sortrec1 = []
sortrec2 = []
for conf in dictrec1:
    new_dict = tilde_map(mu, nu, conf)
    old_conf = SandpileSortConfig(S1.sandpile_struct, conf, S1.perm_group, verts=S1.vertices)
    new_conf = SandpileSortConfig(S2.sandpile_struct, new_dict, S2.perm_group, verts=S2.vertices)
    sortrec1 = sortrec1 + [old_conf]
    sortrec2 = sortrec2 + [new_conf]

# Compute the k_delay for both cases
for i in range(len(sortrec1)):
    [delay1, wtopp1, record1] = sortrec1[i].k_delay_test()
    [delay2, wtopp2, record2] = sortrec2[i].k_delay_test()
    print("Configurations: \t{}\t->\t{}".format(sortrec1[i].sandpile_config, sortrec2[i].sandpile_config))
    print("Delay value:\t{}\t->\t{}".format(delay1, delay2))
    print("Toppl word: \t{}\n\t\t{}".format(wtopp1, wtopp2))
    print("Record word:\t{}\n\t\t{}\n".format(record1, record2))