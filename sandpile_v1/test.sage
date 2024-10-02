### Definiamo un grafo

G = graphs.CompleteGraph(10)
H = graphs.CompleteMultipartiteGraph([4,5,6])

S = Sandpile2(G,0)

list_rec = S.recurrent_conf(check = True)
print(list_rec[1])

for i in list_rec[1]:
    S.current_conf = i
    if not S.is_recurrent():
        print("NOOOO")
    else:
        print("YESSS")

