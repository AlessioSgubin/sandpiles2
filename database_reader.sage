#   This code reads informations in the database of computed sorted recurrents
#   
#   We are following this naming scheme:
#       -   CHC_XX_YY    : a chain of cliques, sorted sandpile with XX = a and YY = n
#       -   CHI_XX_YY    : a chain of independents, sorted sandpile with XX = a and YY = n
#       -   LPC_XX_YY    : a loop of cliques, sorted sandpile with XX = a and YY = n
#       -   LPI_XX_YY    : a loop of independents, sorted sandpile with XX = a and YY = n

import pickle
import time
import os


dir_path = 'database'               # Directory/folder path
sandp_files = []                    # List to store files

for file_path in os.listdir(dir_path):                      # Iterate directory
    if os.path.isfile(os.path.join(dir_path, file_path)):   # Check if current file_path is a file
        sandp_files.append(file_path)

sandp_files.sort()                  # Sort names in the list

for name in sandp_files:
    namefile = str(dir_path + "/" + name)
    S = SortedSandpile.load(namefile)

    cell_list = S.specific_opt[1]
    card_cell = S.specific_opt[2]
    dict_cell = {}                  # Dictionary "cell -> list vertices"
    index_next_cell = 1
    for cell in cell_list:
        dict_cell = dict_cell | {cell:[index_next_cell + j for j in range(abs(card_cell[cell]))]}
        index_next_cell += abs(card_cell[cell])

    order = []                      # Reading order for qt-Polynomials
    for cell in cell_list:                                    
        order = order + dict_cell[cell]

    qt_poly = S.qt_Polynomial(ordered = order)
    q_poly = qt_poly(t = 1)
    t_poly = qt_poly(q = 1)
    num = qt_poly(q = 1, t = 1)
    print("Sorted Sandpile: {}".format(name))
    print("qt-poly: \t\t\t{}".format(qt_poly))
    print(" q-poly: \t\t\t{}".format(q_poly))
    print(" t-poly: \t\t\t{}".format(t_poly))
    print(" number: \t\t\t{}\n".format(num))

