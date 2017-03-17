#!/usr/bin/python2
import sys
inf = sys.argv[1]

mass = (16, 1, 1)
mol = [ [] for i in range(3) ]

def next_mol(inf):
    infp = open(inf,"r")
    mol_x = []
    atom_num = 3
    count = 0
    for i in infp:
        dat = i.split( )
        if len(dat) == 3:
            mol_x.append([float(x) for x in dat])
            count += 1
        if count == atom_num:
            yield mol_x
            mol_x = []
            count = 0
        

def torqueM(mol):
    m = [ 0.0 for i in range(6)]
    for i in range(3):
        m[0] += mol[i][1]**2 * mass[i] + mol[i][2]**2 * mass[i]
        m[1] += mol[i][0]**2 * mass[i] + mol[i][2]**2 * mass[i]
        m[2] += mol[i][0]**2 * mass[i] + mol[i][1]**2 * mass[i]
        m[3] -= mol[i][0] * mol[i][1] * mass[i]
        m[4] -= mol[i][0] * mol[i][2] * mass[i]
        m[5] -= mol[i][1] * mol[i][2] * mass[i]
    return m

for mol in next_mol(inf):
    tM = torqueM(mol) 
    tfm = '%5.2f ' * 5 + '%5.2f\n'
    print(tfm%tuple(tM))

