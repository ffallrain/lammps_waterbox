#!/usr/bin/python
import fpdb,gmx_top
import sys,os
import matplotlib.pyplot as plt

nonbondfile = '/opt/gromacs502/share/gromacs/top/amber99sb.ff/ffnonbonded.itp'
rtpfile = '/opt/gromacs502/share/gromacs/top/amber99sb.ff/aminoacids.rtp'
waterfile = '/opt/gromacs502/share/gromacs/top/amber99sb.ff/tip3p.itp'

infile = sys.argv[1]

def next_frame(infile):
    lines = list()
    for line in open(infile):
        if "MODEL" in line:
            lines = list()
            index = int(line.split()[1])
            lines.append(line)
        elif "ENDMDL" in line:
            lines.append(line)
            yield index,lines
        else:
            lines.append(line)

ofp = open('tip3p_potential.data','w')
gmxtop = gmx_top.gmxtop(nonbondfile,rtpfile,waterfile)
index_l = list()
vdw_l = list()
chg_l = list()
total_l = list()
for index,frame in next_frame(infile):
    pdb = fpdb.fPDB(frame)
    pdb.load_ff_params(gmxtop)
    water1 = pdb.topology.residues[0]
    water2 = pdb.topology.residues[1]
    vdw,chg = fpdb.potential_resi(water1,water2)
    ofp.write("%d\t%g\t%g\t%g\n"%(index,vdw,chg,vdw+chg))
    index_l.append(index)
    vdw_l.append(vdw)
    chg_l.append(chg)
    total_l.append(vdw+chg)
ofp.close()
    
data_sets = (vdw_l,chg_l,total_l)
labels = ('vDW','Charge','Total')
plt.figure(figsize = (14,4) )
plt.title("Potential in TIP3p")
for i in range(len(data_sets)):
    label = labels[i]
    data = data_sets[i]
    ax = plt.subplot( 1,len(data_sets),i+1 )
    ax.plot(index_l,data,label = label,color = 'g')
    ax.legend(loc='best')
    ax.set_xlabel("Index of frame")
    ax.set_ylabel("%s (kcal/mol)"%label)
plt.show()
    
    


