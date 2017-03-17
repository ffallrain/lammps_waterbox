#!/usr/bin/python
import sys,os

infile = sys.argv[1]

prefix = "./debug2pdb/"
outfile = prefix+"debug2pdb.data"

if not os.path.isdir(prefix):
    os.mkdir(prefix)

dbl_wat_str = '''REMARK 
ATOM      1  OW  SOL     1    %8.3f%8.3f%8.3f   0.0  0.00
ATOM      2  HW1 SOL     1    %8.3f%8.3f%8.3f   0.0  0.00
ATOM      3  HW2 SOL     1    %8.3f%8.3f%8.3f   0.0  0.00
ATOM      4  OW  SOL     2    %8.3f%8.3f%8.3f   0.0  0.00
ATOM      5  HW1 SOL     2    %8.3f%8.3f%8.3f   0.0  0.00
ATOM      6  HW2 SOL     2    %8.3f%8.3f%8.3f   0.0  0.00
END
'''

log_ofp = open(outfile,'w')
i = 0
for line in open(infile):
    items = line.split()
    if i % 8 == 0:
        file_index = int(items[0])
        temp = float(items[1])
        temp_body = float(items[2])
        potential = float(items[3])
        total = float(items[4])
    elif i % 8 == 1:
        oax = float(items[0])
        oay = float(items[1])
        oaz = float(items[2])
    elif i % 8 == 2:
        hax = float(items[0])
        hay = float(items[1])
        haz = float(items[2])
    elif i % 8 == 3:
        Hax = float(items[0])
        Hay = float(items[1])
        Haz = float(items[2])
    elif i % 8 == 4:
        obx = float(items[0])
        oby = float(items[1])
        obz = float(items[2])
    elif i % 8 == 5:
        hbx = float(items[0])
        hby = float(items[1])
        hbz = float(items[2])
    elif i % 8 == 6:
        Hbx = float(items[0])
        Hby = float(items[1])
        Hbz = float(items[2])
    elif i % 8 == 7:
        forcex = float(items[0])
        forcey = float(items[1])
        forcez = float(items[2])
        torquex = float(items[3])
        torquey = float(items[4])
        torquez = float(items[5])
        ofp = open( "%s/%d.pdb"%(prefix,file_index),'w' )
        ofp.write( dbl_wat_str%( oax,oay,oaz, hax,hay,haz, Hax,Hay,Haz, obx,oby,obz, hbx,hby,hbz, Hbx,Hby,Hbz ) )
        ofp.close()
        log_ofp.write("%d\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n"%(file_index,temp_body,potential,total,forcex,forcey,forcez,torquex,torquey,torquez) )
    i += 1
log_ofp.close()
print "Done."
    
