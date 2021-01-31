import sys
"""
Get the transition numbers for OMEGA extraction of ionization cross sections. NOTE: Assumes elastic collisions are not included. Needs modification if they are present (one extra transition per level).
"""

# 15960 for w0
# 369 for w0

def gettransnos(trmno, total, iptrm):
    ntr = 0
    ntrmn = 0
    for i in range(trmno-1):
        total -= 1
        iptrm -= 1
        ntr += total

    iptrm -= 1    
    ntrmn = ntr + iptrm
    total -= 1
    ntr += total
    return ntrmn, ntr

def writetranstofile(trmno,total,iptrm):
    file = open('multistgsig.txt', 'a')
    file.write(' '.join([str(i) for i in gettransnos(trmno,total,iptrm)])+" "+str(trmno)+" 1\n")
    file.close()
