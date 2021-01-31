
#avalues = {'5 6 2 50 3 6 2 1': 74900000, '11 6 3 7 9 6 2 1': 374000000, \
#        '3 6 1 7 5 6 0 3': 126000000}


def compareavals(nistaval):
    oic = open('oic')

    line = oic.readline()
    while "K   LV    T 2S+1    L   2J" not in line:
        line = oic.readline()
    
    uptrans = []
    lowtrans = []
    for trans in nistaval.keys():
        upper = ' '.join(trans.split()[:4])
        lower = ' '.join(trans.split()[-4:])
        uptrans.append(upper)
        lowtrans.append(lower)

    tofind = list(dict.fromkeys(uptrans+lowtrans))
    lvlnums = {}

    while len(tofind) != 0:
        line = oic.readline()
        linesplit = line[16:].split()
        qnums = ' '.join([linesplit[2],str(abs(int(linesplit[0]))),\
             linesplit[1],linesplit[3]])
        for trans in tofind:
            if qnums == trans:
                lvlnums[qnums] = line[5:10].strip()
                tofind.pop(tofind.index(qnums))

    line = oic.readline()
    while "CF   LV    W   CF   LV" not in line:
        line = oic.readline()

    findaval = {i:nistaval[i] for i in nistaval.keys()}

    def check(line, upper, lower, autoaval):
        split = line.split()
        for uplvl in upper:
            if uplvl.split()[3]+lvlnums[uplvl] or \
               lvlnums[uplvl] in split[:2]:
                lowlvl = lower[upper.index(uplvl)]
                if (lowlvl.split()[3]+lvlnums[lowlvl] in split[2:5]) or \
                    (lvlnums[lowlvl] in split[2:5]):
                    print("LOWLVL",lowlvl)
                    print("SPLIT",split[2:5])

                    autoaval[upper.pop(upper.index(uplvl)) + ' ' + \
                            lower.pop(lower.index(lowlvl))] = \
                            abs(float(split[-3]))
                    return 

    autoaval = {}
    while len(uptrans) != 0:
        line = oic.readline()
        check(line, uptrans, lowtrans, autoaval)
#        print("UPTRANS",uptrans)

    assert len(lowtrans) == 0
    oic.close()
    print("AUTOAVAL",autoaval)

    diff = 0
    for trans in nistaval.keys():
        diff += ((float(nistaval[trans])-float(autoaval[trans]))/ \
                float(nistaval[trans]))**2
    return diff

#print("AUTOAVALS",autoaval)
#print("DIFF",aval_diff(autoaval,avalues))
