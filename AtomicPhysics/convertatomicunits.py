
def convert(num, unit):
    return num*unit

def cmtoev(num):
    return convert(num, 1.2398*10**-4)

def evtocm(num):
    return convert(num, 8065.817)

def cmtoryd(num):
    return convert(num, 9.1122*10**-6)

def evtoryd(num):
    return convert(num, 1.0/13.606)

def rydtoev(num):
    return convert(num, 13.606)

print(evtocm(0.5))
