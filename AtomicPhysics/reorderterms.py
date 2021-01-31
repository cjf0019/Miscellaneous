import numpy as np
file = open('omaddt')

file.readline()
file.readline()
line = file.readline().split()
terms = []

for i in range(15960):
    one = int(line.pop(0))
    two = int(line.pop(0))
    terms += [[abs(one), two]]

termsfile = open('TERMSreorderedtest')
firstline = termsfile.readline()

#j = termsfile.readline()
#newterms = [j]
#for i in termsfile:
#    if np.round(float(i.split()[-1]),5) == np.round(float(j.split()[-1]),5):
#        if 
writefile = open('TERMSreordered', 'w')
writefile.write(firstline)

#j = termsfile.readline()
#writefile.write(j)
reorder = 0

#for term in terms[:100]:
#    i = termsfile.readline()
#    isplit = i.split()
#    termsterms = [int(k) for k in isplit[:2]]
#    if termsterms == term:
#        writefile.write(i)
#    else:
#        print('FALSE')
#        if reorder == 0:
#            j = i
#            reorder = 1
#        else:
#            writefile.write(i)
#            writefile.write(j)
#            reorder = 0


termlines = []

for termline in termsfile:
    termlines.append(termline)
    

reordered = [] 

for i in range(len(terms)):
    match = 0
    termno = 0
    while match == 0:
        if [int(k) for k in termlines[termno].split()[:2]] == terms[i]:
            reordered.append(termlines.pop(termno))
            match = 1
        else:
            termno += 1

print(termlines)    

writefile.write(''.join(reordered))
writefile.write(termlines.pop(0))
writefile.close()
