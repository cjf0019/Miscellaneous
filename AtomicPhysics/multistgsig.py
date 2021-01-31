import re
import subprocess
from gettransnums import *

termlist = open('termlist.txt')
for i in termlist:
    writetranstofile(int(i),15960,369)
#    writetranstofile(int(i),18960, 393)

inputfile = open('multistgsig.txt')
inputtext = inputfile.read()
newntr = [line.split() for line in inputtext.split("\n")]
template = open('bpytemplate').read()

for i in newntr:
    if i == []:
       continue
    print(i)
    bpy = open('bpy'+i[-2],'w')
    bpy.write(template+'\n' + 'python pystgsig.py ../omaddt '+' '.join(i)+"\n")
    bpy.close()
