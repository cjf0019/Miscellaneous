file = open('TERMS')
header = file.readline()
newfile = open('TERMSfixed', 'a')
newfile.write(header)
for i in file:
   if len(i.split()) == 5: 
      new = i[:2]+' '+i[2:]
      newfile.write(new)
   else:
      newfile.write(i)
newfile.close()
