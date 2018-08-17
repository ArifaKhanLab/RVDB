import sys
sys.path.append('F:\\UPDATE_SCRIPTS_LOGS')

homedir=sys.argv[1]
datetag=sys.argv[2]
currentvs=sys.argv[3]
wdir=homedir+'\\'+'RVDBv'+currentvs
rvdbfilename=wdir+'\\'+sys.argv[4]
clstrfilename=wdir+'\\'+sys.argv[5]

inf=open(clstrfilename)
creps=set([])
for line in inf:
    if line.endswith('*crep*\n'):
        acc=line.split('|')[2]
        creps.add(acc)
inf.close()
inf=open(rvdbfilename)
crvdbfilename=rvdbfilename.replace('U-','C-')
outf=open(crvdbfilename,'w')
match=False
c=0
d=0
for i,line in enumerate(inf):
    if line.startswith('>acc'):
        acc=line.split('|')[2]
        if acc in creps:
            match=True
            c+=1
        else:
            if '|TPA|' in line:
                match=True
            else:
                match=False
    if match:
        outf.write(line)
        d+=1
print crvdbfilename+' created with '+str(c)+' entries and '+str(d)+' lines'
