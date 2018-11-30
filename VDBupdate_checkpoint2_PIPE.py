import sys

homedir=sys.argv[1]
datetag=sys.argv[2]
currentvs=sys.argv[3]
relnotes=sys.argv[4]
wdir=homedir+'\\RVDBv'+currentvs
gbdir=wdir+'\\GenBank_raw_data_'+datetag
logdir=gbdir+'\\'+'log'
##### Put release file for update here #####
gbrelfilename=logdir+'\\'+relnotes

import os
import re
import time
import datetime

gb_divtypes=['ENV','HTC','INV','MAM','PLN','PRI','ROD','TSA','VRL','VRT']
basefilenamefinder=re.compile('gb[a-z]{3}[0-9]+.')
gbtypefinder=re.compile('gb[a-z]{3}')
fnnumfinder=re.compile('[0-9]+.')
gb_update_counts=dict()
##### This is the first of three output files. This prints the real-time progress of the unzipping, file per file,
##### and also prints out the total file count and seqs count per division,each time a file is unzipped
outf=open(logdir+'\\'+'RVDBv'+currentvs+'_checkpt2a.log','w')
unzipfns=[]

for f,fn in enumerate(os.listdir(gbdir)):
    if 'features' in fn:
        ###Only want to look at fasta file, not features file
        continue
    if not fn.endswith(datetag):
        continue
    unzipfns.append(fn)
    ts=time.time()
    st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
    outf.write('Initiating scan of GenBank file '+fn+' '+st+'\n')
    try:
        basefilename=basefilenamefinder.findall(fn)[0]
    except IndexError:
        ###This is for TPA sequences - no official release notes for them, so we have to count them separately
        continue
    gbtype=gbtypefinder.findall(basefilename)[0][2:]
    fnnum=fnnumfinder.findall(basefilename)[0][:-1]
    print gbtype+'\t'+str(fnnum)+'\t'+'('+basefilename[:-1]+')'
    try:
        existing_gbtype=gb_update_counts[gbtype]
    except KeyError:
        existing_gbtype=dict()
    existing_fnnum=0
    inf=open(gbdir+'\\'+fn)
    for i,line in enumerate(inf):
        if line.startswith('>acc'):
            existing_fnnum+=1
    existing_gbtype[fnnum]=existing_fnnum
    gb_update_counts[gbtype]=existing_gbtype
    inf.close()
    outf.write(str(existing_fnnum)+' '+'sequences scanned\n')
    outf.write('Files scanned so far: '+'\n')
    for gbtype in gb_update_counts.keys():
        fns=gb_update_counts[gbtype]
        fns=[int(ff) for ff in fns]
        fns=sorted(fns)
        fns=[str(ff) for ff in fns]
        fns=sorted(fns)
        outf.write('\t'+gbtype+' '+','.join(fns)+'\n')
    outf.write('Running total of sequences from files scanned: '+'\n')
    for gbtype in gb_update_counts.keys():
        total=0
        existing_gbtype=gb_update_counts[gbtype]
        for thisfn in existing_gbtype.keys():
            total+=int(existing_gbtype[thisfn])
        outf.write('\t'+gbtype+' '+str(total)+'\n')
    ts=time.time()
    st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
    outf.write('Completed scan of GenBank file '+fn+' '+st+'\n\n')


##################################################################################################
########## Beginning of execute control commands                                       ###########
##################################################################################################

##### This is the second log file. This file prints out a summary of the seqs unzipped
##### with just the file name tag and # seqs (division+fnnum    *count*, e.x. GBENV12     82356)
outf=open(logdir+'\\'+currentvs+'_checkpt2b.log','w')
cdict2=dict()
for div in sorted(gb_update_counts.keys()):
    bydiv=gb_update_counts[div]
    for fnnum in sorted(bydiv.keys()):
	c=bydiv[fnnum]
	fntag=div.upper()+fnnum
	outf.write(fntag+'\t'+str(c)+'\n')
	cdict2[fntag]=c
outf.close()

##### Collects counts per division / per file (checkpt2b.log info)
inf=open(logdir+'\\'+currentvs+'_checkpt2b.log')
cdict2=dict()
for i,line in enumerate(inf):
    fntag,c=line.strip().split('\t')
    cdict2[fntag]=c
inf.close()

##### Collects official release counts, at the same level (i.e. file  #seqs)
##### so as to compare to unzipped totals
cdict1=dict()
inf=open(gbrelfilename)
intext=inf.read().split('2.2.6')[3].split('Division     Entries    Bases')[1].split('2.2.7')[0].strip()
inf.close()
for line in intext.split('\n'):
    fntag,c,bases=line.strip().split()
    cdict1[fntag]=c

##### This is the third log file. This file prints out any discrepancy in seq counts
##### between official release notes (cdict1) and unzipped files (cdict2)
##### Note that, since the "for" loop reads in keys from cdict2, this procedure will not
##### verify whether all the files are present, only whether files present (i.e. unzipped)
##### contain the full, correct number of seqs
outf=open(logdir+'\\'+currentvs+'_checkpt2c.log','w')
for fntag in sorted(cdict2.keys()):
    c2=cdict2[fntag]
    try:
        c1=cdict1[fntag]
        if not c1==c2:
            outf.write(fntag+'\t'+str(c1)+'\t'+str(c2)+'\n')
    except KeyError:
        outf.write(fntag+'\t'+'dne'+'\t'+str(c2)+'\n')
outf.close()

##### This is the fourth log file. This file prints out a total amount of entries per division fo release notes vs. unzipped + counted in-house
bydiv1=dict()
bydiv2=dict()
import re
divfinder=re.compile('[A-Z]+')
for fntag in sorted(cdict1.keys()):
    div=divfinder.findall(fntag)[0]
    c=cdict1[fntag]
    c=int(c)
    try:
        existing=bydiv1[div]
    except KeyError:
        existing=0
    existing+=c
    bydiv1[div]=existing
for fntag in sorted(cdict2.keys()):
    div=divfinder.findall(fntag)[0]
    c=cdict2[fntag]
    c=int(c)
    try:
        existing=bydiv2[div]
    except KeyError:
        existing=0
    existing+=c
    bydiv2[div]=existing
outf=open(logdir+'\\'+currentvs+'_checkpt2d.log','w')
for div in sorted(bydiv1.keys()):
    c1=bydiv1[div]
    try:
        c2=bydiv2[div]
    except KeyError:
        continue
    outf.write(div+'\t'+str(c1)+'\t'+str(c2)+'\n')
outf.close()

for unzipfn in unzipfns:
    print unzipfn

