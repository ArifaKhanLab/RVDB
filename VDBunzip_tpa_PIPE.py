import gzip
import sys

homedir=sys.argv[1]
datetag=sys.argv[2]
currentvs=sys.argv[3]
wdir=homedir+'\\'+'RVDBv'+currentvs
tpadir=wdir+'\\'+'TPA_raw_data_dec.2017'
print tpadir
tpafiletag=sys.argv[4]
import os
logf=open(tpadir+'\\log\\'+'unzip_log_tpa.txt','w')
logf.write('raw data download directory: '+wdir+'\n')
logf.write('tpafiletag: '+tpafiletag+'\n')
logf.write('datetag: '+datetag+'\n')
logf.write('begin unzipped TPA filenames\n')
unzipfilenames=[]
for fn in os.listdir(tpadir):
    if not tpafiletag in fn or not fn.endswith('.gz'):
        continue
    if not 'con_' in fn:
        continue
    print fn
    fh = gzip.open(tpadir+'\\'+fn, 'rb')
    unzipfn=fn.replace('.gz','.seq.'+datetag)
    unzipfilenames.append(unzipfn)
    unzipf=open(tpadir+'\\'+unzipfn,'w')
    seq=''
    c=0
    for i,line in enumerate(fh):
        if line.startswith('>'):
            if not i==0:
                unzipf.write(header+'\n'+seq+'\n')
                seq=''
            sl=line.strip().split()
            acc=sl[0].split('>')[1]
            desc=' '.join(sl[1:])
            header='>acc|TPA|'+acc+'|'+desc
        else:
            seq+=line.strip()
        if i-c==1000000:
            print i
            c=i
    unzipf.write(header+'\n'+seq+'\n')
    unzipf.close()

for unzipfn in unzipfilenames:
    logf.write(unzipfn)
    print unzipfn
