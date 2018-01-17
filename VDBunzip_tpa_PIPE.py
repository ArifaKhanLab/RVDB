import gzip
import sys

homedir=sys.argv[1]
datetag=sys.argv[2]
currentvs=sys.argv[3]
wdir=homedir+'\\'+'RVDBv'+currentvs
tpadir=wdir+'\\'+'TPA_raw_data_dec.2017'
tpafiletag=sys.argv[4]
import os
logf=open(tpadir+'\\log\\'+'unzip_log_tpa.txt','w')
logf.write('raw data download directory: '+wdir+'\n')
logf.write('tpafiletag: '+tpafiletag+'\n')
logf.write('datetag: '+datetag+'\n')
logf.write('begin unzipped TPA filenames\n')
unzipfilenames=[]
for fn in os.listdir(tpadir):
    if tpafiletag in fn or not fn.endswith('.gz'):
        continue
    fh = gzip.open(wdir+'\\'+fn, 'rb')
    unzipfn=fn.replace('.gz','.seq.'+datetag)
    logf.write(unzipfn+'\n')
    unzipfilenames.append(unzipfn)
    unzipf=open(wdir+'\\'+unzipfn,'w')
    seq=''
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
    unzipf.write(header+'\n'+seq+'\n')
    unzipf.close()

for unzipfn in unzipfilenames:
    logf.write(unzipfn)
    print unzipfn
