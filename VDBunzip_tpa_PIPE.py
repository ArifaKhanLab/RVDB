import gzip
import sys

wdir=sys.argv[1]
datetag=sys.argv[2]
tpafiletag=sys.argv[3]
import os
logf=open(wdir+'\\TPA_raw_data_'+datetag+'\\log\\'+'unzip_log_tpa.txt','w')
logf.write('raw data download directory: '+wdir+'\n')
logf.write('gzfiletag: '+gzfiletag+'\n')
logf.write('datetag: '+datetag+'\n')
logf.write('begin unzipped TPA filenames\n')
unzipfilenames=[]
for fn in os.listdir(wdir):
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
