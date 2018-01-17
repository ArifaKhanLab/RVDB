import sys
import gzip
homedir=sys.argv[1]
date=sys.argv[2]
currentvs=sys.argv[3]
wdir=homedir+'\\'+'RVDBv'+currentvs
refseqdir=wdir+'\\RefSeq_raw_data_'+date
zipfilenames=sys.argv[4:-1]
zipfilenames=[refseqdir+'\\'+zipfilename for zipfilename in zipfilenames]
unzipfn=refseqdir+'\\'+sys.argv[-1]

unzipf=open(unzipfn,'w')
for zipfilename in zipfilenames:
    fh = gzip.open(zipfilename, 'rb')
    for i,line in enumerate(fh):
        unzipf.write(line)
    fh.close()
unzipf.close()

print unzipfn
