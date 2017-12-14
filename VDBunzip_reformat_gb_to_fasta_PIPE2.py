import sys
sys.path.append('C:\\Python27\\Lib\\site-packages')
from Bio import SeqIO
import gzip


#t=[]
#for line in open(sys.argv[2], 'rb').readlines():
 #   t.append(line.strip())
#print t

wdir=sys.argv[1]
datetag=sys.argv[2]
gzfiletag=sys.argv[3]
gbdir=wdir+'\\'+'GenBank_raw_data_'+datetag
import os
logf=open(gbdir+'\\log\\'+'unzip_log.txt','w')
logf.write('raw data download directory: '+gbdir+'\n')
logf.write('gzfiletag: '+gzfiletag+'\n')
logf.write('datetag: '+datetag+'\n')
logf.write('begin unzipped filenames\n')
unzipfilenames=[]
neighboraccs_filename=sys.stdin.readlines()[0].strip()
##print neighboraccs_filename
inf=open(neighboraccs_filename)
neighbors=set(inf.read().strip().split('\n'))
inf.close()
for fn in os.listdir(gbdir):
    if not fn.startswith(gzfiletag) or not fn.endswith('.gz'):
        continue
    fh = gzip.open(gbdir+'\\'+fn, 'rb')
    unzipfn=fn.replace('.gz','.'+datetag)
    logf.write(unzipfn+'\n')
    unzipfilenames.append(unzipfn)
    unzipf=open(gbdir+'\\'+unzipfn,'w')
    for g,gb_record in enumerate(SeqIO.parse(fh,'genbank')):
        acc = gb_record.annotations['accessions'][0]
        organism = gb_record.annotations['organism']
        tax_line = ("; ").join(gb_record.annotations['taxonomy'])
        feat = gb_record.features
        disc = gb_record.description
        seqlen = len(gb_record)
        seq = gb_record.seq.tostring()
        acc = gb_record.id
        source = gb_record.annotations['source']
        key = gb_record.annotations['keywords']
        file_type = gb_record.annotations['data_file_division']
        date = gb_record.annotations['date']
        if acc.split('.')[0] in neighbors:
            db='NEIGHBOR'
        else:
            db='GENBANK'
        header='>acc|'+db+'|'+acc+'|'+disc+'|'+organism+'|'+file_type+'|'+ date
        unzipf.write(header+'\n'+seq+'\n')
    unzipf.close()

##logf.write('end unzipped filenames\n')
##logf.close()

##for unzipfn in unzipfilenames:
##    print unzipfn

       
