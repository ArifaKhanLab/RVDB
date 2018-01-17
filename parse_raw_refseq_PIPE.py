import sys
homedir=sys.argv[1]
scriptdir=homedir+'\\'+'UPDATE_SCRIPTS_LOGS'
datetag=sys.argv[2]
currentvs=sys.argv[3]
wdir=homedir+'\\RVDBv'+currentvs
refseqdir=wdir+'\\RefSeq_raw_data_'+datetag
logdir=refseqdir+'\\'+'log'
refseq_viral_zipfilenames=sys.argv[4:]
refseq_viral_zipfilenames=[refseqdir+'\\'+refseq_viral_zipfilename for refseq_viral_zipfilename in refseq_viral_zipfilenames]
logdir=refseqdir+'\\'+'log'

####################################################################################################
######## Simple function to reformat a raw NCBI Viral fast file to have 2 lines /entry      ########
######## and output two files, phage and eukaryotic viral + archaea                         ########
####################################################################################################
def parse_raw_refseqviral(refseq_viral_unzipfilename):
##    print 'parsing '+raw_ncbi_viral_filename+' into phage and eukaryotic viral + archael sequences'
    outf1_name=refseqdir+'\\'+'viral.genomic.eukviral.fasta'
    outf2_name=refseqdir+'\\'+'viral.genomic.phage.fasta'
    logf.write('Refseq viral unzipfilename for eukaryotic entries: '+outf1_name+'\n')
    logf.write('Refseq viral unzipfilename for phage entries: '+outf2_name+'\n')
    outf1=open(outf1_name,'w')
    outf2=open(outf2_name,'w')
    phage_names_fn=logdir+'\\'+'phage_kws.txt'
    phage_names_inf=open(phage_names_fn)
    phage_names=phage_names_inf.read().strip().split('\n')
    phage_names_inf.close()
    p=0
    e=0
    eukviralaccs=[]
    phageaccs=[]
    inf=open(refseq_viral_unzipfilename)
    phage=False
    seq=''
    for i,line in enumerate(inf):
        if i==0:
            header=line.strip()
            acc=header.split('|')[2]
            for phage_name in phage_names:
                if phage_name in header.lower():
                    phage=True
                    break
            continue
        if line.startswith('>'):
            if phage:
                outf2.write(header+'\n'+seq+'\n')
                p+=1
                phageaccs.append(acc)
                phage=False
            else:
                outf1.write(header+'\n'+seq+'\n')
                e+=1
                eukviralaccs.append(acc)
                phage=False
            header=line.strip()
            acc=header.split('|')[2]
            seq=''
            for phage_name in phage_names:
                if phage_name in header.lower():
                    phage=True
                    break
        else:
            seq+=line.strip()
    if phage:
        outf2.write(header+'\n'+seq+'\n')
        p+=1
    else:
        outf1.write(header+'\n'+seq+'\n')
        e+=1
    inf.close()
##    print "read and reformatted "+str(e)+" eukaryotic viral + archaeal viral sequences"
##    print "read and reformatted "+str(p)+" phage sequences"
    outf1.close()
    outf2.close()
    outf3=open(refseqdir+'\\'+'viral.genomic.eukviral.accs.txt','w')
    outf4=open(refseqdir+'\\'+'viral.genomic.phage.accs.txt','w')
    outf3.write('\n'.join(eukviralaccs))
    outf4.write('\n'.join(phageaccs))
    logf.write(str(len(eukviralaccs))+' eukaryotic viral entries and '+str(len(phageaccs))+' phage entries'+'\n')
    outf3.close()
    outf4.close()


import gzip


logf=open(logdir+'\\'+'unzip_log_refseqviral.txt','w')
logf.write('RefSeq viral zipfilenames: '+','.join(refseq_viral_zipfilenames)+'\n')
refseq_viral_unzipfilename=refseqdir+'\\'+'viral.genomic.fna'
unzipf=open(refseq_viral_unzipfilename,'w')
logf.write('RefSeq viral unzipfilename: '+refseq_viral_unzipfilename+'\n')

for fn in refseq_viral_zipfilenames:
    fh = gzip.open(fn, 'rb')
    seq=''
    for i,line in enumerate(fh):
        if line.startswith('>'):
            if not i==0:
                unzipf.write(header+'\n'+seq+'\n')
                seq=''
            sl=line.strip().split()
            acc=sl[0].split('>')[1]
            desc=' '.join(sl[1:])
            header='>acc|REFSEQ|'+acc+'|'+desc
        else:
            seq+=line.strip()
    fh.close()
    unzipf.write(header+'\n'+seq+'\n')
unzipf.close()

parse_raw_refseqviral(refseq_viral_unzipfilename)
logf.close()
