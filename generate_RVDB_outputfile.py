import sys
homedir=sys.argv[1]
date=sys.argv[2]
currentvs=sys.argv[3]
wdir=homedir+'\\'+'RVDBv'+currentvs

def collect_SEMRqualifying_filenames(wdir,date,currentvs):
    qfilenames=[]
    gb_pass_dir=wdir+'\\'+'GenBank_raw_data'+date+'\\'+'negkw_out_'+date
    refseq_pass_dir=wdir+'\\'+'RefSeq_raw_data'+date
    tpa_pass_dir=wdir+'\\'+'TPA_raw_data'+date+'\\'+'negkw_out_'+date
    import os
    qfilenames.append(refseq_pass_dir+'\\'+'viral.genomic.eukviral.fasta')
    for fn in os.listdir(gb_pass_dir):
        if fn.endswith('OK.fasta') or fn.endswith('VRL.fasta'):
            qfilenames.append(gb_pass_dir+'\\'+fn)
    for fn in os.listdir(tpa_pass_dir):
        if fn.endswith('OK.fasta'):
            qfilenames.append(tpa_pass_dir+'\\'+fn)
    outf=open(wdir+'\\'+'U-RVDBv'+currentvs+'.fasta','w')
    for fn in qfilenames:
        inf=open(fn)
        for line in inf:
            outf.write(line.strip()+'\n')
        inf.close()
    outf.close()
