import sys
sys.path.append('E:\\TOOLBOX')
from sequence_record_functions_PIPE import get_accs_flatfile as getaccs
from sequence_record_functions_PIPE import get_accs_fastafile as getaccsfa
from sequence_record_functions_PIPE import get_filenames as getfns
homedir=sys.argv[1]
datetag=sys.argv[2]
currentvs=sys.argv[3]
removeaccsfn=sys.argv[4]
wdir=homedir+'\\'+'RVDBv'+currentvs
gbdir=wdir+'\\'+'GenBank_raw_data_'+datetag
refseqdir=wdir+'\\'+'RefSeq_raw_data_'+datetag
tpadir=wdir+'\\'+'TPA_raw_data_'+datetag
gb_negkwdir=gbdir+'\\'+'negkw_out_'+datetag
tpa_negkwdir=tpadir+'\\'+'negkw_out_'+datetag
dupaccsfn=refseqdir+'\\'+'refseq_viral_originalaccs.txt'

postags=['OK','VRL']
negtags=['FLAG','AMB','headers']
gbfns=getfns(gb_negkwdir,postags,negtags)
tpafns=getfns(tpa_negkwdir,postags,negtags)
allfns=[]
allfns.append(refseqdir+'\\'+'viral.genomic.eukviral.fasta')
allfns.extend(gbfns)
allfns.extend(tpafns)
removeaccs=getaccs(wdir+'\\'+removeaccsfn)
dupaccs=getaccs(dupaccsfn)
outf=open(wdir+'\\'+'U-RVDBv'+currentvs+'.fasta','w')
c=0
match=False
written=set([])
for fn in allfns:
    inf=open(fn)
    for line in inf:
        if line.startswith('>acc'):
            acc=line.split('|')[2].split('.')[0]
            if acc in removeaccs or acc in dupaccs or acc in written:
                match=False
            else:
                match=True
                written.add(acc)
                c+=1
        if match:
            outf.write(line.strip()+'\n')
    inf.close()
outf.close()
