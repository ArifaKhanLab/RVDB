import sys
sys.path.append('E:\\TOOLBOX')
from sequence_record_functions_PIPE import get_accs_flatfile as getaccs
from sequence_record_functions_PIPE import get_accs_fastafile as getaccsfa
from sequence_record_functions_PIPE import get_filenames as getfns
homedir=sys.argv[1]
datetag=sys.argv[2]
currentvs=sys.argv[3]
previous_urvdb_filename=sys.argv[4]
##homedir='E:'
##datetag='apr.2018'
##currentvs='13.0'
wdir=homedir+'\\'+'RVDBv'+currentvs
refseqdir=wdir+'\\'+'RefSeq_raw_data_'+datetag
gbdir=wdir+'\\'+'GenBank_raw_data_'+datetag
tpadir=wdir+'\\'+'TPA_raw_data_'+datetag
gb_negkwdir=gbdir+'\\'+'negkw_out_'+datetag
tpa_negkwdir=tpadir+'\\'+'negkw_out_'+datetag
dupaccsfn=refseqdir+'\\'+'refseq_viral_originalaccs.txt'
dupaccs=getaccs(dupaccsfn)


def write_update_accs_outfile(postags,negtags,accstype):
    print "writing out + accessions, those that are RefSeq eukaryotic or pass SEM-R_PIPE screen"
    outf=open(wdir+'\\'+'RVDBv'+currentvs+'_accs'+accstype+'.txt','w')
    if accstype=='OK':
        refseq_accs=getaccs(refseqdir+'\\'+'viral.genomic.eukviral.accs.txt')
        outf.write('\n'.join(refseq_accs)+'\n')
    print "finished collecting RefSeq Viral accessions"
    gbfns=getfns(gb_negkwdir,postags,negtags)
    for gbfn in gbfns:
        print gbfn
        gbaccsfa=getaccsfa(gbfn)
        gbaccs=[]
        for gbacc in gbaccsfa:
            if not gbacc.split('.')[0] in dupaccs:
                gbaccs.append(gbacc)
        outf.write('\n'.join(list(gbaccs))+'\n')
    print "finished collecting GenBank+ accessions"
    tpafns=getfns(tpa_negkwdir,postags,negtags)
    for tpafn in tpafns:
        tpaaccs=getaccsfa(tpafn)
        outf.write('\n'.join(list(tpaaccs))+'\n')
    print "finished collecting TPA+ accessions"
    outf.close()

postags=['OK','VRL']
negtags=['FLAG','headers']
accstype='OK'
write_update_accs_outfile(postags,negtags,accstype)
postags=['AMB']
negtags=['FLAG','headers']
accstype='AMB'
write_update_accs_outfile(postags,negtags,accstype)

def make_review_sheets(previous_urvdb_filename):
    oldaccs=getaccsfa(previous_urvdb_filename)
    newaccs=getaccs(wdir+'\\'+'RVDBv'+currentvs+'_accsOK.txt')
    d1=set(oldaccs).difference(set(newaccs))
    d2=set(newaccs).difference(set(oldaccs))
    i1=set(oldaccs).intersection(set(newaccs))
    d1out=[]
    inf=open(previous_urvdb_filename)
    print "writing out headers for entries present in: "+previous_urvdb_filename+" but not in update v"+currentvs
    for line in inf:
        if line.startswith('>acc'):
            sl=line.strip().split('|')
            acc=sl[2]
            if acc in d1:
                d1out.append(sl[1:])
    inf.close()
    outf=open(wdir+'\\'+'RVDBv'+currentvs+'.missing.csv','wb')
    import csv
    writer=csv.writer(outf)
    d1out.insert(0,['SOURCE','ACCESSION','DESCRIPTION'])
    writer.writerows(d1out)
    outf.close()
    d2out=[]
    inf=open(refseqdir+'\\'+'viral.genomic.eukviral.fasta')
    match=False
    for line in inf:
        if line.startswith('>acc'):
            sl=line.strip().split('|')
            acc=sl[2]
            if acc in d2:
                d2out.append(sl[1:])
    inf.close()
    postags=['OK','VRL']
    negtags=['FLAG','AMB','headers']
    gbfns=getfns(gb_negkwdir,postags,negtags)
    tpafns=getfns(tpa_negkwdir,postags,negtags)
    readfns=[]
    readfns.extend(list(set(gbfns)))
    readfns.extend(list(set(tpafns)))
    print "writing out headers for entries present in update v"+currentvs+" that were not present in the previous version, "+previous_urvdb_filename
    for fn in readfns:
        print fn
        inf=open(fn)
        for line in inf:
            if line.startswith('>acc'):
                sl=line.strip().split('|')
                acc=sl[2]
                if acc in d2:
                    d2out.append(sl[1:])
        inf.close()
    outf=open(wdir+'\\'+'RVDBv'+currentvs+'.new.csv','wb')
    writer=csv.writer(outf)
    d2out.insert(0,['SOURCE','ACCESSION','DESCRIPTION'])
    writer.writerows(d2out)
    outf.close()

make_review_sheets(previous_urvdb_filename)
