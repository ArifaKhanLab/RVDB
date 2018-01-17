import csv
import sys

def get_duplicate_acc(clustering_filename):
    inf=open(clustering_filename)
    dups=dict()
    mappable=False
    for i,line in enumerate(inf):
        if line.startswith('>Cluster'):
            continue
        if '*crep*' in line:
            if 'C_' in line:
                mappable=True
                rs_len=int(line.split()[1].split('nt')[0])
                rs_acc=line.split('|')[2]
                rs_desc=line.split('|')[3]
            else:
                mappable=False
        else:
            if mappable:
                s_len=int(line.split()[1].split('nt')[0])
                s_acc=line.split('|')[2]
                s_desc=line.split('|')[3]
                identity=float(line.split('%')[0].split('/')[1])
                if s_len==rs_len and s_desc==rs_desc and identity==float(100.00):
                    dups[rs_acc]=s_acc

def get_refseq_viral_neighbors_dict(rsn_filename):
    inf=open(rsn_filename)
    reader=csv.reader(inf)
    neighbors_dict=dict()
    for r,row in enumerate(reader):
        if r<=1:
            continue
        rsaccs=row[0].strip().split(',')
        nacc=row[1].strip().split('.')[0]
        for rsacc in rsaccs:
            try:
                existing=neighbors_dict[rsacc]
            except KeyError:
                existing=[]
            existing.append(nacc)
            neighbors_dict[rsacc]=existing
    del reader
    inf.close()
    return neighbors_dict

def get_gb_comments(comments_dict,gb_filename):
##    print "comments dict is of length: "+str(len(comments_dict))
    inf=open(gb_filename)
    records=inf.read().strip().split('//\n')
    comments=[]
    c=0
    for record in records:
        try:
            comment=record.split('COMMENT')[1].split('FEATURES')[0].strip()
        except IndexError:
            comment='NA'
        comments.append(comment)
    rs_accs=[]
    for record in records:
        try:
            rs_acc=record.split('ACCESSION')[1].split('\n')[0].strip()
        except IndexError:
            c+=1
            rs_acc='NA'+str(c)
        rs_accs.append(rs_acc)
    for r,rs_acc in enumerate(rs_accs):
        rs_acc_s=rs_acc.strip().split()
        for rs_acc_ss in rs_acc_s:
            comments_dict[rs_acc_ss]=comments[r]
    inf.close()
    return comments_dict

def make_duplicates_dict(neighbors_dict,comments_dict):
    dups=dict()
    for rs_acc in comments_dict.keys():
        comment=comments_dict[rs_acc]
        try:
            neighbor_accs=neighbors_dict[rs_acc]
        except KeyError:
            continue
        for neighbor_acc in neighbor_accs:
            if neighbor_acc in comment:
                dups[rs_acc]=neighbor_acc
    return dups

def extend_duplicates_dict(comments_dict):
    import re
    accfinder=re.compile('[A-Z]+[0-9]+[A-Z0-9]*')
    dups2=dict()
    for rs_acc in comments_dict.keys():
        comment=comments_dict[rs_acc]
        accs=accfinder.findall(comment)
        if len(accs)==0:
            continue
        acc=accs[0]
        dups2[rs_acc]=acc
    return dups2


homedir=sys.argv[1]
date=sys.argv[2]
currentvs=sys.argv[3]
wdir=homedir+'\\RVDBv'+currentvs
refseqdir=wdir+'\\RefSeq_raw_data_'+date
sys.path.append(homedir+'\\TOOLBOX')
from sequence_record_functions_PIPE import get_gis_flatfile as getaccs
from sequence_record_functions_PIPE import collate_dict_allvalues as colldict
logdir=refseqdir+'\\'+'log'
logf=open(logdir+'\\'+'mapping_refseq_genbank.txt','w')
import os
for fn in os.listdir(refseqdir):
    if fn.startswith('refseqviral_neighbors_mapping') and fn.endswith('.csv'):
        nmapfilename=refseqdir+'\\'+fn
        break
refseq_gbfilenames=[]
for fn in os.listdir(refseqdir):
    if fn.startswith('viral.genomic_'):
        refseq_gbfilenames.append(refseqdir+'\\'+fn)
neighbors_dict=get_refseq_viral_neighbors_dict(nmapfilename)
neighbors=colldict(neighbors_dict)
comments_dict=dict()
for refseq_gbfilename in refseq_gbfilenames:
    comments_dict=get_gb_comments(comments_dict,refseq_gbfilename)
dups=make_duplicates_dict(neighbors_dict,comments_dict)
dups2=extend_duplicates_dict(comments_dict)
dups3=dict()
import sys
rvdb_rs_accs=getaccs(refseqdir+'\\'+'viral.genomic.eukviral.accs.txt')
rvdb_rs_accs=[rs_acc.split('.')[0] for rs_acc in rvdb_rs_accs]
logf.write(str(len(rvdb_rs_accs))+' refseq viral accessions'+'\n')
logf.write(str(len(neighbors))+' neighbor accessions'+'\n')
logf.write(str(len(dups))+' refseq viral accessions mapped to original entries using neighbors annotation'+'\n')
logf.write(str(len(dups2))+' refseq viral accessions mapped to original entries using GenBank metadata (comments section)'+'\n')
logf.write(str(len(set(dups.keys()).intersection(set(dups2.keys()))))+' refseq viral accessions in common between the two above'+'\n')
for rs_acc in rvdb_rs_accs:
    try:
        dups3[rs_acc]=dups[rs_acc]
    except KeyError:
        try:
            dups3[rs_acc]=dups2[rs_acc]
        except KeyError:
            logf.write(rs_acc+' could not be mapped'+'\n')
            continue
logf.write(str(len(dups3))+' of '+str(len(rvdb_rs_accs))+' refseq viral accessions mapped to original entries in total'+'\n')
outfilename=refseqdir+'\\'+'refseq_viral_originalaccs.txt'
outf=open(outfilename,'w')
outf.write('\n'.join(sorted(dups3.values())))
outf.close()
logf.close()
outfilename=refseqdir+'\\'+'neighbor_accs.txt'
print outfilename
outf=open(outfilename,'w')
outf.write('\n'.join(neighbors))
outf.close()
        
            
        
