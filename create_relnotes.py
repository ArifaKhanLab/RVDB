# -*- coding: cp1252 -*-
import sys
from collections import defaultdict
homedir=sys.argv[1]
date=sys.argv[2]
currentvs=sys.argv[3]
wdir=homedir+'\\'+'RVDBv'+currentvs
relday=sys.argv[4]
relmon=sys.argv[5]
relyr=sys.argv[6]
gbmon=sys.argv[7]
gbmon=gbmon[0].upper()+gbmon[1:].lower()
gbyr=sys.argv[8]
rsmon=sys.argv[9]
rsmon=rsmon[0].upper()+rsmon[1:].lower()
rsyr=sys.argv[10]
rvdbfilename=wdir+'\\'+sys.argv[11]
try:
    relmon=int(relmon)
except ValueError:
    relmon=relmon[0].upper()+relmon[1:].lower()
    mondict=dict(zip(['Jan','Feb','March','April','May','June','July','August','Sept','Oct','Nov','Dec'],['1','2','3','4','5','6','7','8','9','10','11','12']))
    relmon=mondict[relmon]
##homedir='F:'
##currentvs='12.1'
##date='dec.2017'
##wdir=homedir+'\\'+'RVDBv'+currentvs
##relday='5'
##relmon='Feb'
##relyr='2018'
##gbmon='dec'
##gbyr='2017'
##rsmon='nov'
##rsyr='2017'
##rvdbfilename='U-RVDBv12.1.fasta'
divs=['VRL','ENV','HTC','INV','MAM','PLN','PRI','ROD','VRT','REFSEQ','NEIGHBOR','TPA']
### -*- coding: cp1252 -*-

def summary_stats(rvdbfilename):
    stats=defaultdict(int)
    r=0
    c=0
    l=0
    inf=open(rvdbfilename)
    for line in inf:
        l+=1
        if line.startswith('>acc'):
            c+=1
            if c-r==10000:
                print str(c)+' entries in '+rvdbfilename+' scanned'
                r=c
            sl=line.split('|')
            source=sl[1]
            if source=='REFSEQ' or source=='NEIGHBOR' or source=='TPA':
                stats[source]+=1
            if source=='GENBANK':
                for div in divs[:-3]:
                    if '|'+div+'|' in line:
                        stats[div]+=1
    inf.close()
    d=0
    for s in stats.keys():
        d+=stats[s]
    print 'number of sequences in '+rvdbfilename+': '+str(c)
    print 'number of sequences in '+rvdbfilename+' characterized by source/division: '+str(d)+' (should be same as the first line)'
    print 'number of lines in '+rvdbfilename+': '+str(l)+' (should be twice the first line)'
    return [stats,c]

def write_relnotes(stats,numentries,rvdbfilename,relday,relmon,relyr,gbmon,gbyr,rsmon,rsyr):
    outf=open(rvdbfilename.replace('.fasta','.releasenotes.txt'),'w')
    if 'U-RVDB' in rvdbfilename:
        clustered=False
        tag='U-RVDB'
    else:
        if 'C-RVDB' in rvdbfilename:
            clustered=True
            tag='C-RVDB'
        else:
            print 'the input RVDB file cannot be classified as unclustered or clustered - please check'
    if clustered:
        outf.write(tag+': '+'Clustered')
    else:
        outf.write(tag+': '+'Unclustered')
    outf.write(' Reference Viral Database\n\n\n\n')
    outf.write('RELEASENOTES '+currentvs+'\n\n')
    outf.write('Version '+currentvs+' of the ')
    if clustered:
        outf.write('Clustered')
    else:
        outf.write('Unclustered')
    outf.write(" Reference Viral Database ("+tag+"v"+currentvs+")"+"\n"+"was released on "+relmon+"/"+relday+"/"+relyr+" by Arifa Khan's Group at CBER, FDA. This update is based on the "+gbmon+", "+gbyr+" GenBank release, and the "+rsmon+", "+rsyr+" RefSeq viral and genome neighbors download.\n\n\n")
    outf.write('FEATURES OF '+tag+'\n')
    outf.write('----------------------------------------------------------------\n\n')
    outf.write('Size: \n')
    outf.write('Number of sequences/records: '+str(numentries)+'\n')
    outf.write('Number of bases: '+'\n')
    outf.write('Checksum/MD5 number:\n\n\n')
    outf.write('SUMMARY STATISTICS\n')
    outf.write('-----------------------------------------------------------------\n')
    outf.write('Genbank/Other Division	        Number of Sequences\n')
    outf.write('VRL	……………………………………………………………………… '+str(stats['VRL'])+'\n')
    outf.write('ENV	……………………………………………………………………… '+str(stats['ENV'])+'\n')
    outf.write('HTC	……………………………………………………………………… '+str(stats['HTC'])+'\n')
    outf.write('INV	……………………………………………………………………… '+str(stats['INV'])+'\n')
    outf.write('MAM	……………………………………………………………………… '+str(stats['MAM'])+'\n')
    outf.write('PLN	……………………………………………………………………… '+str(stats['PLN'])+'\n')
    outf.write('PRI	……………………………………………………………………… '+str(stats['PRI'])+'\n')
    outf.write('ROD	……………………………………………………………………… '+str(stats['ROD'])+'\n')
    outf.write('VRT	……………………………………………………………………… '+str(stats['VRT'])+'\n')
    outf.write('TPA	……………………………………………………………………… '+str(stats['TPA'])+'\n\n')
    outf.write('NCBI RefSeq Viral	…………………………… '+str(stats['REFSEQ'])+'\n')
    outf.write('NCBI Genome Neighbors   …………………………… '+str(stats['NEIGHBOR'])+'\n')
    outf.write('DISCLAIMER\n')
    outf.write('-----------------------------------------------------------------\n')
    outf.write('This material is not intended to be used as a regulatory standard.')
    outf.close()

stats,numentries=summary_stats(rvdbfilename)
write_relnotes(stats,numentries,rvdbfilename,relday,relmon,relyr,gbmon,gbyr,rsmon,rsyr)
