import csv
import re
import sys
homedir=sys.argv[1]
datetag=sys.argv[2]
currentvs=sys.argv[3]
wdir=homedir+'\\'+'RVDBv'+currentvs
urvdb_filename=wdir+'\\'+'U-RVDBv'+currentvs+'.fasta'
raw_clustering_infilename=wdir+'\\'+sys.argv[4]


###########################################################################################
########## First function for parsing a raw CD-HIT-EST clustering output file   ###########
########## repairs headers that were curtailed during the clustering run        ###########
########## and promotes either RefSeq Viral or neighbor sequence (respectively) ###########
########## to cluster representative status                                     ###########
###########################################################################################
def write_promoted_clusters(raw_clustering_infilename,descs):
    print "Writing cluster-representative output file using: "+raw_clustering_infilename+" as input"
    inf=open(raw_clustering_infilename)
    clusters=inf.read().split('>Cluster')[1:]
    inf.close()
    creps=[]
    endline_examiner=re.compile('[\Wa-zA-Z0-9_:]*\.\.\.')
    for cluster in clusters:
        lines=cluster.strip().split('\n')
        nl=len(lines)
        hasrs=False
        hasncbi=False
        crepindex='NA'
        outlines=[]
        for i,line in enumerate(lines):
            line=line.strip()
            if i==0:
                cnum=line
                continue
            sdesc=line.split('|')
            acc=sdesc[2]
            desc=descs[acc]
            sl=line.split('|')
            endline_slot=sl[-1]
            try:
                endline_chars=endline_examiner.findall(endline_slot)[0]
            except IndexError:
                print line
            try:
                line=line.replace(endline_chars,desc+'|...')
            except TypeError:
                print line
##            line=line.replace(endline,'|'+desc+'|...')
            lines[i]=line
            if '|REFSEQ|' in line:
                hasrs=True
                crepindex=i
                if line.endswith('*'):
                    crep=line.replace('*','*crep*')
                else:
                    crep=line+' *crep*'
            if not hasrs and '|NEIGHBOR|' in line:
                if not 'patent' in line:
                    hasncbi=True
                    crepindex=i
                    if line.endswith('*'):
                        crep=line.replace('*','*crep*')
                    else:
                        crep=line+' *crep*'
            if not hasrs and not hasncbi and line.endswith('... *'):
                crepindex=i
                if line.endswith('*'):
                    crep=line.replace('*','*crep*')
                else:
                    crep=line+' *crep*'
        try:
            outlines=[crep]
        except UnboundLocalError:
            print line
            continue
        try:
            outlines.extend(lines[1:crepindex])
        except TypeError:
            print outlines
            print lines
            print crepindex
        try:
            outlines.extend(lines[crepindex+1:])
        except IndexError:
            print ''.join(lines)
            print '\n\n'
        creps.append('>Cluster '+cnum+'\n'+'\n'.join(outlines))
    promoted_clustering_outfilename=raw_clustering_infilename.replace('.clstr','.promoted.clstr')
    outf=open(promoted_clustering_outfilename,'w')
    outf.write('\n'.join(creps)+'\n')
    outf.close()

###########################################################################################
##### Aux function #1 for first CD-HIT-EST parsing script ("write_promoted_clusters"  #####
##### collects descriptions for all U-RVDB (of the present update) sequences          #####
def collect_rvdb_descriptions(urvdb_filename):
    descs=dict()
    inf=open(urvdb_filename)
    for line in inf:
        if line.startswith('>acc'):
            sl=line.strip().split('|')
            acc=sl[2]
            desc='|'.join(line.strip().split('|')[3:])
            descs[acc]=desc
    inf.close()
    return descs

##############################################################################################
########## Second function for parsing a raw CD-HIT-EST clustering output file      ##########
########## writes a second final output file containing just creps                  ##########
##############################################################################################
def write_crep_outfile(promoted_clustering_filename):
    print "Writing creps outfile using: "+promoted_clustering_filename+" as input"
    creps_outfilename=promoted_clustering_filename.replace('.promoted.clstr','.creps.clstr')
    outf=open(creps_outfilename,'w')
    inf=open(promoted_clustering_filename)
    clusters=inf.read().strip().split('>Cluster')[1:]
    for cluster in clusters:
        cnum=cluster.strip().split('\n')[0].strip()
        lines=cluster.strip().split('\n')[1:]
        for line in lines:
            if '*crep*' in line:
                outf.write('>Cluster '+cnum+'\n'+line+'\n')
    outf.close()

descs=collect_rvdb_descriptions(urvdb_filename)
write_promoted_clusters(raw_clustering_infilename,descs)
