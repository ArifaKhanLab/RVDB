import csv
import re
wdir='F:\\VDBv10.1'

##### Given a file with one gi per row, gives you the gis
def get_gis_flatfile(flatfilename):
##    print "Getting gene identifiers from flat file: "+flatfilename
    inf=open(flatfilename)
    gis=set(inf.read().strip().split('\n'))
    inf.close()
    return gis

##### Given a fasta file, gives you the gis
def get_gis_fastafile(fasta_filename):
##    print "Getting gene identifiers from fasta file: "+fasta_filename
    gis=set([])
    inf=open(fasta_filename)
    for i,line in enumerate(inf):
        if line.startswith('>'):
            gi=line.split('|')[1]
            gis.add(gi)
    inf.close()
    return gis

##### Joins two or more sets
def join_sets(sets):
    outset=[]
    for s in sets:
        outset.extend(list(s))
    outset=set(outset)
    return outset

##### make a dictionary with a simple 1:1 mapping, where key and value are separated by the parameter 'delimiter'
def make_dict(mapping_filename,delimiter,cumulative):
##    print 'Creating dictionary from mapping in '+mapping_filename
    d=dict()
    inf=open(mapping_filename)
    for i,line in enumerate(inf):
        sl=line.strip().split(delimiter)
        key=sl[0].strip()
        value=sl[1].strip()
        if cumulative:
            try:
                existing=d[key]
            except KeyError:
                existing=[]
            existing.append(value)
            d[key]=existing
        else:
            d[key]=value
    inf.close()
    return d

##### Given a headers file, gives you the gis
def get_gis_headersfile(headers_filename):
    gis=set([])
    inf=open(headers_filename)
    for i,line in enumerate(inf):
        if line.startswith('>'):
            gi=line.strip().split('|')[1]
            gis.add(gi)
    inf.close()
    return gis

#################################################################################
######## Simple function to retrieve gene identifiers from a RefSeq file ########
#################################################################################
def get_refseq_jan2016_gis(refseq_filename):
##    print "Fetching RefSeq (Jan 2016 version) gis for non-phage viral entries"
    inf=open(refseq_filename)
    ncbigis=set([])
    for i,line in enumerate(inf):
        if line.startswith('>gi'):
            gi=line.split('|')[1]
            ncbigis.add(gi)
    inf.close()
    return rsgis


################################################################################
######## Simple function to retrieve descriptions from a VDB fasta file ########
################################################################################
def get_vdb_desc(vdb_filename):
##    print "Getting VDB headers dict using "+vdb_filename
    descs=dict()
    inf=open(vdb_filename)
    for i,line in enumerate(inf):
        if line.startswith('>gi'):
            gi=line.split('|')[1]
            desc='|'.join(line.strip().split('|')[4:])
            descs[gi]=desc
    return descs

####################################################################################################
######## Simple function to retrieve gi and descriptions from a VDB descriptions quick file ########
####################################################################################################
def get_vdb_desc_quick(vdb_desc_filename):
##    print "Getting VDB headers dict quickly using "+vdb_desc_filename
    inf=open(vdb_desc_filename)
    descs=dict()
    for i,line in enumerate(inf):
        gi,desc=line.strip().split('\t')
        descs[gi]=desc
    inf.close()
    return descs

####################################################################################################
######## Simple function to reformat a raw NCBI Viral fast file to have 2 lines /entry      ########
######## and output two files, phage and eukaryotic viral + archaea                         ########
####################################################################################################
def parse_raw_NCBI_Viral(raw_ncbi_viral_filename):
##    print 'parsing '+raw_ncbi_viral_filename+' into phage and eukaryotic viral + archael sequences'
    outf1_name=raw_ncbi_viral_filename.replace('.fasta','.eukviral.fasta')
    outf2_name=raw_ncbi_viral_filename.replace('.fasta','.phage.fasta')
    outf1=open(outf1_name,'w')
    outf2=open(outf2_name,'w')
    phage_names_inf=open('F:\\NCBI_viral\\phage_kws.txt')
    phage_names=phage_names_inf.read().strip().split('\n')
    phage_names_inf.close()
    inf=open(raw_ncbi_viral_filename)
    phage=False
    p=0
    e=0
    for i,line in enumerate(inf):
        if i==0:
            entry=line
            continue
        if line.startswith('>'):
            if phage:
                outf2.write(entry)
                p+=1
                phage=False
            else:
                outf1.write(entry)
                e+=1
            entry=line
            header=line.lower()
            for phage_name in phage_names:
                if phage_name in header:
                    phage=True
                    break
        else:
            entry+=line
    if phage:
        outf2.write(entry)
        p+=1
    else:
        outf1.write(entry)
        e+=1
##    print "read and reformatted "+str(e)+" eukaryotic viral + archaeal viral sequences"
##    print "read and reformatted "+str(p)+" phage sequences"
    outf1.close()
    outf2.close()
    
##################################################################################################
########## Simple function to extract clusters from clustering file                    ###########
##################################################################################################
def get_clusters(cluster_infilename):
##    print "Retrieving clusters from file: "+cluster_infilename
    inf=open(cluster_infilename)
    clusters=inf.read().strip().split('>Cluster')[1:]
    inf.close()
    return clusters
    
##################################################################################################
########## Simple function to extract headers from .fasta file                         ###########
##################################################################################################
def get_headers(fasta_filename,filterset):
##    print "Retrieving headers from file: "+fasta_filename
    headers=[]
    inf=open(fasta_filename)
    for i,line in enumerate(inf):
        if line.startswith('>'):
            gi=line.strip().split('|')[1]
            if len(filterset)>0:
                if not gi in filterset:
                    continue
            headers.append(line)
    inf.close()
    return headers

def read_file_byblock(entries,infilename,blockparsefn,delim):
##    print "Reading "+infilename+" by block"
    inf=open(infilename)
    nlines=100000
    block=''
    n=0
    lastblock=False
    for i,line in enumerate(inf):
        if n<nlines:
            block+=line
        if n>=nlines:
            entries,block=blockparsefn(entries,block,delim,lastblock)
            n=len(block.split('\n'))
        n+=1
    lastblock=True
    entries,block=blockparsefn(entries,block,delim,lastblock)
    return entries

####### Sub-function for read_file_byblock when parsing GenBank raw entries
####### ... passed as 3rd argument to read_file_byblock master function
def parse_genbank_block(entries,gb_block,delim,lastblock):
    records=gb_block.strip().split(delim)
    finalrecind=len(records)-1
    regexes=genbank_record_regexes()
##    print "initiating block parse with "+str(finalrecind+1)+" records"
    for r,record in enumerate(records):
        if not lastblock:
            if r==finalrecind:
                gb_block=record
                return [entries,gb_block]
            else:
                gi=record.split('VERSION')[1].split('\n')[0].split('GI:')[1].strip()
                entries[gi]=parse_genbank_record(record,regexes)
        else:
            gi=record.split('VERSION')[1].split('\n')[0].split('GI:')[1].strip()
            entries[gi]=parse_genbank_record(record,regexes)
            gb_block=''
    return [entries,gb_block]

####### Sub-function for parse_genbank_block, parses individual GenBank record
def parse_genbank_record(gb_record,regexes):
    byrec=dict()
    features=[]
    featuretitlefinder=regexes[0]
    try:
        featureinfo=gb_record.strip().split('\nFEATURES')[1].split('ORIGIN')[0].strip()
        featureinfo='\n'.join(featureinfo.strip().split('\n')[1:])
    except IndexError:
        byrec[gi]=features
        return byrec
    gi=gb_record.split('VERSION')[1].split('\n')[0].split('GI:')[1].strip()
    featuretitlelines=featuretitlefinder.findall(featureinfo)
    find=len(featuretitlelines)-1
    for f,featuretitleline in enumerate(featuretitlelines):
        try:
            start,end=featuretitleline.strip().split()[1].split('..')
        except ValueError:
            continue
        foundfeatures=[]
        if f<find:
            featureblock=featureinfo.split(featuretitleline)[1].split(featuretitlelines[f+1])[0].strip()
        else:
            featureblock=featureinfo.split(featuretitleline)[1].strip()
        featurelines=featureblock.strip().split('=')[1:]
        for featureline in featurelines:
            sd=simplify_description(featureline)
            viral=False
            viral,poskws=pos_screen(sd)
            if viral:
                try:
                    feature=featureline.split('"')[1].strip()
                except IndexError:
                    feature=featureline.strip()
                feature=' '.join(feature.split())
                foundfeatures.append(feature)
        features.append(['; '.join(foundfeatures),[start,end]])
    byrec['feature']=features
    return byrec

##### Enter regular expression for detecting "anatomy" of the GenBank raw entries
def genbank_record_regexes():
    print "Initiating GenBank record regular expression list"
    import re
    featuretitlefinder=re.compile('     [a-z]+_*[a-z]+\s+[\<\>]*[0-9]+..[\<\>]*[0-9]+')
    regexes=[featuretitlefinder]
    return regexes

####### Prints out information parsed from raw GenBank entries file
def print_out_gbinfo(entries,infotypes,outfilename):
    print "Printing results file with GenBank entry information: "+outfilename
    outf=open(outfilename,'w')
    gis=entries.keys()
    for gi in gis:
        entry=entries[gi]
        outf.write('>gi|'+gi+'\n')
        for infotyp in infotypes:
            outf.write(infotyp.upper()+'\n')
            infos=entry[infotyp]
            if infotyp=='feature':
                outf.write(str(sdict[gi][1])+'\n')
                for info in infos:
                    start,end=info[1]
                    desc=info[0]
                    if not desc=='':
                        outf.write(desc+';'+start+'-'+end+'\n')
        outf.write('\n')

##def iseven(integer):
##    if integer/2==float(integer)/2:
##        return True
##    else:
##        return False

########## Takes a standard GenBank entry description and simplifies the format - all non-alphanumeric characters replaced by a single whitespace character each (' '),
########## underscore is replaced by a single whitespace  (' '), the description is made lowercase, then a single whitespace is added to either side
########## is added to the beginning and end
def simplify_description(description):
    description_s=re.sub(r'\W+',' ',description)
    description_s=description_s.replace('_',' ')
    description_s=description_s.lower()
    description_s=description_s.strip()
    description_s=' '+description_s+' '
    return description_s

####### Collects sequences, sequence lengths from a given fastafile and gimask input
def fetch_seqs(fastafilename,gimask):
    print "Fetching sequences for "+fastafilename+" using a GI mask of length: "+str(len(gimask))
    sdict=dict()
    inf=open(fastafilename)
    match=False
    for i,line in enumerate(inf):
        if line.startswith('>gi'):
            gi=line.strip().split('|')[1]
            if gi in gimask:
                match=True
            else:
                match=False
        else:
            if match:
                seq=line.strip()
                slen=len(seq)
                sdict[gi]=[seq,slen]
    return sdict

####### Collects gis or accs from a given input file, either a list of identifiers
####### .... or a fasta file
def collect_gis(infilename,infiletype,idtype):
    ids=set([])
    inf=open(infilename)
    if infiletype=='simple':
        ids=set(inf.read().strip().split('\n'))
    if infiletype=='fasta':
        for i,line in enumerate(inf):
            if line.startswith('>gi'):
                if idtype=='gi':
                    thisid=line.strip().split('|')[1]
                if idtype=='ac':
                    thisid=line.strip().split('|')[3]
                ids.add(thisid)
    if infiletype=='clstr':
        clusters=inf.read().strip().split('>Cluster')[1:]
        for cluster in clusters:
            csplit=cluster.strip().split('|')
            if idtype=='gi':
                thisid=csplit[1].strip()
            if idtype=='ac':
                thisid=csplit[3].strip()
            ids.add(thisid)
    inf.close()
    return ids

####### Discovers which GIs in the input GenBank entries feature file have stretches of sequence without feature annotation
def lacking_features(infilename,outfilename):
    inf=open(infilename)
    entries=inf.read().strip().split('>gi|')[1:]
    numfinder=re.compile('[0-9]+')
    outf=open(outfilename,'w')
    for entry in entries:
        slen=int(entry.strip().split('\n')[2])
        gaps=[[1,slen]]
        features=entry.strip().split('\n')[3:]
        for feature in features:
            feature=feature.strip().split(';')[-1]
            fstart,fend=numfinder.findall(feature)
            fstart=int(fstart)
            fend=int(fend)
            gaps=update_gaps(gaps,fstart,fend)
        if len(gaps)>0:
            outf.write('>gi|'+entry)
    outf.close()

####### Given a stretch of sequence (e.g. for an additional feature), updates remaining gap or areas not covered by any such stretches so far
def update_gaps(gaps,fstart,fend):
    newgaps=[]
    for gap in gaps:
        gapstart,gapend=gap
        if gapstart<fstart:
            newgaps.append([gapstart,fstart])
        if gapend>fend:
            newgaps.append([fend,gapend])
##    print 'NOW LOOK'
##    print gaps
##    print newgaps
    return newgaps

##### Given a file in multi-block fasta format, reformats with sequence in a single line only #####
def reformat_fasta_1lineseq(in_fastafile_name,dbtag):
    out_fastafile_name=in_fastafile_name.replace('.fna','_formatted.fna')
    outf=open(out_fastafile_name,'w')
    inf=open(in_fastafile_name)
    seq=''
    for i,line in enumerate(inf):
        if line.startswith('>'+dbtag):
            if not seq=='':
                outf.write(header+seq+'\n')
            header=line
            seq=''
        else:
            seq+=line.strip()
    outf.write(header+seq)
    inf.close()
    outf.close()
##
####### Creates a dictionary from a file with two entities per line #####
##def make_dict(mapfile,reverse):
##    d=dict()
##    inf=open(mapfile)
##    for i,line in enumerate(inf):
##        sl=line.split()
##        if reverse:
##            d[sl[1].strip()]=sl[0].strip()
##        else:
##            d[sl[0].strip()]=sl[1].strip()
##    inf.close()
##    return d

##### Calculates length of all sequences in input fasta file and writes output
def calculate_seqlength(infastafilename):
    outfilename=infastafilename.replace('.fasta','_lens.txt')
    l=0
    outlens=[]
    inf=open(infastafilename)
    for i,line in enumerate(inf):
        if line.startswith('>gi'):
            if not l==0:
                outlens.append([gi,l])
            gi=line.strip().split('|')[1].strip()
        else:
            l+=len(line.strip())
    outlens.append([gi,l])
    inf.close()
    import operator
    order=operator.itemgetter(1)
    outlens=sorted(outlens,key=order,reverse=True)
    outf=open(outfilename,'w')
    for outlen in outlens:
        outf.write(outlen[0]+' '+str(outlen[1])+'\n')
    outf.close()

##### Reads in as input the output of the function above, storing sequence lengths
##### as a dictionary
def fetch_seqlength(seqlen_filename,filtergis):
    lens=dict()
    inf=open(seqlen_filename)
    for i,line in enumerate(inf):
        gi,l=line.strip().split()
        if len(filtergis)>0:
            if not gi in filtergis:
                continue
        lens[gi]=int(l)
    inf.close()
    return lens

##### Removes all sequences under a given length, writes two output files:
##### One for the removed sequences, the other for the remainder sequences
def extract_short_seqs(inputfastafilename,lens,dist_cutoff):
    print 'Extracting short sequences (under '+str(dist_cutoff)+' b.p.) from file: '+inputfastafilename
    out_shortseqs_filename=inputfastafilename.replace('.fasta','_shortseqs_'+str(dist_cutoff)+'.fasta')
    out_shortseqsremainder_filename=inputfastafilename.replace('.fasta','_shortseqs_'+str(dist_cutoff)+'_'+'remainder.fasta')
    outf1=open(out_shortseqs_filename,'w')
    outf2=open(out_shortseqsremainder_filename,'w')
    match=False
    inf=open(inputfastafilename)
    for i,line in enumerate(inf):
        if line.startswith('>gi'):
            gi=line.strip().split('|')[1]
            seqlen=lens[gi]
            if seqlen<=dist_cutoff:
                match=True
            else:
                match=False
        if match:
            outf1.write(line)
        else:
            outf2.write(line)
    inf.close()
    outf1.close()
    outf2.close()

##### Removes all sequences with the given set of keywords, writes two output files:
##### One for removed sequences, the other for the remainder sequences
def extract_kw_seqs(inputfastafilename,kws,tag):
    print 'Extracting sequences with '+tag+'-related keywords from file: '+inputfastafilename+'kws: '+', '.join(kws)
    out_withkw_filename=inputfastafilename.replace('.fasta','_'+tag+'+'+'.fasta')
    out_withoutkw_filename=inputfastafilename.replace('.fasta','_'+tag+'-'+'.fasta')
    outf1=open(out_withkw_filename,'w')
    outf2=open(out_withoutkw_filename,'w')
    match=False
    inf=open(inputfastafilename)
    for i,line in enumerate(inf):
        if line.startswith('>gi'):
            gi=line.strip().split('|')[1]
            desc=' '.join(line.strip().split('|')[4:])
            ds=simplify_description(desc)
            match=False
            for kw in kws:
                if kw in ds:
                    match=True
                    break
        if match:
            outf1.write(line)
        else:
            outf2.write(line)
    inf.close()
    outf1.close()
    outf2.close()

def invert_dict(dict1):
    dict2=dict()
    for key1 in dict1.keys():
        values1=dict1[key1]
        for val1 in values1:
            try:
                existing=dict2[val1]
            except KeyError:
                existing=[]
            existing.append(key1)
            dict2[val1]=existing
    return dict2

def combine_dicts(dict1,dict2,islist):
    print 'Combining dictionaries of length '+str(len(dict1))+' and length '+str(len(dict2))
    dict3=dict()
    keys1=dict1.keys()
    keys2=dict2.keys()
    keys3=[]
    keys3.extend(keys1)
    keys3.extend(keys2)
    keys3=list(set(keys3))
    for key in keys3:
        entries3=[]
        try:
            entries1=dict1[key]
            if islist:
                entries3.extend(entries1)
            else:
                entries3.append(entries1)
        except KeyError:
            a=1
        try:
            entries2=dict2[key]
            if islist:
                entries3.extend(entries2)
            else:
                entries3.append(entries2)
        except KeyError:
            a=1
        entries3=list(set(entries3))
        dict3[key]=entries3
    return dict3

def combine_sets(set1,set2):
    set3=[]
    set3.extend(list(set1))
    set3.extend(list(set2))
    set3=set(set3)
    return set3

def combine_listoflists(lol1,lol2):
    found=set([])
    lol3=[]
    for pair in lol1:
        opair=list(pair)
        opair=sorted(opair)
        tag='_'.join(sorted(opair))
        if tag in found:
            continue
        found.add(tag)
        lol3.append(pair)
    for pair in lol2:
        opair=list(pair)
        opair=sorted(opair)
        tag='_'.join(sorted(opair))
        if tag in found:
            continue
        found.add(tag)
        lol3.append(pair)
    return lol3

def combine_dict_listoflists(dict1,dict2):
    print 'Combining dictionaries of length '+str(len(dict1))+' and length '+str(len(dict2))
    dict3=dict()
    keys1=dict1.keys()
    keys2=dict2.keys()
    keys3=[]
    keys3.extend(keys1)
    keys3.extend(keys2)
    keys3=list(set(keys3))
    for key in keys3:
        try:
            lol1=dict1[key]
        except KeyError:
            lol1=[]
        try:
            lol2=dict2[key]
        except KeyError:
            lol2=[]
        lol3=combine_listoflists(lol1,lol2)
        dict3[key]=lol3
    return dict3

def collate_dict_allvalues(d1):
    vals=[]
    for key in d1.keys():
        vals.extend(d1[key])
    vals=set(vals)
    return vals    

######### Collects sequences and sequence lengths for all sequences in gimask
##infilename='F:\\VDBv9.1\\Clustering\\'+'VDBv9.1.creps.retro.clstr'
##infiletype='clstr'
##idtype='gi'
##gimask1=collect_gis(infilename,infiletype,idtype)
##infilename='F:\\VDBv9.1\\Clustering\\'+'VDBv9.1.creps.retro.clstr'
##infiletype='clstr'
##idtype='gi'
##gimask2=collect_gis(infilename,infiletype,idtype)
##gimask=[]
##gimask.extend(list(gimask1))
##gimask.extend(list(gimask2))
##gimask=set(gimask)
##del gimask1
##del gimask2
##fastafilename='F:\\VDBv9.1\\'+'rVDBv9.1.creps+rhabdo_033116s.fasta'
##sdict=fetch_seqs(fastafilename,gimask)
######### Collects only those features / GI that have positive keywords in them
##import sys
##sys.path.append('F:\\TOOLBOX')
##from rvdbkeywords import *
##wdir='F:\\VDBv9.1\\GenBank_format_entries_retro'
##delim='\n//\n'
##import os
##entries=dict()
##for fn in os.listdir(wdir):
##    if 'parsed' in fn:
##        continue
##    if not 'VDBv9.1.creps.retro' in fn:
##        continue
##    print fn
##    gbrecs_filename=wdir+'\\'+fn
##    entries=read_file_byblock(entries,gbrecs_filename,parse_genbank_block,delim)
####outfilename=wdir+'\\'+'VDBv9.1.creps.retro.featuresparsed2.txt'
####infotypes=['feature']
####print_out_gbinfo(entries,infotypes,outfilename)
##
####infilename='F:\\VDBv9.1\\GenBank_format_entries_viral\\VDBv9.1.creps.notinRodney90k.viral.featuresparsed2.txt'
####outfilename='F:\\VDBv9.1\\GenBank_format_entries_viral\\VDBv9.1.creps.notinRodney90k.viral.featuresparsed2_withflanking.txt'
####lacking_features(infilename,outfilename)

########## Block of code below is to format a fasta file in which the sequence covers multiple lines (length 60 or other)
########## as a file in which sequence is in one line only
##in_fastafile_name='F:\\REFSEQ\\viral.1.1.genomic_july2016.fna'
##dbtag='gi'
##reformat_fasta_1lineseq(in_fastafile_name,dbtag)

####### Block of code below is to extract sequences 100 bp or less
##inputfastafilename=wdir+'\\'+'VDBv10.0.raw.nr.fasta'
##lens_filename=inputfastafilename.replace('.fasta','_lens.txt')
##lens=fetch_seqlength(lens_filename)
##dist_cutoff=100
##extract_short_seqs(inputfastafilename,lens,dist_cutoff)
####inputfastafilename=wdir+'\\'+'VDBv10.0.creps.fasta')
####extract_short_seqs(inputfastafilename,lens,dist_cutoff)
##
####### Block of code below is to extract LTR only sequences
##inputfastafilename=wdir+'\\'+'VDBv10.0.raw.nr.fasta'
####lens=fetch_seqlength(lens_filename)
##dist_cutoff=1000
##extract_short_seqs(inputfastafilename,lens,dist_cutoff)
##kws=["terminal repeat"," ltr "," 5' end "," 3' end "]
##tag='LTRonly'
##inputfastafilename=wdir+'\\'+'VDBv10.0.raw.nr_shortseqs_1000.fasta'
##extract_kw_seqs(inputfastafilename,kws,tag)
##inputfastafilename='F:\\VDBv10.0\\VDBv10.0.creps.fasta'
##extract_short_seqs(inputfastafilename,lens,dist_cutoff)
##inputfastafilename='F:\\VDBv10.0\\VDBv10.0.creps_shortseqs_1000.fasta'
##extract_kw_seqs(inputfastafilename,kws,tag)

##vdb_filename='F:\\VDBv10.1\\VDBv10.1.raw.fasta'
##descs=get_vdb_desc(vdb_filename)

##raw_ncbi_viral_filename='F:\\NCBI_viral\\NCBI_viral_05222017_raw.fasta'
##parse_raw_NCBI_Viral(raw_ncbi_viral_filename)
