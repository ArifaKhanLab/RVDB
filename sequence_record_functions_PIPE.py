import csv
import re

##### Given a file with one gi per row, gives you the accs
def get_accs_flatfile(flatfilename):
##    print "Getting gene identifiers from flat file: "+flatfilename
    inf=open(flatfilename)
    accs=set(inf.read().strip().split('\n'))
    inf.close()
    return accs

##### Given a fasta file, gives you the accs
def get_accs_fastafile(fasta_filename):
##    print "Getting gene identifiers from fasta file: "+fasta_filename
    accs=set([])
    inf=open(fasta_filename)
    for i,line in enumerate(inf):
        if line.startswith('>'):
            acc=line.split('|')[1]
            accs.add(acc)
    inf.close()
    return accs

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

################################################################################
######## Simple function to retrieve descriptions from a VDB fasta file ########
################################################################################
def get_vdb_desc(vdb_filename):
##    print "Getting VDB headers dict using "+vdb_filename
    descs=dict()
    inf=open(vdb_filename)
    for i,line in enumerate(inf):
        if line.startswith('>gi'):
            acc=line.split('|')[1]
            desc='|'.join(line.strip().split('|')[4:])
            descs[acc]=desc
    return descs

#################################################################################
######## Simple function to retrieve gene identifiers from a RefSeq file ########
#################################################################################
def get_refseq_accs(refseq_filename):
##    print "Fetching RefSeq accs for non-phage viral entries"
    inf=open(refseq_filename)
    rsaccs=set([])
    for i,line in enumerate(inf):
        if line.startswith('>acc'):
            acc=line.split('|')[1]
            rsaccs.add(acc)
    inf.close()
    return rsaccs

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
def fetch_seqs(fastafilename,filterset):
    print "Fetching sequences for "+fastafilename+" using a acc mask of length: "+str(len(filterset))
    sdict=dict()
    inf=open(fastafilename)
    match=False
    for i,line in enumerate(inf):
        if line.startswith('>acc'):
            acc=line.strip().split('|')[1]
            if acc in filterset:
                match=True
            else:
                match=False
        else:
            if match:
                seq=line.strip()
                slen=len(seq)
                sdict[acc]=[seq,slen]
    return sdict

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
            acc=line.strip().split('|')[1].strip()
        else:
            l+=len(line.strip())
    outlens.append([acc,l])
    inf.close()
    import operator
    order=operator.itemgetter(1)
    outlens=sorted(outlens,key=order,reverse=True)
    outf=open(outfilename,'w')
    for outlen in outlens:
        outf.write(outlen[0]+' '+str(outlen[1])+'\n')
    outf.close()

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

def collate_dict_allvalues(d1):
    vals=[]
    for key in d1.keys():
        vals.extend(d1[key])
    vals=set(vals)
    return vals    

