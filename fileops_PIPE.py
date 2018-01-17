


def segment_file(infilename,nrows,recstart,recend):
    a=0
    b=0
    filetype='.'+infilename.split('.')[-1]
    inf=open(infilename)
    segmented_filenames=[]
##  print 'Counting number of rows from input file: '+infilename
    for i,line in enumerate(inf):
        a=1
    totalrows=i
    inf.close()
    starts=range(0,totalrows,nrows)
    for s,start in enumerate(starts):
        inf=open(infilename)
        outfilename=infilename.replace(filetype,'_'+str(start)+'_'+str(start+nrows)+filetype)
        segmented_filenames.append(outfilename)
##      Writing output file
        outf=open(outfilename,'w')
        r=0
        for i,line in enumerate(inf):
            if s==0:
                if i<start:
                    continue
            else:
                if i<firstnewrow:
                    continue
            r+=1
            outf.write(line)
            if line.startswith(recend):
                if r>nrows:
                    firstnewrow=i+1
                    break
        inf.close()
        outf.close()
    return segmented_filenames

import sys

homedir=sys.argv[1]
date=sys.argv[2]
currentvs=sys.argv[3]
wdir=homedir+'\\'+'RVDBv'+currentvs
refseqdir=wdir+'\\'+'RefSeq_raw_data_'+date
input_filename=refseqdir+'\\'+'viral.genomic.gbff'
filetype=sys.argv[4]
if filetype=='gbff':
    recstart='LOCUS'
    recend='//\n'
nrows=int(sys.argv[5])
segmented_filenames=segment_file(input_filename,nrows,recstart,recend)
for segmented_filename in segmented_filenames:
    print segmented_filename
    
                    
                
    
