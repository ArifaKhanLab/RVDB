import sqlite3

## Initiaties db
def create_new_db(dbname,tablenames,tablefields,tabletypes,pkeyindex):
# Connecting to the database file
    print wdir+'/'+dbname
    conn = sqlite3.connect(wdir+'/'+dbname)
    c = conn.cursor()
# Creating SQLite tables according to tablenames (each initialized with 1 column as in tablefields, data type as in tabletypes)
    for t,tablename in enumerate(tablenames):
        fieldname=tablefields[t]
        ftype=tabletypes[t]
        if t==pkeyindex-1:
            c.execute('CREATE TABLE {tn} ({f} {ft} PRIMARY KEY)'\
                .format(tn=tablename, f=fieldname, ft=ftype))
        else:
            c.execute('CREATE TABLE {tn} ({f} {ft})'\
                .format(tn=tablename, f=fieldname, ft=ftype))
    conn.commit()
    conn.close()


## Adds columns
def alter_db(dbname,tablename,fieldnames,fieldtypes):
# Connecting to the database file
    conn = sqlite3.connect(wdir+'\\'+dbname)
    c = conn.cursor()
# Adding columns to the SQLite table passed as "tablename"; these columns are named "fieldnames" and are of type "fieldtypes"
    for f,fieldname in enumerate(fieldnames):
        ftype=fieldtypes[f]
        c.execute("ALTER TABLE {tn} ADD COLUMN '{f}' {ft}"\
                .format(tn=tablename, f=fieldname, ft=ftype))
    conn.commit()
    conn.close()

## Adds rows
def build_db(rvdb_dir,rvdb_filename):
# Connecting to the database file
    conn = sqlite3.connect(wdir+'\\'+dbname)
    c = conn.cursor()
# Opening the RVDB file, to be able to population the sqlite database
    inf=open(rvdb_dir+'\\'+rvdb_filename)
    c1=0
    c2=0
    seqlen='NA'
    chardict=getchar(rvdb_dir,rvdb_filename)
    for line in inf:
        if line.startswith('>acc'):
            if not seqlen=='NA':
                c.execute("INSERT INTO "+str(maintable)+" ("+str(field1)+", "+str(field2)+", "+str(field3)+", "+str(field4)+", "+str(field5)+", "+str(field6)+", "+str(field7)+", "+str(field8)
                +") VALUES (?,?,?,?,?,?,?,?)", (acc,source,description,seqlen,organism,division,date,biocategory))
                c1+=1
                if c1-c2==100000:
                    print str(c1)+' entries from '+rvdb_filename+' parsed into sqlite3 format'
                    c2=c1
            sl=line.strip().split('|')
            source=sl[1]
            acc=sl[2]
            seqlen=0
            biocategory=chardict[acc]
            if source=='REFSEQ' or source == 'TPA':
                description=' '.join(sl[3:])
                organism='NA'
                division='NA'
                date='NA'
            else:
                description=' '.join(sl[3:-3])
                organism=sl[-3]
                division=sl[-2]
                date=sl[-1]
        else:
            seqlen+=len(line.strip())
    inf.close()
    conn.commit()
    conn.close()

###############################################################################################################
########## Retrieves characeterized sequences                                                       ###########
########## (i.e. output of write_characterization_output from RVDB_characterization.py script)      ###########
########## this is needed for the build_db function above                                           ###########
###############################################################################################################
def getchar(wdir,rvdb_filename):
    ex_fn=wdir+'\\'+rvdb_filename+'.EX.headers.txt'
    ex_inf=open(ex_fn)
    enrv_fn=wdir+'\\'+rvdb_filename+'.ENRV.headers.txt'
    enrv_inf=open(enrv_fn)
    erv_fn=wdir+'\\'+rvdb_filename+'.ERV.headers.txt'
    erv_inf=open(erv_fn)
    retro_fn=wdir+'\\'+rvdb_filename+'.LTR-RETRO.headers.txt'
    retro_inf=open(retro_fn)
    other_fn=wdir+'\\'+rvdb_filename+'.UNASSIGNED.headers.txt'
    other_inf=open(other_fn)
    cdict2=dict()
    ex_heads=ex_inf.read().strip().split('\n')
    ex_inf.close()
    enrv_heads=enrv_inf.read().strip().split('\n')
    enrv_inf.close()
    erv_heads=erv_inf.read().strip().split('\n')
    erv_inf.close()
    retro_heads=retro_inf.read().strip().split('\n')
    retro_inf.close()
    other_heads=other_inf.read().strip().split('\n')
    other_inf.close()
    for ex_head in ex_heads:
        acc=ex_head.split('|')[2]
        cdict2[acc]='EX'
    for enrv_head in enrv_heads:
        acc=enrv_head.split('|')[2]
        cdict2[acc]='ENRV'
    for erv_head in erv_heads:
        acc=erv_head.split('|')[2]
        cdict2[acc]='ERV'
    for retro_head in retro_heads:
        acc=retro_head.split('|')[2]
        cdict2[acc]='LTR-RETRO'
    for other_head in other_heads:
        acc=other_head.split('|')[2]
        cdict2[acc]='OTHER'
    return cdict2

def parse_sysargvs(sysargvs):
    s_args=sysargvs.split()
    conditions=dict()
    for s_arg in s_args:
        field,vals=s_arg.split('=')
        vals=s_arg.split('[')[1].split(']')[0].split(',')
        conditions[field]=vals
    return conditions

## "("+str(field1)+", "+str(field2)+", "+str(field3)+", "+str(field4)+", "+str(field5)+", "+str(field6)+")"

def query_db(dbname,outfilename,conditions,boolcond):
# Connecting to the database file
    conn = sqlite3.connect(wdir+'\\'+dbname)
    c = conn.cursor()
    executestring="SELECT "+"*"+" FROM "+str(maintable)
    conditions_string=[]
    for field in conditions.keys():
        field_conditions=conditions[field]
        if field=='description':
            for field_condition in field_conditions:
                conditions_string.append(field+" LIKE "+"'"+"%"+field_condition+"%""'")
        if field=='source':
            for field_condition in field_conditions:
                conditions_string.append(field+"="+"'"+field_condition+"'")
    boolcond=" "+boolcond.upper()+" "
    executestring+="\nWHERE "+boolcond.join(conditions_string)
    print executestring
    c.execute(executestring)
    rows=c.fetchall()
    outf=open(wdir+'\\'+outfilename,'w')
    for row in rows:
        print row
        outf.write('|'.join(row)+'\n')
    outf.close()

### Collect arguments for the directory of the update; also adds path to scripts for update
import sys
sys.path.append('E:/UPDATE_SCRIPTS_LOGS')
homedir=sys.argv[1]
datetag=sys.argv[2]
currentvs=sys.argv[3]
wdir=homedir+'/'+'RVDBv'+currentvs
rvdb_filename=sys.argv[4]

### Initialize sqlite database
dbname=rvdb_filename.replace('.fasta','.sqlite.db') # name of the sqlite db file to be created
maintable = 'rvdb'  # name of the first table, which contains RVDB information
pkeytable = 'accessions'  # name of the second table to be created, which contains RVDB accessions
tablenames=[maintable,pkeytable]
pkeyindex=2 # states that the second table contains the PRIMARY KEYS
field1 = 'accs' # name of the column containing the PRIMARY KEY values - the accession, with version number (*.[1-9])
tablefields=[field1,field1]
tabletypes=['TEXT','TEXT']

### Create sqlite database for RVDB; "create_new_db" creates the tables, alter_db adds fields
print dbname
create_new_db(dbname,tablenames,tablefields,tabletypes,pkeyindex)
field2 = 'source' # name of the column containing the provenance information - GENBANK, REFSEQ, NEIGHBOR, TPA
field3 = 'description'  # name of the column containing the entry description
field4 = 'seqlen' # name of the column containing the sequence length
field5 = 'organism' # name of the column containing the organism (absent from REFSEQ entries)
field6 = 'division' # name of the column containing the GenBank division (absent from REFSEQ and TPA entries)
field7 = 'date' # name of the column containing the date of submission of the entry (absent from REFSEQ entries because REFSEQ are based on prior entries in GenBank)
field8 = 'biocategory' # name of the column containing the level 1 biological category as obtained from RVDB_characterization.py (E:\\UPDATE_SCRIPTS_LOGS)
fieldnames=[field2,field3,field4,field5,field6,field7,field8]
fieldtypes=['TEXT']*7
alter_db(dbname,maintable,fieldnames,fieldtypes)

### Populates the RVDB sqlite database with values, taking rvdb_filename as an input (typically U-RVBv$currentvs.fasta, but can also be C-RVDBv$currentvs.fasta)
build_db(wdir,rvdb_filename)


