##################################################################################################
########## Beginning of key[words,regexs,rules]                                        ###########
##################################################################################################

import re

##################################################################################################
########## Positive keywords, positive regular expressions.                            ###########
########## Keywords < 5 characters in length receive flanking ' ' (single whitespace)  ###########
########## on either side, with two exceptions: 'virus' and 'viral'.                   ###########
##################################################################################################
poskws=[' copia ',' env ',' erv ',' gag ',' gipsy ',' herv ','harlequin',' iap ',' ltr ',' mdg3 ',' mdg1 ',
        ' merv ',' pol ',' rt gene ',' sire ',' vlp ','blastopia','capsid',' delta element','dutpase',
        'endogenous','envelop','erranti','gypsy element','gypsy like','gypsy type',' iap ltr ','insertion element',
        'integrase','intracisternal a particle',' long terminal repeat',' ltr ','micropia','polycomb response element',
        'polyprotein','provir','replicase','retrovir','retroelement','retrotranspos','reverse transc','spuma',' ty element',
        ' ty insertion','viral','virales','viridae','virion','viroid','virus']

########## Matches HERVs of the following form: HERV, HERVE, HERV-I,  HERV I, HERV9, HERV 9, HERV-FC
pos_herv_finder=re.compile('\s{1}herv\s{0,1}[a-z0-9]{0,2}\s{1}')
########## Matches HMLs of the following form: hml1, hml 2
pos_hml_finder=re.compile('\s{1}hml\s{0,1}[0-9]{1}\s{1}')
########## Matches Ty elements of the following form: ty1, ty 2
pos_ty_finder=re.compile('\s{1}ty\s{0,1}[0-9]{1}\s{1}')
def positive_rules(description_s):
    ds=description_s
    viral=[]
    if 'transpos' in ds:
        if ' bel ' in ds:
            viral.append('transpos + bel')
        if ' pao ' in ds:
            viral.append('transpos + pao')
        if ' roo ' in ds:
            viral.append('transpos + roo')
        if 'morgane' in ds:
            viral.append('transpos + morgane')
    return viral

##################################################################################################
########## Negative keywords, positive regular expressions.                            ###########
########## Keywords < 5 characters in length receive flanking ' ' (single whitespace). ###########
##################################################################################################
negkws=[' 12s ', ' 16s ', ' 18s ', ' 23s ',' 26s ', ' 28s ', ' alu ', 'mammalian genes involved in viral infection',
        ' method ',' methods ',' moab ','oncogene',' p sine',' phage ',' microphage ','pseudogene','receptor',
        ' rip tad ',' rt pcr ',' rtpcr ',' scfv ',' sine ','anaphagefferens transposon insertion flanking sequence',
        'antibod','antisense iap','antivir','apolipoprotein','apoptosis inhibitor','archae','artificial sequences',
        'atpase','bacteri',' bark of madame vinous',' beta lactamase','binding protein','caudovirales','cellular',
        'chloroplast','choretroides','composition and methods','compositions and assays',
        ' elk retro ',' entry protein','enveloping','express immune library','flanking',' foxy transposable element',
        'gammaretroviral integration into nucleosomal target dna ',' iap associated factor',' iap like protein',' iap nucleobase',
        'immunoglobulin','inclusion protein','inhibition of therapy','inhibitor of apoptosis','integration 1 site',
        'integration junction','integration site','integron','interacting protein','interactive','interactor',
        'intergenic spacer','inverted repeat','involved in phenomena of','leishmania repetitive element','macrophage',
        'markers for detecting genetic polymorphism','method for evaluating','microviridae','mitochondri',
        ' motor protein like ',' mouse retroviral tagged','myovirid','neutralizing influenza virus','non ltr','nonltr',
        'oncogene','oncoprotein',' phage ','phagemid','phgmid','piggybac',' poly a specific ribonuclease','polyadenylation',
        'preintegration site','primer walking',' probe set ','proviral insertion site','pseudogene',' pspr rt mr ',
        'putative condensation factor','receptor','region flanking','retroviral mutagenesis library',
        'rhamnosyltransferase','ribosom',' rs3 rt ',' rt raphanus sativus cdna 5 mrna sequence',' sire guyss ',
        'tandem repeats and reverse transcriptase primer binding site ','taphrinaceae','telomerase','tobacco mosaic virus induced library',
        'ubiquitin specific peptidase','vector',' viral and viral associated mirnas and uses thereof',
        'viral genomic survey of stool from south east asian children with acute flaccid paralysis',
        'metagenome genomic survey sequence','metagenome genomic genomic survey sequence',' virus resistance',
        'antifreeze','associating',' au element','boolean','cheilotoma musciformis','coatomer','cytochrome',
        'endogenous cellulase','enhancer factor',' entry mediator','envelope membrane','fertilization envelope',
        ' gut membrane polyprotein','harbinger','herves','histone polyprotein','hivep2',' iap antagonist','insertion site',
        'instability elements','interacting','junction',' line retro ',' line 1 ',' line1 ',' line like ',' line type ',
        ' ltr centromeric satellite',' ltr race ','microphage',' non autonomous','nuclear envelope',' outer envelope',
        'pegasus','plastid envelope','recombination protein','recombination site ','resistance','retrotransposon line ',' seed library','suppressor of ty',
        ' taco ','tigger',' tip100 ',' tn env ',' trna gag ',' trna leu ','ubiquitin polyprotein',' viral associated','viralicidic',' virus activating protease',
        ' virus transformed','vitelline envelope',' line element','escherichia','phosphate starvation inducible protein',
        ' t4 like major capsid protein',' solo ltr ','spumarius','subtilisin like','retroposed','deleted element',
        'pheidole harlequina','integrin associated protein','intestinal alkaline phosphatase','vitalline coat protein',
        'islet amyloid polypeptide precursor',
        'satellite','integrin associated protein',' ltr 1 interferon regulatory factor',' pol rfamide','plasmid replicase',
        " 3' flank "," 5' flank ",' pol cytoplasmic male sterility fragment',' late cornified envelope',' mlgapc ',' g3pdh ','transposon etn ','retrotransposon line',
        'retroposed','philaenus spumarius','atelopus spumarius','isolate pol ','isolate ltr ','human iap homolog','glucagon polyprotein','m gag 384',
        'polyprotein encoding gonadotropin releasing hormone','lysidyl aminotransferase','cyclophilin a variant iap ','isolate 208 pol ','translocation breakpoint region',
        'intestinal alkaline phosphatase','islet amyloid polypeptide','allobates spumaponens','synthetic linker','m pol ex 1169','sftpb jm ltr',
        'ribulosebisphosphate carboxylase form II polyprotein precursor','euglossa retroviridis','diabrotica virgifera virgifera endogenous endoglucanase',
        'clone hhrii pol','endogenous vascular elastase',' shell coat protein',' ias virus ',' egg envelope',' t4 like ',' t7 like ',' t4like',' t7like',
        'integration region','virus susceptibility','uncultured marine microorganism','caenorhabditis elegans chromosome unknown clone','viral responsive protein',
        ' iap binding motif ',' pcr product',' anti viral ',' env capsular antigen protein',' ecotype ty ',' zea mays cultivar ty','cornified envelope',
        'preintegration sequence','virion associated','endogenous antisense rna ','crossover point ','intetration site ','chicken dna that hybridizes']



########## Matches pol genes of the form: pol i, polii, pol a, pol 1, pol alpha
neg_pol_finder=re.compile('\s{1}pol\s{0,1}[a,b,d,e,h,i,k,1,2]\s{1}|\s{1}pol\s{0,1}ii\s{1}|\s{1}pol\s{0,1}iii\s{1}\s{1}pol\s{0,1}alpha\s{1}|\s{1}pol\s{0,1}beta\s{1}|\s{1}pol\s{0,1}delta\s{1}|\s{1}pol\s{0,1}epsilon\s{1}|\s{1}pol\s{0,1}eta\s{1}|\s{1}pol\s{0,1}ii\s{1}|\s{1}polymerase\s{0,1}[a,b,d,e,h,i,k,1,2]\s{1}|\s{1}polymerase\s{0,1}alpha\s{1}|\s{1}polymerase\s{0,1}beta\s{1}|\s{1}polymerase\s{0,1}delta\s{1}|\s{1}polymerase\s{0,1}epsilon\s{1}|\s{1}polymerase\s{0,1}eta\s{1}|\s{1}polymerase\s{0,1}ii\s{1}|\s{1}pol\s{0,1}mu\s{1}|\s{1}polymerase\s{0,1}mu\s{1}|\s{1}pol\s{0,1}kappa\s{1}|\s{1}polymerase\s{0,1}kappa\s{1}|\s{1}pol\s{0,1}[i,ii,iii][a,b,c,d,e,f,g,h,i,k]\s{1}|\s{1}polymerase\s{0,1}[i,ii,iii][a,b,c,d,e,f,g,h,i,k]\s{1}')
neg_bacter_finder=re.compile('\s[a-z]+bacter\s')
neg_contig_finder=re.compile('\scontig[0-9]*\s')
neg_unculturedvirus_finder=re.compile('\suncultured\svirus\sehp\s[a-z]{1}[0-9]{1}[a-z]{1}\s')
neg_isolatety_finder=re.compile('\sisolate\sty[1-9]*\s')
########## All of the specific combinations of words that serve as negative criteria (i.e. flagging a sequence as non-viral).
########## These negative rules serve as highly-specific neg selection criteria when the true and false positives are biologically
########## very similar (e.x. flagging DNA transposons from retrotransposons).
########## Please note that the minimum length for not requiring a flanking whitespace in a keyword is 3, rather than 5, for the secondary conditionals
########## , i.e. those additional keywords to 'tranpos'
def negative_rules(description_s):
    ds=description_s
    nonviral=[]
    ambiguous=False
    if 'virus' in ds or 'viral' in ds:
        if 'challenged' in ds or 'induced' in ds or 'injected' in ds:
            nonviral.append('virus/viral + challenged/induced/injected')
            ambiguous=True
    if 'infection' in ds or 'infected' in ds:
        nonviral.append('infected/infection')
    if 'transposon insertion' in ds and 'sequence tagged site' in ds:
        nonviral.append('transposon insertion + sequence tagged site')
    if 'phage' in ds:
        if not ' ogre ' in ds and not 'anaphagefferens' in ds:
            nonviral.append('phage - ogre/anaphagefferens')
    if 'transpos' in ds:
        if 'mariner' in ds:
            nonviral.append('transpos + mariner')
        if 'pogo' in ds:
            nonviral.append('transpos + pogo')
        if ' tgm ' in ds:
            nonviral.append('transpos + tgm')
        if ' fot ' in ds:
            nonviral.append('transpos + fot')
        if 'abr1' in ds:
            nonviral.append('transpos + abr1')
        if 'idle' in ds:
            nonviral.append('transpos + idle')
        if ' tam ' in ds:
            nonviral.append('transpos + tam')
        if ' hat ' in ds:
            nonviral.append('transpos + hat')
        if 'helitron' in ds:
            nonviral.append('transpos + helitron')
        if ' tc1 ' in ds or ' tc 1 ' in ds:
            nonviral.append('transpos + tc1/tc 1')
        if 'cacta' in ds:
            nonviral.append('transpos + cacta')
        if 'krak' in ds:
            nonviral.append('transpos + krak')
        if ' sx ' in ds:
            nonviral.append('transpos + sx')
        if ' hop ' in ds:
            nonviral.append('transpos + hop')
        if 'mutator' in ds and 'like' in ds:
            nonviral.append('transpos + mutator + like')
        if ' tart ' in ds:
            nonviral.append('transpos + tart')
        if 'hsmar' in ds:
            nonviral.append('transpos + hsmar')
        if 'talisker' in ds:
            nonviral.append('transpos + talisker')
        if 'mite jura' in ds:
            nonviral.append('transpos + mite jura')
        if 'mite arran' in ds:
            nonviral.append('transpos + mite arran')
        if 'tpn1' in ds:
            nonviral.append('transpos + tpn1')
        if ' l1p ' in ds and ' ma2 ' in ds:
            nonviral.append('transpos + l1p + ma2')
        if ' pif ' in ds and 'like' in ds:
            nonviral.append('transpos + pif + like')
        if 'tourist' in ds and ' zm ' in ds:
            nonviral.append('transpos + tourist + zm')
        if 'scooter' in ds:
            nonviral.append('transpos + scooter')
        if 'tat1' in ds or 'tat 1' in ds:
            nonviral.append('transpos + tat1/tat 1')
        if 'caspar' in ds:
            nonviral.append('transpos + caspar')
        if ' tip ' in ds:
            nonviral.append('transpos + tip')
        if ' ac ' in ds:
            if ' ds ' in ds or 'element' in ds:
                nonviral.append('transpos + ac + ds/element')
        if 'activ' in ds:
            nonviral.append('transpos + activ')
        if 'mudr' in ds:
            nonviral.append('transpos + mudr')
        if 'au element' in ds:
            nonviral.append('transpos + au element')
        if 'bilbo' in ds:
            nonviral.append('transpos + bilbo')
        if 'revolver' in ds:
            nonviral.append('transpos + revolver')
        if 'pong like' in ds:
            nonviral.append('transpos + pong like')
        if 'minos' in ds:
            nonviral.append('transpos + minos')
        if ' cr1 ' in ds or ' cr 1 ' in ds:
            nonviral.append('transpos + cr1/cr 1')
        if 'collect2' in ds:
            nonviral.append('transpos + collect2')
        if 'tpo1' in ds:
            nonviral.append('transpos + tpo1')
        if 'albatross' in ds:
            nonviral.append('transpos + albatgross')
        if ' tgm ' in ds:
            nonviral.append('transpos + tgm')
        if 'en 1' in ds:
            nonviral.append('transpos + en 1')
        if 'tam3' in ds:
            nonviral.append('transpos + tam3')
        if 'pptn' in ds:
            nonviral.append('transpos + pptn')
        if ' tn ' in ds:
            nonviral.append('transpos + tn')
        if 'ping' in ds:
            nonviral.append('transpos + ping')
        if 'line1' in ds or 'line 1' in ds or 'line like' in ds or 'line type' in ds or ' l1 ' in ds:
            nonviral.append('transpos + line1/line 1/line like/ line type/l1')
        if ' sine' in ds:
            nonviral.append('transpos + sine')
    if 'retroelement' in ds or 'retrotranspos' in ds or 'reverse transc' in ds:
        if 'line1' in ds or 'line 1' in ds or 'line like' in ds or 'line type' in ds or ' l1 ' in ds:
            nonviral.append('retroelement/retrotranspos/reverse transc + line1/line 1/line like/ line type/l1')
    if 'baculovir' in ds and ' iap ' in ds and 'repeat' in ds:
        nonviral.append('baculovir + iap + repeat')
    if ' clc ' in ds and ' ltr ' in ds and ' lactuca sativa ' in ds:
        nonviral.append('clc + ltr + lactuca sativa')
    if ' fst ' in ds:
        if ' ty 1 ' in ds or ' ty1 ' in ds:
            nonviral.append('fst + ty1 / ty 1')
    if ' bos taurus' in ds and ' v region' in ds:
        nonviral.append('bos taurus + v region')
    if 'rt gene' in ds and not ' rt gene' in ds:
        nonviral.append('rt gene')
    if ' fst ' in ds:
        if ' ty1 ' in ds or ' ty 1 ' in ds:
            nonviral.append('fst + ty1/ty 1')
    if 'rt gene' in ds:
        if not ' rt gene' in ds:
            nonviral.append('rt gene')
    if ' etn ' in ds:
        if ' ltr ' in ds:
            nonviral.append('etn + ltr')
    if ' coat protein' in ds or ' core protein' in ds:
        if not 'virus' in ds and not 'viral' in ds:
            nonviral.append('coat/core protein - virus/viral')
    if ' ig ' in ds:
        if ' heavy chain ' in ds or ' light chain ' in ds:
            nonviral.append('ig + heavy/light chain')
    if ' l1 ' in ds:
        if not ' l1 gene ' in ds:
            nonviral.append('l1')
    if 'envelope protein' in ds:
        if not 'virus' in ds and not 'viral' in ds:
            if 'envelope protein odv' in ds:
                nonviral.append('envelope protein odv')
    if ' ty1 ' in ds or ' ty 1 ' in ds or ' ty2 ' in ds or ' ty 2 ' in ds or ' ty3 ' in ds or ' ty 3 ' in ds:
        if 'tabanus yao' in ds:
            nonviral.append('tabanus yap ty')
    if ' env ' in ds:
        if 'uncultured trypanosome' in ds:
            nonviral.append('uncultured trypanosome isolate *-env')
    if ' pol ' in ds:
        if 'mycoplasma synoviae' in ds:
            nonviral.append('mycoplasma synoviae isolate */POL/*')
        if 'aurelia sp 1 sensu dawson' in ds:
            nonviral.append('aurelia sp 1 sensu dawson, pol')
    if ' ty2 ' in ds or ' ty3 ' in ds  or ' ty4 ' in ds:
        if ' gs gene ' in ds or ' nak gene ' in ds or ' pgd gene ' in ds or ' idh gene ' in ds or ' cad gene ' in ds or ' gln gene ' in ds or ' tpi gene ' in ds:
            nonviral.append('ty2/ty3/ty4 + gs/nak/pgd/idh/cad/gln/tpi')
    if 'pacifastacus gambelii isolate' in ds:
        if ' pol ' in ds:
            nonviral.append('pacifastacus gambelii isolate + pol ')
    return [nonviral,ambiguous]
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

##################################################################################################
########## End of key[words,regexs,rules]                                              ###########
##################################################################################################

##################################################################################################
##################################################################################################

##################################################################################################
########## Beginning of execute code                                                   ###########
##################################################################################################

def run_screen(screentype,date,istpa,dupaccs,neighaccs):
    if screentype=='poskw':
        fastafilenames=[]
        if not istpa:
            for fn in os.listdir(gbdir):
                if fn.startswith('gb') and fn.endswith('.seq.'+date):
                    fastafilenames.append(fn)
        else:
            for fn in os.listdir(tpadir):
                if 'tpa_cu' in fn and fn.endswith('.seq.'+date):
                    fastafilenames.append(fn)
        run_pos_screen_allgb(date,fastafilenames,istpa,neighaccs)
    if screentype=='sizemirna':
        pscreen_filenames=[]
        if not istpa:
            for fn in os.listdir(poskw_outdir):
                if fn.startswith('gb') and fn.endswith('seq.'+date+'.pscreen.fasta'):
                    pscreen_filenames.append(fn)
        else:
            for fn in os.listdir(tpadir):
                if 'tpa_cu' in fn and fn.endswith('seq.'+date+'.pscreen.fasta'):
                    fastafilenames.append(fn)
        run_size_mirna_screen(date,pscreen_filenames,istpa)
    if screentype=='negkw':
        sizemirna_filenames=[]
        if not istpa:
            for fn in os.listdir(sizemirna_outdir):
                if fn.startswith('gb') and fn.endswith('seq.'+date+'.sizemirnascreenOK.fasta'):
                    sizemirna_filenames.append(fn)
        else:
            for fn in os.listdir(tpadir):
                if 'tpa_cu' in fn and fn.endswith('seq.'+date+'.sizemirnascreenOK.fasta'):
                    fastafilenames.append(fn)
        run_neg_screen(date,sizemirna_filenames,istpa,dupaccs)

def run_pos_screen_allgb(date,fasta_filenames,istpa,neighaccs):
    outheadersfn='headers_pscreen.txt'
    outheadersfile=open(poskw_outdir+'\\'+outheadersfn,'w')
    logf=open(logdir+'\\'+'poskw_screen_log.txt','w')
    poskw_outfns=[]
    for fn in fasta_filenames:
        fn=fn.strip()
        inf=open(gbdir+'\\'+fn)
        poskw_outfn=fn.replace('seq.'+date,'seq.'+date+'.pscreen.fasta')
        poskw_outf=open(poskw_outdir+'\\'+poskw_outfn,'w')
        for i,line in enumerate(inf):
            if line.startswith('>acc'):
                if istpa:
                    description=line.split('|')[3]
                else:
                    description=' '.join([subdesc.strip() for subdesc in line.split('|')[3:-2]])
                ds=simplify_description(description)
                viral,d_poskws=pos_screen(ds)
                if '|VRL|' in line:
                    viral=True
                if viral:
                    outheadersfile.write(line.strip()+'@'+', '.join(d_poskws)+'\n')
            if viral:
                poskw_outf.write(line)
        poskw_outf.close()
        poskw_outfns.append(poskw_outfn)
    for poskw_outfn in poskw_outfns:
        logf.write(poskw_outfn+'\n')
        print poskw_outfn
    logf.close()
    outheadersfile.close()

def run_size_mirna_screen(date,pscreen_filenames,istpa):
    outheadersfn1='headers_sizemirnaOK.txt'
    outheadersfn2='headers_sizemirnaFLAG.txt'
    outheadersfile1=open(sizemirna_outdir+'\\'+outheadersfn1,'w')
    outheadersfile2=open(sizemirna_outdir+'\\'+outheadersfn2,'w')
    logf=open(logdir+'\\'+'sizemirna_screen_log.txt','w')
    sizemirna_outfns=[]
    for fn in pscreen_filenames:
        fn=fn.strip()
        inf=open(poskw_outdir+'\\'+fn)
        entries=inf.read().strip().split('>acc')[1:]
        size_mirna_screen_fn1=fn.replace('.pscreen.fasta','.sizemirnascreenOK.fasta')
        sizemirna_outfns.append(size_mirna_screen_fn1)
        size_mirna_outf1=open(sizemirna_outdir+'\\'+size_mirna_screen_fn1,'w')
        size_mirna_screen_fn2=fn.replace('.pscreen.fasta','.sizemirnascreenFLAG.fasta')
        size_mirna_outf2=open(sizemirna_outdir+'\\'+size_mirna_screen_fn2,'w')
        for entry in entries:
            fail=False
            header,seq=entry.strip().split('\n')
            if istpa:
                description=header.split('|')[3].strip()
            else:
                description=' '.join([subdesc.strip() for subdesc in header.split('|')[3:-2]])
            sdesc=simplify_description(description)
            if ' mirna ' in sdesc or 'microrna' in sdesc or 'micro rna' in sdesc:
                fail=True
            if len(seq.strip())<=50:
                fail=True
            if not fail:
                outheadersfile1.write(header+'\n')
                size_mirna_outf1.write('>acc'+entry)
            else:
                outheadersfile2.write(header+'\n')
                size_mirna_outf2.write('>acc'+entry)
        size_mirna_outf1.close()
        size_mirna_outf2.close()
    for sizemirna_outfn in sizemirna_outfns:
        logf.write(sizemirna_outfn+'\n')
        print sizemirna_outfn
    logf.close()
    outheadersfile1.close()
    outheadersfile2.close()

def run_neg_screen(date,sizemirna_filenames,istpa,dupaccs):
##    print 'Flagging all remaining entries that contain non-viral kws/regexes/rules as non-viral, from update: '+currentvs
    outheadersfn1='headers_nscreenVRL_'+currentvs+'.txt'
    outheadersfn2='headers_nscreenOK_'+currentvs+'.txt'
    outheadersfn3='headers_nscreenFLAG_'+currentvs+'.txt'
    outheadersfn4='headers_nscreenAMB_'+currentvs+'.txt'
    outheadersfn5='headers_nscreenDUP_'+currentvs+'.txt'
    outheadersfile1=open(negkw_outdir+'\\'+outheadersfn1,'w')
    outheadersfile2=open(negkw_outdir+'\\'+outheadersfn2,'w')
    outheadersfile3=open(negkw_outdir+'\\'+outheadersfn3,'w')
    outheadersfile4=open(negkw_outdir+'\\'+outheadersfn4,'w')
    outheadersfile5=open(negkw_outdir+'\\'+outheadersfn5,'w')
    logf=open(logdir+'\\'+'negkw_screen_log.txt','w')
    for fn in sizemirna_filenames:
        fn=fn.strip()
        inf=open(sizemirna_outdir+'\\'+fn)
        entries=inf.read().strip().split('>acc')[1:]
        if 'gbvrl' in fn:
            vrl=True
        else:
            vrl=False
        nscreen_fn1=fn.replace('.sizemirnascreenOK.fasta','.nscreenVRL_.fasta')
        nscreen_f1=open(negkw_outdir+'\\'+nscreen_fn1,'w')
        nscreen_fn2=fn.replace('.sizemirnascreenOK.fasta','.nscreenOK_.fasta')
        nscreen_f2=open(negkw_outdir+'\\'+nscreen_fn2,'w')
        nscreen_fn3=fn.replace('.sizemirnascreenOK.fasta','.nscreenFLAG_.fasta')
        nscreen_f3=open(negkw_outdir+'\\'+nscreen_fn3,'w')
        nscreen_fn4=fn.replace('.sizemirnascreenOK.fasta','.nscreenAMB_.fasta')
        nscreen_f4=open(negkw_outdir+'\\'+nscreen_fn4,'w')
        nscreen_fn5=fn.replace('.sizemirnascreenOK.fasta','.nscreenDUP_.fasta')
        nscreen_f5=open(negkw_outdir+'\\'+nscreen_fn5,'w')
        logf.write(nscreen_fn1+'\n')
        logf.write(nscreen_fn2+'\n')
        logf.write(nscreen_fn3+'\n')
        logf.write(nscreen_fn4+'\n')
        logf.write(nscreen_fn5+'\n')
        for entry in entries:
            entry='>acc'+entry
            viral=True
            header,seq=entry.strip().split('\n')
            acc=header.split('|')[2]
            if istpa:
                description=header.split('|')[3].strip()
            else:
                description=' '.join([subdesc.strip() for subdesc in header.split('|')[3:-2]])
            sdesc=simplify_description(description)
            viral,d_negkws,ambiguous=neg_screen(sdesc)
            if vrl:
                outheadersfile1.write(header+'@'+', '.join(d_negkws)+'\n')
                nscreen_f1.write(entry)
            else:
                if viral:
                    outheadersfile2.write(header+'@'+', '.join(d_negkws)+'\n')
                    nscreen_f2.write(entry)
                else:
                    if ambiguous:
                        outheadersfile4.write(header+'@'+', '.join(d_negkws)+'\n')
                        nscreen_f4.write(entry)
                    else:
                        if acc.split('.')[0] in dupaccs:
                            outheadersfile5.write(header+'@'+', '.join(d_negkws)+'\n')
                            nscreen_f5.write(entry)
                        else:
                            outheadersfile3.write(header+'@'+', '.join(d_negkws)+'\n')
                            nscreen_f3.write(entry)
        nscreen_f1.close()
        nscreen_f2.close()
        nscreen_f3.close()
        nscreen_f4.close()
    logf.close()
    outheadersfile1.close()
    outheadersfile2.close()
    outheadersfile3.close()
        
def pos_screen(description_s):
    viral=False
    ds=description_s
    d_poskws=[]
    d_negkws=[]
    for poskw in poskws:
        if poskw in ds:
            d_poskws.append(poskw)
    hervs=pos_herv_finder.findall(ds)
    hmls=pos_hml_finder.findall(ds)
    tys=pos_ty_finder.findall(ds)
    rules_viral=positive_rules(ds)
    if len(hervs)>0:
        hervs=list(set(hervs))
        d_poskws.extend(hervs)
    if len(hmls)>0:
        hmls=list(set(hmls))
        d_poskws.extend(hmls)
    if len(tys)>0:
        tys=list(set(tys))
        d_poskws.extend(tys)
    if len(rules_viral)>0:
        rules_viral=list(set(rules_viral))
        d_poskws.extend(rules_viral)
    if len(d_poskws)>0:
        viral=True
        d_poskws=sorted(d_poskws)
    return [viral,d_poskws]

def neg_screen(description_s):
    viral=True
    ds=description_s
    d_negkws=[]
    ambiguous=False
    for negkw in negkws:
        if negkw in ds:
            d_negkws.append(negkw)
    pols=neg_pol_finder.findall(ds)
    if len(pols)>0:
        pols=list(set(pols))
        d_negkws.extend(pols)
    bacters=neg_bacter_finder.findall(ds)
    if len(bacters)>0:
        bacters=list(set(bacters))
        d_negkws.extend(bacters)
    contigs=neg_contig_finder.findall(ds)
    if len(contigs)>0:
        contigs=list(set(contigs))
        d_negkws.extend(contigs)
    uncultured_viruses=neg_unculturedvirus_finder.findall(ds)
    isolate_tys=neg_isolatety_finder.findall(ds)
    if len(isolate_tys)>0:
        isolate_tys=list(set(isolate_tys))
        d_negkws.extend(isolate_tys)
    if len(uncultured_viruses)>0:
        uncultured_viruses=list(set(uncultured_viruses))
        d_negkws.extend(uncultured_viruses)
    d_negkws=sorted(d_negkws)
    nonviralrules,ambiguous=negative_rules(ds)
    if len(nonviralrules)>0:
        nonviralrules=sorted(nonviralrules)
        d_negkws=sorted(d_negkws)
        d_negkws.extend(nonviralrules)
    if len(d_negkws)>0:
        viral=False
    return [viral,d_negkws,ambiguous]

##################################################################################################
########## Beginning of execute control commands                                       ###########
##################################################################################################

import os
import sys
homedir=sys.argv[1]
date=sys.argv[2]
currentvs=sys.argv[3]
screentype=sys.argv[4]
seqtype=sys.argv[5]
if seqtype=='tpa':
    istpa=True
else:
    istpa=False
wdir=homedir+'\\RVDBv'+currentvs
gbdir=wdir+'\\GenBank_raw_data_'+date
poskw_outdir=gbdir+'\\poskw_out_'+date
sizemirna_outdir=gbdir+'\\sizemirna_out_'+date
negkw_outdir=gbdir+'\\negkw_out_'+date
scriptdir=gbdir+'\\scripts'
logdir=gbdir+'\\log'
refseqdir=wdir+'\\RefSeq_raw_data_'+date
dupacc_filename=refseqdir+'\\'+'refseq_viral_originalaccs.txt'
dupacc_infile=open(dupacc_filename)
dupaccs=set(dupacc_infile.read().strip().split('\n'))
dupacc_infile.close()
neighacc_filename=refseqdir+'\\'+'neighbor_accs.txt'
neighacc_infile=open(neighacc_filename)
neighaccs=set(neighacc_infile.read().strip().split('\n'))
neighacc_infile.close()
run_screen(screentype,date,istpa,dupaccs,neighaccs)

