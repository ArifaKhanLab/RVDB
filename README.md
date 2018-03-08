## RVDB - The Reference Virus Database (RVDB).

The RVDB (Reference Virus DataBase) is a collection of all known viral, viral-related, and viral-like sequences, except bacterial viruses. The RVDB is updated directly from the NCBI RefSeq, NCBI GenBank and NCBI TPA sequences collections. The files contained in this directory contain instructions for performing an update of the RVDB, as well as characterizing the resulting sequences by the following categories: exogenous viral, endogenous nonretroviral, endogenous retroviral, LTR-retrotransposon, and unassigned viral genes / fragments. Further instructions for performing the RVDB update are provided in the Supplementary Methods S1 document, while instructions for characterization are provided in the Supplementary Methods S2 document. The “RVDB_characterization.py” python script is used for the characterization; all other python scripts are used for updating RVDB. The manuscript describing RVDB is currently under review (Norman Goodacre, Aisha Aljanahi, Subhiksha Nandakumar, Mike Mikailov, Arifa S. Khan. A reference viral database (RVDB) to enhance bioinformatics analysis of high throughput sequencing for novel virus detection).


## Instructions for Characterizing RVDB Sequences.docx

This contains instructions for characterizing sequences in an RVDB .fasta file, using the RVDB_characterization.py script. 

## RVDB_characterization.py

Divides an input RVDB .fasta file's sequences into 5 groups: exogenous viral, endogenous nonretroviral, endogenous retroviral, LTR-retrotransposon, and viral gene/fragment. 

## Instructions for Updating RVDB.docx

This contains instructions for running the semi-automated RVDB update pipeline. This pipeline consists of 3 main parts: 1) downloading raw entries from GenBank, RefSeq, and TPA ftp sites, 2) parsing the files, assimilating meta-data, and running filters and semantic screens to select viral, viral-related, and viral-like sequences and filter out host cellular and phage sequences, 3) manual review, clustering, and database (i.e. .fasta) file generation. 

## parse_raw_refseq_PIPE.py 

Takes the two RefSeq viral files and outputs a eukaryotic viral fasta file formatted with two lines per entry (header and sequences), as well as a phage file (same format. Only the eukaryotic viral file will be incorporated into the RVDB update. 

## multiple_gzunzip_PIPE.py

Unzips the two GenBank-format flat files (.gbff) for refseq viral and combines them into a single file. 

## fileops_PIPE.py 

Splits the combined GenBank-format flat file into multiple files, so that each can be read into Python efficiently. 

## rs_acc_mapping_PIPE.py 

Using the GenBank-format flat file metadata for RefSeq viral, finds the duplicate entries’ accessions (original entries, upon which RefSeq viral entries were based. During the unzipping and reformatting of GenBank division files, these duplicates are not included. Also identifies which sequences are RefSeq viral neighbors, so that this can be incorporated into the headers of sequences from GenBank division files as they are unzipped. 

## VDBunzip_reformat_gb_to_fasta_PIPE.py 

Unzips the GenBank division files, labels sequences that are RefSeq viral neighbors during the unzipping. Please note that a modified form of the GenBank Scanner.py script (found in Biopython, typically in Python sub-folder: Lib/site-packages/Bio/GenBank) should be used, in place of the original. This is described below (last .py file). 

## VDBunzip_tpa_PIPE.py 

Unzips the TPA files. TPA files are simply zipped .fasta files. TPA files do not have any associated metadata, so the Scanner.py script is not needed to reformat them. 

## VDBupdate_checkpoint2_PIPE.py 

Runs checkpoint2, which verifies that all GenBank division files were unzipped and reformatted as .fasta files successfully, using the most recent official GenBank release notes as a reference (see "Instructions for Updating RVDB.docx" for information on where to access these release notes.

## SEM-R_PIPE.py 

Semantic screening for the RVDB. Contains three parts: a positive semantic screen to pull in potentially-viral sequences, a size/mirna screen to remove mirnas and sequences <50 b.p. in length, and a negative semantic screen to remove all false-positive sequences pulled in by the positive screen. 

## Scanner.py
This is a modified version of the Scanner.py script that can simply be copied and pasted into the same directory as the original Scanner.py script (overwrite previous). To avoid having to change additional lines in calling scripts to accommodate a different name for Scanner.py, the name of Scanner.py was not changed. The modifications should not hinder any existing functionality of the script, so it can safely be used in place of the original. The modified Scanner.py contains a try/except block in two places, to correct for an error that occurred in a small number of entries, related to the ‘Structured Comments’ metadata. In a very small number of cases, it was noted that the standard Scanner.py script tries to extract the ‘Structured Comments’ metadata, when this metadata is in fact not present. 

## create_U-RVDB.py

Takes all .fasta files containing sequences that are either RefSeq eukaryotic viral or GenBank/TPA that have passed all three parts of the SEM-R screen, and combines them into one file - the unclustered RVDB, named U-RVDBv$version.fasta. 

## create_C-RVDB.py

Takes the cluster repressentatives from the CD-HIT-EST clustering output (.clstr file) and uses them as a filter to extract a subset of the unclustered RVDB, resulting in the clustered RVDB, named C-RVDBv$version.fasta. 