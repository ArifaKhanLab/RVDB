# RVDB
A nonredundant, reference viral database (RVDB)
README

We have chosen to apply the Creative Commons Attribution 3.0 
Unsupported License to this version of the software.

Reference Viral Database (RVDB) is developed by Arifa Khan's group at CBER, FDA for enhancing virus
detection using next-generation sequencing (NGS) technologies. U-RVDB is the unclustered database 
and C-RVDB is clustered to reduce redundancy and retain sequence diversity, representing all viral 
sequences regardless of size. The initial version of C-RVDB (10.2) was released on 12/28/2016. The 
initial version of U-RVDBv10.2 was released on 5/5/2017. The release notes for current versions of 
both U-RVDBv and C-RVDBv can be accessed in DOWNLOADS section. Updated versions will be provided as 
they became available. 


Contributors: Norman Goodacre, Subhiksha Nandakumar, and Arifa S. Khan

Questions/Comments:
If you have any questions or comments regarding rVDB please contact Arifa Khan (Arifa.Khan@fda.hhs.gov).



RVDB USAGE
--------------------------------------------------------------------------------------------------

RVDB can be downloaded in fasta format from the DOWNLOADS section of the HIVE page. 
Fasta file of RVDB can be formatted as a Blast searchable database using makeblastdb script from the 
blast suite of tools. 
Instructions are available at: https://www.ncbi.nlm.nih.gov/books/NBK279688/ 
Once formatted, submit queries of interest using BLASTN or TBLASTX.



THE INPUT
---------------------------------------------------------------------------------------------------

Input is a user-supplied sequences of any type in fasta format. RVDB can also 
be queried with large NGS datasets.

RVDB-prot and RVDB-prot-HMM were developed by Thomas Bigot in Marc Eloit’s Pathogen Discovery group in 
collaboration with Center of Bioinformatics, Biostatistics and Integrative Biology (C3BI) at Institut Pasteur, 
for enhancing virus detection using NGS technologies. They are updated after each 
new release of the nucleotidic database. The version number of the proteic databases follows the one of the 
original nucleic database. More information is provided by the link for the proteic databases at 
Institut Pasteur  (http://rvdb-prot.pasteur.fr/).

