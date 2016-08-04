GFF3 to EMBL convertion script
==============================

script to attempt to convert GFF3 and fasta to legal EMBL format suitable for 
ENA submission.

Based on documentation from http://www.insdc.org/files/feature_table.html, and
ftp://ftp.ebi.ac.uk/pub/databases/embl/doc/usrman.txt.

The output can be validated using the ENA flat file validator "embl-client_10_09_2015.jar". For an up-to-date ENA flat file validator, please visit http://www.ebi.ac.uk/ena/software/flat-file-validator.

# PREREQUISITE

**Python 2.7**, **biopython** and the **bcbio-gff** python packages.

In order to install biopython and bcbio-gff please use the following steps:

**Mac OS X:**

 Intall pip the python package manager:
 >sudo easy_install pip 

 Install biopython using pip:
 >pip install biopython

 Install bcbio-gff using pip:
 >pip install bcbio-gff
