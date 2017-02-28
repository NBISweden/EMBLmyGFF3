GFF3 to EMBL convertion script
==============================

Software to convert GFF3 and fasta to legal EMBL format suitable for 
ENA submission.

Based on documentation from http://www.insdc.org/files/feature_table.html, http://www.ebi.ac.uk/ena/WebFeat/ and
ftp://ftp.ebi.ac.uk/pub/databases/embl/doc/usrman.txt.

The output can be validated using the ENA flat file validator distributed by EMBL. Please visit http://www.ebi.ac.uk/ena/software/flat-file-validator and/or https://github.com/enasequence/sequencetools for more information.

## VERSION 
**GFF2EMBL.1.0.0**

This is the first version released (1 March 2017). 

## PREREQUISITE

**Python 2.7**, **biopython** and the **bcbio-gff** python packages.

In order to install biopython and bcbio-gff please use the following steps:

**Mac OS X:**

 Intall pip the python package manager:
 >sudo easy_install pip 

 Install biopython using pip:
 >pip install biopython==1.67

 Install bcbio-gff using pip:
 >pip install bcbio-gff=0.6.4

## Check the installation

 Executing:
 >./GFF2EMBL.py
 
 or
 
 >./GFF2EMBL.py -h
 
will display some help.
 

## Example

### Simple case

A correct **gff3 file** and the **genome in fasta format** that has been used to produce the gff file are the only things mandatory.
Test data from the Drosophila melanogaster species are located in the example folder.

 Executing:
 >./GFF2EMBL.py example/dmel_chr4.gff3 example/dmel_chr4.fa
 
 should prompt you a remind that even if it's not mandatory, some information are always nice to add to an EMBL file.
 To skip that step just press ENTER.
 
 Then as the taxonomy option has not ben filled in this case and it's a mandatory field, a prompt ask you top chose between these 15 values:
  - MUS	Mus musculus
  - PHG	Bacteriophage
  - UNC	Unclassified
  - ENV	Environmental Sample
  - FUN	Fungal
  - VRT	Other Vertebrate
  - HUM	Human
  - INV	Invertebrate
  - TGN	Transgenic
  - ROD	Other Rodent
  - SYN	Synthetic
  - PLN	Plant
  - MAM	Other Mammal
  - PRO	Prokaryote
  - VRL	Viral
 
 type:
 >INV

 and press ENTER.

 That's it ! The result will be printed to STDOUT.
 
 In order to write the result in the desired file use this command:
 
 >./GFF2EMBL.py example/dmel_chr4.gff3 example/dmel_chr4.fa > result.embl

### Use through a bash script

In order to help its use, especially when you want to fill many optional information, you could write a quick bash script. We provide an example of such script here **example/script.sh**.

In order to use it move in the example folder, then launch the script:
./script.sh



