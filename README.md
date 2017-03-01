GFF3 to EMBL convertion script
==============================
**GFF3+FASTA => EMBL format**
Software to convert GFF3 and fasta to legal EMBL format suitable for ENA submission.

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
 
## USAGE

### FOREWORD


To get a valid EMBL flat file suitable for submission you have to check that all mandatory metadata are correct,  where necessary fill the information needed to be sure that the software is aware about all information needed.

####**/!\/!\/!\ In order to submit an embl file to [ENA](http://www.ebi.ac.uk/ena) you will need a project ID provided by EMBL. Please visit the [EMBL web site](http://www.ebi.ac.uk/ena/support/genome-submission-faq) to learn how to obtain a project ID. The project ID must be provided to the software through the -a or --accession in order to get a valid embl file for submission !**

Please add the  **--project_id PRJXXXX** parameter in any of these cases when you have your EMBL project ID in order to have a correct EMBL file for submission. If you don't have yet this information you can add it later in the PR line. You don't need this information if you don't plan to submit the data.

### Simple case (Suitable for common submissions)

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
 
 >./GFF2EMBL.py example/dmel_chr4.gff3 example/dmel_chr4.fa -o result.embl

When you use the software in its simpliest way, some default values are assumed.
When you must control the parameters ? (i.e Parameter paragraph)

  - When your data class are **part of this list**: Patent (PAT), Expressed Sequence Tag (EST), Genome Survey Sequence (GSS), High Thoughput CDNA sequencing (HTC), High Thoughput Genome sequencing (HTG), Mass Genome Annotation (MGA), Whole Genome Shotgun (WGS), Transcriptome Shotgun AssEMBLy (TSA), Sequence Tagged Site (STS). Otherwise by default we use the Standard class (STD)
  - When the topology of your sequence **is not** linear. You will have to set it to "circular"
  - When the molecule type **is not** "genomic DNA". The possible value are: genomic RNA", "mRNA", "tRNA", "rRNA", "other RNA", "other DNA", "transcribed RNA", "viral cRNA", "unassigned DNA", "unassigned RNA"
  - When your organism do not use the Standard genetic code (table 1). Please visit this [NCBI genetic code] (https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi) page for more information.
  - When you want to add extra information.

### Complete case (Suitable when default values are not adapted)

 >./GFF2EMBL.py example/dmel_chr4.gff3 example/dmel_chr4.fa --data_class STD --topology linear --molecule_type "genomic DNA" --table 1 -o result.embl

### Comprehensive case (When you want to create an EMBL file with all the possible information)

 >./GFF2EMBL.py example/dmel_chr4.gff3 example/dmel_chr4.fa --data_class STD --topology linear --molecule_type "genomic DNA" --table 1 -o result.embl

### Use through a bash script

In order to help its use, especially when you want to fill many optional information, you could write a quick bash script. We provide an example of such script here **example/script.sh**.

In order to use it move in the example folder, then launch the script:
./script.sh

## PARAMETER

The software can work directly from the annotation file in gff3 format and the genome file in fasta format. 
To produce a proper EMBL output file the software actually needs some information, only the _taxonomy_ information will be asked to the user, all others inportant and mandatory information will be filled with default values. If these options don't reflect your data, you must inform the tool using the corresponding options. Here is a list of such options:

  - **--data_class** the default value is **STD** *(This option is used to set up the 5th token of the ID line.)*
  - **--topology** the default value is **linear** *(This option is used to set up the Topology that is the 3th token of the ID line.)*
  - **--molecule_type** the default value is **genomic DNA** *(This option is used to set up the Molecule type that is the 4th token of the ID line.)*
  - **--table** the defalut value is **1** *(This option is used to set up the translation table qualifier transl_table of the CDS features.)*
  - **--project_id** the defalut value is **Unknown** *(This option is used to set up the PR line.)*
  
Some fields of the EMBL output are optional and are no used by default. If you want to fill them, you will have to inform the tool with the corresponding options. Please use the software help to get a comprehensive list of the available options.


