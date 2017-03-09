GFF3 to EMBL convertion script
==============================
**GFF3+FASTA => EMBL format**

Software to convert GFF3 and fasta to legal EMBL format suitable for [ENA](http://www.ebi.ac.uk/ena) submission.

Based on documentation from http://www.insdc.org/files/feature_table.html, http://www.ebi.ac.uk/ena/WebFeat/ and
ftp://ftp.ebi.ac.uk/pub/databases/embl/doc/usrman.txt.

The output can be validated using the ENA flat file validator distributed by EMBL. Please visit http://www.ebi.ac.uk/ena/software/flat-file-validator and/or https://github.com/enasequence/sequencetools for more information.

## INDEX

[Version](#version)</br>
[Prerequisite](#prerequisite)

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

A correct **gff3 file** and the **genome in fasta format** that has been used to produce the gff file are the mandatory input files.
Then, in order to get a valid EMBL flat file suitable for submission you have to fill carefully all mandatory metadata.

**/!\ Please be aware that a __project ID__ and an __accession number__ are mandatory for a submission to [ENA](http://www.ebi.ac.uk/ena). You don't need this information if you don't plan to submit the data. If you don't have yet this information you can add it later by replacing the corresponding fields. Please visit the [EMBL web site](http://www.ebi.ac.uk/ena/support/genome-submission-faq) to learn how to obtain a __project ID__ and an __accession number__.**

Test data from the Drosophila melanogaster species are located in the example folder.

### Simple case (Suitable for common submissions)

Executing:

 >./GFF2EMBL.py example/dmel_chr4.gff3 example/dmel_chr4.fa
 
Will prompt you to fill one by one the mandatory information needed to produce a proper EMBL file.
Most of time a default value is suggested. Once the software has all the information it need, it will process the input files and will print the result to STDOUT.
 
In order to write the result in the desired file use the **-o** option:
 
 >./GFF2EMBL.py example/dmel_chr4.gff3 example/dmel_chr4.fa -o result.embl

### Complete case (Minimum requirement to launch the software and avoid any prompt)

 >./GFF2EMBL.py example/dmel_chr4.gff3 example/dmel_chr4.fa --data_class STD --topology linear --molecule_type "genomic DNA" --table 1  --taxonomy INV --project_id PRJXXXX -o result.embl

### Advanced case (When you want add more information than those mandatory: e.g publication)

 >./GFF2EMBL.py example/dmel_chr4.gff3 example/dmel_chr4.fa --data_class STD --topology linear --molecule_type "genomic DNA" --table 1  --taxonomy INV --project_id Unknown --author 'author for the reference' --rt 'reference title' --rl 'Some journal' -o result.embl

### Use through a bash script

You may prefer to launch the software through a bash script especially when you want to fill many information, so we provide an example of such script here **example/script.sh**.

In order to use it move in the example folder, then launch the script:
./script.sh

## PARAMETER

Some parameters are mandatory and some others are not. Here is a list of all parameters available. 
You can also find a comprehensive help about the different parameters using the software help command.

positional arguments:
  gff_file              input gff-file
  fasta                 input fasta sequence
  
Arguments related to the EMBL format to check carrefully:

  - -p , --project_id     Project ID. The defalut value is **Unknown** *(This option is used to set up the PR line.)*
  - -r , --table          Translation table. The defalut value is **1** *(This option is used to set up the translation table qualifier transl_table of the CDS features.)* Please visit this [NCBI genetic code] (https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi) page for more information.
  - -s , --species        Sample Species, formatted as 'Genus species (english name)'. Default value: **Genus species (english name)**. This option is used to set up the OS line.
  - -t , --topology       Sequence topology. The default value is **linear** *(This option is used to set up the Topology that is the 3th token of the ID line.)*
  - -d , --data_class     Data class of the sample. The default value is **STD** *(This option is used to set up the 5th token of the ID line.)*
  - -m , --molecule_type  Molecule type of the sample. -the default value is **genomic DNA**
  - -a , --accession      Accession number(s) for the entry. Default value: **UNKNOWN** . This option is used to set up the accession number of the AC line and the first token of the ID line as well as the prefix of the locus_tag qualifier.    
  - -x , --taxonomy       Source taxonomy. No default value. This option is used to set the taxonomic division within ID line (6th token).
  
optional arguments related to the software:

  - -h, --help            Show this help message and exit
  - -v, --verbose         increase verbosity
  - -q, --quiet           decrease verbosity
  - --shame               Suppress the shameless plug
  - -z, --gzip            Gzip output file
  - -o , --output         Output filename.

optional arguments related to the EMBL format:

  - -c , --created        Creation time of the original entry. The default value is the **date of the day**.
  - -g , --organelle      Sample organelle. No default value.
  - -k , --keyword        Keywords for the entry. No default value.
  - -l , --classification Organism classification. The default value is **Life** 
  - --rc                  Reference Comment. No default value.
  - --rx                  Reference cross-reference. No default value.
  - --rg                  Reference Group, the working groups/consortia that produced the record. No default value.
  - --ra , --author       Author for the reference. No default value.
  - --rt                  Reference Title. No default value.
  - --rl                  Reference publishing location. No default value.
  - --translate           Include translation in CDS features. Not activated by default.
  - --version             Sequence version number. The default value is **1** 
  
