GFF3 to EMBL conversion tool
==============================

<img src="EMBLmyGFF3.png" width="500" height="150" />

Software to convert GFF3 and fasta to legal EMBL format suitable for [ENA](http://www.ebi.ac.uk/ena) submission.

Based on documentation from http://www.insdc.org/files/feature_table.html, http://www.ebi.ac.uk/ena/WebFeat/ and
ftp://ftp.ebi.ac.uk/pub/databases/embl/doc/usrman.txt.

The output can be validated using the ENA flat file validator distributed by EMBL. Please visit http://www.ebi.ac.uk/ena/software/flat-file-validator and/or https://github.com/enasequence/sequencetools for more information.

## INDEX

[Prerequisite](#prerequisite)</br>
[Installation](#installation)</br>
[Usage](#usage)</br>
&nbsp;&nbsp;&nbsp;[Foreword](#foreword)</br>
&nbsp;&nbsp;&nbsp;[Simple case](#simple-case)</br>
&nbsp;&nbsp;&nbsp;[Complete case](#complete-case)</br>
&nbsp;&nbsp;&nbsp;[Advanced case](#advanced-case)</br>
&nbsp;&nbsp;&nbsp;[Use through a bash script](#use-through-a-bash-script)</br>
[Parameter](#parameter)</br>
[Mapping](#mapping)</br>
&nbsp;&nbsp;&nbsp;[Feature type](#feature-type)</br>
&nbsp;&nbsp;&nbsp;[GFF3 Attribute to EMBL qualifier](#gff3-attribute-to-embl-qualifier)</br>
&nbsp;&nbsp;&nbsp;[Other](#other)</br>
[Known issues](#known-issues)</br>
[Citation](#citation)
[Author](#author)

## PREREQUISITE

**Python 2.7**, **biopython 1.67** and the **bcbio-gff 0.6.4** python packages.

In order to install biopython and bcbio-gff please use the following steps:

**Mac OS X / LINUX:**

 Intall pip the python package manager:
 >sudo easy_install pip 

 Install biopython using pip:
 >pip install biopython==1.67

 Install bcbio-gff using pip:
 >pip install bcbio-gff==0.6.4

## Installation

Clone the repository:

 >git clone https://github.com/NBISweden/EMBLmyGFF3.git
 
Move into the folder:
 >cd EMBLmyGFF3/
 
 Executing:
 >./EMBLmyGFF3.py
 
 or
 
 >./EMBLmyGFF3.py -h
 
will display some help.

## USAGE

### FOREWORD

A correct **GFF3 file** and the **genome in FASTA format** that has been used to produce the GFF3 file are the mandatory input files.
Then, in order to get a valid EMBL flat file suitable for submission you have to fill carefully all mandatory metadata.

**/!\ Please be aware that a *project ID* and an *accession number* are mandatory for a submission to [ENA](http://www.ebi.ac.uk/ena). You don't need this information if you don't plan to submit the data (In case you just want an EMBL-like flat file for other purposes). If you don't have yet those information you can add them later by replacing the corresponding fields. Please visit the [EMBL web site](http://www.ebi.ac.uk/ena/support/genome-submission-faq) to learn how to obtain a *project ID* and an *accession number*.**

Test data from the Drosophila melanogaster species are located in the example folder.

### Simple case 

 >./EMBLmyGFF3.py example/maker.gff3 example/maker.fa
 
Will prompt you to fill one by one the mandatory information needed to produce a proper EMBL file.
Once the software has all the information it needs, it will process the input files and will print the result to STDOUT.
 
In order to write the result in the desired file use the **-o** option:
 
 >./EMBLmyGFF3.py example/maker.gff3 example/maker.fa -o result.embl

### Complete case 

Minimum requirement to launch the software and avoid any prompt.

 >./EMBLmyGFF3.py example/maker.gff3 example/maker.fa --data_class STD --topology linear --molecule_type 'genomic DNA' --table 1  --species 'Drosophila melanogaster (fly)' --taxonomy INV --accession ERSXXXXXXX --project_id PRJXXXXXXX --rg MYGROUP -o result.embl

### Advanced case 

When you want add more information than those mandatory: e.g publication.

 >./EMBLmyGFF3.py example/maker.gff3 example/maker.fa --data_class STD --topology linear --molecule_type "genomic DNA" --table 1  --species 'Drosophila melanogaster (fly)' --taxonomy INV --accession ERSXXXXXXX --project_id PRJXXXXXXX --rg MYGROUP --author 'author for the reference' --rt 'reference title' --rl 'Some journal' -o result.embl

### Use through a bash script

You may prefer to launch the software through a bash script especially when you want to fill many information, so we provide examples of such script here **example/maker_example.sh**.

In order to use it move in the example folder, then launch the script:
>./maker_example.sh

## PARAMETER

Some parameters are mandatory and some others are not. Here is a list of all parameters available. 
You can also find a comprehensive help about the different parameters using the software help command.

**positional arguments:**

| Parameter | Description |
| --- | --- |
|  gff_file|              Input gff-file.|
|  fasta|                 Input fasta sequence.|
  
**Mandatory Arguments related to the EMBL format to check carrefully:**

| Parameter | Description |
| --- | --- |
|  -p , --project_id     |Project ID. The defalut value is **UNKNOWN** *(This option is used to set up the PR line.)*|
|  -r , --table          |Translation table. No default value. *(This option is used to set up the translation table qualifier| transl_table of the CDS features.)* Please visit this [NCBI genetic code](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi) page for more information.|
|  -s , --species        |Sample Species, formatted as 'Genus species (english name)'. No default value. This option is used to set up the OS line.|
|  -t , --topology       |Sequence topology. No default value. *(This option is used to set up the Topology that is the 3th token of the ID line.)*|
|  -d , --data_class     |Data class of the sample. No default value. *(This option is used to set up the 5th token of the ID line.)*
|  -m , --molecule_type  |Molecule type of the sample. No default value.|
|  -a , --accession      |Accession number(s) for the entry. Default value: **UNKNOWN** . This option is used to set up the accession number of the AC line and the first token of the ID line as well as the prefix of the locus_tag qualifier.|    
|  -x , --taxonomy       |Source taxonomy. No default value. This option is used to set the taxonomic division within ID line (6th token).|
|  --rg                  |Reference Group, the working groups/consortia that produced the record. No default value.|
  
**Optional arguments related to the software:**

| Parameter | Description |
| --- | --- |
|  -h, --help            |Show this help message and exit.|
|  -v, --verbose         |Increase verbosity.|
|  -q, --quiet           |Decrease verbosity.|
|  --shame               |Suppress the shameless plug.|
|  -z, --gzip            |Gzip output file.|
|  -o , --output         |Output filename.|

**Optional arguments related to the EMBL format:**

| Parameter | Description |
| --- | --- |
|  -c , --created|        Creation time of the original entry. The default value is the **date of the day**.|
|  -g , --organelle|      Sample organelle. No default value.|
|  -k , --keyword|        Keywords for the entry. No default value.|
|  -l , --classification| Organism classification. The default value is **Life**.|
|  --rc|                  Reference Comment. No default value.|
|  --rx|                  Reference cross-reference. No default value.|
|  --ra , --author|       Author for the reference. No default value.|
|  --rt|                  Reference Title. No default value.|
|  --rl|                  Reference publishing location. No default value.|
|  --translate|           Include translation in CDS features. Not activated by default.|
|  --version|             Sequence version number. The default value is **1**.|
|  --keep_duplicates|       Do not remove duplicate features during the process.|
|  --interleave_genes|    Print gene features with interleaved mRNA and CDS features.|

## MAPPING

The challenge for a correct conversion is the correct mapping between the feature types described in the 3th column as well as the different attribute’s tags of the 9th column of the GFF3 file and the corresponding EMBL features and qualifiers.</br>
If you figure out that a feature type or an attribute's tag is not mapped to the corresponding EMBL features or qualifiers you would like, you will have to modify the corresponding information in the mapping files.</br>
The software will skip the unknown feature types (Non EMBL feature types that are not mapped against an EMBL feature type) and the unknown qualifiers (Non EMBL qualifiers that are not mapped against an EMBL qualifier) and will inform you during the conversion process. If you want to include them within the output, you can add the information needed in the corresponding mapping file.


### Feature type

The EMBL format accepts 52 different feature types whereas the GFF3 is constrained to be a Sequence Ontology term or accession number (3th column of the GFF3), but nevertheless this constitutes 2278 terms in version 2.5.3 of the Sequence Ontology.

The file handling the proper mapping is ***translation_gff_feature_to_embl_feature.json***

**example:**

  >"three_prime_UTR": {</br>
  &nbsp;&nbsp;"target": "3'UTR"</br>
  }
 
This will map the **three_prime_UTR** feature type from the 3th column of the GFF3 file to the **3'UTR** EMBL feature type.
**When the feature type from the GFF3 is identical to the EMBL feature no need to specify any target.** If a target is needed and you didn't specified it, the tool will throw a warning message during the process.
 
You can decide which features will be printed in the output using the **remove** parameter:

  >"exon": {</br>
  &nbsp;&nbsp;"remove": true</br>
  }

Like that no exon feature will be display in the output.
 
### GFF3 Attribute to EMBL qualifier

The embl format accepts 98 different qualifiers where the corresponding attribute tag types in the 9th column of the GFF3 are unlimited.
The file handling the proper mapping is ***translation_gff_attribute_to_embl_qualifier.json***

**example:**

  >"Dbxref": {</br>
  &nbsp;&nbsp;"source description": "A database cross reference.",</br>
  &nbsp;&nbsp;"target": "db_xref",</br>
  &nbsp;&nbsp;"dev comment": ""</br>
  },
 
This will map the **Dbxref** attribute's tag from the 9th columm of the GFF3 file to the **db_xref** embl qualifier.

### Other

The **source** (2nd column) as well as the **score** (6th column) from the GFF3 file can also be handled through the ***translation_gff_other_to_embl_qualifier.json*** mapping file.

  >"source": {</br>&nbsp;&nbsp;
  "source description": "The source is a free text qualifier intended to describe the algorithm or operating procedure that generated this
                           feature. Typically this is the name of a piece of software, such as Genescan or a database name, such
                           as Genbank. In effect, the source is used to extend the feature ontology by adding a qualifier to the type 
                           creating a new composite type that is a subclass of the type in the type column.", </br>
  &nbsp;&nbsp;"target": "note",</br>
  &nbsp;&nbsp;"prefix": "source:",</br>
  &nbsp;&nbsp;"dev comment": "EMBL qualifiers tend to be more specific than this, so very hard to create a good mapping."</br>
  },
 
This will map the **source** from the 2nd columm of the GFF3 file to the **note** embl qualifier. 

**/!\\** Please notice the *prefix* allows to add information dowstream the source value wihtin the qualifier (Upstream information is also possible using *suffix*).</br>
e.g: The source value is "Prokka":</br> 
Within the embl file, instead to get **note="Prokka"**, here we will get **note="source:Prokka"**

## KNOWN ISSUES

**biopython version**
There's a bug between bcbio-gff 0.6.4 and biopython 1.68 though, so use biopython 1.67.

If you have several version of biopython or bcbio-gff on your computer it is possible that an incompatible version is used by default which will lead to an execution error. To check the real version used during the execution you can use this command:
 
 >python -c "import Bio; from BCBio import GFF; print 'biopython version: '+Bio.\_\_version\_\_; print 'bcbio-gff version: '+GFF.\_\_version\_\_"

**Duplicated Features**</br>
Features that have the same key (feature type) and location as another feature are considered as duplicates and aren't allowed by the EMBL database. So they are remove during the process. If you don't plan to submit the file to ENA and you wish to keep these features, use the *--keep_duplicates* option. 

## CITATION

[![DOI](https://zenodo.org/badge/64933149.svg)](https://zenodo.org/badge/latestdoi/64933149)

## AUTHOR

Martin Norling<sup>1,2</sup>, Niclas Jareborg<sup>1,3</sup>, Jacques Dainat<sup>1,2</sup>

<sup>1</sup>National Bioinformatics Infrastructure Sweden (NBIS), SciLifeLab, Uppsala Biomedicinska Centrum (BMC), Husargatan 3, S-751 23 Uppsala, SWEDEN.</br>
<sup>2</sup>IMBIM - Department of Medical Biochemistry and Microbiology, Box 582, S-751 23 Uppsala, SWEDEN.</br>
<sup>3</sup>Department of Biochemisty and Biophys-ics, Stockholm University / SciLifeLab, Box 1031, S-171 21 Solna, SWEDEN.
