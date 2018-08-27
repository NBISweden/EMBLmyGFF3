#!/bin/bash

########################################################
# Script example to simplify the use of many options #
########################################################

#PATH to the FASTA file used to produce the annotation
GENOME=`dirname "$0"`"/augustus.fa"

#PATH to the ANNOTATION in gff3 FORMAT
ANNOTATION=`dirname "$0"`"/augustus.gff3"

#PROJECT name registered on EMBL
PROJECT="17285"

#Locus tag registered on EMBL
LOCUS_TAG="MYLOCUSTAG"

# species name
SPECIES="Drosophila melanogaster"

# Taxonomy
TAXONOMY="INV"

#The working groups/consortia that produced the record. No default value
REFERENCE_GROUP="MyGroup"

#Translation table
TABLE="1"

#Molecule type of the sample.
MOLECULE="genomic DNA"

myCommand="EMBLmyGFF3 --rg $REFERENCE_GROUP -i $LOCUS_TAG -p $PROJECT -m \"$MOLECULE\" -r $TABLE -t linear -s \"$SPECIES\" -x $TAXONOMY -o EMBLmyGFF3-augustus-example.embl $ANNOTATION $GENOME $@"
echo -e "Running the following command:\n$myCommand"

#execute the command
eval $myCommand
