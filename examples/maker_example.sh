#!/bin/bash

########################################################
# Script example to simplify the use of many options #
########################################################

#FASTA file used to produce the annotation
GENOME="maker.fa"

# ANNOATION in gff3 FORMAT
ANNOTATION="maker.gff3"

#PROJECT name registered on EMBL
PROJECT="17285"

#Locus tag registered on EMBL
LOCUS_TAG="MY_LOCUS_TAG"

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

myCommand="EMBLmyGFF3 --rg $REFERENCE_GROUP -i $LOCUS_TAG -p $PROJECT -m \"$MOLECULE\" -r $TABLE -t linear -s \"$SPECIES\" -x $TAXONOMY -o EMBLmyGFF3-maker-test.embl $ANNOTATION $GENOME"
echo -e "Running the following command:\n$myCommand"

#execute the command
eval $myCommand
