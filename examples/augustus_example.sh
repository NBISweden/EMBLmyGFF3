#!/bin/bash

########################################################
# Script example to simplify the use of many options #
########################################################

#FASTA file used to produce the annotation
GENOME="augustus.fa"

# ANNOATION in gff3 FORMAT
ANNOTATION="augustus.gff3"

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

# Converter script location
EMBLmyGFF3="../EMBLmyGFF3.py"

myCommand="$EMBLmyGFF3 --rg $REFERENCE_GROUP -i $LOCUS_TAG -p $PROJECT -m \"$MOLECULE\" -r $TABLE -t linear -s \"$SPECIES\" -x $TAXONOMY $ANNOTATION $GENOME $@ > augustus.embl"
echo -e "Running the following command:\n$myCommand"

#execute the command
eval $myCommand
