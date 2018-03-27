#!/bin/bash

########################################################
# Script example to simplify the use of many options #
########################################################

#PATH to the FASTA file used to produce the annotation
GENOME=`dirname "$0"`"/prokka.fa"

#PATH to the ANNOTATION in gff3 FORMAT
ANNOTATION=`dirname "$0"`"/prokka.gff3"

#PROJECT name registered on EMBL
PROJECT="17285"

#Locus tag registered on EMBL
LOCUS_TAG="MY_LOCUS_TAG"

# species name
SPECIES="escherichia coli"

# Taxonomy
TAXONOMY="PRO"

#The working groups/consortia that produced the record. No default value
REFERENCE_GROUP="MyGroup"

#Translation table
TABLE="11"

#Strain
STRAIN="K-12"


#Topology.
TOPOLOGY="circular"

#Molecule type of the sample.
MOLECULE="genomic DNA"

myCommand="EMBLmyGFF3 --rg $REFERENCE_GROUP -i $LOCUS_TAG -p $PROJECT -t $TOPOLOGY -m \"$MOLECULE\" -r $TABLE -t linear --strain $STRAIN -s \"$SPECIES\" -q -x $TAXONOMY -o EMBLmyGFF3-prokka-example.embl $ANNOTATION $GENOME $@"
echo -e "Running the following command:\n$myCommand"

#execute the command
eval $myCommand
