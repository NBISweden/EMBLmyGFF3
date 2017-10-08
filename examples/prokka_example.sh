#!/bin/bash

########################################################
# Script example to simplify the use of many options #
########################################################

#FASTA file used to produce the annotation
GENOME="prokka.fa"

# ANNOATION in gff3 FORMAT
ANNOTATION="prokka.gff3"

#PROJECT name registered on EMBL
PROJECT="17285"

#Sample accession registered on EMBL
ACCESSION="MY_SAMPLE_ACCESSION"

# species name
SPECIES="escherichia coli"

# Taxonomy
TAXONOMY="PRO"

# The Organism Classification (can be retrieved here http://www.ncbi.nlm.nih.gov/Taxonomy/)
CLASSIFICATION=""

#The working groups/consortia that produced the record. No default value
REFERENCE_GROUP="MyGroup"

#Translation table
TABLE="11"

#Data class of the sample.
CLASS="STD"

#Topology.
TOPOLOGY="circular"

#Molecule type of the sample.
MOLECULE="genomic DNA"

# Converter script location
EMBLmyGFF3="../EMBLmyGFF3.py"

echo -e "Running the following command:\n$EMBLmyGFF3 --rg $REFERENCE_GROUP -a $ACCESSION -p $PROJECT -t $TOPOLOGY -l \"$CLASSIFICATION\" -m \"$MOLECULE\" -d $CLASS -r $TABLE -t linear -s \"$SPECIES\" -q -x $TAXONOMY $ANNOTATION $GENOME $@ > prokka.embl"

$EMBLmyGFF3 --rg $REFERENCE_GROUP -a $ACCESSION -p $PROJECT -t $TOPOLOGY -l "$CLASSIFICATION" -m "$MOLECULE" -d $CLASS -r $TABLE -t linear -s "$SPECIES" -q -x $TAXONOMY $ANNOTATION $GENOME $@ > prokka.embl
