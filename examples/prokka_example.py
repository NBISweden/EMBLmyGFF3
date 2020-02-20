#!/usr/bin/env python3

#############################
# Test case on prokka gff3  #
#############################

import subprocess
import os
import sys

def fill_path(file):
	path =  os.path.realpath(__file__)
	tail = path.rsplit('/',1)
	return tail[0]+"/"+file

def main():
	#PATH to the FASTA file used to produce the annotation
	GENOME="prokka.fa"

	#PATH to the ANNOTATION in gff3 FORMAT
	ANNOTATION="prokka.gff3"

	#PROJECT name registered on EMBL
	PROJECT="17285"

	#Locus tag registered on EMBL
	LOCUS_TAG="MYLOCUSTAG"

	# species name
	SPECIES="escherichia coli"

	# Taxonomy
	TAXONOMY="PRO"

	#The working groups/consortia that produced the record. No default value
	REFERENCE_GROUP="MyGroup"

	#Translation table
	TABLE="11"

	#Topology.
	TOPOLOGY="circular"

	#Molecule type of the sample.
	MOLECULE="genomic DNA"

	#Strain
	STRAIN="K-12"

	#Create the command
	command = "EMBLmyGFF3 --rg REFERENCE_GROUP -i "+LOCUS_TAG+" -p "+PROJECT+" -m \""+MOLECULE+"\" -r "+TABLE+" -t "+TOPOLOGY+" --strain \""+STRAIN+"\" -s \""+SPECIES+"\" -x "+TAXONOMY+" -o EMBLmyGFF3-prokka-example.embl "+fill_path(ANNOTATION)+" "+fill_path(GENOME)
	print("Running the following command: "+command)

	#Execute the command
	process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
	stdout, stderr = process.communicate()
	print(stdout.decode('utf-8'))
	sys.exit(process.returncode)

if __name__ == '__main__':
	main()
