#!/bin/sh

# ---------------------------------------------------------------------------------------------------------------------
#
#                                   Center for Computational Science
#										http://www.ccs.miami.edu/
#                             			  University of Miami
#
#   This software is a "University of Miami Work" under the terms of the United States Copyright Act.
#   Please cite the author(s) in any work or product based on this material.
#
#   OBJECTIVE:
#	The purpose of this script is to collect the list of bacterial strain names to use in the selection of a
#	Conditional Population.
#
#   NOTES:
#   Please see the dependencies section below for the required libraries (if any).
#
#   DEPENDENCIES:
#
#       • None. Uses standard python foundation classes.
#
#   The above libraries & modules are required. You can check the modules currently installed in your
#   system by running: python -c "help('modules')"
#
#   USAGE:
#   Run the program with the "--help" flag to see usage instructions.
#
#
#
#	AUTHOR:	Camilo Valdes (cvaldes3@med.miami.edu)
#			Computational Biology and Bioinformatics Group
#			Center for Computational Science (CCS), University of Miami
#
#
#   COLLABORATORS:
#
#           Jennifer Clarke (jclarke3@unl.edu)
#           Director, Computational Sciences Initiative
#           Associate Professor Department of Food Science and Technology & Department of Statistics
#           University of Nebraska-Lincoln
#
#           Bertrand Clarke (bclarke3@unl.edu)
#           Professor and Department Head, Department of Statistics
#           University of Nebraska-Lincoln
#
#
# ---------------------------------------------------------------------------------------------------------------------


echo ""
echo "[" `date '+%m/%d/%y %H:%M:%S'` "] Starting..."

BACTERIAL_GENOMES_DIR='/Path/To/Bacteria/Genome/Directory'

OUTPUT_FILE='/Path/To/Output/File/Directory/list_bacterial_strain_sequence_names.txt'

rm -fR $OUTPUT_FILE

for GENOME_DIR in `ls $BACTERIAL_GENOMES_DIR`
{
	echo "[" `date '+%m/%d/%y %H:%M:%S'` "] Processing "$GENOME_DIR

	FULL_PATH_GENOME_DIR=$BACTERIAL_GENOMES_DIR'/'$GENOME_DIR

	cd $FULL_PATH_GENOME_DIR

	cat *.fna | grep ">" | sed 's/| /|/' >> $OUTPUT_FILE

}


echo
echo "[" `date '+%m/%d/%y %H:%M:%S'` "] Done."
echo