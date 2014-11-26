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
#	The purpose of this script is to split a conditional population file into sub-files that contains a
#	maximum of 50 entries each.  This constraint is due to MetaSim - it cannot load a large database all
#	at once.  A starting conditional population of 600 generates a file of about 2GB, and MetaSim is just
#	incapable of loading it.  So we have to split the large file into files that it can handle -- 50 has
#	been the magic number so far.
#
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

BASE_DIR='/Path/To/Simulations/Directory/Random_Picks'

CONDITIONAL_POPULATION_TO_SPLIT=$BASE_DIR'/Conditional_Population/conditional_population.txt'

OUTPUT_DIR=$BASE_DIR'/Conditional_Population'

# Number of lines per output file --> MetaSim cannot handle more than this many sequences at a time
NUMBER_OF_LINES='50'

PREFIX='conditional_population_'
SUFFIX_LENGTH='1'

# --------------------------------------------- Script ---------------------------------------------


cd $BASE_DIR'/Conditional_Population'

split -a $SUFFIX_LENGTH -l $NUMBER_OF_LINES $CONDITIONAL_POPULATION_TO_SPLIT $PREFIX



echo
echo "[" `date '+%m/%d/%y %H:%M:%S'` "] Done."
echo