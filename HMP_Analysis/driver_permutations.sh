#!/bin/sh
#BSUB -J JOB_NAME
#BSUB -o /Path/To/Log/Directory/permutations/mid_vaginal/2_log.txt
#BSUB -W 168:00
#BSUB -u "EMAIL_GOES_HERE"
#BSUB -x
#BSUB -n 16
#BSUB -R "span[ptile=16]"

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
#	The purpose of this script is to serve as a template for the permutation runs. LSF jobs are dispatched with this
#	script, which in turns calls the 'permutations.sh' script to execute the mutation & alignment commands in each
#	node.
#
#
#   NOTES:
#   Each HMP sample group has its own "driver" which allows it to execute independently from other samples.
#
#
#   DEPENDENCIES:
#
#       • None.
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

PROB=2

# Location of the driving scripts
RUN_DIR='/Path/To/permutations/mid_vaginal/drivers'

# Define the output directory
OUTPUT_DIR='/Path/To/permutations/mid_vaginal/results/reads_permuted@'$PROB

# We'll create our output directory
mkdir -p $OUTPUT_DIR

# .. And here we go!
for ((i = 207; i <= 350; i++))
{
	echo "Creating permutation $i"

	INSTANCE="$OUTPUT_DIR/perm_$i"

	mkdir -p "$INSTANCE"

	$RUN_DIR/permutations.sh $INSTANCE $PROB
}

# ----------------------------------------------------- Summaries -----------------------------------------------------

SUM_FILE="$OUTPUT_DIR/permutation_summary@$PROB.txt"

touch $SUM_FILE

DATE=`date '+%m/%d/%y %H:%M:%S'`

echo $DATE >> $SUM_FILE

FILE_HEADER=`printf "Permutation\tAlignment Counts\tPercentage"`
echo $FILE_HEADER >> $SUM_FILE

for i in `find $OUTPUT_DIR -name "flagStats.txt"`
do

	PERM=`echo $i | cut -d'/' -f2 | cut -d'_' -f2`

	ALN_COUNT=`grep "mapped (" $i | cut -d' ' -f1`
	ALN_PERCT=`grep "mapped (" $i | cut -d' ' -f5 | sed 's/(//' | sed 's/\:nan\%)//'`

	ROW=`printf "$PERM\t$ALN_COUNT\t$ALN_PERCT"`
	echo $ROW >> $SUM_FILE

done
