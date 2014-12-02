#!/bin/sh
#BSUB -J JOB_NAME-MERGE_COUNTS
#BSUB -o /Path/To/Log/Directory/permutations/mid_vaginal/merger_counts.txt
#BSUB -W 168:00
#BSUB -u "EMAIL_GOES_HERE"
#BSUB -x
#BSUB -n 16
#BSUB -R "span[ptile=16]"
#BSUB -M 33554432

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
#	The purpose of this program is to drive the fractional-counts jobs for the HMP data.
#
#   NOTES:
#   Please see the dependencies section below for the required libraries (if any).
#
#   DEPENDENCIES:
#
#       • None. Standard bash shell script.
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

COUNT_BASE_DIR='/Path/To/permutations/mid_vaginal/analysis/counts'

APP_PATH='/Path/To/Project/software'

AVERAGE_HIT_DIR=$COUNT_BASE_DIR'/average_hits'

mkdir -p $AVERAGE_HIT_DIR


echo ""

for Q in 2 4 6 8 10 12 14 16 18 20 22 24 26 28 30
{
	echo ""
	echo "[" `date '+%m/%d/%y %H:%M:%S'` "] Q"$Q

	Q_BASE_DIR=$COUNT_BASE_DIR'/'$Q

	# ------------------------------------------------------ Counts ------------------------------------------------------

	echo "[" `date '+%m/%d/%y %H:%M:%S'` "] Q"$Q" — Processing Counts..."

	FILE_TO_REMOVE_2=$Q_BASE_DIR'/'$Q'_counts.txt'

	rm -f $FILE_TO_REMOVE_2

	OUTPUT_FILE_WITH_MERGED_COUNTS=$Q_BASE_DIR'/'$Q'-file_list_with_counts.txt'

	rm -f $OUTPUT_FILE_WITH_MERGED_COUNTS

	FILE_PATTERN_TO_LOOK_FOR=$Q_BASE_DIR'/'$Q'_counts-*.txt'

	ls $FILE_PATTERN_TO_LOOK_FOR >> $OUTPUT_FILE_WITH_MERGED_COUNTS

	$APP_PATH/countMerger.pl -i $OUTPUT_FILE_WITH_MERGED_COUNTS -o $Q_BASE_DIR -p $Q --header --type counts

}

echo ""
echo ""
echo "[" `date '+%m/%d/%y %H:%M:%S'` "] Fin."
echo ""