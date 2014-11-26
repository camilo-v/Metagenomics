#!/bin/sh
#BSUB -J JOB_NAME-MERGE_READ_HITS
#BSUB -o /Path/To/Log/Directory/permutations/mid_vaginal/merger_read_hit_distribution.txt
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
#	The purpose of this program is to drive the collection scripts for the fractional counts and obtain the average number
#	of genomic hits (to the reference database) for each NGS read.
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

COUNT_BASE_DIR='/Path/To/metagenomics/permutations/mid_vaginal/analysis/counts'

APP_PATH='/Path/To/Project/software'

AVERAGE_HIT_DIR=$COUNT_BASE_DIR'/average_hits'

mkdir -p $AVERAGE_HIT_DIR


echo ""

for Q in 2 4 6 8 10 12 14 16 18 20 22 24 26 28 30
{
	echo ""
	echo "[" `date '+%m/%d/%y %H:%M:%S'` "] Q"$Q

	Q_BASE_DIR=$COUNT_BASE_DIR'/'$Q

	#------------------------------------------------------------------------------------------------
	#
	#	Read-Hit Distribution
	#

	echo "[" `date '+%m/%d/%y %H:%M:%S'` "] Q"$Q" — Processing Hit Distributions..."

	FILE_TO_REMOVE_DISTRIBS=$Q_BASE_DIR'/'$Q'_read_hit_distribution.txt'
	rm -f $FILE_TO_REMOVE_DISTRIBS

	mv $Q_BASE_DIR'/'$Q'_read_hit_distribution-1-25.txt' $Q_BASE_DIR'/'$Q'_read_hit_distribution-001-25.txt'
	mv $Q_BASE_DIR'/'$Q'_read_hit_distribution-26-50.txt' $Q_BASE_DIR'/'$Q'_read_hit_distribution-026-50.txt'
	mv $Q_BASE_DIR'/'$Q'_read_hit_distribution-51-75.txt' $Q_BASE_DIR'/'$Q'_read_hit_distribution-051-75.txt'
	mv $Q_BASE_DIR'/'$Q'_read_hit_distribution-76-100.txt' $Q_BASE_DIR'/'$Q'_read_hit_distribution-076-100.txt'

	# First batch of merges -- alignments are too big
	for BATCH in 0 1 2 3
	{
		echo "[" `date '+%m/%d/%y %H:%M:%S'` "] Q"$Q" — Batch "$BATCH

		BATCH_OUTPUT_FILE_WITH_MERGED_DISTS=$Q_BASE_DIR'/'$Q'-batch_'$BATCH'-file_list_with_distributions.txt'

		rm -f $BATCH_OUTPUT_FILE_WITH_MERGED_DISTS

		BATCH_FILE_PATTERN=$Q_BASE_DIR'/'$Q'_read_hit_distribution-'$BATCH'*.txt'

		# Aggregate the files we need into a list file
		ls $BATCH_FILE_PATTERN >> $BATCH_OUTPUT_FILE_WITH_MERGED_DISTS

		TYPE='read_hit_distribution-batch_'$BATCH

		$APP_PATH/countMerger.pl -i $BATCH_OUTPUT_FILE_WITH_MERGED_DISTS -o $Q_BASE_DIR -p $Q --omit --type $TYPE
	}

	# A-B Merge

	echo "[" `date '+%m/%d/%y %H:%M:%S'` "] A-B Merge..."

	#-----------------------------------------------------------------
	#
	#	A
	#
	OUTPUT_FILE_FOR_LAST_MERGE_A=$Q_BASE_DIR'/'$Q'-final_merge_file_list_for_hit_distribution_A.txt'

	rm -f $OUTPUT_FILE_FOR_LAST_MERGE_A

	PATERN_FOR_LAST_MERGE_A=$Q_BASE_DIR'/'$Q'_read_hit_distribution-batch_[0-1].txt'

	ls $PATERN_FOR_LAST_MERGE_A >> $OUTPUT_FILE_FOR_LAST_MERGE_A

	TYPE_FOR_LAST_MERGE_A='read_hit_distribution_A'

	$APP_PATH/countMerger.pl -i $OUTPUT_FILE_FOR_LAST_MERGE_A -o $Q_BASE_DIR -p $Q --omit --type $TYPE_FOR_LAST_MERGE_A --header

	#-----------------------------------------------------------------
	#
	#	B
	#
	OUTPUT_FILE_FOR_LAST_MERGE_B=$Q_BASE_DIR'/'$Q'-final_merge_file_list_for_hit_distribution_B.txt'

	rm -f $OUTPUT_FILE_FOR_LAST_MERGE_B

	PATERN_FOR_LAST_MERGE_B=$Q_BASE_DIR'/'$Q'_read_hit_distribution-batch_[2-3].txt'

	ls $PATERN_FOR_LAST_MERGE_B >> $OUTPUT_FILE_FOR_LAST_MERGE_B

	TYPE_FOR_LAST_MERGE_B='read_hit_distribution_B'

	$APP_PATH/countMerger.pl -i $OUTPUT_FILE_FOR_LAST_MERGE_B -o $Q_BASE_DIR -p $Q --omit --type $TYPE_FOR_LAST_MERGE_B --header

	#-----------------------------------------------------------------
	#
	#	FINAL
	#

	echo "[" `date '+%m/%d/%y %H:%M:%S'` "] Final Merge..."

	OUTPUT_FILE_FOR_LAST_MERGE_FINAL=$Q_BASE_DIR'/'$Q'-final_merge_file_list_for_hit_distribution_FINAL.txt'

	rm -f $OUTPUT_FILE_FOR_LAST_MERGE_FINAL

	PATERN_FOR_LAST_MERGE_FINAL=$Q_BASE_DIR'/'$Q'_read_hit_distribution_[AB].txt'

	ls $PATERN_FOR_LAST_MERGE_FINAL >> $OUTPUT_FILE_FOR_LAST_MERGE_FINAL

	TYPE_FOR_LAST_MERGE_FINAL='read_hit_distribution_FINAL'

	$APP_PATH/countMerger.pl -i $OUTPUT_FILE_FOR_LAST_MERGE_FINAL -o $Q_BASE_DIR -p $Q --omit --type $TYPE_FOR_LAST_MERGE_FINAL --header


	#-----------------------------------------------------------------
	#
	#	Averages
	#	We are only interested in the mean number of hits for each read, which are in the last column
	#

	echo "[" `date '+%m/%d/%y %H:%M:%S'` "] Extracting Averages..."

	FINAL_FILE_WITH_HITS=$Q_BASE_DIR'/'$Q'_read_hit_distribution_FINAL.txt'

	OUTPUT_FILE_WITH_AVERAGE_HITS=$Q_BASE_DIR'/'$Q'_average_hits.txt'

	COLUMN_TO_CUT='352'

	cut -f $COLUMN_TO_CUT $FINAL_FILE_WITH_HITS > $OUTPUT_FILE_WITH_AVERAGE_HITS

	cp $OUTPUT_FILE_WITH_AVERAGE_HITS $AVERAGE_HIT_DIR

}

echo ""
echo ""
echo "[" `date '+%m/%d/%y %H:%M:%S'` "] Fin."
echo ""