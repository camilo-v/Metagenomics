#!/bin/sh
#BSUB -J JOB_NAME
#BSUB -o /Path/To/Log/Directory/group_COND_POPULATION_GRP_TEMPLATE_log.txt
#BSUB -W 168:00
#BSUB -N
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
#	The purpose of this script is serve as a template for the permutation runs. LSF jobs are dispatched with this
#	script, which in turns calls the 'permutations.sh' script to execute the mutation & alignment commands.
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

# ----------------------------------------------------- Params --------------------------------------------------------

echo ""
echo "[" `date '+%m/%d/%y %H:%M:%S'` "] Starting..."
echo "[" `date '+%m/%d/%y %H:%M:%S'` "]"

# Location of the permutation script
PERMUTATION_APP_PATH='/Path/To/Permutation/scripts'

GROUP_NUMBER='COND_POPULATION_GRP'

# Metagenomic Sample
SAMPLE_BASE='/Path/To/MetaSim/Samples/group_'$GROUP_NUMBER'/sample'

FASTA_READS_1=$SAMPLE_BASE'/group_'$GROUP_NUMBER'_1.fa'
FASTA_READS_2=$SAMPLE_BASE'/group_'$GROUP_NUMBER'_2.fa'


# ------------------------------------------------ Probability Loop ---------------------------------------------------

PROBABILITIES=( 'TEMPLATE' )


for PROB in "${PROBABILITIES[@]}"
{
	echo "[" `date '+%m/%d/%y %H:%M:%S'` "]"
	echo "[" `date '+%m/%d/%y %H:%M:%S'` "]"
	echo "[" `date '+%m/%d/%y %H:%M:%S'` "] Probabiltiy: "$PROB
	echo "[" `date '+%m/%d/%y %H:%M:%S'` "]"

	# Main Output Directory
	OUTPUT_DIR='/Path/To/Output/Directory/metagenomic_simulations/groups/group_'$GROUP_NUMBER'/permutations/results/reads_permuted@'$PROB
	mkdir -p $OUTPUT_DIR

	# ----------------------------------------------- Permutation Loop -----------------------------------------------

	for ((i = 1; i <= 200; i++))
	{
		echo "[" `date '+%m/%d/%y %H:%M:%S'` "]"
		echo "[" `date '+%m/%d/%y %H:%M:%S'` "] Creating Read Permutation: "$i
		echo "[" `date '+%m/%d/%y %H:%M:%S'` "]"

		INSTANCE="$OUTPUT_DIR/perm_$i"

		mkdir -p "$INSTANCE"

		$PERMUTATION_APP_PATH/permutations.sh $INSTANCE $PROB $FASTA_READS_1 $FASTA_READS_2
	}

}


echo "[" `date '+%m/%d/%y %H:%M:%S'` "] "
echo "[" `date '+%m/%d/%y %H:%M:%S'` "] Done."
echo "[" `date '+%m/%d/%y %H:%M:%S'` "] "