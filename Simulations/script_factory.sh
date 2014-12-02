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
#	The purpose of this script is create the necessary worker scripts to run the permutaion jobs in a distributed
#	cluster.  The factory script will take the template script and replace the adequate parameters so that each group
#	job is submitted to its own computational node, and executes independetly from the other group jobs.
#
#   NOTES:
#   Please see the dependencies section below for the required.
#
#   DEPENDENCIES:
#
#       • None.
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
echo "[" `date '+%m/%d/%y %H:%M:%S'` "]"

OUTPUT_PATH='/Path/To/Output/Directory/Simulations/Scripts/permutations/scratch_runs/drivers'

TEMPLATE_SCRIPT_PATH='/Path/To/Template/Script/Simulations/Scripts/permutations/scratch_runs/drivers/driver_permutations_group_TEMPLATE.sh'

GROUPS_TO_PROCESS=( '0' '1' '2' '3' '4' '5' '6' '7' '8' '9' )

PROBABILITIES=( '0' '2' '4' '6' '8' '10' '12' '14' '16' '18' '20' '22' '24' '26' '28' '30' )

for GROUP in "${GROUPS_TO_PROCESS[@]}"
{
	echo "[" `date '+%m/%d/%y %H:%M:%S'` "] Creating Scripts for Group: "$GROUP

	for PROB in "${PROBABILITIES[@]}"
	{
		echo "[" `date '+%m/%d/%y %H:%M:%S'` "] Creating "$PROB" Script..."

		OUTPUT_DIR=$OUTPUT_PATH'/group_'$GROUP

		mkdir -p $OUTPUT_DIR

		SCRIPT_NAME=$OUTPUT_DIR'/driver_permutations_group_'$GROUP'_'$PROB'.sh'

		# Remove the old copy of the script
		rm -f $SCRIPT_NAME

		# Create the new copy of the script
		touch $SCRIPT_NAME

		JOB_NAME='group_'$GROUP'_'$PROB

		NEW_SCRIPT=`sed s/TEMPLATE/$PROB/g $TEMPLATE_SCRIPT_PATH | sed s/COND_POPULATION_GRP/$GROUP/g | sed s/JOB_NAME/$JOB_NAME/g`

		echo "$NEW_SCRIPT" >> $SCRIPT_NAME

		chmod ug+rwx $SCRIPT_NAME
	}

}


echo ""
echo "[" `date '+%m/%d/%y %H:%M:%S'` "]"
echo "[" `date '+%m/%d/%y %H:%M:%S'` "] Done."
echo "[" `date '+%m/%d/%y %H:%M:%S'` "]"
echo ""