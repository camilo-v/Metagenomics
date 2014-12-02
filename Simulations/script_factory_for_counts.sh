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
#	The purpose of this script is to create the counting scripts from that will get submitted to the Pegasus2 cluster
#	(via LSF).  Each script is created from the "count_rundown_TEMPLATE.sh" blueprint script.  The scheme here is to copy
#	the template script 'N' number of times and replace markers in the template with group specific variables and
#	parameters — with 'N' being the number of of subgroups from the conditional population plus the number of mutation
#	rate probabilities we wish to test.
#
#
#   NOTES:
#   Please see the dependencies section below for the required libraries (if any).
#
#   DEPENDENCIES:
#
#       • None. Uses standard modules.
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

OUTPUT_PATH='/Path/to/Simulations/Scripts/counts/counting_scripts'

TEMPLATE_SCRIPT_PATH='/Paht/to/Simulations/Scripts/counts/count_rundown_TEMPLATE.sh'

GROUPS_TO_PROCESS=( '0' '1' '2' '3' '4' '5' '6' '7' '8' '9' '10' )

PROBABILITIES=( '0' '2' '4' '6' '8' '10' '12' '14' '16' '18' '20' '22' '24' '26' '28' '30' )

for GROUP in "${GROUPS_TO_PROCESS[@]}"
{
	echo "[" `date '+%m/%d/%y %H:%M:%S'` "] Creating Scripts for Group: "$GROUP

	for PROB in "${PROBABILITIES[@]}"
	{
		echo "[" `date '+%m/%d/%y %H:%M:%S'` "] Creating "$PROB" Script..."

		OUTPUT_DIR=$OUTPUT_PATH'/group_'$GROUP

		mkdir -p $OUTPUT_DIR

		SCRIPT_NAME=$OUTPUT_DIR'/driver_count_rundown_group_'$GROUP'_prob_'$PROB'.sh'

		# Remove the old copy of the script
		rm -f $SCRIPT_NAME

		# Create the new copy of the script
		touch $SCRIPT_NAME

		JOB_NAME='group_'$GROUP'_'$PROB'_count'

		NEW_SCRIPT=`sed s/TEMPLATE/$PROB/g $TEMPLATE_SCRIPT_PATH | sed s/COND_POPULATION_GRP/$GROUP/g | sed s/JOB_NAME/$JOB_NAME/g`

		echo "$NEW_SCRIPT" >> $SCRIPT_NAME

		chmod ug+rwx $SCRIPT_NAME

	}
}

echo ""
echo "[" `date '+%m/%d/%y %H:%M:%S'` "] Done."
echo "[" `date '+%m/%d/%y %H:%M:%S'` "]"
echo ""