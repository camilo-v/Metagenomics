#!/bin/sh
#BSUB -J JOB_NAME
#BSUB -o /Path/To/LSF/Log/Directory/simulations/group_COND_POPULATION_GRP_TEMPLATE_count_rundown_log.txt
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
#	The purpose of this script is to serve as a template for the fractional-counts factory script.  The factory script
#	will create a stand alone scripts that will execute the counting script, "perm_rundown_simulations.pl" for each of
#	the subgroup jobs.  Each of these scripts is a submission job for LSF in the Pegasus2 cluster.
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


PROB='TEMPLATE'

GROUP_NUMBER='COND_POPULATION_GRP'

# 	Path to framework (software)
PROJECT_DIR='/Path/To/metagenomic_simulations/scripts'

#	Path to FASTA of all IMG database
REF_SEQS='/Path/To/references/simulations/Conditional_Population.fa'

#	Reference Names — from FASTA file ">" header records
REF_NAMES='/Path/To/references/simulations/references.txt'

# 	Output Directory for sum counts
OUTPUT_DIR='/Path/To/metagenomic_simulations/groups/group_'$GROUP_NUMBER'/analysis/counts/'$PROB

mkdir -p $OUTPUT_DIR

# 	Directory for permutations
PERMUTATIONS_BASE='/Path/To/metagenomic_simulations/groups/group_'$GROUP_NUMBER'/permutations/results'

PERMUTATIONS_DIR=$PERMUTATIONS_BASE'/reads_permuted@'$PROB

# Use a MAPQ filter to filter-out low quality reads
MAPQ='10'

#	Build & Analyze
$PROJECT_DIR/perm_rundown_simulations.pl -i $REF_NAMES -o $OUTPUT_DIR -p $PERMUTATIONS_DIR -f $REF_SEQS --mapq $MAPQ --prefix $PROB --read_hits --low 1 --high 25
$PROJECT_DIR/perm_rundown_simulations.pl -i $REF_NAMES -o $OUTPUT_DIR -p $PERMUTATIONS_DIR -f $REF_SEQS --mapq $MAPQ --prefix $PROB --read_hits --low 26 --high 50
$PROJECT_DIR/perm_rundown_simulations.pl -i $REF_NAMES -o $OUTPUT_DIR -p $PERMUTATIONS_DIR -f $REF_SEQS --mapq $MAPQ --prefix $PROB --read_hits --low 51 --high 75
$PROJECT_DIR/perm_rundown_simulations.pl -i $REF_NAMES -o $OUTPUT_DIR -p $PERMUTATIONS_DIR -f $REF_SEQS --mapq $MAPQ --prefix $PROB --read_hits --low 76 --high 100
$PROJECT_DIR/perm_rundown_simulations.pl -i $REF_NAMES -o $OUTPUT_DIR -p $PERMUTATIONS_DIR -f $REF_SEQS --mapq $MAPQ --prefix $PROB --read_hits --low 101 --high 125
$PROJECT_DIR/perm_rundown_simulations.pl -i $REF_NAMES -o $OUTPUT_DIR -p $PERMUTATIONS_DIR -f $REF_SEQS --mapq $MAPQ --prefix $PROB --read_hits --low 126 --high 150
$PROJECT_DIR/perm_rundown_simulations.pl -i $REF_NAMES -o $OUTPUT_DIR -p $PERMUTATIONS_DIR -f $REF_SEQS --mapq $MAPQ --prefix $PROB --read_hits --low 151 --high 175
$PROJECT_DIR/perm_rundown_simulations.pl -i $REF_NAMES -o $OUTPUT_DIR -p $PERMUTATIONS_DIR -f $REF_SEQS --mapq $MAPQ --prefix $PROB --read_hits --low 176 --high 200