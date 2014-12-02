#!/bin/sh
#BSUB -J JOB_NAME
#BSUB -o /Path/To/Log/File/For/permutations/mid_vaginal/2_count_rundown_log.txt
#BSUB -W 168:00
#BSUB -u "EMAIL_GOES_HERE"
#BSUB -x
#BSUB -n 16
#BSUB -R "span[ptile=16]"
#BSUB -M 33554432

# ---------------------------------------------------------------------------------------------------------------------
#
#                                   		Center for Computational Science
#												http://www.ccs.miami.edu/
#                             			  			University of Miami
#
#   This software is a "University of Miami Work" under the terms of the United States Copyright Act.
#   Please cite the author(s) in any work or product based on this material.
#
#   OBJECTIVE:
#	The purpose of this program is to obtain the fractional counts for the HMP dataset.  This script calls the
#	"perm_rundown.pl" Perl script that actually obtains the counts.  The 'low' and 'high' parameters in the Perl script
#	demarcate the lower and upper bounds to slice the permutation directories by. Due to the size of the alignment
#	files, processing all the permutations at once is not feassible, so they are split into batches that can then
#	be easily processed.  The boundaries to split by, are identified by 'low' and 'high' — they are interval markers.
#
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

PROB='PROBABILITY_VALUE'

# 	Path to framework (software)
PROJECT_DIR='/Path/To/Project/software'

#	Path to FASTA of all IMG database
REF_SEQS='/Paht/To/annotations/img/4/bacteria/fasta/img_bacteria.fa'

#	Reference Names — from FASTA file ">" header records
REF_NAMES='/Path/To/annotations/img/4/bacteria/metadata/references.txt'

# 	Output Directory for sum counts
OUTPUT_DIR='/Path/to/metagenomics/permutations/mid_vaginal/analysis/counts/'$PROB

mkdir -p $OUTPUT_DIR

# 	Directory for permutations
PERMUTATIONS_BASE='/Path/to/metagenomics/permutations/mid_vaginal/results'

PERMUTATIONS_DIR=$PERMUTATIONS_BASE'/reads_permuted@'$PROB

# Use a MAPQ filter to filter-out low quality reads
MAPQ='10'

#	Build & Analyze
$PROJECT_DIR/perm_rundown.pl -i $REF_NAMES -o $OUTPUT_DIR -p $PERMUTATIONS_DIR -f $REF_SEQS --mapq $MAPQ --prefix $PROB --read_hits --low 1 --high 25
$PROJECT_DIR/perm_rundown.pl -i $REF_NAMES -o $OUTPUT_DIR -p $PERMUTATIONS_DIR -f $REF_SEQS --mapq $MAPQ --prefix $PROB --read_hits --low 26 --high 50
$PROJECT_DIR/perm_rundown.pl -i $REF_NAMES -o $OUTPUT_DIR -p $PERMUTATIONS_DIR -f $REF_SEQS --mapq $MAPQ --prefix $PROB --read_hits --low 51 --high 75
$PROJECT_DIR/perm_rundown.pl -i $REF_NAMES -o $OUTPUT_DIR -p $PERMUTATIONS_DIR -f $REF_SEQS --mapq $MAPQ --prefix $PROB --read_hits --low 76 --high 100
$PROJECT_DIR/perm_rundown.pl -i $REF_NAMES -o $OUTPUT_DIR -p $PERMUTATIONS_DIR -f $REF_SEQS --mapq $MAPQ --prefix $PROB --read_hits --low 101 --high 125
$PROJECT_DIR/perm_rundown.pl -i $REF_NAMES -o $OUTPUT_DIR -p $PERMUTATIONS_DIR -f $REF_SEQS --mapq $MAPQ --prefix $PROB --read_hits --low 126 --high 150
$PROJECT_DIR/perm_rundown.pl -i $REF_NAMES -o $OUTPUT_DIR -p $PERMUTATIONS_DIR -f $REF_SEQS --mapq $MAPQ --prefix $PROB --read_hits --low 151 --high 175
$PROJECT_DIR/perm_rundown.pl -i $REF_NAMES -o $OUTPUT_DIR -p $PERMUTATIONS_DIR -f $REF_SEQS --mapq $MAPQ --prefix $PROB --read_hits --low 176 --high 200
$PROJECT_DIR/perm_rundown.pl -i $REF_NAMES -o $OUTPUT_DIR -p $PERMUTATIONS_DIR -f $REF_SEQS --mapq $MAPQ --prefix $PROB --read_hits --low 201 --high 225
$PROJECT_DIR/perm_rundown.pl -i $REF_NAMES -o $OUTPUT_DIR -p $PERMUTATIONS_DIR -f $REF_SEQS --mapq $MAPQ --prefix $PROB --read_hits --low 226 --high 250
$PROJECT_DIR/perm_rundown.pl -i $REF_NAMES -o $OUTPUT_DIR -p $PERMUTATIONS_DIR -f $REF_SEQS --mapq $MAPQ --prefix $PROB --read_hits --low 251 --high 275
$PROJECT_DIR/perm_rundown.pl -i $REF_NAMES -o $OUTPUT_DIR -p $PERMUTATIONS_DIR -f $REF_SEQS --mapq $MAPQ --prefix $PROB --read_hits --low 276 --high 300
$PROJECT_DIR/perm_rundown.pl -i $REF_NAMES -o $OUTPUT_DIR -p $PERMUTATIONS_DIR -f $REF_SEQS --mapq $MAPQ --prefix $PROB --read_hits --low 301 --high 325
$PROJECT_DIR/perm_rundown.pl -i $REF_NAMES -o $OUTPUT_DIR -p $PERMUTATIONS_DIR -f $REF_SEQS --mapq $MAPQ --prefix $PROB --read_hits --low 326 --high 350