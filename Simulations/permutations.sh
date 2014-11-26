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
#	The purpose of this script is to call the mutation script and to drive the alignment of the metagenomic sample to
#	the reference database.
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

PROJECT_DIR='/Path/To/Project/software'

OUTPUT_DIR_READS=$1

PROB=$2

FASTA_READS_1=$3
FASTA_READS_2=$4


PROC_THREADS='32'

# ------------------------------------------------ Read Permutation ---------------------------------------------------
#
# Permute the Reads (paired-ends) IFF its not the baseline (0) run
#

if [ $PROB != 0 ]; then
	$PROJECT_DIR/sequencePermuter.pl -p $PROB -t $PROC_THREADS -i $FASTA_READS_1 -o $OUTPUT_DIR_READS
	mv $OUTPUT_DIR_READS/permutedSequences.fasta $OUTPUT_DIR_READS/permutedSequences-1.fasta

	$PROJECT_DIR/sequencePermuter.pl -p $PROB -t $PROC_THREADS -i $FASTA_READS_2 -o $OUTPUT_DIR_READS
	mv $OUTPUT_DIR_READS/permutedSequences.fasta $OUTPUT_DIR_READS/permutedSequences-2.fasta
fi


# ---------------------------------------------------- Alignment ------------------------------------------------------

#	SAM Header Location
# HEADER='/projects/bioinf/ref/annotations/img/4/bacteria/metadata/header.sam'

# 	Reads
READS_1="reads1File"
READS_2="reads2File"

# 	If this is the baseline run, we'll just use the pristine reads
if [ $PROB == 0 ]; then
	READS_1=$FASTA_READS_1
	READS_2=$FASTA_READS_2
fi

#	If this is not the baseline run, we'll use the permuted reads
if [ $PROB != 0 ]; then
	READS_1=$OUTPUT_DIR_READS/permutedSequences-1.fasta
	READS_2=$OUTPUT_DIR_READS/permutedSequences-2.fasta
fi

# 	Output Directory for alignment run

OUTDIR=$1/bowtie
mkdir -p $OUTDIR

# 	Bowtie Index
BOWTIE_INDEX='/Path/To/Bowtie2/Index/Files/Conditional_Population'


# ----------------------------------------------------- Bowtie --------------------------------------------------------

# Location of the Bowtie app executable
BOWTIE2_APP_PATH='/Path/To/Bowtie2/Aligner/Installation/bowtie/v2/2.1.0'

# Filename for the SAM/BAM files
FILENAME="alignments"

# Local Alignemnt and Gap Alignment Parameters for bowtie
SEED_LENGTH=20
SEED_INTERVAL="S,1,0.50"

ALIGNMENT_SUMMARY_FILE=$OUTDIR'/alignment_summary.txt'

$BOWTIE2_APP_PATH/bowtie2 --local -D 20 -R 3 -N 0 -L $SEED_LENGTH -i $SEED_INTERVAL -a --time -f -x $BOWTIE_INDEX -1 $READS_1 -2 $READS_2 -S $OUTDIR/$FILENAME.sam 2> $ALIGNMENT_SUMMARY_FILE


# ---------------------------------------------------- Samtools ------------------------------------------------------

SAMTOOLS_APP_PATH='/Path/To/SamTools/Installation/samtools/0.1.19'

# Convert the Bowtie SAM output to standard BAM
/nethome/bioinfo/tools/samtools/0.1.19/samtools view -b -S -o $OUTDIR/$FILENAME.bam $OUTDIR/$FILENAME.sam

FLAGSTATS_FILE=$OUTDIR'/flagStats.txt'
REF_COUNTS_FILE=$OUTDIR'/reference_Counts.txt'

# Sort the alignments and get basic alignment statistics
$SAMTOOLS_APP_PATH/samtools sort $OUTDIR/$FILENAME.bam $OUTDIR/$FILENAME.sorted
$SAMTOOLS_APP_PATH/samtools index $OUTDIR/$FILENAME.sorted.bam
$SAMTOOLS_APP_PATH/samtools flagstat $OUTDIR/$FILENAME.sorted.bam > $FLAGSTATS_FILE
$SAMTOOLS_APP_PATH/samtools idxstats $OUTDIR/$FILENAME.sorted.bam > $REF_COUNTS_FILE


# ----------------------------------------------------- Cleanup -------------------------------------------------------

rm -fR $OUTDIR/$FILENAME.sam
rm -fR $OUTDIR/$FILENAME.bam

# --------------------------------------------------- End of Line -----------------------------------------------------
#
# Have a nice day!
# :)
#
