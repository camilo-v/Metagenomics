#!/bin/sh


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
#	The purpose of this script is to create a single permutation of the mutated NGS reads (the sample).  This script is
#	called N-number of times by "driver_permutations_group_TEMPLATE.sh", where "N" is the number of permutations desired.
#	This script drives the bulk of the work: it aligns the NGS reads at different mutation rates: from '0' (no mutations)
#	all the way to '.30' in intervals of '0.02'.
#	The main task of this script is to mutate the reads at a given mutation rate ($PROB) by calling the perl script
#	"sequencePermuter.pl"; it then aligns the output (the newly mutated NGS sample) to the reference conditional
#	population database of strains using bowtie2.
#	Samtools is used in the end to sort and index the alignments, as well as to obtain basic stats ('flagstat').
#
#
#   NOTES:
#   The script uses the bowtie2 DNA aligner for the alignment step.  A bowtie2 index must have been previously created
#	with the necessary conditional population stratins. Also, see the dependencies section below for the required
#	binaries.
#	Note that at the time this scrip was developed, Bowtie2 did not support large indices (it does so now).  So the app
#	that created the index, "bowtie2-build", was limitted to creating indices from sequence databases smaller than
#	4 GB, a 32-bit limitation.  This is the main reason why this script calls bowtie multiple times: so that it can
#	align the reads to mutliple slices of the full index. The SAM/BAM results are then merged using samtools.
#	Please see the paper supplementals for a detailed explanation.
#
#
#   DEPENDENCIES:
#
#       • sequencePermuter.pl — a Perl script that performs the point mutattions at the requested rate.
#		• bowtie2 — a tool for aligning sequencing reads to reference sequences.
#		• samtools — a suite of programs for interacting with high-throughput sequencing data
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


# ------------------------------------------------------- Setup -------------------------------------------------------

PROJECT_DIR='/Path/To/Project/software'

OUTPUT_DIR_READS=$1
PROB=$2

FASTA_READS_1='/Path/To/data/hmp/samples/mid_vaginal/SRS015072/SRS015072-1.fa'
FASTA_READS_2='/Path/To/data/hmp/samples/mid_vaginal/SRS015072/SRS015072-2.fa'

PROC_THREADS=64

echo ""
echo "[" `date '+%m/%d/%y %H:%M:%S'` "]"
echo "[" `date '+%m/%d/%y %H:%M:%S'` "] Starting..."
echo "[" `date '+%m/%d/%y %H:%M:%S'` "]"

#
# Permute the Reads (paired-ends) IFF its not the baseline (0) run
#
if [ $PROB != 0 ]; then
	$PROJECT_DIR/sequencePermuter.pl -p $PROB -t $PROC_THREADS -i $FASTA_READS_1 -o $OUTPUT_DIR_READS
	mv $OUTPUT_DIR_READS/permutedSequences.fasta $OUTPUT_DIR_READS/permutedSequences-1.fasta

	$PROJECT_DIR/sequencePermuter.pl -p $PROB -t $PROC_THREADS -i $FASTA_READS_2 -o $OUTPUT_DIR_READS
	mv $OUTPUT_DIR_READS/permutedSequences.fasta $OUTPUT_DIR_READS/permutedSequences-2.fasta
fi


# --------------------------------------------------- DNA Alignment ---------------------------------------------------

#	SAM Header Location
HEADER='/Path/To/annotations/img/4/bacteria/metadata/header.sam'

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

# 	Output Directory
OUTDIR=$1/bowtie
mkdir -p $OUTDIR

# 	Index - 1.A
FILENAME_1_A="alignments_1_A"
INDEX_1_A=/projects/bioinf/ref/annotations/img/4/bacteria/index/1/A/img_bacteria_1_A

# 	Index - 1.B.1
FILENAME_1_B_1="alignments_1_B_1"
INDEX_1_B_1=/projects/bioinf/ref/annotations/img/4/bacteria/index/1/B/1/img_bacteria_1_B_1

# 	Index - 1.B.2
FILENAME_1_B_2="alignments_1_B_2"
INDEX_1_B_2=/projects/bioinf/ref/annotations/img/4/bacteria/index/1/B/2/img_bacteria_1_B_2


# 	Index - 2.A
FILENAME_2_A="alignments_2_A"
INDEX_2_A=/projects/bioinf/ref/annotations/img/4/bacteria/index/2/A/img_bacteria_2_A

# 	Index - 2.B.1
FILENAME_2_B_1="alignments_2_B_1"
INDEX_2_B_1=/projects/bioinf/ref/annotations/img/4/bacteria/index/2/B/1/img_bacteria_2_B_1

# 	Index - 2.B.2
FILENAME_2_B_2="alignments_2_B_2"
INDEX_2_B_2=/projects/bioinf/ref/annotations/img/4/bacteria/index/2/B/2/img_bacteria_2_B_2


# ------------------------------------------------------- Bowtie ------------------------------------------------------

# 	Parameters for bowtie
SEED_LENGTH=20
SEED_INTERVAL="S,1,0.50"

/nethome/bioinfo/tools/bowtie/v2/2.1.0/bowtie2 --local -D 20 -R 3 -N 0 -L $SEED_LENGTH -i $SEED_INTERVAL -a --time -f -x $INDEX_1_A -1 $READS_1 -2 $READS_2 -S $OUTDIR/$FILENAME_1_A.sam
/nethome/bioinfo/tools/bowtie/v2/2.1.0/bowtie2 --local -D 20 -R 3 -N 0 -L $SEED_LENGTH -i $SEED_INTERVAL -a --time -f -x $INDEX_1_B_1 -1 $READS_1 -2 $READS_2 -S $OUTDIR/$FILENAME_1_B_1.sam
/nethome/bioinfo/tools/bowtie/v2/2.1.0/bowtie2 --local -D 20 -R 3 -N 0 -L $SEED_LENGTH -i $SEED_INTERVAL -a --time -f -x $INDEX_1_B_2 -1 $READS_1 -2 $READS_2 -S $OUTDIR/$FILENAME_1_B_2.sam

/nethome/bioinfo/tools/bowtie/v2/2.1.0/bowtie2 --local -D 20 -R 3 -N 0 -L $SEED_LENGTH -i $SEED_INTERVAL -a --time -f -x $INDEX_2_A -1 $READS_1 -2 $READS_2 -S $OUTDIR/$FILENAME_2_A.sam
/nethome/bioinfo/tools/bowtie/v2/2.1.0/bowtie2 --local -D 20 -R 3 -N 0 -L $SEED_LENGTH -i $SEED_INTERVAL -a --time -f -x $INDEX_2_B_1 -1 $READS_1 -2 $READS_2 -S $OUTDIR/$FILENAME_2_B_1.sam
/nethome/bioinfo/tools/bowtie/v2/2.1.0/bowtie2 --local -D 20 -R 3 -N 0 -L $SEED_LENGTH -i $SEED_INTERVAL -a --time -f -x $INDEX_2_B_2 -1 $READS_1 -2 $READS_2 -S $OUTDIR/$FILENAME_2_B_2.sam


# ------------------------------------------------------- Samtools -----------------------------------------------------

# 	Convert the results for INDEX_1_A, INDEX_1_B_1, INDEX_1_B_2
/nethome/bioinfo/tools/samtools/0.1.19/samtools view -b -S -o $OUTDIR/$FILENAME_1_A.bam $OUTDIR/$FILENAME_1_A.sam
/nethome/bioinfo/tools/samtools/0.1.19/samtools view -b -S -o $OUTDIR/$FILENAME_1_B_1.bam $OUTDIR/$FILENAME_1_B_1.sam
/nethome/bioinfo/tools/samtools/0.1.19/samtools view -b -S -o $OUTDIR/$FILENAME_1_B_2.bam $OUTDIR/$FILENAME_1_B_2.sam

# 	Rehead the results for INDEX_1_A
/nethome/bioinfo/tools/samtools/0.1.19/samtools reheader $HEADER $OUTDIR/$FILENAME_1_A.bam > $OUTDIR/$FILENAME_1_A.reheaded.bam
/nethome/bioinfo/tools/samtools/0.1.19/samtools reheader $HEADER $OUTDIR/$FILENAME_1_B_1.bam > $OUTDIR/$FILENAME_1_B_1.reheaded.bam
/nethome/bioinfo/tools/samtools/0.1.19/samtools reheader $HEADER $OUTDIR/$FILENAME_1_B_2.bam > $OUTDIR/$FILENAME_1_B_2.reheaded.bam


# 	Convert the results for INDEX_2_A, INDEX_2_B_1, INDEX_2_B_2
/nethome/bioinfo/tools/samtools/0.1.19/samtools view -b -S -o $OUTDIR/$FILENAME_2_A.bam $OUTDIR/$FILENAME_2_A.sam
/nethome/bioinfo/tools/samtools/0.1.19/samtools view -b -S -o $OUTDIR/$FILENAME_2_B_1.bam $OUTDIR/$FILENAME_2_B_1.sam
/nethome/bioinfo/tools/samtools/0.1.19/samtools view -b -S -o $OUTDIR/$FILENAME_2_B_2.bam $OUTDIR/$FILENAME_2_B_2.sam

# 	Rehead the results for the 2nd half
/nethome/bioinfo/tools/samtools/0.1.19/samtools reheader $HEADER $OUTDIR/$FILENAME_2_A.bam > $OUTDIR/$FILENAME_2_A.reheaded.bam
/nethome/bioinfo/tools/samtools/0.1.19/samtools reheader $HEADER $OUTDIR/$FILENAME_2_B_1.bam > $OUTDIR/$FILENAME_2_B_1.reheaded.bam
/nethome/bioinfo/tools/samtools/0.1.19/samtools reheader $HEADER $OUTDIR/$FILENAME_2_B_2.bam > $OUTDIR/$FILENAME_2_B_2.reheaded.bam


# ------------------------------------------------- SAM/BAM Processing ------------------------------------------------

# 	Merge the BAM files
/nethome/bioinfo/tools/samtools/0.1.19/samtools merge -f $OUTDIR/alignments-1.bam $OUTDIR/$FILENAME_1_A.reheaded.bam $OUTDIR/$FILENAME_1_B_1.reheaded.bam $OUTDIR/$FILENAME_1_B_2.reheaded.bam
/nethome/bioinfo/tools/samtools/0.1.19/samtools merge -f $OUTDIR/alignments-2.bam $OUTDIR/$FILENAME_2_A.reheaded.bam $OUTDIR/$FILENAME_2_B_1.reheaded.bam $OUTDIR/$FILENAME_2_B_2.reheaded.bam
/nethome/bioinfo/tools/samtools/0.1.19/samtools merge -f $OUTDIR/alignments.bam $OUTDIR/alignments-1.bam $OUTDIR/alignments-2.bam


# Sort the main Alignment
/nethome/bioinfo/tools/samtools/0.1.19/samtools sort $OUTDIR/alignments.bam $OUTDIR/alignments.sorted
/nethome/bioinfo/tools/samtools/0.1.19/samtools index $OUTDIR/alignments.sorted.bam
/nethome/bioinfo/tools/samtools/0.1.19/samtools flagstat $OUTDIR/alignments.sorted.bam > $OUTDIR/flagStats.txt
/nethome/bioinfo/tools/samtools/0.1.19/samtools idxstats $OUTDIR/alignments.sorted.bam > $OUTDIR/reference_Counts.txt

# Sort the child alignments (cuts of main reference)
/nethome/bioinfo/tools/samtools/0.1.19/samtools sort $OUTDIR/$FILENAME_1_A.bam $OUTDIR/$FILENAME_1_A.sorted
/nethome/bioinfo/tools/samtools/0.1.19/samtools index $OUTDIR/$FILENAME_1_A.sorted.bam
/nethome/bioinfo/tools/samtools/0.1.19/samtools flagstat $OUTDIR/$FILENAME_1_A.sorted.bam > $OUTDIR/flagStats-$FILENAME_1_A.txt
/nethome/bioinfo/tools/samtools/0.1.19/samtools idxstats $OUTDIR/$FILENAME_1_A.sorted.bam > $OUTDIR/reference_Counts-$FILENAME_1_A.txt

/nethome/bioinfo/tools/samtools/0.1.19/samtools sort $OUTDIR/$FILENAME_1_B_1.bam $OUTDIR/$FILENAME_1_B_1.sorted
/nethome/bioinfo/tools/samtools/0.1.19/samtools index $OUTDIR/$FILENAME_1_B_1.sorted.bam
/nethome/bioinfo/tools/samtools/0.1.19/samtools flagstat $OUTDIR/$FILENAME_1_B_1.sorted.bam > $OUTDIR/flagStats-$FILENAME_1_B_1.txt
/nethome/bioinfo/tools/samtools/0.1.19/samtools idxstats $OUTDIR/$FILENAME_1_B_1.sorted.bam > $OUTDIR/reference_Counts-$FILENAME_1_B_1.txt

/nethome/bioinfo/tools/samtools/0.1.19/samtools sort $OUTDIR/$FILENAME_1_B_2.bam $OUTDIR/$FILENAME_1_B_2.sorted
/nethome/bioinfo/tools/samtools/0.1.19/samtools index $OUTDIR/$FILENAME_1_B_2.sorted.bam
/nethome/bioinfo/tools/samtools/0.1.19/samtools flagstat $OUTDIR/$FILENAME_1_B_2.sorted.bam > $OUTDIR/flagStats-$FILENAME_1_B_2.txt
/nethome/bioinfo/tools/samtools/0.1.19/samtools idxstats $OUTDIR/$FILENAME_1_B_2.sorted.bam > $OUTDIR/reference_Counts-$FILENAME_1_B_2.txt


/nethome/bioinfo/tools/samtools/0.1.19/samtools sort $OUTDIR/$FILENAME_2_A.bam $OUTDIR/$FILENAME_2_A.sorted
/nethome/bioinfo/tools/samtools/0.1.19/samtools index $OUTDIR/$FILENAME_2_A.sorted.bam
/nethome/bioinfo/tools/samtools/0.1.19/samtools flagstat $OUTDIR/$FILENAME_2_A.sorted.bam > $OUTDIR/flagStats-$FILENAME_2_A.txt
/nethome/bioinfo/tools/samtools/0.1.19/samtools idxstats $OUTDIR/$FILENAME_2_A.sorted.bam > $OUTDIR/reference_Counts-$FILENAME_2_A.txt

/nethome/bioinfo/tools/samtools/0.1.19/samtools sort $OUTDIR/$FILENAME_2_B_1.bam $OUTDIR/$FILENAME_2_B_1.sorted
/nethome/bioinfo/tools/samtools/0.1.19/samtools index $OUTDIR/$FILENAME_2_B_1.sorted.bam
/nethome/bioinfo/tools/samtools/0.1.19/samtools flagstat $OUTDIR/$FILENAME_2_B_1.sorted.bam > $OUTDIR/flagStats-$FILENAME_2_B_1.txt
/nethome/bioinfo/tools/samtools/0.1.19/samtools idxstats $OUTDIR/$FILENAME_2_B_1.sorted.bam > $OUTDIR/reference_Counts-$FILENAME_2_B_1.txt

/nethome/bioinfo/tools/samtools/0.1.19/samtools sort $OUTDIR/$FILENAME_2_B_2.bam $OUTDIR/$FILENAME_2_B_2.sorted
/nethome/bioinfo/tools/samtools/0.1.19/samtools index $OUTDIR/$FILENAME_2_B_2.sorted.bam
/nethome/bioinfo/tools/samtools/0.1.19/samtools flagstat $OUTDIR/$FILENAME_2_B_2.sorted.bam > $OUTDIR/flagStats-$FILENAME_2_B_2.txt
/nethome/bioinfo/tools/samtools/0.1.19/samtools idxstats $OUTDIR/$FILENAME_2_B_2.sorted.bam > $OUTDIR/reference_Counts-$FILENAME_2_B_2.txt


# ------------------------------------------------------ Clean Up -----------------------------------------------------

rm -fR $OUTDIR/$FILENAME_1_A.sam
rm -fR $OUTDIR/$FILENAME_1_A.bam
rm -fR $OUTDIR/$FILENAME_1_B_1.sam
rm -fR $OUTDIR/$FILENAME_1_B_1.bam
rm -fR $OUTDIR/$FILENAME_1_B_2.sam
rm -fR $OUTDIR/$FILENAME_1_B_2.bam

rm -fR $OUTDIR/$FILENAME_2_A.sam
rm -fR $OUTDIR/$FILENAME_2_A.bam
rm -fR $OUTDIR/$FILENAME_2_B_1.sam
rm -fR $OUTDIR/$FILENAME_2_B_1.bam
rm -fR $OUTDIR/$FILENAME_2_B_2.sam
rm -fR $OUTDIR/$FILENAME_2_B_2.bam

rm -fR $OUTDIR/alignments-1.bam
rm -fR $OUTDIR/alignments-2.bam

rm -fR $OUTDIR/$FILENAME_1_A.reheaded.bam
rm -fR $OUTDIR/$FILENAME_1_B_1.reheaded.bam
rm -fR $OUTDIR/$FILENAME_1_B_2.reheaded.bam
rm -fR $OUTDIR/$FILENAME_2_A.reheaded.bam
rm -fR $OUTDIR/$FILENAME_2_B_1.reheaded.bam
rm -fR $OUTDIR/$FILENAME_2_B_2.reheaded.bam

rm -fR $OUTDIR/alignments.bam
rm -fR $OUTDIR/alignments-MAPQ30.bam

# -------------------------------------------------------- End --------------------------------------------------------

echo "[" `date '+%m/%d/%y %H:%M:%S'` "]"
echo "[" `date '+%m/%d/%y %H:%M:%S'` "]"
echo "[" `date '+%m/%d/%y %H:%M:%S'` "] Done."
echo "[" `date '+%m/%d/%y %H:%M:%S'` "]"
echo ""