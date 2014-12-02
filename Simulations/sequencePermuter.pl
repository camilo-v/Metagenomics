#!/usr/bin/perl -w

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
#	The purpose of this program is to permute a particular sequence using a uniform distribution.  The script mutates
#	the reads by walking along the length of each mate read, randomly drawing from a uniform distribution, and
#	comparing the result of the draw to the requested mutation rate.  If the result is less thant the mutation rate,
#	then the nucleotide is "mutated".
#	Nucleotides are randomly mutated with other (non-self) nucleotides by randomly drawing from the remaining three (3)
#	nucleotides, and Chargaff's base-pairing rules are not observed in the mutation step.  Once all the NGS reads
#	have been processed, they are saved to disk in FASTQ format.
#
#
#   NOTES:
#   A paired-end fragment consists of two (2) mates: the left read and the right read.  Independent mutataion steps are
#	carried out for each mate in the fragment.
#
#
#   DEPENDENCIES:
#
#       • Getopt::Long
#		• Bio::Perl
#       • Cwd
#       • IO::Pipe
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

# Perl Modules
use strict;
use Getopt::Long;
use Term::ANSIColor qw(:constants);
use Benchmark;
use Cwd;
use Bio::Perl;
use Bio::Seq;
use Bio::SeqIO;
use Bio::Seq::Quality;
use IO::Pipe;

# ------------------------------------------------------ Main ---------------------------------------------------------

my $file;               # Input File to be disturbed (FASTA format, reads or genomes)
my $outputDir;			# Output Directory
my $numThreads;         # Number of processor threads to use
my $probability;        # Target Probability of event (mutation) occurying


GetOptions ( 'i=s' => \$file, 'o:s' => \$outputDir, 't:s' => \$numThreads, 'p:i' => \$probability );

if( ! defined( $file) )
{ print( RED, "\n\nWARNING! Missing Files ", RESET ); print( ".. $!\n\n" ); exit(1); }

chomp($file);

my $currentWD = getcwd();       # Get a pointer to the Current Working Directory

# Default to the current working directory if none provided
if( !defined( $outputDir ) )
{ $outputDir = $currentWD; }

# Default to using only the main thread
if( !defined( $numThreads ) )
{ $numThreads = 1; }

# Default probability is 0.02%
if( !defined( $probability ) )
{ $probability = 0.02; }

$probability = ( $probability / 100 );

print( YELLOW, "\nProb: $probability", RESET );


$| = 1;  # Flush STDOUT

# -------------------------------------------------- Output Files -----------------------------------------------------

my $fout1 = $outputDir . "/" . "details.txt";
unless( open( OUTFILE, ">$fout1" ) ) { print "File $fout1 does not exist"; exit; }


# ------------------------------------------------- Sequenc Loader ----------------------------------------------------
#
#	Load sequences from file into a BioPerl object
#

print( BOLD, GREEN, "\n--\n", RESET );
print( BOLD, "Loading file(s) ... \n", RESET );

# Bioperl object to read-in the sequences.
my $seq_in = Bio::SeqIO->new('-file' => "<$file", '-format' => 'fasta');

#
# Array to hold sequences (makes it easier to loop & manipulate)
#
my @seq_array = ();

print( BOLD, GREEN, "\n--\n", RESET );
print( BOLD, "Caching sequences ... \n", RESET );

# Put the sequences into an array
while( my $seq = $seq_in->next_seq() )
{ push( @seq_array, $seq ); }

print( BOLD, GREEN, "\n--\n", RESET );
print( BOLD, "Sorting sequences ... \n", RESET );

#
# Sort the sequences by length
#
@seq_array = sort { $a->length <=> $b->length } @seq_array;

print( BOLD, GREEN, "\n--\n", RESET );
print( BOLD, "Gathering descriptive stats ... \n", RESET );

# Descriptive stats on the sequences
my $totalLength = 0;
my $sequenceCount = 0;

foreach my $seq ( @seq_array )
{
    $totalLength += $seq->length();
    $sequenceCount++;
}

my $mean_length = ( $totalLength / $sequenceCount );
my $median_length = $seq_array[ $sequenceCount / 2 ]->length();



print OUTFILE ( "Details for:\n$file" );
print OUTFILE ( "\n\nOutput Dir:\n$outputDir\n\n" );

print( CYAN, "\n Number of Sequences in input file: ", RESET );
print( " " . &addCommas($sequenceCount) . "" );
print OUTFILE ( "\n Number of Sequences in input file: " . &addCommas($sequenceCount) . "" );

print( CYAN, "\n Total Length of Sequences: ", RESET );
print( " " . &addCommas($totalLength) . "" );
print OUTFILE ( "\n Total Length of Sequences: "  . &addCommas($totalLength) . "" );

print( CYAN, "\n Mean Length of Sequences: ", RESET );
print( " " . &addCommas($mean_length) . "" );
print OUTFILE ( "\n Mean Length of Sequences: "  . &addCommas($mean_length) . "" );

print( CYAN, "\n Median Length of Sequences: ", RESET );
print( " " . &addCommas($median_length) . "\n" );
print OUTFILE ( "\n Median Length of Sequences: "  . &addCommas($median_length) . "" );


print( BOLD, GREEN, "\n--\n", RESET );
print( BOLD, "Permuting Sequences ... \n", RESET );


# Output FASTA file with permuted sequences
my $outputFastaFile = $outputDir . "/" . "permutedSequences.fasta";
my $seqout = Bio::SeqIO->new( -file   => ">$outputFastaFile", -format => 'fasta' );


# Hash with Swapping Table
my %swapTable = (   'A' => [ 'T', 'G', 'C' ],
                    'T' => [ 'A', 'G', 'C' ],
                    'G' => [ 'A', 'T', 'C' ],
                    'C' => [ 'A', 'T', 'G' ],
                    'N' => [ 'N', 'N', 'N' ],
                    'M' => [ 'M', 'M', 'M' ]
                );



my $processedReads = 0;


# Perl Block for parallel code
{
	my $objects = int(@seq_array / $numThreads);

	my @subsets = map { [splice @seq_array, 0, $objects] } 0 .. $numThreads-1;

	push @{ $subsets[$_] }, $seq_array[$_] for 0 .. $#seq_array;

	my $pipe = new IO::Pipe; 	# Data channel pipe

	my @workers = ();        	# Array to store worker thread PID's

	for my $i ( 0 .. $numThreads - 1 )
	{
		my $pid = fork();

		if( $pid )
		{
			# Parent process, $pid is the child's process ID
			push( @workers, $pid );
		}
		elsif( $pid == 0 )
		{
			$pipe->writer();
			$pipe->autoflush(1);

			for my $seq ( @{ $subsets[ $i ] } )
			{
				# Do some work on the sequence objects
				my $sequenceID      = $seq->id();
                my $sequence        = $seq->seq();
                my $sequenceLength  = $seq->length();
                #my $sequenceQual	= $seq->qual();

                my $p = $probability;
                my $n = $sequenceLength;
                my $k = sprintf( "%.3f", ( $p * $n ) );

                my $x;

                my @nucleotideArray = split(//,"$sequence");

                my $numberOfSwaps = 0;

                my $newPermutedSequence = "";

                foreach my $nt ( @nucleotideArray )
                {
                    if( defined( $swapTable{ $nt } ) )
                    {
                        my @swapArr     = @{ $swapTable{ $nt } };
                        my $toSwapWith  = $swapArr[ int( rand( 3 ) ) ];

                        my $x = sprintf( "%.3f", ( rand() ) );

                        if( $x < $p )
                        {
                            $newPermutedSequence = $newPermutedSequence . $toSwapWith;
                            $numberOfSwaps++;
                        }

                        if( $x >= $p )
                        {
                            $newPermutedSequence = $newPermutedSequence . $nt;
                        }

                    } # defined nucleotide

                } # Loop through the nucleotide array

                my $descriptionString = "| p: $p, n: $n, k: $k, swaps: $numberOfSwaps";

                # We send back the results along with some differentiating ID
				print $pipe ( "1;$newPermutedSequence;$sequenceID;$descriptionString" . "\n" );

			}

			exit(0);

		} # End Worker Thread
		else
		{
			die("Can't fork worker: $!\n");
		}

	} # End thread-spawn loop

	# Read data back from the workers as it becomes available
	$pipe->reader();

	while( my $dataReceived = <$pipe> )
	{
		# Data is available in $dataReceived variable
		chomp( $dataReceived );

        # The data read back from the child processes is available in the following array
        my @permSeqData = split(/\;/,$dataReceived);

        if( $permSeqData[ 0 ] == 1 )
        {
            $processedReads++;

            my $permutedSeq     = $permSeqData[ 1 ];
            my $permutedID      = $permSeqData[ 2 ];
            my $permutedDesc    = $permSeqData[ 3 ];

			# Artificial Qualities
			my $qual = '0';
			my $qual_string = "";
			for my $n ( 0 .. length( $permutedSeq ) )
			{ $qual_string = $qual_string . $qual; }

            my $new_seq_obj = Bio::Seq->new(-seq => "$permutedSeq", -display_id => "$permutedID", -desc => "$permutedDesc" );

            # Write it out to the FASTA file
            $seqout->write_seq( $new_seq_obj );
        }

	}

	# Reap the worker threads, i.e., wait for all the processes to end
	foreach( @workers )
	{ waitpid( $_, 0 ); }

}

print( CYAN, "\n Number of Sequences permuted: ", RESET );
print( " " . &addCommas($processedReads) . "" );
print OUTFILE ( "\n Number of Sequences permuted: " . &addCommas($processedReads) . "" );

print("\n\n");

# Re-flush STDOUT
$| = 1;

close(OUTFILE);


# ----------------------------------------------- Functions & Methods -------------------------------------------------


#	addCommas
#	Simple method to format a number with commas
#
sub addCommas
{
    local($_) = shift;
    1 while s/^(-?\d+)(\d{3})/$1,$2/;
    return $_;
}


# Exit the program with a normal condition flag
exit( 0 );