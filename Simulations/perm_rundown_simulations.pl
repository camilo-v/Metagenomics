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
#	The purpose of this program is to cycle through the permutation directories and summarize the various
#   alignment hits (counts) for each reference sequence present in the SAM/BAM files.  This program is called from each
#	of the "count_rundown_TEMPLATE.sh" scripts.
#
#
#   NOTES:
#   The program takes a range of permutation directories to process, rather than trying to 'slurp' all the permutation
#	directories at once.  Please see the dependencies section below for the required modules.
#
#
#   DEPENDENCIES:
#
#       • Bio::Perl
#		• Bio::Seq
#       • File::Slurp
#       • Getopt::Long
#		• Memory::Usage
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
use File::Slurp qw(read_dir);
use Bio::Perl;
use Bio::Seq;
use Bio::SeqIO;
use Memory::Usage;

# ------------------------------------------------------ Main ---------------------------------------------------------

my $file;               # Input File (reference sequence names)
my $outputDir;			# Output Directory
my $permutationsDir;    # Directory holding permutations
my $refSeqs;            # FASTA reference sequences
my $prefix;             # Prefix string for output files
my $read_hits = 0;      # Flag to print the Read-Hit distribution (time consuming)
my $low;                # Low bound for directory slurping
my $high;               # High bound for directory slurping
my $mapq;               # MAPQ threshold to use for filtering

GetOptions ( 'i=s' => \$file, 'o:s' => \$outputDir, 'p=s' => \$permutationsDir, 'f=s' => \$refSeqs, 'prefix:s' => \$prefix,
                'read_hits' => \$read_hits, 'low:i' => \$low, 'high:i' => \$high, 'mapq:i' => \$mapq );

if( ! defined( $file)  || ! defined($permutationsDir) ) {
    print( RED, "\n\nWARNING! Missing Files ", RESET ); print( ".. $!\n\n" ); exit(1);
}

chomp($file);
chomp($permutationsDir);
chomp($refSeqs);
chomp($low);
chomp($high);
chomp($mapq);

if( defined( $prefix ) )
{ chomp( $prefix ); }

my $currentWD = getcwd();       # Get a pointer to the Current Working Directory

# Default to the current working directory if none provided
if( !defined( $outputDir ) )
{ $outputDir = $currentWD; }

$| = 1;  # Flush STDOUT

# ------------------------------------------------ Output Files ---------------------------------------------------------

my $fout1 = $outputDir . "/" . $prefix . "_counts-" . $low . "-" . $high . ".txt";
unless( open( OUTFILE, ">$fout1" ) ) { print "File $fout1 does not exist"; exit; }

if( $read_hits ) {
    my $fout2 = $outputDir . "/" . $prefix . "_read_hit_distribution-" . $low . "-" . $high . ".txt";
    unless( open( OUTFILE2, ">$fout2" ) ) { print "File $fout2 does not exist"; exit; }
}

# ------------------------------------------ Reference Name File Loader --------------------------------------------------

print( BOLD, GREEN, "\n--\n", RESET );
print( BOLD, "Loading file(s) ... ", RESET );

unless( open( INFILE, $file ) )
{ die( "\n\nUnable to open Input File: $file.  Please verify.\n\n" ); }

my %referenceNames = ();

my $numberOfLines = 0;

while( my $line = <INFILE> ) {

    chomp($line);
	$line   =~ s/\r//g;
    my @lineArray = split(/ /,$line);

	my $aRefName = $lineArray[ 0 ];	# Modified for MetaSim and Bowtie2 Indexing

    $referenceNames{ $aRefName } = $lineArray[ 0 ];

    $numberOfLines++;
}

close(INFILE);

print( CYAN, "\n Number of Reference Names in input file: ", RESET );
print( " " . &addCommas($numberOfLines) . "\n" );


# ---------------------------------------------- Permutation Processor ---------------------------------------------------

print( BOLD, GREEN, "\n--\n", RESET );
print( BOLD, "Looking at permutations in dir: " . $permutationsDir . "", RESET );

# Hash for the reference counts,
# KEY: reference ID
# VAL: Nested Hash --> KEY: Permutation Number, VAL: Array with counts
my %referenceCounts = ();

my %permutationColumns = ();

my $root = $permutationsDir;

my @alignmentFiles = (  "alignments.sorted.bam" );

# Hash to keep track of how many genomes a read maps to
# KEY: Read ID
# VAL: Count of how many genomes KEY maps to
my %readHitDistribution = ();

for my $dir ( grep { -d "$root/$_" } read_dir( $root ) ) {

	# Trying to load all the reads for all the permutations into memory is not optimal.  So we split the
	# permutation directories into batches — 'low' and 'high' demarcate the permutation interval to use.
    my $dirNumber = $dir;
    $dirNumber =~ s/perm_//g;

    if( ( $low <= $dirNumber ) && ( $dirNumber <= $high ) ) {

        $permutationColumns{ $dir } = 1;

        my $bowtie_dir  = $permutationsDir . "/" . $dir . "/bowtie";

        foreach my $aFile ( @alignmentFiles ) {

            my $aln_file = $bowtie_dir . "/$aFile";

            print( GREEN, "\n  Processing $aln_file ...", RESET );

            my @alignments  = `samtools view $aln_file`;
            my $alnCount    = @alignments;

            print( "\n  - Dealing with $alnCount ..." );

            my $totalRefHits = 0;

            my $readsPassingFilter = 0;

            foreach my $alignment ( @alignments ) {

                my @alnFields    = split(/\t/,$alignment);
                my $referenceHit = $alnFields[ 2 ];
                my $readName     = $alnFields[ 0 ];
                my $readQuality  = $alnFields[ 4 ];

                # We are only interested in "high" quality alignments -- for some definition of "high"
                # Here, our definition is "high enough" and it is set with the --mapq parameter
                if ( $readQuality > $mapq ) {

                    $readsPassingFilter++;

                    if( defined( $referenceNames{ $referenceHit } ) ) {

                        # Number of read hits a reference has
                        # Reference --> No. Reads
                        push( @{ $referenceCounts{ $referenceHit }{ $dir } }, $readName );

                        # Read --> No. of References
                        my $count = 1;

                        if( defined( $readHitDistribution{ $readName } ) ) {
                            my $oldCount = $readHitDistribution{ $readName }{ $dir };
                            $oldCount++;
                            $readHitDistribution{ $readName }{ $dir } = $oldCount;
                        }

                        if( ! defined( $readHitDistribution{ $readName } ) ) {
                            $readHitDistribution{ $readName }{ $dir } = $count;
                        }
                    }

                } # End MAPQ filter

            } # End Alignment loop

            print( "\n  - Passing $mapq MAPQ Filter: $readsPassingFilter ..." );

            print( "\n" );
        }

    }

}


# --------------------------------------------- Fractional Counts ---------------------------------------------
#
#   Maps a read to the number (int) of genomes it is mapped to within a given permutation.
#

if( $read_hits ) {

    print( BOLD, GREEN, "\n--\n", RESET );
    print( BOLD, "Generating Read-Hit Distribution...\n", RESET );

    my %readGenomicIncidence = ();

    for my $aRead ( sort{ $a cmp $b || $a <=> $b } keys %readHitDistribution ) {

        print OUTFILE2 ( "$aRead\t" );

        my $cumulativeHit = 0;

        my $numberOfPermutations += scalar keys %permutationColumns;

        for my $aCol ( sort{ $a cmp $b || $a <=> $b } keys %permutationColumns ) {

            if( defined( $readHitDistribution{ $aRead }{ $aCol } ) ) {

                my $numberOfGenomeHits = $readHitDistribution{ $aRead }{ $aCol };

                $cumulativeHit += $numberOfGenomeHits;

                $readGenomicIncidence{ $aRead }{ $aCol } = $numberOfGenomeHits;

                print OUTFILE2 ( "$numberOfGenomeHits\t" );
            }

            if( ! defined( $readHitDistribution{ $aRead }{ $aCol } ) ) {
                print OUTFILE2 ( "0\t" );
            }
        }

        my $averageNumberOfHitsAccrossPermutations = $cumulativeHit / $numberOfPermutations;

        print OUTFILE2 ( "$averageNumberOfHitsAccrossPermutations\n" );
    }

    close(OUTFILE2);
}



# ------------------------------------------- Sequence Lengths -------------------------------------------
#
#	Reference Sequence File Loader (used for sequence lengths in reporting files)
#

print( BOLD, GREEN, "\n--\n", RESET );
print( BOLD, "Loading FASTA reference sequences ... ", RESET );

my $seq_in = Bio::SeqIO->new( -format => 'fasta', -file => $refSeqs );

my $refsMatched = 0;

my %seqLengths = ();

while( my $seq = $seq_in->next_seq() ) {

    if( defined( $referenceNames{ $seq->id } ) ) {
        $refsMatched++;
        $seqLengths{ $seq->id } = $seq->length;
    }
}

my $sizeOfGenSeq = keys( %seqLengths );

print( CYAN, "\n Sequences matched with Reference List: ", RESET );
print( " " . &addCommas($sizeOfGenSeq) . "(" . &addCommas($refsMatched) . ")\n" );



# ------------------------------------------- Output File Reports -------------------------------------------
#
#	Write the reference count data structures to the output files
#

print( BOLD, GREEN, "\n--\n", RESET );
print( BOLD, "Creating reports...\n", RESET );

print( BOLD, "Printing Column Headers...\n", RESET );

my $headerCounter = 0;

# First loop is just to pickup the correct column headers
for my $referenceSeq ( sort{ $a cmp $b || $a <=> $b } keys %referenceCounts ) {

    if( $headerCounter == 1 ) {
        last; }

    print OUTFILE ( "Reference Sequence" . "\t" . "Length" . "" );

    for my $aCol ( sort{ $a cmp $b || $a <=> $b } keys %permutationColumns ) {
        print OUTFILE ( "\t" . $aCol . "" );
    }

    $headerCounter++;
}

print( BOLD, "Printing Data Columns (reference sequence rows)...", RESET );

# Second loop is to fill out the data columns
for my $referenceSeq ( sort{ $a cmp $b || $a <=> $b } keys %referenceNames ) {

    print OUTFILE ( "\n" . $referenceNames{ $referenceSeq } . "\t" );

    my $refLength = 0;

    if( defined( $seqLengths{ $referenceSeq } ) ) {

        $refLength = $seqLengths{ $referenceSeq };
        print OUTFILE ( "" . $refLength . "" );
    }

    # $aCol --> is a perm_* directory
    for my $aCol ( sort{ $a cmp $b || $a <=> $b } keys %permutationColumns ) {

        my $hits = 0;

        if( defined( $referenceCounts{ $referenceSeq } ) ) {

            # $referenceSeq --> is a $referenceHit name from samtools field 2
            my %permutations = %{ $referenceCounts{ $referenceSeq } };

            if( defined( $permutations{ $aCol } ) ) {

                my @arrayOfHits  = @{ $permutations{ $aCol } };

                foreach my $aHit ( @arrayOfHits ) {

                	my $numberOfGenomesReadMapsTo = $readHitDistribution{ $aHit }{ $aCol };

                	$hits += ( 1 / $numberOfGenomesReadMapsTo );
                }
            }
        }

        print OUTFILE ( "\t" .$hits . "" );
    }

}

print("\n\n");

close(OUTFILE);

$| = 1; # Re-flush STDOUT

# ---------------------------------------------- Functions & Methods -------------------------------------------------


#	addCommas()
#	Simple method to format a number with commas
#
#
sub addCommas
{
    local($_) = shift;
    1 while s/^(-?\d+)(\d{3})/$1,$2/;
    return $_;
}


# Exit the program with a normal condition flag
exit( 0 );