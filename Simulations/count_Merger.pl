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
#	The purpose of this program is to paste (in parallel) a set of columns from a set of tab-delimited files.  The
#	purpose of this program is similar to the Unix 'paste' utility.  The difference is that this program pads empty
#	columns to 'square' the output file — something the 'paste' utility does not do.
#
#
#   NOTES:
#   Please see the dependencies section below for the required libraries (if any).
#
#   DEPENDENCIES:
#
#       • None. Uses standard Perl modules.
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
use POSIX qw/strftime/;

#############################################################################################################
#
#													Main
#
#-----------------------------------------------------------------------------------------------------------
#
#	Preliminaries
#

my $bundle = "0.1";
my $build  = "A0001";

my $file;               # Input File list (the files that will be merged)
my $outputDir;			# Output Directory
my $prefix;             # Prefix for output file
my $last_column = 0;    # Omit the last column in the merged files (read-hit distribution files have it)
my $fileHeader = 0;     # Flag to skip the file header
my $file_type;          # Files to merged can be either 'counts' or 'hit distributions'

GetOptions (    'i=s' => \$file, 'o:s' => \$outputDir, 'p=s' => \$prefix, 'type=s' => \$file_type,
                'omit' => \$last_column, 'header' => \$fileHeader );

if( ! defined( $file) )
{ print( RED, "\n\nWARNING! Missing Files ", RESET ); print( ".. $!\n\n" ); exit(1); }

chomp($file);

my $currentWD = getcwd();       # Get a pointer to the Current Working Directory

# Default to the current working directory if none provided
if( !defined( $outputDir ) )
{ $outputDir = $currentWD; }

if( !defined( $prefix ) )
{ $prefix = "prefixDefault"; }

if( !defined( $file_type ) )
{ $file_type = "typeDefault"; }

my $startingPosition = 1;

my $flagForFileTypeOfCounts = 0;

if( defined( $file_type ) ) {
	if( $file_type eq "counts" ) {
		$startingPosition = 2;
		$flagForFileTypeOfCounts = 1;
	}
}

$| = 1;  # Flush STDOUT

print( "\n" );
print( "[" . strftime('%d-%b-%Y %H:%M:%S',localtime()) . "] Starting...\n" );

#-----------------------------------------------------------------------------------------------------------
#
#	Output Files with aggregated counts
#

print( "[" . strftime('%d-%b-%Y %H:%M:%S',localtime()) . "] Creating output files...\n" );

my $fout1 = $outputDir . "/" . $prefix . "_" . $file_type . ".txt";
unless( open( OUTFILE, ">$fout1" ) ) { print "File $fout1 does not exist"; exit; }

#-----------------------------------------------------------------------------------------------------------
#
#	File Loader
#

print( "[" . strftime('%d-%b-%Y %H:%M:%S',localtime()) . "] Loading file list...\n" );

unless( open( INFILE, $file ) )
{ die( "\n\nUnable to open List Input File: $file.  Please verify.\n\n" ); }

# Maps a file ID (int) to the path of the file to merge
# KEY: File ID
# VAL: Path of file in KEY
my %filesToAggregate = ();

my $numberOfLines = 0;

while( my $line = <INFILE> )
{
    chomp($line);
	$line   =~ s/\r//g;
    my @lineArray = split(/\t/,$line);

    my $filePath = $lineArray[ 0 ];

    $filesToAggregate{ $numberOfLines } = $filePath;

    $numberOfLines++;
}

close(INFILE);

print( "[" . strftime('%d-%b-%Y %H:%M:%S',localtime()) . "] Number of lines in file list for Q" . $prefix . ":" );
print( " " . &addCommas($numberOfLines) . "\n" );

# Maps a fileName to its contents (nested array of lines)
# KEY: File Name (string)
# VAL: Nested array of lines
#my %fileContents = ();

# Maps a FileID to the FileName (column number in resulting merged file)
my %fileToID = ();

# Maps a Sequence Name to a nested hash of its fileName and array of permutations
my %referencePermutationsCounts = ();

# Maps a Sequence Name to its length
my %referenceSequenceLength = ();

my $totalNumberOfColumns = 0;

for my $fileID ( sort{ $a cmp $b || $a <=> $b } keys %filesToAggregate )
{
    my $filePath = $filesToAggregate{ $fileID };

    my @fileIDArr = split(/\//,$filePath);
    my $fileName = $fileIDArr[ $#fileIDArr ];  # The last item in the array is the file name

    # Read the file contents
    unless( open( INFILE2, $filePath ) )
    { die( "\n\nUnable to open Input File in loop: $filePath.  Please verify.\n\n" ); }

#    my @fileLineArray = ();

    my $numberOfLinesInFile = 1;

    if( $fileHeader ) { # All this causes is to move the INFILE2 pointer ahead by one cycle
        my $fileHeader = <INFILE2>;
    }

    my $flagForLastColRead = 1;

    print( "[" . strftime('%d-%b-%Y %H:%M:%S',localtime()) . "] Reading in file $fileID...\n" );

    while( my $line = <INFILE2> )
    {
        chomp($line);
        $line   =~ s/\r//g;
        my @lineArray = split(/\t/,$line);

        my $sequenceName = $lineArray[ 0 ];
        my $sequenceLength = $lineArray[ 1 ];

        my @permutationsArray = ();

        my $lastColumnIndex = $#lineArray;

        # Read-Hit distribution files have a last column with the average number of hits
        if( $last_column ) {
            $lastColumnIndex = $#lineArray - 1;
        }

        if( $flagForLastColRead ) {
            $totalNumberOfColumns += $lastColumnIndex;
            $flagForLastColRead = 0;
        }

        for my $n ( $startingPosition .. $lastColumnIndex ) {
            push( @permutationsArray, $lineArray[ $n ] );
        }

#        push( @fileLineArray, $line );

        $numberOfLinesInFile++;

        $referencePermutationsCounts{ $sequenceName }{ $fileName } = \@permutationsArray;

        $referenceSequenceLength{ $sequenceName } = $sequenceLength;
    }

#    my $numberOfLinesInFile = @fileLineArray;

    print( "[" . strftime('%d-%b-%Y %H:%M:%S',localtime()) . "] $fileID) $fileName: " . &addCommas($numberOfLinesInFile) . " lines\n" );

#    $fileContents{ $fileName } = \@fileLineArray;

    $fileToID{ $fileID } = $fileName;

    close(INFILE2);
}

print( "[" . strftime('%d-%b-%Y %H:%M:%S',localtime()) . "] Total Columns: $totalNumberOfColumns\n" );

print( "[" . strftime('%d-%b-%Y %H:%M:%S',localtime()) . "] Preparing Output Files...\n" );

my $numberOfColumnsToPrint += scalar keys %fileToID;

# Print the output file header
print OUTFILE ( "Reference Sequence\t" );

if( $flagForFileTypeOfCounts ) {
	print OUTFILE ( "Reference Length\t" );
}

for my $n ( 1 .. $totalNumberOfColumns ) {
    print OUTFILE ( "perm_" . $n . "\t" );
}

print OUTFILE ( "Average Hits\n" );

print( "[" . strftime('%d-%b-%Y %H:%M:%S',localtime()) . "] Generating Reports...\n" );

for my $reference ( sort{ $a cmp $b || $a <=> $b } keys %referencePermutationsCounts )
{
    print OUTFILE ( $reference . "\t" );

    if( $flagForFileTypeOfCounts ) {
		print OUTFILE ( $referenceSequenceLength{ $reference } . "\t" );
    }

    my %innerHash = %{ $referencePermutationsCounts{ $reference } };

    my $totalHitCount = 0;

    my $colsPopulated = 0;

    # This loop merges through all the batches
    for my $key ( sort{ $a cmp $b || $a <=> $b } keys %innerHash )
    {
        my @innerArray       = @{ $innerHash{ $key } };
        my $sizeOfInnerArray = @innerArray;

        for my $n ( 0 .. $#innerArray ) {

            my $aCount = 0;

            if( defined( $innerArray[ $n ] ) ) {

                $aCount = $innerArray[ $n ];

                print OUTFILE ( $aCount . "\t" );

                $colsPopulated++;

            }

            $totalHitCount += $aCount;
        }
    }

    my $colsToPad = ( $totalNumberOfColumns - $colsPopulated ) - 1; # adjust for read name

    for my $i ( 0 .. $colsToPad ) {
        print OUTFILE ( "0\t" );
    }


    my $averageNumberOfHits = $totalHitCount / $totalNumberOfColumns;

    if( $last_column ) {
        print OUTFILE ( "$averageNumberOfHits\t" );
    }

    print OUTFILE ( "\n" );

    # Because the Hash is pretty big, and we only have 32GB to work with (only!?)
    # we'll need to delete the key-value pair after it has been printed
    # This ensures that we'll be able to make it through the hash-inner-hash structures without
    # having the kernel kill our process.
    delete( $referencePermutationsCounts{ $reference } );
}

print( "[" . strftime('%d-%b-%Y %H:%M:%S',localtime()) . "] Done.\n" );


close(OUTFILE);

print("\n\n");

# Re-flush STDOUT
$| = 1;



#-----------------------------------------------------------------------------------------------------------
#
#                                               End of Line
#
#############################################################################################################
#
#	addCommas() Method
#
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