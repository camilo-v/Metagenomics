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
#	The purpose of this program is to summarize the pvalues of a set of baseline strains (HMP) from the
#   different genomic sequences.
#
#   NOTES:
#   Please see the dependencies section below for the required libraries (if any).
#
#
#   DEPENDENCIES:
#
#       • Cwd.
#       • List::Util
#       • Benchmark
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
use List::Util ();

# ------------------------------------------------------ Main ---------------------------------------------------------

my $file;               # Input File with P-values
my $outputDir;			# Output Directory
my $columnToUse;		# The Column in which the pvalue is found at
my $delimiter;			# Delimiter pattern to split file with (default is TAB)

GetOptions ( 'i=s' => \$file, 'o:s' => \$outputDir, 'c:i' => \$columnToUse, 'd:s' => \$delimiter );

if( ! defined( $file) )
{ print( RED, "\n\nWARNING! Missing Files ", RESET ); print( ".. $!\n\n" ); exit(1); }

chomp($file);

my $currentWD = getcwd();       # Get a pointer to the Current Working Directory

# Default to the current working directory if none provided
if( !defined( $outputDir ) ) {
	$outputDir = $currentWD; }

if( !defined( $columnToUse ) ) {
	$columnToUse = 1; }

if( !defined( $delimiter ) ){
	$delimiter = "\t";
}

if( defined( $delimiter ) ) {
	chomp($delimiter);
	$delimiter = "" . $delimiter . "";
}

$| = 1;  # Flush STDOUT

# -------------------------------------------------- Output Files -----------------------------------------------------

my $fout1 = $outputDir . "/" . "pvals-composite.txt";
unless( open( OUTFILE, ">$fout1" ) ) { print "File $fout1 does not exist"; exit; }


# --------------------------------------------------- File Loader -----------------------------------------------------

print( BOLD, GREEN, "\n--\n", RESET );
print( BOLD, "Loading file(s) ... ", RESET );

unless( open( INFILE, $file ) )
{ die( "\n\nUnable to open Input File: $file.  Please verify.\n\n" ); }

my %pvals = ();

my $numberOfLines = 0;

my $header = <INFILE>;

while( my $line = <INFILE> )
{
    chomp($line);
	$line   =~ s/\r//g;

	#
	# WARNING!
	# This line is tricky, when using a tab (\t) delimiter the split function should be called with
	# enclosing //, but when its a comma, it shouldn't.... just stick to TABS and DO NOT USE COMMAS!!!
	# The Strain names have commas in them, so avoid them as a delimiter.
	#
    my @lineArray = split(/$delimiter/,$line);

    my $strain_line = $lineArray[ 0 ];
    my $pval   		= $lineArray[ $columnToUse ];

    $strain_line =~ s/"//g;

    my @strainLineArray = split(/ /,$strain_line);

    my $strain = "";

    if( defined ($strainLineArray[ 1 ]) && defined($strainLineArray[ 2 ]) && defined($strainLineArray[ 3 ]) ) {
		$strain = $strainLineArray[ 1 ] . " " . $strainLineArray[ 2 ] . " " . $strainLineArray[ 3 ] . "";
    }

	if( $strain ne "" ) {
    	push( @{ $pvals{ $strain } }, $pval);
    }

    $numberOfLines++;
}

close(INFILE);

print( CYAN, "\n\n\n Number of lines processed in input file: ", RESET );
print( " " . &addCommas($numberOfLines) . "\n" );

my $numberOfStrains += scalar keys %pvals;

print( CYAN, "\n Number of aggregated strains: ", RESET );
print( " " . &addCommas($numberOfStrains) . "\n" );

for my $strain ( sort{ $a cmp $b || $a <=> $b } keys %pvals )
{
    my @pvalArray = @{ $pvals{ $strain } };
    my $numOfPvals = @pvalArray;

    my $sum = 0;

    foreach my $n ( @pvalArray )
    { $sum += $n; }

    my $avgPval = $sum / $numOfPvals;

    my $minPval = List::Util::min( @pvalArray );
    my $maxPval = List::Util::max( @pvalArray );

    print( "\n$strain " );
    print( GREEN, " $numOfPvals ($sum)", RESET );
    print( MAGENTA, " $avgPval", RESET );

    print OUTFILE ( "\n$strain\t$minPval" );
}

close(OUTFILE);

print("\n\n");

# Re-flush STDOUT
$| = 1;



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