#!/usr/bin/python
# coding: utf-8

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
#	The purpose of this program is to create a properly formatted MetaSim profile using GI accessions rather than strain
#   names.  Because of the disparity in naming conventions, and the potential for name space collisions, straight strain
#   names are not useful (and potentially problematic) as their names vary too much — that the names of the strains contain
#   many characters that create lots of problems and edge cases for a regular expression.  It is easier, and more accurate,
#   to simply use GI accession numbers for the MetaSim profiles — using these removes any ambiguity as to what sequence
#   is in the database.
#
#   NOTES:
#   Please see the dependencies section below for the required libraries (if any).
#
#   DEPENDENCIES:
#
#       • None. Uses standard python foundation classes.
#
#   The above libraries & modules (if any) are required. You can check the modules currently installed in your
#   system by running: python -c "help('modules')"
#
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

# 	Python Modules
import os, sys
import argparse
import time
import math
import csv
import re

# ------------------------------------------------------ Main ---------------------------------------------------------

sys.stdout.flush()

print( "[ " + time.strftime('%d-%b-%Y %H:%M:%S',time.localtime()) + " ]" + " Python Starting" + "" )

# 	Pick up the command line arguments
parser = argparse.ArgumentParser()
parser.add_argument( "--mapping-file", dest="mapping_file", required=True, type=str, help="Mapping file" )
parser.add_argument( "--input-profile", dest="input_profile", required=True, type=str, help="Input Profile to clean" )
parser.add_argument( "--output-file", dest="output_file", required=False, type=str, help="Affix for Output files" )
parser.add_argument( "--output-dir", dest="output_dir", required=False, type=str, help="Output Directory" )
args = parser.parse_args()

# 	Variables initialized from the command line arguments
filePathForInputFile    = args.mapping_file
inputProfile            = args.input_profile
output_file             = args.output_file if args.output_file is not None else "default_out.txt"
output_directory        = args.output_dir if args.output_dir is not None else "./out"

# Remove any whitespaces around the file paths
filePathForInputFile.strip()

# -------------------------------------------- Output Files & Directories ---------------------------------------------

# We'll check if the output directory exists — either the default (current) or the requested one
if not os.path.exists( output_directory ):
    print( "[ " + time.strftime('%d-%b-%Y %H:%M:%S',time.localtime()) + " ]" + " Output directory does not exist.  Creating..." + "" )
    os.makedirs( output_directory )

outputFileWithCleanProfile = output_directory + "/" + output_file

outputFileWithNameGIMappingOnly = output_directory + "/mappingFile-name_gi.txt"

# ----------------------------------------------- Mapping File Loading ------------------------------------------------
#
#   The gimmick here is to load the mapping file and create a dictionary of Strain_Name to GI accession number
#

dictionaryOfStrainNamesToGI = {}

print( "[ " + time.strftime('%d-%b-%Y %H:%M:%S',time.localtime()) + " ]" + " Loading Mapping File..." + "" )

numberOfLinesInMappingFile = 0

with open( filePathForInputFile,'r' ) as INFILE:

    reader = csv.reader( INFILE, delimiter='|' )
    writer = csv.writer( open( outputFileWithNameGIMappingOnly, "wb" ), delimiter='\t', lineterminator="\n" )

    try:
        for row_line in reader:

            accession_gi   = str(row_line[ 1 ])
            strain_info    = str(row_line[ 4 ])

            # Clean up the fields we need before we store in the dictionary
            # strain_name  = strain_info.split( "," )[ 0 ]
            # strain_name.lstrip()

            strain_name = strain_info.lstrip()

            dictionaryOfStrainNamesToGI[ strain_name ] = accession_gi

            writer.writerow( [ accession_gi, strain_name ] )

            numberOfLinesInMappingFile += 1

    except csv.Error as e:
        sys.exit( "File %s, line %d: %s" % ( aFile, reader.line_num, e ) )

print( "[ " + time.strftime('%d-%b-%Y %H:%M:%S',time.localtime()) + " ]" + " Number of Lines in mapping file: " + '{:0,.0f}'.format(numberOfLinesInMappingFile) + "" )

print( "[ " + time.strftime('%d-%b-%Y %H:%M:%S',time.localtime()) + " ] " )


# ---------------------------------------------- "Dirty" Profile Loading ----------------------------------------------
#
#   We have a profile from the R script... but the problem is that MetaSim cannot match the strain names.  So we'll load
#   the 'dirty' profiles and change them so that they only contain 'gi' accession references.
#

print( "[ " + time.strftime('%d-%b-%Y %H:%M:%S',time.localtime()) + " ]" + " Loading Profile: " + inputProfile )

numberOfLinesInProfile  = 0
numberOfMatchedStrains  = 0
numberOfStrainsNotFound = 0

with open( inputProfile,'r' ) as INFILE2:

    reader = csv.reader( INFILE2, delimiter='\t' )
    writer = csv.writer( open( outputFileWithCleanProfile, "wb" ), delimiter='\t', lineterminator="\n" )

    try:
        for row_line in reader:

            abundance       = str(row_line[ 0 ])
            key_identifier  = str(row_line[ 1 ])
            strain_info     = str(row_line[ 2 ])

            # Clean the strain name so that we can match it in the GI mapping file
            # strain_name_clean = strain_info.replace( "_", " " )
            strain_name_clean = strain_info
            strain_name_clean = re.sub( " uid[0-9]*", "", strain_name_clean)

            # The profile we are reading in contains the name of the Base-Class Strain.  It does not contain entries for
            # its constituent sequences (chrs, plasmids).  So the job here is to read-in the base strain, extract its
            # abundance value, and then search the Dictionary with GI values to match the strain name.  Once we find a
            # match, we can write out all the matches — ideally these are the constituent sequences.

            print( "\nSearching for: ]" + strain_name_clean + "[")

            for key, value in dictionaryOfStrainNamesToGI.iteritems():  # Iterate on both keys and values

                if key in strain_name_clean:

                    print( "Match: " + strain_name_clean + " -- KEY: " + key + ", VAL: " + value )

                    writer.writerow( [ abundance, "gi", value ] )

                    numberOfMatchedStrains += 1

                else:
                    numberOfStrainsNotFound += 1

            numberOfLinesInProfile += 1

    except csv.Error as e:
        sys.exit( "File %s, line %d: %s" % ( aFile, reader.line_num, e ) )

print( "[ " + time.strftime('%d-%b-%Y %H:%M:%S',time.localtime()) + " ]" + " Number of Lines in profile: " + '{:0,.0f}'.format(numberOfLinesInProfile) + "" )
print( "[ " + time.strftime('%d-%b-%Y %H:%M:%S',time.localtime()) + " ]" + " Number of Matched Strains: " + '{:0,.0f}'.format(numberOfMatchedStrains) + "" )
print( "[ " + time.strftime('%d-%b-%Y %H:%M:%S',time.localtime()) + " ]" + " Number of Not Matched: " + '{:0,.0f}'.format(numberOfStrainsNotFound) + "" )

print( "[ " + time.strftime('%d-%b-%Y %H:%M:%S',time.localtime()) + " ] " )


# ---------------------------------------------------- End of Line ----------------------------------------------------

print( "[ " + time.strftime('%d-%b-%Y %H:%M:%S',time.localtime()) + " ] " )
print( "[ " + time.strftime('%d-%b-%Y %H:%M:%S',time.localtime()) + " ]" + " Done." + "\n" )

sys.exit()
