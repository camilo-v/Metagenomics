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
#	The purpose of this program is to select a random number of bacterial genomes and split them into a given number
#   of groups.  The original version of this script randomly selected the groups at the Genus level; the current
#   version selects the groups at the Strain level.  An additional parameter was added to define the starting number
#   of the "Conditional Population", that is, not the overall Population taken from the total number of Strains, but
#   the subset of the Population to be used as a "Conditional Population" — a subset that for all intents and
#   purposes is treated as the starting Population.
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

# 	Python Modules
import os, sys
import collections
import random
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
parser.add_argument( "--strains", required=True, type=str, help="Input File with Strains" )
parser.add_argument( "--groups",  required=True, type=int, help="Number of groups to split by" )
parser.add_argument( "--population", required=True, type=int, help="Conditional Population to start with" )
parser.add_argument( "--affix", required=False, type=str, help="Affix for Output files" )
parser.add_argument( "--out", required=False, type=str, help="Output Directory" )
args = parser.parse_args()

# 	Variables initialized from the command line arguments
filePathForInputFile    = args.strains
numberOfGroups          = args.groups
coditionalPopulation    = args.population
affixForOutputFile      = args.affix if args.affix is not None else "default_out"
output_directory        = args.out if args.out is not None else "./out"

# Remove any whitespaces around the file paths
filePathForInputFile.strip()

# -------------------------------------------- Output Files & Directories ---------------------------------------------

# We'll check if the output directory exists — either the default (current) or the requested one
if not os.path.exists( output_directory ):
    print( "[ " + time.strftime('%d-%b-%Y %H:%M:%S',time.localtime()) + " ]" + " Output directory does not exist.  Creating..." + "" )
    os.makedirs( output_directory )

# --------------------------------------------------- File Loading ----------------------------------------------------

print( "[ " + time.strftime('%d-%b-%Y %H:%M:%S',time.localtime()) + " ]" + " Loading Files..." + "" )

arrayOfAllStrains = []

numberOfLinesInFile = 0

with open( filePathForInputFile,'r' ) as INFILE:

    reader = csv.reader( INFILE, delimiter='|' )

    try:
        for row_line in reader:

            strain_info = str(row_line[ 4 ])

            # Only use 'complete genomes'
            regexp  = re.compile(r'complete genome')

            # And by 'complete genomes' we mean NO contigs
            regexp2 = re.compile(r'contig')

            if( ( regexp.search(strain_info) ) and (not regexp2.search(strain_info) ) ):

                strain = strain_info.strip()

                arrayOfAllStrains.append( strain )

            numberOfLinesInFile += 1

    except csv.Error as e:
        sys.exit( "File %s, line %d: %s" % ( aFile, reader.line_num, e ) )

print( "[ " + time.strftime('%d-%b-%Y %H:%M:%S',time.localtime()) + " ]" + " Number of Lines in file: " + '{:0,.0f}'.format(numberOfLinesInFile) + "" )

print( "[ " + time.strftime('%d-%b-%Y %H:%M:%S',time.localtime()) + " ] " )


# ------------------------------------------- Conditional Population ------------------------------------------

print( "[ " + time.strftime('%d-%b-%Y %H:%M:%S',time.localtime()) + " ] Number of Strains: " + str( len( arrayOfAllStrains ) ) )
print( "[ " + time.strftime('%d-%b-%Y %H:%M:%S',time.localtime()) + " ] Conditional Population: " + str( coditionalPopulation ) )
print( "[ " + time.strftime('%d-%b-%Y %H:%M:%S',time.localtime()) + " ] Number of Groups: " + str( numberOfGroups ) )
print( "[ " + time.strftime('%d-%b-%Y %H:%M:%S',time.localtime()) + " ] " )

#   Shuffle the Strain list so that no input bias is present
random.shuffle( arrayOfAllStrains )

#   Randomly Select a "N" number of Strains from the overall list.  "N" value is given by the "population" CLI argument
arrayOfRandomConditionalPopulation = []

print( "[ " + time.strftime('%d-%b-%Y %H:%M:%S',time.localtime()) + " ] Selecting Conditional Population..." )

# Output file for the Conditional Population
outputFileForConditionalPopulation = output_directory + "/conditional_population.txt"

writerForConditionalPopulation = csv.writer( open( outputFileForConditionalPopulation, "wb" ), delimiter='\t' )

for populationNumber in range( 0, coditionalPopulation ):

    randomStrain = random.choice( arrayOfAllStrains )
    arrayOfRandomConditionalPopulation.append( randomStrain )
    arrayOfAllStrains.remove( randomStrain )
    random.shuffle( arrayOfAllStrains )

    writerForConditionalPopulation.writerow( [ randomStrain ] )


print( "[ " + time.strftime('%d-%b-%Y %H:%M:%S',time.localtime()) + " ] Conditional Population Size: " + str( len( arrayOfRandomConditionalPopulation ) ) )
print( "[ " + time.strftime('%d-%b-%Y %H:%M:%S',time.localtime()) + " ] No. of Strains after selection: " + str( len( arrayOfAllStrains ) ) )
print( "[ " + time.strftime('%d-%b-%Y %H:%M:%S',time.localtime()) + " ] " )

# ------------------------------------------ Random Group Selection -------------------------------------------

strainsPerGroup = int( math.floor( len( arrayOfRandomConditionalPopulation ) / numberOfGroups ) )

print( "[ " + time.strftime('%d-%b-%Y %H:%M:%S',time.localtime()) + " ] Creating Groups..." )
print( "[ " + time.strftime('%d-%b-%Y %H:%M:%S',time.localtime()) + " ] Strains per Group: " + str( strainsPerGroup ) )


for groupNumber in range( 0, numberOfGroups ):

    print( "[ " + time.strftime('%d-%b-%Y %H:%M:%S',time.localtime()) + " ] Creating Group: " + str(groupNumber) )

    random.shuffle( arrayOfRandomConditionalPopulation )    # Re-shuffle the order before starting a new group

    # Groups are writen to output files
    outputFile = output_directory + "/" + affixForOutputFile + "-Group_" + str(groupNumber) + ".txt"

    writer = csv.writer( open( outputFile, "wb" ), delimiter='\t' )

    for element in range( 0, strainsPerGroup ):

        randomStrainFromConditionalPopulation = random.choice( arrayOfRandomConditionalPopulation )
        arrayOfRandomConditionalPopulation.remove( randomStrainFromConditionalPopulation )

        writer.writerow( [ randomStrainFromConditionalPopulation ] )


print( "[ " + time.strftime('%d-%b-%Y %H:%M:%S',time.localtime()) + " ] " )

# ------------------------------------------------ End of Line ------------------------------------------------

print( "[ " + time.strftime('%d-%b-%Y %H:%M:%S',time.localtime()) + " ]" + " Done." + "\n" )

sys.exit()
