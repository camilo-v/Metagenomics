
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
#	The purpose of this program is to assign a set of abundance values to a list of bacterial strains.  The list of
#	strains is a subgroup which was randomly selected from a larger list of bacterial strains known as the
#	"Conditional Population".  The size of the conditional population is 60, and each subgroup has size 60 and should
#	not contain any strains present in the other groups.
#
#	USAGE:
#	The input file 'fileToProcess' is a text-file with the name of the strains.  It has one row per strain with one column,
#	and the biological classification (Genus, Species, etc.) is separated with spaces.
#
#   NOTES:
#   Please see the dependencies section below for the required libraries (if any).
#
#   DEPENDENCIES:
#
#       • reshape2
#       • ggplot2
#       • plyr
#       • untb
#       • vegan
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

# 	Import any necessary libraries here

library(reshape2)
library(ggplot2)
library(plyr)
library(untb)
library(vegan)

options( digits = 15 ) # Enable printing of large numbers. R defaults to "rounding" for display.

# ---------------------------------------------------- Functions ------------------------------------------------------

#
# Function to implement a 'cbind.fill' procedure
# @params	List of data.frames to cbind together.  They could contain different lengths (no. of rows)
# @returns 	A new data.frams with combined columns
# @author	http://stackoverflow.com/questions/7962267/cbind-a-df-with-an-empty-df-cbind-fill
#
cbind.fill <- function(...){
    nm <- list(...)
    nm <- lapply(nm, as.matrix)
    n <- max(sapply(nm, nrow))
    do.call(cbind, lapply(nm, function (x)
        rbind(x, matrix(, n-nrow(x), ncol(x)))))
}

# ------------------------------------------------------ Main ---------------------------------------------------------

# 	Set the current working directory
setwd("/Path/To/Working/Directory/")

print( paste( "[", format(Sys.time(), "%m/%d/%y %H:%M:%S"),"] R Starting... ", sep="") )

# We previously created 10 Groups, and these are contained in 10 files

pathToInputFilesForGroups = "/Path/To/Output/Directory/From/Random/Selections/Random_Picks"

# i.e., Number of Individuals — total mass of Genome Copy Numbers
numberOfSequencingReads = 150000

# the rare species advantage term
c_rare_species = 5

for( groupNumber in 0:9 ) {

	fileToProcess = paste( pathToInputFilesForGroups, "/random_strains_picks-Group_", groupNumber, ".txt", sep="" )

	print(fileToProcess)

	strainTable = read.delim( fileToProcess, header=FALSE )

	numberOfStrains = nrow( strainTable )

	print( numberOfStrains )

	alpha = fishers.alpha( numberOfSequencingReads, numberOfStrains, give=FALSE )

	print( paste( "[", format(Sys.time(), "%m/%d/%y %H:%M:%S"),"] Alpha: ", alpha, sep="") )

	abundancesForStrains = fisher.ecosystem( numberOfSequencingReads, numberOfStrains, numberOfSequencingReads )

	# What type of Object does 'fisher.ecosystem' return
	print(class(abundancesForStrains))

	# Cast the 'fisher.ecosystem' data structure to a more manageable DataFrame
	dataFrameForAbundances = as.data.frame(abundancesForStrains)

	abundanceDataFrame = cbind.fill( strainTable, dataFrameForAbundances )

	colnames(abundanceDataFrame) = c("Strain", "Strain Number","Abundance")

	# print( abundanceDataFrame )

	# Basic output file for record-keeping.  We'll generate another output file that is formatted specifically for MetaSim
	outputFile = paste( pathToInputFilesForGroups, "/Abundances/abundance_strains-Group_", groupNumber, ".txt", sep="" )
	write.table( abundanceDataFrame, file=outputFile, append=FALSE, row.names=FALSE, sep="\t" )

	#
	# MetaSim Profile output file
	# format for profile is <abundance_value> <key_identifier> <key_value>
	#

	metaSimProfile = data.frame( matrix(0, ncol = 1, nrow = 60) )

	for( i in 1:60 ) {
	  metaSimProfile[ i, ] = "name"
	 }

	colnames(metaSimProfile) = c("key_identifier")

	metaSimProfileAsDF = cbind(abundanceDataFrame[,3], metaSimProfile[,1], abundanceDataFrame[,1] )

	# print( metaSimProfileAsDF )

	outputFileMetaSim = paste( pathToInputFilesForGroups, "/Profiles/profileForMetaSim-Group_", groupNumber, ".txt", sep="" )

	# There are some "N/A" present from the previous "cbind.fill" step.  We'll replace them with Zeros so that we have
	# a properly constructed MetaSim profile.
	write.table( metaSimProfileAsDF, file=outputFileMetaSim, append=FALSE, row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t", na="0" )
}

print( paste( "[", format(Sys.time(), "%m/%d/%y %H:%M:%S"),"] FIN.", sep="") )
