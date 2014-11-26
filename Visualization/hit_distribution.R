
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
#	The purpose of this program is to create a Box Plot with GGPLot2 to display the average number of hits that a
#	read has with increasing values of the mutation rate Q.
#
#   NOTES:
#   Please see the dependencies section below for the required libraries (if any).
#
#   DEPENDENCIES:
#
#       • reshape2
#       • ggplot2
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

options( digits = 15 ) # Enable printing of large numbers. R defaults to "rounding" for display.


# ----------------------------------------------------- Setup ---------------------------------------------------------

print( paste( "[", format(Sys.time(), "%m/%d/%y %H:%M:%S"),"] R Starting... ", sep="") )

working_directory="/Path/To/counts/read_hit_distribution/group_NUMBER"

ggplot_output_dir="pdf"

setwd( file.path( working_directory ) )

# ------------------------------------------------------ Main ---------------------------------------------------------

groupNumber = "0"

# Y-axis display limits
yAxis_lower = -0
yAxis_upper = 2

averageCounts = read.csv( "merged_hits.txt", sep="\t", header=TRUE )

averageCounts_melted <- melt( averageCounts )

ggplot(averageCounts_melted, aes(x=variable, y=value, fill=variable )) + geom_boxplot() + xlab("Mutation Rate Q") +
				ylab("Average Number of Reference Mappings") + ylim( yAxis_lower, yAxis_upper ) +
				ggtitle( paste( "Group ", groupNumber, sep="" ) )


# Create the output directory for GGPlot if it does not exist.
dir.create(file.path(working_directory, ggplot_output_dir), showWarnings = FALSE)

ggplot_output_file= paste( ggplot_output_dir, "/group_", groupNumber, "-average_read_mappings-fig-1_B.pdf", sep="" )

ggsave( ggplot_output_file )


print( paste( "[", format(Sys.time(), "%m/%d/%y %H:%M:%S"),"] Done.", sep="") )

