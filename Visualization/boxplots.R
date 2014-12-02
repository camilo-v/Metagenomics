

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
#	The purpose of this program is to generate BoxPlots for the different abundance & length classes from the
#	confusion matrices.  Each class will have two (2) box plots, one based on the true positive rate (10 values) and
#	the other based on the true negative rate (10 values).
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


# ------------------------------------------------------ Main ---------------------------------------------------------

# 	Set the current working directory -- the Dir in which the results are stored in.
setwd("/Path/To/Simulations/Confusion_Matrix/Matrices_Figures/BoxPlots")

stratType="Abundance"

fileVersion="1"

fileToProcess = paste( "data-", stratType, "/conditionals-v", fileVersion, ".txt", sep="" )

expressionTable = read.csv( fileToProcess, sep="\t", header=TRUE, check.names = FALSE )

#	Melt the data for box_plot format
expressionTableMelted = melt( expressionTable, id.vars = "Group", measure.vars = c(2:11))

# 	Quick print to STDOUT for visual inspection
expressionTableMelted


# ----------------------------------------------------- GGPlot --------------------------------------------------------

xLables = c( "Overall TP", "Overall TN", "Class 1 TPR", "Class 1 TNR", "Class 2 TPR", "Class 2 TNR", "Class 3 TPR", "Class 3 TNR", "Class 4 TPR", "Class 4 TNR" )

boxplot =   ggplot( expressionTableMelted, aes( x=variable, y=value, fill=variable ) ) +
				xlab( paste( "", stratType, " Classes", sep="" ) ) + ylab( "Positive / Negative Rates" ) +
				geom_boxplot( outlier.colour = "red" ) +
				ggtitle( paste( "Accuracy Rates by ", stratType, " Classes" , sep="" ) ) +
				guides( fill=guide_legend(title=NULL) ) +
				theme( text = element_text(size=8), axis.text.x = element_text(angle=50, vjust=0.6) ) +
				scale_fill_manual( name="Rates", labels = c("True Positive", "True Negative"), values = c( "#1bb941", "#f57670", "#1bb941", "#f57670", "#1bb941", "#f57670", "#1bb941", "#f57670", "#1bb941", "#f57670" ) ) +
				scale_x_discrete( labels=xLables )


# 	If you're using TextMate with the 'R' plugin, then you can quickly view the plot
boxplot

outputdirectory = "/Path/To/Simulations/Confusion_Matrix/Matrices_Figures/PDFs_Boxplots/"

ggsave( paste( outputdirectory, "/boxplot_conditional_by_", stratType, "-rates-v", fileVersion, ".pdf", sep="" ) )

