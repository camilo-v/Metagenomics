
# ---------------------------------------------------------------------------------------------------------------------
#
#                                   		Center for Computational Science
#												http://www.ccs.miami.edu/
#                             			  			University of Miami
#
#
#   This software is a "University of Miami Work" under the terms of the United States Copyright Act.
#   Please cite the author(s) in any work or product based on this material.
#
#   OBJECTIVE:
#	The purpose of this program is to drive calculate the pvalues for the IMG references using the HMP data.
#
#   NOTES:
#   Please see the dependencies section below for the required libraries (if any).
#
#   DEPENDENCIES:
#
#       • multtest
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

rm(list=ls())
setwd("/Path/To/metagenomics/permutations/mid_vaginal/analysis/p_vals/")

# set mutation rate, number of mutated data sets desired, and number of sample reads
qval<-2
nperm<-350
totsampcnt<-495256*2
rlength<-100

# ------------------------------------------------- Base Counts -------------------------------------------------------
#
#	Base "0" Counts
# 	We read in the number of reads that align when no mutations are done.  Note that we do not actually need to do
#	this "nperm" times, since the results should not vary if when no mutations are being done.  The DNA aligner should
#	be deterministic — luckily for us, Bowtie2 is.
#	The Base counts for each group do not change with the choice of q, so it is possible to compute once, save in a
#	workspace, and then load it back into memory.
#

datall<-read.csv('/Path/To/metagenomics/permutations/mid_vaginal/analysis/counts/0/0_counts.txt', sep="\t", header=TRUE, row.names=1, allowEscapes=T, colClasses=c('character',rep('numeric',nperm+1)))

data<-matrix(349,dim(datall)[1],2)
data[,1]<-datall[,1] # get genome lengths

# take mean of null runs (different run results due to variability of aligner)
data[,2]<-apply(datall[,2:dim(datall)[2]],1,mean)
totrefcnt<-sum(data[,1])
nref<-dim(data)[1]
rm(datall)
save.image('0count.RData')

# adjust hypergeometric based p-values by BH, resort to match original list
library(multtest)

# ------------------------------------------------ Mutated Counts ------------------------------------------------------
#
#	Q Counts
#	(Q) is some value from 2..30
#	Read in permutation based counts, assuming nperms for each reference genome
#

filename<-paste('/Path/To/permutations/mid_vaginal/analysis/counts/',qval,'/',qval,'_counts.txt',sep="")

datap<-read.csv(file=filename, sep="\t", header=TRUE, row.names=1, allowEscapes=T, colClasses=c('character',rep('numeric',nperm+1)))


# ----------------------------------------- Multiple Testing Correction -----------------------------------------------
#
#	Apply single-step WY correction as described in Meinshausen et al. (2011)
#	Remember by default sorting is done ascending
#

# convert raw permutation count data to adjusted count data using read and genome lengths
pperm<-apply(datap[,2:(nperm+1)],2,function(x) (x*rlength)/data[,1])

# Adjust permutation based p-values by single-step WY, note rank matrix is nperm*nref
# note if even one reference file is 'weird' this goes very badly wrong
maxpp<-apply(pperm,2,max)

tmp1<-cbind(t(matrix(rep(maxpp,nref),nperm,nref)),data[,2]*rlength/data[,1])

rank<-apply(tmp1,1,function(x) order(x, decreasing=T))

ppvalWF<-apply(rank,2,function(x) (which(x==(nperm+1))-1)/(nperm+1))

outfile2<-paste('perm',qval,'.permWY1.csv',sep="")


# ----------------------------------------------------- Warning -------------------------------------------------------
#
#	No comma-separated output files.  Tab-delimited only.  The strains names have commas and cause all kinds of
#	problems.
#
write.table(cbind(ppvalWF,datap), file=outfile2, sep="\t")
