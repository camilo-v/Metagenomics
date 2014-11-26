
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
#	The purpose of this script is to compute pvalues for the references in the database using the base (0th) counts and
#	the counts at the set value of "q".
#
#   NOTES:
#   Please see the dependencies section below for the required libraries (if any).
#
#   DEPENDENCIES:
#
#       • multtest.
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

rm(list=ls())
setwd("/Path/To/metagenomic_simulations/groups/group_GRP/analysis/pvalues/TEMPLATE")

# set mutation rate, number of mutated data sets desired, and number of sample reads
qval<-TEMPLATE
nperm<-200
totsampcnt<-150000*2
rlength<-100

# read in number of reads that align when no mutations are done
# we don't actually need to do this nperm times, since the results should not vary if
# no mutations are being done

print( paste( "[", format(Sys.time(), "%m/%d/%y %H:%M:%S"),"] ", sep="") )
print( paste( "[", format(Sys.time(), "%m/%d/%y %H:%M:%S"),"] Starting... ", sep="") )

######################################################################################################
#
#	Base "0" Counts
#	this does not change with the choice of q, so it is possible to compute once, save
#	in a workspace, and then simply load
#

datall<-read.csv('/Path/To/metagenomic_simulations/groups/group_GRP/analysis/counts/0/0_counts_clean.txt', sep="\t", header=TRUE, row.names=1, allowEscapes=T, colClasses=c('character',rep('numeric',nperm+1)))

print( paste( "[", format(Sys.time(), "%m/%d/%y %H:%M:%S"),"] Null Counts (0th Counts)", sep="") )

print( "/Path/To/metagenomic_simulations/groups/group_0/analysis/counts/0/0_counts_clean.txt" )

data<-matrix(199,dim(datall)[1],2)
data[,1]<-datall[,1] # get genome lengths

data[,1]

# take mean of null runs (different run results due to variability of aligner)
data[,2]<-apply(datall[,2:dim(datall)[2]],1,mean)
totrefcnt<-sum(data[,1])
nref<-dim(data)[1]
rm(datall)
save.image('0count.RData')

# adjust hypergeometric based p-values by BH, resort to match original list
library(multtest)

######################################################################################################
#
#	Q Counts
#	(Q) is some value from 2..30
#	Read in permutation based counts, assuming nperms for each reference genome
#

filename<-paste('/Path/To/metagenomic_simulations/groups/group_GRP/analysis/counts/',qval,'/',qval,'_counts_clean.txt', sep="")

datap<-read.csv(file=filename, sep="\t", header=TRUE, row.names=1, allowEscapes=T, colClasses=c('character',rep('numeric',nperm+1)))

print( paste( "[", format(Sys.time(), "%m/%d/%y %H:%M:%S"),"] Q Counts", sep="") )
print( filename )


######################################################################################################
#
#	Apply single-step WY correction as described in Meinshausen et al. (2011)
#	Remember by default sorting is done ascending
#

# convert raw permutation count data to adjusted count data using read and genome lengths

pperm<-apply(datap[,2:(nperm+1)],2,function(x) (x*rlength)/data[,1])

quant<-seq(0.90,1,0.01)
tmp<-apply(pperm,2,function(x) quantile(x,quant))
plot(quant,apply(tmp,1,mean))

# Adjust permutation based p-values by single-step WY, note rank matrix is nperm*nref
# note if even one reference file is 'weird' this goes very badly wrong
# maxpp<-apply(pperm,2,max)

valueForQuantile = 0.95
maxpp<-apply(pperm,2,function(x) quantile(x,valueForQuantile))  # The Non-Permuted Has to beat the MAX!

tmp1<-cbind(t(matrix(rep(maxpp,nref),nperm,nref)),data[,2]*rlength/data[,1])

rank<-apply(tmp1,1,function(x) order(x, decreasing=T))

ppvalWF<-apply(rank,2,function(x) (which(x==(nperm+1))-1)/(nperm+1))

outfile2<-paste('perm',qval,'.permWY1-', valueForQuantile, '.txt',sep="")

#
#	NO COMMA SEPARATED OUTPUT FILES, THE STRAIN NAMES HAVE COMMAS AND ANNIHILATE ANY PARSING EFFORTS.
#

#write.csv(cbind(ppvaladj$rawp,ppvaladj$adjp,ppvalWF,datap), file=outfile2)
#write.csv(cbind(ppvalWF,datap), file=outfile2)

write.table(cbind(ppvalWF,datap), file=outfile2, sep="\t")
