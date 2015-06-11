#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use File::Temp;
use Env qw/TMPDIR/;
use Statistics::R;
use File::Basename;

$SIG{INT} = sub{exit}; #Handle ungraceful exits with CTRL-C.

my $tmpdir = File::Temp->newdir(
    "tmpDirItaggerStatsXXXXXXX",
    DIR => $TMPDIR."/",
    CLEANUP => 0 #Let itagger.pl manage the removal of temp files.
);

sub END{};

my $usage=<<'ENDHERE';
NAME:
itaggerStats.pl

PURPOSE:
Generates stacked barplots in pdf format
from multiple taxonomic summary tables 
typically generated with Qiime.

INPUT:
--infile <string> : OTU table

OUTPUT:
--outdir <string> : one pdf file contaning barplot graph.

NOTES:

BUGS/LIMITATIONS:

AUTHOR/SUPPORT:
Julien Tremblay - jtremblay@lbl.gov
ENDHERE

## OPTIONS
my ($help, $infile, $outdir);
my $verbose = 0;

GetOptions(
	'infile=s' 	=> \$infile,
	'outdir=s'	=> \$outdir,
    'help' 		=> \$help
);
if ($help) { print $usage; exit; }

## VALIDATE
die("--infile please provide one infile\n") unless($infile);
die("--outdir please provide an output directory\n") unless($outdir);

## MAIN
my $R = Statistics::R->new();
$R->startR;

my $hist_string= "";
my @array;
my $row = 2;
my $col = 2;

## R COMMANDS
#$R->run("options(stringsAsFactors = FALSE)");
#$R->run("data <- read.table('".$outfile_table."', sep='\t', header=F, skip=1)");
$R->run("colors = c(
		'#00FF00',
		'#FF8080',
		'#FF00FF',
		'#0000FF',
		'#00CCFF',
		'#CCFFFF',
		'#CCFFCC',
		'#99CCFF',
		'#CC99FF',
		'#FFCC99',
		'#3366FF',
		'#33CCCC',
		'#99CC00',
		'#FF99CC',
		'#FFCC00',
		'#FF9900',
		'#FF6600',
		'#666699',
		'#969696',
		'#003366',
		'#339966',
		'#003300',
		'#333300',
		'#FFFF99',
		'#993300',
		'#993366',
		'#333399',
		'#333333',
		'#000001',
		'#FFFFFF'
	)"
);

my $curr_string = '
	#These libraries are required for the analyses performed below.

	library(gplots)			#used in heat map generation
	library(vegan)			#used in AnoSim
	library(mcpHill)		#for use with multiple adjusted pvalues for alpha diversity 
	library(gclus)			#required for Heatmap of dissimilarity matrix
	library(ade4)			#required for the DISSIMILARITY MATRIX part
	library(cluster) 		#required for the DISSIMILARITY MATRIX part
	library(packfor) 		#required for the RDA analysis
	library(FD) 			#required for hierarchical clustering
	library(RColorBrewer)	#required for hierarchical clustering
	library(labdsv) 		#required for hierarchical clustering
	library(ape)			#required for PCA
	library(BiodiversityR)	#required for PCA
	library(ecodist)		#required for NMDS
	library(ellipse)		#required for NMDS
	library(MASS)			#required for NMDS
	library(indicspecies) 	#for indicator species analysis

	####################	
	#### FUNCTIONS  ####
	####################
	## These functions are necessary for the analyses below.
	
	## function for splitting the consensus lineage term by semicolon and parsing to the desired level of phylogeny.
	extract.name.level = function(x, level){
	  a=c(unlist(strsplit(x,";")),"Other")
	  paste(a[1:min(level,length(a))],collapse=";")
	}

	
	## function takes a transposed otu table with the Consensus lienage removed, a "level" 1-7, and the vector containing the consensus lineage names, usually (taxa.names) 
	otu2taxonomy = function(x, level, taxa=NULL){
	  if(is.null(taxa)){
	    taxa = colnames(x)
	  }
	  #taxa is usually supplied in the function call.
	  
	  if(length(taxa)!=dim(x)[2]){
	    print("ERROR: taxonomy should have the same length
	  as the number of columns in OTU table")
	  return;
	  }
	  #checks dimensions are the same
	  
	  level.names = sapply(as.character(taxa), function(x) extract.name.level(x,level=level))
	  #uses the extract.level.name function
	  
	  t(apply(x, 1, function(y) tapply(y,level.names,sum)))
	  #adds the levels together?
	}

	
	## this function takes a transposed otu table (columns are otus and rows are samples)
	## and a minimum number of read counts and minimum number of samples and removes all OTUs 
	## that do not have at least teh min read counts in at least min number of samples. These
	## values are derived from the threshold script.
	measurableOTUs = function(otu_table, numofcounts, numofsamples){
	  otu_kept <- {} 
	  for (i in 1:dim(otu_table)[2]) {
	    current_otu <- otu_table[,c(i)]
	    if( sum(current_otu>numofcounts) > numofsamples) {
	      otu_kept <- c(otu_kept, colnames(otu_table)[i])
	    }    
	  }
	  otu_table_out <- otu_table[,colnames(otu_table) %in% otu_kept]
	  return(otu_table_out)
	}
	
	## Used in the HEATMAP DISSIMILARITY PORTION OF THE SCRIPT
	hcoplot <- function(tree, diss, k, title=paste("Reordered dendrogram from",deparse(tree$call),sep="\n")) 
	{
	  require(gclus)
	  gr <- cutree(tree, k=k)
	  tor <- reorder.hclust(tree, diss)
	  plot(tor, hang=-1, xlab=paste(length(gr),"sites"), sub=paste(k,"groups"), main=title)
	  so <- gr[tor$order]
	  gro <- numeric(k)
	  for (i in 1:k) {
	    gro[i] <- so[1]
	    if (i<k) so = so[so!=gro[i]]
	  }
	  rect.hclust(tor, k=k, border=gro+1, cluster=gr)
	  legend("topright", paste("Group",1:k), pch=22, col=2:(k+1), bty="n")
	}
	
	## Used in the dissimilarity measure part of the script.
	coldiss <- function(D, nc = 4, byrank = TRUE, diag = FALSE)
	{
	  require(gclus)
	  
	  if (max(D)>1) D <- D/max(D)
	  
	  if (byrank) {
	    spe.color = dmat.color(1-D, cm.colors(nc))
	  }
	  else {
	    spe.color = dmat.color(1-D, byrank=FALSE, cm.colors(nc))
	  }
	  
	  spe.o = order.single(1-D)
	  speo.color = spe.color[spe.o,spe.o]
	  
	  op = par(mfrow=c(1,2), pty="s",cex=0.4,cex.main=2,cex.lab=2)
	  
	  if (diag) {
	    plotcolors(spe.color, rlabels=attributes(D)$Labels, 
	               main="Dissimilarity Matrix", 
	               dlabels=attributes(D)$Labels)
	    plotcolors(speo.color, rlabels=attributes(D)$Labels[spe.o], 
	               main="Ordered Dissimilarity Matrix", 
	               dlabels=attributes(D)$Labels[spe.o])
	  }
	  else {
	    plotcolors(spe.color, rlabels=attributes(D)$Labels, 
	               main="Dissimilarity Matrix")
	    plotcolors(speo.color, rlabels=attributes(D)$Labels[spe.o], 
	               main="Ordered Dissimilarity Matrix")
	  }
	  
	  par(op)
	}
	
	## Used in the dissimilarity measure part of the script.
	coldiss2 <- function(D, nc = 4, byrank = TRUE, diag = FALSE)
	{
	  require(gclus)
	  
	  if (max(D)>1) D <- D/max(D)
	  
	  if (byrank) {
	    spe.color = dmat.color(1-D, cm.colors(nc))
	  }
	  else {
	    spe.color = dmat.color(1-D, byrank=FALSE, cm.colors(nc))
	  }
	  
	  spe.o = order.single(1-D)
	  speo.color = spe.color[spe.o,spe.o]
	  
	  op = par(mfrow=c(1,1), pty="s",cex=0.6,cex.main=2,cex.lab=2)
	  
	  if (diag) {
	  #  plotcolors(spe.color, rlabels=attributes(D)$Labels, 
	   #            main="Dissimilarity Matrix", 
	    #           dlabels=attributes(D)$Labels)
	    plotcolors(speo.color, rlabels=attributes(D)$Labels[spe.o], 
	               main="Ordered Dissimilarity Matrix", 
	               dlabels=attributes(D)$Labels[spe.o])
	  }
	  else {
	    #plotcolors(spe.color, rlabels=attributes(D)$Labels, 
	     #          main="Dissimilarity Matrix")
	    plotcolors(speo.color, rlabels=attributes(D)$Labels[spe.o], 
	               main="Ordered Dissimilarity Matrix")
	  }
	  
	  par(op)
	}
	
	## this is necessary to specify the evplot function
	## this is another function required by the RDA analysis
	evplot = function(ev) {
	  # Broken stick model (MacArthur 1957)
	  n = length(ev)
	  bsm = data.frame(j=seq(1:n), p=0)
	  bsm$p[1] = 1/n
	  for (i in 2:n) {
	    bsm$p[i] = bsm$p[i-1] + (1/(n + 1 - i))
	  }
	  bsm$p = 100*bsm$p/n
	  # Plot eigenvalues and % of variation for each axis
	  op = par(mfrow=c(2,1))
	  barplot(ev, main="Eigenvalues", col="bisque", las=2)
	  abline(h=mean(ev), col="red")
	  legend("topright", "Average eigenvalue", lwd=1, col=2, bty="n")
	  barplot(t(cbind(100*ev/sum(ev), bsm$p[n:1])), beside=T, 
	          main="% variation", col=c("bisque",2), las=2)
	  legend("topright", c("% eigenvalue", "Broken stick model"), 
	         pch=15, col=c("bisque",2), bty="n")
	  par(op)
	}

	######################
	####### MAIN #########
	######################

	setwd("'.$outdir.'")
	# Wherever you want files opened from and written to.
	
	
	## IMPORT OTU TABLE
	
	## This portion of the script opens the OTU table as a table and creates several permutations 
	## of this table: 
	## com                             the raw table, with rownames and header: 91733x142
	## raw_com                         the backup file: 91733 x 142
	## com_noRDP                       removes the consensus lineage column: 91733 x 141
	## t_com_noRDP                     the transposed table, rows are samples, columns are 
	##                                 OTUs 141 x 91733
	## t_com_noRDP_1000                interim file before rarifying that has only samples
	##                                 with >1000 reads. 70 x 91733
	## t_com_noRDP_rare                all samples rarefied to 1000 reads 70 x 19688
	## t_com_noRDP_rare_log            log 10 transformed rarefied file. 70 x 19688
	## t_com_noRDP_measurable          as t_com_noRDP, but with OTUs removed that do not 
	##                                 have at least x min reads,in at least y samples,
	##                                 +samples with less than 1000 reads removed. 70 x 2106 
	## t_com_noRDP_measurable_log      as t_com_noRDP_measureable, but log10 transformed: 70 x 2106  
	## t_com_noRDP_measurable_rare     as t_com_noRDP_measurable, but with samples rarefied
	##                                 to 1000 reads per sample 69 x 2106
	## t_com_noRDP_measurable_rare_log as t_com_noRDP_measurable_rare, but with log10 transformed 
	##                                 readcounts 69 x 2106
	
	## Using sep = "\t" removes the need for eliminating spaces in the colnames in 
	## excel before importing to R.
	com <- read.table("otu_table_full.txt", header=TRUE,row.names=1,sep="\t")
	
	## Use the last line to remove Epibuffer and CleanBuffer.
	#com<-com[,-c(134,135)]
	
	## This removes the extra information in chloro and mito annotations : k__chloro_tritcum > k__chloro
	com$Consensus.lineage <- gsub("((\\\\w)__(\\\\w*))_(\\\\w*);", "\\\\1;", com$Consensus.lineage, perl=T)
	
	## Saves the raw OTU table to a new variable for later use
	raw_com <- com
	
	## Saves the OTU ID and RDP lineage for future use in a table called taxa_names.
	taxa_names <- as.data.frame(com[,dim(com)[2]])
	OTU_ID<-rownames(com)
	rownames(taxa_names)<- OTU_ID
	colnames(taxa_names) <- c("Consensus.lineage")
	taxa_names<-cbind(taxa_names,OTU_ID)
	taxa.names = taxa_names$Consensus.lineage
	
	## Removes consensus lineage column
	com_noRDP = as.matrix(com[,-dim(com)[2]])
	
	## *** OPTIONAL *** #
	## com_noRDP_RA = scale(com_noRDP, center=F, scale=colSums(com_noRDP))
	## Can be used to get relative abundance instead of raw read counts; *** would need to change input files below.
	
	## This transposes the com table, so that rows are samples and columns are OTUs.
	t_com_noRDP <- t(com_noRDP)
	
	## Converts to a dataframe.
	t_com_noRDP<-as.data.frame(t_com_noRDP)
	
	## Gives the sorted readcounts for all samples. This can be used to determine a reasonable
	## threshold for the number of reads in the rarefaction below. A value between 1000 and 10,000
	## seems reasonable.
	sort(rowSums(t_com_noRDP))
	
	## Removes all samples with less than the specified number of reads, which is necessary 
	## for the rarefaction that takes place next.
	t_com_noRDP_1000 <- t_com_noRDP[rowSums(t_com_noRDP)>1000,]
	
	## Produces the rarefied OTU table with all samples rarefied to the speceified number of reads
	## NOTE: Can rarefy to a value just below the minimum read count in remaining samples.
	## NOTE: Takes a long time with large files. (with 70 samples and 91000 OTUs ~20 minutes)
	t_com_noRDP_rare <- rrarefy(t_com_noRDP_1000,2500)
	
	## an alternative for reprocessing the same files later is to read in the rare file:
	write.table (t_com_noRDP_rare, file = "t_com_noRDP_rare.txt", col.names = T, row.names = T, 
	             quote=F,sep = "\t")
	t_com_noRDP_rare<- read.table("t_com_noRDP_rare.txt", header=TRUE, row.names=1,sep="\t")
	
	## NOTE***May need to execute the following line if the colnames of t_com_noRDP_rare 
	## have an "X" preceding the otu ids.
	colnames(t_com_noRDP_rare) <- gsub("X(\\\\d*)","\\\\1",colnames(t_com_noRDP_rare),perl=T)
	
	## Removes any OTU that does not have at least 1 read counts. Many OTUs
	## will have 0 reads after rarifying.
	t_com_noRDP_rare <- t_com_noRDP_rare[,colSums(t_com_noRDP_rare)>0]
	
	## Converts from double/numeric matrix to dataframe, necessary for downstream analysis.
	t_com_noRDP_rare <- as.data.frame(t_com_noRDP_rare)
	
	## gets the log transform, necessary for some parts of the script.
	t_com_noRDP_rare_log <- log(t_com_noRDP_rare+1)
	
	## NOTE:For creating a measureable OTU table, the appropriate threshold values should be 
	## determined from an an analysis of replicate samples with the OTU_Threshold.R script.
	
	## NOTE* Define minreadcount and minsamplenumber before proceeding!
	minreadcount=7 # DEFINE!
	minsamplenumber=3 # DEFINE!
	
	## Executes the function that removes all OTUs without at least minreacount number of reads 
	## in at least minsamplenumber of samples.
	t_com_noRDP_measureable <- measurableOTUs(t_com_noRDP,minreadcount,minsamplenumber)
	
	## Gives the sorted readcounts for remaining samples. Set the threshold for rarefying to 
	## just above the minimum value you feel is reasonable. 
	sort(rowSums(t_com_noRDP_measureable))
	
	## Removes all samples with less than the specified number of reads, which is necessary 
	## For the rarefaction that takes place next.
	t_com_noRDP_measureable <- t_com_noRDP_measureable[rowSums(t_com_noRDP_measureable)>5000,]
	
	## Log 10 transform of the measurable OTU file.
	t_com_noRDP_measureable_log <- log(t_com_noRDP_measureable+1)
	
	## Produces the rarefied OTU table with all samples rarefied to the specified number of reads
	## NOTE: Can rarefy to a value just below the min read count in remaining samples.
	## NOTE:takes a long time with large files.
	t_com_noRDP_measureable_rare <- rrarefy(t_com_noRDP_measureable,9000)
	
	## Removes any OTU that doe not have at least 1 read counts. Usuall at this stage no OTUs
	t_com_noRDP_measureable_rare <- t_com_noRDP_measureable_rare[,colSums(t_com_noRDP_measureable_rare)>0]
	## will have 0 reads after rarifying.
	
	## Converts to a data frame, necessary for downstream analysis.
	t_com_noRDP_measureable_rare <- as.data.frame(t_com_noRDP_measureable_rare)
	
	## The associated log file.
	t_com_noRDP_measureable_rare_log <- log(t_com_noRDP_measureable_rare+1)
';

push(@array, $curr_string);

my $outfile_graph = $outdir."/graph.pdf";

$R->run('pdf("'.$outfile_graph.'" , width=7, height=7)');
$R->run('par(mfrow=c('.$row.','.$col.'))'); 
foreach my $cmd (@array){
	$R->run($cmd);
	#print $cmd."\n";
}
$R->run('dev.off()');

system "rm -rf ".$tmpdir; 
1;
exit;
