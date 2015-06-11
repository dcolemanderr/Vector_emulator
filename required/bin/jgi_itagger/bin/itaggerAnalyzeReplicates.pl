#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use File::Temp;
use Data::Dumper;
use Env qw/TMPDIR/;
use Statistics::R;
use File::Basename;
use iTagger::FastaDb;
use Statistics::Descriptive;

$SIG{INT} = sub{exit}; #Handle ungraceful exits with CTRL-C.

my $usage=<<'ENDHERE';
NAME:
itaggerAnalyzeReplicates.pl

PURPOSE:
Compare replicates of OTUs (at the moment, only replicates 1 and 2) and
plot the OTU abundance of rep 1 and 2. Determine R2 correlation and do
a progressive drop-out analysis. Return value at which OTUs shows a R2
square of 0.87. To do that you MUST have these three fields in your
mapping file: replicate, barcode, rep_common_name.

INPUT:
--infile_otu_table <string>  : File having cluster (OTUs) representative.
                               Headers must be properly formatted and have
                               barcodes information.
--infile_mapping <string>    : Mapping file having a field called barcodes
--threshold                  : Abundance threshold. Default = 0.87
--add_duplicates <string>    : Optional. will write sums of duplicates to a
                               separate file.

OUTPUT:
--outfile_graph <string>     : one pdf file contaning plots
--outfile_freq <string>      : outfile for OTU frequency table of all OTUs
--outfile_otu_table <string> : OTU table on which --threshold has been applied.
--outdir <string>            : outdir to write OTU frquency tables for each barcodes.

NOTES:

BUGS/LIMITATIONS:

AUTHOR/SUPPORT:
Julien Tremblay - jtremblay@lbl.gov
ENDHERE

## OPTIONS
my ($help, $infile_otu_table, $infile_mapping, $outfile_graph, $outfile_freq, $outdir, $threshold, $outfile_otu_table, $add_duplicates);
my $verbose = 0;

GetOptions(
	'infile_otu_table=s'	=> \$infile_otu_table,
	'outfile_otu_table=s'	=> \$outfile_otu_table,
	'infile_mapping=s'		=> \$infile_mapping,
	'threshold=f'			=> \$threshold,
	'outfile_graph=s'		=> \$outfile_graph,
	'outfile_freq=s'		=> \$outfile_freq,
	'outdir=s'				=> \$outdir,
	'add_duplicates'		=> \$add_duplicates,
    'help' 					=> \$help
);
if ($help) { print $usage; exit; }

## VALIDATE
die("--infile_otu_table please provide at least one infile\n") 			unless($infile_otu_table);
die("--outfile_otu_table please provide at least one outfile\n") 		unless($outfile_otu_table);
die("--infile_mapping please provide at least one infile\n") 			unless($infile_mapping);
die("--outfile please provide an outfile for tabular values\n") 		unless($outfile_graph);
die("--outfile_freq please provide an outfile for tabular values\n") 	unless($outfile_freq);
die("--outdir please provide an outdir for tabular values\n") 			unless($outdir);

$threshold = 0.87 unless($threshold);
my $threshold_OTU;
my $r_squared_table = $outdir."/r-squared_table.txt";
## MAIN
no warnings 'numeric';
no warnings 'redefine';

my $tmpdir = File::Temp->newdir(
    "tmpDirItaggerAnalyzeReplicatesXXXXXX",
    DIR => $TMPDIR."/",
    CLEANUP => 0 #Let itagger.pl manage the removal of temp files.
);

my $legend_cex = 0.2;
my $row = 3;
my $col = 3;
my $cex = 0.7;
my $format = "c(1,10,20,30,40)";

my $hist_string= "";
my @array;
my $id;

my $tmp_infile = $tmpdir."/tmpinfile.txt";
my $number_of_samples = 0;

my %hash_by_barcodes;
my %hash_all;

my $pos_replicate;
my $pos_barcode;
my $pos_common_id;
my $pos_sample_id;

# Parse mapping file. Look for fields "replicate", "barcode" , "SampleID" and "rep_common_name".

my $counter = 0;
open(IN, "<".$infile_mapping) or die "Can't open file ".$infile_mapping."\n";
while(<IN>){
	chomp;

	next if($counter != 0 && substr($_, 0,1) eq "#");

	# First find the position in the spreadsheet of the required
	# fields.
	if($counter == 0){
		my @row = split(/\t/, $_);
		my $header_pos = 0;
		foreach my $el(@row){
			if($el eq "replicate"){
				$pos_replicate = $header_pos;
			}

			if($el eq "barcode"){
				$pos_barcode = $header_pos;
			}
			
			if($el eq "rep_common_name"){
				$pos_common_id = $header_pos;
			}
	
			if($el =~ m/SampleID/i){
				$pos_sample_id = $header_pos;
			}
	
			$header_pos++;
		}	
	
		$counter++;

	# Then extract the value and store it in data structure.
	}else{
		my @row = split(/\t/, $_);
		my $replicate = $row[$pos_replicate];
		my $barcode = $row[$pos_barcode];
		my $common_id = $row[$pos_common_id];
		my $sample_id = $row[$pos_sample_id];

		#print "replicate:\t".$replicate."\n";
		#print "barcode:\t".$barcode."\n";
		#print "common_id:\t".$id."\n";
		#print "sample_id:\t".$sample_id."\n";

		# Here populate two hash tables. 1st the hash_by_barcodes and 2nd hash_all		
		$hash_by_barcodes{$common_id}{$replicate} = $sample_id;
		$hash_all{$common_id}{$replicate} = $sample_id;

		$counter++;
	}		
}
close(IN);

# First loop through all replicates, regardless of their barcodes.
# Don't forget that we assume that all sequences in input fasta file are unique.
# First, do hash_by_barcodes and then pass a second time for hash_all.

my @tables;
# open file for writing. One outfile per common id.
push(@tables, $outfile_freq);
open(OUT_A, ">".$outfile_freq) or die "Can't open file ".$outfile_freq."\n";
my %curr_hash;
for my $common_id ( keys %hash_all ) {
	
	for my $replicate (keys % {$hash_all{$common_id}}){
		my $curr_sampleID =  $hash_all{$common_id}{$replicate};
		
		# Here loop through the otu table and match the position that match the sample name.
		# Then print it to file.
		open(OTU, "<".$infile_otu_table) or die "Can't open file ".$infile_otu_table."\n";
		my $pos;
		while(<OTU>){
			chomp;
			if($_ =~ m/#/){
				if($_ =~ m/#OTU ID/){
					my @row = split(/\t/, $_);
					shift(@row); pop(@row);
	
					my $header_counter = 0;
					my $found = 0;
					foreach my $el (@row){
						#print $el."\n";
						if($el eq $curr_sampleID){
							$pos = $header_counter;
							$found = 1;
						}
						$header_counter++;
					}
					print "Warning: Did not find sampleID in OTU table header...\n" if($found == 0);	
				}else{
					next;
				}
			}else{
				# Now extract the counts since we know the position of the curr. sample.
				if(defined($pos)){
					my @row = split(/\t/, $_);
					my $cluster_id = shift(@row); 
					my $lineage = pop(@row);
					#$curr_hash{$curr_sampleID}{$cluster_id} = $row[$pos];
					#$curr_hash{$cluster_id}{$curr_sampleID} = $row[$pos];
					if(exists $curr_hash{$cluster_id}{$replicate}){
						$curr_hash{$cluster_id}{$replicate} = $row[$pos] + $curr_hash{$cluster_id}{$replicate};
					}else{
						$curr_hash{$cluster_id}{$replicate} = $row[$pos];
					}
				}
			}
		}
		close(OTU);
	}
}
	
#print Dumper(\%curr_hash);	

# Now that we looped through all replicates of all common id, print the output to a file.
# Print header
print OUT_A "#OTU_ID";
my %names;
for my $cluster_id ( sort {$a<=>$b}  keys %curr_hash ) {
	for my $replicate ( sort {$a<=>$b} keys % {$curr_hash{$cluster_id}}){
		if(exists $names{$replicate}){
			}else{
			print OUT_A "\t".$replicate;
			$names{$replicate} = $replicate;
		}
	}
}
print OUT_A "\n";

# print counts.
for my $cluster_id ( sort {$a<=>$b} keys %curr_hash ) {
	print OUT_A $cluster_id;
	for my $replicate (sort {$a<=>$b} keys %{ $curr_hash{$cluster_id} } ){
		print OUT_A "\t".$curr_hash{$cluster_id}{$replicate};
	}
	print OUT_A "\n";
}

close(OUT_A);

# Next loop through clusters/OTUs and bin sequences based on their replicate.
# Don't forget that we assume that all sequences in input fasta file are unique.
# First, do hash_by_barcodes and then pass a second time for hash_all.
# Also print all results to a single table. (master_table.txt).
my $temp_table = $tmpdir."/temp_table.txt";
open(MASTER_TEMP, ">".$temp_table) or die "Can't open ".$temp_table;
print MASTER_TEMP "#OTU_ID\tRep1\tRep2\n";

my $max = 0;
for my $common_id ( keys %hash_by_barcodes ) {
	
	# open file for writing. One outfile per common id.
	open(OUT_B, ">".$outdir."/".$common_id.".tab") or die "Can't open file ".$outdir."/".$common_id.".tab";
	push(@tables, $outdir."/".$common_id.".tab");
	my %curr_hash;

	for my $replicate (keys % {$hash_by_barcodes{$common_id}}){
		my $curr_sampleID =  $hash_by_barcodes{$common_id}{$replicate};
		
		# Here loop through the otu table and match the position that match the sample name.
		# Then print it to file.
		open(OTU, "<".$infile_otu_table) or die "Can't open file ".$infile_otu_table."\n";
		my $pos;
		while(<OTU>){
			chomp;
			if($_ =~ m/#/){
				if($_ =~ m/#OTU ID/){
					my @row = split(/\t/, $_);
					shift(@row); pop(@row);
	
					my $header_counter = 0;
					my $found = 0;
					foreach my $el (@row){
						#print $el."\n";
						if($el eq $curr_sampleID){
							$pos = $header_counter;
							$found = 1;
						}
						$header_counter++;
					}
					print "Warning: Did not find sampleID ".$curr_sampleID." in OTU table header...\n" if($found == 0);	
				}else{
					next;
				}
			}else{
				# Now extract the counts since we know the position of the curr. sample.
				if(defined($pos)){
					my @row = split(/\t/, $_);
					my $cluster_id = shift(@row); 
					my $lineage = pop(@row);
					#$curr_hash{$curr_sampleID}{$cluster_id} = $row[$pos];
					$curr_hash{$cluster_id}{$curr_sampleID} = $row[$pos];
					$max = $row[$pos] if($row[$pos] > $max);
				}
			}
		}
		close(OTU);
	}
	
	# Now that we looped through all replicates of all common id, print the output to a file.
	print OUT_B "#OTU_ID";
	my %names;
	for my $cluster_id ( sort {$a<=>$b} keys %curr_hash ) {
		for my $curr_sampleID (keys % {$curr_hash{$cluster_id}}){
			if(exists $names{$curr_sampleID}){

			}else{
				print OUT_B "\t".$curr_sampleID;
				$names{$curr_sampleID} = $curr_sampleID;
			}
		}
	}
	print OUT_B "\n";

	for my $cluster_id ( sort {$a<=>$b} keys %curr_hash ) {
		print OUT_B $cluster_id;
		print MASTER_TEMP $cluster_id;

		for my $curr_sampleID (sort keys %{ $curr_hash{$cluster_id} } ){
			print OUT_B "\t".$curr_hash{$cluster_id}{$curr_sampleID};
			print MASTER_TEMP "\t".$curr_hash{$cluster_id}{$curr_sampleID};
		}
		print OUT_B "\n";
		print MASTER_TEMP "\n";
	}
	close(OUT_B);

}
close(MASTER_TEMP);

## Then loop again through table to remove lines where there are only 2 elements.
my $master_table = $outdir."/master_table.txt";
open(MASTER_TEMP, "<".$temp_table) or die "Can't open ".$temp_table;
open(MASTER, ">".$master_table) or die "Can't open ".$master_table;
print MASTER "#OTU_ID\tRep1\tRep2\n";
while(<MASTER_TEMP>){
	chomp;
	if($_ =~ m/#/){
		next;
	}

	my @row = split(/\t/, $_);
	next if(@row < 3);
	print MASTER $_."\n";
}
close(MASTER);

push(@tables, $master_table);
#
##print Dumper(\%hash_y);	

## BUILD GRAPH WITH R
my $log_max = log($max)/log(10);
my $R = Statistics::R->new();
$R->startR;
$R->run("options(stringsAsFactors = FALSE)");
$R->run("library(gplots)");

if($add_duplicates){
	pop(@tables);shift(@tables);
	my @string;
	$string[0] = "#OTU table having added duplicates";
	$string[1] = "#OTU ID\t";
	my $added_otu_table = $outdir."/added_duplicates.tab";
	open(ADD, ">".$added_otu_table) or die "Can't open ".$added_otu_table."\n";
	my $table_counter = 0;
	foreach my $table (@tables){
	
		# Test if table have at least 2 replicates. Length = 3 including sampleID.
		open(TABLE, "<".$table) or die "Can't open file ".$table."\n";
		my $counter = 0;
		my $row_length;
		while(<TABLE>){
			chomp;
			if($_ =~ m/#/){
				my @row = split(/\t/, $_);
				# print header.
				$string[1] .= $row[1]."\t";
				next;
			}
	
			my @row = split(/\t/, $_);
			$row_length = @row;
			if($table_counter == 0){
				# If no duplicates, print only value.
				if($row_length < 3){
					$string[$counter+2] .= $row[0]."\t".$row[1]."\t";
				}else{
					my $sum = $row[1] + $row[2];
					$string[$counter+2] .= $row[0]."\t".$sum."\t";
				}

			}else{
				# If no duplicates, print only value.
				if($row_length < 3){
					$string[$counter+2] .= $row[1]."\t";
				}else{
					my $sum = $row[1] + $row[2];
					$string[$counter+2] .= $sum."\t";
				}
			}
			$counter++;
		}
		$table_counter++;
	}

	# Finally print to file.
	foreach my $line (@string){
		$line = substr($line, 0, -1);
		print ADD $line."\n";
	}
	close(ADD);
	exit;
}


my $file_counter = 0;
foreach my $table (@tables){
	
	# Test if table have at least 2 replicates. Length = 3 including sampleID.
	open(TABLE, "<".$table) or die "Can't open file ".$table."\n";
	my $counter = 0;
	my $row_length;
	while(<TABLE>){
		chomp;
		if($counter == 1){
			my @row = split(/\t/, $_);
			$row_length = @row;
		}
		$counter++;
	}
	next if($row_length < 3);	

	my $curr_string = '
		
		## Functions START ##
		rsqr <- function(x) {
			## creates the new otu table with only otu\'s with more than the i reads (all done at log scale)
			com_test <- com_log_threshold[(com_log_threshold[,1] > log_table[x]) & (com_log_threshold[,2] > log_table[x]),]

			## fits the linear regression model
			fit <- lm(com_test[,1]~com_test[,2])
			
			## pulls out the R-squared statistic
			rsq_value <- summary(fit)[[8]]
		
			## puts the r-squared value in a pair, along with the number i (min number of sequences allowed)
			pair <- c(x, rsq_value)
			
			return(pair)
		}
		## Functions END ##

		data <- read.table("'.$table.'")
		rep1 <- (data[,2]+1)
		rep2 <- (data[,3]+1)
		maxRep1 <- max(rep1)
		maxRep2 <- max(rep2)
		log_rep1 <- log10(data[,2]+1)
		log_rep2 <- log10(data[,3]+1)
		data <- cbind(data, log_rep1)
		data <- cbind(data, log_rep2)
		com_log_threshold <- data[,4:5]
		
		y <- 50
		numiter <- seq(1:y)
		log_table <- log10(numiter)
		rsq_table <- matrix(c(0,0),nrow=1)
		
		## prepares the blank rsqtable
		com_log_threshold <- as.data.frame(com_log_threshold)
		
		## for 1 to 50... for instance.
		for (i in 1:y) {
			## all iterations 1 through y.
			pair<- rsqr(i) ## Call rsqr function here... Will determine R2 values at each log threshold
			## adds the returned pair to the rsq_table
			rsq_table <- rbind(rsq_table,pair)
		}
		
		## Plot tables
		## Before plotting, find the OTU count threshold that correspond to >= 0.87
		THRESHOLD='.$threshold.'
		myThreshold <- 0
		foundThreshold <- 0		

		for (j in 1:y){
			if(rsq_table[j,2] >= THRESHOLD){
				myThreshold <- rsq_table[j,1]
				foundThreshold <- 1
				break
			}
		}
		
		## Plot the R2 table itself
		gplots:::textplot(rsq_table)
		
		## this plots the R-squared values vs. count threshold with R value drawn at specified cutoff.
		plot(
			rsq_table[,1],
			rsq_table[,2],
			col="black",
			main="",
			axes=F,
			xlab="",
			ylab="",
			ylim=c(0,1),
			pch=20,
			cex=0.7
		)
		box(lwd=0.5)
		axis(side=1, las=1,   cex.axis=0.5, lwd=0.5)
		axis(side=2, las=1.5, cex.axis=0.5, lwd=0.5)
		mtext("R-squared value", side=2,line=2, cex=0.6)
		mtext("OTU count thresholds", side=1,line=2, cex=0.6)
		mtext("'.basename($table).'", side=3, cex=0.6)
		if(foundThreshold != 0){
			abline(v=myThreshold, col="red")
		}
		usr <- par("usr") ## get user coordinates
		par(usr = c(0, 1, 0, 1)) ## new relative user coordinates
		if(foundThreshold != 0){
			myString1 = paste("R-squared >", THRESHOLD, sep = " ", collapse = NULL)
			myString2 = paste("OTU threshold value = ", myThreshold, sep = " ", collapse = NULL)
			text(0.7, 0.2, labels=c(myString1), col="black", cex=0.6)
			text(0.7, 0.3, labels=c(myString2),col="black", cex=0.6)
		}
		par(usr = usr) ## restore original user coordinates
	
		plot(
			log_rep1, 
			log_rep2, 
			col="black",
			#log="xy",
			xlab="", 
			ylab="",
			main="",
			xlim=c(0, '.$log_max.'),
			ylim=c(0, '.$log_max.'),
			axes=F,
			pch=20,
			cex=0.35
		)
		box(lwd=0.5)
		axis(side=1, las=1,   cex.axis=0.5, lwd=0.5)
		axis(side=2, las=1.5, cex.axis=0.5, lwd=0.5)
		mtext("OTUs replicate 2", side=2,line=2, cex=0.6)
		mtext("OTUs replicate 1", side=1,line=2, cex=0.6)
		mtext("'.basename($table).'", side=3, cex=0.6)
		abline(h=log10(myThreshold),lwd=1,col="red")
		abline(v=log10(myThreshold),lwd=1,col="red")
	';
	
	if($file_counter == 0){
		$curr_string .= '
			write.table(rsq_table, "'.$r_squared_table.'", sep="\t")
		';
		#print $curr_string."\n";
	}
	$file_counter++;

	push(@array, $curr_string);
}

# Finally combine all plots on one plot.

$R->run('pdf("'.$outfile_graph.'" , width=7, height=7)');
$R->run('par(mfrow=c('.$row.','.$col.'))'); 
foreach my $cmd (@array){
	$R->run($cmd);
	#print $cmd."\n";
}
$R->run('dev.off()');

## Find OTU threshold based on r-squared table. Then 
## Generate a new otu table based on that threshold
#$threshold_OTU = findThreshold($r_squared_table, $threshold);
#
#my $cmd = "itaggerFilterObsTable.pl";
#$cmd .= " --otu_table ".$infile_otu_table;
#$cmd .= " --otu_table_out ".$outfile_otu_table;
#$cmd .= " --frequency 2";
#$cmd .= " --threshold ".$threshold_OTU;
##$cmd .= " --threshold 2";
#system($cmd);

## FIND OTU THREHSOLD
sub findThreshold{
	my ($file, $THRESHOLD) = @_;

	my $threshold = 0;

	open(IN, "<".$file) or die "Can't open file ".$file."\n";
	while(<IN>){
		chomp;
		#next if($. == 1 || $. == 2);
		my @row = split(/\t/, $_);
		if($row[2] >= $THRESHOLD){
			$threshold = $row[1];	
			last;
		}
	}
	close(IN);
	return $threshold;
}

## REMOVE TEMP FILES
sub END{
	system("rm ".$tmpdir." -rf");
}
1;
exit;
