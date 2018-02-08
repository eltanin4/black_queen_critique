#!/usr/local/bin/perl -w

#this script should take as input a file includes organism with their EC numbers and as output generate a matrix of competition scores.
#Input: organism index and a list of ECS separate with space in each line i.e.: 1 4.2.1.10 2.6.1.11 1.1.1.22 2.7.7.18 1.1.1.31 4.2.1.19 


#------------------------main------------------------------------
use strict;
use lib ('./'); #This should point to where you keep NetSeed.pm
use NetSeed;
use File::Copy;
use File::Path qw(remove_tree);


my $in_ec_file = $ARGV[0]; #file of the ec numebrs per organism

my $networkFile = "network.txt";
my $enzymes_labels="enzymes_labels.txt"; #for each ec lable hold ec number (total 2283 enzymes): i.e. 5|1.1.1.10
my $reactions_file="reactions_parsed.txt"; #for each reaction lable holds changes between metabolits lables i.e. 1|=| |1:3|1:4
my $ec_reac_mapping="ec_reac_mapping.txt"; #for each ec lable holds the reaction lable i.e. 5|1: 1131  10|2: 1995 1355

my $paramFile = "parameters.txt";
my $maxSpecies = 500;
my %ec2LabelHash; #for each ec number holds ec label
my %ecLabel2ReactionHash; #hash of arrays for each ec label holds array of reaction
my %reaction2MetabolitesHash; #hash of array for each reaction holds array of metabolitse changes

my $seedFile = "seeds.txt";
my $genomes_nums_used = "genomes_nums_used.txt";


my $run_dir = "Runs/";
mkdir ($run_dir);

print "Reading data files\n";
&readFiles();
print "DONE reading input files\n";
print "Calculate seeds\n";
&runLoopNetSeed();
print "DONE calculate seeds\n";
print "Run NetCmpt\n";
&runSynergy();

print "Run NetCmpt finished\n";
my $tmp = "copy emo_data.txt ..\\Competition_matrix.txt";
system($tmp);

chdir ("..\\");
remove_tree($run_dir);

print "DONE ALL\n";
#-------------------------------------------sub--------------------------------
sub readFiles
{
	my ($line,$ecLabel,$nReaction,$allReaction,$reac,$from, $to,$direction);
	my (@reactionsArr);
	open EcLABLES, $enzymes_labels  or die "cannot open file $enzymes_labels  for reading $!";
	print "Reading $enzymes_labels\n";
	while ($line = <EcLABLES>) 
	{
		chomp $line;
		if ($line=~m/(\d+)\|(.+)/){$ec2LabelHash{$2} = $1;}
	}
	
	close EcLABLES;
	
	print "Reading $ec_reac_mapping\n";
	open EcRECTION, $ec_reac_mapping  or die "cannot open file $ec_reac_mapping  for reading $!";
	while ($line = <EcRECTION>) 
	{
		chomp $line;
		if ($line=~m/(\d+)\|(\d+)\:\s(.+)/)
		{
			$ecLabel = $1;
			$nReaction  = $2;
			$allReaction = $3;
			@reactionsArr = split (/\s/, $3);
			if (@reactionsArr != $nReaction) {print "error in parse $ec_reac_mapping ec $ecLabel\n";}
			
			$ecLabel2ReactionHash{$1} = [@reactionsArr];
			
		
		}
		
	}
	
	close EcRECTION;
	
	print "Reading $reactions_file\n";
	open RectMETABOLITES, $reactions_file  or die "cannot open file $reactions_file  for reading $!";
	while ($line = <RectMETABOLITES>) 
	{
		chomp $line;
		my @changes;
		my @vals = split(/\|/, $line);
		$reac = $vals[0];
		$direction =  $vals[1];
		$from = $vals[3];
		$to  = $vals[4];
		
		my @valsFrom = split (/:/, $from);
		my @metaboFromArr = split (/,/, $valsFrom[1]);
		my @valsTo= split (/\:/, $to);
		my @metaboToArr = split (/,/, $valsTo[1]);
		for (my $i = 0 ;$i<scalar(@metaboFromArr);$i++)
		{
			for (my $j = 0;$j<scalar(@metaboToArr); $j++)
			{				
				push @changes, "$metaboFromArr[$i]\t$metaboToArr[$j]";
				if ($direction eq "=") { push @changes, "$metaboToArr[$j]\t$metaboFromArr[$i]";}		
			}
		}
		
		$reaction2MetabolitesHash{$reac} =[@changes]; 
		
	}
	
	close RectMETABOLITES;
		
}


sub runLoopNetSeed
{

	chdir ($run_dir);
	$in_ec_file = "../" . $in_ec_file;
	open IN, $in_ec_file or die "Can't find $in_ec_file: $!\n";
	open (OUT,">$genomes_nums_used"); #to be used in synergic
	open (SEED,">$seedFile"); #to be used in synergic
	my @species;
	my %seeds;
	my %line;
	my ($ec,$label);
	my $org_index = 0;
	my $org;
	while (<IN>)
	{
		chomp $_;
		my @vals=split(/\s+/,$_);
		if (scalar(@vals)>0)
		{
			$org_index++;
			$org = $vals[0];
			if (!($org eq $org_index))
			{			
				#print "index of org $org changed to $org_index\n";
				$org = $org_index;
			}
			push @species,$org;
			
		}
		my $counter=scalar(@vals)-1;
		#$line{$org}= "$org|$counter: ";
		my $ecNum = 0;
		$line{$org}= "";
		$seeds{$org} = "$org|";
		for (my $id=1;$id<scalar(@vals);++$id){
			$ec = $vals[$id];		
			if (exists $ec2LabelHash{$ec})
			{
				$label = $ec2LabelHash{$ec};
				$line{$org}.= "$label ";
				$ecNum++;
			}
			else {print "$ec  of $org  is missing\n;";}
		}
		$line{$org} = "$org|$ecNum: ".$line{$org};
		&runNetSeed($line{$org},$networkFile,\%seeds,$org);		
		$line{$org}.= "\n";
	}
	
	
	close IN;
	print OUT "number_of_organisms=" . scalar(@species) . "\n";
	foreach my $species(@species){print OUT "$line{$species}";}
	close OUT;
	print SEED "number_of_organisms=" . scalar(@species) . "\n";
	foreach my $species(@species){print SEED "$seeds{$species}";}
	close SEED;
	if (scalar(@species)>$maxSpecies)
	{
		print "Numebr of species must be less then $maxSpecies \n";
		exit;
	}

}


sub runNetSeed
{
	my $inStr = shift;
	my $networkFile = shift;
	my $ref_seedHash = shift;
	my $org = shift;
	open (OUTFILE,">$networkFile") or die "unable to open $networkFile for writing\n";
	my @lablesArr = split (/\s/, $inStr);

	my (@reactionArr,@changes);
	for (my $i = 1; $i <scalar(@lablesArr);++$i)
	{
		my $label = $lablesArr[$i];
		if (exists $ecLabel2ReactionHash{$label})
		{
			
			
			@reactionArr = @{$ecLabel2ReactionHash{$label}};
			my $size  = scalar(@reactionArr);
			
			for (my $r = 0; $r<scalar(@reactionArr); $r++)
			{
				
				if (exists $reaction2MetabolitesHash{$reactionArr[$r]})
				{
					@changes = @{$reaction2MetabolitesHash{$reactionArr[$r]}};
					for (my $c = 0; $c<scalar(@changes); $c++){print OUTFILE "$changes[$c]\n"; }
				}
				else{print "reaction of lable $label is missing in createNet function\n";}
			}
			
			
		}
		else{print "Label $label is missing in createNet function\n";}
		
	}
	
	close OUTFILE;
	&calacSeed($networkFile,$ref_seedHash, $org); # running NetSeedPerl
	
}



sub calacSeed
{	
	my $networkFile = shift;
	my $ref_seedHash = shift;
	my $org = shift;
	open (my $NetFH, "< $networkFile") or die "Can't find $networkFile: $!\n";
	my $sep = "\t";
	my $giantComponent = 1;
	my $minimal = 0;

	# --- -- --- -- --- -- Analyze Network -- --- -- --- -- --- -- --- #

	CalculateSeeds($NetFH,$sep,$giantComponent,$minimal);

	# Grab the resulting data
	my %Seeds = Seeds();
	my %GroupedSeeds = GroupedSeeds();
	my %NonSeeds = NonSeeds();
	my %Ignored = IgnoredNodes();
	my %Nodes = AllNodes();
	my %Edges = AllEdges();
	my $nElements = NumElements();

	# --- -- --- -- --- -- --- --- --- --- --- -- --- -- --- -- --- -- #


	# --- -- --- -- --- -- Write Out Files -- --- -- --- -- --- -- --- #

	# Write out the files
	my $prefix = "example/";
	my $suffix = '.dat';

	unless (-d $prefix) {
  	mkdir $prefix; #change mkdir??
	}
	print "writing seeds to $prefix\n";

	# Seeds
	my $seedLink = "${prefix}seed${suffix}";
	my $numOfSeeds = keys( %Seeds );
	$ref_seedHash->{$org} .= "$numOfSeeds: ";
	open(seedFH, "> ${seedLink}") or print "<br />error opening ${seedLink}<br />";
	foreach my $seed (keys %Seeds) { print seedFH "$seed\t$Seeds{$seed}\n"; $ref_seedHash->{$org} .= "$seed "; }
	$ref_seedHash->{$org} .= "\n";
	close(seedFH);

	# Grouped Seeds
	my $groupedLink = "${prefix}seed_grouped${suffix}";
	open(groupFH, "> ${groupedLink}") or print "<br />error opening ${groupedLink}<br />";
	foreach my $seed (keys %GroupedSeeds) { print groupFH "$seed\n"; }
	close(groupFH);

	 # Non Seed Nodes
	my $nonLink = "${prefix}non_seed${suffix}";
	open(nonFH, "> ${nonLink}") or print "<br />error opening ${nonLink}<br />";
	foreach my $node (keys %NonSeeds) { print nonFH "$node\n"; }
	close(nonFH);

	# Pruned Nodes
	my $prunedLink = "${prefix}ignored${suffix}";
	open(prunedFH, "> ${prunedLink}") or print "<br />error opening ${prunedLink}<br />";
	foreach my $node (keys %Ignored) { print prunedFH "$node\n"; }
	close(prunedFH);
	
}

sub runSynergy
{
	
	open (PARAM_OUTFILE, ">$paramFile") or die "Failed PARAM_OUTFILE!!!\n";
	
	print PARAM_OUTFILE "enzymes_labels=../enzymes_labels.txt\n";
	print PARAM_OUTFILE "reactions_labels=../reactions_labels.txt\n";
	print PARAM_OUTFILE "metabolites_labels=../compounds_labels.txt\n";
	print PARAM_OUTFILE "reactions=../reactions_parsed.txt\n";
	print PARAM_OUTFILE "network_matrix=ReactionsMatrix.txt\n";
	print PARAM_OUTFILE "network_matrix_enzymes=EnzymesMatrix.txt\n";
	print PARAM_OUTFILE "ec_reac_mapping=../ec_reac_mapping.txt\n";
	print PARAM_OUTFILE "reac_ec_mapping=../reac_ec_mapping.txt\n";
	print PARAM_OUTFILE "genomes=$genomes_nums_used\n";
	print PARAM_OUTFILE "biomass=../biomass_vector.txt\n";
	print PARAM_OUTFILE "seeds=$seedFile\n";
	print PARAM_OUTFILE "create_reactions_matrix=Y\n";
	print PARAM_OUTFILE "print_file_name_pre=output=\n";
	print PARAM_OUTFILE "m_print_log=Y\n";
	print PARAM_OUTFILE "m_print_results=Y\n";
	print PARAM_OUTFILE "m_print_metabolism=Y\n";
	print PARAM_OUTFILE "m_print_nets=Y\n";
	print PARAM_OUTFILE "m_print_every=20\n";
	
	close(PARAM_OUTFILE);

	# RUN THE PROGRAM
	print "Running the NetCmpt program\n";
	my $out = `..//NetCmpt.exe`;
	
}