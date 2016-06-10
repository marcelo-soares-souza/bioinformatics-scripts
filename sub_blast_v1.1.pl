#!/usr/bin/perl -w
#-------------------------------------------------------------------
# Split query (fasta) into pieces, Blast them all and concat output
# 
# Requirements:
# - fasta_from_list_v3.1.pl  
# 
# Inputs: 
# - parameters (through sub_blast.conf)
#
# Outputs:
# - Blast results
# 
# Lucas S. de Brito 2016/05/16
#
# Usage: 
# Configure config file <sub_blast.conf>
# Run: perl sub_blast.pl {config file}

# I have large RAM and prefere data loaded into memory => 1
# I have low   RAM and prefere data loaded into files  => 0
$ram = 1;   # Include into config, and update present script


#===================================================================
# CORE
#-------------------------------------------------------------------
# 1. Initial stuff
iniStuff();

# 2. Variables
vars();

# 3. Read config file
configScript();

# 4. Print chosen parameters into LOG
printParameters();

# 5. Create directories
createDirs();

# 6. Build command lines
buildCmdLines();
# 6.1. SUB set_blast_options => build blast options command line

# 7. Run Blast
runBlast();
# 7.1. SUB runBlastParallel => run Blast parallel in more than 2 threads

# Final stuff
$date = `date`;
print LOG "\n\nParallel Blast ' sub_blast.pl ' finished at " . $date;
print "\n\nParallel Blast ' sub_blast.pl ' finished :) at " . $date . "\n";
#===================================================================


#===================================================================
# SUBROUTINES
#-------------------------------------------------------------------

#-------------------------------------------------------------------
# 1. SUB iniStuff => Initial configuration
sub iniStuff
{
	$date = `date`;
	print "\nStarting script at " . $date;

	# Expected date format: "Wed Oct 14 12:48:33 BRT 2015"
	($temp,$mon,$mday,$fullhour,$temp,$year) = split(' ',$date);
	$fullhour =~ s/:/h/;     # replace first ':'
	$fullhour =~ s/:/min/;   # replace second ':'
	$starttime = $year.$mon.$mday."-".$fullhour;

	# Log file
	$log_file = "sub_blast_".$starttime.".log";
	$log_file2 = "temp.log2";
	open (LOG, ">$log_file") or die $!;
	print LOG "\nStarting script at " . $date;

	use YAPE::Regex::Explain;
}

#-------------------------------------------------------------------
# 2. SUB vars => Variables initialization etc
sub vars
{
	# Some default values:
	$path_to_softwares = "/usr/local/bin/";
	$path_to_results   = ".";   # Default
	$log_file          = "sub_blast.log";

	if ($ARGV[0])
	{
		$conf_file = $ARGV[0];
	}
	else
	{
		$conf_file = "sub_blast.conf";
	}

	# Script config file
	$config_file = $path_to_results."/".$conf_file;
}

#-------------------------------------------------------------------
# 3. SUB configScript => Load config file parameters
sub configScript {

	# Filling the hash with recomended values
	%config = 
	(
		NUM_SEQS    => '1',
		BLAST_THREADS => '1',
		BLAST         => 'blastn',
		FFL         => 'fasta_from_list_v3.1.pl',

		PATH_TO_SOFTWARES => '/usr/local/bin/',
		PATH_TO_RESULTS   => '.',

		BLAST_EVALUE         => '1e-5',
		BLAST_PERC_IDENTITY  => '0.70',
		BLAST_OUTPUT_FORMAT  => '6 qseqid sseqid qstart qend sstart send evalue length pident qcovs',

	);

	# Loading config file
	print "\nReading config file... ";
	open(CONFIG, "$config_file");

	while (my $line = <CONFIG>) 
	{
		chop($line);

		# se não for linha comentada e nem linha vazia
		if ($line =~ /^[^\#]/ && !($line =~ /^\s+$/)) 
		{
			my ($key, $value) = split("=", $line);
			$value =~ s/^\s+//g;
			$value =~ s/\s+$//g;
			$key   =~ s/^\s+//g;
			$key   =~ s/\s+$//g;

			if ( $value ne "" )
			{
				$config{$key} = $value;
			}
		}
	}
	close(CONFIG);
	print "done.";
}

#-------------------------------------------------------------------
# 4. SUB printparameters => Print parameters into LOG file
sub printParameters
{

	# Print at log the parameters.................................
	print LOG "\nDefined parameters\n";
	print LOG "-----------------------------\n";
	foreach my $k (sort keys %config)
	{
		print LOG "$k\: $config{$k}\n"
	}
	print LOG "-----------------------------\n";
}

#-------------------------------------------------------------------
# 5. SUB createDirs => Create directories
sub createDirs
{
	print "\nCreating directories... ";
	# Base directory
	system("mkdir $path_to_results/$starttime");
	# For blast
	system("mkdir $path_to_results/$starttime/blast/");
	for ($i6=1;$i6<=$config{NUM_SEQS};$i6++)
	{
		system("mkdir $path_to_results/$starttime/blast/gen" . $i6);
	}
	# other ***
	print "done.";
}

#-------------------------------------------------------------------
# 6. SUB buildCmdLines => Build command lines
sub buildCmdLines
{
	# path to software
	if ( $config{PATH_TO_SOFTWARES} ne ""  )
	{
		$path_to_softwares = $config{PATH_TO_SOFTWARES}; 
	}
	if ( $path_to_softwares !~ /\/$/ ) 
	{
		$path_to_softwares = $path_to_softwares . "\/"
	}

	# path to results
	if ( $config{PATH_TO_RESULTS} ne ""  )
	{
		$path_to_results = $config{PATH_TO_RESULTS}; 
	}
	if ( $path_to_results !~ /\/$/ ) 
	{
		$path_to_results = $path_to_results . "\/"
	}

	# number of blast threads
	if ( $config{BLAST_THREADS} ne ""  )
	{
		$blast_threads = $config{BLAST_THREADS}; 
	}

	# blast
	# Build options
	set_blast_options();

	# build blast command line for >=1 threads
	@cmd_blast = {};	
	for ($i9=1;$i9<=$config{NUM_SEQS};$i9++)
	{
		for ($i14=1;$i14<=$blast_threads;$i14++)
		{
			$fasta_temp = "temp_".$i14.".fasta";
			$cmd_blast{$i9}{$i14} = $config{BLAST}." ".$fasta_temp." ".$blast_options." 2> ".$log_file2." > ".$log_file2;
		}
	}
}

#-------------------------------------------------------------------
# 6.1. SUB set_blast_options => build blast options command line
sub set_blast_options 
{
	if ($config{"TRF_MASKED_FILE"} eq "yes")
	{
		$trf_m = " -m"; 
	}
	else 
	{ 
		if ($config{"TRF_MASKED_FILE"} eq "no")
		{
			$trf_m = ""; 
		}
		else
		{
			die "ERROR: TRF_MASKED_FILE must be \"yes\" or \"no\"."; 
		}
	}

	if ($config{"TRF_FLANKING_SEQUENCE"} eq "yes")
	{
		$trf_f = " -f"; 
	}
	else 
	{
		if ($config{"TRF_FLANKING_SEQUENCE"} eq "no")
		{
			$trf_f = ""; 
		}
		else
		{ 
			die "ERROR: TRF_FLANKING_SEQUENCE must be \"yes\" or \"no\"."; 
		}
	} 

	if ($config{"TRF_DATA_FILE"} eq "yes")
	{ 
		$trf_d = " -d"; 
	}
	else 
	{
		if ($config{"TRF_DATA_FILE"} eq "no")
		{ 
			$trf_d = ""; 
		}
		else
		{ 
			die "ERROR: TRF_DATA_FILE must be \"yes\" or \"no\"."; 
		}
	}

	if ($config{"TRF_SUPPRESS_HTML"} eq "yes")
	{ 
		$trf_h = " -h"; 
	}
	else 
	{
		if ($config{"TRF_SUPPRESS_HTML"} eq "no")
		{ 
			$trf_h = ""; 
		}
		else
		{ 
			die "ERROR: TRF_SUPPRESS_HTML must be \"yes\" or \"no\"."; 
		}
	}

	if ($config{"TRF_NGS"} eq "yes")
	{ 
		$trf_ngs = " -ngs"; 
	}
	else 
	{
		if ($config{"TRF_NGS"} eq "no")
		{ 
			$trf_ngs = ""; 
		}
		else
		{ 
			die "ERROR: TRF_NGS must be \"yes\" or \"no\"."; 
		}
	}

	$blast_options = " " . $config{"TRF_MATCH_INDEL_WEIGHT"} . " " . $config{"TRF_MISMATCH_INDEL_WEIGHT"} . " " . $config{"TRF_DELTA_INDEL_WEIGHT"} . " " . $config{"TRF_PM"} . " " . $config{"TRF_PI"} . " " . $config{"TRF_MIN_SCORE"} . " " . $config{"TRF_MAX_PERIOD"} . $trf_m . $trf_f . $trf_d . $trf_h . $trf_ngs;
}

#-------------------------------------------------------------------
# 7.1. SUB runBlastParallel => run Blast in 1 or more threads
sub runBlastParallel
{

	###################################################################################################
	### UNDER CONSTRUCTION INI

	# ------------------------------------------------------------------
	# *** model for parallel processing of multifastas *** (unfinished)
	# ------------------------------------------------------------------

	# split genome file and run blast in all temp files at the same time

	# creating list from fasta
	$fasta_list = $fasta_file;
	$fasta_list =~ s/.fasta$/.list/;
	$cmd_less = "less " . $fasta_path . $fasta_file . " | grep \">\" > " . $fasta_list . "";
	system( $cmd_less );

	# loading and counting list
	@HLIST = {};
	open ( LIST, "$fasta_list" );
	$cnt_list = 1;
	foreach $line (<LIST>)
	{
		chomp $line;
		$HLIST{$cnt_list} = $line;
		$cnt_list++;
	}
	$num_fasta = $cnt_list - 1;
	print LOG "\n$fasta_file have $num_fasta fastas\n";

	# split fasta file in up to $blast_threads temp fastas
	# create temp file names
	@fasta_temp = {};
	for ($i2=1;$i2<=$blast_threads;$i2++)
	{
		if ( $i2 <= $num_fasta )
		{
			$list_temp{$i2}  = "temp_" . $i2 . ".list";
			$LIST_TEMP      = "LIST_TEMP_" . $i2;
			open ( $LIST_TEMP, ">$list_temp{$i2}" ) or die $!;
		}
	}
	# Filling list temp file lists
	$thread_control = 1;
	for ($i3=1;$i3<=$num_fasta;$i3++)
	{
		print ".";
		$fasta_temp{$thread_control} = "temp_" . $thread_control . ".fasta";            
		$LIST_TEMP       = "LIST_TEMP_" . $thread_control;
		# printing fasta names into files
		# ---------------------------------------------------------------------------------------
		print $LIST_TEMP "".$HLIST{$i3}."\n";   # $LIST_TEMP = filename; $HLIST{$i3} = fasta name
		# ---------------------------------------------------------------------------------------
		if ( $thread_control < $blast_threads )
		{
			$thread_control++;
		}
		else
		{
			$thread_control = 1;
		}
	}

	# closing temp file lists
	for ($i4=1;$i4<=$blast_threads;$i4++)
	{
		if ( $i4 <= $num_fasta )
		{
			$temp = "TEMP_" . $i4;
			close( $temp );
		}
	}

	# create temp file fastas
	for ($i12=1;$i12<=$blast_threads;$i12++)
	{
		$out_temp = $list_temp{$i12};
		$out_temp =~ s/.list$/.fasta/;
		$cmd_ffl = "perl ".$path_to_softwares.$config{FFL}." ".$list_temp{$i12}." ".$fasta_path.$fasta_file." ".$out_temp;
		print "=";
		system( $cmd_ffl );
	}

	# RUN BLAST in all fastas
	# run all at same time and save pid numbers
	@childpid = {};
	for ($i13=1;$i13<=$blast_threads;$i13++)
	{
		print ">";
		print "\n \=\> You can follow this blast analysis in another terminal by\n \=\> ' tail -f ".$cmd_cd."/".$log_file2." '\n";
		$childpid{$i1}{$i13} = fork() or exec("$cmd_blast{$i1}{$i13}");
	}
	# wait all pids to finish
	for ($i14=1;$i14<=$blast_threads;$i14++)
	{
		print LOG "\nwaiting PID $childpid{$i1}{$i14} from $cmd_blast{$i1}{$i14}\n";
		waitpid( $childpid{$i1}{$i14},0 );
	}
	print LOG "... genome ".$i1." done.";
	print "\n\n*********************************************************************************************";   # debug
	print "\n*** Remember to delete list and fasta temporary from Parallel blast at the end of scripts ***";     # debug
	print "\n*********************************************************************************************\n";   # debug

	#------------------------------------#
	# RUNNIG WELL SO FAR, BUT INCOMPLETE #
	#------------------------------------#
	# NEXT STEPS:
	# => JOIN/RENAME FILES FROM BLAST TO RECOGNIZE THEM AS ONE BLAST RESULT
	# => DELETE TEMPORATY FILES
	# ----------------------------------

	### UNDER CONSTRUCTION END
	###################################################################################################

	# end things
	chdir($path_to_results) or die "Cannot change dir to \'$cmd_cd\', error $!";
	print "genome ".$i1." done.\n";
	print LOG "... genome ".$i1." done.";
}
