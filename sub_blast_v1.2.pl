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
runBlastSubfiles();
# 7.1. SUB runBlastSubfiles => run Blast sub files

# Final stuff
$date = `date`;
print LOG "\n\nSub Files Blast ' sub_blast.pl ' finished at " . $date;
print "\n\nSub Files Blast ' sub_blast.pl ' finished :) at " . $date . "\n";
#===================================================================


#===================================================================
# SUBROUTINES
#-------------------------------------------------------------------

#-------------------------------------------------------------------
# 1. SUB iniStuff => Initial configuration
sub iniStuff
{
	my $padded_num = '000';
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

	#use YAPE::Regex::Explain;
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
		BLAST_SUBFILES => '1',
		BLAST          => 'blastn',
		FFL            => 'fasta_from_list_v3.1.pl',

		PATH_TO_SOFTWARES => '/usr/local/bin/',
		PATH_TO_RESULTS   => '.',

		BLAST_EVALUE         => '1e-5',
		BLAST_PERC_IDENTITY  => '0.70',
		#BLAST_OUTPUT_FORMAT  => '6 qseqid sseqid qstart qend sstart send evalue length pident qcovs',

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
	system("mkdir $path_to_results/$starttime/blast/gen");
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

	# number of blast sub files
	if ( $config{BLAST_SUBFILES} ne ""  )
	{
		$blast_subfiles = $config{BLAST_SUBFILES}; 
	}

	# blast
	# Build options
	set_blast_options();

	# build blast command line for >=1 sub files
	@cmd_blast = {};	
	for ($i14=1;$i14<=$blast_subfiles;$i14++)
	{
		$padded_num = sprintf ("%0${num_digits}d", $i14 );
		$fasta_temp = "temp_".$padded_num.".fasta";
		$blast_out = "$path_to_results/$starttime/blast/gen/" . $fasta_temp;
		$blast_out =~ s/\.fasta/$config{"BLAST_OUT_SUFFIX"}/;	
		$cmd_blast{$i14} = $config{BLAST}." -query " . $fasta_temp . " -out " . $blast_out . " -db " . $config{"BLAST_DATABASES_ADDRESS"} . $blast_options . " 2> ".$log_file2." > ".$log_file2;
	}
}

#-------------------------------------------------------------------
# 6.1. SUB set_blast_options => build blast options command line
sub set_blast_options 
{
	$blast_options = "";

	if($config{"BLAST_REWARD"} ne "disable"){
		$blast_options = $blast_options . " -reward " . $config{"BLAST_REWARD"};
	}

	if($config{"BLAST_PENALTY"} ne "disable"){
		$blast_options = $blast_options . " -penalty " . $config{"BLAST_PENALTY"};
	}
	
	if($config{"BLAST_GAPOPEN"} ne "disable"){
		$blast_options = $blast_options . " -gapopen " . $config{"BLAST_GAPOPEN"};
	}
	
	if($config{"BLAST_GAPEXTEND"} ne "disable"){
		$blast_options = $blast_options . " -gapextend " . $config{"BLAST_GAPEXTEND"};
	}
	
	if($config{"BLAST_DUST"} ne "disable"){
		$blast_options = $blast_options . " -dust " . $config{"BLAST_DUST"};
	}

	if($config{"BLAST_SOFT_MASKING"} ne "disable"){
		$blast_options = $blast_options . " -soft_masking " . $config{"BLAST_SOFT_MASKING"};
	}

	if($config{"BLAST_EVALUE"} ne "disable"){
		$blast_options = $blast_options . " -evalue " . $config{"BLAST_EVALUE"};
	}
	
	if($config{"BLAST_SEARCHSP"} ne "disable"){
		$blast_options = $blast_options . " -searchsp " . $config{"BLAST_SEARCHSP"};
	}

	if($config{"BLAST_NUM_THREADS"} ne "disable"){
		$blast_options = $blast_options . " -num_threads " . $config{"BLAST_NUM_THREADS"};
	}

	if($config{"BLAST_PERC_IDENTITY"} ne "disable"){
		$blast_options = $blast_options . " -pident " . $config{"BLAST_PERC_IDENT"};
	}
	
	if($config{"BLAST_NUM_ALIGNMENTS"} ne "disable"){
		$blast_options = $blast_options . " -num_alignments " . $config{"BLAST_NUM_ALIGNMENTS"};
	}

	if($config{"BLAST_OUTPUT_FORMAT"} ne "disable"){
		$blast_options = $blast_options . " -outfmt " . $config{"BLAST_OUTPUT_FORMAT"};
	}

		

}

#-------------------------------------------------------------------
# 7.1. SUB runBlastSubfiles => run Blast in 1 or more sub files
sub runBlastSubfiles
{

	###################################################################################################
	### UNDER CONSTRUCTION INI

	# ------------------------------------------------------------------
	# *** model for sub files processing of multifastas *** (unfinished)
	# ------------------------------------------------------------------

	# RUN blast at $path_to_results/<date>/blast/<gen$i>/
	print "\nStarting Blast at genome ...";
	print LOG "\nStarting Blast at genome ...";
	$cmd_cd = $path_to_results . $starttime . '/blast/gen';
	chdir( $cmd_cd ) or die "Cannot change dir to \'$cmd_cd\', error $!";

	# Setting gemome fasta file
	$fasta = "GEN_FILE";
	$fasta_path = $path_to_results;
	$fasta_file = $config{$fasta};


	# Split genome file and run blast in all temp files at the same time:
	$fasta_list = $fasta_file;
	$fasta_list =~ s/.fasta$/.list/;
	$num_digits = length($blast_subfiles);

	if($config{"FLAG_READS"} eq '0'){
		# creating list from fasta
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
	
		# split fasta file in up to $blast_subfiles temp fastas
		# create temp file names
		@fasta_temp = {};
		for ($i2=1;$i2<=$blast_subfiles;$i2++)
		{
			if ( $i2 <= $num_fasta )
			{
				$padded_num = sprintf ("%0${num_digits}d", $i2 );
				$list_temp{$i2}  = "temp_" . $padded_num . ".list";
				$LIST_TEMP      = "LIST_TEMP_" . $padded_num;
				open ( $LIST_TEMP, ">$list_temp{$i2}" ) or die $!;
			}
		}
		# Filling list temp file lists
		$subfile_control = 1;
		for ($i3=1;$i3<=$num_fasta;$i3++)
		{
			print ".";
			$padded_num = sprintf ("%0${num_digits}d", $subfile_control );
			$fasta_temp{$subfile_control} = "temp_" . $padded_num . ".fasta";            
			$LIST_TEMP       = "LIST_TEMP_" . $subfile_control;
			# printing fasta names into files
			# ---------------------------------------------------------------------------------------
			print $LIST_TEMP "".$HLIST{$i3}."\n";   # $LIST_TEMP = filename; $HLIST{$i3} = fasta name
			# ---------------------------------------------------------------------------------------
			if ( $subfile_control < $blast_subfiles )
			{
				$subfile_control++;
			}
			else
			{
				$subfile_control = 1;
			}
		}
	
		# closing temp file lists
		for ($i4=1;$i4<=$blast_subfiles;$i4++)
		{
			if ( $i4 <= $num_fasta )
			{
				$padded_num = sprintf ("%0${num_digits}d", $i4 );
				$temp = "LIST_TEMP_" . $padded_num;
				close( $temp );
			}
		}
	

		# create temp file fastas:
		for ($i12=1;$i12<=$blast_subfiles;$i12++)
		{
			$out_temp = $list_temp{$i12};
			$out_temp =~ s/.list$/.fasta/;
			$cmd_ffl = "perl ".$path_to_softwares.$config{FFL}." ".$list_temp{$i12}." ".$fasta_path.$fasta_file." ".$out_temp;
			print "=";
			system( $cmd_ffl );
		}
	}
	else{
		# create temp file fastas:
		$cmd_wc = "wc -l $fasta_path.$fasta_file";
		$fasta_file_size = `$cmd_wc`;
		$temp_file_size = $fasta_file_size / $blast_subfiles;
		$cmd_split_reads = "split -a $num_digits -d -l $temp_file_size --additional-suffix=.fasta $fasta_file temp_";
		system($cmd_split_reads);
	}

	# RUN BLAST in all fastas
	# run all at same time and save pid numbers
	@childpid = {};
	for ($i13=1;$i13<=$blast_subfiles;$i13++)
	{
		print ">";
		print "\n \=\> You can follow this blast analysis in another terminal by\n \=\> ' tail -f ".$cmd_cd."/".$log_file2." '\n";
		print "$cmd_blast{$i13}\n";
		system("$cmd_blast{$i13}");
	}

	# Concatenate BLAST results
	$blast_final = $fasta_file;
	$blast_final =~ s/\.fasta/$config{"BLAST_OUT_SUFFIX"}/;
	for ($i15=1;$i15<=$blast_subfiles;$i15++)
	{
		$padded_num = sprintf ("%0${num_digits}d", $i15 );
		$fasta_temp = "temp_".$padded_num.".fasta";
		$blast_out = "$path_to_results/$starttime/blast/gen/" . $fasta_temp;
		$blast_out =~ s/\.fasta/$config{"BLAST_OUT_SUFFIX"}/;	
		system("cat $blast_out >> $blast_final");
	}

	print LOG "... genome done.";
	print "\n\n*********************************************************************************************";   # debug
	print "\n*** Remember to delete list and fasta temporary from blast sub files at the end of scripts ***";     # debug
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
	print "genome done.\n";
	print LOG "... genome done.";
}
