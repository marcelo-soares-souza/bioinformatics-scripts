#!/usr/bin/perl -w
# Quick script to count fasta sequences frequency and lenght
#
# Edu :) 20161018
#
# Usage: perl <script_name>.pl <fasta_file_name>
#===========================================================

if ( !$ARGV[0] )
{       # input files basic control
  die "\nUsage: perl <script_name <fasta_file_name>\n\n";
}

$date = `date`;

$fasta = $ARGV[0];   # capture fastq file name
$out   = $fasta . ".freq.csv";

# open ( FA, "$fasta" ) or die "\nCannot open fasta file\n";

open (FA, "gunzip -c $fasta | ") or die "\nCannot open fasta file\n";
open ( OUT, ">$out" ) or die "\nCannot open output\n";

print "\nStarted file $fasta at $date\n";

my %HSEQ  = ();

while ( my $line = <FA> )
{
  chomp $line;
  if ( $line !~ /^>.*$/ )
  {
    $HSEQ{$line}++;
  }
}

print OUT "Sequence,Freq,Len\n";

foreach $k (sort keys %HSEQ)
{
  $len  = length($k);
  print OUT "$k\,$HSEQ{$k}\,$len\n";
}

$date = `date`;
print "Finished file $fasta at $date\n";
print "Output file name: ' $out '\n\n";

