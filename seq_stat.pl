#!/usr/bin/perl -w
use strict;
my ($len,$total,$contigs,$GC,$Ns,$gaps)=(0,0,1,0,0,0);
my @x; while(<>){
        if(/^[\>\@]/){
                if($len>0){
                        $total+=$len;
                        $contigs ++;
                        push @x,$len;
                }
                $len=0;
        }
        else{
                s/\s//g;
                $len+=length($_);
        $GC++ while $_ =~ /[GC]/gi;
        $Ns++ while $_ =~ /N/gi;
        $gaps++ while $_ =~ /N[ACTG]/gi;
        }
}
if ($len>0){
        $total+=$len;
        push @x,$len;
}
@x=sort{$b<=>$a} @x;
my $max_value = $x[0];
my ($count,$half)=(0,0);
for (my $j=0;$j<@x;$j++){
        $count+=$x[$j];
        if (($count>=$total/2)&&($half==0)){
                print "Number of sequences\t$contigs scaffolds\n";
                print "Total length\t$total bp\n";
                print "Max scaffold length\t$max_value bp\n";
                print "N50\t$x[$j] bp\n";
                $half=$x[$j]
        }elsif ($count>=$total*0.9){
                print "N90\t$x[$j] bp\n";
                my $roundGC = sprintf ("%.2f",($GC / $total) * 100);
                print "GC content\t$roundGC%\n";
                print "# of gaps\t$gaps\n";
                my $roundN = sprintf ("%.2f",($Ns / $total) * 100);
                print "% of genome in gaps\t$roundN\n";
                exit;
        }
}
