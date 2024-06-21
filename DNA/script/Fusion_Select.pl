#!/usr/bin/perl -w
use strict;
use warnings;
no strict 'refs';

die "$0 fusions.txt samplename " if @ARGV < 2 ;

open FUS,'<',"$ARGV[0]" or die "$!\n";

while (<FUS>){
	chomp;
	my $line = $_;
	if (/^Est_Type/){
		print "sampleName\tGene1\tGene2\tBreak1\tBreak2\tFrequency\n";
	}
	else{
		my @tmp = split /\t/,$line;
		my $gene1 = $tmp[1];
		my $gene2 = $tmp[2];
		my $break1 = $tmp[3];
		my $break2 = $tmp[4];
		my $support1 = $tmp[5];
		my $support2 = $tmp[6];
		my $total_depth = $tmp[16];
		my $freq_tmp1 = ($support1 * 4 / ($total_depth + $support1 * 3)) * 100 ;
		my $freq_tmp2 = ($support2 * 4 / ($total_depth + $support2 * 3)) * 100 ;
		my $frequency;
		$freq_tmp1 >= $freq_tmp2 ? ($frequency = (sprintf "%.2f",$freq_tmp1) . '%') : ($frequency = (sprintf "%.2f",$freq_tmp2) . '%');
#		$freq_tmp1 >= $freq_tmp2 ? print "$freq_tmp1\n" : print "$freq_tmp2\n";
		if (grep /RET/,@tmp){
			print "$ARGV[1]\t$gene1\t$gene2\t$break1\t$break2\t$frequency\n";
		}
	}
}
close(FUS);
