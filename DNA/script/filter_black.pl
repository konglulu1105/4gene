#!/usr/bin/perl -w 
use strict;

die "$0 annotation.txt black.list sample_black.txt\n" if @ARGV < 3 ;

open ANNO,'<',"$ARGV[0]" or die "can not open file annotation.txt\n";
open LIST,'<', "$ARGV[1]" or die "can not open file black.list\n";
open BLACK,'>',"$ARGV[2]" or die "can not open file sample_black.txt\n";
my @black_list ;
while (<LIST>){
	chomp;
	my $line = $_;
	my @tmp = split /\t/,$line;
	my $key = join "\t",@tmp[0..5];
	push @black_list,$key;
}
while (<ANNO>){
	chomp;
	my $line = $_;
	if (/^SampleName/){
		print "$line\n";
		print BLACK "$line\n";
	}
	else{
		my @tmp = split /\t/,$line;
		my $key = join "\t",@tmp[1,2,3,4,5,10];
		if (grep /^$key$/,@black_list){
			print BLACK "$line\n"
		}
		else{
			print "$line\n";
		}
	}
}
