#!/usr/bin/perl -w
use strict;
use warnings;

die "$0 sample_varscan_snp.vcf sample_varscan_indel.vcf " if @ARGV < 2 ;

open SNP,'<',"$ARGV[0]" or die "$!\n";
open INDEL,'<',"$ARGV[1]" or die "$!\n";
my $snp_number = 0;
my $indel_number = 0;

while (<SNP>){
    chomp;
    my $line = $_;
    if(/^#/){
        print "$line\n";
    }
    if (/^chr/){
	$snp_number++;
        print "$line\n";
    }
}

while (<INDEL>){
    chomp;
    my $line = $_;
    if (/^#/){
        next;
    }
    else{
        $indel_number++;
        print "$line\n";
    }
}

if ($snp_number == 0 && $indel_number == 0 ){
        print STDERR "Warning：There are no effective mutation sites！\n";
}	
