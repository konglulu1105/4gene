#!/usr/bin/perl -w
use strict;
use warnings;

die "$0 sample_varscan_vep.txt sample_varscan.vcf Canonical.txt" if @ARGV < 3 ;

open VEP,'<',"$ARGV[0]" or die "$!";
open VCF,'<',"$ARGV[1]" or die "$!";
open CANO,'<',"$ARGV[2]" or die "$!";
my @canonical;
my %vcf;
my %vep;
my %trans_count;

while (<CANO>){
    chomp;
    push @canonical,$_;
}
close (CANO);
while (<VCF>){
    chomp;
    unless (/^#/){
        my $line = $_;
        my @tmp = split /\t/,$line;
        my $chr = $tmp[0];
        my $start = $tmp[1];
        my $ref = $tmp[3];
        my $alt = $tmp[4];
        my @info = split /:/,$tmp[9];
        if (length($ref) == length($alt)){
            my $key = join "\t",@tmp[0,1,3,4];
	        if (@info == 14){
                my $value = join "\t",@info[0,3,5,6];
		        $vcf{$key} = $value;
	        }
            if (@info == 5){
                my $gt = $info[0];
                my $dp = $info[2];
                my $ad = (split /\,/,$info[1])[1];
		my $tmp = ($ad/$dp) * 100;
                my $freq =  (sprintf "%.2f",$tmp) . '%';
                my $value = join "\t",($gt,$dp,$ad,$freq);
                $vcf{$key} = $value;
            }
        }
        elsif(length($ref) < length($alt)){
            $ref = "-";
            $alt = substr ($alt,1);
            my $key = join "\t",($chr,$start,$ref,$alt);
            if (@info == 14){
                my $value = join "\t",@info[0,3,5,6];
                $vcf{$key} = $value;
            }
            if (@info == 5){
                my $gt = $info[0];
                my $dp = $info[2];
                my $ad = (split /\,/,$info[1])[1];
		my $tmp = ($ad/$dp) * 100;
                my $freq = (sprintf "%.2f",$tmp) . '%';
                my $value = join "\t",($gt,$dp,$ad,$freq);
                $vcf{$key} = $value;
            }
        }    
        else{
            $ref = substr ($ref,1);
            $alt = "-";
	        $start = $start + 1;
            my $key = join "\t",($chr,$start,$ref,$alt);
            if (@info == 14){
                my $value = join "\t",@info[0,3,5,6];
                $vcf{$key} = $value;
            }
            if (@info == 5){
                my $gt = $info[0];
                my $dp = $info[2];
                my $ad = (split /\,/,$info[1])[1];
		my $tmp = ($ad/$dp) * 100;
                my $freq = (sprintf "%.2f",$tmp) . '%';
                my $value = join "\t",($gt,$dp,$ad,$freq);
                $vcf{$key} = $value;
            }
        }
    }
}
close (VCF);
while (<VEP>){
    chomp;
    my @tmp = split /\t/,$_;
    my $head = join "\t",@tmp[0..4];
    my $tail = join "\t",@tmp[5..$#tmp];
    if (/^Chr/){
        print "$head\t"."GT\tDP\tAD\tMutation_Frequency\t"."$tail\n";
    }
    else{
        my $transcript = (split /\./,$tmp[8])[0];
        if ($transcript ne '' && grep /^$transcript$/,@canonical){
            my $key_vcf = join "\t",@tmp[0,1,3,4];
            my $key_vep = join ("\t",@tmp[0,1,3,4],$transcript);
            $vep{$key_vep} = join ("\t",$head,$vcf{$key_vcf},$tail);
        }
        else{
            next;
        }
    }
}
close (VEP);
foreach my $value (values %vep){
    print "$value\n";
}
