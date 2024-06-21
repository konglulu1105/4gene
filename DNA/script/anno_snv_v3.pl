#!/usr/bin/perl -w
use strict;
use warnings;

die "$0  ZC_hotspot.txt sample_combined_ignore.vcf samplename" if @ARGV < 3 ;


open VCF,'<',"$ARGV[1]" or die "$!\n";

sub snv_anno{
    my $snv_id = shift;
    my @return_part1;
    my @return_part2;
    open ANNO,'<',"$ARGV[0]" or die "$!\n";
    while (<ANNO>){
        chomp;
        my $line = $_;
        my @tmp_part1 = (split /\t/,$line)[0,1,3,4];
        my @tmp_part1_re = (split /\t/,$line)[0..4];
        my @tmp_part2 = (split /\t/,$line)[5..11];
        if (@tmp_part1 ~~ @$snv_id){
            @return_part1 = @tmp_part1_re;
            @return_part2 = @tmp_part2;
            last;
#            return (\@tmp_part1,\@tmp_part2);
        }
        else{
            next;
        }
    }
    return (\@return_part1,\@return_part2);
}
print "Sample_Name\tChr\tPosition_Start\tPosition_End\tReference\tMutation\tDP\tAD\tMutation_Frequency\tGene\tHGVSc\tHGVSp\tExon\tMutation_Type\tdbsnp\tcosmic\tStatus\n";
while (<VCF>){
    chomp;
    my $line =$_;
    if ($line =~ /^#/){
        next;
    }
    else{
        my $status;
        my @tmp_part1 = (split /\t/,$line)[0,1,3,4];
        my ($anno_part1,$anno_part2)= &snv_anno(\@tmp_part1);
        if (@$anno_part1 && @$anno_part2){
            my $startP = @$anno_part1[1];
            my $ref = @$anno_part1[3];
            my $alt = @$anno_part1[4];
            my $tmp_7 = (split /\t/,$line)[7];
            my $tmp_9 = (split /\t/,$line)[9];
            my $adp = (split /;/,$tmp_7)[0];
            $adp =~ tr/ADP=//d;
            my $ad = (split /:/,$tmp_9)[5];
            my $freq = (split /:/,$tmp_9)[6];
            (my $freq_tmp = $freq)=~ s/%//g;
            if ($startP == 140453136 && $ref eq 'A' && $alt eq 'T' && $freq_tmp >= 0.6){
            $status = "Positive";
            }
            elsif ($startP == 1295250 && $ref eq 'G' && $alt eq 'A' && $freq_tmp >= 0.7){
            $status = "Positive";
            }
            elsif ($startP == 1295228 && $ref eq 'G' && $alt eq 'A' && $freq_tmp >= 0.8){
            $status = "Positive";
            }

     #       if ($startP == 105246551 && $ref eq 'C' && $alt eq 'T' && $freq_tmp >= 0.5){
     #       $status = "Positive";
     #       }
     #       elsif ($startP == 140453136 && $ref eq 'A' && $alt eq 'T' && $freq_tmp >= 0.7){
     #       $status = "Positive";
     #       }
     #       elsif ($startP == 178952085 && $ref eq 'A' && $alt eq 'G' && $freq_tmp >= 0.5){
     #       $status = "Positive";
     #       }
     #       elsif ($startP == 178952085 && $ref eq 'A' && $alt eq 'T' && $freq_tmp >= 0.7){
     #       $status = "Positive";
     #       }
     #       elsif ($startP == 1295250 && $ref eq 'G' && $alt eq 'A' && $freq_tmp >= 0.9){
    #        $status = "Positive";
   #         }
  #          elsif ($startP == 1295228 && $ref eq 'G' && $alt eq 'A' && $freq_tmp >= 0.4){
 #           $status = "Positive";
#            }
            else{$status = "Negative";}
            my $output_line = (join "\t", @$anno_part1) . "\t" . (join "\t",($adp,$ad,$freq)) . "\t" . (join "\t",@$anno_part2);
#            if ($status eq "Positive"){
            print "$ARGV[2]\t$output_line\t$status\n";
#            }
        }
        else{
            next;
        }
    }
}
