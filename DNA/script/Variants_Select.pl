#!/usr/bin/perl -w
use warnings;

die "$0 vep.txt hotspots.txt samplename" if @ARGV < 3 ;

open VEP,'<',"$ARGV[0]" or die "$!";
open HOTSPOT,'<',"$ARGV[1]" or die "$!";

my @hotspot;

my @low_type = ("splice_region_variant",
		"intron_variant",
		"start_retained_variant",
		"stop_retained_variant",
		"synonymous_variant",
		"3_prime_UTR_variant",
		"5_prime_UTR_variant",	
		"upstream_gene_variant",
		"downstream_gene_variant",
		"intergenic_variant",
		"non_coding_transcript_variant");
while (<HOTSPOT>){
	chomp;
	my $line = $_;
	unless (/^#/){
		my @tmp = split /\t/,$line;
		my $chr = $tmp[0];
        	my $start = $tmp[1];
        	my $ref = $tmp[3];
        	my $alt = $tmp[4];
		if (length($ref) == length($alt)){
            		my $key = join "\t",($chr,$start,$ref,$alt);
			push @hotspot,$key;
        	}
        	elsif(length($ref) < length($alt)){
            		$ref = "-";
            		$alt = substr ($alt,1);
            		my $key = join "\t",($chr,$start,$ref,$alt);
			push @hotspot,$key;
        	}
        	else{
            		$ref = substr ($ref,1);
            		$alt = "-";
            		$start = $start + 1;
            		my $key = join "\t",($chr,$start,$ref,$alt);
			push @hotspot,$key;
        	}

	}
}
close (HOTSPOT);
LABEL1: while (<VEP>){
	chomp;
        my $line = $_;
	if (/^Chr/){
		print "SampleName\t$line\n";
	}
	else{
		my @tmp = split /\t/,$line;
                my $key = join "\t",@tmp[0,1,3,4];
		my $depth = $tmp[6];
		my $variant_supporting = $tmp[7];
		my $Allele_Frequency = $tmp[8];
		my $gene = $tmp[9];
		chop $Allele_Frequency;
		my $af = $tmp[18];
		my $eas_af = $tmp[19];
		my $gnomAD_AF = $tmp[20];
		my $gnomAD_EAS_AF = $tmp[21];
		my @Mutation_Type = split /,/,$tmp[15];
		foreach my $hotspot (@hotspot){
			if ($hotspot eq $key && $depth >= 50 && $variant_supporting >= 8 && $Allele_Frequency >= 0.1){
				print "$ARGV[2]\t$line\n";
				next LABEL1;
			}
		}
		if($depth >= 50 && $variant_supporting >= 10 && $Allele_Frequency >= 0.1 && ($eas_af eq '.' || $eas_af <= 0.001) && ($gnomAD_EAS_AF eq '.' || $gnomAD_EAS_AF <= 0.001) && ($af eq '.' || $af <= 0.001) && ($gnomAD_AF eq '.' || $gnomAD_AF <= 0.001)){
#			if ($gene eq 'TERT' && grep /upstream_gene_variant/,@Mutation_Type){
#				print "$ARGV[2]\t$line\n";
#			}
#			else{
				foreach my $type (@Mutation_Type){
					unless(grep /^$type$/,@low_type){
						print "$ARGV[2]\t$line\n";
						next LABEL1;
					}
				}
#			}
		}
	}
}
