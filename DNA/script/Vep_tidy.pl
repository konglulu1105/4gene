#!/usr/bin/perl -w
use strict;
use warnings;

die "$0 sample.combined.vepout.txt" if @ARGV < 1 ;

open VEPOUT,'<',"$ARGV[0]" or die "$!";
my $tittle = "Chr\tPosition_Start\tPosition_End\tReference\tMutation\tGene\tHGVSc\tHGVSp\tTranscript\tCanonical\tExon|Intron\tMutation_Type\tdbSNP\tCosmic\t1000G_AF\t1000G_EAS_AF\tgnomAD_AF\tgnomAD_EAS_AF\tSIFT_score\tSIFT\tPolyPhen_score\tPolyPhen\tCLIN_SIG";
print "$tittle\n";
my %AA_hash = (
'Ala' => 'A',
'Arg' => 'R',
'Asp' => 'D',
'Cys' => 'C',
'Gln' => 'Q',
'Glu' => 'E',
'His' => 'H',
'Ile' => 'I',
'Gly' => 'G',
'Asn' => 'N',
'Leu' => 'L',
'Lys' => 'K',
'Met' => 'M',
'Phe' => 'F',
'Pro' => 'P',
'Ser' => 'S',
'Thr' => 'T',
'Trp' => 'W',
'Tyr' => 'Y',
'Val' => 'V'
);
my ($Chr,$Start,$End,$Ref,$Alt,$Transcript_ID,$Mutation_type,$dbSNP,$COSMIC);
my %head;
while (<VEPOUT>){
	chomp;
	my %Extra = ();
        my @dbSNP = ();
        my @COSMIC= ();
	my $line = $_;
	if (/^##/){
		next;	
	}
	elsif (/^#Uploaded_variation/){
		my @tmp = split (/\t/,$line);
		foreach my $value (0..$#tmp){
			$head{$tmp[$value]} = $value;
		}
	}
	else{
		my @element = split (/\t/,$line);
		if (exists $head{'#Uploaded_variation'}){
			my $index = $head{'#Uploaded_variation'};
			my @tmp = split (/_/,$element[$index]);
			$Chr = $tmp[0];
			$Ref = ((split (/\//,$tmp[-1]))[0]);
		}
		if(exists $head{'Location'}){
			my $index = $head{'Location'};
			my @tmp = split (/:/,$element[$index]);
			$Start = ((split (/-/,$tmp[1]))[0]);
			if ((split (/-/,$tmp[1]))[1]){
				$End = ((split (/-/,$tmp[1]))[1]);
			}
			else{
				$End = $Start;
			}
		}
		if(exists $head{'Allele'}){
			my $index = $head{'Allele'};
			$Alt = $element[$index];
		}
		if(exists $head{'Feature'}){
			my $index =  $head{'Feature'};
			$Transcript_ID = $element[$index];
		}
		if(exists $head{'Consequence'}){
			my $index = $head{'Consequence'};
			$Mutation_type = $element[$index];
		}
		if (exists $head{'Existing_variation'}){
			my $index = $head{'Existing_variation'};
			if ($element[$index] eq "-"){
				$dbSNP = ".";
				$COSMIC = ".";
			}
			else{
				my @tmp = split (/,/,$element[$index]);
				foreach (@tmp){
					if (/^rs/){
						push (@dbSNP,$_);
					}
					if (/^COSM/){
						push (@COSMIC,$_);
					}
				}
				$dbSNP = join('|',@dbSNP);
				$COSMIC = join('|',@COSMIC);
			}	
		}
		if (exists $head{'Extra'}){
			my $index = $head{'Extra'};
			my @tmp = split (/;/,$element[$index]);
			foreach my $tmp (@tmp){
				$Extra{((split(/=/,$tmp))[0])} = ((split(/=/,$tmp))[1]);
			}
			unless(exists $Extra{'AF'}){
				$Extra{'AF'} = '.'
			}
			unless(exists $Extra{'EAS_AF'}){
				$Extra{'EAS_AF'} = '.';
			}
			unless(exists$Extra{'gnomAD_AF'}){
				$Extra{'gnomAD_AF'} = '.';
			}
			unless (exists$Extra{'gnomAD_EAS_AF'}){
				$Extra{'gnomAD_EAS_AF'} = '.';
			}
			if (exists $Extra{'SIFT'}){
				$Extra{'SIFT'} =~ /^(\w+)\((\S+)\)/;
				$Extra{'SIFT'} = $1;
				$Extra{'SIFT_score'} = $2;
			}
			else {
				$Extra{'SIFT'} = '.';
				$Extra{'SIFT_score'} = '.';
			}
			if (exists $Extra{'PolyPhen'}){
				$Extra{'PolyPhen'} =~ /^(\w+)\((\S+)\)/;
				$Extra{'PolyPhen'} = $1;	
				$Extra{'PolyPhen_score'}=$2;
			}
			else{
				$Extra{'PolyPhen'} = '.';
				$Extra{'PolyPhen_score'}= '.';
			}
			unless (exists $Extra{'CANONICAL'}){
				$Extra{'CANONICAL'} = '.';
			}
			unless (exists $Extra{'CLIN_SIG'}){
				$Extra{'CLIN_SIG'}= '.';
			}
			if (!exists$Extra{'EXON'} && !exists$Extra{'INTRON'}){
				$Extra{'EXON_INTRON'} = '.';
			}
			elsif (exists$Extra{'EXON'}){
				$Extra{'EXON_INTRON'} = "EXON".$Extra{'EXON'};
			}
			else{
				$Extra{'EXON_INTRON'} = "INTRON".$Extra{'INTRON'};
			}
			unless (exists$Extra{'SYMBOL'}){
				$Extra{'SYMBOL'} = '.';
			}
			if (exists$Extra{'HGVSp'}){
				if ($Extra{'HGVSp'} =~ /^\S+:p\.([A-Z][a-z]{2})[0-9]+(\S{3})$/){
					my $aa_before;
					my $aa_after;
					if (exists $AA_hash{$1}){
						$aa_before = $AA_hash{$1};
					}
					else{
						$aa_before=$1;
					}
					if ($2 eq '%3D'){
						$aa_after = '=';
					}
					elsif(exists$AA_hash{$2}){
						$aa_after = $AA_hash{$2};
					}
					else {
						 $aa_after = $2;	
					}
					$Extra{'HGVSp'} =~ s/^(\S+:p\.)[A-Z][a-z]{2}([0-9]+)\S{3}$/$1$aa_before$2$aa_after/g;
				}
			}
			else{
				$Extra{'HGVSp'} = '.';
			}
			unless (exists$Extra{'HGVSc'}){
				$Extra{'HGVSc'} = '.';
			}								
		}
	print "$Chr\t$Start\t$End\t$Ref\t$Alt\t$Extra{'SYMBOL'}\t$Extra{'HGVSc'}\t$Extra{'HGVSp'}\t$Transcript_ID\t$Extra{'CANONICAL'}\t$Extra{'EXON_INTRON'}\t$Mutation_type\t$dbSNP\t$COSMIC\t$Extra{'AF'}\t$Extra{'EAS_AF'}\t$Extra{'gnomAD_AF'}\t$Extra{'gnomAD_EAS_AF'}\t$Extra{'SIFT_score'}\t$Extra{'SIFT'}\t$Extra{'PolyPhen_score'}\t$Extra{'PolyPhen'}\t$Extra{'CLIN_SIG'}\n";
	}
}
