#!/nextcode/sge_software/anaconda2/bin/perl -w
use warnings;

die "$0 vep.vcf annovar.txt combined.vcf" if @ARGV < 3 ;

open VEP,'<',"$ARGV[0]" or die "$!";
open ANNOVAR,'<',"$ARGV[1]" or die "$!";
open VCF,'<',"$ARGV[2]" or die "$!";

my %annovar;
my %vcf;
my $super_dups_site;
my $rmsk;
while (<ANNOVAR>){
	chomp;
	my $line = $_;
	if (/Chr/){
		my @tmp = split /\t/,$line;
		($super_dups_site) = grep{$tmp[$_] eq 'simpleRepeat'} 0..$#tmp;
		($rmsk) = grep{$tmp[$_] eq 'rmsk'} 0..$#tmp;
	}
	else{
		my @tmp = split /\t/,$line;
		my $key = join "\t",@tmp[0,1,3,4];
		$annovar{$key} = $tmp[$super_dups_site] . '|' . $tmp[$rmsk];
	}
}
close (ANNOVAR);

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
		my $adf = $info[12];
                my $adr = $info[13];
                my $strand = sprintf "%.2f",(($adf * 100)/($adf + $adr));
#                my $value = join "\t",(@info[0,3,5,6],$strand);
                        $vcf{$key} = $strand;
                }
            	if (@info == 5){
#                my $gt = $info[0];
#               my $dp = $info[2];
		my $strand = '.';
#                my $ad = (split /\,/,$info[1])[1];
#                my $tmp = ($ad/$dp) * 100;
#                my $freq =  (sprintf "%.2f",$tmp) . '%';
#                my $value = join "\t",($gt,$dp,$ad,$freq,$strand);
                $vcf{$key} = $strand;
            }
        }
        elsif(length($ref) < length($alt)){
            $ref = "-";
            $alt = substr ($alt,1);
            my $key = join "\t",($chr,$start,$ref,$alt);
            if (@info == 14){
		my $adf = $info[12];
                my $adr = $info[13];
                my $strand = sprintf "%.2f",(($adf * 100)/($adf + $adr));
#                my $value = join "\t",(@info[0,3,5,6],$strand);
                $vcf{$key} = $strand;
            }
            if (@info == 5){
#                my $gt = $info[0];
#                my $dp = $info[2];
		my $strand = '.';
#                my $ad = (split /\,/,$info[1])[1];
#                my $tmp = ($ad/$dp) * 100;
#                my $freq = (sprintf "%.2f",$tmp) . '%';
#                my $value = join "\t",($gt,$dp,$ad,$freq,$strand);
                $vcf{$key} = $strand;
            }
        }
        else{
            $ref = substr ($ref,1);
            $alt = "-";
                $start = $start + 1;
            my $key = join "\t",($chr,$start,$ref,$alt);
            if (@info == 14){
		my $adf = $info[12];
                my $adr = $info[13];
                my $strand = sprintf "%.2f",(($adf * 100)/($adf + $adr));
#                my $value = join "\t",(@info[0,3,5,6],$strand);
                $vcf{$key} = $strand;
            }
            if (@info == 5){
#                my $gt = $info[0];
#                my $dp = $info[2];
		my $strand = '.';
#                my $ad = (split /\,/,$info[1])[1];
#                my $tmp = ($ad/$dp) * 100;
#                my $freq = (sprintf "%.2f",$tmp) . '%';
#                my $value = join "\t",($gt,$dp,$ad,$freq,$strand);
                $vcf{$key} = $strand;
            }
        }
    }
}
close (VCF);

while (<VEP>){
        chomp;
        my $line = $_;
        if (/Chr/){
                print "$line\tRepeat\tStrand\n";
        }
        else{
                my @tmp = split /\t/,$line;
                my $key = join "\t",@tmp[0,1,3,4];
                if (grep /^$key$/,keys %annovar ){
			if (grep /^$key$/,keys %vcf){
				print "$line\t$annovar{$key}\t$vcf{$key}\n";
			}
			else{
				print "$line\t$annovar{$key}"."\t\." ."\n";
			}
		}
		else{
			if (grep /^$key$/,keys %vcf){
				print "$line\t\.\t$vcf{$key}\n";
			}
			else{
				print "$line" . "\t\." ."\n";
			}
		}
	
        }
}
