#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use File::Basename;
use Cwd ('abs_path','getcwd');

my ($help,$sample,$read1,$read2,$path,$conf);

GetOptions (
	"help|h"   =>  \$help,
	"sample|s=s" =>  \$sample,
	"read1|1=s"  =>  \$read1,
	"read2|2=s"  =>  \$read2,
	"path|p=s"   =>  \$path,
	"conf|c=s"   =>  \$conf,
);

my $INFO = <<LINES;
Usage:
	perl $0 -1 R1_fastq -2 R2_fastq -s samplename -c config -p pathway 

Options:
    -help|h     print help information
    -read1|1    whole path of read1
    -read2|2    whole path of read2
    -sample|s   sample name
    -conf|c     config file
    -path|p     analysis output pathway
LINES

die $INFO if ($help);
die $INFO if ( (! $sample) || (! $read1) || (!$read2) || (!$path) ||(! $conf));
$path = abs_path($path);
$read1 = abs_path($read1);
$read2 = abs_path($read2);
#################### software and database ####################################
my $fastqc = &cfg("fastqc");
my $java = &cfg ("java");
my $perl = &cfg("perl");
my $python = &cfg("python");
my $reseqtools = &cfg("reseqtools");
my $trimmomatic = &cfg("trimmomatic");
my $bwa = &cfg("bwa");
my $samtools = &cfg ("samtools");
my $bed = &cfg("bed");
my $ref = &cfg("ref");
my $adapters_path = &cfg("adapters_path");
my $fusion_anno_txt = &cfg("fusion_anno_txt");
my $gatk3 = &cfg("GATK3");
my $qc = &cfg("qc");
my $sam_filter = &cfg("sam_filter");
my $fusion_anno = &cfg("fusion_anno");
my $clear = &cfg("clear");
my $lsf = &cfg("lsf");
my $check = &cfg("check");
##############################################################################
sub cfg {
    my $key_need = shift;
        open CONFIG,'<',"$conf" or die "can not open config file:$!\n";
            while (<CONFIG>){
                chomp;
                if (/^#/){
                        next;
                }
                else{
                        my $key = ((split)[0]);
                        my $value = ((split)[1]);
                        if ($key  eq  $key_need){
                        return $value;
                        }
                }
            }
}

sub dir_structure {
            if (-e "QC_Bam"){
                unless (-e "QC_Bam/Raw"){
                    mkdir "QC_Bam/Raw";
                }
                unless (-e "QC_Bam/Clean"){
                    mkdir "QC_Bam/Clean";
                }
                unless (-e "QC_Bam/Bam"){
                    mkdir "QC_Bam/Bam";
                }
            }
            else{
                mkdir "QC_Bam";
                mkdir "QC_Bam/Raw";
                mkdir "QC_Bam/Clean";
                mkdir "QC_Bam/Bam";
            }
	unless (-e "TMP"){
                mkdir "TMP";
        }
	unless (-e "Run_scripts"){
        	mkdir "Run_scripts";
	}
        unless (-e "results"){
                mkdir "results";
        }
}

sub trim {
        open TRIM,">","$path/Run_scripts/$sample\_trim.sh";
        my $basename_r1 = basename($read1);
        my $basename_r2 = basename($read2);

        unless (-e "$path/QC_Bam/Raw/$basename_r1"){
                system ("ln -s $read1 $path/QC_Bam/Raw/");
        }
        unless (-e "$path/QC_Bam/Raw/$basename_r2"){
                system ("ln -s $read2 $path/QC_Bam/Raw/");
        }
                my $cmd1 = "$fastqc ".
                                "-o $path/TMP/ ".
                                "-t 4 ".
                                "$path/QC_Bam/Raw/$basename_r1 ".
                                "$path/QC_Bam/Raw/$basename_r2 " .
                                "2> $path/TMP/$sample\_fastqc.log";
                my $cmd2 = "$java -jar $trimmomatic PE ".
                                "-phred33 ".
                                "-threads 4 ".
                                "$path/QC_Bam/Raw/$basename_r1 ".
                                "$path/QC_Bam/Raw/$basename_r2 ".
                                "$path/QC_Bam/Clean/$sample\_clean_R1.fq.gz ".
                                "$path/TMP/$sample\_unpaired_R1.fq.gz ".
                                "$path/QC_Bam/Clean/$sample\_clean_R2.fq.gz ".
                                "$path/TMP/$sample\_unpaired_R2.fq.gz ".
                                "ILLUMINACLIP:$adapters_path:2:30:10:1:true LEADING:15 TRAILING:15 SLIDINGWINDOW:30:25 MINLEN:50 ".
                                "2> $path/TMP/$sample\_trim.log";
                my $cmd3 = "$reseqtools Fqtools stat " .
                         "-InFq $path/QC_Bam/Raw/$basename_r1 " .
                         "-InFq $path/QC_Bam/Raw/$basename_r2 " .
                         "-OutStat $path/TMP/$sample\_fqstat.txt ";
                my $cmd4 = "$reseqtools Fqtools stat " .
                         "-InFq $path/QC_Bam/Clean/$sample\_clean_R1.fq.gz " .
                         "-InFq $path/QC_Bam/Clean/$sample\_clean_R2.fq.gz " .
                         "-OutStat $path/TMP/$sample\_clean_fqstat.txt ";
                print ALL "$cmd1\n$cmd2\n$cmd3\n$cmd4\n";
                print TRIM "$cmd1\n$cmd2\n$cmd3\n$cmd4\n";
                close (TRIM);
}
sub bwa{
        open BAM,'>',"$path/Run_scripts/$sample\_bwa.sh";
        my $cmd1 = "$bwa mem ".
                        "-t 4 ".
			"-Y -M ".
                        "-R '\@RG\\tID:$sample\\tSM:$sample\\tPL:illumina\\tLB:library\\tPU:PE\' ".
                        "$ref ".
                        "$path/QC_Bam/Clean/$sample\_clean_R1.fq.gz ".
                        "$path/QC_Bam/Clean/$sample\_clean_R2.fq.gz ".
#                        "| $samtools view -bS -F 4 > $path/QC_Bam/Bam/$sample.bam ".
			"> $path/QC_Bam/Bam/$sample.bam ".
                        "2> $path/TMP/$sample\_bwa.log";
        my $cmd2 = "$samtools sort ".
                        "--threads 5 ".
                        "$path/QC_Bam/Bam/$sample.bam ".
                        "-o $path/QC_Bam/Bam/$sample\_sort.bam ";
        my $cmd3 = "$samtools  index ".
                        "$path/QC_Bam/Bam/$sample\_sort.bam ";
        print BAM "$cmd1\n$cmd2\n$cmd3\n";
        print ALL "$cmd1\n$cmd2\n$cmd3\n";
        close (BAM);
}
sub filter{
	open FILTER,'>',"$path/Run_scripts/$sample\_filter.sh";
	my $cmd1 = "$samtools view -h -F 780 $path/QC_Bam/Bam/$sample\_sort.bam -o $path/QC_Bam/Bam/$sample\_sort.sam";
	my $cmd2 = "$python $sam_filter " .
			"$bed " .
			"$path/QC_Bam/Bam/$sample\_sort.sam > ".
			"$path/TMP/$sample\_reads.stats ";
	my $cmd3 = "$samtools view -b -F 780 -o $path/QC_Bam/Bam/$sample\_sort_filter.bam $path/QC_Bam/Bam/$sample\_sort.sam.filter";
	my $cmd4 = "$samtools index $path/QC_Bam/Bam/$sample\_sort_filter.bam ";
	print FILTER "$cmd1\n$cmd2\n$cmd3\n$cmd4\n";
	print ALL "$cmd1\n$cmd2\n$cmd3\n$cmd4\n";
        close (FILTER);
}

#sub support{
#	open SUPPORT,'>',"$path/Run_scripts/$sample\_support.sh";
#	my $cmd1 = "$java -jar $gatk3 -T DepthOfCoverage " .
#			"-mmq 55 " .
#			"--start 50 " .
#			"--stop 500000 " .
#			"--filter_mismatching_base_and_quals " .
#			"-o $path/TMP/$sample\_DepthOfCoverage " .
#			"-R $ref ".
#			"-I $path/QC_Bam/Bam/$sample\_sort_filter.bam " .
#			"-L $bed ";
#	print SUPPORT "$cmd1\n";
#        print ALL "$cmd1\n";
#        close (SUPPORT);	
#}
sub anno{
	my $cmd1 = "$python $fusion_anno ".
			"$path/TMP/$sample\_reads.stats ".
			"$fusion_anno_txt ".
			"$sample > ".
			"$path/TMP/$sample\_fusion_raw.xls ";
	print ALL "$cmd1\n";
}
sub stat{
	my $cmd1 = "$samtools flagstat " .
                        "$path/QC_Bam/Bam/$sample\_sort.bam " .
                        " > $path/TMP/$sample\_sort_flagstat.txt ";
	my $cmd2 = "$python $qc " .
			"$path/TMP/$sample\_fqstat.txt ".
			"$path/TMP/$sample\_clean_fqstat.txt ".
			"$path/TMP/$sample\_sort_flagstat.txt ".
			"$path/TMP/$sample\_reads.stats ".
			"$sample > ".
			"$path/results/$sample\_QC.xls ";

        my $cmd3 = "$python $check " .
                        "$path/TMP/$sample\_fusion_raw.xls ".
                        "$path/results/$sample\_QC.xls ".
                        "$path/results/$sample\_fusion.xls ";

	print ALL "$cmd1\n$cmd2\n$cmd3\n";
}
sub clear{
#	open CLEAR,'>',"$path/Run_scripts/$sample\_clear.sh";
	my @file = ("$path/Run_scripts/$sample\.sh",
			"$path/QC_Bam/Bam/$sample.bam"
			);
	my $cmd1 = "$perl $clear @file";
	print ALL "$cmd1\n";	
}

sub run{
	if (-e "$path/Run_scripts/$sample\.sh"){
		system ("$python $lsf $path/Run_scripts/$sample\.sh $sample 8");	
	}
}

&dir_structure();
open ALL,">","$path/Run_scripts/$sample\.sh" or die "can not open all.sh file:$!\n";
&trim();
&bwa();
&filter();
#&support();
&anno();
&stat();
#&clear();
&run();
close (ALL);
