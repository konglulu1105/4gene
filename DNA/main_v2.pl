#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use File::Basename;
use Cwd ('abs_path','getcwd');

my ($help,$sample,$read1,$read2,$path,$conf,$run);

GetOptions (
    "help|h"    => \$help,
    "sample|s=s"    => \$sample,
    "read1|1=s"    => \$read1,
    "read2|2=s"   => \$read2,
    "path|p=s" => \$path,
    "conf|c=s"    => \$conf,
    "runall|r!" => \$run,
);

my $INFO = <<LINES;
Usage: 
    perl $0 -1 R1_fastq -2 R2_fastq -s samplename -c config -p pathway -r 

Options:
    -help|h     print help information
    -read1|1	whole path of read1
    -read2|2	whole path of read2
    -sample|s     sample name
    -conf|c     config file
    -path|p	analysis output pathway
    -run|r   run all modules
LINES

die $INFO if ($help);
die $INFO if ( (! $sample) || (! $read1) || (!$read2) || (!$path) ||(! $conf));
$path = abs_path($path);
$read1 = abs_path($read1);
$read2 = abs_path($read2);
############# software and database ##################################################
my $fastqc = &cfg("fastqc");
my $java = &cfg ("java");
my $perl = &cfg("perl");
my $python = &cfg("python");
my $rscript = &cfg("Rscript");
my $trimmomatic = &cfg("trimmomatic");
my $reseqtools =  &cfg("reseqtools");
my $bwa = &cfg("bwa");
my $samtools = &cfg ("samtools");
my $picard = &cfg("picard");
my $gatk3 = &cfg("GATK3");
my $gatk4 = &cfg("GATK4");
my $annovar = &cfg("annovar");
my $varscan = &cfg("varscan");
my $vep = &cfg("vep");
my $factera = &cfg("factera");
my $fusion_select = &cfg("fusion_select");
my $variants_combined = &cfg("variants_combined");
my $variants_combined_ignore = &cfg("variants_combined_ignore");
my $vep_tidy = &cfg("vep_tidy");
my $canonical_select = &cfg("canonical_select");
my $filter_black = &cfg("filter_black");
my $bedtools = &cfg("bedtools");
my $variants_select = &cfg("variants_select");
my $ref = &cfg("ref");
my $adapters_path = &cfg("adapters_path");
my $bed = &cfg("bed");
my $annovar_database = &cfg("annovar_database");
my $exons_bed = &cfg("exons_bed");
my $hg19_2bit = &cfg("hg19_2bit");
my $canonical = &cfg("canonical");
my $hotspots = &cfg("hotspots");
my $qc_stat = &cfg("qc_stat");
my $add_info = &cfg("add_info");
my $black_list = &cfg("black_list");
my $tert_anno = &cfg("tert_anno");
my $clear = &cfg("clear");
my $lsf = &cfg("lsf");
my $target_bed = &cfg("target_bed");
my $anno_snv = &cfg("anno_snv");
my $zc_hotspots = &cfg("zc_hotspots");
####################################################################################
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
	    if (-e "Varscan"){
		    unless (-e "Varscan/VCF"){
			mkdir "Varscan/VCF";
		    }
		    unless(-e "Varscan/Anno"){
			mkdir "Varscan/Anno";
		    }
	     }
	    else {
		    mkdir "Varscan";
		    mkdir "Varscan/VCF";
		    mkdir "Varscan/Anno";
	    }
	unless (-e "Run_scripts"){
		mkdir "Run_scripts";
	}
	unless (-e "TMP"){
		mkdir "TMP";
	}
	unless (-e "results"){
		mkdir "results";
	}
	unless (-e "Fusion"){
		mkdir "Fusion";
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
				"-R '\@RG\\tID:$sample\\tSM:$sample\\tPL:illumina\\tLB:library\\tPU:PE\' ".
				"$ref ".
				"$path/QC_Bam/Clean/$sample\_clean_R1.fq.gz ".
				"$path/QC_Bam/Clean/$sample\_clean_R2.fq.gz ". 
#				"| $samtools view -bS -F 4 > $path/QC_Bam/Bam/$sample.bam ".
				"> $path/QC_Bam/Bam/$sample.bam ".
				"2> $path/TMP/$sample\_bwa.log";
        my $cmd2 = "$samtools sort ".
			"--threads 5 ". 
			"$path/QC_Bam/Bam/$sample.bam ".
			"-o $path/QC_Bam/Bam/$sample\_sort.bam";
        my $cmd3 = "$samtools index ".
			"$path/QC_Bam/Bam/$sample\_sort.bam";
	print BAM "$cmd1\n$cmd2\n$cmd3\n";
	print ALL "$cmd1\n$cmd2\n$cmd3\n";
	close (BAM);
}

sub gatk4{
	open GATK4,'>',"$path/Run_scripts/$sample\_gatk4.sh";
	my $cmd1 = "$gatk4 HaplotypeCaller ".
                        "-R $ref ".
                        "-L $target_bed ".
                        "-I $path/QC_Bam/Bam/$sample\_sort.bam ".
                        "-O $path/Varscan/VCF/$sample\_hc.vcf ";
	print GATK4 "$cmd1\n";
	print ALL "$cmd1\n";
	close (GATK4);	
}
sub read_stat {
	open STAT,'>',"$path/Run_scripts/$sample\_stat.sh" or die "$!\n";
	my $cmd1 = "$java -jar $gatk3 " .
			"-T DepthOfCoverage " .
			"--start 50 " .
			"--stop 300000 " .
			"--filter_mismatching_base_and_quals " .
			"-R $ref " .
			"-o $path/TMP/$sample\_sort_DepthOfCoverage " .
			"-I $path/QC_Bam/Bam/$sample\_sort.bam " .
			"-L $bed " .
			"-ct 100 -ct 200 -ct 500 ";
	my $cmd2 = "$java -jar $gatk3 " .
                        "-T DepthOfCoverage " .
			"--start 50 " .
                        "--stop 20000 " .
			"--filter_mismatching_base_and_quals " .
                        "-R $ref " .
                        "-o $path/TMP/$sample\_sort_dedup_DepthOfCoverage " .
                        "-I $path/QC_Bam/Bam/$sample\_sort_dedup.bam " .
                        "-L $bed " .
                        "-ct 5 -ct 50 -ct 500 ";
	my $cmd3 = "$picard CollectInsertSizeMetrics " .
			"-I $path/QC_Bam/Bam/$sample\_sort.bam " .
			"-O $path/TMP/$sample\_insert_size_metrics.txt " .
			"-H $path/TMP/$sample\_insert_size_histogram.pdf " .
			"-M 0.5 ";
	my $cmd4 = "$bedtools intersect " .
			"-wa " .
			"-a $path/QC_Bam/Bam/$sample\_sort.bam " .
			"-b $bed " .
			"> $path/QC_Bam/Bam/$sample\_sort_cover.bam ";
	my $cmd5 = "$java -jar $picard MarkDuplicates ".
                        "I=$path/QC_Bam/Bam/$sample\_sort_cover.bam ".
                        "O=$path/QC_Bam/Bam/$sample\_sort_cover_dedup.bam ".
                        "M=$path/TMP/$sample\_sort_cover_dedup_metrics.txt ".
                        "CREATE_INDEX=true ";
	my $cmd6 = "$samtools flagstat " .
			"$path/QC_Bam/Bam/$sample\_sort.bam " .
			" > $path/TMP/$sample\_sort_flagstat.txt ";
	my $cmd7 = "$samtools flagstat " .
			"$path/QC_Bam/Bam/$sample\_sort_cover.bam " .
			" > $path/TMP/$sample\_sort_cover_flagstat.txt ";
	my $cmd8 = "$perl $qc_stat " .
			"$path/TMP/$sample\_fqstat.txt " .
			"$path/TMP/$sample\_sort_flagstat.txt " .
			"$path/TMP/$sample\_sort_cover_flagstat.txt " .
			"$path/TMP/$sample\_insert_size_metrics.txt " .
#			"$path/TMP/$sample\_sort_cover_dedup_metrics.txt " .
#			"$path/TMP/$sample\_sort_dedup_metrics.txt " .
			"$path/TMP/$sample\_sort_DepthOfCoverage.sample_summary " .
			"$path/TMP/$sample\_clean_fqstat.txt " .
			"$path/TMP/$sample\_sort_DepthOfCoverage ".
			"$sample " .
			"> $path/results/$sample\_QC.xls " ;

	print STAT "$cmd1\n$cmd3\n$cmd4\n$cmd6\n$cmd7\n$cmd8\n";
	print ALL "$cmd1\n$cmd3\n$cmd4\n$cmd6\n$cmd7\n$cmd8\n";
	close (STAT);
}

sub call_varscan {
        open CALL,'>',"$path/Run_scripts/$sample\_varscan.sh";
        my $cmd1 = "$samtools mpileup ".
			"-d 50000 ". 
			"-f $ref ".
			"-l $bed ".
			"$path/QC_Bam/Bam/$sample\_sort_dedup.bam > ". 
			"$path/QC_Bam/Bam/$sample\_sort_dedup.mpileup ";
        my $cmd2 = "$samtools mpileup ".
			"-d 0 ".
			"-A -B -x ".
			"-f $ref ".
			"-l $target_bed ".
			"$path/QC_Bam/Bam/$sample\_sort.bam > ".
			"$path/QC_Bam/Bam/$sample\_sort.mpileup";
	my $cmd3 = "$samtools mpileup ".
			"-x -d 50000 ".
			"-f $ref ".
			"-l $target_bed ".
			"$path/QC_Bam/Bam/$sample\_sort.bam > ".
			"$path/QC_Bam/Bam/$sample\_sort_ignore.mpileup";
        my $cmd4 = "$java -jar $varscan mpileup2snp ".
			"$path/QC_Bam/Bam/$sample\_sort.mpileup ". 
			"--min-coverage 8 ".
			"--min-reads2 1 ". 
			"--min-var-freq 0.0001 ". 
			"--min-freq-for-hom 0.90 ". 
			"--p-value 0.99 " .
			"--strand-filter 0 " .
			"--output-vcf 1 > ". 
			"$path/Varscan/VCF/$sample\_varscan_snp.vcf 2> ". 
			"$path/Varscan/VCF/$sample\_varscan_snp.err ";
	my $cmd5 = "$java -jar $varscan mpileup2snp ".
			"$path/QC_Bam/Bam/$sample\_sort_ignore.mpileup ".
			"--min-coverage 8 ".
			"--min-reads2 1 ". 
                        "--min-var-freq 0.0001 ". 
                        "--min-freq-for-hom 0.90 ". 
                        "--p-value 0.99 " .
                        "--strand-filter 0 " .
                        "--output-vcf 1 > ".
			"$path/Varscan/VCF/$sample\_varscan_snp_ignore.vcf 2> ".
			"$path/Varscan/VCF/$sample\_varscan_snp_ignore.err ";
        my $cmd6 = "$java -jar $varscan mpileup2indel ". 
			"$path/QC_Bam/Bam/$sample\_sort.mpileup ".
			"--min-coverage 8 ".
			"--min-reads2 1 " .
			"--min-var-freq 0.0001 ". 
			"--min-freq-for-hom 0.90 ".
			"--p-value 0.99 ".
			"--strand-filter 0 ".
			"--output-vcf 1 > ".
			"$path/Varscan/VCF/$sample\_varscan_indel.vcf 2> ".
			"$path/Varscan/VCF/$sample\_varscan_indel.err";
	my $cmd7 = "$java -jar $varscan mpileup2indel ".
                        "$path/QC_Bam/Bam/$sample\_sort_ignore.mpileup ".
                        "--min-coverage 8 ".
                        "--min-reads2 1 " .
                        "--min-var-freq 0.0001 ".
                        "--min-freq-for-hom 0.90 ".
                        "--p-value 0.99 ".
                        "--strand-filter 0 ".
                        "--output-vcf 1 > ".
                        "$path/Varscan/VCF/$sample\_varscan_indel_ignore.vcf 2> ".
                        "$path/Varscan/VCF/$sample\_varscan_indel_ignore.err";
        my $cmd8 = "$java -jar $varscan mpileup2snp ".
			"$path/QC_Bam/Bam/$sample\_sort_dedup.mpileup ". 
			"--min-coverage 20 ".
			"--min-reads2 5 ". 
			"--min-var-freq 0.005 ".
			"--min-freq-for-hom 0.90 ". 
			"--p-value 0.99 ".
			"--strand-filter 0 ".
			"--output-vcf 1 > ".
			"$path/Varscan/VCF/$sample\_dedup_varscan_snp.vcf 2> ".
			"$path/Varscan/VCF/$sample\_dedup_varscan_snp.err";
        my $cmd9 = "$java -jar $varscan mpileup2indel ". 
			"$path/QC_Bam/Bam/$sample\_sort_dedup.mpileup ".
			"--min-coverage 20 ".
			"--min-reads2 5 ". 
			"--min-var-freq 0.005 ". 
			"--min-freq-for-hom 0.90 ".
			"--p-value 0.99 ". 
			"--strand-filter 0 ".
			"--output-vcf 1 > ".
			"$path/Varscan/VCF/$sample\_dedup_varscan_indel.vcf 2> ".
			"$path/Varscan/VCF/$sample\_dedup_varscan_indel.err";
        my $cmd10 = "$perl $variants_combined ".
			"$path/Varscan/VCF/$sample\_dedup_varscan_snp.vcf ".
			"$path/Varscan/VCF/$sample\_dedup_varscan_indel.vcf > ".
			"$path/Varscan/VCF/$sample\_dedup_combined.vcf";
        my $cmd11 = "$perl $variants_combined ".
                        "$path/Varscan/VCF/$sample\_varscan_snp.vcf ".
                        "$path/Varscan/VCF/$sample\_varscan_indel.vcf ".
#			"$path/Varscan/VCF/$sample\_hc.vcf ".
                        "> $path/Varscan/VCF/$sample\_combined.vcf";
	my $cmd12 = "$perl $variants_combined_ignore ".
			"$path/Varscan/VCF/$sample\_varscan_snp_ignore.vcf ".
			"$path/Varscan/VCF/$sample\_varscan_indel_ignore.vcf > ".
			"$path/Varscan/VCF/$sample\_combined_ignore.vcf ";
#	print CALL "$cmd2\n$cmd3\n$cmd4\n$cmd5\n$cmd6\n$cmd7\n$cmd11\n$cmd12\n";
#	print ALL "$cmd2\n$cmd3\n$cmd4\n$cmd5\n$cmd6\n$cmd7\n$cmd11\n$cmd12\n";
	print CALL "$cmd2\n$cmd4\n$cmd6\n$cmd11\n";
	print ALL "$cmd2\n$cmd4\n$cmd6\n$cmd11\n";
	close (CALL); 
}

sub anno {
    open ANNO,'>',"$path/Run_scripts/$sample\_anno.sh";
    my $cmd1 = "$perl $vep ".
		"-i $path/Varscan/VCF/$sample\_combined.vcf ".
		"-o  $path/Varscan/Anno/$sample\_vep.txt ".
		"--assembly GRCh37 ".
		"--fasta $ref ".
		"--force_overwrite ". 
		"--cache ".
		"--offline ".
		"--everything ".
		"--refseq ".
		"--canonical";
    my $cmd2 = "$perl $annovar/table_annovar.pl ".
                "$path/Varscan/VCF/$sample\_combined.vcf ".
                "$annovar_database ".
                "-buildver hg19 ".
                "-protocol refGene,snp142,1000g2015aug_all,1000g2015aug_eas,simpleRepeat,rmsk ".
                "-operation g,f,f,f,r,r ".
                "-nastring . ".
                "-remove  ".
                "-vcfinput ".
                "--outfile $path/Varscan/Anno/$sample ";
    my $cmd3 = "$perl $vep_tidy ".
		"$path/Varscan/Anno/$sample\_vep.txt > ".
		"$path/Varscan/Anno/$sample\_vep_tidy.txt";
    my $cmd4 = "$perl $canonical_select ".
		"$path/Varscan/Anno/$sample\_vep_tidy.txt ".
		"$path/Varscan/VCF/$sample\_combined.vcf ".
		"$canonical > ".
		"$path/Varscan/Anno/$sample\_vep_raw.txt";
    my $cmd5 = "$perl $add_info ".
                "$path/Varscan/Anno/$sample\_vep_raw.txt ".
                "$path/Varscan/Anno/$sample.hg19_multianno.txt ".
		"$path/Varscan/VCF/$sample\_combined_ignore.vcf > ".
                "$path/Varscan/Anno/$sample\_vep_paste.txt";
    my $cmd6 = "$python $tert_anno " .
                "$path/Varscan/Anno/$sample\_vep_paste.txt " .
                " > $path/Varscan/Anno/$sample\_vep_total.txt ";
    my $cmd7 = "$perl $variants_select ".
		"$path/Varscan/Anno/$sample\_vep_total.txt ".
		"$hotspots ".
		"$sample > ".
		"$path/Varscan/Anno/$sample\_somatic_raw.xls";
    my $cmd8 = "$perl $filter_black ".
		"$path/Varscan/Anno/$sample\_somatic_raw.xls ".
		"$black_list ".
		"$path/results/$sample\_black.xls ".
		"> $path/results/$sample\_somatic.xls ";
    my $cmd9 = "$perl $anno_snv ".
                "$zc_hotspots " .
                "$path/Varscan/VCF/$sample\_combined.vcf  " .
                "$sample > " .
                "$path/results/$sample\_somatic.xls ";
#    print ANNO "$cmd1\n$cmd2\n$cmd3\n$cmd4\n$cmd5\n$cmd6\n$cmd7\n$cmd8\n";
#    print ALL "$cmd1\n$cmd2\n$cmd3\n$cmd4\n$cmd5\n$cmd6\n$cmd7\n$cmd8\n";
    print  ANNO "$cmd9\n";
    print ALL "$cmd9\n";
    close (ANNO);
}
sub clear{
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
#&gatk4();
&call_varscan();
&anno();
&read_stat();
close (ALL);
if ($run){
	&run();
}
