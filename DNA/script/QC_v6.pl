#!/usr/bin/perl -w
#
use List::Util qw/sum/;
die "$0 fqstat.txt flagstat.txt sort_cover_flagstat.lxt insert_size_metrics.txt  sort_DepthOfCoverage.sample_summary clean_fqstat.txt sort_DepthOfCoverage samplename" if @ARGV  < 6;

#open FQSTAT,'<',"$ARGV[0]" or die "$!";
#open FLAGS,'<',"$ARGV[1]" or die "$!";
#open COVERFLAGS,'<',"$ARGV[2]" or die "$!";
open INSERTSIZE,'<',"$ARGV[3]" or die "$!";
#open COVERMET,'<',"$ARGV[4]" or die "$!";
#open MET,'<',"$ARGV[5]" or die "$!";
open DEPTH,'<',"$ARGV[4]" or die "$!";

my $sample = $ARGV[7];
my $raw_data_mb;
my $clean_data_mb;
my $raw_q30_ratio;
my $clean_q30_ratio;
#my $map_ratio;
my $insert_bp;
#my $covered_complex;
#my $lib_complex;
my $capture_ratio;
my $mean_depth;
my $median_depth;
my $depth50;
my $raw_reads;
my $clean_reads;
my $raw_gc_content;
my $clean_gc_content;
###############################################
my ($raw_q30,$raw_basenum,$readnum,$gccontent) = extract($ARGV[0]);
my ($clean_q30,$clean_basenum,$clean_readnum,$clean_gccontent) = extract($ARGV[5]);
sub extract{
        my $fqstat_file = shift;
        my @q30;
        my @basenum;
		my @readnum;
		my @gccontent;
        open FQSTAT,'<',"$fqstat_file" or die "$!";
        while (<FQSTAT>){
                chomp;
                my $line = $_;
                if (/^#BaseQ:20--30.*>Q30:\s(.*)%/){
            	    push @q30,$1;
                }
                if (/^#ReadNum:.*BaseNum:\s(\d+)/){
                    push @basenum,$1;
                }
				if (/^#ReadNum:\s(\d+)/){
					push @readnum,$1;
				}
				if (/^#GC%:\s(\d+)/){
					push @gccontent,$1; 
				}
        }
        return (\@q30,\@basenum,\@readnum,\@gccontent);
}
###############################################
sub average{
        my $aver = 0;
        map{$aver += $_} @_;
        $aver /=($#_+1)
}
my $tmp0 = average(@$raw_q30);
my $tmp00 = average(@$clean_q30);
$raw_q30_ratio = sprintf "%.2f",$tmp0 ;
$clean_q30_ratio = sprintf "%.2f",$tmp00 ;
my $tmp1 = (sum (@$raw_basenum))/1000000 ;
my $tmp11 = (sum (@$clean_basenum))/1000000 ;
$raw_data_mb = sprintf "%.0f",$tmp1 ;
$clean_data_mb = sprintf "%.0f",$tmp11 ;
my $tmp2 = sum(@$readnum);
my $tmp22 = sum(@$clean_readnum);
my $tmp222 = (sum (@$clean_readnum))/1000000 ;
$raw_reads = sprintf "%.0f",$tmp2 ;
$clean_reads = sprintf "%.0f",$tmp22;
$clean_reads_mb = sprintf "%.0f",$tmp222;
my $tmp3 = average(@$gccontent);
my $tmp33 = average(@$clean_gccontent);
$raw_gc_content = sprintf "%.2f",$tmp3;
$clean_gc_content = sprintf "%.2f",$tmp33;
###############################################
sub mapreads{
	my $flags = shift;
	my $total_num;
#	my $mate_num;
#	my $singleton_num;
	my $mapped_num;
	open FLAGS,'<',"$flags" or die "$!"; 
	while (<FLAGS>){
        	chomp;
        	my $line = $_;
		if (/(\d+)\s\+\s\d+\sin\stotal\s\(QC-passed\sreads\s\+\sQC-failed\sreads\)/){
			$total_num = $1;	
		}
#		elsif (/(\d+)\s\+\s\d+\swith\sitself\sand\smate\smapped/){
#			$mate_num = $1;
#		}
#		elsif (/(\d+)\s\+\s\d+\ssingletons*/){
#			$singleton_num = $1;
#		}
		elsif (/(\d+)\s\+\s\d+\smapped*/){
			$mapped_num = $1;	
		}
		else{
			next;
		}
	}		
#	my $total_map_num = $mate_num + $singleton_num ; 
	my $ratio = sprintf "%.2f",($mapped_num * 100/$total_num);
	return ($ratio,$mapped_num);
	close (FLAGS);
}
my ($map_ratio,$total_map_num) = mapreads($ARGV[1]);
my $capture_total_map_num = (mapreads($ARGV[2]))[1];
my $tmp_capture_ratio = ($capture_total_map_num/$total_map_num) * 100 ;
$capture_ratio =  sprintf "%.2f",$tmp_capture_ratio ;
##############insertsize##########################################
my @insertsize = <INSERTSIZE>;
my $tmp_insert_bp = (split /\t/,$insertsize[7])[0];
$insert_bp = sprintf "%.f",$tmp_insert_bp ;
############Duplicate#############################################
#my @duplicate = <MET>;
#my $tmp4 = ((split /\t/,$duplicate[7])[7]);
#$lib_complex = sprintf "%.3f",$tmp4 ; 
#################################################################
#my @cover_dup = <COVERMET>;
#my $tmp5 = ((split /\t/,$cover_dup[7])[7]); 
#$covered_complex = sprintf "%.3f",$tmp5 ;
#################################################################
my @depth = <DEPTH>;
my $tmp6 = (split /\t/,$depth[1])[2];
$mean_depth = sprintf "%.f",$tmp6;
my $tmp7 = (split /\t/,$depth[1])[4];
$median_depth = sprintf "%.f",$tmp7;
my $tmp8 = ((split /\t/,$depth[1])[8])/100;
$depth50 = sprintf "%.3f",$tmp8;
##################################################################
sub greater{
	my $mean_depth = shift;
	my $depth_base = shift;
	my $line_depth = shift;
#	my $pct20_depth = 0.2 * $mean_depth;
	my $n_total = 0;
	my $n = 0;
	open DEPTH_BASE,'<',"$depth_base" or die "$!\n";
	while (<DEPTH_BASE>){
		chomp;
		my $line = $_;
		if ($line =~ /^Locus/){
			next;
		}
		else{
			$n_total++;
			my $depth_sample = (split /\t/,$line)[3];
			if ($depth_sample >= $line_depth){
				$n++;
			}
		}
	}
	close (DEPTH_BASE);
	my $tmp = ($n/$n_total)*100;
	my $greater = sprintf "%.2f",$tmp;
	return ($greater);
}
my $greater_1x = greater($mean_depth,$ARGV[6],1);
my $greater_pct20 = greater($mean_depth,$ARGV[6],0.2 * $mean_depth);
my $greater_min = greater($mean_depth,$ARGV[6],2000);
my $greater_mean = greater($mean_depth,$ARGV[6],$mean_depth);
my @warning;
my $status;
my $join_warn;
#if ($clean_reads_mb < 0.45){
#	push @warning,"CleanReads.Mb=$clean_reads_mb";
#}
if ($raw_q30_ratio < 80){
	push @warning,"Raw.Q30.Ratio=$raw_q30_ratio";
}
if ($map_ratio < 90){
	push @warning,"Map.Ratio=$map_ratio";
}
#if ($capture_ratio < 80){
#	push @warning,"On_target_ratio=$capture_ratio";
#}
if ($mean_depth < 20000){
	push @warning,"Mean.Depth=$mean_depth";
}
if ($greater_pct20 < 70){
	push @warning,"Greater_20Pct=$greater_pct20";
}
if (@warning >= 1){
	$join_warn = join ';',@warning;
	$status = 'Fail';
}
else{
	$join_warn = 'NA';
	$status = 'Pass';
}
#print "Sample.name\tRawData.Mb\tCleanData.Mb\tRaw.Q30.Ratio\tClean.Q30.Ratio\tMap.Ratio\tInsert.Bp\tOn_target_ratio\tMean.Depth\tMedian.Depth\tDep500X.Ratio\t>=0.2X\twarning\tstatus\n";
print "Sample_Name\tRaw_Data(Mb)\tRaw_Reads\tClean_Data(Mb)\tClean_Reads\tRaw_Q30\tClean_Q30\tGC_Content\tMapped_Ratio\tInsert_Size\tOn_Target_Ratio\tMean_Depth\tMedian_Depth\tGreater_1X\tGreater_Min\tGreater_20Pct\tGreater_Mean\tWarning\tStatus\n";
#print "$sample\t$raw_data_mb\t$clean_data_mb\t$raw_q30_ratio\t$clean_q30_ratio\t$map_ratio\t$insert_bp\t$capture_ratio\t$mean_depth\t$median_depth\t$depth50\t$pct20_depth\t$join_warn\t$status\n";
print "$sample\t$raw_data_mb\t$raw_reads\t$clean_data_mb\t$clean_reads\t$raw_q30_ratio\t$clean_q30_ratio\t$clean_gc_content\t$map_ratio\t$insert_bp\t$capture_ratio\t$mean_depth\t$median_depth\t$greater_1x\t$greater_min\t$greater_pct20\t$greater_mean\t$join_warn\t$status\n";
