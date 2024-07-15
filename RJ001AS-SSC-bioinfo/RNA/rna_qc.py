#!/usr/bin/env python3
from statistics import mean
import os
import sys
import argparse
from decimal import Decimal
import re
import numpy as np

REAL_DIR = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(os.path.dirname(REAL_DIR), "common"))
import tools

def parse_arguments():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("--sample-name", required=True, dest="sample_name", help="sample name")
    parser.add_argument("--fqstat-raw", dest="fqstat_raw", required=True,
                        help="fqstat for raw data")
    parser.add_argument("--fqstat-clean", dest="fqstat_clean", required=True,
                        help="fqstat for clean data")
    parser.add_argument("--flagstat", dest="flagstat",
                        required=True, help="flagstat")
    parser.add_argument("--reads-stat", required=True,
                        dest="reads_stat", help="reads stat")
    parser.add_argument("--out-file", dest="out_file", required=True, help="out file")
    return parser.parse_args()

def main():
    args = parse_arguments()
    sample_name = args.sample_name
    fqstat_raw = args.fqstat_raw
    fqstat_clean = args.fqstat_clean
    flagstat = args.flagstat
    reads_stat = args.reads_stat
    out_file = args.out_file
    qc_summary = RNAQcSummary(fqstat_raw=fqstat_raw, 
                              fqstat_clean=fqstat_clean,
                              flagstat=flagstat,
                              reads_stat=reads_stat,
                              sample_name=sample_name)
    qc_summary.write_qc_summary_file(out_file)

class RNAQcSummary(object):
    headers = ["Sample_Name", "Raw_Data(Mb)", "Raw_Reads", "Clean_Data(Mb)", "Clean_Reads", "Raw_Q30", "Clean_Q30", "GC_Content",
               "Mapped_Num", "Mapped_Ratio", "ITGB7", "TBP", "PUM1", "POP4", "BLOC1S2", "IC_total_reads", "Warning", "Status"]
    gene_exon_list = ('ITGB7_E5-E6', 'TBP_E6-E7', 'PUM1_E21-E22',
                'POP4_E6-E7', 'BLOC1S2_E5-E6')

    warning_threshold = {
        "IC_num": 2,
        "Raw_Q30": 80,
        "Map_Ratio": 90
    }

    def __init__(self, fqstat_raw, fqstat_clean, flagstat, reads_stat, sample_name) -> None:
        self.sample_name = sample_name
        self.raw_q30, self.raw_basenum, self.readnum, self.gccontent = self.parse_fqstat_file(
            fqstat_raw)
        self.clean_q30, self.clean_basenum, self.clean_readnum, self.clean_gccontent = self.parse_fqstat_file(
            fqstat_clean)
        self.map_ratio, self.total_map_num = self.load_flagstat_file(flagstat)
        self.cov_dict = self.inter_ref_gene(reads_stat)
        ic_num = 0
        self.ic_read = 0
        cont_num = 0
        self.status = ""
        self.gene_cov_dict = dict()
        self.warning_msg = []
        for gene_exon in self.gene_exon_list:
            if gene_exon in self.cov_dict:
                if self.clean_basenum == 0:
                    gene_num = int(self.cov_dict[gene_exon])
                else:
                    gene_num = Decimal(
                        int(self.cov_dict[gene_exon]) * 50 / self.clean_basenum).quantize(Decimal("1."))
                self.gene_cov_dict[gene_exon.split("_")[0]]= str(gene_num)
                if int(gene_num) > 0:
                    ic_num += 1
                    self.ic_read += int(gene_num)
                    if int(gene_num) >= 500:
                        cont_num = cont_num + 1
            else:
                self.gene_cov_dict[gene_exon.split("_")[0]] = '.'
        if cont_num < self.warning_threshold["IC_num"]:
            self.warning_msg.append(f"IC_num={cont_num}")
        if self.raw_q30 < self.warning_threshold['Raw_Q30']:
            self.warning_msg.append(f"Raw_Q30={self.raw_q30}")
        if self.map_ratio < self.warning_threshold["Map_Ratio"]:
            self.warning_msg.append(f"Map_Ratio={self.map_ratio}")
        if len(self.warning_msg) == 0:
            self.status = "Pass"
        else:
            self.status = "Fail"

    def inter_ref_gene(self, reads_stat):
        cov_dict = {}
        with open(reads_stat, 'r') as finp:
            for line in finp.readlines():
                line = line.strip("\n")
                gene = line.split("\t")[0]
                cov_depth = line.split("\t")[1]
                if gene in self.gene_exon_list:
                    cov_depth = str(Decimal(cov_depth).quantize(Decimal("1.")))
                    cov_dict[gene] = cov_depth
                else:
                    continue
        return cov_dict
    
    def write_qc_summary_file(self, out_file):
        fout = open(out_file, 'w')
        fout.write("\t".join(self.headers) + "\n")
        summary_dict = dict()
        summary_dict["Sample_Name"] = self.sample_name
        summary_dict["Raw_Data(Mb)"] = self.raw_basenum
        summary_dict["Raw_Reads"] = self.readnum
        summary_dict["Clean_Data(Mb)"] = self.clean_basenum
        summary_dict["Clean_Reads"] = self.clean_readnum
        summary_dict["Raw_Q30"] = self.raw_q30
        summary_dict["Clean_Q30"] = self.clean_q30
        summary_dict["GC_Content"] = self.clean_gccontent
        summary_dict["Mapped_Num"] = self.total_map_num
        summary_dict["Mapped_Ratio"] = self.map_ratio
        summary_dict["ITGB7"] = self.gene_cov_dict["ITGB7"]
        summary_dict["TBP"] = self.gene_cov_dict["TBP"]
        summary_dict["PUM1"] = self.gene_cov_dict["PUM1"]
        summary_dict["POP4"] = self.gene_cov_dict["POP4"]
        summary_dict["BLOC1S2"] = self.gene_cov_dict["BLOC1S2"]
        summary_dict["IC_total_reads"]= self.ic_read
        summary_dict["Warning"] = ";".join(self.warning_msg)
        summary_dict["Status"]=  self.status
        fout.write("\t".join( [str(summary_dict[h]) for h in self.headers]) + "\n")
        fout.close()
    
    @staticmethod
    def format_float(float_digit):
        format_digit = Decimal(float_digit).quantize(Decimal("0.00"))
        return format_digit

    def load_flagstat_file(self, flagstat_file):
        """parse flagstat file, that is generated by command `samtools flagstat`

        Args:
            flagstat_file (string): path to sample's flagstat file

        Returns:
            ratio(float): ratio of mapped reads
            mapped_num(int): number of mapped reads
        """
        total_num = 0
        mapped_num = 0
        with open(flagstat_file, 'r') as finp:
            for line in finp.readlines():
                line = line.strip("\n")
                rest = re.match(
                    r'(\d+)\s\+\s\d+\sin\stotal\s\(QC-passed\sreads\s\+\sQC-failed\sreads\)', line)
                if rest is not None:
                    total_num = int(rest.group(1))
                rest = re.match(r'(\d+)\s\+\s\d+\smapped*', line)
                if rest is not None:
                    mapped_num = int(rest.group(1))
        ratio = self.format_float(round(float(mapped_num * 100) / float(total_num), 2))
        return ratio, mapped_num

    def parse_fqstat_file(self, fqstat_file):
        """parse fqstat file generated by reseqtools to get several metrics value for fastq
            the fqstat file can be for raw fastq or cleaned fastq

        Args:
            fqstat_file (string): path to fqstat 

        Returns:
            q30_mean(decimal): _description_
        """
        q30 = list()
        basenum = list()
        readnum = list()
        gccontent = list()
        with open(fqstat_file, 'r') as finp:
            for line in finp.readlines():
                line = line.strip("\n")
                rest = re.match(r'^#BaseQ:20--30.*>Q30:\s(.*)%', line)
                if rest is not None:
                    q30.append(float(rest.group(1)))
                rest = re.match(r'^#ReadNum:.*BaseNum:\s(\d+)', line)
                if rest is not None:
                    basenum.append(float(rest.group(1)))
                rest = re.match(r'^#ReadNum:\s(\d+)', line)
                if rest is not None:
                    readnum.append(int((rest.group(1))))
                rest = re.match(r'^#GC%:\s(\d+)', line)
                if rest is not None:
                    gccontent.append(float(rest.group(1)))
        q30_mean = self.get_mean(q30)
        basenum_sum = np.sum(basenum)
        basenum_sum = Decimal(basenum_sum/1000000).quantize(Decimal("1."))
        readnum_sum = np.sum(readnum)
        gccontent_mean = self.get_mean(gccontent)
        return q30_mean, basenum_sum, readnum_sum, gccontent_mean

    @staticmethod
    def get_mean(value_list):
        total = sum(value_list)
        length = len(value_list)
        mean  = Decimal(total/length).quantize(Decimal("0.00"))
        return mean



if __name__ == '__main__':
    main()