#!/usr/bin/env python3
import os
import sys
import argparse
from decimal import Decimal
import re

REAL_DIR = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(os.path.dirname(REAL_DIR), "common"))
import tools


def parse_arguments():
    """Create an instance of the ArgumentParser class and assign it to the variable parser
    parses the arguments and returns the parsed arguments.
    """
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("--fusion-anno", dest="fusion_anno", required=True,
                        help="fusion anno")
    parser.add_argument("--fusion-reads", required=True,
                        dest="fusion_reads", help="reads stat")
    parser.add_argument("--sample-name", required=True,
                        dest="sample_name", help="sample name")
    parser.add_argument("--out-file", required=True,
                        dest="out_file", help="out file")
    return parser.parse_args()


def main():
    args = parse_arguments()
    fusion_anno = args.fusion_anno
    fusion_reads = args.fusion_reads
    sample_name = args.sample_name
    out_file = args.out_file
    output_result(fusion_anno, fusion_reads, sample_name, out_file)

def load_fusion_reads(anno_file):
    """take in the path to the annotation file as an argument
    store the information in an dictionary

    Args:
        anno_file (string): path to the annotation file.

    Returns:
        dict[string, list]: dict containing the fusion reads information
    """
    anno_dict = {}
    with open(anno_file) as finp:
        for i in finp.readlines():
            i = i.strip('\n')
            fusion_id = i.split('\t')[0]
            fusion_gene = i.split('\t')[1]
            break_point1 = i.split('\t')[2]
            break_point2 = i.split('\t')[3]
            tran1 = i.split('\t')[4]
            tran2 = i.split('\t')[5]
            anno_dict[fusion_id] = [fusion_gene, tran1,
                                    tran2, break_point1, break_point2]
    return anno_dict


def load_reads_stat_file(reads_stat_file):
    reads_stat_list = list()
    with open(reads_stat_file) as finp:
        for line in finp.readlines():
            line = line.strip('\n')
            entries = line.split("\t")
            fusion_id = entries[0]
            reads_num = int(entries[1])
            reads_stat_list.append((fusion_id, reads_num))
    return reads_stat_list


def output_result(anno_file, reads_stat_file, sample_name, out_file):
    """a function to annotate the reads stat file that contain fusion with prefined annotation file 
    to get the format output file that contains list of annotated fusion.

    Args:
        anno_file (string): path to annotation database file, predefined in the config file
        reads_stat_file (string): read stats got from script <PIPELINE_ROOT>/RNA/sam_filter.py
        sample_name (string): name of sample
        out_file (string): output formated file.
    """
    # Load fusion reads from the annotation file and store them in a dictionary called anno_dict.
    anno_dict = load_fusion_reads(anno_file)
    reads_stat_list = load_reads_stat_file(reads_stat_file)
    # Define the headers for the output file as a list of column names.
    headers = ["Sample_Name", "Gene", "Transcript_1", "Transcript_2",
               "Exon_Rank_1", "Exon_Rank_2", "Support_Reads", "Status"]
    fout = open(out_file, 'w')
    fout.write("\t".join(headers) + "\n")
    # Iterate through each fusion_id and reads_num in the reads_stat_list.
    for fusion_id,  reads_num in reads_stat_list:
        out_dict = dict()
        # check if reads num is greater than or equal to 8 and if fusion_id exists in anno_dict
        if reads_num >= 8 and fusion_id in anno_dict:
            out_dict["Sample_Name"] = sample_name
            out_dict["Gene"] = anno_dict[fusion_id][0]
            out_dict["Transcript_1"] = anno_dict[fusion_id][1]
            out_dict["Transcript_2"] = anno_dict[fusion_id][2]
            out_dict["Exon_Rank_1"] = anno_dict[fusion_id][3]
            out_dict["Exon_Rank_2"] = anno_dict[fusion_id][4]
            out_dict["Support_Reads"] = reads_num
            out_dict["Status"] = "NA"
            fout.write("\t".join([str(out_dict[h]) for h in headers]) + "\n")
    fout.close()


if __name__ == '__main__':
    main()